import numpy as np
import os
import rasterio
from rasterio.windows import from_bounds,Window
from tqdm import trange
from rasterio.merge import merge
import glob
import ctypes
import math
from pathlib import Path
try:
    from .def_resize import resize_clip
except ImportError:  # Allow running this file directly.
    from def_resize import resize_clip

MEMORY_LIMIT_RATIO = 0.8
WORKING_ARRAY_COUNT = 12


def _get_memory_info():
    """Return (total_bytes, available_bytes, used_ratio)."""
    try:
        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        memory_status = MEMORYSTATUSEX()
        memory_status.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory_status))
        total = int(memory_status.ullTotalPhys)
        available = int(memory_status.ullAvailPhys)
        used_ratio = float(memory_status.dwMemoryLoad) / 100
        return total, available, used_ratio
    except Exception:
        pass

    try:
        import psutil

        memory_status = psutil.virtual_memory()
        return int(memory_status.total), int(memory_status.available), float(memory_status.percent) / 100
    except Exception:
        return None, None, None


def _estimate_calculation_bytes(width, height, datasets):
    dtype_size = 8
    for dataset in datasets:
        dtype_size = max(dtype_size, np.dtype(dataset.dtypes[0]).itemsize)
    return int(width * height * dtype_size * WORKING_ARRAY_COUNT)


def _format_gb(byte_count):
    if byte_count is None:
        return "unknown"
    return f"{byte_count / 1024 ** 3:.2f} GB"


def _resolve_paths(fs_path, output_path):
    if output_path.lower().endswith((".tif", ".tiff")):
        result_path = output_path
        output_dir = os.path.dirname(output_path)
    else:
        output_dir = output_path
        result_path = os.path.join(output_dir, "hw.tif")

    dem_candidates = [
        os.path.join(output_dir, "dem.tif"),
        os.path.join(fs_path, "dem.tif"),
    ]
    dem_path = next((path for path in dem_candidates if os.path.exists(path)), dem_candidates[0])
    return dem_path, result_path


def _read_inputs(srcs, window=None, b=30):
    uca = srcs["uca"].read(1, window=window) + 1
    slope = srcs["slope"].read(1, window=window)
    soil_depth = srcs["soil_depth"].read(1, window=window)
    ks = srcs["ks"].read(1, window=window)
    pr = srcs["rainfall"].read(1, window=window) / 1000 / 365
    eva = srcs["Eva"].read(1, window=window) / 1000 / 365
    ks = np.power(10, ks) / 100
    return uca, slope, soil_depth, ks, pr, eva


def _iter_block_windows(row_off, col_off, height, width, max_pixels):
    block_size = max(1, int(math.sqrt(max_pixels)))
    row_end = row_off + height
    col_end = col_off + width

    for row_s in range(row_off, row_end, block_size):
        for col_s in range(col_off, col_end, block_size):
            row_e = min(row_s + block_size, row_end)
            col_e = min(col_s + block_size, col_end)
            yield Window(col_s, row_s, col_e - col_s, row_e - row_s)


def _calculate_block_pixels(available_memory, datasets):
    dtype_size = 8
    for dataset in datasets:
        dtype_size = max(dtype_size, np.dtype(dataset.dtypes[0]).itemsize)

    bytes_per_pixel = dtype_size * WORKING_ARRAY_COUNT
    if available_memory is None:
        target_bytes = 512 * 1024 ** 2
    else:
        target_bytes = min(max(available_memory * 0.25, 32 * 1024 ** 2), 512 * 1024 ** 2)
    return max(1, int(target_bytes / bytes_per_pixel))

# 计算初始地下水位，自动检测内存不足时的块大小，并进行分块计算

def Cal_water_table(uca, slope, soil_depth, ks, pr, eva, b=30):
    """
    直接使用数组计算稳态地下水位深度

    参数:
        uca_array: 汇流累积面积数组 dem
        slope_array: 坡度数组
        soil_thickness_array: 土壤厚度数组 1m.2m,3m
        ks_array: 饱和导水率数组  #
        pr_array: 前期降雨数组  cm/day
        eva_array: 蒸散发数组  cm/day
        b: 土壤单元宽度，默认30米

    返回:
        M_staz: 稳态地下水位深度数组
    注意：输入的降雨和蒸散发数据单位保持统一，如mm/h或mm/day
    """
    # 确保坡度为弧度制
    if np.nanmax(np.abs(slope)) > 3:
        slope_radians = np.radians(slope)
    else:
        slope_radians = slope

    # eff_infl = 0.8*(pr - eva)
    eff_infl = 0.8 * pr
    eff_infl = np.nan_to_num(eff_infl, nan=0, posinf=0, neginf=0)

    denominator = ks * b
    denominator = np.nan_to_num(denominator, nan=1e-6, posinf=1e6, neginf=1e-6)

    # 防止除零错误
    safe_denominator = np.where(np.abs(denominator) < 1e-6, 1e-6, denominator)

    # 计算Ikb
    Ikb = np.divide(eff_infl, safe_denominator)
    Ikb = np.nan_to_num(Ikb, nan=0, posinf=0, neginf=0)

    # 防止土壤厚度为零
    safe_soil_thickness = np.where(np.abs(soil_depth) < 1e-6, 1e-6, soil_depth)
    soil_invalid = soil_depth < 1

    # 土壤厚度经验赋值
    soil_invalid = soil_depth < 1
    safe_soil_thickness[soil_invalid] = 1

    # 防止坡度计算错误
    cos_slope = np.cos(slope_radians)
    sin_slope = np.sin(slope_radians)
    cos_slope = np.where(np.abs(cos_slope) < 1e-6, 1e-6, cos_slope)
    sin_slope = np.where(np.abs(sin_slope) < 1e-6, 1e-6, sin_slope)

    # 计算M_staz
    eps = 0.0001
    M_staz = (Ikb * uca * b * b) / ((safe_soil_thickness * cos_slope * sin_slope) + eps)  #

    # 数值异常处理
    # 土壤厚度为0就设置hw为0
    # ks为nodata就设置hw为1
    ks_invalid = np.isnan(ks)  # ks 是 NaN 的地方
    M_staz[ks_invalid] = 1
    soil_invalid = np.isnan(soil_depth)  # ks 是 NaN 的地方
    M_staz[soil_invalid] = 1
    # 处理数值异常
    M_staz = np.nan_to_num(M_staz, nan=0, posinf=1, neginf=0)
    # 限制M_staz的范围在0-1之间
    M_staz[ks_invalid] = 1
    M_staz[soil_invalid] = 1
    M_staz = np.nan_to_num(M_staz, nan=0, posinf=1, neginf=0)
    M_staz = np.clip(M_staz, 0, 1)
    # 对于坡度过小的地方（小于约2.87度），设为0
    M_staz[slope_radians < 0.05] = 0

    return M_staz

def merge_blocks(temp_path, output_path):
    """
    将分块的tif拼接为一个完整tif
    参数:
        temp_path (str): 存放分块结果的目录
        output_path (str): 输出文件路径
    """
    # 找到所有分块tif
    search_path = os.path.join(temp_path, "hw_*.tif")
    tif_files = sorted(glob.glob(search_path))

    if len(tif_files) == 0:
        raise FileNotFoundError(f"在 {temp_path} 中没有找到 hw_*.tif 文件")

    # 打开所有分块
    src_files_to_mosaic = [rasterio.open(fp) for fp in tif_files]
    mosaic, out_trans = merge(src_files_to_mosaic)

    # 使用第一个分块的元信息作为模板
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "dtype": "float32",
        "count": 1
    })
    # 写出结果
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic.astype("float32"))
    for src in src_files_to_mosaic:
        src.close()
    print(f"拼接完成，结果保存到 {output_path}")


# def calculate_initial_hw_old(fs_path, output_path,b=30):
#     """
#     计算初始地下水位深度的封装函数
#
#     参数:
#         tif_files (list): 输入TIFF文件路径列表
#         output_path (str): 输出文件路径
#         b (float/int): 土壤单元宽度，默认30
#
#     返回:
#         M_staz: 计算结果数组，同时保存到文件
#     """
#     # region 计算hw
#     with rasterio.open(os.path.join(output_path, "dem.tif")) as src:
#         global_meta = src.meta.copy()
#         global_transform = src.transform
#         bounds = src.bounds
#         minx, miny, maxx, maxy = bounds
#     # 数据读取
#     uca = rasterio.open(os.path.join(fs_path, 'uca.tif')).read(1)*b*b + b*b  # 数值要转化为汇流面积
#     slope = rasterio.open(os.path.join(fs_path, 'slope.tif')).read(1)
#     soil_depth = rasterio.open(os.path.join(fs_path, 'soil_depth.tif')).read(1)
#     ks = rasterio.open(os.path.join(fs_path, 'ks.tif')).read(1)
#     pr = rasterio.open(os.path.join(fs_path, 'pr.tif')).read(1) / 1000 / 365  # 换算为mm/year->m/year-> m/day
#     eva = rasterio.open(os.path.join(fs_path, 'eva.tif')).read(1) / 1000 / 365  # 换算为mm/year->m/year-> m/day
#     ks = np.power(10, ks) / 100  # 10^ks → cm/day → m/day
#
#     ## 定义长和宽的分块数量
#     h, w = uca.shape
#     window = from_bounds(minx, miny, maxx, maxy, transform=global_transform)
#     # 实际像素范围
#     col_off = int(window.col_off);
#     row_off = int(window.row_off)
#     width = int(window.width);
#     height = int(window.height)
#     ## 使用分块对数据进行处理
#     fenkuai_id = 2  # 每一维度分成几个部分
#     width_block = int(width / 2)
#     height_block = int(height / 2)
#     for i in trange(fenkuai_id * fenkuai_id):
#         if i == 0:
#             row_s = row_off
#             row_e = row_s + height_block
#             col_s = col_off
#             col_e = col_s + width_block
#         elif i == 1:
#             row_s = row_off
#             row_e = row_s + height_block
#             col_s = col_off + width_block
#             col_e = col_s + width_block
#         elif i == 2:
#             row_s = row_off + height_block
#             row_e = row_s + height_block
#             col_s = col_off
#             col_e = col_s + width_block
#         elif i == 3:
#             row_s = row_off + height_block
#             row_e = row_s + height_block
#             col_s = col_off + width_block
#             col_e = col_s + width_block
#
#         box_window = Window(col_s, row_s, col_e - col_s, row_e - row_s)
#         block_transform = rasterio.windows.transform(box_window, global_transform)
#
#         uca_blk = uca[row_s:row_e, col_s:col_e]
#         slope_blk = slope[row_s:row_e, col_s:col_e]
#         soil_depth_blk = soil_depth[row_s:row_e, col_s:col_e]
#         ks_blk = ks[row_s:row_e, col_s:col_e]
#         pr_blk = pr[row_s:row_e, col_s:col_e]
#         eva_blk = eva[row_s:row_e, col_s:col_e]
#
#         hw_arr = Cal_water_table(uca_blk, slope_blk, soil_depth_blk, ks_blk, pr_blk, eva_blk)
#         # del uca_blk,slope_blk,soil_depth_blk,ks_blk,pr_blk,eva_blk
#         # 准备输出
#         meta = global_meta.copy()
#         meta.update({
#             'dtype': 'float32',
#             'count': 1,
#             'height': uca_blk.shape[0],
#             'width': uca_blk.shape[1],
#             'transform': block_transform  # 关键！
#         })
#         with rasterio.open(os.path.join(fs_path, f"hw_{i + 1}.tif"), 'w',
#                            **meta) as dst:
#             dst.write(hw_arr.astype(np.float32), 1)
#     del uca, slope, soil_depth, ks, pr, eva, hw_arr, uca_blk, slope_blk, soil_depth_blk, ks_blk, pr_blk, eva_blk
#     # endregion
#     merge_blocks(fs_path, output_path)
#     print("All done!")


def calculate_initial_hw(output_path, b=30):
    """
    Calculate initial groundwater table and automatically switch to block
    processing when memory pressure is high or the full calculation is too large.
    """
    fs_path = os.path.join(output_path, "temp")
    dem_path, result_path = _resolve_paths(fs_path, output_path)
    result_dir = os.path.dirname(result_path)
    if result_dir:
        os.makedirs(result_dir, exist_ok=True)

    input_files = {
        "uca": os.path.join(fs_path, "uca.tif"),
        "slope": os.path.join(fs_path, "slope.tif"),
        "soil_depth": os.path.join(fs_path, "soil_depth.tif"),
        "ks": os.path.join(fs_path, "ks.tif"),
        "rainfall": os.path.join(fs_path, "rainfall.tif"),
        "Eva": os.path.join(fs_path, "Eva.tif"),
    }

    with rasterio.open(dem_path) as dem_src:
        global_meta = dem_src.meta.copy()
        global_transform = dem_src.transform
        bounds = dem_src.bounds
        minx, miny, maxx, maxy = bounds

    srcs = {name: rasterio.open(path) for name, path in input_files.items()}
    try:
        window = from_bounds(minx, miny, maxx, maxy, transform=global_transform)
        col_off = max(0, int(round(window.col_off)))
        row_off = max(0, int(round(window.row_off)))
        width = min(int(round(window.width)), global_meta["width"] - col_off)
        height = min(int(round(window.height)), global_meta["height"] - row_off)
        base_window = Window(col_off, row_off, width, height)
        result_transform = rasterio.windows.transform(base_window, global_transform)

        datasets = list(srcs.values())
        _, available_memory, used_ratio = _get_memory_info()
        estimated_bytes = _estimate_calculation_bytes(width, height, datasets)
        memory_is_high = used_ratio is not None and used_ratio >= MEMORY_LIMIT_RATIO
        calculation_is_large = (
            available_memory is not None
            and estimated_bytes >= available_memory * MEMORY_LIMIT_RATIO
        )
        use_blocks = memory_is_high or calculation_is_large

        if used_ratio is None:
            print("Memory check unavailable; using estimated calculation size only.")
        else:
            print(
                f"Memory used: {used_ratio:.1%}; available: {_format_gb(available_memory)}; "
                f"estimated full calculation: {_format_gb(estimated_bytes)}."
            )

        out_meta = global_meta.copy()
        out_meta.update({
            "dtype": "float32",
            "count": 1,
            "height": height,
            "width": width,
            "transform": result_transform,
        })

        if not use_blocks:
            print("Memory is sufficient; running full-array calculation.")
            uca, slope, soil_depth, ks, pr, eva = _read_inputs(srcs, window=base_window, b=b)
            hw_arr = Cal_water_table(uca, slope, soil_depth, ks, pr, eva, b=b)
            with rasterio.open(result_path, "w", **out_meta) as dst:
                dst.write(hw_arr.astype(np.float32), 1)
            print(f"All done! Result saved to {result_path}")

            # 赋值栅格重采样，输出到临时文件夹
            files_to_resample = ["hw.tif"]
            files_to_resample = ["hw.tif"]
            file_paths = [os.path.join(output_path, f) for f in files_to_resample]
            resize_clip(os.path.join(output_path, "dem.tif"), file_paths,
                        fs_path, resample_method=2)
            return None

        max_pixels = _calculate_block_pixels(available_memory, datasets)
        block_windows = list(_iter_block_windows(row_off, col_off, height, width, max_pixels))
        print(
            f"Memory ratio is over {MEMORY_LIMIT_RATIO:.0%} or estimated calculation is large; "
            f"running block calculation with {len(block_windows)} blocks."
        )

        with rasterio.open(result_path, "w", **out_meta) as dst:
            for block_index in trange(len(block_windows)):
                src_window = block_windows[block_index]
                dst_window = Window(
                    src_window.col_off - col_off,
                    src_window.row_off - row_off,
                    src_window.width,
                    src_window.height,
                )

                uca_blk, slope_blk, soil_depth_blk, ks_blk, pr_blk, eva_blk = _read_inputs(
                    srcs, window=src_window, b=b
                )
                hw_arr = Cal_water_table(
                    uca_blk, slope_blk, soil_depth_blk, ks_blk, pr_blk, eva_blk, b=b
                )
                dst.write(hw_arr.astype(np.float32), 1, window=dst_window)

                del uca_blk, slope_blk, soil_depth_blk, ks_blk, pr_blk, eva_blk, hw_arr

        print(f"All done! Result saved to {result_path}")

        # 赋值栅格重采样，输出到临时文件夹
        files_to_resample = ["hw.tif"]
        files_to_resample = ["hw.tif"]
        file_paths = [os.path.join(output_path, f) for f in files_to_resample]
        resize_clip(os.path.join(output_path, "dem.tif"), file_paths,
                    fs_path, resample_method=2)
        return None
    finally:
        for src in srcs.values():
            src.close()



if __name__ == '__main__':
    project_root = Path(__file__).resolve().parents[3]
    out_path = project_root / "example" / "output"
    calculate_initial_hw(str(out_path))
