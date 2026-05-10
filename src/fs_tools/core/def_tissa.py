import os
import glob
import numpy as np
from osgeo import gdal, gdalconst
import scipy.special as ss
import scipy as sc
import rasterio
from rasterio.windows import from_bounds,Window
from tqdm import trange
from rasterio.merge import merge
import ctypes
import math
from pathlib import Path

MEMORY_LIMIT_RATIO = 0.8
WORKING_ARRAY_COUNT = 60


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


def _format_gb(byte_count):
    if byte_count is None:
        return "unknown"
    return f"{byte_count / 1024 ** 3:.2f} GB"


def _estimate_calculation_bytes(width, height, datasets):
    dtype_size = 8
    for dataset in datasets:
        dtype_size = max(dtype_size, np.dtype(dataset.dtypes[0]).itemsize)
    return int(width * height * dtype_size * WORKING_ARRAY_COUNT)


def _calculate_block_pixels(available_memory, datasets):
    dtype_size = 8
    for dataset in datasets:
        dtype_size = max(dtype_size, np.dtype(dataset.dtypes[0]).itemsize)

    bytes_per_pixel = dtype_size * WORKING_ARRAY_COUNT
    if available_memory is None:
        target_bytes = 256 * 1024 ** 2
    else:
        target_bytes = min(max(available_memory * 0.2, 32 * 1024 ** 2), 512 * 1024 ** 2)
    return max(1, int(target_bytes / bytes_per_pixel))


def _iter_block_windows(row_off, col_off, height, width, max_pixels):
    block_size = max(1, int(math.sqrt(max_pixels)))
    row_end = row_off + height
    col_end = col_off + width

    for row_s in range(row_off, row_end, block_size):
        for col_s in range(col_off, col_end, block_size):
            row_e = min(row_s + block_size, row_end)
            col_e = min(col_s + block_size, col_end)
            yield Window(col_s, row_s, col_e - col_s, row_e - row_s)


def _resolve_output_paths(out_path, uncertainty_out_path=None):
    if out_path.lower().endswith((".tif", ".tiff")):
        fs_out_path = out_path
        output_dir = os.path.dirname(out_path)
    else:
        output_dir = out_path
        fs_out_path = os.path.join(output_dir, "fs.tif")

    if uncertainty_out_path is None:
        root, ext = os.path.splitext(fs_out_path)
        uncertainty_out_path = f"{root}_std{ext or '.tif'}"

    return fs_out_path, uncertainty_out_path


def _read_fs_inputs(srcs, window=None):
    Cr = srcs["Cr"].read(1, window=window)
    Cs = srcs["Cs"].read(1, window=window)
    muw = srcs["muw"].read(1, window=window)
    soil_depth = srcs["soil_depth"].read(1, window=window)
    suw = srcs["suw"].read(1, window=window)
    hw = srcs["hw"].read(1, window=window)
    slope = srcs["slope"].read(1, window=window)
    phi = srcs["phi"].read(1, window=window)
    qt = srcs["qt"].read(1, window=window)

    soil_depth[soil_depth < 1] = 1
    nodata = srcs["qt"].nodata
    if nodata is not None:
        qt[qt == nodata] = 0
    qt = np.nan_to_num(qt, nan=0)
    if np.nanmax(slope) > 3:
        slope = np.deg2rad(slope)

    return Cr, Cs, muw, soil_depth, suw, hw, slope, phi, qt

def merge_blocks(temp_path, output_path,name_fs):
    """
    将分块的tif拼接为一个完整tif

    参数:
        temp_path (str): 存放分块结果的目录
        output_path (str): 输出文件路径
    """
    # 找到所有分块tif
    search_path = os.path.join(temp_path, name_fs)
    tif_files = sorted(glob.glob(search_path))

    if len(tif_files) == 0:
        raise FileNotFoundError(f"在 {temp_path} 中没有找到 fs_*.tif 文件")

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

    # 关闭文件
    for src in src_files_to_mosaic:
        src.close()

    print(f"拼接完成，结果保存到 {output_path}")

def Calculate_fs(Cr,Cs,muw,soil_depth,suw,hw,slope,phi,qt,rw=9810):
    """
    数据易发性分析，各参数已经确定完毕
    fs_arr = Calculate_fs(Cr_blk, Cs_blk, muw_blk,soil_depth_blk, suw_blk, hw_blk, slope_blk,phi_blk,qt_blk,rw=rw)
    """
    np.seterr(divide="ignore", invalid="ignore")
    tan_phi = np.tan(phi)
    cos_slope = np.cos(slope)
    sin_slope = np.sin(slope)
    numerator = Cs + Cr + (qt+muw * soil_depth + (suw - rw - muw) *hw * soil_depth) * cos_slope * cos_slope* tan_phi
    denominator = (qt+muw * soil_depth + (suw - muw) * hw * soil_depth) * cos_slope * sin_slope
    fs = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, 10, dtype=np.float64),
        where=denominator != 0
    )
    fs[np.isinf(fs)] = 10
    fs[np.isnan(fs)] = 10
    fs[slope< 0.087] = 10
    fs = np.clip(fs, 0, 10)

    return fs


def Calculate_fs_uncertainty(Cr,Cs,gm,D,gs,Hw,beta,phi,qt,gw=9810,cv=0.3):
    """
    此程序确定fs的置信边界
    """
    # FS 对各变量的偏导数

    np.seterr(divide="ignore", invalid="ignore")
    common = (qt+gm*D+(gs-gm)*Hw*D)*np.sin(beta)*np.cos(beta)
    eps = 0.001
    numerator = (Cr + Cs + (qt + gm * D + (gs - gw - gm) * Hw * D) * (np.cos(beta)) ** 2 * np.tan(phi))
    common = np.where(np.abs(common) < eps, eps, common)  # 避免除0
    Dcr = 1 / common
    Dcs = 1 / common
    Dqt = ((np.cos(beta)) ** 2 * np.tan(phi) * common
           - np.sin(beta) * np.cos(beta) * numerator) / (common ** 2)
    Dgm = (((np.cos(beta)) ** 2 * np.tan(phi) * (D - Hw * D)) * common
           - np.sin(beta) * np.cos(beta) * (D - Hw * D) * numerator) / (common ** 2)
    Dgs = (((np.cos(beta)) ** 2 * np.tan(phi) * (Hw * D)) * common
           - np.sin(beta) * np.cos(beta) * (Hw * D) * numerator) / (common ** 2)
    DD = (((np.cos(beta)) ** 2 * np.tan(phi) * (gm + (gs - gw - gm) * Hw)) * common
          - np.sin(beta) * np.cos(beta) * (gm + (gs - gm) * Hw) * numerator) / (common ** 2)
    DHw = (((gs - gw - gm) * D * (np.cos(beta)) ** 2 * np.tan(phi)) * common
           - (gs - gm) * D * np.sin(beta) * np.cos(beta) * numerator) / (common ** 2)
    Dbeta = ((2 * np.cos(beta) * (-np.sin(beta)) * (qt + gm * D + (gs - gw - gm) * Hw * D) * np.tan(phi)) * common
             - (qt + gm * D + (gs - gm) * Hw * D) * (
                         (np.cos(beta)) ** 2 + (-np.sin(beta)) * np.sin(beta)) * numerator) / (common ** 2)
    Dphi = ((qt + gm * D + (gs - gw - gm) * Hw * D) * (np.cos(beta)) ** 2 * (1 / (np.cos(phi) ** 2))) / common

    # ====== 2. 假设输入参数的标准差 (CV = 0.3 默认) ======
    Cr_SD = cv * Cr
    Cs_SD = cv * Cs
    qt_SD = cv * qt
    gm_SD = cv * gm
    gs_SD = cv * gs
    Hw_SD = cv * Hw
    D_SD = cv * D
    beta_SD = cv * beta
    phi_SD = cv * phi

    # ====== 3. 误差传播公式计算FS_SD ======
    FS_SD = np.sqrt(
        Dcr ** 2 * Cr_SD ** 2 +
        Dcs ** 2 * Cs_SD ** 2 +
        Dqt ** 2 * qt_SD ** 2 +
        Dgm ** 2 * gm_SD ** 2 +
        Dgs ** 2 * gs_SD ** 2 +
        DD ** 2 * D_SD ** 2 +
        DHw ** 2 * Hw_SD ** 2 +
        Dbeta ** 2 * beta_SD ** 2 +
        Dphi ** 2 * phi_SD ** 2
    )
    mask = FS_SD>10
    FS_SD[mask] = 0
    return FS_SD


def run_fs(out_path, uncertainty_out_path=None, cv=0.3):
    """
    Calculate FS and FS uncertainty in one run.

    The function reads the memory status automatically. If current memory usage is
    over 80%, or the estimated full calculation would use more than 80% of
    available memory, it switches to block-window processing.
    """
    fs_path = os.path.join(out_path, "temp")
    fs_out_path, fs_std_out_path = _resolve_output_paths(out_path, uncertainty_out_path)
    fs_out_dir = os.path.dirname(fs_out_path)
    fs_std_out_dir = os.path.dirname(fs_std_out_path)
    if fs_out_dir:
        os.makedirs(fs_out_dir, exist_ok=True)
    if fs_std_out_dir:
        os.makedirs(fs_std_out_dir, exist_ok=True)

    input_files = {
        "Cr": os.path.join(fs_path, "Cr.tif"),
        "Cs": os.path.join(fs_path, "Cs.tif"),
        "muw": os.path.join(fs_path, "muw.tif"),
        "soil_depth": os.path.join(fs_path, "soil_depth.tif"),
        "suw": os.path.join(fs_path, "suw.tif"),
        "hw": os.path.join(fs_path, "hw.tif"),
        "slope": os.path.join(fs_path, "slope.tif"),
        "phi": os.path.join(fs_path, "phi.tif"),
        "qt": os.path.join(fs_path, "qt.tif"),
    }

    with rasterio.open(os.path.join(out_path, "dem.tif")) as dem_src:
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
            print("Memory check unavailable; running block calculation for safety.")
            use_blocks = True
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
            print("Memory is sufficient; running full-array FS and uncertainty calculation.")
            Cr, Cs, muw, soil_depth, suw, hw, slope, phi, qt = _read_fs_inputs(
                srcs, window=base_window
            )
            fs_arr = Calculate_fs(
                Cr=Cr, Cs=Cs, muw=muw, soil_depth=soil_depth,
                suw=suw, hw=hw, slope=slope, phi=phi, qt=qt, rw=9810
            )
            fs_std = Calculate_fs_uncertainty(
                Cr=Cr, Cs=Cs, gm=muw, D=soil_depth,
                gs=suw, Hw=hw, beta=slope, phi=phi, qt=qt, gw=9810, cv=cv
            )
            fs_std = np.nan_to_num(fs_std, nan=0, posinf=0, neginf=0)

            with rasterio.open(fs_out_path, "w", **out_meta) as dst:
                dst.write(fs_arr.astype(np.float32), 1)
            with rasterio.open(fs_std_out_path, "w", **out_meta) as dst:
                dst.write(fs_std.astype(np.float32), 1)

            print(f"FS saved to {fs_out_path}")
            print(f"FS uncertainty saved to {fs_std_out_path}")
            print("All done!")
            return fs_arr, fs_std

        max_pixels = _calculate_block_pixels(available_memory, datasets)
        block_windows = list(_iter_block_windows(row_off, col_off, height, width, max_pixels))
        print(
            f"Memory ratio is over {MEMORY_LIMIT_RATIO:.0%} or estimated calculation is large; "
            f"running block calculation with {len(block_windows)} blocks."
        )

        with rasterio.open(fs_out_path, "w", **out_meta) as fs_dst, \
                rasterio.open(fs_std_out_path, "w", **out_meta) as fs_std_dst:
            for block_index in trange(len(block_windows)):
                src_window = block_windows[block_index]
                dst_window = Window(
                    src_window.col_off - col_off,
                    src_window.row_off - row_off,
                    src_window.width,
                    src_window.height,
                )

                Cr, Cs, muw, soil_depth, suw, hw, slope, phi, qt = _read_fs_inputs(
                    srcs, window=src_window
                )
                fs_arr = Calculate_fs(
                    Cr=Cr, Cs=Cs, muw=muw, soil_depth=soil_depth,
                    suw=suw, hw=hw, slope=slope, phi=phi, qt=qt, rw=9810
                )
                fs_std = Calculate_fs_uncertainty(
                    Cr=Cr, Cs=Cs, gm=muw, D=soil_depth,
                    gs=suw, Hw=hw, beta=slope, phi=phi, qt=qt, gw=9810, cv=cv
                )
                fs_std = np.nan_to_num(fs_std, nan=0, posinf=0, neginf=0)

                fs_dst.write(fs_arr.astype(np.float32), 1, window=dst_window)
                fs_std_dst.write(fs_std.astype(np.float32), 1, window=dst_window)

                del Cr, Cs, muw, soil_depth, suw, hw, slope, phi, qt, fs_arr, fs_std

        print(f"FS saved to {fs_out_path}")
        print(f"FS uncertainty saved to {fs_std_out_path}")
        print("All done!")
        return None
    finally:
        for src in srcs.values():
            src.close()

def cal_pfail(out_path):
    """
    计算pfail与RI
    """
    fs_mean = rasterio.open(os.path.join(out_path, "fs.tif")).read(1)
    fs_std = rasterio.open(os.path.join(out_path, "fs_std.tif")).read(1)
    with rasterio.open(os.path.join(out_path, "dem.tif")) as dem_src:
        global_meta = dem_src.meta.copy()
    out_meta = global_meta.copy()
    out_meta.update({
        "dtype": "float32",
        "count": 1,
    })

    RI = (fs_mean - 1) / fs_std
    output_path_RI = os.path.join(out_path, 'RI.tif')
    with rasterio.open(output_path_RI, "w", **out_meta) as dst:
        dst.write(RI.astype(np.float32), 1)

    x = -1 * (np.log(fs_mean) - 0.5 * np.log((fs_std / fs_mean) ** 2 + 1)) / (np.sqrt(2) * np.sqrt(np.log((fs_std / fs_mean) ** 2 + 1)))
    # 利用对数正态分布的累计分布函数计算Pfail
    # 估计对数正态分布的参数
    mu = np.log(fs_mean) - 0.5 * np.log((fs_std / fs_mean) ** 2 + 1)
    sigma = np.sqrt(np.log((fs_std / fs_mean) ** 2 + 1))
    # 计算Pfail
    Pfail = 0.5 + 0.5 * ss.erf((np.log(1.25) - mu) / (np.sqrt(2) * sigma))
    # 输出Pfail
    output_path_Pfail = os.path.join(out_path, 'pfail.tif')
    with rasterio.open(output_path_Pfail, "w", **out_meta) as dst:
        dst.write(Pfail.astype(np.float32), 1)

    try:
        import scipy as sc
        pf = sc.stats.norm.cdf(x)
        output_path_pf = os.path.join(out_path, 'pf.tif')
        with rasterio.open(output_path_pf, "w", **out_meta) as dst:
            dst.write(pf.astype(np.float32), 1)

    finally:
        print("pf tif need scipy to calculate, if you want to calculate pf, please install scipy and run this function again.")



if __name__ == '__main__':
    project_root = Path(__file__).resolve().parents[3]
    out_path = project_root / "example" / "output"
    run_fs(str(out_path))
