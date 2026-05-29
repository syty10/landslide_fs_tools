import rasterio
from osgeo import gdal
gdal.DontUseExceptions()
import numpy as np
import os
from pathlib import Path
try:
    from .def_resize import resize_clip
except ImportError:  # Allow running this file directly.
    from def_resize import resize_clip



def Cal_soilthickness(output_dir):
    """
    计算土壤厚度
    参数：
        river_dist  : 距河流距离
        twi         : 地形湿度指数
        curvature   : 曲率
        slope_array : 坡度（°）
    返回：
        soil_depth  : 土壤厚度（float32）
    """
    file_list = ["twi.tif", "curve1.tif", "slope.tif", "river_dist.tif"]
    # 读取数据
    fs_path = os.path.join(output_dir, "temp")
    twi = rasterio.open(os.path.join(fs_path, file_list[0])).read(1)
    curvature = rasterio.open(os.path.join(fs_path, file_list[1])).read(1)
    slope = rasterio.open(os.path.join(fs_path, file_list[2])).read(1)
    d_river = rasterio.open(os.path.join(fs_path, file_list[3])).read(1)

    # ---- 定义一个统一的清理函数 ----
    def clean_invalid(arr):
        """将 np.nan 和 np.inf 替换为该数组的最小非NaN值"""
        arr = arr.copy()
        valid_min = np.nanmin(arr[np.isfinite(arr)]) if np.any(np.isfinite(arr)) else 0
        arr[~np.isfinite(arr)] = valid_min
        return arr

    # ---- 清理输入数据 ----
    d_river = clean_invalid(d_river)
    twi = clean_invalid(twi)
    curvature = clean_invalid(curvature)
    slope = clean_invalid(slope)
    # 判断输入slope_array是否为度，如果是则转换为弧度
    if np.nanmax(slope) > 3:
        slope = np.deg2rad(slope)

    # ---- 归一化 ----
    river_max = np.nanmax(d_river) or 1
    twi_max = np.nanmax(twi) or 1
    min_slope, max_slope = np.nanmin(slope), np.nanmax(slope)
    min_curv, max_curv = np.nanmin(curvature), np.nanmax(curvature)

    normalized_slope = (slope - min_slope) / (max_slope - min_slope + 1e-9)
    normalized_curv = (curvature - min_curv) / (max_curv - min_curv + 1e-9)

    # ---- 参数 ----
    a, b, c, d, e = 0.7, 0.6, 0.7, 0.3, 0.9

    # ---- 主计算 ----
    core = 1 - a*normalized_slope - b*(d_river/river_max) + c*normalized_curv + d*(twi/twi_max)
    core = np.clip(core, 0, None)  # 负数设为0，避免出现复数
    soil_depth = core**e + 2

    # ---- 输出 ---- # 宜良区域经验公式，该区域土很厚
    soil_depth = soil_depth.astype(np.float32)+ 10
    # 输出路径
    out_soil_depth = os.path.join(output_dir, "soil_depth.tif")
    meta = rasterio.open(os.path.join(fs_path, file_list[0])).meta
    with rasterio.open(out_soil_depth, 'w', **meta) as dst:
        dst.write(soil_depth, 1)

    # 赋值栅格重采样，输出到临时文件夹
    files_to_resample3 = ["soil_depth.tif"]
    file_paths = [os.path.join(output_dir, f) for f in files_to_resample3]
    resize_clip(os.path.join(output_dir, "dem.tif"), file_paths,
                    fs_path, resample_method=2)


if __name__ == '__main__':
    # 输入路径  # 印度地区经验公式
    project_root = Path(__file__).resolve().parents[3]
    out_path = project_root / "example" / "output"
    Cal_soilthickness(str(out_path))

