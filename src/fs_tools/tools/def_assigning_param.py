import rasterio
from scipy.stats import kstwo
from pathlib import Path

try:
    from .def_resize import resize_clip
except ImportError:  # Allow running this file directly.
    from def_resize import resize_clip
import pandas as pd
import os
import numpy as np
"""
参数赋值：
Cs: 土壤粘聚力
ksat: 土壤饱和水导率
suw: 饱和土壤密度
muw：湿润土壤密度
phi: 土壤内摩擦角
# uca要转换为m²乘900，并且要加1，防止出现为0情况
"""
xlsx_path = Path(__file__).resolve().parent / "soil_texture.xlsx"
soil_df = pd.read_excel(xlsx_path, sheet_name="Sheet1")
g = 9.81  # 重力加速度 m/s²
rho_s = 2.65  # 土粒密度 g/cm³


def assign_params(out_path):
    """
    土壤参数赋值函数phi,Cs,suw,muw,ksat等参数的赋值函数,其中ksat启用开源1km土壤饱和导水率进行赋值
    Cr赋值函数依据土地利用类型赋值
    :param fs_path:
    :param outp_path:
    :return:
    """
    fs_path = os.path.join(out_path, "temp")
    soil_type_path = os.path.join(fs_path, "soil_type.tif")
    dry_soil_path = os.path.join(fs_path, "dry_soil.tif")

    # 读取土壤类型TIF
    with rasterio.open(soil_type_path) as src:
        soil_type = src.read(1)
        ref_profile = src.profile
        soil_nodata = src.nodata

    # 读取干土密度TIF (kg/m³)
    with rasterio.open(dry_soil_path) as src:
        dry_density = src.read(1)
        dry_nodata = src.nodata

    # assign phi
    phi = np.zeros_like(soil_type, dtype=float)
    for _, row in soil_df.iterrows():
        mask = soil_type == row["ID"]  # 假设soil_type.tif存的类别编号对应表中ID
        phi[mask] = np.deg2rad(row["Angle of Internal Friction(phi 度)"])

    phi_invalid_mask = (soil_type == soil_nodata)
    phi[phi_invalid_mask] = 0.436332  # 相当于25°
    profile = ref_profile.copy()
    profile.update(dtype=rasterio.float32, count=1)
    with rasterio.open(os.path.join(out_path, "phi.tif"), "w", **profile) as dst:
        dst.write(phi.astype(np.float32), 1)

    # Cs 依据土壤类型赋值
    Cs = np.zeros_like(soil_type, dtype=float)
    for _, row in soil_df.iterrows():
        mask = soil_type == row["ID"]  # 假设soil_type.tif存的类别编号对应表中ID
        Cs[mask] = row["SOIL COHESION(KPa)"] * 1000

    soil_invalid_mask = (soil_type == soil_nodata)
    Cs[soil_invalid_mask] = 15000  # pa 内聚力空值赋值
    profile = ref_profile.copy()
    profile.update(dtype=rasterio.float32, count=1)
    with rasterio.open(os.path.join(out_path, "Cs.tif"), "w", **profile) as dst:
        dst.write(Cs.astype(np.float32), 1)

    # muw与suw的赋值
    # 读取土壤参数表 (第二个表)
    soil_df2 = pd.read_excel(xlsx_path, sheet_name="Sheet2")
    soil_df2 = soil_df2.dropna(subset=["ID"])  # 删除 ID 为空的行
    soil_dict = soil_df2.set_index("ID")[["Field Capacity(cm3 )", "fraction"]].to_dict("index")
    # dry_density干土密度栅格值乘10方得到对应的kg/立方米    soil_type 土壤类型
    gamma_d = dry_density * 98.1  # kg/立方米->N/立方米=pa
    gamma_w = 9810
    theta_fc = np.zeros_like(soil_type, dtype=float)
    theta_sat = np.zeros_like(soil_type, dtype=float)
    # 依据ID 进行field Capacity(立方厘米)与fraction进行赋值,此处可以依据粘黏土比例进行估算，在此选用VIC模型赋值
    for soil_id, vals in soil_dict.items():
        mask = (soil_type == soil_id)
        theta_fc[mask] = vals["Field Capacity(cm3 )"]
        theta_sat[mask] = vals["fraction"]
    # 无效值替换1
    mask = (soil_type == soil_nodata)
    theta_sat[mask] = 0.45  # 饱和含水率
    theta_fc[mask] = 0.35  # 田间含水率
    # theta_avg = (theta_fc + theta_sat)/2
    suw = gamma_d + gamma_w * theta_sat
    muw = gamma_d + gamma_w * theta_fc
    mask = (dry_density == dry_nodata)
    suw[mask] = np.nan
    muw[mask] = np.nan
    for name, arr in zip(["suw", "muw"], [suw, muw]):
        profile = ref_profile.copy()
        profile.update(dtype=rasterio.float32, count=1)
        with rasterio.open(os.path.join(out_path, f"{name}.tif"), "w", **profile) as dst:
            dst.write(arr.astype(np.float32), 1)

    # ksat赋值，启用开源1km土壤饱和导水率进行赋值
    # current_file = Path(__file__).resolve()
    # src_dir = current_file.parents[2]
    # ksat_path = src_dir / "fs_datas" / "ks.tif"
    # resize_clip(os.path.join(out_path,"dem.tif"), [ksat_path], out_path, resample_method=1)

    # # Uca参数修正,在后期读取时候修正
    # uca_path = os.path.join(fs_path, "uca.tif")
    # out_uca = os.path.join(out_path, "uca1.tif")
    # with rasterio.open(uca_path) as src:
    #     global_meta = src.meta.copy()
    #     uca = src.read(1)
    # uca = (uca + 1) * 900
    # uca_meta = global_meta.copy()
    # with rasterio.open(out_uca, 'w', **uca_meta) as dst:
    #     dst.write(uca, 1)

    # 赋值栅格重采样，输出到临时文件夹
    files_to_resample3 = ["phi.tif", "Cs.tif", "suw.tif", "muw.tif"]
    file_paths = [os.path.join(out_path, f) for f in files_to_resample3]
    resize_clip(os.path.join(out_path, "dem.tif"), file_paths,
                fs_path, resample_method=2)


def def_resample_to_dem(output_dir):
    # 将下载的所有数据重采样到dem的分辨率与范围
    # 创建临时文件夹路径
    temp_path = os.path.join(output_dir, "temp")
    if not os.path.exists(temp_path):
        os.makedirs(temp_path)
    dem_path = os.path.join(output_dir, "dem.tif")
    files_to_resample1 = ["slope.tif", "twi.tif", "uca.tif", "curve1.tif", "dry_soil.tif", "soil_type.tif", "ks.tif",
                         "qt.tif", "Cr.tif", "rainfall.tif", "Eva.tif", "river_dist.tif"]
    file_paths = [os.path.join(output_dir, f) for f in files_to_resample1]
    resize_clip(dem_path, file_paths, temp_path, resample_method=2)
    files_to_resample2 = ["soil_type.tif",]
    file_paths = [os.path.join(output_dir, f) for f in files_to_resample2]
    resize_clip(dem_path, file_paths, temp_path, resample_method=0)






if __name__ == '__main__':
    project_root = Path(__file__).resolve().parents[3]
    output_dir = project_root / "example" / "output"
    # 注意此处fs_data即temp文件夹下面的uca是没有换算乘900的原始数据，计算的时候需要小心
    # 数据重采样至相同分辨率与范围，输出至临时文件夹
    def_resample_to_dem(str(output_dir))
    # 赋值后参数输出至output文件夹，调用的是临时文件夹中的重采样数据
    assign_params(str(output_dir))
