"This module provides functions for downloading data from Google Earth Engine (GEE) using the geemap library."
import os
from pathlib import Path
import rasterio
import numpy as np
import pandas as pd
import ee
import geemap

import whitebox
wbt = whitebox.WhiteboxTools()
wbt.verbose = True
_EE_INITIALIZED = False


def initialize_gee(project=None, authenticate=True):
    """Initialize Google Earth Engine only when a download is requested."""
    global _EE_INITIALIZED

    if _EE_INITIALIZED:
        return

    if authenticate:
        ee.Authenticate()

    # 这里直接写 project id
    project = project or "landslide-495707"

    ee.Initialize(project=project)

    _EE_INITIALIZED = True

"""依据提供的shp文件将fs计算所需要的数据下载到本地

一：地形衍生因子
dem: 哥白尼数据集中的数字高程模型
slope: 依据dem计算得到的坡度数据， 目前单位为度
twi: 依据dem计算得到的湿润指数数据
uca: 依据dem计算得到的上游汇流面积数据, 记得需要乘像元面积
curve: 地形曲率，依据dem计算得到的地形曲率数据，使用的为剖面曲率

二：土壤质地因子
soil_type: 土壤质地
dry_soil: 干土密度
Cs: 土壤粘聚力
ksat: 土壤饱和水导率
suw: 饱和土壤密度
muw：湿润土壤密度


三：植被赋值与水文因子
Cr: 根系作用力 <-依据土地利用赋值得到
qt： 植被附加力<-依据AGB赋值获得
rainfall: 降水量  需要指定时间范围，一般默认为近10年的年平均值
evapotranspiration: 蒸散量  需要指定时间范围，一般默认为近10年的年平均值
river_dist: 河流距离  依据河流数据计算得到
"""

xlsx_path = r".\soil_texture.xlsx"
# soil_df = pd.read_excel(xlsx_path, sheet_name="Sheet1")
# g = 9.81  # 重力加速度 m/s²
# rho_s = 2.65  # 土粒密度 g/cm³
# region  数字高程模型下载，及其衍生因子计算
def get_zfactor(dem_path):
    with rasterio.open(dem_path) as src:
        bounds = src.bounds
        lat = (bounds.top + bounds.bottom) / 2

    zfactor = 1 / (111320 * np.cos(np.deg2rad(lat)))
    return zfactor


def download_dem(shp_path, output_dir):
    """
    地形及其衍生因子下载函数
    :param shp_path:
    :param output_dir:
    :return:
    """
    initialize_gee()
    roi = geemap.shp_to_ee(shp_path)
    dem = ee.Image("NASA/NASADEM_HGT/001").select("elevation").clip(roi)
    # 如果输出文件夹中已经存在 dem.tif 文件，则跳过下载
    if os.path.exists(os.path.join(output_dir, "dem.tif")):
        print("dem.tif already exists, skipping download.")
        return
    else:
        print("Downloading dem.tif...")
        geemap.download_ee_image(
            image=dem,
            filename=os.path.join(output_dir, "dem.tif"),
            region=roi.geometry(),
            crs="EPSG:4326",
        )



def dem_factors(dem_path, output_dir):
    """
    依据dem计算地形衍生因子，
    :param dem_path:
    :param output_dir:uca,slope,twi,curve
    :return:
    """
    # 1. 读取dem，并进行填洼处理
    wbt.breach_depressions(dem_path, os.path.join(output_dir, "dem2.tif"))
    # 2. 读取uca，并进行填洼处理
    wbt.d_inf_flow_accumulation(
        os.path.join(output_dir, "dem2.tif"),
        os.path.join(output_dir, "uca.tif"),
        out_type="Specific Contributing Area",
        log=False
    )

    # 3. 计算坡度，WetnessIndex 要求 slope 用 degrees,不要在gee中计算坡度，插值函数有问题
    wbt.slope(
        os.path.join(output_dir, "dem2.tif"),
        os.path.join(output_dir, "slope.tif"),
        units="degrees"
    )

    # 4. 计算twi
    wbt.wetness_index(
        os.path.join(output_dir, "uca.tif"),
        os.path.join(output_dir, "slope.tif"),
        os.path.join(output_dir, "twi.tif")
    )

    # 5. 计算曲率
    # zfactor = get_zfactor(os.path.join(output_dir, "dem2.tif"))
    # print(f"zfactor为{zfactor}")
    # wbt.profile_curvature(
    #     os.path.join(output_dir, "dem2.tif"),
    #     os.path.join(output_dir, "curve.tif"),
    #     zfactor = zfactor
    # )
    wbt.slope(
        os.path.join(output_dir, "dem2.tif"),
        os.path.join(output_dir, "slope1.tif"),
    )
    wbt.slope(
        os.path.join(output_dir, "slope1.tif"),
        os.path.join(output_dir, "curve.tif"),
    )
    zfactor = get_zfactor(os.path.join(output_dir, "dem2.tif"))
    wbt.profile_curvature(
        os.path.join(output_dir, "dem2.tif"),
        os.path.join(output_dir, "curve1.tif"),
        zfactor = zfactor
    )

def download_soil(shp_path, output_dir):
    """
    土壤质地因子下载函数
    :param shp_path:
    :param output_dir:
    :return:
    """
    initialize_gee()
    roi = geemap.shp_to_ee(shp_path)
    # dry_soil
    dataset = ee.Image("OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02")
    # 选择 b0 波段 表层0-5m土厚，分辨率为250m
    selectedBand = dataset.select('b0')
    geemap.download_ee_image(
        image=selectedBand,
        filename=os.path.join(output_dir, "dry_soil.tif"),
        region=roi.geometry(),
        crs="EPSG:4326",  # 常用经纬度坐标系
        scale=30,  # 分辨率 (m)
    )

    # soil_type
    dataset = ee.Image("projects/sat-io/open-datasets/FAO/HWSD_V2_SMU")
    soil_type = dataset.select("TEXTURE_USDA").toInt16()

    geemap.download_ee_image(
        image=soil_type,
        filename=os.path.join(output_dir, "soil_type.tif"),
        region=roi.geometry(),
        crs="EPSG:4326",
        scale=1000,
    )

    # Ksat下载，全球ksat已经上传到GEE，分辨率为1km
    dataset = ee.Image("projects/landslide-495707/assets/ksat")
    geemap.download_ee_image(
        image=dataset,
        filename=os.path.join(output_dir, "ks.tif"),
        region=roi.geometry(),
        crs="EPSG:4326",
    )


def download_veg_hydro(shp_path, output_dir, start_date="2012-01-01", end_date="2022-01-01"):
    """
    植被赋值与水文因子下载函数
    :param shp_path:
    :param output_dir:
    :return:
    """
    initialize_gee()
    # AGB qt
    roi = geemap.shp_to_ee(shp_path)
    dataset = ee.ImageCollection("projects/sat-io/open-datasets/ESA/ESA_CCI_AGB").filterDate(f"2021-01-01", f"2022-01-01")
    img = ee.Image(dataset.first())  # 每年 1 幅图像（AGB+SD 两个波段）
    agb = img.select('AGB')
    # 将 Mg/ha 转换为 Pa
    agb_2021_pa = agb.multiply(0.1 * 9.81)
    geemap.download_ee_image(
        image=agb_2021_pa,
        filename=os.path.join(output_dir, "qt.tif"),
        scale=100,  # 100 m 原生分辨率
        region=roi.geometry(),
        crs="EPSG:4326"
    )

    # Cr
    dataset = ee.ImageCollection("projects/sat-io/open-datasets/landcover/ESRI_Global-LULC_10m_TS")
    lulc2021 = dataset.filter(ee.Filter.stringContains('system:index', '_2021'))
    lulc2021_mosaic = lulc2021.mosaic()
    in_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    out_list = [0, 10000, 600, 2000, 0, 100, 0, 0, 5000, 2000, 10000]
    kp_img = lulc2021_mosaic.remap(in_list, out_list).rename("kp")
    geemap.download_ee_image(
        image=kp_img,
        filename=os.path.join(output_dir, "Cr.tif"),
        region=roi.geometry(),
        crs="EPSG:4326",
        scale=30,
        dtype="int16"
    )

    # rainfall
    chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")  # type: ignore
    import datetime
    current_year = datetime.datetime.now().year
    years = list(range(2009, current_year))
    annual_means = []
    # 计算每年影像的 ROI 平均降水
    for year in range(2009, 2025):
        start = ee.Date.fromYMD(year, 1, 1)
        end = ee.Date.fromYMD(year, 12, 31)
        # 年累计降水
        yearly_img = chirps.filterDate(start, end).sum().rename('precipitation').clip(roi)
        annual_means.append(yearly_img)
    # 将每年计算结果转成 FeatureCollection
    multi_year_img = ee.ImageCollection(annual_means).mean().rename('mean_precip')
    geemap.download_ee_image(
        image=multi_year_img,
        filename=os.path.join(output_dir, "rainfall.tif"),
        region=roi.geometry(),
        crs="EPSG:4326",
        scale=5000,
    )

    # Evapotranspiration
    mod16 = ee.ImageCollection("MODIS/006/MOD16A2")
    annual_images = []
    for year in range(2009, 2023):
        start = ee.Date.fromYMD(year, 1, 1)
        end = ee.Date.fromYMD(year, 12, 31)
        yearly_img = mod16.filterDate(start, end).select("ET").sum().multiply(0.1).rename('Evapotranspiration').clip(
            roi)
        annual_images.append(yearly_img)

    multi_year_img = ee.ImageCollection(annual_images).mean().rename('Evapotranspiration_mean')
    geemap.download_ee_image(
        image=multi_year_img,
        filename=os.path.join(output_dir, "Eva.tif"),
        region=roi.geometry(),
        crs="EPSG:4326",
        scale=500,
    )

    # river distance
    rivers = ee.FeatureCollection("WWF/HydroSHEDS/v1/FreeFlowingRivers").filterBounds(roi)
    dist_river = rivers.distance(searchRadius=50000, maxError=30).clip(roi).rename("dist_to_river_m")
    geemap.download_ee_image(
        image=dist_river,
        filename=os.path.join(output_dir, "river_dist.tif"),
        region=roi.geometry() if isinstance(roi, ee.featurecollection.FeatureCollection) else roi,
        scale=30,
        crs="EPSG:4326"
    )



def download_factors(shp_path, output_dir):
    # dem栅格数据下载
    download_dem(shp_path, output_dir)
    # dem相关因子计算
    dem_factors(os.path.join(output_dir, "dem.tif"), output_dir)
    # 土壤相关因子下载
    download_soil(shp_path, output_dir)
    download_veg_hydro(shp_path, output_dir)



if __name__ == '__main__':
    project_root = Path(__file__).resolve().parents[3]
    shp_path = project_root / "example" / "shp" / "demo.shp"
    output_dir = project_root / "example" / "output"
    download_factors(str(shp_path), str(output_dir))
