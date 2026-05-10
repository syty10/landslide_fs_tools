import os

from osgeo import gdal


def resize_clip(ref_tif, file_paths, output_folder, resample_method=1):
    """
    Resample and clip rasters to match a reference raster.

    Parameters
    ----------
    ref_tif : str
        Reference raster. Its extent, projection, width, and height are used.
    file_paths : list[str]
        Rasters to resample.
    output_folder : str
        Folder where resampled rasters are written.
    resample_method : int
        0 = nearest, 1 = bilinear, 2 = cubic.
    """
    ref_ds = gdal.Open(str(ref_tif))
    if ref_ds is None:
        raise FileNotFoundError(f"Cannot open reference raster: {ref_tif}")

    ref_proj = ref_ds.GetProjection()
    ref_trans = ref_ds.GetGeoTransform()
    ref_width = ref_ds.RasterXSize
    ref_height = ref_ds.RasterYSize

    os.makedirs(output_folder, exist_ok=True)

    if resample_method == 0:
        resample_alg = gdal.GRA_NearestNeighbour
    elif resample_method == 1:
        resample_alg = gdal.GRA_Bilinear
    elif resample_method == 2:
        resample_alg = gdal.GRA_Cubic
    else:
        raise ValueError("resample_method must be 0, 1, or 2")

    output_bounds = [
        ref_trans[0],
        ref_trans[3] + ref_trans[5] * ref_height,
        ref_trans[0] + ref_trans[1] * ref_width,
        ref_trans[3],
    ]

    for path in file_paths:
        src_ds = gdal.Open(str(path))
        if src_ds is None:
            print(f"Cannot open: {path}")
            continue

        filename = os.path.basename(path)
        out_path = os.path.join(output_folder, filename)
        gdal.Warp(
            out_path,
            src_ds,
            format="GTiff",
            width=ref_width,
            height=ref_height,
            outputBounds=output_bounds,
            dstSRS=ref_proj,
            resampleAlg=resample_alg,
        )
        print(f"{filename} resampled.")
