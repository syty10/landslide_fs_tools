# landslide_fs_tools

`landslide_fs_tools` is a Python toolkit for calculating landslide factor of
safety (FS) from a study-area shapefile. It downloads and prepares the required
terrain, soil, vegetation, and hydrological factors, then runs the FS workflow
through a single high-level function.

## Workflow

The integrated pipeline runs these steps:

1. Download input factors from Google Earth Engine.
2. Resample all factors to the DEM grid.
3. Assign soil and vegetation parameters.
4. Calculate soil depth.
5. Calculate initial groundwater table height.
6. Calculate FS and FS uncertainty rasters.

The final output is `fs.tif`; important intermediate rasters such as
`soil_depth.tif`, `hw.tif`, and `fs_std.tif` are also saved.

## Installation

Create and activate a Python environment, then install dependencies:

```bash
pip install -r requirements.txt
```

GDAL and rasterio can be sensitive to platform-specific binary dependencies. If
`pip` installation fails, installing them with conda first is usually more
reliable:

```bash
conda install -c conda-forge gdal rasterio
pip install -r requirements.txt
```

When running directly from the repository, make sure `src` is on your Python
path. The demo notebook does this automatically.

## Google Earth Engine Setup

The download stage uses Google Earth Engine through `earthengine-api` and
`geemap`. Authenticate before running the full workflow:

```bash
earthengine authenticate
```

If your Earth Engine account requires an explicit project ID, set `EE_PROJECT`:

```bash
export EE_PROJECT="your-earth-engine-project-id"
```

On Windows PowerShell:

```powershell
$env:EE_PROJECT = "your-earth-engine-project-id"
```

Do not commit local credential files such as `client_secret.json`.

## Quick Start

```python
from pathlib import Path
import sys

project_root = Path.cwd()
sys.path.insert(0, str(project_root / "src"))

from fs_tools import calculate_fs_from_shp

result = calculate_fs_from_shp(
    shp_path=project_root / "example" / "shp" / "demo.shp",
    output_dir=project_root / "example" / "output_demo",
    download=True,
    reuse_existing=True,
)

print(result.fs)
```

To reuse already downloaded factors and run only the calculation stages, set:

```python
calculate_fs_from_shp(
    "example/shp/demo.shp",
    "example/output_demo",
    download=False,
    reuse_existing=True,
)
```

## Demo Notebook

Open the notebook below for a complete example with maps of the final FS raster
and intermediate outputs:

```text
example/demo_fs_workflow.ipynb
```

The notebook uses `example/shp/demo.shp` and writes generated rasters to
`example/output_demo/`.

## Main API

```python
calculate_fs_from_shp(
    shp_path,
    output_dir,
    download=True,
    reuse_existing=True,
    b=30,
    cv=0.3,
    calculate_probability=False,
)
```

Parameters:

- `shp_path`: input study-area shapefile.
- `output_dir`: folder for downloaded factors, intermediate rasters, and final
  outputs.
- `download`: download missing raw factors from Google Earth Engine.
- `reuse_existing`: skip a step when all required files for that step already
  exist.
- `b`: soil unit width used for the initial groundwater calculation.
- `cv`: coefficient of variation used for FS uncertainty.
- `calculate_probability`: additionally calculate `RI.tif`, `pfail.tif`, and
  `pf.tif`.

## Output Files

Typical output files include:

- `fs.tif`: landslide factor of safety.
- `fs_std.tif`: propagated FS uncertainty.
- `soil_depth.tif`: calculated soil depth.
- `hw.tif`: initial groundwater table height.
- `temp/`: resampled intermediate rasters used by the calculation.

## Repository Layout

```text
src/fs_tools/
  pipeline.py              High-level integrated workflow
  core/                    Core raster calculations
  tools/                   Download, resampling, and parameter assignment tools
example/
  shp/                     Example study-area shapefiles
  demo_fs_workflow.ipynb   End-to-end notebook demo
```

## License

This project is released under the MIT License. See `LICENSE` for details.
