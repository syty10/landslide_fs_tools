"""Integrated FS calculation workflow.

The public entry point is :func:`calculate_fs_from_shp`. It keeps the existing
validated calculation functions intact and coordinates them in the required
order from a single shapefile input.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


RAW_FACTOR_FILES = (
    "dem.tif",
    "slope.tif",
    "twi.tif",
    "uca.tif",
    "curve1.tif",
    "dry_soil.tif",
    "soil_type.tif",
    "ks.tif",
    "qt.tif",
    "Cr.tif",
    "rainfall.tif",
    "Eva.tif",
    "river_dist.tif",
)

TEMP_FACTOR_FILES = (
    "slope.tif",
    "twi.tif",
    "uca.tif",
    "curve1.tif",
    "dry_soil.tif",
    "soil_type.tif",
    "ks.tif",
    "qt.tif",
    "Cr.tif",
    "rainfall.tif",
    "Eva.tif",
    "river_dist.tif",
)

TEMP_PARAMETER_FILES = ("phi.tif", "Cs.tif", "suw.tif", "muw.tif")


@dataclass(frozen=True)
class FsPipelineResult:
    """Paths produced by the integrated FS workflow."""

    output_dir: Path
    temp_dir: Path
    fs: Path
    fs_std: Path
    soil_depth: Path
    hw: Path


def _paths_exist(paths: list[Path] | tuple[Path, ...]) -> bool:
    return all(path.exists() for path in paths)


def _should_run(paths: list[Path] | tuple[Path, ...], reuse_existing: bool) -> bool:
    return not reuse_existing or not _paths_exist(paths)


def calculate_fs_from_shp(
    shp_path: str | Path,
    output_dir: str | Path,
    *,
    download: bool = True,
    reuse_existing: bool = True,
    b: int | float = 30,
    cv: float = 0.3,
    calculate_probability: bool = False,
) -> FsPipelineResult:
    """Run the complete FS workflow from one shapefile.

    Parameters
    ----------
    shp_path : str | Path
        Input shapefile used as the study area.
    output_dir : str | Path
        Folder for downloaded factors, intermediate rasters, and final outputs.
    download : bool, default True
        Download raw factors with ``def_download.py`` when they are missing.
        If all raw factors already exist and ``reuse_existing`` is True, the
        download step is skipped automatically.
    reuse_existing : bool, default True
        Reuse existing rasters when all files required by a step are present.
    b : int | float, default 30
        Soil unit width used by the initial groundwater calculation.
    cv : float, default 0.3
        Coefficient of variation used for FS uncertainty.
    calculate_probability : bool, default False
        Also calculate RI, pfail, and pf rasters after FS is complete.
    """
    shp_path = Path(shp_path).resolve()
    output_dir = Path(output_dir).resolve()
    temp_dir = output_dir / "temp"

    if not shp_path.exists():
        raise FileNotFoundError(f"Shapefile not found: {shp_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    raw_factor_paths = tuple(output_dir / name for name in RAW_FACTOR_FILES)
    if _should_run(raw_factor_paths, reuse_existing):
        if not download:
            missing = [path.name for path in raw_factor_paths if not path.exists()]
            raise FileNotFoundError(
                "Raw factor rasters are missing and download=False: "
                + ", ".join(missing)
            )
        from .tools.def_download import download_factors

        download_factors(str(shp_path), str(output_dir))
        missing_after_download = [path.name for path in raw_factor_paths if not path.exists()]
        if missing_after_download:
            raise FileNotFoundError(
                "Download finished, but these raw factor rasters are still missing: "
                + ", ".join(missing_after_download)
            )
    else:
        print("Raw factor rasters already exist; skipping download.")

    temp_factor_paths = tuple(temp_dir / name for name in TEMP_FACTOR_FILES)
    if _should_run(temp_factor_paths, reuse_existing):
        from .tools.def_assigning_param import def_resample_to_dem

        def_resample_to_dem(str(output_dir))
    else:
        print("Resampled factor rasters already exist; skipping resampling.")

    temp_parameter_paths = tuple(temp_dir / name for name in TEMP_PARAMETER_FILES)
    if _should_run(temp_parameter_paths, reuse_existing):
        from .tools.def_assigning_param import assign_params

        assign_params(str(output_dir))
    else:
        print("Assigned parameter rasters already exist; skipping assignment.")

    soil_depth_paths = (output_dir / "soil_depth.tif", temp_dir / "soil_depth.tif")
    if _should_run(soil_depth_paths, reuse_existing):
        from .core.def_calsoildepth import Cal_soilthickness

        Cal_soilthickness(str(output_dir))
    else:
        print("Soil depth raster already exists; skipping soil-depth calculation.")

    hw_paths = (output_dir / "hw.tif", temp_dir / "hw.tif")
    if _should_run(hw_paths, reuse_existing):
        from .core.def_initialhw import calculate_initial_hw

        calculate_initial_hw(str(output_dir), b=b)
    else:
        print("Initial groundwater raster already exists; skipping hw calculation.")

    fs_paths = (output_dir / "fs.tif", output_dir / "fs_std.tif")
    if _should_run(fs_paths, reuse_existing):
        from .core.def_tissa import run_fs

        run_fs(str(output_dir), cv=cv)
    else:
        print("FS rasters already exist; skipping FS calculation.")

    if calculate_probability:
        probability_paths = (
            output_dir / "RI.tif",
            output_dir / "pfail.tif",
            output_dir / "pf.tif",
        )
        if _should_run(probability_paths, reuse_existing):
            from .core.def_tissa import cal_pfail

            cal_pfail(str(output_dir))
        else:
            print("Probability rasters already exist; skipping probability calculation.")

    return FsPipelineResult(
        output_dir=output_dir,
        temp_dir=temp_dir,
        fs=output_dir / "fs.tif",
        fs_std=output_dir / "fs_std.tif",
        soil_depth=output_dir / "soil_depth.tif",
        hw=output_dir / "hw.tif",
    )
