import rasterio
import xarray as xr
import numpy as np
import ndpyramid
from pathlib import Path
import warnings


def geotiff_to_zarr(
    input_tif, 
    output_zarr, 
    chunks=None,
    levels=4,
    pyramid_method="coarsen",
    max_memory_mb=1000
):
    """
    Convert a multi-band GeoTIFF to a multi-resolution Zarr pyramid (NGFF-compatible).
    Fully supports any number of bands and preserves band names.
    """
    if chunks is None:
        chunks = {"band": 1, "y": 1024, "x": 1024}

    try:
        # ------------------------------
        # 1. Open GeoTIFF + Read Metadata
        # ------------------------------
        with rasterio.open(input_tif) as src:
            bands = src.count
            height = src.height
            width = src.width
            transform = src.transform
            crs = src.crs
            nodata = src.nodata
            dtype = src.dtypes[0]

            # Extract band names
            band_names = [
                src.descriptions[i] if src.descriptions[i] else f"band_{i+1}"
                for i in range(bands)
            ]

            # Estimate memory usage
            estimated_mb = (bands * height * width * np.dtype(dtype).itemsize) / 1024**2
            if estimated_mb > max_memory_mb:
                warnings.warn(
                    f"Large file (~{estimated_mb:.1f} MB). "
                    "Using Dask is recommended for optimal performance."
                )

            # Read full raster (band, y, x)
            data = src.read()

            # Build spatial coordinates
            if transform != rasterio.Affine.identity():
                cols = np.arange(width)
                rows = np.arange(height)
                xs, _ = rasterio.transform.xy(transform, np.zeros(width), cols)
                _, ys = rasterio.transform.xy(transform, rows, np.zeros(height))
                x_coords = np.array(xs)
                y_coords = np.array(ys)
            else:
                x_coords = np.arange(width)
                y_coords = np.arange(height)

            # Metadata
            metadata = {
                "source_file": str(Path(input_tif).name),
                "crs": crs.to_string() if crs else None,
                "transform": transform.to_gdal(),
                "nodata": nodata,
            }

        # ---------------------------------------
        # 2. Create Xarray Dataset with band names
        # ---------------------------------------
        da = xr.DataArray(
            data,
            dims=("band", "y", "x"),
            coords={
                "band": np.arange(1, bands + 1),
                "band_name": ("band", band_names),
                "x": x_coords,
                "y": y_coords,
            },
            attrs=metadata,
        )

        ds = xr.Dataset({"raster": da})

        # --------------------------
        # 3. Chunk if Dask available
        # --------------------------
        try:
            import dask
            ds = ds.chunk(chunks)
            print("✓ Dask available — chunking enabled.")
        except ImportError:
            print("⚠ Dask NOT installed — processing in memory.")

        print(f"\nBuilding multiscale pyramid ({levels} levels)...")

        # -------------------------
        # 4. Build multiscale pyramid
        # -------------------------
        if pyramid_method == "coarsen":
            pyramid = ndpyramid.pyramid_coarsen(
                ds, 
                factors=[2] * levels, 
                dims=["y", "x"], 
                boundary="trim"
            )
        else:
            pyramid = ndpyramid.pyramid_reproject(
                ds,
                levels=levels,
                resampling="average",
            )

        # ---------------------
        # 5. Write to Zarr store
        # ---------------------
        print(f"\nWriting to Zarr store: {output_zarr}")
        pyramid.to_zarr(output_zarr, mode="w", consolidated=True)

        print("\n✓ Zarr pyramid created successfully.")
        print(f"  Bands: {bands}")
        print(f"  Band names: {band_names}")
        print(f"  Dimensions: {height} × {width}")
        print(f"  Chunks: {chunks}")

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        raise



def addNGIF(output_zarr):
    """
    Add NGFF multiscales metadata for CarbonPlan Maps / Viv / Vizarr compatibility.
    """
    import zarr

    root = zarr.open(output_zarr, mode="r+")

    # Detect pyramid levels (0, 1, 2, 3, ...)
    level_keys = sorted(
        [k for k in root.group_keys() if k.isdigit()],
        key=int
    )

    multiscales = [{
        "version": "0.4",
        "name": "raster",
        "type": "pyramid",
        "axes": [
            {"name": "band", "type": "channel"},
            {"name": "y",    "type": "space"},
            {"name": "x",    "type": "space"},
        ],
        "datasets": [{"path": k} for k in level_keys],
        "pixels_per_tile": 256,
    }]

    root.attrs["multiscales"] = multiscales
    zarr.consolidate_metadata(output_zarr)

    print("\n✓ NGFF metadata added.")
    print(f"✓ Levels detected: {level_keys}\n")



# -------------------------
# 6. RUN DIRECTLY FROM HERE
# -------------------------
if __name__ == "__main__":

    # >>>>>> EDIT THESE TWO PATHS <<<<<<

    INPUT_TIF  = r"C:/Users/Batman/Desktop/Tech/BTP/Data/flood_probability_2023-09-25_88.19823706960149_24.314907220056885.tif"
    OUTPUT_ZARR = r"C:/Users/Batman/Desktop/Tech/BTP/Data/flood_pyramid_multiband.zarr"

    # Convert
    geotiff_to_zarr(
        input_tif=INPUT_TIF,
        output_zarr=OUTPUT_ZARR,
        chunks={"band": 1, "y": 1024, "x": 1024},
        levels=6,
        pyramid_method="coarsen"
    )

    # Add NGFF metadata
    addNGIF(OUTPUT_ZARR)

    print("✓ Completed.")
