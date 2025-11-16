import rasterio
import xarray as xr
import numpy as np
import ndpyramid
from rasterio.enums import Resampling
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
    Convert a GeoTIFF to a multi-resolution Zarr pyramid.
    
    Parameters
    ----------
    input_tif : str or Path
        Path to input GeoTIFF file
    output_zarr : str or Path
        Path to output Zarr store
    chunks : dict, optional
        Chunk sizes for each dimension. Default: {"band": 1, "y": 1024, "x": 1024}
    levels : int, optional
        Number of pyramid levels to create. Default: 4
    pyramid_method : str, optional
        Pyramid creation method: 'coarsen' (mean downsampling) or 'reproject' (supports custom resampling)
        Default: 'coarsen'
    max_memory_mb : int, optional
        Approximate max memory to use in MB for windowed reads. Default: 1000
    """
    if chunks is None:
        chunks = {"band": 1, "y": 1024, "x": 1024}
    
    try:
        # Open GeoTIFF with rasterio
        with rasterio.open(input_tif) as src:
            bands = src.count
            height = src.height
            width = src.width
            transform = src.transform
            crs = src.crs
            nodata = src.nodata
            dtype = src.dtypes[0]
            
            # Check if georeferenced
            if transform == rasterio.Affine.identity():
                warnings.warn(
                    "Dataset is not georeferenced. Using pixel coordinates."
                )
            
            # Estimate memory usage
            estimated_mb = (bands * height * width * np.dtype(dtype).itemsize) / (1024**2)
            
            # Read data (with memory consideration)
            if estimated_mb > max_memory_mb:
                warnings.warn(
                    f"Large file detected ({estimated_mb:.1f} MB). "
                    f"Consider using windowed reads or installing dask for better memory management."
                )
            
            # Read all bands
            data = src.read()
            
            # Build coordinate arrays
            if transform != rasterio.Affine.identity():
                # Georeferenced: use actual coordinates
                cols = np.arange(width)
                rows = np.arange(height)
                xs, _ = rasterio.transform.xy(transform, np.zeros(width), cols, offset='center')
                _, ys = rasterio.transform.xy(transform, rows, np.zeros(height), offset='center')
                x_coords = np.array(xs)
                y_coords = np.array(ys)
            else:
                # Not georeferenced: use pixel indices
                x_coords = np.arange(width)
                y_coords = np.arange(height)
            
            # Collect metadata
            metadata = {
                "source_file": str(Path(input_tif).name),
            }
            
            if transform != rasterio.Affine.identity():
                metadata["transform"] = transform.to_gdal()
            
            if crs:
                metadata["crs"] = crs.to_string()
            
            if nodata is not None:
                metadata["nodata"] = nodata
            
            # Add band names if available
            band_descriptions = [src.descriptions[i] or f"band_{i+1}" 
                               for i in range(bands)]
        
        # Create xarray DataArray
        da = xr.DataArray(
            data,
            dims=("band", "y", "x"),
            coords={
                "band": np.arange(1, bands + 1),
                "x": x_coords,
                "y": y_coords,
            },
            attrs=metadata
        )
        
        # Add band descriptions as coordinate
        da = da.assign_coords(band_name=("band", band_descriptions))
        
        # Wrap into dataset (don't chunk yet - will chunk during write)
        ds = xr.Dataset({"raster": da})
        
        print(f"Building {levels}-level pyramid using '{pyramid_method}' method...")
        
        # Check if dask is available
        try:
            import dask
            HAS_DASK = True
            print("Using dask for chunked operations")
        except ImportError:
            HAS_DASK = False
            print("Dask not available - using in-memory operations")
        
        # Build multi-resolution pyramid
        if HAS_DASK:
            # Chunk before pyramid building if dask is available
            ds = ds.chunk(chunks)
        
        # Choose pyramid method
        if pyramid_method == "coarsen":
            # pyramid_coarsen uses mean downsampling (no custom resampling in v0.3.x)
            pyramid = ndpyramid.pyramid_coarsen(
                ds,
                factors=[2] * levels,
                dims=["y", "x"],
                boundary="trim"
            )
        elif pyramid_method == "reproject" and crs and transform != rasterio.Affine.identity():
            # pyramid_reproject supports custom resampling but requires georeferenced data
            pyramid = ndpyramid.pyramid_reproject(
                ds,
                levels=levels,
                resampling="average"  # Options: 'average', 'bilinear', 'nearest', etc.
            )
        else:
            # Default to coarsen
            if pyramid_method == "reproject":
                warnings.warn("Reproject method requires georeferenced data. Using coarsen instead.")
            pyramid = ndpyramid.pyramid_coarsen(
                ds,
                factors=[2] * levels,
                dims=["y", "x"],
                boundary="trim"
            )
        
        # Write to Zarr with compression
        print(f"Writing to Zarr store: {output_zarr}")
        
        # DataTree structures handle encoding differently
        # Simply write without custom encoding (uses defaults)
        pyramid.to_zarr(
            output_zarr, 
            mode="w", 
            consolidated=True
        )
        
        print(f"✓ Zarr pyramid successfully created at: {output_zarr}")
        print(f"  Levels: {levels}")
        print(f"  Bands: {bands}")
        print(f"  Dimensions: {height} x {width}")
        print(f"  Chunk size: {chunks}")
        
    except rasterio.errors.RasterioIOError as e:
        print(f"✗ Error reading GeoTIFF: {e}")
        raise
    except Exception as e:
        print(f"✗ Error during conversion: {e}")
        raise


def inspect_zarr_pyramid(zarr_path):
    """
    Inspect a Zarr pyramid and print information about its structure.
    
    Parameters
    ----------
    zarr_path : str or Path
        Path to Zarr store
    """
    import zarr
    
    store = zarr.open(zarr_path, mode='r')
    print(f"\n=== Zarr Pyramid Structure: {zarr_path} ===")
    print(store.tree())
    
    # Load with xarray
    ds = xr.open_zarr(zarr_path, consolidated=True)
    print(f"\n=== Dataset Info ===")
    print(ds)



def addNGIF(output_zarr):
    import zarr
    root = zarr.open(output_zarr, mode="r+")

    # Count pyramid levels
    level_keys = sorted([k for k in root.group_keys() if k.isdigit()], key=int)

    multiscales = [{
        "version": "0.4",
        "name": "raster",
        "type": "pyramid",
        "axes": [
            {"name": "band", "type": "other"},
            {"name": "y", "type": "space"},
            {"name": "x", "type": "space"},
        ],
        "datasets": [{"path": k} for k in level_keys],
        "pixels_per_tile": 256,      # REQUIRED BY CARBONPLAN
    }]

    root.attrs["multiscales"] = multiscales

    # Consolidate metadata
    zarr.consolidate_metadata(output_zarr)

    print("✓ Added NGFF / multiscales metadata")
    print(f"✓ Pyramid levels detected: {level_keys}")


if __name__ == "__main__":
    # Example usage
    OUTPUT_ZARR = "output.zarr"
    INPUT_TIF = "sample.tif"
    geotiff_to_zarr(
        input_tif=INPUT_TIF,
        output_zarr=OUTPUT_ZARR,
        chunks={"band": 1, "y": 512, "x": 512},
        levels=4,
        pyramid_method="coarsen"  # Options: "coarsen" (mean) or "reproject" (custom resampling)
    )
    addNGIF(OUTPUT_ZARR)
    
    # Optional: inspect the output
    # inspect_zarr_pyramid("output.zarr")
