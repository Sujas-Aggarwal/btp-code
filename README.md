# Scalable Cloud-Optimised Geospatial Dashboard

**[B.Tech Project](https://github.com/Sujas-Aggarwal/btp-code)**  
Department of Civil Engineering, IIT Delhi  
November 2025  
By Sujas Kumar

***

## Overview

This repository contains code and resources for the B.Tech project on developing a scalable, cloud-optimised geospatial dashboard. The focus is on efficient handling of large flood inundation raster datasets, conversion to cloud-ready formats, and foundational steps toward real-time geospatial visualization.

### Core Components

- **DummyZarrVisualizer.html**  
  Simple HTML file for visualizing Zarr-format raster tiles (demo for Rajasthan region).

- **TiffToZarr.py**  
  Python script for converting GeoTIFF files to multi-band, multi-resolution Zarr arrays. Supports extraction of spatial metadata, chunking, and pyramid generation for scalable storage and access.

- **server.py**  
  CORS-enabled HTTP server for serving local raster/Zarr assets; helps in benchmarking and debugging visualization and data access.

***

## Repository Structure

- [`DummyZarrVisualizer.html`](https://github.com/Sujas-Aggarwal/btp-code/blob/master/DummyZarrVisualizer.html) — Simple visualizer for Zarr outputs
- [`TiffToZarr.py`](https://github.com/Sujas-Aggarwal/btp-code/blob/master/TiffToZarr.py) — GeoTIFF to Zarr conversion pipeline (multiband/chunked)
- [`server.py`](https://github.com/Sujas-Aggarwal/btp-code/blob/master/server.py) — Local dev server for visualizations

***

## Getting Started

1. **Clone this repository:**
   ```bash
   git clone https://github.com/Sujas-Aggarwal/btp-code.git
   cd btp-code
   ```
2. **Install the required Python packages:**
   ```bash
   pip install rasterio xarray numpy ndpyramid zarr flask
   ```
3. **Convert GeoTIFF to Zarr:**
   Execute `TiffToZarr.py` with your input file as described in the script's docstring.
4. **Visualize Zarr output:**
   Open `DummyZarrVisualizer.html` in your browser and follow comments for rudimentary visualization.
5. **Run local server (optional):**
   ```bash
   python server.py
   ```
   This starts a simple HTTP server for your Zarr or raster data assets.

***

## Data & Methods

- **Data sources:** Flood inundation GeoTIFF datasets derived from Sentinel-1A SAR data (see project documentation for data provenance).
- **Processing:** Conversion to Zarr format includes chunking, metadata extraction (CRS, transform, band labels), and multiscale pyramid generation.
- **Visualization:** The basic HTML visualizer provides a starting point for viewing Zarr raster tiles; further dashboard development may follow.

***

## License, Contact & Acknowledgement

This code is for academic use as part of B.Tech project requirements.  
Supervisor: Prof. Manabendra Saharia (IIT Delhi)  
Contact: [ce1221047@iitd.ac.in](mailto:ce1221047@iitd.ac.in)  
GitHub: [@Sujas-Aggarwal](https://github.com/Sujas-Aggarwal)

***

© 2025 Sujas Kumar | IIT Delhi B.Tech Project
