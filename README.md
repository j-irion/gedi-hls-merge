# GEDI–HLS Joiner (30 m)

This script builds a **per-shot, per-month** dataset by aligning **GEDI L4A** above-ground biomass density (AGBD) measurements with **HLS v2.0** optical reflectance from **both** Sentinel-2 (HLSS30) and Landsat-8/9 (HLSL30).

For each month you choose, it:
1. Downloads GEDI L4A granules, keeps only **high-quality shots** inside your AOI.
2. Finds HLS (L30+S30) scenes over a **rolling window** around that month.
3. Masks clouds/shadows (HLS QA/Fmask).
4. Samples HLS bands **at the GEDI shot locations** and takes the **nan-median** across all clear looks in the window.
5. Writes a CSV with one row per GEDI shot per month and these columns:

`footprint_id, lon, lat, agbd, agbd_se, beam, shot_number, ym, Blue, Green, Red, NIR, SWIR1, SWIR2`

> Notes  
> • `agbd` is in **Mg/ha** (from GEDI L4A).  
> • HLS surface reflectance is typically **scaled** (e.g., scale factor 0.0001). This script preserves the stored integer values; apply scaling downstream if needed.

---

## Requirements
- A working GDAL/rasterio environment

Install Python deps (example):
```bash
pip install earthaccess rasterio pystac-client planetary-computer pandas numpy h5py netCDF4 shapely tqdm requests python-dateutil
```

### Authenticate to NASA Earthdata
Create `~/.netrc` with your Earthdata Login (https://urs.earthdata.nasa.gov):
```
machine urs.earthdata.nasa.gov
  login YOUR_USERNAME
  password YOUR_PASSWORD
```
Then the script will log in via `earthaccess.login(strategy="netrc")`.

---

## Quick start

Run with defaults (same AOI/month as in the code, July 2022, ±7-day window, up to 200 HLS items/month):
```bash
python gedi_hls_join.py
```

Choose a different area (bbox), months, and output location:
```bash
python gedi_hls_join.py   --aoi-bbox -123.8,45.0,-122.9,45.8   --year 2021 --months 6,7,8   --workdir ./work_pnw   --out-csv biomass_join_2021_summer.csv
```

Use all months and allow clear water pixels:
```bash
python gedi_hls_join.py   --months all   --accept-water
```

Pass a GeoJSON AOI (first feature used) or a WKT polygon:
```bash
python gedi_hls_join.py --aoi-geojson ./aoi.geojson
python gedi_hls_join.py --aoi-wkt "POLYGON((-150 64,-150 65,-148 65,-148 64,-150 64))"
```

---

## CLI options

**AOI (choose one; priority is bbox > geojson > wkt):**
- `--aoi-bbox minLon,minLat,maxLon,maxLat`  
  Example: `--aoi-bbox -97.5,46.5,-90.0,49.5`
- `--aoi-geojson /path/to/aoi.geojson`  
  Accepts Feature, FeatureCollection, or bare Geometry.
- `--aoi-wkt "WKT_STRING"`

**Time selection**
- `--year INT` (default: `2022`)
- `--months STR` months list like `"6,7,8"` or `"all"` (default: `"7"`)

**Rolling window & limits**
- `--rolling-days INT` ± days around each month for HLS median (default: `7`)
- `--max-hls-per-month INT` cap on HLS scenes (L30+S30 combined) per month (default: `200`)
- `--limit-gedi-per-month INT` cap on GEDI granules per month (default: `40`; use `--limit-gedi-per-month 0` to read all)

**Workspace & outputs**
- `--workdir PATH` root folder for caches & outputs (default: `"."`)  
  Creates:
  - `workdir/gedi_cache/`
  - `workdir/hls_cache_hls_v2/`
- `--out-csv PATH` output CSV; if relative, it’s placed under `--workdir` (default: `gedi_hls_monthly_30m.csv`)

**Product & quality**
- `--gedi-collection-id STR` GEDI L4A concept ID (default: `C2237824918-ORNL_CLOUD`)
- `--accept-water` keep pixels flagged as water (off by default)
- `--accept-adjacent` keep pixels adjacent to cloud/shadow (off)
- `--accept-snow` keep snow/ice pixels (off)

---

## What the output means

- **Alignment:** Each row is a GEDI footprint aligned to the **nearest HLS pixel** (via sampling in the tile CRS).
- **Missing values:** If no clear HLS pixel is available for a shot in the window, the band values are `NaN` for that month.
- **Bands:**  
  - Blue = B02  
  - Green = B03  
  - Red = B04  
  - NIR = **B08** (S2) or **B05** (Landsat)  
  - SWIR1 = **B11** (S2) or **B06** (Landsat)  
  - SWIR2 = **B12** (S2) or **B07** (Landsat)

---

## Setting the working directory

Use `--workdir` to control **where caches and outputs go**. Example:
```bash
python gedi_hls_join.py --workdir /data/biomass_runs/alaska_2020 --out-csv alaska_2020.csv
```
You can safely delete `gedi_cache/` and `hls_cache_hls_v2/` afterward to free space.

---

## Performance & disk tips

- HLS tiles can be large; the cache for one month over a big AOI can be **tens of GB**.  
  Reduce with `--max-hls-per-month`, a smaller AOI, or a shorter `--rolling-days`.
- Increase robustness: keep defaults; the script retries CMR queries and skips unreadable assets gracefully.

---

## Troubleshooting

- **Auth errors:** Check your `~/.netrc` or run `earthaccess.login()` once in a Python shell to verify.
- **RasterIO/GDAL errors:** Ensure `rasterio` can open remote TIFFs; a working GDAL install is required.
- **All NaNs for bands:** Too cloudy, small/empty AOI, or strict masks. Try enabling `--accept-water` (for water bodies), or increasing `--rolling-days` / `--max-hls-per-month`.

---

## Data sources & attribution

- **GEDI L4A AGBD** (ORNL DAAC, Mg/ha).  
- **Harmonized Landsat–Sentinel-2 (HLS) v2.0**: HLSS30 (Sentinel-2) + HLSL30 (Landsat-8/9).  
Please follow the original data providers’ citation and license terms when publishing results.
