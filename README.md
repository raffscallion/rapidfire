# rapidfire

**R**elatively **A**ccurate **P**articulate **I**nformation **D**erived **F**rom **I**nformation **R**etrieved **E**asily

`rapidfire` estimates daily surface-level PM2.5 concentrations by fusing routine and temporary air
quality monitors, low-cost sensors, satellite aerosol retrievals, and
meteorological model output using a random forest regression. It can be used for
operational near-real-time (previous day) estimates or for retrospective analysis.

---

## Data Sources

| Source | Variable | Function |
|---|---|---|
| [AirNow](https://www.airnow.gov) file archive | Hourly PM2.5 (air agency monitors) | `airnow_acquire()` |
| [AIRSIS / WRCC](https://github.com/MazamaScience/AirMonitor) (via AirMonitor) | Hourly PM2.5 (temporary monitors) | `temp_monitors_acquire()` |
| [PurpleAir](https://www.purpleair.com) API | Estimated PM2.5 (low-cost sensors) | `pa_acquire_realtime()`, `pa_preprocess()` |
| [HRRR-Smoke](https://rapidrefresh.noaa.gov/hrrr/) (NOAA S3) | Boundary layer height, temperature, winds, humidity, smoke mass density, precipitation | `hrrr_acquire()` |
| [MAIAC AOD](https://lpdaac.usgs.gov/products/mcd19a2v061/) MCD19A2 (NASA EarthData) | 470 nm aerosol optical depth | `maiac_acquire()`, `maiac_preprocess()` |

---

## Installation

```r
# install.packages("pak")
pak::pkg_install("raffscallion/rapidfire")
```

### Prerequisites

- **NASA EarthData account** — required for MAIAC satellite data.
  Register at <https://urs.earthdata.nasa.gov>. Store credentials as environment
  variables `EARTHDATA_USER` and `EARTHDATA_PASSWORD` (e.g., in `.Renviron`).
- **PurpleAir API key** — required for sensor data acquisition.
  Request a key at <https://community.purpleair.com/t/creating-api-keys/3951>.
  Store as `PURPLEAIR_READ_KEY` in `.Renviron`. Alternately, a source of PurpleAir data.

---

## Workflow

The package supports three main stages: **model development**, **daily
prediction**, and **validation**. Each stage reads from and writes
to a set of directories passed as a named `paths` list. You can change the paths to 
suit your environment, but the names must be as listed here.

### 1. Set up paths

```r
paths <- list(
  airnow          = "./data/airnow/",
  temp_monitors   = "./data/temp_monitors/",
  monitors        = "./data/monitors/",
  hrrr            = "./data/hrrr/",
  maiac           = "./data/maiac/",
  purpleair       = "./data/purpleair/",
  output          = "./output/",
  cross_validation = "./output/cross_validation/"
)
```

### 2. Acquire daily inputs

Run these steps for each date you want to process. For operational use, schedule
them daily (e.g., with `cronR` or `taskscheduleR`).

```r
library(rapidfire)

date <- "2024-11-15"

# Load the grid definition as created by rapidfire::grid_state()
grid <- terra::rast("./grid/prediction_grid.tif")

# Regulatory monitors (AirNow)
airnow_acquire(date, output_path = paths[["airnow"]])

# Temporary monitors (AIRSIS / WRCC)
temp_monitors_acquire(date, output_path = paths[["temp_monitors"]])

# Combine monitor sources and optionally clip to an extent
# Here I loaded a pregenerated extent of 100 km around CA saved as a SpatVector
extent <- terra::readRDS("./extents/ca_100km_buffer.RDS")
monitors_combine(date,
  input_path_airnow        = paths[["airnow"]],
  input_path_temp_monitors = paths[["temp_monitors"]],
  output_path              = paths[["monitors"]],
  extent                   = extent)

# HRRR-Smoke meteorology (raw GRIB2 files deleted after processing)
hrrr_acquire(date, outpath = paths[["hrrr"]], extent_grid = grid,
             delete_files = TRUE)

# MAIAC AOD satellite data
maiac_files <- maiac_acquire(date,
  user     = Sys.getenv("EARTHDATA_USER"),
  password = Sys.getenv("EARTHDATA_PASSWORD"),
  outpath  = paths[["maiac"]])
maiac_preprocess(date, input_path = paths[["maiac"]],
                 output_path = paths[["maiac"]])
fs::file_delete(maiac_files)  # remove intermediate HDF files

# PurpleAir low-cost sensors
sensors <- readRDS("./data/purpleair/purpleair_sensors.RDS")
pa_preprocess(date,
  duckdb_path = "./data/purpleair/purpleair.duckdb",
  sensors     = sensors,
  output_path = paths[["purpleair"]])
```

For PurpleAir, a local DuckDB database must first be populated using
`pa_acquire_realtime()` on a schedule (e.g., hourly). The sensor metadata
file can be built once with `pa_find_sensors()`.

### 3. Develop a model

Train a model over a date range. This only needs to be re-run periodically
(e.g., monthly or seasonally) as conditions change.

```r
model_list <- develop_model(
  dt1   = "2024-10-01",
  dt2   = "2024-11-30",
  paths = paths
)
saveRDS(model_list, "./models/rapidfire_model.RDS")
```

`develop_model()` returns a list containing the fitted random forest model and
the kriging variogram models needed for prediction. These are saved together so
that `predict_grid()` and `daily_cross_validate()` have everything they need in
one file.

### 4. Generate daily predictions

```r
result <- predict_grid(
  dt         = date,
  grid_file  = "./grid/prediction_grid.tif",
  model_file = "./models/rapidfire_model.RDS",
  paths      = paths
)
```

Output is a multi-layer GeoTIFF (`rapidfire_YYYY-MM-DD.tif`) containing the
PM2.5 prediction (`PM25_RF`, µg/m³) alongside the log-scale prediction and all
input predictor layers.

### 5. Validate predictions

Leave-one-out cross-validation can be run alongside daily predictions to track
ongoing model performance:

```r
cv_results <- daily_cross_validate(
  dt         = date,
  model_file = "./models/rapidfire_model.RDS",
  paths      = paths
)
```

Results are saved as `rapidfire_l-o-o_YYYY-MM-DD.csv`.

---

## Function Reference

### Data acquisition

| Function | Description |
|---|---|
| `airnow_acquire()` | Download AirNow 24-hr PM2.5 data |
| `temp_monitors_acquire()` | Download AIRSIS/WRCC temporary monitor data |
| `pa_find_sensors()` | Get PurpleAir sensor metadata for a bounding box |
| `pa_acquire_realtime()` | Fetch current PurpleAir readings into a DuckDB database |
| `hrrr_acquire()` | Download and preprocess HRRR-Smoke data |
| `maiac_acquire()` | Download MAIAC AOD granules from NASA EarthData |
| `maiac_preprocess()` | Merge MAIAC tiles and fill spatial gaps |

### Data preparation

| Function | Description |
|---|---|
| `monitors_combine()` | Combine AirNow and temporary monitor data |
| `pa_preprocess()` | Compute 24-hr PurpleAir averages from a DuckDB database |

### Modeling

| Function | Description |
|---|---|
| `develop_model()` | Train a random forest model over a date range |
| `predict_grid()` | Generate a gridded PM2.5 surface for a single day |
| `daily_cross_validate()` | Leave-one-out cross-validation at monitor locations |

### Utilities

| Function | Description |
|---|---|
| `monitors_variogram()` | Fit a single-day variogram to monitor data |
| `monitors_variogram_pooled()` | Fit a pooled multi-day variogram |
| `monitors_krige_grid()` | Krige PM2.5 values onto a prediction grid |
| `monitors_krige_points()` | Krige PM2.5 values at point locations |
| `monitors_split()` | Split monitor data into train/test sets |
| `grid_state()` | Create an output modeling grid based on a spatial extent (such as a state) |
| `pa_check_api_key()` | Verify a PurpleAir API key |
| `hrrr_stack()` | Load preprocessed HRRR variables into a stacked SpatRaster |

---

## Citation

If you use `rapidfire` in your research, please cite:

> Raffuse, S., O'Neill, S., and Schmidt, R.: A model for rapid PM2.5 exposure
> estimates in wildfire conditions using routinely available data: rapidfire
> v0.1.3, *Geosci. Model Dev.*, 17, 381–397,
> https://doi.org/10.5194/gmd-17-381-2024, 2024.

---

## License

GPL (>= 3)
