
# Functions for handling monitor data, updated to use terra rather than sp

# Using updated AirMonitor functions, retrieve data from all of AirNow, AIRSIS, and WRCC
monitors_acquire <- function(start, end, states = "CA", output_path = "./processed_data/monitors/",
                             timezone = "America/Los_Angeles", crs = "EPSG:3395") {
  
  start <- as.Date(start)
  end <- as.Date(end)
  #end <- end + 1

  # Load data from the AirMonitor archive
  mon <- AirMonitor::monitor_load(start, end + 1, epaPreference = "epa_aqs",
                                  timezone = timezone)

  # Subset to the desired states  
  if (!is.na(states)) {
    mon <- AirMonitor::monitor_filter(mon, stateCode %in% states)
  }
  
  # Calculated 24-hr average PM2.5 (and log PM2.5), make spatial, and convert projections
  spat <- monitors_preprocess(mon, crs = crs, timezone = timezone)

  # Export as daily RDS files
  dates <- seq.Date(start, end, by = "1 day")
  export_rds <- function(date) {
    one_day <- terra::subset(spat, spat$Day == date)
    fname <- fs::path_join(c(output_path, 
                             paste0("monitors_processed_", 
                                    strftime(date, "%Y-%m-%d"), ".RDS")))
    terra::saveRDS(one_day, fname)
    return(fname)
  }
  files <- purrr::map_chr(dates, export_rds, .progress = TRUE)
  
}

# Convert a mts_monitor object to a spatvector in a standard projection. Also calculate
# 24-hr average and log.
monitors_preprocess <- function(mon, crs = "EPSG:3395", timezone = timezone) {

  df <- mon$data |>
    tidyr::pivot_longer(names_to = "monitorID", values_to = "PM25", -datetime)
    #tidyr::gather("monitorID", "PM25", -datetime)
  meta <- mon$meta |>
    select(monitorID=deviceDeploymentID, longitude, latitude, stateCode, deviceType)
  
  df <- left_join(df, meta, by = "monitorID") |>
    filter(!is.na(PM25)) |>
    mutate(LocalTime = lubridate::with_tz(datetime, timezone),
           Day = lubridate::floor_date(LocalTime, "days")) |>
    summarise(PM25 = mean(PM25, na.rm = FALSE),
              Hours = n(),
              .by = c(monitorID, Day)) |>
    filter(Hours >= 16,
           PM25 > 0) |>
    mutate(PM25_log = log(PM25))
  
  # Attach lat/lon from metadata
  sites <- meta |>
    select(monitorID, longitude, latitude) |>
    distinct()
  df <- df |>
    inner_join(sites, by = "monitorID")
  
  if (nrow(df) == 0) {
    "No monitor data for this range..."
    return(NULL)
  }

  # make spatial and project to desired coordinates
  spat <- terra::vect(df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  spat <- terra::project(spat, crs)
  
}

# Create variogram from spatvectors
monitors_variogram <- function(mon, cutoff = NULL, width = 15000) {
  
  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  df <- terra::as.data.frame(mon, geom = "XY")
  vgm_mod <- gstat::vgm(model = "Gau", nugget = 0.02)
  vgm_data <- gstat::variogram(PM25_log ~  1, locations = ~x+y, data = df,
                               cutoff = cutoff, width = width)
  v <- gstat::fit.variogram(vgm_data, vgm_mod)
  
}

# Given input measurements, output locations, and a variogram, predict values at the
# output locations
monitors_krige_points <- function(mon, locs, vgm) {
  
  pts_sf <- sf::st_as_sf(mon)
  ok <- gstat::krige(PM25_log ~ 1, locations = pts_sf, newdata = locs,
                     model = vgm, maxdist = 100000)
  
}

#Given input measurements, a prediction grid, and a variogram, predict values on the grid
monitors_krige_grid <- function(mon, grid, vgm) {
  
  ## YOU ARE HERE
  res <- terra::interpolate(grid, vgm)
  
}

test_locs <- tibble(id = 1:3,
                    lon = c(-121.72693, -122.58036, -117.91731),
                    lat = c(38.54058, 37.92791, 33.81088))
test_locs <- sf::st_as_sf(test_locs, coords = c("lon", "lat"), crs = sf::st_crs(4326))
test_locs <- sf::st_transform(test_locs, crs = sf::st_crs(3395))
