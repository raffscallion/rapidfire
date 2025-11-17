
# Functions for handling monitor data, updated to use terra rather than sp

# Combine data from airnow and temporary monitors from AirMonitors, clip to an extent, and
# remove any null data
monitors_combine <- function(date, input_path_airnow, input_path_temp_monitors, 
                             output_path = "./processed_data/monitors/",
                             extent = NULL) {
  
  dt <- as.Date(date)
  airnow_filename <- paste0("airnow_pm25_24hr_", strftime(dt, "%Y-%m-%d"), ".RDS")
  an <- terra::readRDS(fs::path_join(c(input_path_airnow, airnow_filename)))
  an$source <- "airnow"
  
  tempm_filename <- paste0("temp_monitors_pm25_24hr_", strftime(dt, "%Y-%m-%d"), ".RDS")
  tempm <- terra::readRDS(fs::path_join(c(input_path_temp_monitors, tempm_filename)))
  tempm$source <- "airsis_wrcc"

  mon <- rbind(an, tempm)
  
  # If wanted, limit to monitors within the provided extent, which is a SpatVector or SpatExtent
  if (!is.null(extent)) {
    if (!identical(terra::crs(mon), terra::crs(extent))) {
      extent <- terra::project(extent, mon)      
    }
    mon <- mon[extent]
  }

  # Further modeling is not tolerant of missing data, so remove any NAs
  mon <- mon[!is.na(mon$PM25_log),]
  
  
  fname <- paste0("monitors_combined_", strftime(dt, "%Y-%m-%d"), ".RDS")
  fname <- fs::path_join(c(output_path, fname))
  terra::saveRDS(mon, fname)
  mon
  
}


# Create variogram from spatvectors
monitors_variogram <- function(mon, cutoff = NULL, width = 15000) {
  
  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  df <- terra::as.data.frame(mon, geom = "XY")
  # make sure all values are finite
  df <- dplyr::filter(df, is.finite(PM25_log))
  vgm_mod <- gstat::vgm(model = "Gau", nugget = 0.02)
  vgm_data <- gstat::variogram(PM25_log ~  1, locations = ~x+y, data = df,
                               cutoff = cutoff, width = width)
  v <- gstat::fit.variogram(vgm_data, vgm_mod)
  
}

# Given input measurements, output locations, and a variogram, predict values at the
# output locations
monitors_krige_points <- function(mon, locs, vgm, maxdist = Inf) {
  
  pts_sf <- sf::st_as_sf(mon)
  ok <- gstat::krige(PM25_log ~ 1, locations = pts_sf, newdata = locs,
                     model = vgm, maxdist = maxdist)
  
}

#Given input measurements, a prediction grid, and a variogram, predict values on the grid
monitors_krige_grid <- function(mon, grid, vgm, maxdist = Inf) {
  
  df <- terra::as.data.frame(mon, geom = "XY")
  # make sure all values are finite
  df <- dplyr::filter(df, is.finite(PM25_log))
  ok <- gstat::gstat(NULL, "PM25_log", PM25_log ~ 1, locations = ~x+y, model = vgm, 
                     maxdist = maxdist, data = df)
  res <- terra::interpolate(grid, ok)
  
}
 
# test_locs <- tibble(id = 1:3,
#                     lon = c(-121.72693, -122.58036, -117.91731),
#                     lat = c(38.54058, 37.92791, 33.81088))
# test_locs <- sf::st_as_sf(test_locs, coords = c("lon", "lat"), crs = sf::st_crs(4326))
# test_locs <- sf::st_transform(test_locs, crs = sf::st_crs(3395))
