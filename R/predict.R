# Use an existing developed model to predict values at given new locations and dates

# Predict a single day across an output grid
#' Title
#'
#' @param dt 
#' @param grid_file 
#' @param model_file 
#' @param paths 
#'
#' @returns
#' @export
#'
#' @examples
predict_grid <- function(dt, grid_file, model_file, paths) {
  
  # Load template grid
  grid <- terra::rast(grid_file)
  
  # load monitoring data
  monitoring_file <- fs::path_join(c(paths[["monitors"]], 
                                     paste0("monitors_combined_",
                                            strftime(dt, "%Y-%m-%d"),
                                            ".RDS")))
  mon <- terra::readRDS(monitoring_file)
  
  models <- readRDS(model_file)
  
  # load met and smoke data from hrrr
  hrrr <- hrrr_stack(dt, paths[["hrrr"]])
  
  # Convert MASSDEN to log scale
  massdenlog <- log(hrrr$MASSDEN)
  names(massdenlog) <- "MASSDEN_log"
  hrrr <- terra::rast(list(massdenlog, hrrr))
  hrrr <- terra::subset(hrrr, which(names(hrrr)=="MASSDEN"), negate = TRUE)
  
  # interpolate to grid
  # mon_vgm <- monitors_variogram(mon)
  mon_vgm <- models$vgm_mon
  monitor_grid <- monitors_krige_grid(mon, grid, mon_vgm)
  ank <- monitor_grid[[1]]
  names(ank) <- "PM25_log_ANK"

  
  
  # load satellite data
  maiac_file <- fs::path_join(c(paths[["maiac"]],
                                paste0("MAIAC_processed_", strftime(dt, "%Y-%m-%d"), ".tif")))
  maiac <- terra::rast(maiac_file)
  
  # interpolate to grid
  maiac <- maiac_regrid(maiac, grid)
  
  # load sensor data
  ### Sensor data will be coming from CDPH? Need non-cdph way as well
  pa_file <- fs::path_join(c(paths[["purpleair"]],
                             paste0("purpleair_processed_",
                                    strftime(dt, "%Y-%m-%d"),
                                    ".RDS")))
  pa <- terra::readRDS(pa_file)
  
  # interpolate to grid
  #pa_vgm <- monitors_variogram(pa)
  pa_vgm <- models$vgm_pa
  pa_grid <- monitors_krige_grid(pa, grid, pa_vgm, nmax = 100)
  pak <- pa_grid[[1]]
  names(pak) <- "PM25_log_PAK"

  # combine data
  stack <- terra::rast(list(ank, pak, hrrr, maiac))

  # run model
  #model <- readRDS(model_file)$model
  model <- models$model
  res <- terra::predict(stack, model)
  names(res) <- "PM25_log_RF"

  # convert to linear
  linear <- exp(res)
  names(linear) <- "PM25_RF"
  
  # stack output (including input)
  final <- terra::rast(list(linear, res, stack))
  
  # Export file
  outfile <- fs::path_join(c(paths[["output"]], 
                             paste0("rapidfire_", 
                                    strftime(dt, format = "%Y-%m-%d"), ".tif")))
  
  terra::writeRaster(final, outfile, overwrite = TRUE)

  return(final)
  
}