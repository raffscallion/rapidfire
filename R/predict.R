# Use an existing developed model to predict values at given new locations and dates

#' Predict gridded PM2.5 concentrations for a single day
#'
#' Assembles all predictor layers for a given date — kriged regulatory and temprorary monitor
#' estimates (ANK), kriged PurpleAir estimates (PAK), HRRR meteorological variables,
#' and MAIAC AOD — and applies a pre-trained random forest model to generate a
#' gridded PM2.5 surface. Both log-scale and linear-scale predictions are included
#' in the output stack alongside all input predictor layers.
#'
#' @param dt A \code{Date} or date-coercible character string specifying the day to
#'   predict.
#' @param grid_file Path to a GeoTIFF file defining the prediction grid (CRS, extent,
#'   and resolution). Loaded with \code{terra::rast}.
#' @param model_file Path to an RDS file containing a named list with elements:
#'   \describe{
#'     \item{model}{A fitted random forest model (e.g., from \code{\link{develop_model}}).}
#'     \item{vgm_mon}{A \code{variogramModel} for regulatory monitor kriging.}
#'     \item{vgm_pa}{A \code{variogramModel} for PurpleAir sensor kriging
#'       (\code{nmax = 100}).}
#'   }
#' @param paths A named list of directory paths with the following elements:
#'   \describe{
#'     \item{monitors}{Directory containing combined monitor RDS files
#'       (\code{monitors_combined_YYYY-MM-DD.RDS}).}
#'     \item{hrrr}{Directory containing preprocessed HRRR GeoTIFF files.}
#'     \item{maiac}{Directory containing preprocessed MAIAC GeoTIFF files.}
#'     \item{purpleair}{Directory containing preprocessed PurpleAir RDS files.}
#'     \item{output}{Directory where the output GeoTIFF will be saved.}
#'   }
#'
#' @returns A multi-layer \code{SpatRaster} on the \code{grid_file} grid with layers:
#'   \describe{
#'     \item{PM25_RF}{Random forest PM2.5 prediction (µg/m³, linear scale).}
#'     \item{PM25_log_RF}{Random forest PM2.5 prediction (log scale).}
#'     \item{PM25_log_ANK}{Kriged regulatory monitor estimate (log scale).}
#'     \item{PM25_log_PAK}{Kriged PurpleAir estimate (log scale).}
#'     \item{MASSDEN_log, HPBL, RH, UGRD, VGRD, TMP, APCP}{HRRR meteorological
#'       predictors.}
#'     \item{MAIAC_AOD}{MAIAC aerosol optical depth.}
#'   }
#'   The raster is also written to \code{paths[["output"]]} as
#'   \code{rapidfire_YYYY-MM-DD.tif}.
#' @export
#'
#' @seealso \code{\link{hrrr_stack}}, \code{\link{monitors_krige_grid}},
#'   \code{\link{maiac_regrid}}
#'
#' @examples
#' \dontrun{
#' result <- predict_grid(
#'   dt = "2024-11-15",
#'   grid_file = "./grid/prediction_grid.tif",
#'   model_file = "./models/rapidfire_model.RDS",
#'   paths = list(
#'     monitors  = "./processed_data/monitors/",
#'     hrrr      = "./processed_data/HRRR/",
#'     maiac     = "./processed_data/MAIAC/",
#'     purpleair = "./processed_data/purpleair/",
#'     output    = "./output/"
#'   )
#' )
#' }
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