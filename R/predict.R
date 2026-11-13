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
#' @param grid A \code{SpatRaster} defining the prediction grid (CRS, extent, and
#'   resolution).
#' @param mon A \code{SpatVector} of combined regulatory and temporary monitor
#'   observations for the day, with a \code{PM25_log} column.
#' @param models A named list with elements:
#'   \describe{
#'     \item{model}{A fitted random forest model (e.g., from \code{\link{develop_model}}).}
#'     \item{vgm_mon}{A \code{variogramModel} for regulatory monitor kriging.}
#'     \item{vgm_pa}{A \code{variogramModel} for PurpleAir sensor kriging
#'       (\code{nmax = 100}).}
#'   }
#' @param hrrr A \code{SpatRaster} of HRRR meteorological and smoke layers for the
#'   day (e.g., from \code{\link{hrrr_stack}}).
#' @param maiac A \code{SpatRaster} of processed MAIAC AOD for the day. It is
#'   regridded to \code{grid} internally.
#' @param pa A \code{SpatVector} of processed PurpleAir sensors for the day (e.g.,
#'   from \code{\link{pa_preprocess}}).
#'
#' @returns A multi-layer \code{SpatRaster} on the \code{grid} grid with layers:
#'   \describe{
#'     \item{PM25_RF}{Random forest PM2.5 prediction (µg/m³, linear scale).}
#'     \item{PM25_log_RF}{Random forest PM2.5 prediction (log scale).}
#'     \item{PM25_log_ANK}{Kriged regulatory monitor estimate (log scale).}
#'     \item{PM25_log_PAK}{Kriged PurpleAir estimate (log scale).}
#'     \item{MASSDEN_log, HPBL, RH, UGRD, VGRD, TMP, APCP}{HRRR meteorological
#'       predictors.}
#'     \item{MAIAC_AOD}{MAIAC aerosol optical depth.}
#'   }
#'   Use \code{\link{predict_grid_write}} to save the result to disk.
#' @export
#'
#' @seealso \code{\link{predict_grid_write}}, \code{\link{hrrr_stack}},
#'   \code{\link{monitors_krige_grid}}, \code{\link{maiac_regrid}}
#'
#' @examples
#' \dontrun{
#' result <- predict_grid(
#'   dt = "2024-11-15",
#'   grid = terra::rast("./grid/prediction_grid.tif"),
#'   mon = mon,
#'   models = readRDS("./models/rapidfire_model.RDS"),
#'   hrrr = hrrr_stack("2024-11-15", "./processed_data/HRRR/"),
#'   maiac = terra::rast("./processed_data/MAIAC/MAIAC_processed_2024-11-15.tif"),
#'   pa = pa
#' )
#' }
predict_grid <- function(dt, grid, mon, models, hrrr, maiac, pa) {
  
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
  
  # interpolate maiac to grid
  maiac <- maiac_regrid(maiac, grid)
  
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
  
}

#' Write predicted PM2.5 grids to a GeoTIFF
#'
#' Saves the output of \code{\link{predict_grid}} as
#' \code{rapidfire_YYYY-MM-DD.tif}, overwriting any existing file.
#'
#' @param final_grids A \code{SpatRaster} returned by \code{\link{predict_grid}}.
#' @param dt A \code{Date} or date-coercible character string for the predicted day,
#'   used in the file name.
#' @param output_path Directory where the GeoTIFF will be saved.
#'
#' @returns Called for its side effect of writing a file; returns
#'   \code{final_grids} invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' predict_grid_write(result, "2024-11-15", "./output/")
#' }
predict_grid_write <- function(final_grids, dt, output_path) {
  # Export file
  outfile <- fs::path_join(c(output_path,
                             paste0("rapidfire_",
                                    strftime(dt, format = "%Y-%m-%d"), ".tif")))

  terra::writeRaster(final_grids, outfile, overwrite = TRUE)
  
}
