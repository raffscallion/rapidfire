
#' Leave-one-out cross-validation of daily PM2.5 predictions
#'
#' Evaluates model performance for a single day using leave-one-out cross-validation
#' across all regulatory monitors. For each monitor, its observation is withheld and
#' kriging is used to interpolate a PM2.5 estimate at that location from the remaining
#' monitors. Predictor variables (HRRR meteorology, MAIAC AOD, kriged PurpleAir) are
#' then extracted at each monitor location and fed into the pre-trained random forest
#' model to generate a final PM2.5 prediction. Results are saved as a CSV.
#'
#' @param dt A \code{Date} or date-coercible character string specifying the day to
#'   validate.
#' @param model_file Path to an RDS file containing a named list with elements:
#'   \describe{
#'     \item{model}{A fitted random forest model (e.g., from \code{\link{develop_model}}).}
#'     \item{vgm_mon}{A \code{variogramModel} for regulatory monitor kriging.}
#'     \item{vgm_pa}{A \code{variogramModel} for PurpleAir sensor kriging.}
#'   }
#' @param paths A named list of directory paths with the following elements:
#'   \describe{
#'     \item{monitors}{Directory containing combined monitor RDS files
#'       (\code{monitors_combined_YYYY-MM-DD.RDS}).}
#'     \item{hrrr}{Directory containing preprocessed HRRR GeoTIFF files.}
#'     \item{maiac}{Directory containing preprocessed MAIAC GeoTIFF files.}
#'     \item{purpleair}{Directory containing preprocessed PurpleAir RDS files.}
#'     \item{cross_validation}{Directory where the output CSV will be saved.}
#'   }
#'
#' @returns A data frame with one row per monitor containing:
#'   \describe{
#'     \item{monitorID}{Monitor identifier.}
#'     \item{PM25_log}{Observed log-transformed PM2.5.}
#'     \item{PM25_log_ANK}{Leave-one-out kriged estimate from regulatory monitors.}
#'     \item{PM25_log_PAK}{Kriged estimate from PurpleAir sensors.}
#'     \item{MAIAC_AOD}{MAIAC aerosol optical depth extracted at the monitor.}
#'     \item{HPBL, RH, UGRD, VGRD, TMP, MASSDEN_log, APCP}{HRRR meteorological
#'       variables extracted at the monitor.}
#'     \item{PM25_log_RF}{Random forest PM2.5 prediction.}
#'   }
#'   The data frame is also written to \code{paths[["cross_validation"]]} as
#'   \code{rapidfire_l-o-o_YYYY-MM-DD.csv}.
#' @export
#'
#' @seealso \code{\link{monitors_krige_points}}, \code{\link{hrrr_stack}}
#'
#' @examples
#' \dontrun{
#' results <- daily_cross_validate(
#'   dt = "2024-11-15",
#'   model_file = "./models/rapidfire_model.RDS",
#'   paths = list(
#'     monitors       = "./processed_data/monitors/",
#'     hrrr           = "./processed_data/HRRR/",
#'     maiac          = "./processed_data/MAIAC/",
#'     purpleair      = "./processed_data/purpleair/",
#'     cross_validation = "./output/cross_validation/"
#'   )
#' )
#' }
daily_cross_validate <- function(dt, model_file, paths) {
  
  # load monitoring data
  monitoring_file <- fs::path_join(c(paths[["monitors"]], 
                                     paste0("monitors_combined_",
                                            strftime(dt, "%Y-%m-%d"),
                                            ".RDS")))
  mon <- terra::readRDS(monitoring_file)

  # Remove any non-finite values
  mon <- mon[is.finite(mon$PM25_log),]
  
  # For each site, interpolate without it's data
  ids <- mon$monitorID
  
    # load the variogram
  models <- readRDS(model_file)
  mon_vgm <- models$vgm_mon
  
  leave_out_interpolate <- function(id) {
    
    mon_sub <- mon[mon$monitorID != id,]
    mon_one <- mon[mon$monitorID == id,]
    loc <- terra::crds(mon_one)
    
    res <- monitors_krige_points(mon_sub, mon_one, mon_vgm)
    df <- tibble(monitorID = id,
                 PM25_log_ANK = res$var1.pred)
  }
  
  interpolated_monitors <- purrr::map(ids, leave_out_interpolate, 
                                      .progress = "Leave-one-out interpolation") |>
    purrr::list_rbind()
  
  # load met and smoke data from hrrr
  hrrr <- hrrr_stack(dt, paths[["hrrr"]])
  
  # Convert MASSDEN to log scale
  massdenlog <- log(hrrr$MASSDEN)
  names(massdenlog) <- "MASSDEN_log"
  hrrr <- terra::rast(list(massdenlog, hrrr))
  hrrr <- terra::subset(hrrr, which(names(hrrr)=="MASSDEN"), negate = TRUE)
  hrrr_extracted <- terra::extract(hrrr, mon)
  
  # load satellite data
  maiac_file <- fs::path_join(c(paths[["maiac"]],
                                paste0("MAIAC_processed_", strftime(dt, "%Y-%m-%d"), ".tif")))
  maiac <- terra::rast(maiac_file)
  
  # Extract at monitor points
  maiac_crs <- terra::crs(maiac)
  mon_maiac <- terra::project(mon, maiac_crs)
  maiac_extracted <- terra::extract(maiac, mon_maiac)
  
  # load sensor data
  ### Sensor data will be coming from CDPH? Need non-cdph way as well
  pa_file <- fs::path_join(c(paths[["purpleair"]],
                             paste0("purpleair_processed_",
                                    strftime(dt, "%Y-%m-%d"),
                                    ".RDS")))
  pa <- terra::readRDS(pa_file)

  # interpolate to monitor points
  pa_vgm <- models$vgm_pa
  pak <- monitors_krige_points(pa, mon, pa_vgm, nmax = 100)

  # Compile and predict
  model_inputs <- as.data.frame(mon) |>
    select(monitorID, Day, PM25_log, source) |>
    left_join(interpolated_monitors, by = "monitorID")
  model_inputs$PM25_log_PAK <- pak$var1.pred
  model_inputs$MAIAC_AOD <- maiac_extracted$MAIAC_AOD
  model_inputs <- model_inputs |>
    bind_cols(select(hrrr_extracted, -ID))
  
  model_inputs <- model_inputs[complete.cases(model_inputs),]
  
  model <- models$model
  pred <- predict(model, model_inputs)
    
  # combine with observed inputs and return
  results <- model_inputs
  results$PM25_log_RF <- pred
  
  # Export file
  outfile <- fs::path_join(c(paths[["cross_validation"]], 
                             paste0("rapidfire_l-o-o_", 
                                    strftime(dt, format = "%Y-%m-%d"), ".csv")))
  readr::write_csv(results, outfile)
  
  return(results)
}


