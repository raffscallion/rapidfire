
#' Leave-one-out cross-validation of daily PM2.5 predictions
#'
#' Evaluates model performance for a single day using leave-one-out cross-validation
#' across all regulatory monitors. For each monitor, its observation is withheld and
#' kriging is used to interpolate a PM2.5 estimate at that location from the remaining
#' monitors. Predictor variables (HRRR meteorology, MAIAC AOD, kriged PurpleAir) are
#' then extracted at each monitor location and fed into the pre-trained random forest
#' model to generate a final PM2.5 prediction.
#'
#' @param dt A \code{Date} or date-coercible character string specifying the day to
#'   validate.
#' @param mon A \code{SpatVector} of monitor observations for the day, with
#'   \code{monitorID}, \code{Day}, \code{PM25_log}, and \code{source} columns.
#' @param models A named list with elements:
#'   \describe{
#'     \item{model}{A fitted random forest model (e.g., from \code{\link{develop_model}}).}
#'     \item{vgm_mon}{A \code{variogramModel} for regulatory monitor kriging.}
#'     \item{vgm_pa}{A \code{variogramModel} for PurpleAir sensor kriging.}
#'   }
#' @param hrrr A \code{SpatRaster} of HRRR layers for the day (e.g., from
#'   \code{\link{hrrr_stack}}).
#' @param maiac A \code{SpatRaster} of processed MAIAC AOD for the day.
#' @param pa A \code{SpatVector} of processed PurpleAir sensors for the day.
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
#'   Use \code{\link{daily_cross_validate_write}} to save the result as a CSV.
#' @export
#'
#' @seealso \code{\link{daily_cross_validate_write}},
#'   \code{\link{monitors_krige_points}}, \code{\link{hrrr_stack}}
#'
#' @examples
#' \dontrun{
#' results <- daily_cross_validate(
#'   dt = "2024-11-15",
#'   mon = mon,
#'   models = readRDS("./models/rapidfire_model.RDS"),
#'   hrrr = hrrr_stack("2024-11-15", "./processed_data/HRRR/"),
#'   maiac = terra::rast("./processed_data/MAIAC/MAIAC_processed_2024-11-15.tif"),
#'   pa = pa
#' )
#' }
daily_cross_validate <- function(dt, mon, models, hrrr, maiac, pa) {
  
  
  # Remove any non-finite values
  mon <- mon[is.finite(mon$PM25_log),]
  
  # For each site, interpolate without it's data
  ids <- mon$monitorID
  
  # load the variogram
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
  
  # Convert MASSDEN to log scale
  massdenlog <- log(hrrr$MASSDEN)
  names(massdenlog) <- "MASSDEN_log"
  hrrr <- terra::rast(list(massdenlog, hrrr))
  hrrr <- terra::subset(hrrr, which(names(hrrr)=="MASSDEN"), negate = TRUE)
  hrrr_extracted <- terra::extract(hrrr, mon)
  
  # Extract MAIAC at monitor points
  maiac_crs <- terra::crs(maiac)
  mon_maiac <- terra::project(mon, maiac_crs)
  maiac_extracted <- terra::extract(maiac, mon_maiac)
  
  # Get median value to fill NAs
  maiac_median <- terra::global(maiac, fun = median, na.rm = TRUE) %>%
    .$global
  
  # Use median to fill NAs
  maiac_extracted <- maiac_extracted |>
    tidyr::replace_na(list(MAIAC_AOD = maiac_median))
  
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

  return(results)
}

#' Write leave-one-out cross-validation results to CSV
#'
#' Saves the output of \code{\link{daily_cross_validate}} as
#' \code{rapidfire_l-o-o_YYYY-MM-DD.csv}.
#'
#' @param results A data frame returned by \code{\link{daily_cross_validate}}.
#' @param dt A \code{Date} or date-coercible character string for the validated day,
#'   used in the file name.
#' @param output_path Directory where the CSV will be saved.
#'
#' @returns Called for its side effect of writing a file; returns
#'   \code{results} invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' daily_cross_validate_write(results, "2024-11-15", "./output/cross_validation/")
#' }
daily_cross_validate_write <- function(results, dt, output_path) {
  
  outfile <- fs::path_join(c(output_path,
                             paste0("rapidfire_l-o-o_",
                                    strftime(dt, format = "%Y-%m-%d"), ".csv")))
  readr::write_csv(results, outfile)
}
