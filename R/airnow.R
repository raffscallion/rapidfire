
# AirNow data via their file archive
# This seems more reliable than AirSensor, but doesn't include AIRSIS and other temporary data

#' Acquire AirNow daily PM2.5 monitor data
#'
#' Downloads the AirNow daily data file for a given date from the AirNow S3
#' archive, filters to 24-hour average PM2.5 observations, log-transforms the
#' values, and returns a projected \code{SpatVector} saved to disk. This source
#' covers regulatory monitors but does not include temporary monitors from AIRSIS
#' or WRCC (see \code{\link[rapidfire]{temp_monitors_acquire}} for those).
#'
#' @param date A \code{Date} or date-coercible character string specifying the day
#'   to acquire.
#' @param output_path Directory path where the processed RDS file will be saved.
#'   Default is \code{"./processed_data/airnow/"}.
#' @param crs Coordinate reference system string passed to \code{terra::project}.
#'   Default is \code{"EPSG:3395"}.
#'
#' @returns A \code{SpatVector} with one point per monitor and columns
#'   \code{monitorID} (AQS site identifier), \code{Day} (date), and
#'   \code{PM25_log} (log-transformed 24-hour average PM2.5). Also written to
#'   \code{output_path} as \code{airnow_pm25_24hr_YYYY-MM-DD.RDS}.
#' @export
#'
#' @seealso \code{\link{monitors_combine}}
#'
#' @examples
#' \dontrun{
#' airnow_acquire("2024-11-15")
#' }
airnow_acquire <- function(date, output_path = "./processed_data/airnow/",
                           crs = "EPSG:3395") {
  
  an_url <- "https://s3-us-west-1.amazonaws.com//files.airnowtech.org/airnow"
  
  dt <- as.Date(date)
  full_url <- paste(an_url, strftime(dt, "%Y"), strftime(dt, "%Y%m%d"), "daily_data_v2.dat",
                    sep = "/")
  
  airnow <- readr::read_delim(url(full_url), delim = "|", 
                              col_names = c("valid_date", "aqsid", "site_name", 
                                            "parameter_name", "reporting_units", "value",
                                            "averaging_period", "data_source", "aqi",
                                            "aqi_category", "latitude", "longitude",
                                            "full_aqsid"),
                              show_col_types = FALSE) |>
    mutate(Day = as.Date(valid_date, format = "%m/%d/%y")) |>
    filter(parameter_name == "PM2.5-24hr",
           value != -999) |>
    mutate(PM25_log = log(value)) |>
    select(monitorID=full_aqsid, Day, longitude, latitude, PM25_log)
  
  # make spatial and project to desired coordinates
  spat <- terra::vect(airnow, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  spat <- terra::project(spat, crs)
  
  fname <- paste0("airnow_pm25_24hr_", strftime(date, "%Y-%m-%d"), ".RDS")
  fname <- fs::path_join(c(output_path, fname))
  terra::saveRDS(spat, fname)
  spat
  
}