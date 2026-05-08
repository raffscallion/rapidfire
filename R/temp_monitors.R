
# Because the AirMonitor archive is not always up-to-date, we will go directly to airnow.
# However, we still use AirMonitor for the full list of temporary monitors that may not be
# in airnow.

#' Acquire temporary monitor (AIRSIS/WRCC) daily PM2.5 data
#'
#' Downloads AIRSIS and WRCC temporary monitor data via the \pkg{AirMonitor} package,
#' selects the appropriate archive tier based on how recent the requested date is
#' (latest < 10 days, daily < 45 days, annual otherwise), computes 24-hour average
#' PM2.5, and returns a projected \code{SpatVector} saved to disk. Only monitors with
#' at least 16 hours of valid observations are retained.
#'
#' @param date A \code{Date} or date-coercible character string specifying the day
#'   to acquire.
#' @param states An optional character vector of two-letter state codes used to
#'   filter monitors (e.g., \code{c("CA", "OR", "WA")}). If \code{NULL}, all
#'   available monitors are retained.
#' @param timezone Local timezone string used to assign observations to the correct
#'   calendar day. Default is \code{"America/Los_Angeles"}.
#' @param output_path Directory path where the processed RDS file will be saved.
#'   Default is \code{"./processed_data/temp_monitors/"}.
#' @param crs Coordinate reference system string passed to \code{terra::project}.
#'   Default is \code{"EPSG:3395"}.
#'
#' @returns A \code{SpatVector} with one point per monitor and columns
#'   \code{monitorID} (device deployment ID), \code{Day} (date), \code{PM25}
#'   (24-hour mean PM2.5), \code{Hours} (number of hourly observations used),
#'   and \code{PM25_log} (log-transformed PM2.5). Returns \code{NULL} if no
#'   monitors with sufficient data are found. Also written to \code{output_path}
#'   as \code{temp_monitors_pm25_24hr_YYYY-MM-DD.RDS}.
#' @export
#'
#' @seealso \code{\link{monitors_combine}}, \code{\link{airnow_acquire}}
#'
#' @examples
#' \dontrun{
#' temp_monitors_acquire("2024-11-15", states = c("CA", "OR", "NV"))
#' }
temp_monitors_acquire <- function(date, states = NULL, timezone = "America/Los_Angeles",
                                  output_path = "./processed_data/temp_monitors/",
                                  crs = "EPSG:3395") {
  
  dt <- as.Date(date)
  
  if (Sys.Date() - dt < 10) {
    airsis <- AirMonitor::airsis_loadLatest()
    wrcc <- AirMonitor::wrcc_loadLatest()
  } else if (Sys.Date() - dt < 45) {
    airsis <- AirMonitor::airsis_loadDaily()
    wrcc <- AirMonitor::wrcc_loadDaily()
  } else {
    airsis <- AirMonitor::airsis_loadAnnual(lubridate::year(dt))
    wrcc <- AirMonitor::wrcc_loadAnnual(lubridate::year(dt))
  }

  mon <- AirMonitor::monitor_combine(airsis, wrcc) |>
    AirMonitor::monitor_filterDate(startdate = dt, enddate = dt,
                                   timezone = timezone)
  
  if (!is.null(states)) {
    mon <- mon  |>
      AirMonitor::monitor_filter(stateCode %in% states)
  }
  
  spat <- temp_monitors_preprocess(mon, crs, timezone)
  
  fname <- paste0("temp_monitors_pm25_24hr_", strftime(date, "%Y-%m-%d"), ".RDS")
  fname <- fs::path_join(c(output_path, fname))
  terra::saveRDS(spat, fname)
  spat

}

#' Convert an AirMonitor \code{mts_monitor} object to a projected SpatVector
#'
#' Pivots hourly monitor data to long format, joins site metadata, computes daily
#' mean PM2.5 in local time, and filters to days with at least 16 valid hourly
#' observations and positive mean values. Log-transforms the daily mean and returns
#' the result as a projected \code{SpatVector}.
#'
#' @param mon An \code{mts_monitor} object as returned by \pkg{AirMonitor} functions
#'   such as \code{AirMonitor::monitor_combine}.
#' @param crs Coordinate reference system string passed to \code{terra::project}.
#'   Default is \code{"EPSG:3395"}.
#' @param timezone Local timezone string used to assign hourly observations to the
#'   correct calendar day.
#'
#' @returns A \code{SpatVector} with columns \code{monitorID}, \code{Day},
#'   \code{PM25}, \code{Hours}, and \code{PM25_log}, or \code{NULL} if no monitors
#'   meet the data completeness threshold.
#'
#' @seealso \code{\link{temp_monitors_acquire}}
temp_monitors_preprocess <- function(mon, crs = "EPSG:3395", timezone = timezone) {
  
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