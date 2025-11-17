
# Because the AirMonitor archive is not always up-to-date, we will go directly to airnow.
# However, we still use AirMonitor for the full list of temporary monitors that may not be
# in airnow.

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

# Convert a mts_monitor object to a spatvector in a standard projection. Also calculate
# 24-hr average and log.
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