
# AirNow data via their api
# This seems more reliable than AirSensor, but doesn't include AIRSIS and other temporary data


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