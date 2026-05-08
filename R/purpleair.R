#' Take realtime purpleair data from a duckdb file and create 24-hr averages ready for
#' interpolation
#'
#' @param dt
#' @param duckdb_path
#' @param sensors
#' @param location_types
#' @param timezone
#' @param clean_outliers
#' @param crs
#' @param output_path
#'
#' @returns
#' @export
#'
#' @examples
pa_preprocess <- function(dt, duckdb_path, sensors, location_types = 0,
                          timezone = "America/Los_Angeles",
                          clean_outliers = TRUE,
                          crs = "EPSG:3395",
                          output_path = "./processed_data/purpleair/") {

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_path)
  
  # Pull days around the date from the database before converting to correct tz
  start <- as.Date(dt) - 1
  end <- as.Date(dt) + 2
  df <- tbl(con, "pa_realtime") |>
    filter(confidence > 60,
           last_seen > start,
           last_seen < end) |>
    collect()
  
  DBI::dbDisconnect(con)

  # Convert to local time and filter to requested date
  df <- df |>
    mutate(LocalTime = lubridate::with_tz(last_seen, tzone = timezone),
           LocalDay = lubridate::floor_date(LocalTime, "days")) |>
    filter(LocalDay == dt) |>
    select(sensor_index, LocalTime, LocalDay, pm2.5, confidence)
  
  # log-transform and calculate 24-hr average
  # Convert 0 to 0.01 first
  avg <- df |>
    mutate(pm2.5 = if_else(pm2.5 == 0, 0.01, pm2.5),
           PM25_log = log(pm2.5)) |>
    summarise(PM25_log = mean(PM25_log, na.rm = TRUE),
              n = n(),
              .by = sensor_index)
  
  # Join to geolocation and type info from sensors df
  sensors <- sensors |>
    select(sensor_index, latitude, longitude, location_type) |>
    mutate(across(everything(), as.numeric))
  
  avg <- avg |>
    left_join(sensors, by = "sensor_index") |>
    filter(location_type %in% location_types)  

  # filter based on 75% completeness assuming a sampling rate of once per 30 minutes
  avg <- avg |>
    filter(n >= 36) |>
    select(-n)

  # convert to spatial
  pa_sp <- terra::vect(avg, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # project to desired coordinates
  pa_sp <- terra::project(pa_sp, crs)
  
  # clean spatial outliers
  if (clean_outliers) {
    pa_sp <- pa_clean_spatial_outliers(pa_sp)  
  }
  
  fname <- paste0("purpleair_processed_", strftime(dt, "%Y-%m-%d"), ".RDS")
  fname <- fs::path_join(c(output_path, fname))
  terra::saveRDS(pa_sp, fname)
  pa_sp
  
}

# Get the current data for all sensors within a bounding box. Schedule this to build a db
#' Title
#'
#' @param nwlng 
#' @param nwlat 
#' @param selng 
#' @param selat 
#' @param max_age 
#' @param key 
#' @param dbdir 
#'
#' @returns
#' @export
#'
#' @examples
pa_acquire_realtime <- function(nwlng = -124.41, nwlat = 42.01, selng = -114.13, 
                                selat = 32.53, max_age = 3600, key = NULL,
                                dbdir = NULL) {

  if (is.null(key)) {
    key <- Sys.getenv("PURPLEAIR_READ_KEY")
  }
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = "purpleair.duckdb")
  
  url <- "https://api.purpleair.com/v1/sensors/"
  
  query_string <- list(
    api_key = key,
    max_age = max_age,
    fields = "last_seen,confidence,pm2.5",
    nwlat = nwlat,
    nwlng = nwlng,
    selat = selat,
    selng = selng
  )
  
  tryCatch({
    response <- httr::VERB("GET", url, query = query_string,
                           httr::content_type("application/octet-stream"),
                           httr::accept("application/json"))
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  },
  error = function(e) {
    message("An error occurred:\n", e)
  }
  )
  
  df <- as.data.frame(result$data)
  names(df) <- result$fields
  
  df <- df |>
    dplyr::mutate(last_seen = as.POSIXct(last_seen, tz = "UTC",
                                         format = "%FT%TZ"))
  
  recs <- DBI::dbAppendTable(con, "pa_realtime", df)
  message(paste(recs, "records added:", Sys.time()))
  DBI::dbDisconnect(con)
  
}


#' Get all of the PurpleAir sensor indices within a bounding box. See
#' https://api.purpleair.com/#api-sensors-get-sensors-data. Includes indoor and outdoor
#'
#' @param nwlng A northwest longitude for the bounding box
#' @param nwlat A northwest latitude for the bounding box
#' @param selng A southeast longitude for the bounding box
#' @param selat A southeast latitude for the bounding box
#' @param key A PurpleAir READ API key
#'
#' @return A data frame with sensor_index, date_created, last_seen, name,
#'   latitude, and longitude for all sensors within the specified bounding box.
#' @export
#'
#' @examples # CA Bounding Box
#' nwlng <- -124.41
#' nwlat <- 42.01
#' selng <- -114.13
#' selat <- 32.53
#' pa_find_sensors(nwlng, nwlat, selng, selat, key = my_api_key)
pa_find_sensors <- function(nwlng, nwlat, selng, selat, key) {

  query_string <- list(
    api_key = key,
    fields = "name,latitude,longitude,location_type,last_seen,date_created",
    #location_type = "0", # outside
    max_age = "0",
    nwlng = nwlng,
    nwlat = nwlat,
    selng = selng,
    selat = selat
  )

  url <- "https://api.purpleair.com/v1/sensors"

  response <- httr::VERB("GET", url, query = query_string,
                         httr::content_type("application/octet-stream"),
                         httr::accept("application/json"))
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

  if (!is.null(result$error)) {
    stop(result$error, ": ", result$description)
  }

  df <- as.data.frame(result$data)
  names(df) <- result$fields
  df %>%
    mutate(last_seen = lubridate::as_datetime(as.numeric(last_seen)),
           date_created = lubridate::as_datetime(as.numeric(date_created)))

}

#' Get the historical data for a specific PurpleAir sensor for a single local day. This
#' requires a PurpleAir API key with access to historical data, which can be requested.
#' For multiple days, use pa_sensor_history. This function specifies the exact hours
#' needed in local time to avoid additional request costs.
#'
#' @param sensor_index The numeric PurpleAir sensor index for for the sensor to request
#' @param key A PurpleAir API read key with access to historical data
#' @param date The date to retrieve
#' @param timezone The sensor local time zone (default is "America/Los_Angeles")
#'
#' @return A data frame for a single sensor
#' @export
#'
#' @examples
pa_sensor_history_daily <- function(sensor_index, date, timezone = "America/Los_Angeles", key) {
  
  # Convert to UTC datetime
  start_date <- as.POSIXct(date, tz = timezone)
  start_date <- lubridate::with_tz(start_date, tzone = "UTC")
  end_date <- start_date + lubridate::days(1)

  query_string <- list(
    api_key = key,
    start_timestamp = strftime(start_date, format = "%FT%TZ"),
    end_timestamp = strftime(end_date, format = "%FT%TZ"),
    average = 60,
    fields = "pm2.5_atm"
  )

  url <- paste0("https://api.purpleair.com/v1/sensors/",
                sensor_index,
                "/history")
    
  # The purpleair historical api rate limit is one per second and you can only get one
  # sensor at a time (!)
  Sys.sleep(1)
    
  response <- httr::VERB("GET", url, query = query_string,
                         httr::content_type("application/octet-stream"),
                         httr::accept("application/json"))
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

  if (!is.null(result$error)) {
    warning(result$error, ": ", result$description)
    return(NULL)
  }
  
  df <- as.data.frame(result$data)
  if (nrow(df) == 0) {
    return(NULL)
  }
  names(df) <- result$fields
  df %>%
    mutate(sensor_index = sensor_index,
           pm2.5_atm = as.numeric(pm2.5_atm),
           time_stamp = as.POSIXct(time_stamp, format = "%FT%TZ"))


}

# takes a df from pa_sensor_history_daily and a df from pa_find_sensors to construct a
# spatvector of purpleair data for modeling and save as RDS
pa_process_daily <- function(pa_raw, pa_sensors, timezone = "America/Los_Angeles",
                             clean_outliers = TRUE,
                             crs = "EPSG:3395",
                             output_path = "./processed_data/purpleair/") {
  
  # Calculate daily average and convert to log (need to remove zeros)
  pa_daily <- pa_raw |>
    mutate(time_local = lubridate::with_tz(time_stamp, tzone = timezone),
           day_local = lubridate::floor_date(time_local, unit = "days")) |>
    summarise(PM25_daily = mean(pm2.5_atm, na.rm = TRUE),
              n = n(),
              .by = c(sensor_index, day_local)) |>
    filter(PM25_daily > 0,
           n >= 15) |>
    mutate(PM25_log = log(PM25_daily))

  dates <- unique(pa_daily$day_local)
  if (length(dates) != 1) {
    stop("pa_process_daily currently only works on single day data")
  }
  
  # convert to spatial
  pa_df <- pa_daily |>
    left_join(pa_sensors, by = "sensor_index") |>
    mutate(latitude = as.numeric(latitude),
           longitude = as.numeric(longitude)) |>
    select(sensor_index, latitude, longitude, day_local, PM25_daily, PM25_log)
  
  pa_sp <- terra::vect(pa_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # project to desired coordinates
  pa_sp <- terra::project(pa_sp, crs)
  
  # clean spatial outliers
  if (clean_outliers) {
    pa_sp <- pa_clean_spatial_outliers(pa_sp)  
  }
  
  fname <- paste0("purpleair_processed_", strftime(dates, "%Y-%m-%d"), ".RDS")
  fname <- fs::path_join(c(output_path, fname))
  terra::saveRDS(pa_sp, fname)
  pa_sp
  

}

# like purpleair_clean_spatial_outliers, but using the newer spatial libraries and updated
# field names. Some updates on the algorithm. Not convert to log before applying tests.
# Increased search distance from 10km to 20km. New test is z-score. Also fail any
# locations where their are no neighbors within the radius and value of PM25_log > 7
# (~1000 ug/m3)

# New version here assumes that there is only one day in the spatv
pa_clean_spatial_outliers <- function(spatv) {

  region_stats2 <- function(ind, var) {
    mean <- mean(var[ind])
    std <- sd(var[ind])
    n <- length(ind)
    data.frame(mean, std, n)
  }
  
  # Need to convert to sf to use the operator we need
  sfs <- sf::st_as_sf(spatv)
  
  within <- sf::st_is_within_distance(sfs, dist = 20000) # trying double the search distance
  stats <- purrr::map_dfr(within, region_stats2, sfs$PM25_log)
  sfs <- bind_cols(sfs, stats) %>%
    tidyr::replace_na(list(std = 0)) %>%
    mutate(ZScore = (PM25_log - mean) / std,
           Outlier = if_else(abs(ZScore) >= 1.5, 1, 0),
           Outlier = if_else(n == 1 & PM25_log > 7, 1, Outlier))
  good_devices <- filter(sfs, Outlier == 0)
    
  # convert back to terra
  spat <- terra::vect(good_devices)
  
}





#' Check that your PurpleAir API key is working and that it is the correct type
#'
#' @param key A PurpleAir API key
#'
#' @return The type of API key. Either "READ" or "WRITE"
#' @export
#'
#' @examples
pa_check_api_key <- function(key) {
  
  query_string <- list(
    api_key = key
  )
  
  url <- "https://api.purpleair.com/v1/keys"
  
  response <- httr::VERB("GET", url, query = query_string,
                         httr::content_type("application/octet-stream"),
                         httr::accept("application/json"))
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  if (!is.null(result$error)) {
    stop(result$error, ": ", result$description)
  }
  result$api_key_type
  
}





