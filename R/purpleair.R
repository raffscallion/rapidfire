
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

#' Get all of the PurpleAir sensor indices within a bounding box. See
#' https://api.purpleair.com/#api-sensors-get-sensors-data
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
    fields = "name,latitude,longitude,last_seen,date_created",
    location_type = "0", # outside
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

#' Get the historical data for a specific PurpleAir sensor over a time period.
#' This requires a PurpleAir API key with access to historical data, which can
#' be requested.
#'
#' @param sensor_index The numeric PurpleAir sensor index for for the sensor to request
#' @param start_date Date the first date to request
#' @param end_date Date the last date to request
#' @param key A PurpleAir API read key with access to historical data
#'
#' @return
#' @export
#'
#' @examples
pa_sensor_history <- function(sensor_index, start_date, end_date, key) {

  # Convert to date (in case we have a full timestamp)
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date) + 1

  # Current API limits to 180 day windows for hourly averages, so break into
  # chunks if longer than that.
  start_dates <- seq.Date(from = start_date, to = end_date, by = "180 days")
  end_dates <- start_dates + 180
  end_dates[length(end_dates)] <- end_date

  pa_sensor_week <- function(start, end) {
    query_string <- list(
      api_key = key,
      start_timestamp = strftime(start, format = "%FT00:00:00Z"),
      end_timestamp = strftime(end, format = "%FT00:00:00Z"),
      average = 60,
      fields = "pm2.5_atm"
    )

    url <- paste0("https://api.purpleair.com/v1/sensors/",
                  sensor_index,
                  "/history")

    # If this doesn't work - take a 30 minute pause
    for(i in 1:4){
      if (i > 1) {
        cat("Call failed ", i, " times. Pausing 30 minutes before retry.")
        Sys.sleep(30 * 60)
      }
      try({
        response <- httr::VERB("GET", url, query = query_string,
                               httr::content_type("application/octet-stream"),
                               httr::accept("application/json"))
        result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
        break #break/exit the for-loop
      }, silent = FALSE)
    }

    if (!is.null(result$error)) {
      stop(result$error, ": ", result$description)
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

  purrr::map2_dfr(start_dates, end_dates, pa_sensor_week)

}

#' Download all of the PurpleAir sensors within a bounding polygon and a time
#' range, one sensor at a time. This requires a PurpleAir API read key with
#' access to historical data. Note that the PurpleAir API has a nominal limit of
#' 1000 queries per key per day.
#'
#' @param first_date Date First day of data to request
#' @param last_date Date Last day of data to request
#' @param bounding_poly A simple features polygon, such as a state or country
#' @param api_key A PurpleAir API read key with access to historical data
#' @param output_folder The location where the output RDS files should be stored
#' @param starting_point If specified, the rownumber to begin processing.
#'
#' @return
#' @export
#'
#' @examples
pa_batch_download_sensors <- function(first_date, last_date, bounding_poly,
                                      api_key, output_folder,
                                      starting_point = NULL) {

  # Make sure the bounding poly is in the same crs
  bounding_poly <- sf::st_transform(bounding_poly, crs = 4326)

  # Get the bounding box to query the PA API for sensors
  bbox <- sf::st_bbox(bounding_poly)

  # Query
  sensors_in_box <- pa_find_sensors(nwlng = bbox["xmin"], nwlat = bbox["ymax"],
                                    selng = bbox["xmax"], selat = bbox["ymin"],
                                    key = api_key)

  # Limit to within bounding polygon
  sens_sf <- sf::st_as_sf(sensors_in_box, coords = c("longitude", "latitude"),
                          crs = 4326)
  ints <- sf::st_intersects(sens_sf, bounding_poly, sparse = FALSE)
  pa <- sensors_in_box[ints[,1],]

  runs <- pa %>%
    select(sensor_index, start_date=date_created, end_date=last_seen) %>%
    mutate(key = api_key)

  # We save each monitor to RDS in case this fails at any time, so we can pick
  # up where we left off
  save_pa <- function(sensor_index, start_date, end_date, key) {
    df <- pa_sensor_history(sensor_index, start_date, end_date, key)
    filename <- glue::glue("{sensor_index}.RDS")
    saveRDS(df, fs::path_join(c(output_folder, filename)))
  }

  # Limit to the time period of interested (within first and last date)
  first <- as.POSIXct(first_date)
  last <- as.POSIXct(last_date)
  runs2 <- runs %>%
    mutate(start_date = if_else(start_date < first, first, start_date),
           end_date = if_else(end_date > last, last, end_date)) %>%
    filter(start_date < end_date)

  # This may fail after a while do to API limits or other issues. If so, we can
  # continue from the specified starting point
  if (!is.null(starting_point)) {
    runs2 <- slice(runs2, starting_point:nrow(runs2))
  }
  purrr::pwalk(runs2, save_pa, .progress = TRUE)
}


#' Compile downloaded individual sensor information from the PurpleAir API into
#' a single data set covering a specified time period
#'
#' @param start_date Date The first date of the time period
#' @param end_date Date The last date of the time period
#' @param sensors A data frame from \code{\link{pa_find_sensors}}
#' @param sensor_folder The folder location of the sensor data in individual
#'   files acquired from \code{\link{pa_batch_download_sensors}}
#'
#' @return A data frame of PurpleAir data that can be passed to other routines,
#'   specifically \code{\link{purpleair_spatial}}
#' @export
#'
#' @examples
pa_build_dataset <- function(start_date, end_date, bounding_poly, api_key,
                             sensor_folder, sensors = NULL) {

  if (is.null(sensors)) {

    # Make sure the bounding poly is in the same crs
    bounding_poly <- sf::st_transform(bounding_poly, crs = 4326)

    # Get the bounding box to query the PA API for sensors
    bbox <- sf::st_bbox(bounding_poly)

    # Query
    sensors_in_box <- pa_find_sensors(nwlng = bbox["xmin"], nwlat = bbox["ymax"],
                                      selng = bbox["xmax"], selat = bbox["ymin"],
                                      key = api_key)

    # Limit to within bounding polygon
    sens_sf <- sf::st_as_sf(sensors_in_box, coords = c("longitude", "latitude"),
                            crs = 4326)
    ints <- sf::st_intersects(sens_sf, bounding_poly, sparse = FALSE)
    sensors <- sensors_in_box[ints[,1],]

  }

  # Eliminate any sensors that do not have data within the time range
  df <- sensors %>%
    filter(date_created < end_date,
           last_seen > start_date)

  # For each sensor, open the file, get the data
  process_sensor <- function(i) {

    sensor_info <- sensors %>%
      filter(sensor_index == i)

    filename <- fs::path_join(c(sensor_folder, glue::glue("{i}.RDS")))
    try_read <- purrr::possibly(readRDS)
    df <- try_read(filename)
    if (is.null(df)) {
      return(NULL)
    }

    # This is in local time zone, which works for me in CA
    # Calculated daily mean, require at least 15 hours of data
    if (nrow(df) == 0) {
      return(NULL)
    }
    df <- df %>%
      mutate(Date = lubridate::floor_date(time_stamp, "days")) %>%
      filter(Date >= start_date,
             Date <= end_date) %>%
      group_by(Date) %>%
      summarise(pm25_1day = mean(pm2.5_atm, na.rm = TRUE),
                n = n()) %>%
      filter(n >= 15) %>%
      mutate(latitude = as.numeric(sensor_info$latitude),
             longitude = as.numeric(sensor_info$longitude),
             deviceID = i) %>%
      select(deviceID, Date, latitude, longitude, pm25_1day, n)
  }

  purrr::map_dfr(sensors$sensor_index, process_sensor, .progress = TRUE)

}




#' purpleair_spatial
#'
#' Get log of daily mean, select only one sensor per location, convert to
#' SpatialPointsDataFrame, and project to planar coordinates
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
purpleair_spatial <- function(df) {

  df <- df %>%
    filter(pm25_1day > 0) %>%
    mutate(Purple_log = log(pm25_1day)) %>%
    select(Date, deviceID, latitude, longitude, pm25_1day, Purple_log) %>%
    group_by(Date, latitude, longitude) %>%
    summarise(deviceID = first(deviceID),
              pm25_1day = first(pm25_1day),
              Purple_log = first(Purple_log),
              .groups = "drop")

  sp::coordinates(df) <- ~longitude+latitude
  sp::proj4string(df) <- sp::CRS("+init=epsg:4326")
  sp::spTransform(df, sp::CRS("+init=epsg:3395"))

}

#' purpleair_clean_spatial_outliers
#'
#' Scrub spatial outliers by removing and sensor that is >2 sd away from the
#' local median within a 10-km radius
#'
#' @param spdf A SpatialPointsDataFrame of PurpleAir data such as from
#'   \code{\link{purpleair_spatial}}
#'
#' @return The data from \emph{spdf} with spatial outliers removed
#' @export
#'
#' @examples     pa_data <- get_purpleair_daterange(dt1, dt2, states)
#' pas <- purpleair_spatial(pa_data)
#' pa_clean <- purpleair_clean_spatial_outliers(pas)
#'
purpleair_clean_spatial_outliers <- function(spdf) {


  region_stats <- function(ind, var) {
    med <- median(var[ind])
    std <- sd(var[ind])
    n <- length(ind)
    data.frame(med, std, n)
  }

  # Need to convert to sf to use the operator we need
  sfs <- sf::st_as_sf(spdf)

  # Could not get grouping to work, so do this one date at a time manually
  datelist <- unique(sfs$Date)

  calc_daily_stats <- function(dt, sfs, spdf) {
    sfs <- filter(sfs, Date == dt)
    within <- sf::st_is_within_distance(sfs, dist = 10000)
    stats <- purrr::map_dfr(within, region_stats, sfs$pm25_1day)
    sfs <- bind_cols(sfs, stats) %>%
      tidyr::replace_na(list(std = 0)) %>%
      mutate(AbsDiff = abs(pm25_1day - med),
             Outlier = if_else(AbsDiff > 2 * std, 1, 0))

    good_devices <- as.data.frame(sfs) %>%
      filter(Outlier == 0) %>%
      .$deviceID
    clean <- spdf[spdf$deviceID %in% good_devices,]
    clean <- clean[clean$Date == dt,]

  }

  clean_spdf <- purrr::map(datelist, calc_daily_stats, sfs, spdf)
  merged <- do.call(rbind, clean_spdf)
}


# Create binned variograms using cleaned dataset, one day at a time
# pa_vgms <- create_purpleair_variograms(pa_clean, cutoff = 100000)


#' create_purpleair_variograms
#'
#' Create binned variograms using a cleaned PurpleAir data set, one day at a
#' time
#'
#' @param df A cleaned SpatialPointsDataFrame of PurpleAir data, such as from
#'   \code{\link{purpleair_clean_spatial_outliers}}
#' @param cutoff numeric A cutoff distance in meters that is passed to
#'   \code{\link[gstat]{vgm}}
#'
#' @return A list containing daily variograms and model fits
#' @export
#'
#' @examples pa_vg <- create_purpleair_variograms(pa_clean, cutoff = 100000)
create_purpleair_variograms <- function(df, cutoff = NULL) {

  # If fit does not converge, just pass null
  try_fit <- purrr::possibly(gstat::fit.variogram, otherwise = NULL)

  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  vgm_mod <- gstat::vgm(model = "Gau", nugget = 0.02)

  daily_variogram <- function(day, cutoff) {
    df <- df[df$Date == day,]
    vgm_data <- gstat::variogram(Purple_log ~  1, data = df, cutoff = cutoff, width = 15000) %>%
      mutate(Date = day)
    mod_fit <- try_fit(vgm_data, vgm_mod)
    return(list(d=vgm_data, m=mod_fit))
  }

  dates <- unique(df$Date)
  purrr:::map(dates, daily_variogram, cutoff)

}


# Now run OK day-by-day outputting to the airnow test and training sites
# Run Ordinary Kriging on the training data at the test locations using the
# models created in create_airnow_variograms


#' krige_purpleair
#'
#' Performs daily ordinary kriging interpolation at specified locations and
#' dates
#'
#' @param pa_data Cleaned PurpleAir data such as produced by
#'   \code{\link{purpleair_clean_spatial_outliers}}
#' @param out_locs A SpatialPointsDataFrame with monitor data such as from
#'   \code{\link{recast_monitors}}. This provides the locations and dates that
#'   will be interpolated.
#' @param vgms Daily variograms as produced by
#'   \code{\link{create_purpleair_variogram}}
#'
#' @return The data from \emph{out_locs} with interpolated PurpleAir PM2.5 in
#'   log scale attached
#' @export
#'
#' @examples pa_ok <- krige_purpleair(pa_clean, mon, pa_vg)
krige_purpleair <- function(pa_data, out_locs, vgms) {

  dates <- unique(pa_data$Date)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Date)

  process_one_ok <- function(date, out_locs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Date == date) %>%
      distinct() %>%
      .$row
    test <- out_locs[out_locs$Day == date,]
    mod <- vgms[[row]]$m
    train <- pa_data[pa_data$Date == date,]
    print(date)
    ok <- gstat::krige(Purple_log ~ 1, locations = train, newdata = test,
                       model = mod)

    # Attach to measured values
    test$purp.pred <- ok$var1.pred
    test$purp.var <- ok$var1.var
    test

  }

  all <- purrr::map(dates, process_one_ok, out_locs, vgms, rows)
  # Combine together - this looks like it will be slow
  merged <- do.call(rbind, all)

}

krige_purpleair_all <- function(pa_data, out_locs, vgms) {
  dates <- unique(pa_data$Date)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Date)

  process_one_ok <- function(date, out_locs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Date == date) %>%
      distinct() %>%
      .$row
    mod <- vgms[[row]]$m
    train <- pa_data[pa_data$Date == date,]
    print(date)
    ok <- gstat::krige(Purple_log ~ 1, locations = train, newdata = out_locs,
                       model = mod)

    # Attach to measured values
    output <- out_locs@data
    output$Day <- date
    output$PM25_log_PAK <- ok$var1.pred
    output$PM25_log_PAvar <- ok$var1.var
    output

  }

  all <- purrr::map_dfr(dates, process_one_ok, out_locs, vgms, rows)

}



# Krige monitor data for all locations and dates in an input SpatialPointsDataFrame


#' krige_purpleair_sitedates
#'
#' Krige PurpleAir data for all locations and dates in a SpatialPointsDataFrame
#'
#' @param pa_data PurpleAir data as created by
#'   \code{\link{get_purpleair_daterange}} or
#'   \code{\link{create_purpleair_archive}}
#' @param outlocs A SpatialPointsDataFrame of locations (and dates) to predict
#' @param vgms Daily variograms produced by
#'   \code{\link{create_purpleair_variograms}}
#'
#' @return A dataframe with the data from \emph{outlocs} with log PM2.5
#'   (PM25_log_PAK) and variability (PM25_log_PAvar) appended
#' @export
#'
#' @examples pa_ok <- krige_purpleair_sitedates(pa_clean, locations, pa_vg)
krige_purpleair_sitedates <- function(pa_data, outlocs, vgms) {

  dates <- unique(outlocs$Day)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Date)

  process_one_ok <- function(date, pa_data, outlocs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Date == date) %>%
      distinct() %>%
      .$row
    mod <- vgms[[row]]$m
    locs <- pa_data[pa_data$Date == date,]
    outlocs <- outlocs[outlocs$Day == date,]
    print(date)
    ok <- gstat::krige(Purple_log ~ 1, locations = locs, newdata = outlocs,
                       model = mod, nmax = 100) #large maxdist was still slow - trying nmax

    # Attach to measured values
    output <- outlocs@data
    output$PM25_log_PAK <- ok$var1.pred
    output$PM25_log_PAvar <- ok$var1.var
    output

  }
  all <- purrr::map_dfr(dates, process_one_ok, pa_data, outlocs, vgms, rows)
}
