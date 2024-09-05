
#' get_airnow_daterange
#'
#' For the specified date range and US states, downloads AirNow data.
#'
#' @param start Date The earliest date to search
#' @param end Date The latest date to search. Must be in the same calendar year
#'   as \code{start}.
#' @param states character The states to include as a vector of two-character
#'   state abbreviations
#'
#' @return AirNow data for the specified dates and locations as SpatialPointsDataFrame
#' @export
#'
#' @examples dt1 <- as.Date("2018-11-01")
#'   dt2 <- as.Date("2018-11-30")
#'   an_ws <- get_airnow_daterange(dt1, dt2, c("CA", "NV"))
get_airnow_daterange <- function(start, end, states) {

  end <- end + 1

  if (Sys.Date() - end < 45) {
    raw1 <- airnow_loadDaily() |>
      subset_monitors(start, end, states) |>
      recast_monitors()

  }
  if (Sys.Date() - start > 45) {
    raw2 <- airnow_loadAnnual(lubridate::year(start)) |>
      subset_monitors(start, end, states) |>
      recast_monitors()
  }

  if (exists("raw1")) {
    if (exists("raw2")) {
      # Need to merge the data sets. There may be overlap to remove.
      last_date_raw2 <- max(raw2$Day)
      raw1 <- raw1[raw1$Day > last_date_raw2,]
      return(rbind(raw1, raw2))
    } else {
      return(raw1)
    }
  } else {
    return(raw2)
  }
  
}

# Limit monitor data downloaded from AirFire to a specific date range and state list
subset_monitors <- function(ws, start, end, states) {
  ws$meta <- ws$meta |>
    dplyr::filter(stateCode %in% states)
  ws$data <- ws$data |>
    dplyr::filter(datetime >= lubridate::force_tz(start, "America/Los_Angeles"),
                  datetime <= lubridate::force_tz(end, "America/Los_Angeles"))
  return(ws)
}

#' get_airsis_daterange
#' 
#' NOTE!! As of 2024-09-05, the most recent airsis data in the AirFire archive was 2022-09-19
#'
#' For the specified date range and US states, downloads AIRSIS data.
#'
#' @param start Date The earliest date to search
#' @param end Date The latest date to search. Must be in the same calendar year
#'   as \code{start}.
#' @param states character The states to include as a vector of two-character
#'   state abbreviations
#'
#' @return AIRSIS data for the specified dates and locations as A SpatialPointsDataFrame
#' @export
#'
#' @examples dt1 <- as.Date("2018-11-01")
#'   dt2 <- as.Date("2018-11-30")
#'   an_ws <- get_airsis_daterange(dt1, dt2, c("CA", "NV"))
get_airsis_daterange <- function(start, end, states) {

  end <- end + 1

  if (Sys.Date() - end < 45) {
    raw1 <- airsis_loadDaily() |>
      subset_monitors(start, end, states) |>
      recast_monitors()
  }
  if (Sys.Date() - start > 45) {
    raw2 <- airsis_loadAnnual(lubridate::year(start)) |>
      subset_monitors(start, end, states) |>
      recast_monitors()
  }

  if (exists("raw1")) {
    if (exists("raw2")) {
      # Need to merge the data sets. There may be overlap to remove.
      last_date_raw2 <- max(raw2$Day)
      raw1 <- raw1[raw1$Day > last_date_raw2,]
      return(rbind(raw1, raw2))
      
    } else {
      return(raw1)
    }
  } else {
    return(raw2)
  }

}

#' recast_monitors
#'
#' Switch from the compact \emph{ws_monitor} format to a square tibble, get
#' daily mean, calculate Log, convert to spatial, and reproject to planar
#' coordinates
#'
#' @param mon ws_monitor object of AirNow or AIRSIS data
#'
#' @return A SpatialPointsDataFrame
#' @export
#'
#' @examples dt1 <- as.Date("2018-11-01")
#'   dt2 <- as.Date("2018-11-30")
#'   an_ws <- get_airsis_daterange(dt1, dt2, c("CA", "NV")) %>%
#'     recast_monitors()
recast_monitors <- function(mon) {
  df <- mon$data %>%
    tidyr::gather("monitorID", "PM25", -datetime)
  meta <- mon$meta %>%
    select(monitorID, longitude, latitude, stateCode, monitorType)
  df <- left_join(df, meta, by = "monitorID") %>%
    filter(!is.na(PM25)) %>%
    mutate(LocalTime = lubridate::with_tz(datetime, "America/Los_Angeles"),
           Day = lubridate::floor_date(LocalTime, "days")) %>%
    group_by(monitorID, Day) %>%
    summarise(PM25 = mean(PM25, na.rm = FALSE),
              Hours = n(),
              .groups = "drop") %>%
    filter(Hours >= 16,
           PM25 > 0) %>%
    mutate(PM25_log = log(PM25))

  sites <- mon$meta %>%
    select(monitorID, longitude, latitude) %>%
    distinct()

  # edited to inner_join from left_join when joining sites to df
  df <- df %>%
    inner_join(sites, by = "monitorID")

  # make spatial and change to planar coordinates
  sp::coordinates(df) <- ~longitude+latitude
  sp::proj4string(df) <- sp::CRS("+init=epsg:4326")
  df <- sp::spTransform(df, sp::CRS("+init=epsg:3395"))
  return(df)

}

#' create_airnow_variograms
#'
#' Calculates daily spatial variograms for kriging interpolation on monitor data
#'
#' @param df A SpatialPointsDataFrame from \code{\link{recast_monitors}}. Can
#'   include AirNow and/or AIRSIS data.
#' @param cutoff numeric A distance passed to \code{\link[gstat]{variogram}}
#'
#' @return A list of variograms, one for each date in the input data set
#' @export
#'
#' @examples an_vg <- create_airnow_variograms(an_ws)
create_airnow_variograms <- function(df, cutoff = NULL) {

  # If fit does not converge, just pass null
  try_fit <- purrr::possibly(gstat::fit.variogram, otherwise = NULL)

  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  vgm_mod <- gstat::vgm(model = "Gau", nugget = 0.02)

  daily_variogram <- function(day, cutoff) {
    df <- df[df$Day == day,]
    vgm_data <- gstat::variogram(PM25_log ~  1, data = df, cutoff = cutoff, width = 15000) %>%
      mutate(Day = day)
    mod_fit <- try_fit(vgm_data, vgm_mod)
    return(list(d=vgm_data, m=mod_fit))
  }

  dates <- unique(df$Day)
  purrr:::map(dates, daily_variogram, cutoff)

}

#' split_airnow_data
#'
#' Create a training and test data set - test dataset is a fraction of the
#' samples, chosen randomly for each day.
#'
#' @param df A SpatialPointsDataFrame from \code{\link{recast_monitors}}. Can
#'   include AirNow and/or AIRSIS data.
#' @param test_fraction The fraction of data to withold as the test data set.
#'   Default is 0.3.
#' @param seed numeric A starting randomization seed. Default is 1977.
#'
#' @return A list with two SpatialPointsDataFrames, test and train
#' @export
#'
#' @examples mon_split <- split_airnow_data(an_ws, test_fraction = 0.3)
split_airnow_data <- function(df, test_fraction = 0.3, seed = 1977) {

  if (is.numeric(seed)) {
    set.seed(seed)
  }

  n <- nrow(df@data)
  test_ind <- sample(n, size = round(n * test_fraction))
  test <- df[test_ind,]
  train <- df[-test_ind,]
  list(test=test, train=train)

}

# Run Ordinary Kriging on the training data at the test locations using the
# models created in create_airnow_variograms

#' krige_airnow
#'
#' Run ordinary kriging interpolation on the training data at the test locations
#' using the data set created by \code{\link{split_airnow_data}} and the
#' variograms created by \code{\link{create_airnow_variograms}}. This is to
#' create subset variograms that can be used in developing and validating the
#' final RF model. To interpolate at all sites, use
#' \code{\link{krige_airnow_all}}.
#'
#' @param an_data Named list from \code{\link{split_airnow_data}}.
#' @param vgms A list of daily variograms from
#'   \code{\link{create_airnow_variograms}}.
#'
#' @return A data frame of the test data with the kriged values of PM25_log from
#'   AirNow and/or AIRSIS data.
#' @export
#'
#' @examples mon_split <- split_airnow_data(an_ws, test_fraction = 0.3)
#'           an_vg <- create_airnow_variograms(an_ws)
#'           an_ok <- krige_airnow(mon_split, an_vg)
krige_airnow <- function(an_data, vgms) {

  dates <- unique(an_data$train$Day)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Day)

  process_one_ok <- function(date, an_data, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Day == date) %>%
      distinct() %>%
      .$row
    mod <- vgms[[row]]$m
    test <- an_data$test
    test <- test[test$Day == date,]
    train <- an_data$train
    train <- train[train$Day == date,]
    ok <- gstat::krige(PM25_log ~ 1, locations = train, newdata = test,
                       model = mod)

    # Attach to measured values
    test$var1.pred <- ok$var1.pred
    test$var1.var <- ok$var1.var
    test

  }

  all <- purrr::map(dates, process_one_ok, an_data, vgms, rows)
  # Combine together
  merged <- do.call(rbind, all)

}

# We now have values at the test sites that can be used in developing the RF model,
# plus the model itself for use in predicting the full surface

# This returns a data frame

#' krige_airnow_all
#'
#' Runs daily ordinary kriging based on PM25_log data from a monitor SPDF from
#' \code{\link{recast_monitors}} to a set of locations and dates.
#'
#' @param an A SpatialPointsDataFrame from \code{\link{recast_monitors}}. These
#'   are the input locations with date specific data to be interpolated.
#' @param outlocs A SpatialPointsDataFrame from \code{\link{recast_monitors}}.
#'   These are the output locations to be estimated.
#' @param vgms A list of daily variograms from
#'   \code{\link{create_airnow_variograms}}.
#'
#' @return The dataframe portion of \emph{outlocs} with the kriging prediction
#'   and variance added (PM25_log_ANK and PM25_log_var)
#' @export
#'
#' @examples ank_complete <- krige_airnow_all(an_ws, an_ws an_vg)
krige_airnow_all <- function(an, outlocs, vgms) {

  dates <- unique(an$Day)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Day) %>%
    distinct()

  process_one_ok <- function(date, an, outlocs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Day == date) %>%
      distinct() %>%
      .$row
    mod <- vgms[[row]]$m
    locs <- an[an$Day == date,]
    outlocs <- outlocs[outlocs$Day == date,]
    print(date)
    ok <- gstat::krige(PM25_log ~ 1, locations = locs, newdata = outlocs,
                       model = mod)

    # Attach to measured values
    output <- outlocs@data
    output$Day <- date
    output$PM25_log_ANK <- ok$var1.pred
    output$PM25_log_var <- ok$var1.var
    output

  }

  all <- purrr::map_dfr(dates, process_one_ok, an, outlocs, vgms, rows)
}

#' krige_airnow_sitedates
#'
#' Krige monitor data for all locations and dates in an input
#' SpatialPointsDataFrame
#'
#' @param an A SpatialPointsDataFrame as produced by
#'   \code{\link{recast_monitors}}
#' @param outlocs A SpatialPointsDataFrame of locations (and dates) to predict
#' @param vgms Daily variograms as produced by
#'   \code{\link{create_airnow_variograms}}
#'
#' @return A dataframe with the data from \emph{outlocs} with log PM2.5
#'   (PM25_log_ANK) and variability (PM25_log_var) appended
#' @export
#'
#' @examples ank <- krige_airnow_sitedates(mon, locations, an_vg)
krige_airnow_sitedates <- function(an, outlocs, vgms) {

  dates <- unique(outlocs$Day)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Day)

  process_one_ok <- function(date, an, outlocs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Day == date) %>%
      distinct() %>%
      .$row
    mod <- vgms[[row]]$m
    locs <- an[an$Day == date,]
    outlocs <- outlocs[outlocs$Day == date,]
    print(date)
    ok <- gstat::krige(PM25_log ~ 1, locations = locs, newdata = outlocs,
                       model = mod)

    # Attach to measured values
    output <- outlocs@data
    output$PM25_log_ANK <- ok$var1.pred
    output$PM25_log_var <- ok$var1.var
    output

  }
  all <- purrr::map_dfr(dates, process_one_ok, an, outlocs, vgms, rows)
}
