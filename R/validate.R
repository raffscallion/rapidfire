# Validate developed models with some hand-built cross validation

# Gather all AirNow and AirSIS data for a date range and area and split into k
# folds. For each fold, use the remaining data to predict the values. Then,
# collect and summarize.
#'
#'
#' @param dt1 date Start date of the model input data
#' @param dt2 date End date of the model input data
#' @param states character A vector of two-letter state abbreviations
#' @param k numeric Number of folds to split into (default is 10)
#' @param pa_data Optional PurpleAir data in pas-like format. If not provided, a
#'   download will be attempted.
#' @param model Path to a final model object in RDS format, such as produced by
#'   \code{\link{develop_model}}.
#' @param pa_cutoff
#'
#' @return
#' @export
#'
#' @examples
validate <- function(dt1, dt2, states = "CA", k = 10, pa_data = NULL,
                     model, pa_cutoff = 100000) {

  # Extract the model
  mod <- readRDS(model)
  mod <- mod$model$finalModel
  
  # Get and prep AirNow and mobile monitor data
  print("AirNow, AirSIS, and WRCC data...")
  an_ws <- get_monitor_daterange(dt1, dt2, states, "airnow")
  as_ws <- get_monitor_daterange(dt1, dt2, states, "airsis")
  wr_ws <- get_monitor_daterange(dt1, dt2, states, "wrcc")
  mon <- do.call(rbind, c(an_ws, as_ws, wr_ws))

  an_vg <- create_airnow_variograms(mon)

  # Get and prep all of the necessary data once
  print("BlueSky data...")
  bluesky <- bluesky_archive_at_locs(mon)

  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(mon)

  # Prep NARR
  print("Preparing NARR...")
  narr <- narr_at_airnow(mon)

  # PurpleAir
  print("Collecting PurpleAir...")
  if (is.null(pa_data)) {
    pa_data <- get_purpleair_daterange(dt1, dt2, states)
  } else {
    pa_data <- readRDS(pa_data) %>%
      filter(Date >= dt1,
             Date <=dt2)
  }
  pas <- purpleair_spatial(pa_data)
  print("Cleaning PurpleAir spatial outliers...")
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)
  pa_ok <- krige_purpleair_sitedates(pa_clean, mon, pa_vg)

  # Convert to monitoring data to data frame for sampling
  df <- mon@data
  df$lon <- mon@coords[,1]
  df$lat <- mon@coords[,2]

  # Randomize the order, then split into k groups
  samples <- slice_sample(df, prop = 1) %>%
    mutate(Row = row_number(),
           Group = (Row %% k) + 1)
  groups <- 1:k

  predict_group <- function(g) {
    # For each group, run the model
    test <- filter(samples, Group == g)
    train <- setdiff(samples, test)

    sp::coordinates(train) <- ~lon+lat
    sp::proj4string(train) <- sp::proj4string(mon)

    sp::coordinates(test) <- ~lon+lat
    sp::proj4string(test) <- sp::proj4string(train)

    ank <- krige_airnow_sitedates(train, test, an_vg)

    # All the inputs
    # PurpleAir Kriged
    pak_in <- pa_ok %>%
      select(monitorID, Day, PM25_log_PAK)
    # NARR
    narr_in <- narr %>%
      select(-one_of("PM25"), -one_of("Source"), -Hours, -PM25_log)
    # Bluesky
    bluesky_in <- bluesky %>%
      select(monitorID, Day, PM25_bluesky)
    # MAIAC AOD
    maiac_in <- maiac %>%
      select(monitorID, Day, MAIAC_AOD)

    results <- ank %>%
      left_join(pak_in, by = c("monitorID", "Day")) %>%
      left_join(narr_in, by = c("monitorID", "Day")) %>%
      left_join(bluesky_in, by = c("monitorID", "Day")) %>%
      left_join(maiac_in, by = c("monitorID", "Day")) %>%
      # Replace missing and non-finite values with overall median
      mutate(across(where(is.numeric),
                    ~ifelse(is.finite(.x), .x, median(.x, na.rm = TRUE))))

    results$PM25_log_RF <- predict(mod, results)
    results

  }

  purrr::map_dfr(groups, predict_group)

}
