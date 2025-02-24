# Use an existing developed model to predict values at given new locations and dates

#' predict_locs
#'
#' Uses an existing model, such as developed by \code{\link{develop_model}}, to
#' predict values a given locations and dates. Prepares the data necessary to do
#' so. MAIAC, NARR, BlueSky, and (optionally) PurpleAir data should already be
#' downloaded.
#'
#' @param dt1 Date The earliest date of data to predict
#' @param dt2 Date The latest date to predict
#' @param states character A vector of two-character state codes
#' @param model A final model object, such as extracted from the results of
#'   \code{\link{develop_model}}. For example, if the object was called mod,
#'   this would be mod$model$finalModel
#' @param locations A SpatialPointsDataFrame of locations (and dates) to predict
#' @param pa_cutoff A cutoff value passed to \code{\link[gstat]{vgm}} for
#'   interploation of PurpleAir data
#' @param bluesky_special  character Handles three special cases of BlueSky
#'   data. If "2020", processes data prior to October 10, 2020 separately from
#'   that after, as the BlueSky data format changed. If "HAQAST", uses custom
#'   BlueSky-CMAQ output created during the HAQAST campaign (see
#'   \url{https://doi.org/10.1080/10962247.2021.1891994}{O'Neill et al., 2021}).
#'   If "nominal", instead use a nominal placeholder value for BlueSky PM2.5 of
#'   0.1, which will have a near neutral impact on predictions.
#' @param pa_data PurpleAir data in pas-like format.
#'
#' @return A dataframe with all model input values and the resulting predictions
#'   for the locations and dates specified in \emph{locations}
#' @export
#'
#' @examples pred <- predict_locs(dt1, dt2, states = "CA", model = mod$model$finalModel,
#'                                locations = locs, pa_data = pa)
predict_locs <- function(dt1, dt2, states = "CA", model, locations,
                         pa_cutoff = 100000, bluesky_special = NULL,
                         pa_data = NULL) {

  # Get and prep AirNow and mobile monitor data
  print("AirNow, AirSIS, and WRCC data...")
  an_ws <- get_monitor_daterange(dt1, dt2, states, "airnow")
  as_ws <- get_monitor_daterange(dt1, dt2, states, "airsis")
  wr_ws <- get_monitor_daterange(dt1, dt2, states, "wrcc")
  mon <- do.call(rbind, c(an_ws, as_ws, wr_ws))

  an_vg <- create_airnow_variograms(mon)
  ank <- krige_airnow_sitedates(mon, locations, an_vg)

  # Prep Bluesky For 2020, the grid changed midway through october - will need
  # to process them separately. I added this hack to address it
  print("BlueSky data...")
  if (!is.null(bluesky_special)) {
    if (bluesky_special == "2020") {
      bluesky <- bluesky_archive_at_locs(locations)
    } else if (bluesky_special == "HAQAST") {
      #Use the HAQAST S1 CMAQ run instead of the standard Bluesky archive
      stop("Need to implement bluesky_haqast_at_locs")
      bluesky_stack <- stack_haqast_archive(dt1, dt2)
      bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
    } else if (bluesky_special == "nominal") {
      bluesky <- bluesky_nominal(locations)
    } else {
      stop(paste("bluesky_special:", bluesky_special, "not supported"))
    }
  } else {
    bluesky <- bluesky_archive_at_locs(locations)
  }

  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(locations)

  # Prep NARR
  print("Preparing NARR...")
  narr <- narr_at_airnow(locations)

  # Get and prep PurpleAir data For dates prior to April 5, 2019, the data must
  # be acquired by create_purpleair_archive before running this script
  print("PurpleAir data...")
  pa_data <- readRDS(pa_data) %>%
    filter(Date >= dt1,
          Date <=dt2)

  pas <- purpleair_spatial(pa_data)
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)

  # Krige purple air data for all locations in the monitor data set
  pa_ok <- krige_purpleair_sitedates(pa_clean, locations, pa_vg)

  # All the inputs
  # PurpleAir Kriged
  pak_in <- pa_ok %>%
    select(monitorID, Day, PM25_log_PAK)
  # NARR
  narr_in <- narr %>%
    select(-one_of("PM25"), -one_of("Source"))
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

  # How do I require randomForest be loaded here?
  results$PM25_log_RF <- predict(model, results)
  return(results)
}


predict_locs_h <- function(dt1, dt2, states = "CA", model, locations,
                         pa_cutoff = 100000, pa_data = NULL,
                         hrrr_path = "./data/hrrr/processed/") {
  
  # Get and prep AirNow and mobile monitor data
  print("AirNow, AirSIS, and WRCC data...")
  an_ws <- get_monitor_daterange(dt1, dt2, states, "airnow")
  as_ws <- get_monitor_daterange(dt1, dt2, states, "airsis")
  wr_ws <- get_monitor_daterange(dt1, dt2, states, "wrcc")
  mon <- do.call(rbind, c(an_ws, as_ws, wr_ws))
  
  an_vg <- create_airnow_variograms(mon)
  ank <- krige_airnow_sitedates(mon, locations, an_vg)
  
  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(locations)
  
  # Prep HRRR-smoke
  print("Preparing HRRR-smoke")
  hrrr <- hrrr_at_sites(mon, input_path = hrrr_path)
  
  # the data must created by pa_build_dataset
  print("PurpleAir data...")
  pa_data <- readRDS(pa_data) %>%
    filter(Date >= dt1,
           Date <=dt2)
  
  pas <- purpleair_spatial(pa_data)
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)
  
  # Krige purple air data for all locations in the monitor data set
  pa_ok <- krige_purpleair_sitedates(pa_clean, locations, pa_vg)
  
  # All the inputs
  # PurpleAir Kriged
  pak_in <- pa_ok %>%
    select(monitorID, Day, PM25_log_PAK)
  # HRRR
  hrrr_in <- hrrr %>%
    select(-PM25, -Hours, -PM25_log)
  # MAIAC AOD
  maiac_in <- maiac %>%
    select(monitorID, Day, MAIAC_AOD)
  
  results <- ank %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(hrrr_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day")) %>%
    # Replace missing and non-finite values with overall median
    mutate(across(where(is.numeric),
                  ~ifelse(is.finite(.x), .x, median(.x, na.rm = TRUE))))
  
  # How do I require randomForest be loaded here?
  results$PM25_log_RF <- predict(model, results)
  return(results)
}