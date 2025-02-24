# Automation functions (scripts) for developing models over a given period using
# the standard stuff (AirNow/AirSis, PurpleAir, MAIAC AOD, BlueSky, NARR)

#' develop_model
#'
#' Specify start and end dates, and states to train a model. Note that this
#' function requires MAIAC and NARR data be downloaded already. See
#' \code{\link{maiac_download}} and \code{\link{get_narr}}.
#'
#' @param dt1 date Start date of the model input data
#' @param dt2 date End date of the model input data
#' @param states character A vector of two-letter state abbreviations
#' @param pa_cutoff numeric The distance in meters to consider PurpleAir data in
#'   kriging. The default value is 100,000 m (100 km).
#'   \code{\link{create_purpleair_variogram}}.
#' @param seed numeric A number to set the randomization seed for
#'   reproducibility
#' @param bluesky_special character Handles two special cases of BlueSky data.
#'   If "2020", processes data prior to October 10, 2020 separately from that
#'   after, as the BlueSky data format changed. If "HAQAST", uses custom
#'   BlueSky-CMAQ output created during the HAQAST campaign (see
#'   \url{https://doi.org/10.1080/10962247.2021.1891994}{O'Neill et al., 2021}).
#' @param pa_data character Path to an RDS file of pre-retrieved PurpleAir data
#'   from \code{\link{get_purpleair_daterange}}.
#'
#' @return A list with two items. model contains the model as returned by
#'   \code{\link[caret]{train}}. output is a dataframe with the input data and
#'   the modeled predictions.
#' @export
#'
#' @examples dt1 <- as.Date("2018-11-01")
#' dt2 <- as.Date("2018-11-30")
#' pa <- "./data/purpleair/purpleair_2018.RDS"
#' mod_2018_nov_<- develop_model(dt1, dt2, states = "CA", pa_data = pa)
develop_model <- function(dt1, dt2, states, pa_cutoff = 100000, seed = 1977,
                          bluesky_special = NULL, pa_data = NULL,
                          bs_path = "./data/bluesky/archive/") {

  # Get and prep AirNow and mobile monitor data
  print("AirNow, AirSIS, and WRCC data...")
  an_ws <- get_monitor_daterange(dt1, dt2, states, "airnow")
  as_ws <- get_monitor_daterange(dt1, dt2, states, "airsis")
  wr_ws <- get_monitor_daterange(dt1, dt2, states, "wrcc")
  mon <- do.call(rbind, c(an_ws, as_ws, wr_ws))

  an_vg <- create_airnow_variograms(mon)
  # trying a larger fraction of test data to make a more complete RF model
  mon_split <- split_airnow_data(mon, test_fraction = 0.3)
  an_ok <- krige_airnow(mon_split, an_vg)

  ank_complete <- krige_airnow_all(mon, mon, an_vg) %>%
    select(monitorID, Day, PM25_log_ANK, PM25_log_var) %>%
    distinct()

  # Prep Bluesky For 2020, the grid changed midway through october - will need
  # to process them separately. I added this hack to address it
  print("BlueSky data...")
  if (!is.null(bluesky_special)) {
    if (bluesky_special == "2020") {
      bluesky_stack1 <- stack_bluesky_archive(dt1, as.Date("2020-10-09"))
      bluesky_stack2 <- stack_bluesky_archive(as.Date("2020-10-10"), dt2)
      mon1 <- mon[mon@data$Day <= as.POSIXct("2020-10-09"),]
      mon2 <- mon[mon@data$Day > as.POSIXct("2020-10-09"),]
      bluesky1 <- preprocessed_bluesky_at_airnow(bluesky_stack1, mon1)
      bluesky2 <- preprocessed_bluesky_at_airnow(bluesky_stack2, mon2)
      bluesky <- bind_rows(bluesky1, bluesky2)
    } else if (bluesky_special == "HAQAST") {
      #Use the HAQAST S1 CMAQ run instead of the standard Bluesky archive
      bluesky_stack <- stack_haqast_archive(dt1, dt2)
      bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
    } else {
      stop(paste("bluesky_special:", period, "not supported"))
    }
  } else {
    bluesky_stack <- stack_bluesky_archive(dt1, dt2, path = bs_path)
    bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
  }

  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(mon)

  # Prep NARR
  print("Preparing NARR...")
  narr <- narr_at_airnow(mon)

  # Get and prep PurpleAir data For dates prior to April 5, 2019, the data must
  # be acquired by create_purpleair_archive before running this script
  print("PurpleAir data...")
  if (is.null(pa_data)) {
    pa_data <- get_purpleair_daterange(dt1, dt2, states)
  } else {
    pa_data <- readRDS(pa_data) %>%
      filter(Date >= dt1,
             Date <=dt2)
  }
  pas <- purpleair_spatial(pa_data)
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)

  # Krige purple air data for all locations in the monitor data set
  pa_ok <- krige_purpleair(pa_clean, mon, pa_vg)

  # Now need to put altogether in prep for RF modeling
  print("Combining inputs for modeling...")
  train <- mon_split$train@data %>%
    mutate(Split = "Train")
  test <- mon_split$test@data %>%
    mutate(Split = "Test")

  # These have the measured PM25 - the target variable
  model_in <- bind_rows(train, test)

  # All the inputs
  # AirNow Kriged
  ank_in <- an_ok@data %>%
    select(monitorID, Day, PM25_log_ANK = var1.pred)
  # PurpleAir Kriged
  pak_in <- pa_ok@data %>%
    select(monitorID, Day, PM25_log_PAK = purp.pred)
  # NARR
  narr_in <- narr %>%
    select(-PM25, -Hours, -PM25_log)
  # Bluesky
  bluesky_in <- bluesky %>%
    select(monitorID, Day, PM25_bluesky)
  # MAIAC AOD
  maiac_in <- maiac %>%
    select(monitorID, Day, MAIAC_AOD)

  model_in <- model_in %>%
    left_join(ank_in, by = c("monitorID", "Day")) %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(narr_in, by = c("monitorID", "Day")) %>%
    left_join(bluesky_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day"))

  # Check for missing values
  model_in <- model_in %>%
    filter(!is.na(PM25_bluesky))

  # Replace missing and non-finite values with overall median
  model_in <- model_in %>%
    mutate(across(PM25_log_PAK:MAIAC_AOD,
                  ~if_else(is.finite(.x), .x, median(.x, na.rm = TRUE))))

  # Train the model
  set.seed(seed)
  train_control <- caret::trainControl(method = "cv", number = 10)
  tune_grid <- data.frame(mtry = c(2,3,4,5))

  model_in_test <- filter(model_in, Split == "Test") %>%
    filter(!is.na(PM25_log_ANK))

  model <- caret::train(PM25_log ~ PM25_log_ANK + PM25_log_PAK + PM25_bluesky + MAIAC_AOD +
                          air.2m + uwnd.10m + vwnd.10m + rhum.2m + apcp + hpbl,
                          data = model_in_test, tuneGrid = tune_grid, do.trace = 100,
                          ntree = 500, trControl = train_control, method = "rf",
                          importance = TRUE)

  print(model)

  # Now, krige all the monitor sites and make a prediction at each. Then return
  # both the model and the data frame including the inputs and prediction.
  output <- bind_rows(train, test) %>%
    left_join(ank_complete, by = c("monitorID", "Day")) %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(narr_in, by = c("monitorID", "Day")) %>%
    left_join(bluesky_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day")) %>%
    filter(!is.na(PM25_bluesky)) %>%
    filter(if_all(where(is.numeric), is.finite))

  output$PM25_log_RF <- predict(model$finalModel, output)
  list(model = model, output=output)

}

# New version substituting HRRR-smoke for both NARR and BlueSky
develop_model_h <- function(dt1, dt2, states, pa_cutoff = 100000, seed = 1977,
                          pa_data = NULL, maiac_path = "./data/MAIAC/") {
  
  # Get and prep AirNow and mobile monitor data
  print("AirNow, AirSIS, and WRCC data...")
  an_ws <- get_monitor_daterange(dt1, dt2, states, "airnow")
  as_ws <- get_monitor_daterange(dt1, dt2, states, "airsis")
  wr_ws <- get_monitor_daterange(dt1, dt2, states, "wrcc")
  mon <- do.call(rbind, c(an_ws, as_ws, wr_ws))
  
  an_vg <- create_airnow_variograms(mon)
  # trying a larger fraction of test data to make a more complete RF model
  mon_split <- split_airnow_data(mon, test_fraction = 0.3)
  an_ok <- krige_airnow(mon_split, an_vg)
  
  ank_complete <- krige_airnow_all(mon, mon, an_vg) %>%
    select(monitorID, Day, PM25_log_ANK, PM25_log_var) %>%
    distinct()

  # Prep HRRR-smoke
  hrrr <- hrrr_at_sites(mon, input_path = "../rapidfire_project/data/hrrr/processed/")
  
  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(mon, maiac_path = maiac_path)
     
  # Get and prep PurpleAir data For dates prior to April 5, 2019, the data must
  # be acquired by create_purpleair_archive before running this script
  print("PurpleAir data...")
  if (is.null(pa_data)) {
    stop("PurpleAir data must be provided and will not be acquired at runtime")
  } else {
    pa_data <- readRDS(pa_data) %>%
      filter(Date >= dt1,
             Date <=dt2)
  }
  pas <- purpleair_spatial(pa_data)
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)
  
  # Krige purple air data for all locations in the monitor data set
  pa_ok <- krige_purpleair(pa_clean, mon, pa_vg)
  
  # Now need to put altogether in prep for RF modeling
  print("Combining inputs for modeling...")
  train <- mon_split$train@data %>%
    mutate(Split = "Train")
  test <- mon_split$test@data %>%
    mutate(Split = "Test")
  
  # These have the measured PM25 - the target variable
  model_in <- bind_rows(train, test)
  
  # All the inputs
  # AirNow Kriged
  ank_in <- an_ok@data %>%
    select(monitorID, Day, PM25_log_ANK = var1.pred)
  # PurpleAir Kriged
  pak_in <- pa_ok@data %>%
    select(monitorID, Day, PM25_log_PAK = purp.pred)
  # HRRR
  hrrr_in <- hrrr %>%
    select(-PM25, -Hours, -PM25_log)
  # MAIAC AOD
  maiac_in <- maiac %>%
    select(monitorID, Day, MAIAC_AOD)
  
  model_in <- model_in %>%
    left_join(ank_in, by = c("monitorID", "Day")) %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(hrrr_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day"))
  
  # # Check for missing values
  # model_in <- model_in %>%
  #   filter(!is.na(PM25_bluesky))
  
  # Replace missing and non-finite values with overall median
  model_in <- model_in %>%
    mutate(across(PM25_log_PAK:MAIAC_AOD,
                  ~if_else(is.finite(.x), .x, median(.x, na.rm = TRUE))))
  
  # Train the model
  set.seed(seed)
  train_control <- caret::trainControl(method = "cv", number = 10)
  tune_grid <- data.frame(mtry = c(2,3,4,5))
  
  model_in_test <- filter(model_in, Split == "Test") %>%
    filter(!is.na(PM25_log_ANK))
  
  model <- caret::train(PM25_log ~ PM25_log_ANK + PM25_log_PAK + MAIAC_AOD +
                          HPBL + RH + UGRD + VGRD + TMP + MASSDEN,
                        data = model_in_test, tuneGrid = tune_grid, do.trace = 100,
                        ntree = 500, trControl = train_control, method = "rf",
                        importance = TRUE)
  
  print(model)
  
  # Now, krige all the monitor sites and make a prediction at each. Then return
  # both the model and the data frame including the inputs and prediction.
  output <- bind_rows(train, test) %>%
    left_join(ank_complete, by = c("monitorID", "Day")) %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(hrrr_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day")) %>%
    filter(if_all(where(is.numeric), is.finite))
  
  output$PM25_log_RF <- predict(model$finalModel, output)
  list(model = model, output=output)
  
}

