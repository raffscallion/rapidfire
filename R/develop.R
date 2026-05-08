# New functions for developing models

# Start with daily processed data on disk

#' Develop a random forest model for PM2.5 estimation
#'
#' Trains a random forest model for estimating surface-level PM2.5 over a date range using
#' preprocessed inputs from regulatory and temporary monitors, PurpleAir sensors, HRRR
#' meteorology, and MAIAC AOD. Monitors are split into training and test sets; kriged
#' estimates from the training monitors and PurpleAir sensors are computed at the test
#' locations along with co-located HRRR and MAIAC values. A random forest is then trained
#' using 10-fold cross-validation with \code{mtry} tuning. The function returns the fitted
#' model and variogram models needed for prediction and validation.
#'
#' @param dt1 A \code{Date} or date-coercible character string giving the first day of the
#'   training period (inclusive).
#' @param dt2 A \code{Date} or date-coercible character string giving the last day of the
#'   training period (inclusive).
#' @param paths A named list of directory paths with the following elements:
#'   \describe{
#'     \item{monitors}{Directory containing combined monitor RDS files
#'       (\code{monitors_combined_YYYY-MM-DD.RDS}).}
#'     \item{purpleair}{Directory containing preprocessed PurpleAir RDS files.}
#'     \item{hrrr}{Directory containing preprocessed HRRR GeoTIFF files.}
#'     \item{maiac}{Directory containing preprocessed MAIAC GeoTIFF files.}
#'   }
#' @param test_extent An optional \code{SpatVector} or \code{SpatExtent} used to restrict
#'   monitors to a specific region before model training. If \code{NULL}, all monitors are
#'   used.
#' @param pa_nmax Maximum number of nearest PurpleAir neighbors to use when kriging sensor
#'   estimates to monitor locations. Default is \code{100}.
#' @param seed Integer random seed for monitor splitting and model training. Default is
#'   \code{1977}.
#'
#' @returns A named list with three elements:
#'   \describe{
#'     \item{model}{The final fitted \code{randomForest} model object.}
#'     \item{vgm_mon}{A pooled \code{variogramModel} fitted to regulatory monitor
#'       data, for use in \code{\link{predict_grid}} and
#'       \code{\link{daily_cross_validate}}.}
#'     \item{vgm_pa}{A pooled \code{variogramModel} fitted to PurpleAir sensor
#'       data, for use in \code{\link{predict_grid}} and
#'       \code{\link{daily_cross_validate}}.}
#'   }
#' @export
#'
#' @seealso \code{\link{predict_grid}}, \code{\link{daily_cross_validate}},
#'   \code{\link{monitors_variogram_pooled}}, \code{\link{monitors_split}}
#'
#' @examples
#' \dontrun{
#' model_list <- develop_model(
#'   dt1 = "2024-10-01",
#'   dt2 = "2024-11-30",
#'   paths = list(
#'     monitors  = "./processed_data/monitors/",
#'     purpleair = "./processed_data/purpleair/",
#'     hrrr      = "./processed_data/HRRR/",
#'     maiac     = "./processed_data/MAIAC/"
#'   )
#' )
#' saveRDS(model_list, "./models/rapidfire_model.RDS")
#' }
develop_model <- function(dt1, dt2, paths, test_extent = NULL, pa_nmax = 100, seed = 1977) {
  
  dates <- seq.Date(as.Date(dt1), as.Date(dt2), by = "1 day")
  
  
  # load monitoring data
  load_monitoring <- function(dt) {
    monitoring_file <- fs::path_join(c(paths[["monitors"]], 
                                       paste0("monitors_combined_",
                                              strftime(dt, "%Y-%m-%d"),
                                              ".RDS")))
    mon <- terra::readRDS(monitoring_file)
  }
  mons <- purrr::map(dates, load_monitoring)  
  mon <- do.call(rbind, mons)
  
  # Remove any sites outside of the test extent (e.g., California)
  if (!is.null(test_extent)) {
    crs <- terra::crs(mon)
    test_extent <- terra::project(test_extent, crs)
    mon <- terra::intersect(mon, test_extent)
  }

  # Remove any non-finite values
  mon <- mon[is.finite(mon$PM25_log),]
  
  # Try a pooled variogram using all data
  vgm <- monitors_variogram_pooled(mon)
  
  # Split to test and training
  mon_split <- monitors_split(mon, test_fraction = 0.3, seed = seed)

  # Get the daily interpolated monitor values at the test locations using the train locations
  daily_interpolation <- function(dt, train_locs, test_locs, vgm, nmax) {
    test <- test_locs[test_locs$Day == dt,]
    train <- train_locs[train_locs$Day == dt,]
    res <- monitors_krige_points(train, test, vgm, nmax = nmax)
    coords <- sf::st_coordinates(res)
    res <- res |>
      mutate(Day = dt,
             x = coords[,1],
             y = coords[,2]) |>
      sf::st_drop_geometry()
      
  }
  dates <- unique(mon_split$test$Day)
  mon_krig <- purrr::map(dates, \(x) daily_interpolation(x, mon_split$train, 
                                                         mon_split$test, vgm, nmax = 50)) |>
    purrr::list_rbind() |>
    rename(PM25_log_ANK=var1.pred) |>
    select(-var1.var)

  # Load purpleair sensor data
  load_sensors <- function(dt) {
    pa_file <- fs::path_join(c(paths[["purpleair"]],
                               paste0("purpleair_processed_",
                                      strftime(dt, "%Y-%m-%d"),
                                      ".RDS")))
    pa <- terra::readRDS(pa_file)
    pa$Day <- dt
    pa
  }
  pas <- purrr::map(dates, load_sensors)
  pa <- do.call(rbind, pas)
  
  # remove any non-finite values
  pa <- pa[is.finite(pa$PM25_log),]
  
  # Create a pooled variogram for sensors
  vgm_pa <- monitors_variogram_pooled(pa)
  
  # get daily interpolated sensor values at the test locations
  pa_krig <- purrr::map(dates, \(x) daily_interpolation(x, pa, mon_split$test, vgm_pa,
                                                        nmax = pa_nmax)) |>
    purrr::list_rbind() |>
    rename(PM25_log_PAK=var1.pred) |>
    select(-var1.var)

  # load HRRR, get values at test locations, and log-transform MASSDEN
  load_hrrr <- function(dt, test_locs) {
    hrrr <- hrrr_stack(dt, paths[["hrrr"]])
    test <- test_locs[test_locs$Day == dt,]
    coords <- terra::crds(test)
    extracted <- terra::extract(hrrr, test)
    extracted <- extracted |>
      mutate(Day = dt,
             MASSDEN_log = log(MASSDEN),
             x = coords[,1],
             y = coords[,2]) |>
      select(-ID)
    
  }
  hrrr <- purrr::map(dates, \(x) load_hrrr(x, mon_split$test)) |>
    purrr::list_rbind()

  # Load MAIAC and get values at test locations
  load_maiac <- function(dt, test_locs) {
    maiac_file <- fs::path_join(c(paths[["maiac"]],
                                  paste0("MAIAC_processed_", strftime(dt, "%Y-%m-%d"), ".tif")))
    maiac <- terra::rast(maiac_file)
    maiac <- terra::project(maiac, terra::crs(test_locs))
    test <- test_locs[test_locs$Day == dt,]
    coords <- terra::crds(test)
    extracted <- terra::extract(maiac, test)
    extracted <- extracted |>
      mutate(Day = dt,
             x = coords[,1],
             y = coords[,2]) |>
      select(-ID)
  }
  maiac <- purrr::map(dates, \(x) load_maiac(x, mon_split$test)) |>
    purrr::list_rbind()
  
  # Combine all data at test locations
  model_inputs <- mon_krig |>
    full_join(pa_krig, by = c("x", "y", "Day")) |>
    full_join(hrrr, by = c("x", "y", "Day")) |>
    full_join(maiac, by = c("x", "y", "Day"))
  
  # Join with actual monitored PM2.5
  m_coords <- terra::crds(mon_split$test)
  m <- as.data.frame(mon_split$test) |>
    select(monitorID, Day, PM25_log) |>
    mutate(x = m_coords[,1],
           y = m_coords[,2])
  
  final_input <- m |>
    full_join(model_inputs, by = c("x", "y", "Day"))
  
  # Replace any missing and non-finite values with overall median
  final_input <- final_input %>%
    mutate(across(PM25_log:MAIAC_AOD,
                  ~if_else(is.finite(.x), .x, median(.x, na.rm = TRUE))))
  
  # Train the model
  set.seed(seed)
  train_control <- caret::trainControl(method = "cv", number = 10)
  tune_grid <- data.frame(mtry = c(2,3,4,5))
  
  model <- caret::train(PM25_log ~ PM25_log_ANK + PM25_log_PAK + MAIAC_AOD +
                          HPBL + RH + UGRD + VGRD + TMP + MASSDEN_log + APCP,
                        data = final_input, tuneGrid = tune_grid, do.trace = 100,
                        ntree = 500, trControl = train_control, method = "rf",
                        importance = TRUE)
  
  print(model)
  list(model = model$finalModel, vgm_mon = vgm, vgm_pa = vgm_pa)
  
}
