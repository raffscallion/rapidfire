
# Functions for handling monitor data, updated to use terra rather than sp.
# Spatial functions here also work on purpleair sensor data.

#' Combine AirNow and temporary monitor data for a single day
#'
#' Reads preprocessed AirNow and temporary monitor (AIRSIS/WRCC) RDS files for a
#' given date, combines them into a single \code{SpatVector}, optionally clips to a
#' spatial extent, and removes any records with missing \code{PM25_log} values.
#' A \code{source} column is added to identify the origin of each record.
#'
#' @param date A \code{Date} or date-coercible character string specifying the day
#'   to combine.
#' @param input_path_airnow Directory containing preprocessed AirNow RDS files
#'   named \code{airnow_pm25_24hr_YYYY-MM-DD.RDS}.
#' @param input_path_temp_monitors Directory containing preprocessed temporary
#'   monitor RDS files named \code{temp_monitors_pm25_24hr_YYYY-MM-DD.RDS}.
#' @param output_path Directory where the combined RDS file will be saved.
#'   Default is \code{"./processed_data/monitors/"}.
#' @param extent An optional \code{SpatVector} or \code{SpatExtent} used to
#'   spatially clip monitors. If the CRS differs from the monitor data it is
#'   reprojected automatically. If \code{NULL}, all monitors are retained.
#'
#' @returns A \code{SpatVector} of combined monitor observations with a
#'   \code{source} column (\code{"airnow"} or \code{"airsis_wrcc"}), with all
#'   \code{NA} \code{PM25_log} records removed. Also written to \code{output_path}
#'   as \code{monitors_combined_YYYY-MM-DD.RDS}.
#' @export
#'
#' @examples
#' \dontrun{
#' monitors_combine(
#'   "2024-11-15",
#'   input_path_airnow = "./processed_data/airnow/",
#'   input_path_temp_monitors = "./processed_data/temp_monitors/",
#'   extent = grid
#' )
#' }
monitors_combine <- function(date, input_path_airnow, input_path_temp_monitors, 
                             output_path = "./processed_data/monitors/",
                             extent = NULL) {

  dt <- as.Date(date)
  airnow_filename <- paste0("airnow_pm25_24hr_", strftime(dt, "%Y-%m-%d"), ".RDS")
  an <- terra::readRDS(fs::path_join(c(input_path_airnow, airnow_filename)))
  an$source <- "airnow"
  
  tempm_filename <- paste0("temp_monitors_pm25_24hr_", strftime(dt, "%Y-%m-%d"), ".RDS")
  tempm <- terra::readRDS(fs::path_join(c(input_path_temp_monitors, tempm_filename)))
  if (!is.null(tempm)) {
    tempm$source <- "airsis_wrcc"  
  }
  
  mon <- rbind(an, tempm)
  
  # If wanted, limit to monitors within the provided extent, which is a SpatVector or SpatExtent
  if (!is.null(extent)) {
    if (!identical(terra::crs(mon), terra::crs(extent))) {
      extent <- terra::project(extent, mon)      
    }
    mon <- mon[extent]
  }

  # Further modeling is not tolerant of missing data, so remove any NAs
  mon <- mon[!is.na(mon$PM25_log),]
  
  
  fname <- paste0("monitors_combined_", strftime(dt, "%Y-%m-%d"), ".RDS")
  fname <- fs::path_join(c(output_path, fname))
  terra::saveRDS(mon, fname)
  mon
  
}


#' Fit a variogram model to monitor PM2.5 data
#'
#' Computes an empirical variogram from log-transformed PM2.5 values in a
#' \code{SpatVector} of monitors and fits a Gaussian model with a nugget of 0.02.
#' The fitted model is intended for use with \code{\link{monitors_krige_grid}} or
#' \code{\link{monitors_krige_points}}.
#'
#' @param mon A \code{SpatVector} of monitor data with a \code{PM25_log} column
#'   (e.g., from \code{\link{monitors_combine}}).
#' @param cutoff Maximum distance (in CRS map units) for variogram estimation. If
#'   \code{NULL}, \code{gstat} defaults to half the diagonal of the bounding box.
#' @param width Bin width (lag size) for the empirical variogram, in CRS map units.
#'   Default is \code{15000} (15 km for a metric CRS).
#'
#' @returns A \code{variogramModel} object fitted with a Gaussian model and nugget
#'   of 0.02.
#' @export
#'
#' @examples
#' \dontrun{
#' mon <- monitors_combine("2024-11-15", input_path_airnow = "./processed_data/airnow/",
#'                         input_path_temp_monitors = "./processed_data/temp_monitors/")
#' vgm <- monitors_variogram(mon)
#' }
monitors_variogram <- function(mon, cutoff = NULL, width = 15000) {
  
  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  df <- terra::as.data.frame(mon, geom = "XY")
  # make sure all values are finite
  df <- dplyr::filter(df, is.finite(PM25_log))
  vgm_mod <- gstat::vgm(model = "Gau", nugget = 0.02)
  vgm_data <- gstat::variogram(PM25_log ~  1, locations = ~x+y, data = df,
                               cutoff = cutoff, width = width)
  v <- gstat::fit.variogram(vgm_data, vgm_mod)
  
}

#' Fit a pooled variogram model across multiple days of monitor data
#'
#' Computes a pooled empirical variogram from log-transformed PM2.5 values across
#' multiple days. The best-fitting model is selected automatically from Exponential,
#' Spherical, Gaussian, and Matérn candidates, with a nugget of 0.02.
#'
#' @param mon A \code{SpatVector} of multi-day monitor data with \code{PM25_log}
#'   and \code{Day} columns (e.g., combined output of multiple
#'   \code{\link{monitors_combine}} calls).
#'
#' @returns A \code{variogramModel} object with the best-fitting model type
#'   selected from Exponential, Spherical, Gaussian, and Matérn.
#' @export
#'
#' @seealso \code{\link{monitors_variogram}}
#'
#' @examples
#' \dontrun{
#' vgm <- monitors_variogram_pooled(mon_multiday)
#' }
monitors_variogram_pooled <- function(mon) {
  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  df <- terra::as.data.frame(mon, geom = "XY")
  # make sure all values are finite
  df <- dplyr::filter(df, is.finite(PM25_log))
  # Let the data pick the best model
  vgm_mod <- gstat::vgm(model =  c("Exp", "Sph", "Gau", "Mat"), nugget = 0.02)
  vgm_data <- gstat::variogram(PM25_log ~  Day, locations = ~x+y, data = df,
                                dX = 0)
  v <- gstat::fit.variogram(vgm_data, vgm_mod)

}

#' Krige PM2.5 values at point locations
#'
#' Performs ordinary kriging using a fitted variogram model to predict log-transformed
#' PM2.5 values at arbitrary point locations. Duplicate geometries are resolved by
#' taking the median \code{PM25_log} value, and non-finite values are removed before
#' kriging.
#'
#' @param mon A \code{SpatVector} of monitor observations with a \code{PM25_log}
#'   column used as kriging input.
#' @param locs Prediction locations as a \code{SpatVector} or \code{sf} object.
#' @param vgm A \code{variogramModel} object (e.g., from
#'   \code{\link{monitors_variogram}} or \code{\link{monitors_variogram_pooled}}).
#' @param maxdist Maximum search radius (in CRS map units) for the local kriging
#'   neighborhood. Default is \code{Inf} (global kriging).
#' @param nmax Maximum number of nearest neighbors to include in kriging. Default
#'   is \code{Inf}.
#'
#' @returns An \code{sf} object at \code{locs} with columns \code{var1.pred}
#'   (kriged \code{PM25_log} predictions) and \code{var1.var} (prediction variance).
#'
#' @seealso \code{\link{monitors_krige_grid}}, \code{\link{monitors_variogram}}
monitors_krige_points <- function(mon, locs, vgm, maxdist = Inf, nmax = Inf) {
  pts_sf <- sf::st_as_sf(mon)
  if (class(locs) == "SpatVector") {
    locs <- sf::st_as_sf(locs)
  }
  
  # Make sure all values are finite and remove duplicates to avoid covariance matrix errors
  pts_sf <- filter(pts_sf, is.finite(PM25_log))
  pts_sf <- pts_sf |>
    summarise(PM25_log = median(PM25_log, na.rm = TRUE),
              .by = geometry)
  
  ok <- gstat::krige(PM25_log ~ 1, locations = pts_sf, newdata = locs,
                     model = vgm, maxdist = maxdist, nmax = nmax)
  
}

#' Krige PM2.5 values onto a prediction grid
#'
#' Performs ordinary kriging using a fitted variogram model to predict log-transformed
#' PM2.5 values at every cell of a \code{SpatRaster} grid. Duplicate coordinates are
#' resolved by taking the median \code{PM25_log} value, and non-finite values are
#' removed before kriging.
#'
#' @param mon A \code{SpatVector} of monitor observations with a \code{PM25_log}
#'   column used as kriging input.
#' @param grid A \code{SpatRaster} defining the target prediction grid (CRS, extent,
#'   and resolution).
#' @param vgm A \code{variogramModel} object (e.g., from
#'   \code{\link{monitors_variogram}} or \code{\link{monitors_variogram_pooled}}).
#' @param maxdist Maximum search radius (in CRS map units) for the local kriging
#'   neighborhood. Default is \code{Inf} (global kriging).
#' @param nmax Maximum number of nearest neighbors to include in kriging. Default
#'   is \code{Inf}.
#'
#' @returns A \code{SpatRaster} with kriged \code{PM25_log} predictions interpolated
#'   onto \code{grid}.
#'
#' @seealso \code{\link{monitors_krige_points}}, \code{\link{monitors_variogram}}
monitors_krige_grid <- function(mon, grid, vgm, maxdist = Inf, nmax = Inf) {
  
  df <- terra::as.data.frame(mon, geom = "XY")
  # make sure all values are finite
  df <- dplyr::filter(df, is.finite(PM25_log))
  # Need to remove duplicates to avoid singular covarience matrix errors
  df <- dplyr::summarise(df,
                         PM25_log = median(PM25_log, na.rm = TRUE),
                         .by = c(x, y))
  ok <- gstat::gstat(NULL, "PM25_log", PM25_log ~ 1, locations = ~x+y, model = vgm, 
                     maxdist = maxdist, nmax = nmax, data = df)
  res <- terra::interpolate(grid, ok)
  
}

#' Split monitor data into training and test sets
#'
#' Randomly partitions a \code{SpatVector} of monitor observations into training
#' and test sets by sampling row indices without replacement.
#'
#' @param mon A \code{SpatVector} of monitor data to split.
#' @param test_fraction Fraction of monitors to allocate to the test set.
#'   Default is \code{0.3}.
#' @param seed An optional integer random seed for reproducibility. If \code{NULL},
#'   no seed is set.
#'
#' @returns A named list with two elements:
#'   \describe{
#'     \item{train}{A \code{SpatVector} of training observations.}
#'     \item{test}{A \code{SpatVector} of held-out test observations.}
#'   }
#'
#' @seealso \code{\link{monitors_combine}}
monitors_split <- function(mon, test_fraction = 0.3, seed = NULL) {
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  test_size <- round(length(mon) * test_fraction)
  i <- seq(1, length(mon), by = 1)
  test_i <- sample(i, test_size)
  test <- mon[test_i, ]
  
  train_i <- setdiff(i, test_i)
  train <- mon[train_i, ]
  
  list(test=test, train=train)
  
}

