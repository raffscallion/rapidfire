# Daily diagnostic plots

#' Map daily predictions with monitor differences
#'
#' Plots the predicted PM2.5 surface for a day with monitor locations overlaid,
#' colored by the difference between the leave-one-out random forest prediction
#' and the observed monitor value (prediction minus observation).
#'
#' @param dt Day to plot, as a \code{Date} or a character string in
#'   \code{YYYY-MM-DD} format (used to build file names and the plot title).
#' @param data_path Directory containing the prediction GeoTIFFs written by
#'   \code{\link{predict_grid_write}}.
#' @param validation_path Directory containing the cross-validation CSVs written by
#'   \code{\link{daily_cross_validate_write}}.
#' @param monitor_path Directory containing the combined monitor RDS files
#'   (\code{monitors_combined_YYYY-MM-DD.RDS}).
#' @param state Optional \code{SpatVector} used to crop and mask the prediction
#'   raster and draw an outline. Default is \code{NULL}.
#'
#' @returns A \code{ggplot} object.
#' @export
#'
#' @seealso \code{\link{predict_grid_write}}, \code{\link{daily_cross_validate_write}}
#'
#' @examples
#' \dontrun{
#' daily_map_diagnostic(
#'   "2024-11-15",
#'   data_path = "./output/",
#'   validation_path = "./output/cross_validation/",
#'   monitor_path = "./processed_data/monitors/",
#'   state = california_shapefile
#' )
#' }
daily_map_diagnostic <- function(dt, data_path, validation_path, monitor_path, state = NULL) {
  
  rf_file <- fs::path_join(c(data_path, paste0("rapidfire_", dt, ".tif")))
  rf <- terra::rast(rf_file)
  
  state <- terra::project(state, terra::crs(rf))
  
  if (!is.null(state)) {
    rf <- terra::crop(rf, state, mask = TRUE)  
  }
  
  val_file <- fs::path_join(c(validation_path, paste0("rapidfire_l-o-o_", dt, ".csv")))
  val <- readr::read_csv(val_file)
  
  mon_file <- fs::path_join(c(monitor_path, paste0("monitors_combined_", dt, ".RDS")))
  mon <- terra::readRDS(mon_file)
  monsf <- sf::st_as_sf(mon) |>
    select(monitorID)
  
  val <- monsf |>
    left_join(val, by ="monitorID") |>
    mutate(PM25 = exp(PM25_log),
           PM25_RF = exp(PM25_log_RF),
           Diff = PM25_RF - PM25)

  ints <- terra::vect(val)
  ints <- terra::extract(rf, ints)
  
  val$PM25_RF_full <- ints$PM25_RF
  val <- filter(val, !is.na(PM25_RF_full))

  range <- (max(val$Diff, na.rm = TRUE) - min(val$Diff, na.rm = TRUE))
  zero_pt <- -min(val$Diff, na.rm = TRUE) / range
  
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = rf, ggplot2::aes(fill = PM25_RF)) +
    tidyterra::geom_spatvector(data = state, color = "gray30", fill = NA) +
    ggplot2::geom_sf(data = val, ggplot2::aes(color = Diff), size = 1) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 50), oob = scales::squish, na.value = NA) +
    colorspace::scale_color_continuous_divergingx(palette = "RdBu", mid = zero_pt,
                                                  limits = c(-20, 20), oob = scales::squish) +
    ggplot2::theme_void() +
    ggplot2::labs(fill = expression(PM[2.5]~~mu*g/m^3),
         color = "Difference\nfrom monitor",
         title = dt) +
    ggplot2::theme(legend.position = "inside",
          legend.position.inside = c(0.1, 0.3))
  
}

daily_map <- function(dt, data_path, state = NULL, borders = NULL, aggregate = NULL,
                      limits = NULL) {
  
  rf_file <- fs::path_join(c(data_path, paste0("rapidfire_", dt, ".tif")))
  rf <- terra::rast(rf_file)

  if (!is.null(state)) {
    state <- terra::project(state, terra::crs(rf))
    rf <- terra::crop(rf, state, mask = TRUE)  
  }
  
  g <- ggplot2::ggplot()

  if (!is.null(aggregate)) {
    agg <- terra::project(aggregate, terra::crs(rf))
    rf_agg <- terra::extract(rf$PM25_RF, agg, fun = mean, na.rm = TRUE)
    agg$PM25_RF <- rf_agg$PM25_RF
    agg <- agg[is.finite(agg$PM25_RF),]
    g <- g +
      tidyterra::geom_spatvector(data = agg, ggplot2::aes(fill = PM25_RF),
                                 color = NA)
    
  } else {
    g <- g +
      tidyterra::geom_spatraster(data = rf, ggplot2::aes(fill = PM25_RF))
  }

  if (!is.null(borders)) {
    borders <- terra::project(borders, terra::crs(rf))  
    g <- g +
      tidyterra::geom_spatvector(data = borders, color = "gray30", fill = NA)  
  }
  
  g <- g +
    ggplot2::scale_fill_viridis_c(limits = limits, oob = scales::squish, na.value = NA) +
    ggplot2::theme_void() +
    ggplot2::labs(fill = expression(PM[2.5]~~mu*g/m^3),
                  title = dt) +
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(0.1, 0.3))
  
  
}
