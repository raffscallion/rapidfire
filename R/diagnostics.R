# Daily diagnostic plots

daily_map <- function(dt, data_path, validation_path, monitor_path, state = NULL) {
  
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
    colorspace::scale_color_continuous_divergingx(palette = "RdBu", mid = zero_pt) +
    ggplot2::theme_void() +
    labs(fill = expression(PM[2.5]~~mu*g/m^3),
         color = "Difference\nfrom monitor",
         title = dt) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.1, 0.3))
  
}


