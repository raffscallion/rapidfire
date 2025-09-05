# Functions for dealing with HRRR-smoke data

create_hrrrsmoke_daily <- function(date, outpath, extent_vector,
                                   tz = "America/Los_Angeles",
                                   delete_files = FALSE) {
  
  # determine files needed for 24-hr average over local diurnal time
  start_time <- ISOdate(lubridate::year(date), lubridate::month(date),
                          lubridate::day(date), 0, tz = tz)
  # We pull forcast hour 1 (F001) because precip is an accumulation variable. If we take
  # time 0, it will always be zero. Because of this, we need to start an hour (3600 s) early.
  start_time <- start_time - 3600
  times <- seq(start_time, length.out = 24, by = "1 hour")
  times_utc <- lubridate::with_tz(times, tzone = "UTC")
  dt_str <- format(times_utc, "%Y%m%d")
  hr_str <- format(times_utc, "%H")
  
  # download the files if we don't already have them in the folder
  files <- purrr::map2_chr(dt_str, hr_str, \(x, y) acquire_hrrr(x, y, outpath),
                           .progress = TRUE)

  # Calculate the 24-hr average, one layer at a time for all layers but precip
  cropped_rasters <- purrr::map(files, \(x) subset_hrrr(x, extent_vector))
  layers <- c("HPBL", "RH", "UGRD", "VGRD", "TMP", "MASSDEN")
  for (lay in layers) {
    mean_raster <- calculate_layer_mean(cropped_rasters, layer = lay)
    filename <- glue::glue("{lay}_{strftime(date, '%Y%m%d')}_hrrrsmoke.tif")
    filename <- fs::path_join(c(outpath, "processed", filename))
    terra::writeRaster(mean_raster, filename, overwrite = TRUE)
  }

  # Calculate total precip as a sum of accumulations for each 1-hour forecast
  accum_raster <- calculate_layer_sum(cropped_rasters, layer = "APCP")
  filename <- glue::glue("APCP_{strftime(date, '%Y%m%d')}_hrrrsmoke.tif")
  filename <- fs::path_join(c(outpath, "processed", filename))
  terra::writeRaster(accum_raster, filename, overwrite = TRUE)

  if (delete_files) {
    # Delete the grib files
    fs::file_delete(files)
    # Also delete the idx files
    idx <- paste0(files, ".idx")
    fs::file_delete(idx)
  }
  
}

acquire_hrrr <- function(date_str, model_hr, outpath) {
  # Check for the file in the path
  filename <- glue::glue("{date_str}_hrrr.t{model_hr}z.wrfsfcf01.grib2")
  filepath <- fs::path_join(c(outpath, filename))
  if (!fs::file_exists(filepath)) {
    url <- glue::glue("https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.{date_str}/conus/hrrr.t{model_hr}z.wrfsfcf01.grib2")
    options(timeout = max(1000, getOption("timeout")))
    success <- download.file(url, filepath, mode = "wb")
    url_idx <- paste0(url, ".idx")
    outfile_idx <- paste0(filepath, ".idx")
    success <- download.file(url_idx, outfile_idx)
  } else {
    message(filename, " already acquired")
  }
  return(filepath)

}

subset_hrrr <- function(grib_file, extent_vector = NULL, 
                        layers = c("HPBL:surface", "APCP:surface", "RH:2 m above ground",
                                   "UGRD:10 m above ground", "VGRD:10 m above ground",
                                   "TMP:2 m above ground", "MASSDEN:8 m above ground")) {

  r <- terra::rast(grib_file)

  # grab the coordinates
  coords <- terra::crs(r)

  # get the variables we want
  inds <- purrr::map_int(layers, \(x) which(stringr::str_starts(r@ptr$names, x)))
  r <- r[[inds]]

  # crop to the extent vector supplied
  if (!is.null(extent_vector)) {
    # make sure the vector is in the same coordinates
    extent <- terra::project(extent_vector, coords)
    r <- terra::crop(r, extent, mask = TRUE)
  }

  return(r)

}

calculate_layer_mean <- function(ras_list, layer) {
  
  # Select the named layer by index and reduce to the mean
  indices <- purrr::map(ras_list, \(x) which(stringr::str_starts(x@ptr$names, layer)))
  layers <- purrr::map2(ras_list, indices, \(x, y) x[[y]]) |>
    purrr::reduce(terra::mean)

}

calculate_layer_sum <- function(ras_list, layer) {
  
  # Select the named layer by index and reduce to the mean
  indices <- purrr::map(ras_list, \(x) which(stringr::str_starts(x@ptr$names, layer)))
  # This collapses to a single SpatRaster with just the layer of interest across time
  stack <- purrr::map2(ras_list, indices, \(x, y) x[[y]]) |>
    terra::rast()
  ras <- app(stack, sum)

}

#' hrrrr_at_sites
#'
#' Extract smoke PM2.5, surface temperature, winds, RF, precip, and planetary boundary
#' layer height from pre-downloaded and processed HRRRR-smoke data
#'
#' @param mon A SpatialPointsDataFrame with monitor data such as from
#'   \code{\link{recast_monitors}}
#'
#' @return The data frame from \emph{mon} with the extracted values from the HRRR-smoke
#'   data appended
#' @export
#'
#' @examples hrrrr <- hrrr_at_sites(mon)
hrrr_at_sites <- function(mon, input_path = "./data/hrrr/processed/") {
  
  date <- sort(unique(mon$Day))
  
  # for each date and each variable
  var <- c("HPBL", "RH", "UGRD", "VGRD", "TMP", "MASSDEN")
  cases <- tidyr::expand_grid(date, var)

  one_day <- function(date, var) {
    i <- mon$Day == date
    l = mon[i, ]
    vals <- extract_hrrr(date, var, l, input_path)
    colnames(vals) <- "hrrr_var"
    df <- bind_cols(l@data, vals) |>
      mutate(hrrr_name = var)
    df
  }
  res <- purrr::pmap(cases, one_day) |>
    purrr::list_rbind() |>
    tidyr::pivot_wider(names_from = hrrr_name, values_from = hrrr_var)
  
}


extract_hrrr <- function(date, var, locs, input_path = "./data/hrrr/processed/") {

  filename <- paste(var, strftime(date, "%Y%m%d"), "hrrrsmoke.tif", sep = "_")
  filename <- fs::path_join(c(input_path, filename))
  r <- terra::rast(filename)

  # Pull lon/lat out of the SpatialPointsDF in the same coords as the raster
  proj <- terra::crs(r)
  ltrans <- sp::spTransform(locs, proj)
  ll <- sp::coordinates(ltrans)
  vals <- terra::extract(x = terra::subset(r, 1), ll)
  colnames(vals) <- var
  vals
  
}