# Functions for dealing with HRRR-smoke data

create_hrrrsmoke_daily <- function(date, outpath, extent_vector,
                                   tz = "America/Los_Angeles") {
  
  
  # determine files needed for 24-hr average over local diurnal time
  start_time <- ISOdate(lubridate::year(date), lubridate::month(date),
                        lubridate::day(date), 0, tz = tz)
  times <- seq(start_time, length.out = 24, by = "1 hour")
  times_utc <- lubridate::with_tz(times, tzone = "UTC")
  dt_str <- format(times_utc, "%Y%m%d")
  hr_str <- format(times_utc, "%H")
  
  # download the files if we don't already have them in the folder
  files <- purrr::map2_chr(dt_str, hr_str, \(x, y) acquire_hrrr(x, y, outpath),
                           .progress = TRUE)
  
  # Calculate the 24-hr average, one layer at a time
  cropped_rasters <- purrr::map(files, \(x) subset_hrrr(x, extent_vector))
  layers <- c("HPBL", "APCP", "RH", "UGRD", "VGRD", "TMP", "MASSDEN")
  for (lay in layers) {
    mean_raster <- calculate_layer_mean(cropped_rasters, layer = lay)
    filename <- glue::glue("{lay}_{strftime(date, '%Y%m%d')}_hrrrsmoke.tif")
    filename <- fs::path_join(c(outpath, "processed", filename))
    terra::writeRaster(mean_raster, filename)
  }
  
}

acquire_hrrr <- function(date_str, model_hr, outpath) {
  # Check for the file in the path
  filename <- glue::glue("{date_str}_hrrr.t{model_hr}z.wrfsfcf00.grib2")
  filepath <- fs::path_join(c(outpath, filename))
  if (!fs::file_exists(filepath)) {
    url <- glue::glue("https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.{date_str}/conus/hrrr.t{model_hr}z.wrfsfcf00.grib2")
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