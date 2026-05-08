# Functions for dealing with HRRR-smoke data

#' Download and preprocess HRRR-smoke data for a single day
#'
#' Downloads the 24 HRRR forecast hour-1 (F001) GRIB2 files needed to represent a
#' full local calendar day, subsets and regrids them to \code{extent_grid}, then
#' computes 24-hour means for atmospheric variables (HPBL, RH, UGRD, VGRD, TMP,
#' MASSDEN) and the cumulative precipitation total (APCP). Each variable is written
#' to \code{outpath} as a separate GeoTIFF named
#' \code{VARNAME_YYYYMMDD_hrrrsmoke.tif}.
#'
#' @param date A \code{Date} or date-coercible character string specifying the local
#'   day to acquire.
#' @param outpath Directory path for saving both the raw GRIB2 downloads and the
#'   processed GeoTIFF outputs.
#' @param extent_grid A \code{SpatRaster} defining the target CRS, spatial extent,
#'   and resolution for all output files (e.g., the model prediction grid).
#' @param tz Local timezone string used to identify the 24 UTC hours that correspond
#'   to the requested local calendar day. Default is \code{"America/Los_Angeles"}.
#' @param crs Output coordinate reference system string passed to
#'   \code{terra::project}. Default is \code{"EPSG:3395"}.
#' @param delete_files Logical. If \code{TRUE}, the raw GRIB2 files and their
#'   \code{.idx} index files are deleted after processing. Default is \code{FALSE}.
#'
#' @returns Called for its side effect of writing one GeoTIFF per variable to
#'   \code{outpath}. Invisibly returns \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' hrrr_acquire("2024-11-15", outpath = "./raw_data/HRRR/", extent_grid = grid)
#' }
hrrr_acquire <- function(date, outpath, extent_grid,
                                   tz = "America/Los_Angeles",
                                   crs = "EPSG:3395",
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
  files <- purrr::map2_chr(dt_str, hr_str, \(x, y) hrrr_download(x, y, outpath),
                           .progress = TRUE)

  # Calculate the 24-hr average, one layer at a time for all layers but precip
  cropped_rasters <- purrr::map(files, \(x) hrrr_subset(x, extent_grid))
  layers <- c("HPBL", "RH", "UGRD", "VGRD", "TMP", "MASSDEN")
  for (lay in layers) {
    mean_raster <- calculate_layer_mean(cropped_rasters, layer = lay)
    mean_raster <- terra::project(mean_raster, crs)
    filename <- glue::glue("{lay}_{strftime(date, '%Y%m%d')}_hrrrsmoke.tif")
    filename <- fs::path_join(c(outpath, filename))
    terra::writeRaster(mean_raster, filename, overwrite = TRUE)
  }
  
  # Calculate total precip as a sum of accumulations for each 1-hour forecast
  accum_raster <- calculate_layer_sum(cropped_rasters, layer = "APCP")
  accum_raster <- terra::project(accum_raster, crs)
  filename <- glue::glue("APCP_{strftime(date, '%Y%m%d')}_hrrrsmoke.tif")
  filename <- fs::path_join(c(outpath, filename))
  terra::writeRaster(accum_raster, filename, overwrite = TRUE)
  
  if (delete_files) {
    # Delete the grib files
    fs::file_delete(files)
    # Also delete the idx files
    idx <- paste0(files, ".idx")
    fs::file_delete(idx)
  }
  
}

#' Download a single HRRR-smoke GRIB2 file from the NOAA S3 archive
#'
#' Downloads forecast hour F001 for a given UTC date and model run hour from the
#' NOAA HRRR public S3 bucket. If the file is already present in \code{outpath},
#' the download is skipped. The companion \code{.idx} index file is also downloaded.
#'
#' @param date_str UTC date string in \code{"YYYYMMDD"} format.
#' @param model_hr UTC model run hour string in \code{"HH"} format (e.g.,
#'   \code{"06"}).
#' @param outpath Directory path where GRIB2 files are stored.
#'
#' @returns A character string giving the local file path of the GRIB2 file
#'   (whether newly downloaded or already present).
#'
#' @seealso \code{\link{hrrr_acquire}}
hrrr_download <- function(date_str, model_hr, outpath) {
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

#' Subset and optionally regrid a HRRR GRIB2 file
#'
#' Reads a HRRR GRIB2 file, extracts the specified layers by matching name prefixes,
#' and optionally reprojects and resamples the result to match an \code{extent_grid}.
#'
#' @param grib_file Path to a HRRR GRIB2 file.
#' @param extent_grid A \code{SpatRaster} defining the target CRS, extent, and
#'   resolution. If \code{NULL}, no spatial transformation is applied.
#' @param layers A character vector of GRIB2 layer name prefixes to extract.
#'   Defaults to the standard set used in rapidfire modeling: HPBL, APCP, RH,
#'   UGRD, VGRD, TMP, and MASSDEN at their surface or near-surface HRRR levels.
#'
#' @returns A \code{SpatRaster} containing the requested layers, optionally
#'   reprojected and resampled to match \code{extent_grid}.
#'
#' @seealso \code{\link{hrrr_acquire}}
hrrr_subset <- function(grib_file, extent_grid = NULL, 
                        layers = c("HPBL:surface", "APCP:surface", "RH:2 m above ground",
                                   "UGRD:10 m above ground", "VGRD:10 m above ground",
                                   "TMP:2 m above ground", "MASSDEN:8 m above ground")) {

  r <- terra::rast(grib_file)

  # get the variables we want
  inds <- purrr::map_int(layers, \(x) which(stringr::str_starts(names(r), x)))
  r <- r[[inds]]
  
  # crop and resample to the template extent grid supplied
  if (!is.null(extent_grid)) {
    # reproject to the same coords
    coords <- terra::crs(extent_grid)
    r <- terra::project(r, coords)
    r <- terra::resample(r, extent_grid)
  }
  
  return(r)
  
}


#' Compute pixel-wise mean of a named layer across a list of SpatRasters
#'
#' Extracts a single named layer from each \code{SpatRaster} in \code{ras_list} by
#' matching the layer name prefix, then reduces the stack to a pixel-wise mean.
#'
#' @param ras_list A list of \code{SpatRaster} objects, each containing a layer
#'   whose name starts with \code{layer} (e.g., output of
#'   \code{\link{hrrr_subset}}).
#' @param layer A character string giving the layer name prefix to extract
#'   (e.g., \code{"HPBL"}, \code{"TMP"}).
#'
#' @returns A single-layer \code{SpatRaster} containing the pixel-wise mean of
#'   the named layer across all elements of \code{ras_list}.
#'
#' @seealso \code{\link{calculate_layer_sum}}, \code{\link{hrrr_acquire}}
calculate_layer_mean <- function(ras_list, layer) {
  
  # Select the named layer by index and reduce to the mean
  indices <- purrr::map(ras_list, \(x) which(stringr::str_starts(names(x), layer)))
  layers <- purrr::map2(ras_list, indices, \(x, y) x[[y]]) |>
    purrr::reduce(terra::mean)

}

#' Compute pixel-wise sum of a named layer across a list of SpatRasters
#'
#' Extracts a single named layer from each \code{SpatRaster} in \code{ras_list} by
#' matching the layer name prefix, stacks them, and reduces to a pixel-wise sum.
#' Used primarily to accumulate hourly precipitation into a daily total.
#'
#' @param ras_list A list of \code{SpatRaster} objects, each containing a layer
#'   whose name starts with \code{layer} (e.g., output of
#'   \code{\link{hrrr_subset}}).
#' @param layer A character string giving the layer name prefix to extract
#'   (e.g., \code{"APCP"}).
#'
#' @returns A single-layer \code{SpatRaster} containing the pixel-wise sum of the
#'   named layer across all elements of \code{ras_list}, with the layer renamed to
#'   \code{layer}.
#'
#' @seealso \code{\link{calculate_layer_mean}}, \code{\link{hrrr_acquire}}
calculate_layer_sum <- function(ras_list, layer) {

  # Select the named layer by index and reduce to the mean
  indices <- purrr::map(ras_list, \(x) which(stringr::str_starts(names(x), layer)))
  # This collapses to a single SpatRaster with just the layer of interest across time
  stack <- purrr::map2(ras_list, indices, \(x, y) x[[y]]) |>
    terra::rast()
  ras <- terra::app(stack, sum)
  names(ras) <- layer
  ras

}

#' Load preprocessed HRRR-smoke variables into a stacked SpatRaster
#'
#' Reads the per-variable GeoTIFF files written by \code{\link{hrrr_acquire}} for a
#' given date and stacks them into a single multi-layer \code{SpatRaster}, with each
#' layer named by its variable.
#'
#' @param dt A \code{Date} or date-coercible value specifying the day to load.
#' @param input_path Directory containing preprocessed HRRR GeoTIFF files (output
#'   of \code{\link{hrrr_acquire}}).
#' @param variables Character vector of variable names to load. Files are expected
#'   to follow the naming convention \code{VARNAME_YYYYMMDD_hrrrsmoke.tif}.
#'   Defaults to \code{c("APCP", "MASSDEN", "TMP", "VGRD", "UGRD", "RH", "HPBL")}.
#'
#' @returns A \code{SpatRaster} with one layer per variable, named by
#'   \code{variables}.
#' @export
#'
#' @examples
#' \dontrun{
#' hrrr_stack("2024-11-15", input_path = "./processed_data/HRRR/")
#' }
hrrr_stack <- function(dt, input_path, variables = c("APCP", "MASSDEN", "TMP", "VGRD",
                                                     "UGRD", "RH", "HPBL")) {

  files <- paste(variables, strftime(dt, "%Y%m%d"), "hrrrsmoke.tif", sep = "_")
  files <- fs::path(input_path, files)
  stack <- terra::rast(files)
  names(stack) <- variables
  return(stack)
  
}

