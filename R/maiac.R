

#' Download MAIAC AOD data from the NASA EarthData archive and save as GeoTIFF
#'
#' Queries the NASA CMR API for MAIAC AOD granules (MCD19A2 C6.1) intersecting a
#' bounding box and date. If standard MAIAC data are not yet available, falls back to
#' the Near Real-Time product (MCD19A2N C6.1NRT). Matching HDF granules are downloaded
#' and written to disk as GeoTIFF files.
#'
#' @param dt A \code{Date} or date-coercible character string specifying the local day
#'   to acquire.
#' @param user NASA EarthData username. If \code{NULL}, uses credentials previously
#'   configured by \code{earthdatalogin::edl_netrc}.
#' @param password NASA EarthData password. If \code{NULL}, uses stored credentials.
#' @param outpath Directory path where downloaded GeoTIFF files will be saved.
#'   Default is \code{"./raw_data/MAIAC/"}.
#' @param bounding_box A numeric vector of four coordinates
#'   \code{c(xmin, ymin, xmax, ymax)} in WGS84 (EPSG:4326) defining the spatial
#'   extent of the request. Default covers California.
#' @param tz Time zone string for the requested location, used to convert \code{dt}
#'   to UTC for the CMR query. Default is \code{"America/Los_Angeles"}.
#'
#' @returns A character vector of file paths for the downloaded GeoTIFF files.
#'   Files that fail to download return an error message string in place of a path.
#' @export
#'
#' @examples
#' \dontrun{
#' maiac_acquire(
#'   "2025-09-05",
#'   user = Sys.getenv("EARTHDATA_USER"),
#'   password = Sys.getenv("EARTHDATA_PASSWORD")
#' )
#' }
maiac_acquire <- function(dt, user = NULL, password = NULL, outpath = "./raw_data/MAIAC/",
                          bounding_box = c(-124.5, 32.5, -114, 42.1),
                          tz = "America/Los_Angeles") {
  
  # Construct date query based on time zone
  date_local <- as.POSIXct(dt, tz = tz)
  end_date_local <- date_local + lubridate::days(1)
  
  # Date
  dt_str <- strftime(date_local, tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
  end_dt_str <- strftime(end_date_local, tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
  temporal_str <- paste(dt_str, end_dt_str, sep = ",")
  
  # Spatial query
  bounding_box_str <- paste(bounding_box, collapse = ",")
  
  # set up credentials to EarthData archives
  earthdatalogin::edl_netrc(username = user, password = password)
  #earthdatalogin::edl_netrc()
  
  # First check for the regular MAIAC files (MCD19A2 C6.1). If those aren't ready yet,
  # look for the NRT versions (MCD19A2N C6.1NRT)
  concept_id <- "C2324689816-LPCLOUD"
  
  url <- paste0("https://cmr.earthdata.nasa.gov/search/granules.json?collection_concept_id=",
                concept_id, "&temporal=", temporal_str, "&bounding_box=", bounding_box_str)
  
  result <- RCurl::getURL(url)
  js <- jsonlite::fromJSON(result)
  links <- js$feed$entry$links
  
  if (is.null(links)) {
    warning("No MAIAC data found, looking for NRT version...")
    
    # Look for NRT data
    concept_id <- "C2407807500-LANCEMODIS"
    url <- paste0("https://cmr.earthdata.nasa.gov/search/granules.json?collection_concept_id=",
                  concept_id, "&temporal=", temporal_str, "&bounding_box=", bounding_box_str)
    result <- RCurl::getURL(url)
    js <- jsonlite::fromJSON(result)
    links <- js$feed$entry$links
    
    if (is.null(links)) {
      warning("No MAIAC NRT data found either...")
      return(NULL)
    }
    
    df <- purrr::list_rbind(links)
    df <- df |>
      dplyr::filter(stringr::str_ends(href, ".hdf"))

  } else {
    df <- purrr::list_rbind(links)
    df <- df |>
      dplyr::filter(stringr::str_starts(title, "Download") & stringr::str_ends(title, ".hdf"))
    
  }
  
  urls <- df$href

  # skip any files we already have
  write_new_raster <- purrr::possibly(terra::writeRaster, otherwise = NULL)
  safe_download <- purrr::safely(terra::rast)
  
  download_raster <- function(url) {
    
    results <- safe_download(url) 

    if (is.null(results$error)) {
      raster <- results$result
      fname <- stringr::str_split_1(url, "/")
      fname <- fname[length(fname)]
      fname <- paste0(fs::path_ext_remove(fname), ".tif")
      outfile <- fs::path_join(c(outpath, fname))
      write_new_raster(raster, outfile)
    } else {
      outfile <- results$error$message
    }
    return(outfile)    
  }
  
  # return the list of tif files downloaded (or skipped)
  files <- purrr::map_chr(urls, download_raster, .progress = TRUE)

}

#' Merge MAIAC tiles and fill spatial gaps for a single day
#'
#' Reads all MAIAC GeoTIFF files for a given date from \code{input_path}, extracts the
#' 470 nm AOD layers (\code{Optical_Depth_047}), averages multiple overpasses per tile,
#' mosaics the tiles, and fills remaining gaps using progressive focal windows. The
#' result is saved as a single-layer GeoTIFF in the native MAIAC projection.
#'
#' @param dt A \code{Date} or date-coercible character string specifying the day to
#'   process.
#' @param input_path Directory containing raw MAIAC GeoTIFF files as downloaded by
#'   \code{\link{maiac_acquire}}. Default is \code{"./raw_data/MAIAC/"}.
#' @param output_path Directory where the processed GeoTIFF will be saved.
#'   Default is \code{"./processed_data/MAIAC/"}.
#' @param overwrite Logical. If \code{TRUE}, overwrite an existing output file. If
#'   \code{FALSE} (default), a warning is issued and processing is skipped when the
#'   output file already exists.
#'
#' @returns A \code{SpatRaster} with a single layer named \code{"MAIAC_AOD"} containing
#'   merged, gap-filled 470 nm AOD values in the native MAIAC sinusoidal projection.
#'   The raster is also written to \code{output_path} as
#'   \code{MAIAC_processed_YYYY-MM-DD.tif}.
#' @export
#'
#' @examples
#' \dontrun{
#' maiac_preprocess("2025-09-05")
#' }
maiac_preprocess <- function(dt, input_path = "./raw_data/MAIAC/", 
                             output_path = "./processed_data/MAIAC/",
                             overwrite = FALSE) {
  
  # Merge all passed tiles - confirm they are the same date
  date_str <- paste0("A", strftime(dt, format = "%Y%j"))
  reg <- paste0(".+", date_str, ".+tif")
  files <- fs::dir_ls(input_path, regexp = reg)

  process_tiles <- function(file) {
    r_raw <- terra::rast(file)
    # Get all the layers for AOD at 470 nm (the native for the algorithm)
    layers_needed <- which(strtrim(terra::names(r_raw), 17) == "Optical_Depth_047")
    r_sub <- terra::subset(r_raw, layers_needed)
    # Combine overpasses with mean
    r <- terra::mean(r_sub, na.rm = TRUE)
  }
  tiles <- purrr::map(files, process_tiles)
  
  # Mosiac the tiles
  maiac <- purrr::reduce(tiles, terra::mosaic, fun = "mean")

  # Fill gaps
  maiac <- maiac_fill_gaps_complete(maiac)
  
  # Rename raster layer
  names(maiac) <- "MAIAC_AOD"
  
  # Export file
  outfile <- fs::path_join(c(output_path, 
                             paste0("MAIAC_processed_", 
                                    strftime(dt, format = "%Y-%m-%d"), ".tif")))
  
  if (overwrite) {
    terra::writeRaster(maiac, outfile, overwrite = TRUE)
  } else {
    if (fs::file_exists(outfile)) {
      warning("MAIAC processed file already exists. Skipping.")
    } else {
      terra::writeRaster(maiac, outfile, overwrite = FALSE)  
    }
  }
  
  return(maiac)
  
}

#' Fill gaps in a MAIAC AOD raster using progressive focal smoothing
#'
#' Fills \code{NA} cells in stages using focal mean windows of increasing size
#' (5×5, then 9×9, then 25×25), applying each pass only to cells that remain
#' \code{NA} after the previous pass. Any cells still missing after all three
#' passes are filled with the global median of the raster.
#'
#' @param maiac A \code{SpatRaster} with AOD values, potentially containing
#'   \code{NA} gaps (e.g., from cloud cover).
#'
#' @returns A \code{SpatRaster} of the same extent and resolution as \code{maiac}
#'   with all \code{NA} cells filled.
#'
#' @seealso \code{\link{maiac_preprocess}}
maiac_fill_gaps_complete <- function(maiac) {

  md <- terra::global(maiac, fun = median, na.rm = TRUE) %>%
    .$global

  fill1 <- terra::focal(maiac, w = 5, fun = "mean", na.rm = TRUE, na.policy = "only")
  blanks <- is.na(fill1)
  fill1[blanks] <- NA
  fill2 <- terra::focal(fill1, w = 9, fun = "mean", na.rm = TRUE, na.policy = "only")
  blanks <- is.na(fill2)
  fill2[blanks] <- NA
  fill3 <- terra::focal(fill2, w = 25, fun = "mean", na.rm = TRUE, na.policy = "only")
  blanks <- is.na(fill3)
  med_fill <- blanks * md
  fill3 <- terra::subst(fill3, NA, 0)
  final <- fill3 + med_fill
  

}

#' Reproject and resample a MAIAC raster to match a target grid
#'
#' Reprojects a MAIAC \code{SpatRaster} to the CRS of \code{extent_grid}, resamples
#' it to match the target extent and resolution, and replaces any \code{NA} values
#' introduced during resampling with the raster's global median.
#'
#' @param r A \code{SpatRaster} to reproject and resample (e.g., output of
#'   \code{\link{maiac_preprocess}}).
#' @param extent_grid A \code{SpatRaster} defining the target CRS, extent, and
#'   resolution (e.g., the model prediction grid).
#'
#' @returns A \code{SpatRaster} with the same CRS, extent, and resolution as
#'   \code{extent_grid}, with any \code{NA} cells filled by the global median.
#'
#' @seealso \code{\link{maiac_preprocess}}
maiac_regrid <- function(r, extent_grid) {

  # reproject to the same coords
  coords <- terra::crs(extent_grid)
  r <- terra::project(r, coords)
  r <- terra::resample(r, extent_grid)
  
  # If any NAs result, replace with median
  blanks <- is.na(r)
  md <- terra::global(r, fun = median, na.rm = TRUE) %>%
  .$global
  r <- terra::subst(r, NA, md)
  
}