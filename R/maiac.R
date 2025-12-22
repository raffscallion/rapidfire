

#' Download MAIAC AOD data from NASA EartData archive and save it as tif
#'
#' @param dt character or date - date of MAIAC data to acquire
#' @param user username for NASA EarthData login
#' @param password password for NASA EartData login
#' @param outpath output path for files
#' @param bounding_box a vector of coordinates, x1, y1, x2, y2 to bound the request. The
#'   default bounding box covers California
#' @param tz time zone for the requested location (default is America/Los_Angeles)
#'
#' @returns Returns the list of tif files and downloads files to disk
#' @export
#'
#' @examples maiac_acquire("2025-09-05", user = u, password = p)
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
  
  download_raster <- function(url) {
    raster <- terra::rast(url)
    fname <- stringr::str_split_1(url, "/")
    fname <- fname[length(fname)]
    fname <- paste0(fs::path_ext_remove(fname), ".tif")
    outfile <- fs::path_join(c(outpath, fname))
    write_new_raster(raster, outfile)
    return(outfile)
  }
  
  # return the list of tif files downloaded (or skipped)
  files <- purrr::map_chr(urls, download_raster, .progress = TRUE)

}

# Merge tiles and fill gaps. We keep in native projection to preserve info.
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

# This version fills gaps using surrounding data in stages with increasing window sizes.
# 5x5, then 9x9, then 25x25. Finally filling the remainder with the median value
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

maiac_regrid <- function(r, extent_grid) {
  
  # reproject to the same coords
  coords <- terra::crs(extent_grid)
  r <- terra::project(r, coords)
  r <- terra::resample(r, extent_grid)
  
}



# Not sure if these are still used
safe_cut <- purrr::possibly(`[`, otherwise = NULL)
# Need to handle bad files with a passthrough here


# Read a MAIAC file from disk and flatten
maiac_aod <- function(fname) {
  
  sds <- terra::sds(fname)
  stack <- safe_cut(sds, 1)
  if (!is.null(stack)) {
    # There are multiple overpasses in the file, so we combine them with avg
    r <- terra::mean(stack, na.rm = TRUE)
    return(r)
  } else {
    return(NULL)
  }
  
}

extract_maiac <- function(maiac, locs) {

  sin_proj <- sp::CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  loc_sin <- sp::spTransform(locs, sin_proj)
  vals <- terra::extract(maiac, loc_sin@coords)

}

maiac_one_day <- function(dt, input_path = "./data/MAIAC/",
                          tiles_needed = c("h08v04", "h08v05", "h09v04"),
                          allow_missing = NULL) {

  date_str <- strftime(dt, format = "%Y%j")
  file_prefix <- paste0("MCD19A2.A", date_str)

  all_files <- fs::path_file(fs::dir_ls(input_path))

  pattern <- paste0(file_prefix, ".(?:", paste(tiles_needed, collapse = "|"), ").*hdf")

  files <- stringr::str_subset(all_files, pattern = pattern)
  paths <- paste0(input_path, files)

  # If no MAIAC data for date and allow_missing set, load dummy data with AOD = 0.15
  if (length(files) == 0) {
    if ("MAIAC" %in% allow_missing) {
      warning("No MAIAC data for ", strftime(dt, "%Y-%m-%d"))
      return(NULL)
    } else {
      stop("No MAIAC data for ", strftime(dt, "%Y-%m-%d"))
    }
  }
  
  tiles <- purrr::map(paths, maiac_aod)
  # Remove bad tiles
  tiles <- Filter(Negate(is.null), tiles)
  maiac <- maiac_mosaic(tiles) %>%
    maiac_fill_gaps_complete()

}

#' maiac_at_airnow
#'
#' Extract aerosol optical depth values at point locations on given dates from
#' pre-downloaded MAIAC aerosol tiles (MCD19A2).
#'
#' @param an A SpatialPointsDataFrame with monitor data such as from
#'   \code{\link{recast_monitors}}
#' @param maiac_path The path to the MAIAC data (defaults to "./data/MAIAC/")
#' @param tiles_needed The MODIS tiles required for the region of interest.
#'   Default is c("h08v04", "h08v05", "h09v04") which covers California.
#'
#' @return The data frame from \emph{an} with the extracted values from the
#'   MAIAC data appended
#' @export
#'
#' @examples  maiac <- maiac_at_airnow(mon)
maiac_at_airnow <- function(an, maiac_path = "./data/MAIAC/",
                            tiles_needed = c("h08v04", "h08v05", "h09v04"),
                            allow_missing = NULL) {

  daily_extract <- function(dt) {
    print(dt)
    maiac_aod <- maiac_one_day(dt, input_path = maiac_path, tiles_needed = tiles_needed,
                               allow_missing = allow_missing)
    i <- an$Day == dt
    locs = an[i, ]
    if (is.null(maiac_aod)) {
      df <- locs@data |>
        mutate(MAIAC_AOD = NA)
    } else {
      e <- extract_maiac(maiac_aod, locs)
      df <- locs@data %>%
        mutate(MAIAC_AOD = e$focal_mean)
    }
  }

  dates <- unique(an$Day)
  purrr::map_dfr(dates, daily_extract)

}

maiac_at_grid <- function(start, end, grid, maiac_path = "./data/MAIAC/",
                          tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  daily_extract <- function(dt) {
    print(dt)
    maiac_aod <- maiac_one_day(dt, input_path = maiac_path, tiles_needed = tiles_needed)
    r <- raster::raster(maiac_aod)

    on_grid <- extract_maiac(maiac_aod, grid)

    df <- grid@data %>%
      mutate(MAIAC_AOD = on_grid$mean,
             Day = dt)
  }

  dates <- seq.Date(from = start, to = end, by = "1 day")
  purrr::map_dfr(dates, daily_extract)


}
