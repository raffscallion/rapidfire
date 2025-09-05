

#' Title
#'
#' @param dt 
#' @param user 
#' @param password 
#' @param outpath 
#' @param bounding_box 
#' @param tz 
#'
#' @returns
#' @export
#'
#' @examples
maiac_acquire <- function(dt, user, password, outpath = "./raw_data/MAIAC/",
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
  }
  
  purrr::walk(urls, download_raster, .progress = TRUE)

}


safe_cut <- purrr::possibly(`[`, otherwise = NULL)
# Need to handle bad files with a passthrough here

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

maiac_mosaic <- function(tiles) {

  purrr::reduce(tiles, terra::mosaic, fun = "mean")

}

# Use raster::focal to fill in some missing values by interpolating neighbors,
# Then replace the remaining missing values with the median value
maiac_fill_gaps <- function(maiac, window = 7) {

  md <- terra::global(maiac, fun = median, na.rm = TRUE) %>%
    .$global
  sr <- terra::focal(maiac, w = window, fun = "mean", na.rm = TRUE, na.only = TRUE)

  blank_space <- sr == 0
  fill <- blank_space * md
  filled <- sr + fill

}

# This version fills gaps using surrounding data in stages with increasing window sizes.
# 5x5, then 9x9, then 25x25. Finally filling the remainder with the median value
maiac_fill_gaps_complete <- function(maiac) {

  md <- terra::global(maiac, fun = median, na.rm = TRUE) %>%
    .$global

  fill1 <- terra::focal(maiac, w = 5, fun = "mean", na.rm = TRUE, na.only = TRUE)
  blanks <- fill1 == 0
  fill1[blanks] <- NA
  fill2 <- terra::focal(fill1, w = 9, fun = "mean", na.rm = TRUE, na.only = TRUE)
  blanks <- fill2 == 0
  fill2[blanks] <- NA
  fill3 <- terra::focal(fill2, w = 25, fun = "mean", na.rm = TRUE, na.only = TRUE)
  blanks <- fill3 == 0
  med_fill <- blanks * md
  final <- fill3 + med_fill

}


extract_maiac <- function(maiac, locs) {

  sin_proj <- sp::CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  loc_sin <- sp::spTransform(locs, sin_proj)
  vals <- terra::extract(maiac, loc_sin@coords)

}

maiac_one_day <- function(dt, input_path = "./data/MAIAC/",
                          tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  date_str <- strftime(dt, format = "%Y%j")
  file_prefix <- paste0("MCD19A2.A", date_str)

  all_files <- fs::path_file(fs::dir_ls(input_path))

  pattern <- paste0(file_prefix, ".(?:", paste(tiles_needed, collapse = "|"), ").*hdf")

  files <- stringr::str_subset(all_files, pattern = pattern)
  paths <- paste0(input_path, files)

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
                            tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  daily_extract <- function(dt) {
    print(dt)
    maiac_aod <- maiac_one_day(dt, input_path = maiac_path, tiles_needed = tiles_needed)
    i <- an$Day == dt
    locs = an[i, ]
    e <- extract_maiac(maiac_aod, locs)
    df <- locs@data %>%
      mutate(MAIAC_AOD = e$focal_mean)
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
