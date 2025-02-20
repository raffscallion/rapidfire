
# Data pool location for MAIAC MCD19A2
# https://e4ftl01.cr.usgs.gov/MOTA/MCD19A2.006/


# For a given bounding box, return the MODIS tiles needed for coverage using the
# MODIS siusoidal grid. Definition found at
# https://modis-land.gsfc.nasa.gov/pdf/sn_bound_10deg.txt
modis_tile_pick <- function(llat, ulat, llon, rlon) {

  grid <- read.table("./data/MAIAC/modis_grid/sn_bound_10deg.txt", header = TRUE,
                   row.names = NULL, sep = "", skip = 6, nrow = 648)
  # Get tile for each corner
  ur <- grid %>%
    filter(lon_min < rlon, lon_max > rlon,
           lat_min < ulat, lat_max > ulat)
  lr <- grid %>%
    filter(lon_min < rlon, lon_max > rlon,
           lat_min < llat, lat_max > llat)
  ul <- grid %>%
    filter(lon_min < llon, lon_max > llon,
           lat_min < ulat, lat_max > ulat)
  ll <- grid %>%
    filter(lon_min < llon, lon_max > llon,
           lat_min < llat, lat_max > llat)

  all <- bind_rows(ur, lr, ul, ll) %>%
    distinct()

}

#' maiac_download
#'
#' Download MAIAC collection 6.1 data from the USGS archive
#'
#' @param dt Date to retrieve
#' @param user Username NASA Earthdata
#' @param password Password for NASA Earthdata
#' @param outpath output path to write files
#' @param tiles_needed specific MODIS tiles required for your domain. See
#'   https://modis-land.gsfc.nasa.gov/pdf/sn_bound_10deg.txt
#'
#' @return
#' @export
#'
#' @examples
maiac_download <- function(dt, user, password, outpath = "./data/MAIAC/",
                           tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  base_url <- "https://e4ftl01.cr.usgs.gov/MOTA/MCD19A2.061/"

  # Date
  dt_str <- strftime(dt, format = "%Y.%m.%d")

  # Need to get contents of the folder to create the correct filenames
  folder <- paste0(base_url, dt_str)
  all_files <- rvest::read_html(folder) %>%
    rvest::html_elements("a") %>%
    rvest::html_text()

  # File the filenames that fit the pattern for this date and tile
  dt_str2 <- strftime(dt, format = "%Y%j")
  patterns <- paste("MCD19A2.A", dt_str2, ".", tiles_needed, ".+hdf$",
                    sep = "")
  files_needed <- purrr::map_chr(patterns,
                                 ~stringr::str_subset(all_files, pattern = .x))

  # If the file already exists on disk, just move on to the next one
  write_disk_w_skip <- purrr::possibly(httr::write_disk, NULL)

  # Now, download each file
  download_one <- function(filename, user, pw) {
    fileUrl <- paste0(folder, "/", filename)
    httr::GET(fileUrl, httr::authenticate(user, pw),
              write_disk_w_skip(paste0(outpath, filename)), httr::timeout(60))
  }

  purrr::walk(files_needed, download_one, user, password)

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
