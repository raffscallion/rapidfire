
# Download the NARR data we need
# Surface temp - air.2m
# Surface winds - uwnd.10m and vwnd.10m
# Surface RH - rhum.2m
# Precip rate - apcp
# PBLH - hpbl

get_narr <- function(year, vars = c("air.2m", "uwnd.10m", "vwnd.10m", "rhum.2m",
                                    "apcp", "hpbl"),
                     out_path = "./data/NARR") {

  fnarr <- function(variable, year, path) {
    stop("TODO: Replace NARR download with avaialble code. I believe new package is called RNCEP")
    # Reanalysis::NARRdownloadNetCDF(startYear = year, endYear = year,
    #                                variable = variable, destination = path)
  }

  purrr::walk(vars, fnarr, year = year, path = out_path)

}

extract_narr <- function(date, var, locs, input_path = "./data/NARR/") {

  filename <- paste0(input_path, var, ".", lubridate::year(date), ".nc")
  b <- raster::brick(filename)

  doy <- lubridate::yday(date)
  r <- b[[doy]]
  raster::extract(r, locs)

}

#' narr_at_airnow
#'
#' Extract surface temperature, winds, RF, precip, and planetary boundary layer
#' height from pre-downloaded NARR data
#'
#' @param an_data A SpatialPointsDataFrame with monitor data such as from
#'   \code{\link{recast_monitors}}
#'
#' @return The data frame from \emph{an_data} with the extracted values from the
#'   MAIAC data appended
#' @export
#'
#' @examples narr <- narr_at_airnow(mon)
narr_at_airnow <- function(an_data) {

  date <- unique(an_data$Day)

  # for each date and each variable
  var = c("air.2m", "uwnd.10m", "vwnd.10m", "rhum.2m",
           "apcp", "hpbl")
  cases <- tidyr::expand_grid(date, var)

  one_day <- function(date, var) {
    i <- an_data$Day == date
    l = an_data[i, ]
    vals <- extract_narr(date, var, l)
    df <- l@data
    df$narr_var <- vals
    df$narr_name <- var
    df
  }
  purrr::pmap_dfr(cases, one_day) %>%
    tidyr::pivot_wider(names_from = narr_name, values_from = narr_var)

}


narr_at_grid <- function(start, end, grid) {

  date <- seq.Date(from = start, to = end, by = "1 day")
  # for each date and each variable
  var = c("air.2m", "uwnd.10m", "vwnd.10m", "rhum.2m",
          "apcp", "hpbl")
  cases <- tidyr::expand_grid(date, var)
  one_day <- function(date, var) {
    vals <- extract_narr(date, var, grid)
    df <- grid@data
    df$narr_var <- vals
    df$narr_name <- var
    df$Day <- date
    df
  }
  purrr::pmap_dfr(cases, one_day) %>%
    tidyr::pivot_wider(names_from = narr_name, values_from = narr_var)

}
