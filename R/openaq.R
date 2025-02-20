# Functions for handling PurpleAir data via the OpenAQ archive API. As of early
# 2023, the low-cost sensors in OpenAQ go back to early 2021. Prior to that, the
# AirSensors package can be used.
#

#' Query the OpenAQ API to get all of the low-cost sensors that have reported
#' PM2.5 within a geographic boundary.
#'
#' @param outline_sf A simple features polygon, such as a state, preferably in
#'   an appropriate equal area projection
#' @param search_radius The search radius in meters of each query. The max
#'   allowed by the API is 100,000. However, this will fail for regions with
#'   very dense networks, such as the San Francisco Bay Area. The default is
#'   20,000.
#' @param country_id A two-letter country code. Default is "US". Note that only
#'   one country can be queried at a time.
#'
#' @return A data frame of currently reporting stations
#' @export
#'
#' @examples ca <- sf::st_read("states_hi_res.shp")
#'           ca <- dplyr::filter(ca, STATE_NAME == "California")
#'           ca <- dplyr::select(ca, OBJECTID, STATE_NAME)
#'           proj_string <- "+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
#'           ca_ea <- sf::st_transform(ca, crs = proj_string)
#'           sites <- openaq_find_sites(ca_ea)
#'           date_from <- "2021-02-01"
#'           date_to <- "2021-02-28"
#'           avgs <- openaq_get_averages(sites, date_from, date_to)
openaq_find_sites <- function(outline_sf, search_radius = 20000,
                              country_id = "US") {

  # Create a hex grid to cover the area of outline_sf
  # We would like to end up with circles of radius = search_radius
  # What we need to provide st_make_grid is the diameter of the hexagon.
  # Half of the diameter is the hexagon radius (Ri)
  # D = Ri * 2
  # The circle radius (Rc, search_radius) is related to the hex radius as
  # Ri ~= Rc * 0.866
  # See https://calcresource.com/geom-hexagon.html
  D <- (search_radius * 0.867) * 2
  hex <- sf::st_make_grid(outline_sf, cellsize = c(D, D), square = FALSE)
  hex_int <- sf::st_intersects(hex, outline_sf, sparse = FALSE)
  hex_grid <- hex[hex_int]

  # Centroids of this grid are the locations of the query
  centroids <- sf::st_centroid(hex_grid)

  # Get lat/lons to pass to the API
  centroid_ll <- sf::st_transform(centroids, crs = "+init=epsg:4326")

  # For each coordinate, query the sites within the radius
  query_locations <- function(coord, radius, country_id) {

    coord_string <- paste(round(coord[2], 5), round(coord[1], 5), sep = ",")
    query_string <- list(
      limit = "100000",
      offset = "0",
      sort = "desc",
      parameter_id = "2",
      coordinates = coord_string,
      radius = as.character(radius),
      country = country_id,
      order_by = "location",
      sensorType = "low-cost sensor",
      dumpRaw = "false"
    )

    url <- "https://api.openaq.org/v2/locations"
    response <- httr::VERB("GET", url, query = query_string,
                     httr::content_type("application/octet-stream"),
                     httr::accept("application/json"))
    locs <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    if (is.data.frame(locs$results)) {
      df <- locs$results %>%
        tidyr::unnest(cols= coordinates)
      return(df)
    } else {
      return(NULL)
    }
  }

  sites <- purrr::map(centroid_ll, query_locations, radius = search_radius,
                      country_id = country_id, .progress = "Finding sensors") %>%
    purrr::list_rbind()

  # Since the circles have overlap, remove duplicates
  sites_unique <- unique(sites)

}


#' For the sites found with \code{\link{openaq_find_sites}}, download daily
#' average PM2.5 data from the OpenAQ API
#'
#' @param sites A data frame returned by \code{\link{openaq_find_sites}}. The
#'   relevant fields are id, which is the OpenAQ location id, plus latitude and
#'   longitude
#' @param date_from Earliest date to query in character format "YY-mm-dd"
#' @param date_to Latest date to query in character format "YY-mm-dd"
#'
#' @return A data frame with daily avereage PM2.5
#' @export
#'
#' @examples ca <- sf::st_read("states_hi_res.shp")
#'           ca <- dplyr::filter(ca, STATE_NAME == "California")
#'           ca <- dplyr::select(ca, OBJECTID, STATE_NAME)
#'           proj_string <- "+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
#'           ca_ea <- sf::st_transform(ca, crs = proj_string)
#'           sites <- openaq_find_sites(ca_ea)
#'           date_from <- "2021-02-01"
#'           date_to <- "2021-02-28"
#'           avgs <- openaq_get_averages(sites, date_from, date_to)
openaq_get_averages <- function(sites, date_from, date_to) {

  location_ids <- sites$id

  # Get multiple locations per call
  locs_per <- 40
  batched <- split(location_ids, ceiling(seq_along(location_ids)/locs_per))

  get_averages <- function(locations, date_from, date_to) {

    # Using the averages api, we can get the latest 24-hr averages from these sites
    url <- "https://api.openaq.org/v2/averages"

    # build the location string - needs to look like this: "236295&location=369354"
    loc_string <- paste(locations, collapse = "&location=")

    queryString <- list(
      date_from = date_from,
      date_to = date_to,
      parameter = "pm25",
      country = "US",
      limit = "100000",
      location = loc_string,
      page = "1",
      offset = "0",
      sort = "asc",
      spatial = "location",
      temporal = "day",
      group = "false"
    )

    # build the query string ourselves since we need multiple locations and that
    # doesn't work with httr
    vals <- purrr::imap_chr(queryString, ~paste(.y, .x, sep = "="))
    full_query <- paste(vals, collapse = "&")
    full_url <- paste0(url, "?", full_query)

    response <- httr::VERB("GET", full_url,
                     httr::content_type("application/octet-stream"),
                     httr::accept("application/json"))
    df <- jsonlite::fromJSON(httr::content(response, "text",
                                           encoding = "UTF-8"))$results
    if (is.data.frame(df)) {
      return(df)
    } else {
      return(NULL)
    }

  }

  results <- purrr::map(batched, get_averages, date_from, date_to,
                        .progress = "Getting sensor data") %>%
    purrr::list_rbind()

  # The name field in the averages output maps to the id in the sites data
  site_coords <- sites %>%
    select(id, longitude, latitude) %>%
    mutate(name = as.character(id), .keep = "unused")

  ## export in the correct format:
  ## deviceID, Date, latitude, longitude, pm25_1day
  df <- results %>%
    inner_join(site_coords, by = "name") %>%
    mutate(Date = as.Date(day)) %>%
    select(deviceID=name, Date, latitude, longitude, pm25_1day=average, measurement_count)

}


