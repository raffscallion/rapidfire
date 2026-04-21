# Compare with other products


haqes_acquire <- function(dt, format = "census", product = "PM25_TOT", 
                          include_fips = "all") {
  
  base_url <- "http://air.csiss.gmu.edu/yli/tempdata/ens_data/"
  
  path <- paste0(base_url, format, "/", product)
  
  hrs <- seq(0, 21, by = 3)
  
  #Example: HAQES_NA_v1.0_3h_PM25_TOT_COUNTY.20220102.21z.txt	
  files <- glue::glue("{path}/HAQES_NA_v1.0_3h_{product}_{toupper(format)}.",
                      "{strftime(dt, format = '%Y%m%d')}.",
                      "{sprintf('%02d', hrs)}z.txt")
  
  read_file <- function(f, hr) {
    
    # read the first bit to see where the header ends -  the header is surrounded by `=====`
    safe_read <- purrr::possibly(readLines)
    first_part <- safe_read(f, n = 40)
    if (is.null(first_part)) {
      return(NULL)
    } 
    
    header_ends <- grep("^+=", first_part)[2]
    
    # read in the data
    df <- readr::read_csv(f, skip = header_ends, col_types = "cn", show_col_types = FALSE)
    
    df <- df |>
      mutate(hr = hr)

  }
  
  df <- purrr::map2(files, hrs, \(x, y) read_file(x, y), .progress = TRUE) |>
    purrr::list_rbind()
  
  if (nrow(df) == 0) {
    warning("No data for ", dt)
    return(NULL)
  }

  param <- switch(product,
                  "PM25_TOT" = "pm25",
                  "PM25_OC" = "pm25oc",
                  "PM25_BC" = "pm25bc")
  
  unit <- switch(format,
                 "census" = "TRACT",
                 "county" = "FIPS")
  
  # Limit to selected FIPS
  if (include_fips != "all") {
    
    df <- df |>
      filter(grepl(pattern = paste0("^", include_fips), .data[[unit]]))
    
  }
  
  # Summarize by spatial unit
  df <- df |>
    summarise({{param}} := mean(.data[[param]], na.rm = TRUE),
              .by = {{unit}})
  
}

haqes_archive <- function(start_dt, end_dt, output_path, format = "census", 
                          product = "PM25_TOT", include_fips = "06") {
  
  dts <- seq.Date(as.Date(start_dt), as.Date(end_dt), by = "1 day")
  
  download_haqes <- function(dt, format, product, include_fips, output_path) {
    
    df <- haqes_acquire(dt, format, product, include_fips)
    filename <- paste0("haqes_", strftime(dt, format = "%Y-%m-%d"), ".RDS")
    file_out <- fs::path_join(c(output_path, filename))
    saveRDS(df, file_out)
    
  }
  
  purrr::walk(dts, \(x) download_haqes(x, format, product, include_fips, output_path))
  
}

haqes_compare <- function(dt, rf_path, haqes_path, gis_data) {
  
  haqes_file <- fs::path_join(c(haqes_path, paste0("haqes_", strftime(dt, "%Y-%m-%d"), ".RDS")))
  haqes <- readRDS(haqes_file)
  
  rf_file <- fs::path_join(c(rf_path, paste0("rapidfire_", strftime(dt, "%Y-%m-%d"), ".tif")))
  rf <- terra::rast(rf_file)
  
  crs <- terra::crs(rf)
  if (terra::crs(gis_data) != crs) {
    
    gis_data <- terra::project(gis_data, crs)
  }
  
  ex <- terra::extract(rf, gis_data, fun = mean, layer = "PM25_RF")
  
  tracts <- as.data.frame(gis_data) |>
    select(TRACT=GEOID) |>
    mutate(PM25_RF = ex$value) |>
    left_join(haqes, by = "TRACT") |>
    rename(PM25_HAQES=pm25)
  
  
}

