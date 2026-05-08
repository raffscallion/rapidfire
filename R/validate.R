
# Create a set of monitor predictions by sequential leave-one-out. Then, for each, can run
# prediction at missing location and write out both.

#' Title
#'
#' @param dt 
#' @param model_file 
#' @param paths 
#'
#' @returns
#' @export
#'
#' @examples
daily_cross_validate <- function(dt, model_file, paths) {
  
  # load monitoring data
  monitoring_file <- fs::path_join(c(paths[["monitors"]], 
                                     paste0("monitors_combined_",
                                            strftime(dt, "%Y-%m-%d"),
                                            ".RDS")))
  mon <- terra::readRDS(monitoring_file)

  # Remove any non-finite values
  mon <- mon[is.finite(mon$PM25_log),]
  
  # For each site, interpolate without it's data
  ids <- mon$monitorID
  
    # load the variogram
  models <- readRDS(model_file)
  mon_vgm <- models$vgm_mon
  
  leave_out_interpolate <- function(id) {
    
    mon_sub <- mon[mon$monitorID != id,]
    mon_one <- mon[mon$monitorID == id,]
    loc <- terra::crds(mon_one)
    
    res <- monitors_krige_points(mon_sub, mon_one, mon_vgm)
    df <- tibble(monitorID = id,
                 PM25_log_ANK = res$var1.pred)
  }
  
  interpolated_monitors <- purrr::map(ids, leave_out_interpolate, 
                                      .progress = "Leave-one-out interpolation") |>
    purrr::list_rbind()
  
  # load met and smoke data from hrrr
  hrrr <- hrrr_stack(dt, paths[["hrrr"]])
  
  # Convert MASSDEN to log scale
  massdenlog <- log(hrrr$MASSDEN)
  names(massdenlog) <- "MASSDEN_log"
  hrrr <- terra::rast(list(massdenlog, hrrr))
  hrrr <- terra::subset(hrrr, which(names(hrrr)=="MASSDEN"), negate = TRUE)
  hrrr_extracted <- terra::extract(hrrr, mon)
  
  # load satellite data
  maiac_file <- fs::path_join(c(paths[["maiac"]],
                                paste0("MAIAC_processed_", strftime(dt, "%Y-%m-%d"), ".tif")))
  maiac <- terra::rast(maiac_file)
  
  # interpolate to grid
  maiac <- maiac_regrid(maiac, grid)
  maiac_extracted <- terra::extract(maiac, mon)
  
  # load sensor data
  ### Sensor data will be coming from CDPH? Need non-cdph way as well
  pa_file <- fs::path_join(c(paths[["purpleair"]],
                             paste0("purpleair_processed_",
                                    strftime(dt, "%Y-%m-%d"),
                                    ".RDS")))
  pa <- terra::readRDS(pa_file)

  # interpolate to monitor points
  pa_vgm <- models$vgm_pa
  pak <- monitors_krige_points(pa, mon, pa_vgm, nmax = 100)

  # Compile and predict
  model_inputs <- as.data.frame(mon) |>
    select(monitorID, Day, PM25_log, source) |>
    left_join(interpolated_monitors, by = "monitorID")
  model_inputs$PM25_log_PAK <- pak$var1.pred
  model_inputs$MAIAC_AOD <- maiac_extracted$MAIAC_AOD
  model_inputs <- model_inputs |>
    bind_cols(select(hrrr_extracted, -ID))
  
  model_inputs <- model_inputs[complete.cases(model_inputs),]
  
  model <- models$model
  pred <- predict(model, model_inputs)
    
  # combine with observed inputs and return
  results <- model_inputs
  results$PM25_log_RF <- pred
  
  # Export file
  outfile <- fs::path_join(c(paths[["cross_validation"]], 
                             paste0("rapidfire_l-o-o_", 
                                    strftime(dt, format = "%Y-%m-%d"), ".csv")))
  readr::write_csv(results, outfile)
  
  return(results)
}


