# Functions either directly copied or adapted from PWFSLSmoke or MazamCoreUtils. rapidfire
# used PWFSLSmoke to acquire AirNow and AirSIS data from the AirFire-curated archive. As
# of the deprecation of rgdal and rgeos, PWFSLSmoke is not presently working, so we are
# replicating some functionality here.

#' @keywords AirNow
#' @export
#' @title Load annual AirNow monitoring data
#' @param year Desired year (integer or character representing YYYY).
#' @param parameter Parameter of interest.
#' @param baseUrl Base URL for 'annual' AirNow data files.
#' @param dataDir Local directory containing 'annual' data files.
#' @return A \emph{ws_monitor} object with AirNow data.
#' @description Loads pre-generated .RData files containing annual
#' AirNow data.
#'
#' If \code{dataDir} is defined, data will be loaded from this local
#' dirctory. Otherwise, data will be loaded from the monitoring data repository
#' maintained by PWFSL.
#'
#' The annual files loaded by this function are updated on the 15'th of each
#' month and cover the period from the beginning of the year to the end of the
#' last month.
#'
#' For data during the last 45 days, use \code{airnow_loadDaily()}.
#'
#' For the most recent data, use \code{airnow_loadLatest()}.
#'
#' AirNow parameters include the following:
#' \enumerate{
# #' \item{BARPR}
# #' \item{BC}
# #' \item{CO}
# #' \item{NO}
# #' \item{NO2}
# #' \item{NO2Y}
# #' \item{NO2X}
# #' \item{NOX}
# #' \item{NOOY}
# #' \item{OC}
# #' \item{OZONE}
# #' \item{PM10}
#' \item{PM2.5}
# #' \item{PRECIP}
# #' \item{RHUM}
# #' \item{SO2}
# #' \item{SRAD}
# #' \item{TEMP}
# #' \item{UV-AETH}
# #' \item{WD}
# #' \item{WS}
#' }
#'
#' Available AirNow RData and associated log files can be seen at:
#' \href{https://haze.airfire.org/monitoring/AirNow/RData/}{https://haze.airfire.org/monitoring/AirNow/RData/}
#' @seealso \code{\link{airnow_loadDaily}}
#' @seealso \code{\link{airnow_loadLatest}}
#' @examples
#' \dontrun{
#' # Fail gracefully if any resources are not available
#' try({
#'
#' airnow_loadAnnual(2017) 
#'
#' }, silent = FALSE)
#' }

airnow_loadAnnual <- function(year = NULL,
  parameter = 'PM2.5',
  baseUrl = 'https://haze.airfire.org/monitoring',
  dataDir = NULL) {

# ----- Validate parameters --------------------------------------------------

if ( is.null(year) ) {
stop("Required parameter 'year' is missing.")
}
year <- as.numeric(year)

validParams <- c("PM2.5")
if ( !parameter %in% validParams ) {
paramsString <- paste(validParams, collapse=", ")
stop(paste0("'", parameter,
"' is not a supported parameter. Use 'parameter = ",
paramsString, "'"), call.=FALSE)
}

if ( year < 2016 ) {
stop("PWFSL has no annual data files for AirNow prior to 2016.")
}

# Load the data --------------------------------------------------------------


# Create URL and filename according to the PWFSLSmoke naming scheme
baseUrl <- paste0(baseUrl, '/AirNow/RData/', year)
filename <- paste0("airnow_", parameter, "_", year, ".RData")

ws_monitor <- loadDataFile(filename, baseUrl, dataDir)

return(ws_monitor)

}

#' @keywords AirNow
#' @export
#' @title Load recent AirNow monitoring data
#' @param parameter Parameter of interest.
#' @param baseUrl Base URL for 'daily' AirNow data files.
#' @param dataDir Local directory containing 'daily' data files.
#' @return A \emph{ws_monitor} object with AirNow data.
#' @description Loads pre-generated .RData files containing recent
#' AirNow data.
#'
#' If \code{dataDir} is defined, data will be loaded from this local
#' dirctory. Otherwise, data will be loaded from the monitoring data repository
#' maintained by PWFSL.
#'
#' The daily files loaded by this function are updated once a day, shortly
#' after midnight and contain data for the previous 45 days.
#'
#' For the most recent data, use \code{airnow_loadLatest()}.
#'
#' For data extended more than 45 days into the past, use \code{airnow_loadAnnual()}.
#'
#' AirNow parameters include the following:
#' \enumerate{
# #' \item{BARPR}
# #' \item{BC}
# #' \item{CO}
# #' \item{NO}
# #' \item{NO2}
# #' \item{NO2Y}
# #' \item{NO2X}
# #' \item{NOX}
# #' \item{NOOY}
# #' \item{OC}
# #' \item{OZONE}
# #' \item{PM10}
#' \item{PM2.5}
# #' \item{PRECIP}
# #' \item{RHUM}
# #' \item{SO2}
# #' \item{SRAD}
# #' \item{TEMP}
# #' \item{UV-AETH}
# #' \item{WD}
# #' \item{WS}
#' }
#'
#' Available AirNow RData and associated log files can be seen at:
#' \href{https://haze.airfire.org/monitoring/AirNow/RData/latest/}{https://haze.airfire.org/monitoring/AirNow/RData/latest/}
#' @seealso \code{\link{airnow_loadAnnual}}
#' @seealso \code{\link{airnow_loadLatest}}
#' @examples
#' \dontrun{
#' # Fail gracefully if any resources are not available
#' try({
#'
#' airnow_loadDaily() %>%
#'   monitor_subset(stateCodes=CONUS) %>%
#'   monitor_map()
#'
#' }, silent = FALSE)
#' }

airnow_loadDaily <- function(parameter = 'PM2.5',
                             baseUrl = 'https://haze.airfire.org/monitoring/latest/RData',
                             dataDir = NULL) {

  # Validate parameter
  validParams <- c("PM2.5")
  if ( !parameter %in% validParams ) {
    paramsString <- paste(validParams, collapse=", ")
    stop(paste0("'", parameter,
                "' is not a supported parameter. Use 'parameter = ",
                paramsString, "'"), call.=FALSE)
  }

  # Create filename according to the PWFSLSmoke naming scheme
  filename <- paste0("airnow_", parameter, "_latest45.RData")

  ws_monitor <- loadDataFile(filename, baseUrl, dataDir)

  return(ws_monitor)

}

#' @export
#'
#' @title Load data from URL or local file
#'
#' @param filename Name of the data file to be loaded.
#' @param dataUrl Remote URL directory for data files.
#' @param dataDir Local directory containing data files.
#' @return A data object.
#'
#' @description Loads pre-generated R binary files from a URL or a local
#' directory. This function is intended to be called by other \code{~_load()}
#' functions and can remove internet latencies when local versions of data are
#' available.
#'
#' For this reason, specification of \code{dataDir} always takes precedence over
#' \code{dataUrl}.

loadDataFile <- function(
  filename = NULL,
  dataUrl = NULL,
  dataDir = NULL
) {

  # ----- Validate parameters --------------------------------------------------

  MazamaCoreUtils::stopIfNull(filename)

  if ( is.null(dataUrl) && is.null(dataDir) ) {
    stop("either 'dataUrl' or 'dataDir' must be specified")
  }

  # ----- Load the data --------------------------------------------------------

  try({

    # Always check for dataDir first
    if (
      !is.null(dataDir) &&
      !is.na(dataDir) &&
      dir.exists(path.expand(dataDir))
    ) {

      # Load from a file
      filepath <- file.path(path.expand(dataDir), filename)

      result <- try({
        objectName <- load(filepath)
      }, silent = TRUE)

      if ( "try-error" %in% class(result) ) {
        stop(paste0("data file could not be loaded from: ", filepath))
      } else {
        loadedData <- get(objectName)
      }

    } else {

      # Load from a URL
      filepath <- paste0(dataUrl, '/', filename)

      # Define a 'connection' object so we can close it no matter what happens
      conn <- url(filepath)
      result <- try({
        objectName <- load(conn)
      }, silent = TRUE )
      close(conn)

      if ( "try-error" %in% class(result) ) {
        stop(paste0("data file could not be loaded from: ", filepath))
      } else {
        loadedData <- get(objectName)
      }

    }

  }, silent = TRUE) %>%
    stopOnError(paste0("data file could not be loaded from: ", filepath))

  # ----- Return ---------------------------------------------------------------

  return(loadedData)

}

#' @name stopOnError
#' @export
#' @title Error message generator
#' @param result Return from a \code{try()} block.
#' @param err_msg Custom error message.
#' @param prefix Text string to add in front of the error message.
#' @param maxLength Maximum length of an error message. Error messages
#' beyond this limit will be truncated.
#' @param truncatedLength Length of the output error message.
#' @param call. Logical indicating whether the call should become part of the error message.
#'
#' @return Issues a \code{stop()} with an appropriate error message.
#'
#' @description When writing R code for use in production systems, it is
#' important to enclose chunks of code inside of \code{try()} blocks. This is
#' especially important when processing user input or data obtained from web
#' services which may fail for a variety of reasons. If any problems arise
#' within a \code{try()} block, it is important to generate informative and
#' consistent error messages.
#'
#' Over the years, we have developed our own standard protocol for error handling
#' that is easy to understand, easy to implement, and allows for consistent
#' generation of error messages. To goal is to make it easy for developers to test
#' sections of code that might fail and to create more uniform, more informative
#' error messages than those that might come from deep within the \R execution stack.
#'
#' In addition to the generation of custom error messages, use of \code{prefix}
#' allows for the creation of classes of errors that can be detected and handled
#' appropriately as errors propagate to other functions.
#'
#' @note If logging has been initialized, the customized/modified error message
#' will be logged with \code{logger.error(err_msg)} before issuing
#' \code{stop(err_msg)}.
#'
#' The following examples show how to use this function:
#'
#' \preformatted{
#' library(MazamaCoreUtils)
#'
#' # Arbitrarily deep in the stack we might have:
#'
#' myFunc <- function(x) {
#'   a <- log(x)
#' }
#'
#'
#' # Simple usage
#'
#' userInput <- 10
#' result <- try({
#'   myFunc(x = userInput)
#' }, silent = TRUE)
#' stopOnError(result)
#'
#' userInput <- "ten"
#' result <- try({
#'   myFunc(x = userInput)
#' }, silent = TRUE)
#' stopOnError(result)
#'
#'
#' # More concise code with the '\%>\%' operator
#'
#' try({
#'   myFunc(x = userInput)
#' }, silent = TRUE) \%>\%
#' stopOnError(err_msg = "Unable to process user input")
#'
#' try({
#'   myFunc(x = userInput)
#' }, silent = TRUE) \%>\%
#' stopOnError(prefix = "USER_INPUT_ERROR")
#'
#'
#' # Truncating error message length
#'
#' try({
#'   myFunc(x = userInput)
#' }, silent = TRUE) \%>\%
#' stopOnError(
#'   prefix = "USER_INPUT_ERROR",
#'   maxLength = 40,
#'   truncatedLength = 32
#' )
#'
#' }

stopOnError <- function(
  result,
  err_msg = "",
  prefix = "",
  maxLength = 500,
  truncatedLength = 120,
  call. = FALSE
) {

  if ( "try-error" %in% class(result) ) {

    # Use passed in message or cleaned up version from geterrmessage()
    err_msg <- ifelse(err_msg == "", geterrmessage(), err_msg)

    err_msg <-
      err_msg %>%
      stringr::str_replace("Error : ", "") %>%
      stringr::str_replace("Error: ", "") %>%
      stringr::str_trim()

    if ( prefix != "" )
      err_msg <- paste(prefix, err_msg)

    if ( stringr::str_length(err_msg) > maxLength )
      err_msg <- paste(stringr::str_sub(err_msg, end = truncatedLength), "...")

    # if ( logger.isInitialized() )
    #   logger.error(err_msg)

    stop(err_msg, call. = call.)

  }

}

#' @keywords AIRSIS
#' @export
#' @title Load annual AIRSIS monitoring data
#' @param year Desired year (integer or character representing YYYY).
#' @param parameter Parameter of interest.
#' @param baseUrl Base URL for 'annual' AIRSIS data files.
#' @param dataDir Local directory containing 'annual' data files.
#' @return A \emph{ws_monitor} object with AIRSIS data.
#' @description Loads pre-generated .RData files containing annual
#' AIRSIS data.
#'
#' If \code{dataDir} is defined, data will be loaded from this local
#' dirctory. Otherwise, data will be loaded from the monitoring data repository
#' maintained by PWFSL.
#'
#' The annual files loaded by this function are updated on the 15'th of each
#' month and cover the period from the beginning of the year to the end of the
#' last month.
#'
#' For data during the last 45 days, use \code{airsis_loadDaily()}.
#'
#' For the most recent data, use \code{airsis_loadLatest()}.
#'
#' AIRSIS parameters include the following:
#' \enumerate{
# #' \item{BARPR}
# #' \item{BC}
# #' \item{CO}
# #' \item{NO}
# #' \item{NO2}
# #' \item{NO2Y}
# #' \item{NO2X}
# #' \item{NOX}
# #' \item{NOOY}
# #' \item{OC}
# #' \item{OZONE}
# #' \item{PM10}
#' \item{PM2.5}
# #' \item{PRECIP}
# #' \item{RHUM}
# #' \item{SO2}
# #' \item{SRAD}
# #' \item{TEMP}
# #' \item{UV-AETH}
# #' \item{WD}
# #' \item{WS}
#' }
#'
#' Available AIRSIS RData and associated log files can be seen at:
#' \href{https://haze.airfire.org/monitoring/AIRSIS/RData/}{https://haze.airfire.org/monitoring/AIRSIS/RData/}
#' @seealso \code{\link{airsis_loadDaily}}
#' @seealso \code{\link{airsis_loadLatest}}
#' @examples
#' \dontrun{
#' # Fail gracefully if any resources are not available
#' try({
#'
#' airsis_loadAnnual(2017) 
#'
#' }, silent = FALSE)
#' }

airsis_loadAnnual <- function(year = NULL,
  parameter = 'PM2.5',
  baseUrl = 'https://haze.airfire.org/monitoring',
  dataDir = NULL) {

# ----- Validate parameters --------------------------------------------------

if ( is.null(year) ) {
stop("Required parameter 'year' is missing.")
}
year <- as.numeric(year)

validParams <- c("PM2.5")
if ( !parameter %in% validParams ) {
paramsString <- paste(validParams, collapse=", ")
stop(paste0("'", parameter,
"' is not a supported parameter. Use 'parameter = ",
paramsString, "'"), call.=FALSE)
}

if ( year < 2004 ) {
stop("PWFSL has no annual data files for AIRSIS prior to 2004.")
}

# Load the data --------------------------------------------------------------


# Create URL and filename according to the PWFSLSmoke naming scheme
baseUrl <- paste0(baseUrl, '/AIRSIS/RData/', year)
filename <- paste0("airsis_", parameter, "_", year, ".RData")

ws_monitor <- loadDataFile(filename, baseUrl, dataDir)

return(ws_monitor)

}

#' @keywords AIRSIS
#' @export
#' @title Load recent AIRSIS monitoring data
#' @param parameter Parameter of interest.
#' @param baseUrl Base URL for 'daily' AirNow data files.
#' @param dataDir Local directory containing 'daily' data files.
#' @return A \emph{ws_monitor} object with AIRSIS data.
#' @description Loads pre-generated .RData files containing recent
#' AIRSIS data.
#'
#' If \code{dataDir} is defined, data will be loaded from this local
#' dirctory. Otherwise, data will be loaded from the monitoring data repository
#' maintained by PWFSL.
#'
#' The daily files loaded by this function are updated once a day, shortly
#' after midnight and contain data for the previous 45 days.
#'
#' For the most recent data, use \code{airsis_loadLatest()}.
#'
#' For data extended more than 45 days into the past, use \code{airsis_loadAnnual()}.
#'
#' AIRSIS parameters include the following:
#' \enumerate{
# #' \item{BARPR}
# #' \item{BC}
# #' \item{CO}
# #' \item{NO}
# #' \item{NO2}
# #' \item{NO2Y}
# #' \item{NO2X}
# #' \item{NOX}
# #' \item{NOOY}
# #' \item{OC}
# #' \item{OZONE}
# #' \item{PM10}
#' \item{PM2.5}
# #' \item{PRECIP}
# #' \item{RHUM}
# #' \item{SO2}
# #' \item{SRAD}
# #' \item{TEMP}
# #' \item{UV-AETH}
# #' \item{WD}
# #' \item{WS}
#' }
#'
#' Avaialble AIRSIS RData and associated log files can be seen at:
#' \href{https://haze.airfire.org/monitoring/AIRSIS/RData/latest/}{https://haze.airfire.org/monitoring/AIRSIS/RData/latest/}
#' @seealso \code{\link{airsis_loadAnnual}}
#' @seealso \code{\link{airsis_loadLatest}}
#' @examples
#' \dontrun{
#' # Fail gracefully if any resources are not available
#' try({
#'
#' airsis_loadDaily()
#'
#' }, silent = FALSE)
#' }

airsis_loadDaily <- function(parameter = 'PM2.5',
                             baseUrl = 'https://haze.airfire.org/monitoring/latest/RData',
                             dataDir = NULL) {

  # Validate parameter
  validParams <- c("PM2.5")
  if ( !parameter %in% validParams ) {
    paramsString <- paste(validParams, collapse=", ")
    stop(paste0("'", parameter,
                "' is not a supported parameter. Use 'parameter = ",
                paramsString, "'"), call.=FALSE)
  }

  # Create filename according to the PWFSLSmoke naming scheme
  filename <- paste0("airsis_", parameter, "_latest45.RData")

  ws_monitor <- loadDataFile(filename, baseUrl, dataDir)

  return(ws_monitor)

}