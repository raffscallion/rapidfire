# Just commenting out this old code for now because I may need it to support data I
# delivered using the older grids

# # Tools to create an output grid for predection and pretty pictures. We can also
# # just do the predictions at participant locations, but that requires
# # running kriging and extraction every time a location changed.
# 
# # Make a grid for output predictions from an SPDF
# make_grid <- function(spdf, cellsize) {
#   box <- spdf@bbox
#   xrange <- box[3] - box[1]
#   yrange <- box[4] - box[2]
#   xcells <- ceiling(xrange / cellsize)
#   ycells <- ceiling(yrange / cellsize)
# 
#   x <- seq(1, xcells, by = 1) * cellsize + box[1]
#   y <- seq(1, ycells, by = 1) * cellsize + box[2]
# 
#   df <- data.frame(Id = seq(1, length(x) * length(y), by = 1))
# 
#   grid <- sp::SpatialPointsDataFrame(cbind(rep(x, length(y)), rep(y, each = length(x))),
#                         proj4string = sp::CRS("+init=epsg:3395"), data = df)
#   sp::gridded(grid) <- TRUE
#   return(grid)
# }

#' Creates a grid of the given cellsize (in m) that completely covers the input spatial
#' vector. It uses the standard rapidfire coordinate system of epsg:3395.
#'
#'
#' @param spat_vec A spatial vector that defines the domain for the output grid
#' @param filename The filename for the output grid 
#' @param cellsize The cell size in meters of the output grid (default = 1000 m, 1 km)
#'
#' @returns This function is used for the side-effect of writing the grid to a file
#' @export
#'
#' @examples 
#' \dontrun{
#' grid_state(california_shapefile, "output_grid_1km", cellsize = 1000)
grid_state <- function(spat_vec, filename, cellsize = 1000) {
  
  outline <- terra::project(spat_vec, "epsg:3395")
  box <- terra::ext(outline)
  
  r <- terra::rast(box, resolution = cellsize, crs = terra::crs(outline))
  # initialize all values as NA
  r <- terra::init(r, NA)
  terra::writeRaster(r, filename)
}


# 
# # Make grid for output that covers an SPDF (such as a state or states)
# make_state_grid <- function(outline, cellsize) {
# 
#   # Ensure we are in Mercator
#   outline <- sp::spTransform(outline, sp::CRS("+init=epsg:3395"))
# 
#   # Get the bounding box
#   box <- outline@bbox
# 
#   # create grid cells
#   xrange <- box[3] - box[1]
#   yrange <- box[4] - box[2]
#   xcells <- ceiling(xrange / cellsize)
#   ycells <- ceiling(yrange / cellsize)
# 
#   x <- seq(1, xcells, by = 1) * cellsize + box[1]
#   y <- seq(1, ycells, by = 1) * cellsize + box[2]
# 
#   df <- data.frame(Id = seq(1, length(x) * length(y), by = 1))
# 
#   grid <- sp::SpatialPointsDataFrame(cbind(rep(x, length(y)), rep(y, each = length(x))),
#                                      proj4string = sp::CRS("+init=epsg:3395"), data = df)
# 
#   # filter only to grid cells that fall within state (at least partially)
#   stop("This function is currently inactive while rgeos is replaced")
#   #g <- rgeos::gIntersection(grid, outline, byid = c(TRUE, FALSE))
# 
#   # Convert to spatial pixels
#   sp::gridded(g) <- TRUE
#   return(g)
# 
# }
# 
# # Get ids and lat/lons (in degrees) as a data frame from a grid for making assignements
# grid_to_latlon <- function(grid) {
# 
#   df <- tibble::tibble(monitorID = grid@grid.index)
#   deg <- sp::spTransform(grid, sp::CRS("+init=epsg:4326"))
#   df$longitude <- deg@coords[,1]
#   df$latitude <- deg@coords[,2]
#   return(df)
# 
# }
