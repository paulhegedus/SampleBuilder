#' @title Makes a grid across a polygon boundary
#'
#' @description Takes a path to a polygon layer and creates a
#' grid with cells of a user specified dimension (in meters)
#' across the area. Returns an 'sf' class object.
#'
#' @param poly_layer_path Character, path to the .shp file of the polygon layer
#' of the area of interest.
#' @param cell_size Numeric, the size of the grid cell in meters.
#'
#' @return An 'sf' class grid object
#' @export
make_grid <- function(poly_layer_path,
                      cell_size) {
  ## get the polygon data & put in utm coords
  poly_dat <- get_data(poly_layer_path)
  ## make the grid using meters
  cell_grid <- sf::st_make_grid(poly_dat, cellsize = 61) %>%
    sf::st_as_sf() %>%
    rename_geometry("geometry")
  cell_grid$id <- seq(1:nrow(cell_grid))
  ## put in lat long
  cell_grid_out <- sf::st_transform(cell_grid, 4326)
  return(cell_grid_out)
}

#' @title Get centroids of a grid
#'
#' @description Takes a grid object and returns the centroids of each cell.
#'
#' @param grid_dat An 'sf' object with grid cells for which to get centroid for.
#'
#' @return An 'sf' class grid object
#' @export
get_centroids <- function(grid_dat) {
  utm_epsg <- calcUTMzone(grid_dat)

  ## project if not
  tryCatch({
    grid_dat <- sf::st_transform(grid_dat, utm_epsg)
  },
  warning = function(w) {},
  error = function(e) {
    grid_dat <- sf::st_set_crs(grid_dat, 4326) %>%
      sf::st_transform(utm_epsg)
  })

  ## get centroids
  cell_centroids <- suppressWarnings(sf::st_centroid(grid_dat))

  ## put in lat long
  cell_grid_out <- sf::st_transform(cell_centroids, 4326)

  ## put lat and long coords in file
  temp <- as(cell_grid_out, "Spatial") %>%
    sp::coordinates()
  cell_grid_out <- cbind(cell_grid_out, temp) %>%
    `names<-`(c("id", "long", "lat", "geometry"))

  return(cell_grid_out)
}

rename_geometry <- function(g, name){
  current <- attr(g, "sf_column")
  names(g)[names(g) == current] <- name
  sf::st_geometry(g) <- name
  return(g)
}







