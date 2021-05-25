#' @title Calculate the Universal Transverse Mercator coordinate system zone
#'
#' @description Calculates the UTM zone for the specified data and
#' returns a UTM code for WGS84 (i.e. 32612, 32614). This function is
#' performed in the ManageParms class to store the UTM zone for the farm.
#' A shapefile is passed in with the 'FILE' argument, the spatial data is
#' used to determine the UTM zone mathematically, same if the 'bounds'
#' argument is used. The returned value is the UTM EPSG code for the WGS84 datum.
#'
#' @param FILE Spatial shapefile for determining UTM zone.
#' @return WGS84 UTM code (i.e. 32612, 32614).
#' @source https://apollomapping.com/blog/gtm-finding-a-utm-zone-number-easily
#' @export
calcUTMzone = function(FILE) {
  if (!is.null(FILE)) {
    crs_raster <- suppressWarnings(raster::crs(FILE))
    if (is.na(crs_raster)) {
      sf::st_crs(FILE) <- 4326
    }
    FILE <- sf::st_transform(FILE, 4326)
    # get boundary of file
    bounds <- sf::st_bbox(FILE)
    # calculate zone
    utm_zone <- ceiling((bounds["xmin"] + 180) / 6) %>% as.numeric()
    # check hemisphere
    if (bounds["ymin"] > 0 | bounds["ymax"] > 0) {
      utm_epsg <- paste0(326, utm_zone)
    } else {
      utm_epsg <- paste0(327, utm_zone)
    }
  }
  return(as.integer(utm_epsg))
}
#' @title Import data from path and convert to UTM
#'
#' @description Imports data from specified path. Calculates the UTM
#' code for WGS84 (i.e. 32612, 32614) for the specified data and projects
#' the data. Returns the data projected into the correct UTM zone.
#'
#' @param data_path Path to data for import.
#' @return Data projected in UTM coordinates
#' @export
get_data <- function(data_path) {
  x <- sf::st_read(data_path)
  utm_epsg <- calcUTMzone(x)
  x <- sf::st_transform(x, crs = sf::st_crs(utm_epsg))
  return(x)
}
