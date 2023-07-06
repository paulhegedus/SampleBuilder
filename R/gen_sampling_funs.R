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
    crs <- sf::st_crs(FILE)
    if (is.null(crs) || is.na(crs$proj4string)) {
      crs <- sf::st_crs(FILE, 4326)
    }
    FILE <- sf::st_transform(FILE, 4326)
    # get boundary of file
    bounds <- sf::st_bbox(FILE)
    # calculate zone
    utm_zone <- ceiling((bounds["xmin"] + 180) / 6)
    # check hemisphere
    if (bounds["ymin"] > 0 | bounds["ymax"] > 0) {
      utm_epsg <- paste0("326", utm_zone)
    } else {
      utm_epsg <- paste0("327", utm_zone)
    }
    return(as.integer(utm_epsg))
  }
  return(NULL)
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

#' @title Add sampling info to layers
#'
#' @description Add attributes and fields to layers that you will
#' want to collect in the field. Pass in a data.frame with each row
#' a unique attribute field with the name spacified (e.g. comment), the type of
#' data (e.g. text), and the default value (e.g. foobar).
#'
#' @param layer An 'sf' class object to add information to
#' @param info A data.frame with each row representing a new attribute field.
#' Must have columns called; 'name', 'type', 'default'. Put the name of the
#' attribute field in the rows of the 'name' column, the type of data in the
#' 'type' column (text, integer, real, or logical), and the desired
#' default value in the 'default' column.
#' @return A 'sf' object with the added information added.
#' @export
add_info <- function(layer, info) {
  stopifnot(is.data.frame(info),
            any(grepl("name", names(info))),
            any(grepl("type", names(info))),
            any(grepl("default", names(info)))
            )
  # add fields
  temp <- matrix(NA, nrow = nrow(layer), ncol = nrow(info)) %>%
    as.data.frame()
  names(temp) <- info[, "name"]

  for (i in 1:ncol(temp)) {
    if (info[i, "type"] == "text") {
      temp[, i] <- as.character(info[i, "default"])
    }
    if (info[i, "type"] == "integer") {
      temp[, i] <- as.integer(info[i, "default"])
    }
    if (info[i, "type"] == "real") {
      temp[, i] <- as.numeric(info[i, "default"])
    }
    if (info[i, "type"] == "logical") {
      temp[, i] <- as.logical(info[i, "default"])
    }
  }

  # return points
  layer_update <- cbind(layer, temp)

  return(layer_update)
}



