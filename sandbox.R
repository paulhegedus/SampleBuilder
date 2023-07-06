## sandbox

library(magrittr)
library(sf)
library(ggplot2)

##@@@@@@@@@@@ Functions @@@@@@@@@@@@@@@@@
source("R/gen_sampling_funs.R")
source("R/make_transect_funs.R")

##@@@@@@@@@@@ Parameters @@@@@@@@@@@@@@@@@
## Path
data_path <- "/Users/PaulBriggs/Library/CloudStorage/OneDrive-MontanaStateUniversity/BoxMigratedData/Hegedus/Technician/Misc/ChiefMtn/"
## Specify the name of the boundary shapefile.
boundary_filename <- "ChiefMtnAllotments_Shapefile.shp"
## Specify the name of the line (e.g. road) layer.
line_layer_path_filename <- "ChiefMtn_Allotment_open_roads.shp"
## Specify the name of the polygon layer.
poly_layer_path_filename <- "test_bison_layers.shp"
## Specify the filename to use when saving the transect lines.
transect_filename <- "ChiefMtn_Allotment_open_transects.shp"
## Specify the filename to use when saving the transect points.
transect_points_filename <- "ChiefMtn_Allotment_open_transect_points.shp"


boundary <- sf::st_read(paste0(data_path, "/", boundary_filename)) %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(4326)
line_layer_path <- paste0(data_path, "/", line_layer_path_filename)
poly_layer_path <- paste0(data_path, "/", poly_layer_path_filename)

##@@@@@@@@@@@ make_transects() @@@@@@@@@@@@@@@@@
line_layer_path = line_layer_path
poly_layer_path = poly_layer_path
poly_strat_col = 'strat'
t_number = 20
t_length = 1000
t_size = 10
buddy_t = TRUE
direction = c("positive", "negative")
allow_overlaps = FALSE


## Import line and poly layers (puts in UTM)
line_layer <- get_data(line_layer_path)
if (!is.null(poly_layer_path)) {
  poly_layer <- get_data(poly_layer_path)
}

## Clip line layer to avoid overlapping transects
if (!allow_overlaps) {
  line_layer <- clip_to_valid_sampling_area(line_layer, t_length)
}

## If buddy_t cut t_number in half
if (buddy_t) {
  t_number <- ceiling(t_number / 2)
}

## get_n_strat
if (is.null(poly_layer_path)) {
  line_layer <- get_n(line_layer, t_number)
} else {
  line_layer <- get_strat_n(line_layer, t_number, poly_layer, poly_strat_col)
}

## For each line layer, sample based on n
epsg <- calcUTMzone(line_layer)
sample_result <- sample_line(line_layer, t_size, epsg, buddy_t)
samples <- do.call(rbind, sample_result)

## Create transects
transect_lines <- create_transects(samples, line_layer, t_length, epsg, direction)




transects <- make_transects(
  line_layer_path = line_layer_path,
  poly_layer_path = poly_layer_path,
  poly_strat_col = poly_strat_col,
  t_number = 100,
  t_length = 500,
  t_size = 10,
  buddy_t = TRUE,
  direction = c("positive", "negative"),
  allow_overlaps = FALSE
)


ggplot() +
  # geom_sf(data = boundary) +
  geom_sf(data = get_data(line_layer_path)) +
  geom_sf(data = transects, color="red")


transect_points <- make_transect_pts(
  t_lines = transects,
  t_length = 500,
  t_size = 10
)

sampling_info <- data.frame(
  name = c("date", "comment1", "comment2", "comment3"),
  type = c("character", "character", "character", "character"),
  default = rep(NA, 4)
)
## Code for adding the information to the point layer.
transect_points_wInfo <- add_info(
  transect_points,               ## Shapefile, transect points.
  sampling_info                  ## Dataframe, sampling information.
)


