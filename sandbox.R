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
poly_layer_path_filename <- "bison_habitat_layer.shp"
## Specify the filename to use when saving the transect lines.
transect_filename <- "ChiefMtn_Allotment_open_transects.shp"
## Specify the filename to use when saving the transect points.
transect_points_filename <- "ChiefMtn_Allotment_open_transect_points.shp"


boundary <- sf::st_read(paste0(data_path, "/", boundary_filename)) %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(4326)
line_layer_path <- paste0(data_path, "/", line_layer_path_filename)
poly_layer_path <- paste0(data_path, "/", poly_layer_path_filename)
origin_layer_path <- poly_layer_path
##@@@@@@@@@@@ make_transects() @@@@@@@@@@@@@@@@@
sample_layer_path = poly_layer_path ## path to layer to sample from (polygon or point)
poly_strat_col = 'strat'
t_number = 20
spokes = 8
t_length = 300
t_size = 10
buddy_t = TRUE
direction = c("positive", "negative")
allow_overlaps = FALSE


wagonwheels <- make_wagonwheel_transects(origin_layer_path,
                                         t_length,
                                         spokes,
                                         t_size)
ggplot() +
  # geom_sf(data = boundary) +
  geom_sf(data = get_data(poly_layer_path)) +
  geom_sf(data = wagonwheels, color="red")






