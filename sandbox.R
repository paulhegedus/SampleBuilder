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
poly_centroids <- sf::st_centroid(poly_layer$geometry) %>%
  sf::st_as_sf()
sf::st_geometry(poly_centroids) <- "geometry"
poly_centroids$poly_id <- poly_layer$id
poly_centroids$poly_area <- sf::st_area(poly_layer$geometry) %>% as.numeric()

dists <- sf::st_distance(line_layer, poly_centroids)
poly_centroids$id <- NA
for (i in 1:nrow(poly_centroids)) {
  poly_centroids$id[i] <- line_layer$id[which.min(dists[, i])]
}

lines_close_to_poly <- poly_centroids[, c("poly_area", "id")] %>%
  group_by(id) %>%
  summarize(n_close_polys = n(),
            sum_poly_area =sum(poly_area)) %>%
  sf::st_drop_geometry()
lines_close_to_poly$n_x_area <- lines_close_to_poly$n_close_polys * lines_close_to_poly$sum_poly_area
sub_line_layer <- merge(line_layer, lines_close_to_poly, by="id") # merge(line_layer, temp, by="id", all.x=TRUE)

weights <- sub_line_layer$n_x_area / sum(sub_line_layer$n_x_area)
sub_line_layer$n <- round(weights * t_number)
sub_line_layer$n <- pmax(sub_line_layer$n, 1)
diff_n <- t_number - sum(sub_line_layer$n)
if (diff_n < 0) {
  for (i in 1:nrow(sub_line_layer)) {
    if (sub_line_layer$n[i] > 1) {
      sub_line_layer$n[i] <- sub_line_layer$n[i] - 1
      diff_n <- diff_n + 1
      if (diff_n == 0) {
        break
      }
    }
  }
}



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


