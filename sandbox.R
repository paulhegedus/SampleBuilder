## sandbox

library(magrittr)
library(methods)
library(sf)
library(geosphere)
library(sp)
library(raster)
library(ggplot2)

source("R/gen_sampling_funs.R")
source("R/grid_sampling_funs.R")
source("R/make_transect_funs.R")

data_path <- "/Users/PaulBriggs/Library/CloudStorage/OneDrive-MontanaStateUniversity/BoxMigratedData/Hegedus/Technician/Misc/ChiefMtn/"

line_layer_path <- paste0(data_path, "ChiefMtn_Allotment_roads.shp")
t_number <- 100
t_length <- 1000
t_size <- 10
buddy_t = TRUE
direction = c("positive", "negative")
allow_overlaps = FALSE
stratify_on_length = TRUE

line_layer <- get_data(line_layer_path)

if (buddy_t & (t_number %% 2) != 0) {
  t_number <- t_number + 1
}

line_layer$segment <- 1:nrow(line_layer)
line_layer$n <- get_n(line_layer,
                      t_number,
                      direction,
                      buddy_t,
                      stratify_on_length)

t_lines <- rep(list(NA), nrow(line_layer))


i = 1
line_shp <- line_layer[i, ]

## Make Transect Offsets (in both directions)
# make coordinates
coords <- sp::coordinates(as(line_shp, "Spatial"))[[1]] %>%
  lapply(as.data.frame) %>%
  # do.call(rbind, .) %>%
  lapply(function(x) {x %>% `names<-`(c("x", "y"))})
utm_epsg <- calcUTMzone(line_shp)
# make line offsets
t_sp <- lapply(coords,
               make_line_offsets,
               direction,
               t_length,
               utm_epsg) %>%
  fix_t_sp(direction) %>%
  lapply(sf::st_as_sf)

temp_buff <- sf::st_buffer(line_shp, t_length-1)
for (j in 1:length(t_sp)) {
  t_sp[[j]] <- suppressWarnings(sf::st_difference(t_sp[[j]], temp_buff))
}

n <- line_layer$n[i]

buff_pts <- mapply(get_samp_pts, t_sp, n, utm_epsg)
t_pts <- list(
  buff_pts = lapply(buff_pts, sf::st_as_sf) %>%
    do.call(rbind, .)
)

if (buddy_t) {
  t_pts$buddy_pts <- mapply(
    make_buddies,
    buff_pts,
    t_sp,
    MoreArgs = list(t_size = t_size,
                    utm_epsg = utm_epsg)
  ) %>%
    lapply(sf::st_as_sf) %>%
    do.call(rbind, .)
}

if("buddy_pts" %in% names(t_pts)) {
  bad_buddies <- which(sf::st_is_empty(t_pts$buddy_pts))
  if(length(bad_buddies) > 0) {
    t_pts$buff_pts <- t_pts$buff_pts[-bad_buddies, "x"]
    t_pts$buddy_pts <- t_pts$buddy_pts[-bad_buddies, "x"]
  }
}

## Connect Transects & combine
if(length(t_pts) > 1) {
  if(nrow(t_pts$buddy_pts) != 0 & nrow(t_pts$buff_pts) != 0) {
    t_lines[[i]] <- lapply(t_pts, make_transect_lines, line_shp) %>%
      do.call(rbind, .)
  }
}

bad_buddies <- which(sf::st_is_empty(t_pts$buddy_pts))
if(length(bad_buddies) > 0) {
  t_pts$buff_pts <- t_pts$buff_pts[-bad_buddies, "x"]
  t_pts$buddy_pts <- t_pts$buddy_pts[-bad_buddies, "x"]
}

## Connect Transects & combine
t_lines[[i]] <- lapply(t_pts, make_transect_lines, line_shp) %>%
  do.call(rbind, .)


ggplot() +
  geom_sf(data=line_shp) +
  geom_sf(data=t_sp$pos_t_sp, aes(color='red')) +
  geom_sf(data=t_sp$neg_t_sp, aes(color='red'))

temp_buff <- sf::st_buffer(line_shp, t_length-1)

ggplot() +
  geom_sf(data=temp_buff, color='green') +
  geom_sf(data=line_shp) +
  geom_sf(data=t_sp$pos_t_sp, color='red')

new_offset <- sf::st_difference(t_sp$pos_t_sp, temp_buff)

ggplot() +
  geom_sf(data=temp_buff, color='green') +
  geom_sf(data=line_shp) +
  geom_sf(data=new_offset, color='red')

