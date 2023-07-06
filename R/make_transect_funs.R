#' @title Make transect lines for sampling
#'
#' @description For returning a specified number of transects perpendicular
#' to a user supplied line feature. User options include the number of transects,
#' the length of the transect (in meters), and the size (assumed to be a square) of
#' each transect cell (in meters).
#'
#' Set 'buddy_t' to TRUE to generate a transect less 3 transect sizes away but at
#' least 1 transect size away so that sampling can be performed traveling away from
#' the line feature and back to the line feature.
#'
#' The direction specifies which direction from the line layer transects should be
#' generated. Typically "positive" indicates transects to the North and West while
#' "negative" generates transects to the South and East. Default is both directions
#' from the line layer. Recommended to play with and decide for your needs.
#'
#' While all overlaps are not guaranteed to be removed, setting 'allow_overlaps' to
#' FALSE (default) will restrict sampling to where the risk is minimized. This is done
#' by making a buffer the width of the transects around the line layer and removing
#' sections of the line layer where the buffer overlaps.
#'
#' Note, if the number of transects are smaller than the number of distinct line features
#' then transects will be placed on random line features. If the number of transects is
#' equal to the number of distinct line features then one transect will be placed on each
#' line feature. If the number of transects is greater than the number of line features
#' then at least one transect will be placed on each line feature and remaining transects
#' will be distributed by the length of line features (e.g. longer line features have more
#' transects).
#'
#' @param line_layer_path Character, path to the .shp file of the line layer to build
#' transects from.
#' @param poly_layer_path Optional character, path to the .shp file of the polygon layer to
#' stratify transects on.
#' @param t_number Numeric, the approximate number of transects to generate.
#' @param t_length Numeric, the desired length of transects (in meters).
#' @param t_size Numeric, the size (Y x Y) of the transect cells (in meters).
#' @param buddy_t Logical, default is TRUE,  whether to include a buddy transect
#' for efficient sampling.
#' @param direction Character, default is c("positive", "negative"), the direction from
#' the line layer transects should be built.
#' @param allow_overlaps Logical, default is FALSE, whether to remove areas of line where
#' transects could overlap with other transects.
#' @return A 'sf' object with transect lines and ID numbers in lat long.
#' @export
make_transects <- function(line_layer_path,
                           poly_layer_path = NULL,
                           poly_strat_col = NULL,
                           t_number,
                           t_length,
                           t_size,
                           buddy_t = TRUE,
                           direction = c("positive", "negative"),
                           allow_overlaps = FALSE) {
  # browser()

  ## Parameter checks
  stopifnot(is.character(line_layer_path),
            # is.character(poly_layer_path),
            # is.character(poly_strat_col),
            is.numeric(t_number),
            is.numeric(t_length),
            is.numeric(t_size),
            is.logical(buddy_t),
            is.character(direction),
            !any(!grepl("positive|negative", direction)),
            is.logical(allow_overlaps))

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

  ## Determine number of transects on each line feature
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

  ## Return transect lines
  return(transect_lines)
}

#' @title Make points at center of transect cells
#'
#' @description Lays points out at the center of sampling cells along a transect.
#' User supplies a layer with transect lines, the length of each transect, and the
#' size (assumed to be a square) of the transect cells. Returns a point layer with
#' a transect ID and distance along transect to form a unique point ID.
#'
#' @param t_lines An 'sf' object with transect lines.
#' @param t_length The length of the transects in 't_layer' (in meters).
#' @param t_size The size of the transect cells (in meters), assumed to be
#' square (t_size x t_size).
#' @return An 'sf' object with sampling points along transect lines.
#' @export
make_transect_pts <- function(t_lines,
                              t_length,
                              t_size) {
  utm_epsg <- calcUTMzone(t_lines)
  t_lines_proj <- sf::st_transform(t_lines, utm_epsg)
  point_list <- split(t_lines_proj, seq(nrow(t_lines_proj)))

  # for each transect line
  for (i in 1:length(point_list)) {
    point_list[[i]] <- sf::st_sample(point_list[[i]],
                                     (t_length / (t_size / 2)),
                                     type = "regular") %>%
      sf::st_as_sf() %>%
      sf::st_set_crs(utm_epsg)
    point_list[[i]] <- sf::st_cast(point_list[[i]], 'POINT')
    sf::st_geometry(point_list[[i]]) <- "geometry"

    point_list[[i]] <- point_list[[i]][-seq(1, nrow(point_list[[i]]), 2), ]
    point_list[[i]]$t_id <- rep(t_lines_proj[i, "id"]$id,
                                nrow(point_list[[i]])) %>%
      as.character()
    point_list[[i]]$dist <- seq(t_size,
                                t_length,
                                t_size) %>%
      as.character()
    point_list[[i]]$id <- paste(point_list[[i]]$t_id, point_list[[i]]$dist, sep = ".")
  }
  points_sdf <- do.call(rbind, point_list) %>%
    sf::st_transform(4326)

  return(points_sdf)
}


## clip line layer to a valid sampling area that excludes zones where
## transects coming off of lines could intersect other transects
clip_to_valid_sampling_area <- function(line_layer,
                                        t_length) {
  # make buffer of length of transect around line layer
  buff_layer <- sf::st_buffer(line_layer, t_length)

  # remove overlap to leave valid sampling area
  intersect_list <- suppressWarnings(sf::st_intersects(buff_layer)) %>%
    as.list()
  new_buff_layer <- split(buff_layer, seq(nrow(buff_layer)))
  for (i in 1:nrow(intersect_list)) {
    overlaps <- intersect_list[[i]]
    for (j in 1:length(overlaps)) {
      if (i != overlaps[j]) {
        new_buff_layer[[i]] <- suppressWarnings(
          sf::st_difference(new_buff_layer[[i]],
                            buff_layer[overlaps[j], ])
        )
      }
    }
  }
  one_buff_layer <- lapply(new_buff_layer, function(df) df[, "geometry"]) %>%
    do.call(rbind, .)

  # clip line area by valid sampling area to remove portions of lines
  # where a transect could overlap with another
  new_line_layer <- suppressWarnings(sf::st_intersection(line_layer, one_buff_layer))

  return(new_line_layer)
}


## Get the number of transects for each line feature
get_n <- function(line_layer,
                  t_number) {
  line_layer$len_m <- sf::st_length(line_layer) %>% as.numeric()
  if (nrow(line_layer) == t_number) {
    line_layer$n <- 1
  }
  if (nrow(line_layer) > t_number) {
    line_layer$n <- sample(c(rep(0, nrow(line_layer) - t_number), rep(1, t_number)))
  }
  if (nrow(line_layer) < t_number) {
    weights <- line_layer$len_m / sum(line_layer$len_m)
    line_layer$n <- round(weights * t_number)
    line_layer$n <- pmax(line_layer$n, 1)
    diff_n <- t_number - sum(line_layer$n)
    if (diff_n < 0) {
      for (i in 1:nrow(line_layer)) {
        if (line_layer$n[i] > 1) {
          line_layer$n[i] <- line_layer$n[i] - 1
          diff_n <- diff_n + 1
          if (diff_n == 0) {
            break
          }
        }
      }
    }
  }

  return(line_layer)
}

## Get the number of transects for each line feature
## stratified on how many polygons are closest to
## each line layer.
##
## for example, if a line feature is closest to 8
## polygon features it will get more transects than
## a line layer closest to 1 polygon feature.
##
## if a line feature is not the closest to any
## polygon feature it is not sampled
get_strat_n <- function(line_layer,
                        t_number,
                        poly_layer,
                        strat_col) {
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

  return(sub_line_layer)
}


## sample transect start points on lines
sample_line <- function(line_layer, t_size, epsg, buddy_t) {
  for (i in 1:nrow(line_layer)) {
    samp_pts <- sf::st_sample(
      line_layer$geometry[i],
      line_layer$n[i],
      type = "regular" # "random"
    ) %>%
      sf::st_as_sf()
    is_empty <- sf::st_is_empty(samp_pts)
    samp_pts <- samp_pts[!is_empty, ] %>%
      sf::st_cast('POINT')
    samp_pts <- sf::st_cast(samp_pts, 'POINT')
    sf::st_geometry(samp_pts) <- "geometry"
    # ggplot() + geom_sf(data=line_layer$geometry[i]) + geom_sf(data=samp_pts$geometry)

    # Perform spatial filter to remove points within t_sizes distance
    # min_dist <- t_size * 3
    # repeat {
    #   # Check distance between points and remove any points within t_sizes distance
    #   keep_pts <- vector("logical", nrow(samp_pts))
    #   for (j in 1:nrow(samp_pts)) {
    #     dist <- sf::st_distance(samp_pts[j, ], samp_pts) %>% as.numeric()
    #     dist <- dist[-j]
    #     if (all(dist > min_dist)) {
    #       keep_pts[j] <- TRUE
    #     }
    #   }
    #   samp_pts <- samp_pts[keep_pts, ]
    #   # If the number of points is less than line_layer$n[i], add additional random points
    #   num_points <- nrow(samp_pts)
    #   if (num_points != line_layer$n[i]) {
    #     additional_pts <- sf::st_sample(
    #       line_layer$geometry[i],
    #       abs(nrow(samp_pts) - num_points),
    #       type = "random"
    #     ) %>%
    #       sf::st_as_sf()
    #     is_empty <- st_is_empty(additional_pts)
    #     additional_pts <- additional_pts[!is_empty, ] %>%
    #       sf::st_cast('POINT')
    #     additional_pts <- sf::st_cast(additional_pts, 'POINT')
    #     samp_pts <- rbind(samp_pts, additional_pts)
    #   } else {
    #     break
    #   }
    # }

    if (buddy_t) {
      for (k in 1:nrow(samp_pts)) {
        outer_buff <- sf::st_buffer(samp_pts$geometry[k], t_size*3)
        inner_buff <- sf::st_buffer(samp_pts$geometry[k], t_size)
        buff <- sf::st_difference(outer_buff, inner_buff)
        # ggplot() + geom_sf(data=buff) + geom_sf(data=samp_pts$geometry[k])
        samp_line <- sf::st_intersection(line_layer$geometry[i], buff)
        # ggplot() + geom_sf(data=buff) + geom_sf(data=samp_line) + geom_sf(data=samp_pts$geometry[k])
        buff_pt <- sf::st_sample(
          samp_line,
          1,
          type = "random" # "random"
        ) %>%
          sf::st_as_sf()
        is_empty <- sf::st_is_empty(buff_pt)
        buff_pt <- buff_pt[!is_empty, ] %>%
          sf::st_cast('POINT')
        sf::st_geometry(buff_pt) <- "geometry"
        # ggplot() + geom_sf(data=buff) + geom_sf(data=buff_pt, col='red') + geom_sf(data=samp_line) + geom_sf(data=samp_pts$geometry[k])
        buff_pt$line_id <- line_layer$id[i]
        if (k == 1) {
          buddy_pts <- buff_pt
        } else {
          buddy_pts <- rbind(buddy_pts, buff_pt)
        }
      }
      if (i == 1) {
        buddy_pts_ls <- buddy_pts
      } else {
        buddy_pts_ls <- rbind(buddy_pts_ls, buddy_pts)
      }
    }
    # ggplot() + geom_sf(data=samp_pts$geometry) + geom_sf(data=buddy_pts$geometry, col='red') + geom_sf(data=line_layer$geometry[i])
    samp_pts$line_id <- line_layer$id[i]
    if (i == 1) {
      samp_pts_ls <- samp_pts
    } else {
      samp_pts_ls <- rbind(samp_pts_ls, samp_pts)
    }
  }
  samp_pts_ls$point_id <- 1:nrow(samp_pts_ls)
  if (buddy_t) {
    buddy_pts_ls$point_id <- 1:nrow(buddy_pts_ls)
    output <- list(
      'samp_pts' = samp_pts_ls,
      'buddy_pts' = buddy_pts_ls
    )
  } else {
    output <- list(
      'samp_pts' = samp_pts_ls
    )
  }

  return(output)
}


## make the transect lines for each point
create_transects <- function(samples, line_layer, t_length, epsg, direction) {
  for (i in 1:nrow(samples)) {
    buff <- sf::st_buffer(samples[i, ], t_length)
    sub_line <- line_layer[line_layer$id == buff$line_id, "geometry"] %>%
      sf::st_intersection(buff$geometry)
    # ggplot() + geom_sf(data = buff) + geom_sf(data = sub_line) + geom_sf(data = samples[i, ])

    line_coords <- sf::st_cast(sub_line, "POINT")
    dists <- sf::st_distance(samples[i, ], line_coords) %>% as.numeric()
    line_coords <- line_coords[head(order(dists), 1), ] %>%
      sf::st_coordinates()
    point_coords <- sf::st_coordinates(samples[i, ])

    angle <- atan2(line_coords[1, "Y"] - point_coords[1, "Y"], line_coords[1, "X"] - point_coords[1, "X"])
    angle_degrees <- angle * 180 / pi

    if (any(grepl("positive", direction))) {
      pos_transect <- build_transect(point_coords[1, "X"],
                                     point_coords[1, "Y"],
                                     angle_degrees,
                                     t_length,
                                     t_dir = "positive",
                                     buff$line_id,
                                     buff$point_id,
                                     epsg)
    }
    if (any(grepl("negative", direction))) {
      neg_transect <- build_transect(point_coords[1, "X"],
                                     point_coords[1, "Y"],
                                     angle_degrees,
                                     t_length,
                                     t_dir = "negative",
                                     buff$line_id,
                                     buff$point_id,
                                     epsg)
    }
    # ggplot() +
    #   geom_sf(data = buff) +
    #   geom_sf(data = sub_line) +
    #   geom_sf(data = samples[i, ]) +
    #   geom_sf(data = pos_transect, col='green') +
    #   geom_sf(data = neg_transect, col='red')
    if (any(grepl("positive", direction)) & any(grepl("negative", direction))) {
      transect <- rbind(pos_transect, neg_transect)
    } else {
      transect <- ifelse(any(grepl("positive", direction)), pos_transect, neg_transect)
    }
    if (i == 1) {
      transect_lines <- transect
    } else {
      transect_lines <- rbind(transect_lines, transect)
    }
  }
  transect_lines$id <- 1:nrow(transect_lines)

  return(transect_lines)
}


## build a transect
build_transect <- function(x, y, angle_degrees, t_length, t_dir, line_id, point_id, epsg) {
  dir_deg <- ifelse(t_dir == "positive", 90, -90)

  angle_degrees <- angle_degrees + dir_deg
  angle_radians <- angle_degrees * pi / 180
  x_end <- x + t_length * cos(angle_radians)
  y_end <- y + t_length * sin(angle_radians)

  transect <- sf::st_linestring(matrix(c(x, y,
                                         x_end, y_end),
                                       ncol = 2,
                                       byrow = TRUE)) %>%
    sf::st_cast("MULTILINESTRING") %>%
    sf::st_geometry() %>%
    sf::st_set_crs(epsg) %>%
    sf::st_as_sf()
  sf::st_geometry(transect) <- "geometry"

  transect$line_id <- line_id
  transect$point_id <- point_id
  transect$direction <- t_dir
  transect$azimuth <- angle_degrees
  transect$dist_m <- sf::st_length(transect) %>% as.numeric()

  return(transect)
}


