#' @title Make transect lines for sampling
#'
#' @description For returning a specified number of transects perpendicular
#' to a user supplied line feature. User options include the number of transects,
#' the length of the transect (in meters), and the size (assumed to be a square) of
#' each transect cell (in meters).
#'
#' Set 'buddy_t' to TRUE to generate a transect less than 3 transect sizes away but at
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
#' Note, transect locations are random, so every time this function is run you will
#' receive a new set of transects. Also note, if error occurs saying "n too small"
#' simply rerun the function a couple times.
#'
#' @param line_layer_path Character, path to the .shp file of the line layer to build
#' transects from.
#' @param t_number Numeric, the approximate number of transects to generate.
#' @param t_length Numeric, the desired length of transects (in meters).
#' @param t_size Numeric, the size (Y x Y) of the transect cells (in meters).
#' @param buddy_t Logical, default is TRUE,  whether to include a buddy transect
#' for efficient sampling.
#' @param direction Character, default is c("positive", "negative"), the direction from
#' the line layer transects should be built.
#' @param allow_overlaps Logical, default is FALSE, whether to remove areas of line where
#' transects could overlap with other transects.
#' @param stratify_on_length Logical, default is TRUE, whether to stratify where transects are
#' placed based on the length of transect segments. If FALSE then transects are randomly placed
#' across line segments.
#' @return A 'sf' object with transect lines and ID numbers in lat long.
#' @export
make_transects <- function(line_layer_path,
                           # poly_layer_path,
                           t_number,
                           t_length,
                           t_size,
                           buddy_t = TRUE,
                           direction = c("positive", "negative"),
                           allow_overlaps = FALSE,
                           stratify_on_length = TRUE) {
  # browser()

  ## Parameter checks
  stopifnot(is.character(line_layer_path),
            # is.character(poly_layer_path),
            is.numeric(t_number),
            is.numeric(t_length),
            is.numeric(t_size),
            is.logical(buddy_t),
            is.character(direction),
            !any(!grepl("positive|negative", direction)),
            is.logical(allow_overlaps))

  ## Import line and poly layers (puts in UTM)
  line_layer <- get_data(line_layer_path)
  # poly_layer <- get_data(poly_layer_path)

  ## If buddy_t == TRUE and transect number is odd number
  if (buddy_t & (t_number %% 2) != 0) {
    t_number <- t_number + 1
  }

  ## Clip line layer to avoid overlapping transects
  if (!allow_overlaps) {
    line_layer <- clip_to_valid_sampling_area(line_layer, t_length)
  }

  ## Get the number of transect points for each segment in the line layer
  ## this is based on the distance of the transect segment, whether
  ## transects will be in both directions from the line, and whether
  ## a buddy transect will be used
  line_layer$segment <- 1:nrow(line_layer)
  line_layer$n <- get_n(line_layer,
                        t_number,
                        direction,
                        buddy_t,
                        stratify_on_length)

  ## TODO - make function
  ## For each road in layer make transects
  t_lines <- rep(list(NA), nrow(line_layer))
  for (i in 1:nrow(line_layer)) {
    line_shp <- line_layer[i, ]

    ## Make Transect Offsets (in both directions)
    # make coordinates
    tryCatch({
      coords <- sp::coordinates(as(line_shp, "Spatial"))[[1]] %>%
        lapply(as.data.frame) %>%
        # do.call(rbind, .) %>%
        lapply(function(x) {x %>% `names<-`(c("x", "y"))})
    },
    warning=function(w){
      print(paste0('WARNING making coordinates - ', i))
    },
    error=function(e){
      print(paste0('ERROR making coordinates - ', i))
    })

    utm_epsg <- calcUTMzone(line_shp)
    # make line offsets
    tryCatch({
      t_sp <- lapply(coords,
                     make_line_offsets,
                     direction,
                     t_length,
                     utm_epsg) %>%
        fix_t_sp(direction) %>%
        lapply(sf::st_as_sf)
    },
    warning=function(w){
      print(paste0('WARNING making line offsets - ', i))
    },
    error=function(e){
      print(paste0('ERROR making line offsets - ', i))
    })

    tryCatch({
      temp_buff <- sf::st_buffer(line_shp, t_length-1)
      for (j in 1:length(t_sp)) {
        t_sp[[j]] <- suppressWarnings(sf::st_difference(t_sp[[j]], temp_buff))
      }
    },
    warning=function(w){
      print(paste0('WARNING fixing line offsets - ', i))
    },
    error=function(e){
      print(paste0('ERROR fixing line offsets - ', i))
    })

    ## Sample points along transect offset layers
    n <- line_layer$n[i]
    if (n > 0) {
      # browser()

      tryCatch({
        samp_pts <- mapply(get_samp_pts, t_sp, n, utm_epsg)
        t_pts <- list(
          samp_pts = lapply(samp_pts, sf::st_as_sf) %>%
            do.call(rbind, .)
        )
      },
      warning=function(w){
        print(paste0('WARNING sampling points for transect - ', i))
      },
      error=function(e){
        print(paste0('ERROR sampling points for transect - ', i))
      })

      tryCatch({
        if("samp_pts" %in% names(t_pts)) {
          bad_pts <- which(sf::st_is_empty(t_pts$samp_pts))
          if(length(bad_pts) > 0) {
            samp_pts <- lapply(samp_pts, rm_bad_samp_pts)
            bad_samp_pts <- lapply(samp_pts, length) %>% do.call(c, .)
            if (any(bad_samp_pts == 0)) {
              bad_samp_pts <- which(bad_samp_pts == 0)
              samp_pts[[bad_samp_pts]] <- NULL
              t_sp[[bad_samp_pts]] <- NULL
            }
            t_pts$samp_pts <- t_pts$samp_pts[-bad_pts, "x"]
          }
        }
      },
      warning=function(w){
        print(paste0('WARNING removing empty sampling points - ', i))
      },
      error=function(e){
        print(paste0('ERROR removing empty sampling points - ', i))
      })

      ## If buddy_t == TRUE, make buddy transects
      # make buffer around each point and spsample offset layer within buffer
      ## if buddy system is used
      if (buddy_t) {
        t_pts$buddy_pts <- mapply(
          make_buddies,
          samp_pts,
          t_sp,
          MoreArgs = list(t_size = t_size,
                          utm_epsg = utm_epsg)
        ) %>%
          lapply(sf::st_as_sf) %>%
          do.call(rbind, .)

        # tryCatch({
        #
        # },
        # warning=function(w){
        #   print(paste0('WARNING making buddy points - ', i))
        # },
        # error=function(e){
        #   print(paste0('ERROR making buddy points - ', i))
        # })

        tryCatch({
          if("buddy_pts" %in% names(t_pts)) {
            bad_buddies <- which(sf::st_is_empty(t_pts$buddy_pts))
            if(length(bad_buddies) > 0) {
              t_pts$samp_pts <- t_pts$samp_pts[-bad_buddies, "x"]
              t_pts$buddy_pts <- t_pts$buddy_pts[-bad_buddies, "x"]
            }
          }
        },
        warning=function(w){
          print(paste0('WARNING removing empty buddies - ', i))
        },
        error=function(e){
          print(paste0('ERROR removing empty buddies - ', i))
        })
      }

      ## Connect Transects & combine
      if(length(t_pts) > 1) {
        # browser()
        t_lines[[i]] <- lapply(t_pts, make_transect_lines, line_shp) %>%
          do.call(rbind, .)

        # tryCatch({
        #   if(nrow(t_pts$buddy_pts) != 0 & nrow(t_pts$samp_pts) != 0) {
        #
        #   } else {
        #     next
        #   }
        # },
        # warning=function(w){
        #   print(paste0('WARNING making transect lines - ', i))
        # },
        # error=function(e){
        #   print(paste0('ERROR making transect lines - ', i))
        # })
      } else {
        t_lines[[i]] <- lapply(t_pts, make_transect_lines, line_shp) %>%
          do.call(rbind, .)

        # tryCatch({
        #   if(nrow(t_pts$samp_pts) != 0) {
        #
        #   } else {
        #     next
        #   }
        # },
        # warning=function(w){
        #   print(paste0('WARNING making transect lines - ', i))
        # },
        # error=function(e){
        #   print(paste0('ERROR making transect lines - ', i))
        # })
      }
    }
  } ## end for each segment of line layer

  ## Combine transects from all segments & add attribute to add comment
  tryCatch({
    good_t_lines <- sapply(t_lines, function(x) "sf" %in% class(x))
    t_lines <- t_lines[good_t_lines]
  },
  warning=function(w){
    print(paste0('WARNING removing bad transects - '))
  },
  error=function(e){
    print(paste0('ERROR removing bad transects - '))
  })

  tryCatch({
    transect_lines <- do.call(rbind, t_lines)
    transect_lines$comment <- character(nrow(transect_lines))
  },
  warning=function(w){
    print(paste0('WARNING combining transects - '))
  },
  error=function(e){
    print(paste0('ERROR combining transects - '))
  })

  ## remove transects less than 20% of t_length
  tryCatch({
    transect_lines <- transect_lines[transect_lines$distance >= (t_length * 0.8), ]
  },
  warning=function(w){
    print(paste0('WARNING removing transects less than 80% of t_length - '))
  },
  error=function(e){
    print(paste0('ERROR removing transects less than 80% of t_length - '))
  })

  # Return transect lines
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
    point_list[[i]] <- sp::spsample(as(point_list[[i]], "Spatial"),
                                (t_length / (t_size / 2)),
                                type = "regular",
                                offset = 0) %>%
      sf::st_as_sf()
    #ggplot(t_lines_proj[i, ]) +geom_sf() + geom_sf(data = point_list[[i]] )
    #st_distance( point_list[[i]])
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

#' @title Randomly sample points on a line
#'
#' @description Function for returning a random sample with a specified
#' (approximate) number of points. Must provide a line layer as a 'sf'
#' object class and the EPSG code to define the coordinate reference
#' system of the returned points.
#'
#' @param t_sp Multilinestring object of class 'sf'
#' @param n Approximate number of sample points. See sp::spsample()
#' @param epsg Desired EPSG coordinate reference system code for output
#' @return 'sf' class object with sample points
#' @export
get_samp_pts <- function(t_sp, n, epsg) {
  if (nrow(t_sp) > 0) {
    tryCatch({
      temp <- suppressWarnings(sf::st_sample(
        t_sp,
        n,
        type = "regular" # "random"
      )) %>%
        sf::st_as_sf() %>%
        sf::st_cast(to = "POINT") %>%
        sf::st_set_crs(epsg)
      temp <- temp[!st_is_empty(temp$x), ]
    },
    error=function(e){
      temp <- suppressWarnings(sf::st_sample(
        t_sp,
        n,
        type = "random"
      )) %>%
        sf::st_as_sf() %>%
        sf::st_cast(to = "POINT") %>%
        sf::st_set_crs(epsg)
      temp <- temp[!st_is_empty(temp$x), ]
    })
  } else {
    temp <- st_sf(x = st_sfc(lapply(1, function(x) st_point()))) %>%
      st_set_crs(epsg)
  }
  return(temp)
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

## Given a vector (defined by 2 points) and the distance,
## calculate a new vector that is distance away from the original
segment_shift <- function(x, y, d){
  # calculate vector
  v <- c(x[2] - x[1],y[2] - y[1])
  # normalize vector
  v <- v/sqrt((v[1]^2 + v[2]^2))
  # perpendicular unit vector
  vnp <- c(-v[2], v[1])
  return(list(x = c(x[1] + d*vnp[1], x[2] + d*vnp[1]),
              y = c(y[1] + d*vnp[2], y[2] + d*vnp[2])))
}

## Generate a perpendicular line to another a specified distance away
offset_line <- function(line_coords, d, utm_epsg) {

  # allocate memory
  xn <- numeric((nrow(line_coords) - 1) * 2)
  yn <- numeric((nrow(line_coords) - 1) * 2)
  t_df <- data.frame(x = xn, y = yn)
  # make offset line segment line_coords
  for (i in 1:(nrow(line_coords) - 1)) {
    xs <- c(line_coords$x[i], line_coords$x[i + 1])
    ys <- c(line_coords$y[i], line_coords$y[i + 1])
    new.s <- segment_shift(xs, ys, d)
    t_df[(i - 1) * 2 + 1, "x"] <- new.s$x[1]
    t_df[(i - 1) * 2 + 2, "x"] <- new.s$x[2]
    t_df[(i - 1) * 2 + 1, "y"] <- new.s$y[1]
    t_df[(i - 1) * 2 + 2, "y"] <- new.s$y[2]
  }
  # make spatial
  # plot(t_df$x,t_df$y, col = "red", xlim = c(min(rbind(t_df, line_coords)$x), max(rbind(t_df, line_coords)$x)), ylim = c(min(rbind(t_df, line_coords)$y), max(rbind(t_df, line_coords)$y)))
  # points(line_coords$x, line_coords$y, col = "blue")
  line_sp <- sf::st_linestring(as.matrix(t_df)) %>%
    sf::st_sfc()%>%
    sf::st_set_crs(utm_epsg) %>%
    sf::st_cast("MULTILINESTRING")
  return(line_sp)
}
## Make parallel lines to the line shape
make_line_offsets <- function(coords,
                               direction,
                               t_length,
                               utm_epsg) {
  if (length(direction) > 1) {
    t_sp <- list(
      pos_t_sp = offset_line(coords, t_length, utm_epsg),
      neg_t_sp = offset_line(coords, -t_length, utm_epsg)
    ) %>%
      lapply(sf::st_as_sf)
  } else {
    if (direction == "positive") {
      t_sp <- offset_line(coords, t_length, utm_epsg)
    } else {
      t_sp <- offset_line(coords, -t_length, utm_epsg)
    }
  }
  return(t_sp)
}

# Get the number of transects for each road feature, based on distance of each
# road segment. (samples stratified on length of road segment)
get_n <- function(line_layer,
                  t_number,
                  direction,
                  buddy_t,
                  stratify_on_length) {
  if(stratify_on_length) {
    line_layer$len_m <- sf::st_length(line_layer)
    sum_len <- sum(line_layer$len_m)
    n_frac <- line_layer$len_m / sum_len

    # start with random sampling
    n_vec1 <- runif(nrow(line_layer), 0, ifelse(t_number / nrow(line_layer) < 1, 1, t_number / nrow(line_layer))) %>%
      round()
    # now stratify to add the rest
    n_vec2 <- round(t_number * n_frac) %>%
      as.numeric()
    n_vec <- n_vec1 + n_vec2

    sum_multiplier <- 1
    # if (length(direction) > 1) {
    #   n_vec <- round(n_vec / 2)
    #   sum_multiplier <- sum_multiplier * 2
    # }
    if (buddy_t) {
      n_vec <- round(n_vec / 2)
      sum_multiplier <- sum_multiplier * 2
    }
    sum_n <- sum(n_vec)
    # if there are more ore less sample points than the transect number
    if (sum_n != t_number/sum_multiplier) {
      diff <- sum_n - t_number/sum_multiplier
      diff_frac <- n_vec / sum_n
      diff_n <- round(diff * diff_frac)
      n_vec <- n_vec - diff_n
    }
  } else {
    line_layer$len_m <- sf::st_length(line_layer)

    n_vec <- runif(nrow(line_layer), 0, ifelse(t_number / nrow(line_layer) < 1, 1, t_number / nrow(line_layer))) %>%
      round()

    sum_multiplier <- 1
    if (buddy_t) {
      if (all(n_vec == 1 | n_vec == 0)) {
        rm_vec <- runif(nrow(line_layer), 0, 1) %>% round() %>% as.logical()
        n_vec[rm_vec] <- 0
      } else {
        n_vec <- round(n_vec / 2)
      }
      sum_multiplier <- sum_multiplier * 2
    }
    sum_n <- sum(n_vec)

    if (sum_n != t_number/sum_multiplier) {
      if (sum_n > t_number/sum_multiplier) {
        diff <- sum_n - t_number/sum_multiplier
        diff_n <- sample(which(n_vec > 0), diff, replace = TRUE)
        diff_n <- as.data.frame(diff_n) %>% group_by(diff_n) %>% summarize(plus_n=n())
        n_vec[diff_n$diff_n] <- n_vec[diff_n$diff_n] - diff_n$diff_n
      } else {
        diff <- abs(sum_n - t_number/sum_multiplier)
        diff_n <- sample(1:nrow(line_layer), diff, replace = TRUE)
        diff_n <- as.data.frame(diff_n) %>% group_by(diff_n) %>% summarize(plus_n=n())
        n_vec[diff_n$diff_n] <- n_vec[diff_n$diff_n] + diff_n$diff_n
      }
    }
  }

  return(n_vec)
}

## Make points within one transect size and 3 transect sizes
make_buddies <- function(samp_pts_sub,
                         t_sp_l,
                         t_size,
                         utm_epsg) {
  # browser()

  samp_pts_sub <- sf::st_as_sf(samp_pts_sub)
  samp_pts_l <- rep(list(NA), nrow(samp_pts_sub))
  for (i in 1:length(samp_pts_l)) {
    samp_pts_l[[i]] <- samp_pts_sub[i, ]
  }

  buddy_pts <- lapply(samp_pts_l,
                      make_buddy_pt,
                      t_sp_l,
                      t_size,
                      utm_epsg) %>%
    do.call(rbind, .)

  return(buddy_pts)
}

## make the point on a line that's a buddy to another
make_buddy_pt <- function(samp_pts_j,
                           t_sp_l,
                           t_size,
                           utm_epsg) {
  # browser()

  tryCatch({
    # make buffer around each point in buffer_points
    inner_buff <- sf::st_buffer(samp_pts_j, t_size)
    outer_buff <- sf::st_buffer(samp_pts_j, t_size * 2)
    buff <- sf::st_difference(outer_buff, inner_buff)

    # clip line to buffer
    new_line <- suppressWarnings(sf::st_intersection(t_sp_l, buff))

    # sample within buffer around the point
    samps <- get_samp_pts(new_line, 5, utm_epsg) %>%
      sf::st_cast(to = "POINT")
    buddy_pts_j <- samps[runif(1, 1, length(samps)), ]
    return(buddy_pts_j)
  },
  error=function(e){
    buddy_pts_j <- st_sf(x = st_sfc(lapply(1, function(x) st_point()))) %>%
      st_set_crs(utm_epsg)
    return(buddy_pts_j)
  })
}

## make the lines for transects
make_transect_lines <- function(t_pts,
                                line_shp) {
  # convert into lat long
  t_pts <- sf::st_transform(t_pts, 4326)
  line_shp <- sf::st_transform(line_shp, 4326)
  # get points on line perpendicular to transect pts
  # get shortest dist b/w point and line
  temp <- sf::st_coordinates(t_pts) %>%
    .[, c(1, 2)]
  dists <- geosphere::dist2Line(
    temp,
    as(line_shp, "Spatial")
  )
  # get points on main line perpendicular to transect
  # in lat long
  line_pts <- sp::SpatialPointsDataFrame(
    cbind(dists[, "lon"], dists[, "lat"]),
    as.data.frame(dists)
  ) %>%
    sf::st_as_sf() %>%
    sf::st_set_crs(4326)
  line_pts <- line_pts[, "geometry"]

  # make transect lines between points
  line_list <- rep(list(NA), nrow(line_pts))
  line_sdf <- data.frame(id = 1:length(line_list),
                         geometry = NA)
  for (j in 1:length(line_list)) {
    temp <- sf::st_coordinates(t_pts[j, ]) %>%
                                 .[, c(1, 2)]
    line_list[[j]] <- rbind(sf::st_coordinates(line_pts[j, ]) %>% c(),
                            temp %>% c())
    line_list[[j]] <- sf::st_linestring(line_list[[j]])
  }
  # reproject to lat long

  sfc <- do.call(sf::st_sfc, line_list)
  temp <- sf::st_coordinates(t_pts) %>%
    .[,c(1, 2)]
  t_lines <- sf::st_sf(
    data.frame(id = 1:length(line_list),
               azimuth = geosphere::bearing(sf::st_coordinates(line_pts),
                                            temp),
               geom = sfc)
  ) %>%
    sf::st_set_crs(4326)
  # add distance
  t_lines$distance <- dists[, 1]
  # return transect lines in lat long
  return(t_lines)
}

## Relabel 't_sp' with the direction
fix_t_sp <- function(t_sp, direction) {
  if (length(direction) == 1) {
    t_sp <- lapply(t_sp, sf::st_as_sf) %>%
      do.call(rbind, .)
    t_sp_new <- list(t_sp) %>%
      `names<-`(ifelse(direction == "positive",
                       "pos_t_sp", "neg_t_sp"))

  } else {
    t_sp_new <- rep(list(NULL), length(direction))
    for (j in 1:length(t_sp)) {
      t_sp_new[[1]] <- rbind(t_sp_new[[1]], t_sp[[j]][[1]])
      if (length(t_sp[[j]]) > 1) {
        t_sp_new[[2]] <- rbind(t_sp_new[[2]], t_sp[[j]][[2]])
      }
    }
    names(t_sp_new) <- ifelse(direction == "positive",
                              "pos_t_sp",
                              ifelse(direction == "negative",
                                     "neg_t_sp",
                                     c("pos_t_sp", "neg_t_sp")))
  }
  return(t_sp_new)
}

rm_bad_samp_pts <- function(x) {
  bad_pts <- which(sf::st_is_empty(x))
  if(length(bad_pts) > 0) {
    x <- x[-bad_pts]
  }
  return(x)
}



