#' @title AR_rule
#'
#' @description management model that deploys artificial reefs
#'
#' @param coord_file data frame with red snapper habitat suitability, NR, depth, longitude, and latitude coordinates for each cell/site
#' @param AR_coord_file data frame with AR coordinates (up to 2018)
#' @param ARs_fishery_closure
#' @param SSB_trend spawning biomass trend obtained from estimation model
#' increase", "increase_cont", "decrease", "decrease_cont"
#' @param eff_space effort data (spatial) - for artificial reef placements
#' @param recruit_space recruitment data (spatial) - for artificial reef placements
#' @param mgmt_rule rule on how to implement artificial reefs
#' AR_ref_1 - place artificial reefs anywhere
#' AR_ref_2 - place artificial reefs near historical ones
#' AR_strategy - use SSB trend from estimation model to place artificial reefs
#'    increase - place near historic artificial reefs
#'    increase_cont - place in areas of high effort
#'    decrease - place in areas of high recruitment
#'    decrease_cont - place in areas of high recruitment and no fishing at 10% of artificial reefs
#'

AR_rule <- function(coord_file,
                    AR_coord_file,
                    ARs_fishery_closure,
                    SSB_trend,
                    eff_space,
                    recruit_space,
                    mgmt_rule) {
  cluster_rand <- sample(2:4, size = 1) # amount of clusters per year
  num_rand <- sample(30:60, size = 1) # amount of ARs per year

  # assign amount of ARs per cluster randomly
  AR_cluster_assign <- sample(1:cluster_rand, size = num_rand, replace = TRUE)
  AR_cluster_assign <- as.numeric(table(AR_cluster_assign))

  # convert AR sites
  AR_sites_space <- AR_coord_file
  coordinates(AR_sites_space) <- ~ long + lat
  proj4string(AR_sites_space) <- CRS("+init=epsg:4326")
  AR_sites_space <- spTransform(AR_sites_space, CRS("+init=epsg:4326"))

  # convert spatial extent
  sites_space <- coord_file
  sites_space <- sites_space[, 1:2]
  coordinates(sites_space) <- ~ long + lat
  proj4string(sites_space) <- CRS("+init=epsg:4326")
  sites_space <- spTransform(sites_space, CRS("+init=epsg:4326"))
  gridded(sites_space) <- TRUE


  if (mgmt_rule == "AR_ref_1") { # random regardless of density of ARs
    # assign random coordinates within range per cluster
    new_ARs <- sp::spsample(sites_space, n = cluster_rand + 10, type = "random") # n is "approximate" sample size...
    new_ARs <- as.data.frame(new_ARs@coords)
    # no fishery closures
    ARs_fishery_closure_new <- NULL
  } else if (mgmt_rule == "AR_strategy") { # SSB-based strategy
    if (SSB_trend == "increase") {
      # if SSB trend increase, place artificial reefs near other artificial reefs (same as reference strategy 2)
      # kernel density estimation
      dens_est <- MASS::kde2d(x = AR_sites_space$long, y = AR_sites_space$lat, h = 0.1, n = 50, lims = c(range(coord_file$long), range(coord_file$lat)))
      dens_est_mat <- matrix(NA, nrow = nrow(dens_est$z) * nrow(dens_est$z), ncol = 3)
      zs <- as.vector(dens_est$z)
      dens_est_mat[, 3] <- zs
      dens_est_mat[, 1] <- rep(dens_est$x, nrow(dens_est$z))
      dens_est_mat[, 2] <- rep(dens_est$y, each = ncol(dens_est$z))
      quant_value <- quantile(dens_est_mat[, 3], prob = 0.95)
      dens_est_top <- dens_est_mat[dens_est_mat[, 3] > quant_value, ] # "high" densities of ARs
      # need to keep it within model space
      dens_est_top <- as.data.frame(dens_est_top)
      names(dens_est_top) <- c("x", "y", "z")
      coordinates(dens_est_top) <- ~ x + y
      proj4string(dens_est_top) <- CRS(proj4string(sites_space))

      # where to draw coordinates from
      sel_ARs_space <- c(sp::over(dens_est_top, sites_space))
      sel_ARs_space <- dens_est_top[-which(is.na(sel_ARs_space)), ]

      # choose new AR sites
      samp_num <- 1:nrow(sel_ARs_space)
      samp_num <- sample(samp_num, length(AR_cluster_assign))
      new_ARs <- sel_ARs_space[samp_num, ]
      new_ARs <- as.data.frame(new_ARs@coords)
      colnames(new_ARs) <- c("long", "lat")
      # if there were fishery closures, add them back in
      new_ARs <- rbind(new_ARs, ARs_fishery_closure[, c("long", "lat")])
      # no fishery closures
      ARs_fishery_closure_new <- NULL
    } else if (SSB_trend == "increase_cont") {
      # if SSB trend continues to increase, place artificial reefs in areas of high effort
      eff_space_top <- quantile(eff_space, 0.75)
      AR_eff_space <- which(eff_space >= eff_space_top)

      new_ARs <- sp::spsample(sites_space[AR_eff_space, ], n = cluster_rand + 10, type = "random") # n is "approximate" sample size...
      new_ARs <- as.data.frame(new_ARs@coords)
      colnames(new_ARs) <- c("long", "lat")
      # if there were fishery closures, "add" them back in
      new_ARs <- rbind(new_ARs, ARs_fishery_closure[, c("long", "lat")])
      # no fishery closures
      ARs_fishery_closure_new <- NULL
    } else if (SSB_trend == "decrease") {
      # if SSB trend decrease, place artificial reefs in areas of high recruitment
      recruit_space_top <- quantile(recruit_space, 0.75)
      AR_recruit_space <- which(recruit_space > recruit_space_top)

      new_ARs <- sp::spsample(sites_space[AR_recruit_space, ], n = cluster_rand + 10, type = "random") # n is "approximate" sample size...
      new_ARs <- as.data.frame(new_ARs@coords)

      ARs_fishery_closure_new <- ARs_fishery_closure
    } else if (SSB_trend == "decrease_cont") {
      # if SSB trend continues to decrease, place artificial reefs in areas of high recruitment
      # but also close ARs
      ARs_close_num <- ceiling(0.05 * length(AR_sites_space))

      ARs_close_sites <- 1:length(AR_sites_space)
      ARs_close_site <- sample(ARs_close_num, ARs_close_num, replace = FALSE)

      ARs_fishery_closure_new <- cbind(AR_coord_file[ARs_close_site, ], ARs_close_site)
      names(ARs_fishery_closure_new) <- c("long", "lat", "num")
      ARs_fishery_closure_new <- rbind(ARs_fishery_closure, ARs_fishery_closure_new)

      # add ARs still
      recruit_space_top <- quantile(recruit_space, 0.75)
      AR_recruit_space <- which(recruit_space > recruit_space_top)
      new_ARs <- sp::spsample(sites_space[AR_recruit_space, ], n = cluster_rand + 10, type = "random") # n is "approximate" sample size...
      new_ARs <- as.data.frame(new_ARs@coords)
    }
  } else if (mgmt_rule == "AR_ref_2") {
    # place "anywhere", in areas of historic artificial reefs
    # kernel density estimation - finding areas of high density of artificial reefs
    dens_est <- MASS::kde2d(x = AR_sites_space$long, y = AR_sites_space$lat, h = 0.1, n = 50, lims = c(range(coord_file$long), range(coord_file$lat)))
    dens_est_mat <- matrix(NA, nrow = nrow(dens_est$z) * nrow(dens_est$z), ncol = 3)
    zs <- as.vector(dens_est$z)
    dens_est_mat[, 3] <- zs
    dens_est_mat[, 1] <- rep(dens_est$x, nrow(dens_est$z))
    dens_est_mat[, 2] <- rep(dens_est$y, each = ncol(dens_est$z))
    quant_value <- quantile(dens_est_mat[, 3], prob = 0.95)
    dens_est_top <- dens_est_mat[dens_est_mat[, 3] > quant_value, ] # "high" densities of ARs
    # need to keep it within model space
    dens_est_top <- as.data.frame(dens_est_top)
    names(dens_est_top) <- c("x", "y", "z")
    coordinates(dens_est_top) <- ~ x + y
    proj4string(dens_est_top) <- CRS(proj4string(sites_space))

    # where to draw coordinates from
    sel_ARs_space <- c(sp::over(dens_est_top, sites_space))
    sel_ARs_space <- dens_est_top[-which(is.na(sel_ARs_space)), ]

    # choose new AR sites
    samp_num <- 1:nrow(sel_ARs_space)
    samp_num <- sample(samp_num, length(AR_cluster_assign))
    new_ARs <- sel_ARs_space[samp_num, ]
    new_ARs <- as.data.frame(new_ARs@coords)

    # no fishery closures
    ARs_fishery_closure_new <- NULL
  }

  # create new ARs
  new_ARs_coords <- list()
  for (i in 1:cluster_rand) new_ARs_coords[[i]] <- jitter(rep(c(as.numeric(new_ARs[i, ])), AR_cluster_assign[i]), 0.001)
  new_ARs_coords <- matrix(unlist(new_ARs_coords), ncol = 2, byrow = TRUE)
  new_ARs_coords <- as.data.frame(new_ARs_coords)
  names(new_ARs_coords) <- c("long", "lat")



  ##############
  ## Export ####
  ##############
  output <- NULL
  output <- list(
    new_ARs_coords = new_ARs_coords,
    ARs_fishery_closure = ARs_fishery_closure_new
  )

  return(output)
} # end of AR_rule function
