##########################
## Calculate distance ####
##########################
## PC dist - calculate distance to each site
earth.dist <- function(long1, lat1, long2, lat2) {
  rad <- pi / 180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat / 2))^2 + cos(a1) * cos(b1) * (sin(dlon / 2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.137 # Mean Earth Radius
  d <- R * c
  return(d)
}



###########################
## fit_tmb (TMBhelper) ####
###########################
fit_tmb <- function(obj, fn = obj$fn, gr = obj$gr, startpar = NULL, lower = -Inf,
                    upper = Inf, getsd = TRUE, control = list(
                      eval.max = 10000,
                      iter.max = 10000, trace = 1
                    ), bias.correct = FALSE, bias.correct.control = list(
                      sd = FALSE,
                      split = NULL, nsplit = NULL, vars_to_correct = NULL
                    ),
                    savedir = NULL, loopnum = 2, newtonsteps = 0, n = Inf, getReportCovariance = FALSE,
                    getJointPrecision = FALSE, getHessian = FALSE, quiet = FALSE,
                    start_time_elapsed = as.difftime("0:0:0"), ...) {
  if (is.null(startpar)) {
    startpar <- obj$par
  }
  if (bias.correct == TRUE & is.null(obj$env$random)) {
    message("No random effects detected in TMB model, so overriding user input to `TMBhelper::fit_tmb` to instead specify `bias.correct=FALSE`")
    bias.correct <- FALSE
  }
  if (getReportCovariance == FALSE) {
    message("Note that `getReportCovariance=FALSE` causes an error in `TMB::sdreport` when no ADREPORTed variables are present")
  }
  List <- list(...)
  combine_lists <- function(default, input) {
    output <- default
    for (i in seq_along(input)) {
      if (names(input)[i] %in% names(default)) {
        output[[names(input)[i]]] <- input[[i]]
      } else {
        output <- c(output, input[i])
      }
    }
    return(output)
  }
  BS.control <- list(
    sd = FALSE, split = NULL, nsplit = NULL,
    vars_to_correct = NULL
  )
  BS.control <- combine_lists(default = BS.control, input = bias.correct.control)
  nlminb.control <- list(
    eval.max = 10000, iter.max = 10000,
    trace = 0
  )
  nlminb.control <- combine_lists(
    default = nlminb.control,
    input = control
  )
  start_time <- Sys.time()
  parameter_estimates <- nlminb(
    start = startpar, objective = fn,
    gradient = gr, control = nlminb.control, lower = lower,
    upper = upper
  )
  for (i in seq(2, loopnum, length = max(0, loopnum - 1))) {
    Temp <- parameter_estimates[c("iterations", "evaluations")]
    parameter_estimates <- nlminb(
      start = parameter_estimates$par,
      objective = fn, gradient = gr, control = nlminb.control,
      lower = lower, upper = upper
    )
    parameter_estimates[["iterations"]] <- parameter_estimates[["iterations"]] +
      Temp[["iterations"]]
    parameter_estimates[["evaluations"]] <- parameter_estimates[["evaluations"]] +
      Temp[["evaluations"]]
  }
  for (i in seq_len(newtonsteps)) {
    g <- as.numeric(gr(parameter_estimates$par))
    h <- optimHess(parameter_estimates$par, fn = fn, gr = gr)
    parameter_estimates$par <- parameter_estimates$par - solve(
      h,
      g
    )
    parameter_estimates$objective <- fn(parameter_estimates$par)
  }
  parameter_estimates <- parameter_estimates[c(
    "par",
    "objective", "iterations", "evaluations"
  )]
  parameter_estimates[["time_for_MLE"]] <- Sys.time() -
    start_time
  parameter_estimates[["max_gradient"]] <- max(abs(gr(parameter_estimates$par)))
  parameter_estimates[["Convergence_check"]] <- ifelse(parameter_estimates[["max_gradient"]] <
    1e-04, "There is no evidence that the model is not converged",
  "The model is likely not converged"
  )
  parameter_estimates[["number_of_coefficients"]] <- c(
    Total = length(unlist(obj$env$parameters)),
    Fixed = length(startpar), Random = length(unlist(obj$env$parameters)) -
      length(startpar)
  )
  parameter_estimates[["AIC"]] <- TMBAIC(opt = parameter_estimates)
  if (n != Inf) {
    parameter_estimates[["AICc"]] <- TMBAIC(
      opt = parameter_estimates,
      n = n
    )
    parameter_estimates[["BIC"]] <- TMBAIC(
      opt = parameter_estimates,
      p = log(n)
    )
  }
  parameter_estimates[["diagnostics"]] <- data.frame(
    Param = names(startpar),
    starting_value = startpar, Lower = lower, MLE = parameter_estimates$par,
    Upper = upper, final_gradient = as.vector(gr(parameter_estimates$par))
  )
  if (getsd == TRUE) {
    sd_time <- Sys.time()
    h <- optimHess(parameter_estimates$par, fn = fn, gr = gr)
    if (is.character(try(chol(h), silent = TRUE))) {
      warning("Hessian is not positive definite, so standard errors are not available")
      if (!is.null(savedir)) {
        capture.output(parameter_estimates, file = file.path(
          savedir,
          "parameter_estimates.txt"
        ))
      }
      return(list(opt = parameter_estimates, h = h))
    }
    if (bias.correct == FALSE | is.null(BS.control[["vars_to_correct"]])) {
      if (!is.null(BS.control[["nsplit"]])) {
        if (BS.control[["nsplit"]] == 1) {
          BS.control[["nsplit"]] <- NULL
        }
      }
      parameter_estimates[["SD"]] <- TMB::sdreport(
        obj = obj,
        par.fixed = parameter_estimates$par, hessian.fixed = h,
        bias.correct = bias.correct, bias.correct.control = BS.control[c(
          "sd",
          "split", "nsplit"
        )], getReportCovariance = getReportCovariance,
        getJointPrecision = getJointPrecision, ...
      )
    } else {
      if ("ADreportIndex" %in% names(obj$env)) {
        Which <- as.vector(unlist(obj$env$ADreportIndex()[BS.control[["vars_to_correct"]]]))
      } else {
        parameter_estimates[["SD"]] <- TMB::sdreport(
          obj = obj,
          par.fixed = parameter_estimates$par, hessian.fixed = h,
          bias.correct = FALSE, getReportCovariance = FALSE,
          getJointPrecision = FALSE, ...
        )
        Which <- which(rownames(summary(
          parameter_estimates[["SD"]],
          "report"
        )) %in% BS.control[["vars_to_correct"]])
      }
      if (!is.null(BS.control[["nsplit"]]) && BS.control[["nsplit"]] >
        1) {
        Which <- split(Which, cut(seq_along(Which), BS.control[["nsplit"]]))
      }
      Which <- Which[sapply(Which, FUN = length) > 0]
      if (length(Which) == 0) {
        Which <- NULL
      }
      message(paste0(
        "Bias correcting ", length(Which),
        " derived quantities"
      ))
      parameter_estimates[["SD"]] <- TMB::sdreport(
        obj = obj,
        par.fixed = parameter_estimates$par, hessian.fixed = h,
        bias.correct = TRUE, bias.correct.control = list(
          sd = BS.control[["sd"]],
          split = Which, nsplit = NULL
        ), getReportCovariance = getReportCovariance,
        getJointPrecision = getJointPrecision, ...
      )
    }
    parameter_estimates[["Convergence_check"]] <- ifelse(parameter_estimates$SD$pdHess ==
      TRUE, parameter_estimates[["Convergence_check"]],
    "The model is definitely not converged"
    )
    parameter_estimates[["time_for_sdreport"]] <- Sys.time() -
      sd_time
    if (getHessian == TRUE) {
      parameter_estimates[["hessian"]] <- h
    }
  }
  parameter_estimates[["time_for_run"]] <- Sys.time() -
    start_time + start_time_elapsed
  if (!is.null(savedir)) {
    save(parameter_estimates, file = file.path(savedir, "parameter_estimates.RData"))
    capture.output(parameter_estimates, file = file.path(
      savedir,
      "parameter_estimates.txt"
    ))
  }
  if (quiet == FALSE & parameter_estimates[["Convergence_check"]] !=
    "There is no evidence that the model is not converged") {
    message("#########################")
    message(parameter_estimates[["Convergence_check"]])
    message("#########################")
  }
  return(parameter_estimates)
}


TMBAIC <- function(opt, p = 2, n = Inf) {
  k <- length(opt[["par"]])
  if (all(c("par", "objective") %in% names(opt))) {
    negloglike <- opt[["objective"]]
  }
  if (all(c("par", "value") %in% names(opt))) {
    negloglike <- opt[["value"]]
  }
  Return <- p * k + 2 * negloglike + 2 * k * (k + 1) / (n - k -
    1)
  return(Return)
}



#########################
## convert_AR_coords ####
#########################
#' @title convert_AR_coords
#'
#' @description convert artificial reef coordinates to site numbers and amount per site
#'
#' @param nsites number of cells/sites in model
#' @param AR_coord_file artifical reefs coordinates
#' @param coord_file data frame with red snapper habitat suitability, NR, depth, longitude, and latitude coordinates for each cell/site
#'

convert_AR_coords <-
  function(nsites,
           AR_coord_file,
           coord_file) {
    # match artificial reef coordinates with spatial model cells
    AR_comp <- hutilscpp::match_nrst_haversine(AR_coord_file$long, AR_coord_file$lat,
      coord_file$long, coord_file$lat,
      cartesian_R = 1
    )
    AR_comp <- coord_file[as.numeric(unlist(AR_comp[, 1])), 1:2]
    # get amount of artificial reefs per spatial cell
    AR_comp2 <- group_by(AR_comp, long, lat) %>%
      count(long, lat)
    AR_comp2 <- as.data.frame(AR_comp2)
    # integrate artificial reef numbers into coord_file
    AR_comp3 <- hutilscpp::match_nrst_haversine(coord_file$long, coord_file$lat,
      AR_comp2$long, AR_comp2$lat,
      cartesian_R = 1
    )
    AR_comp3 <- as.numeric(unlist(AR_comp3[, 1]))
    AR_comp3 <- unique(AR_comp3)
    AR_num <- rep(0, nsites)
    AR_num[AR_comp3] <- AR_comp2[AR_comp3, 3]

    return(AR_num)
  } # end of convert_AR_coords



#################
## plot_grid ####
#################
#' @title plot_grid
#'
#' @description function to plot spatial model (includes grid of sites, artificial and natural reef locations, landing sites, and sea depth in meters)
#'
#' @param coord_file data frame with depth, longitude, and latitude coordinates for each cell/site
#' @param AR_sites data frame with longitude and latitude coordinates of artificial reefs
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param sea_col color for the coord_file (color of the ocean)
#' @param xlim_map x axis ranges for overall map (longitude); c(long1, long2)
#' @param ylim_map y axis ranges for overall map (latitude); c(lat1, lat2)
#' @param maxdepth maximum depth in map for plotting
#' @param save_plot save the plot to local directory? TRUE/FALSE
#' @param plot_legend
#' @param scenario_name
#' @param savedir directory name to save plot to, only works if save_plot = TRUE

plot_grid <-
  function(coord_file,
           AR_sites,
           landing_site_file,
           sea_col,
           xlim_map,
           ylim_map,
           maxdepth,
           save_plot = FALSE,
           plot_legend = TRUE,
           scenario_name,
           savedir) {
    coord_file_sp <- coord_file[, c("long", "lat", "depth")]
    colnames(coord_file_sp) <- c("x", "y", "depth")
    coordinates(coord_file_sp) <- ~ x + y
    proj4string(coord_file_sp) <- CRS("+init=epsg:4326")
    coord_file_sp <- spTransform(coord_file_sp, CRS("+init=epsg:4326"))
    gridded(coord_file_sp) <- TRUE
    r <- raster(coord_file_sp)
    projection(r) <- CRS("+init=epsg:4326")
    newmap <- getMap(resolution = "high")

    NR_long <- coord_file[which(coord_file$HB > 0), "long"]
    NR_lat <- coord_file[which(coord_file$HB > 0), "lat"]

    AR_long <- coord_file[AR_sites, "long"]
    AR_lat <- coord_file[AR_sites, "lat"]

    ## plot
    if (save_plot == TRUE) jpeg(filename = paste0(savedir, "/", scenario_name, "_map.jpeg"), pointsize = 16, res = 300, units = "in", width = 8, height = 8)
    if (plot_legend == FALSE) par(oma = c(0, 0, 0, 0))
    sp::plot(r, col = sea_col, legend = FALSE)
    sp::plot(newmap, add = TRUE, xlim = xlim_map, ylim = ylim_map, asp = 1, col = "gray95")
    if (plot_legend == TRUE) {
      sp::plot(r,
        add = TRUE, legend.only = TRUE, breaks = ceiling(seq(0, maxdepth, length = length(sea_col))), col = sea_col, legend.width = 1, legend.shrink = 0.5,
        axis.args = list(
          at = seq(0, maxdepth, 50),
          labels = seq(0, maxdepth, 50),
          cex.axis = 0.7
        )
      )
    }
    points(NR_long, NR_lat, cex = 1.5, pch = 0)
    points(AR_long, AR_lat, cex = 0.9, pch = 17)
    points(landing_site_file$long, landing_site_file$lat, cex = 1.3, pch = 13, col = "red")
    if (plot_legend) legend("bottomright", legend = c("natural reef", "artificial reef", "landing site"), pch = c(0, 17, 13), col = c("black", "black", "red"), bg = "white")
    if (save_plot == TRUE) dev.off()
  } # end of plot_grid function



#####################
## res_site_plot ####
#####################
#' @title res_site_plot
#'
#' @description function to plot results of models spatially
#'
#' @param coord_file data frame with depth, longitude, and latitude coordinates for each cell/site
#' @param AR_sites data frame with longitude and latitude coordinates of artificial reefs
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param res_col color for the coord_file (color of the ocean)

res_site_plot <-
  function(coord_file,
           AR_sites,
           res,
           res_col,
           cuts,
           save_plot = FALSE,
           scenario_name,
           savedir,
           plot_legend,
           landing_site_file = NULL,
           plot_landingsites = FALSE) {
    coord_file_res <- cbind(coord_file[, c("long", "lat")], res)
    colnames(coord_file_res) <- c("x", "y", "output")
    coordinates(coord_file_res) <- ~ x + y
    proj4string(coord_file_res) <- CRS("+init=epsg:4326")
    coord_file_res <- spTransform(coord_file_res, CRS("+init=epsg:4326"))
    gridded(coord_file_res) <- TRUE
    r <- raster(coord_file_res)
    projection(r) <- CRS("+init=epsg:4326")

    NR_long <- coord_file[which(coord_file$HB > 0), "long"]
    NR_lat <- coord_file[which(coord_file$HB > 0), "lat"]

    AR_long <- coord_file[AR_sites, "long"]
    AR_lat <- coord_file[AR_sites, "lat"]

    ## plot
    if (save_plot == TRUE) jpeg(filename = paste0(savedir, "/", scenario_name, "_siteoutput.jpeg"), pointsize = 16, res = 300, units = "in", width = 8, height = 8)
    sp::plot(r, col = res_col, breaks = cuts, asp = 1.5, box = FALSE, legend = FALSE, axes = FALSE)
    points(NR_long, NR_lat, cex = 1.5, pch = 0)
    points(AR_long, AR_lat, cex = 1, pch = 17)
    if (plot_landingsites) points(landing_site_file$long, landing_site_file$lat, cex = 1.3, pch = 13, col = "red")
    if (plot_legend) legend("bottomright", legend = c("natural reef", "artificial reef"), pch = c(0, 17), col = c("black", "black"), bg = "white")
    if (save_plot == TRUE) dev.off()
  } # end of res_site_plot function



#####################
## AddInsetMapZS ####
#####################
AddInsetMapZS <- function(p, col = c("#D8D8D8", "#BFA76F"), main.label = list(label = NA, adj = NULL), sub.label = list(label = NA, adj = NULL), loc = c("bottomleft", "topleft", "topright", "bottomright"), inset = 0.02, width = NULL, bg = c("#FFFFFFE7", "#FFFFFFCC")) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  loc <- match.arg(loc)
  if (!inherits(p, c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    stop("polygon 'p' is the incorrect class")
  }
  usr <- graphics::par("usr")
  crds <- cbind(c(usr[1:2], usr[2:1], usr[1]), c(rep(
    usr[3],
    2
  ), rep(usr[4], 2), usr[3]))
  b <- SpatialPolygons(list(Polygons(list(Polygon(crds)), "bbox")),
    proj4string = crs(p)
  )
  b@bbox <- matrix(c(-87, -83.2333, 29, 30.4), ncol = 2, byrow = TRUE) ## LC modified this
  if (length(rgeos::gIntersection(p, b)) == 0) {
    stop("user coordinates of the plotting region do not intersect polygon")
  }
  ext <- extent(rgeos::gUnion(p, b))
  if (is.null(width)) {
    dx <- 0.2 * diff(usr[1:2])
  } else {
    dx <- width * (diff(usr[1:2]) / graphics::par("pin")[1])
  }
  dy <- dx * (diff(ext[3:4]) / diff(ext[1:2]))
  if (length(inset) == 1) {
    inset <- rep(inset, 2)
  }
  padx <- inset[1] * diff(usr[1:2])
  pady <- inset[2] * diff(usr[3:4])
  if (loc == "bottomleft") {
    loc <- c(usr[1] + padx, usr[3] + pady)
  } else if (loc == "topleft") {
    loc <- c(usr[1] + padx, usr[4] - pady - dy)
  } else if (loc == "topright") {
    loc <- c(usr[2] - padx - dx, usr[4] - pady - dy)
  } else if (loc == "bottomright") {
    loc <- c(usr[2] - padx - dx, usr[3] + pady)
  }
  graphics::rect(loc[1], loc[2], loc[1] + dx, loc[2] + dy,
    col = bg[1], border = NA
  )
  plt <- c(graphics::grconvertX(
    c(loc[1], loc[1] + dx), "user",
    "nfc"
  ), graphics::grconvertY(
    c(loc[2], loc[2] + dy),
    "user", "nfc"
  ))
  graphics::par(plt = plt, bg = bg[2], new = TRUE)
  xlim <- range(ext[1:2])
  ylim <- range(ext[3:4])
  graphics::plot.window(xlim = xlim, ylim = ylim)
  plot(p, col = col[1], border = NA, lwd = 0.25, add = TRUE)
  plot(extent(b), col = col[2], border = "#090909", lwd = 3, add = TRUE)
  plot(p, col = NA, border = "#090909", lwd = 0.25, add = TRUE)
  if (!is.na(main.label[[1]])) {
    x <- coordinates(rgeos::gUnaryUnion(p))[1, ]
    text(x[1], x[2],
      labels = main.label[[1]], adj = main.label$adj,
      cex = 0.7, font = 2
    )
  }
  if (!is.na(sub.label[[1]])) {
    x <- coordinates(rgeos::gUnaryUnion(b))[1, ]
    text(x[1], x[2],
      labels = sub.label[[1]], adj = sub.label$adj,
      cex = 0.6
    )
  }
  graphics::box(lwd = 0.5)
  invisible(NULL)
}



################
## col2rgbA ####
################
col2rgbA <- function(color, transparency) {
  rgb(t(col2rgb(color)) / 255, alpha = transparency)
}



##################
## calc_error ####
##################
#' @title calc_error
#'
#' @description function to calculate median absolute relative error, median relative error, and errors within 50% and 90%
#'
#' @param true values from model that is being compared to, i.e. "true" estimates of metric
#' @param est values being compared to "true" values, "estimated"
#' @param nsites number of sites in model

calc_error <- function(true, est) {
  list <- list(abs_error = NA, rel_error = NA, MRE = NA, MARE = NA)
  # relative error
  list$rel_error[i] <- (est[i] - true[i]) / true[i]

  # median absolute relative error
  list$MARE <- median(abs(list$rel_error), na.rm = TRUE)
  # median relative error
  list$MRE <- median(list$rel_error, na.rm = TRUE)
  
  return(list)
}
