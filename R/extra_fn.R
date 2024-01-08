##########################
## Calculate distance ####
##########################
# code from N. Fisch
## PC dist - calculate distance to each site
earth_dist <- function(long1, lat1, long2, lat2) {
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


#########################
## convert_AR_coords ####
#########################
convert_ar_coords <- function(nsites,
                              ar_coord_file,
                              coord_file) {
  match_nrst <- rep(NA, nrow(coord_file))
  new_ar_site <- rep(NA, nrow(ar_coord_file))
  for (i in 1:nrow(ar_coord_file)) {
    for (j in 1:nrow(coord_file)) {
      match_nrst[j] <- earth_dist(
        ar_coord_file$long[i], ar_coord_file$lat[i],
        coord_file$long[j], coord_file$lat[j]
      )
    }
    new_ar_site[i] <- which.min(match_nrst)
  }
  # sum for ar composition
  ar_comp <- table(new_ar_site)
  ar_sites <- as.numeric(names(ar_comp)) # site numbers and coordinates
  ar_comp2 <- rep(0, nsites)
  ar_comp2[ar_sites] <- as.numeric(ar_comp)
  ar_coords <- coord_file[ar_sites, 1:2]

  # output
  output <- list(
    ar_comp = ar_comp2,
    ar_sites = ar_sites,
    ar_coords = ar_coords
  )
  return(output)
}


###############
## wtd_var ####
###############
wtd_var <- function(x, weights=NULL, normwt=FALSE, na.rm=TRUE,
                    method = c('unbiased', 'ML'))
  ## By Benjamin Tyner <btyner@gmail.com> 2017-0-12
{
  method <- match.arg(method)
  if(! length(weights)) {
    if(na.rm) x <- x[!is.na(x)]
    return(var(x))
  }
  
  if(na.rm) {
    s       <- !is.na(x + weights)
    x       <- x[s]
    weights <- weights[s]
  }

  if(normwt)
    weights <- weights * length(x) / sum(weights)

  if(normwt || method == 'ML')
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))

  # the remainder is for the special case of unbiased frequency weights
  sw  <- sum(weights)
  if(sw <= 1)
      warning("only one effective observation; variance estimate undefined")
  
  xbar <- sum(weights * x) / sw
  sum(weights*((x - xbar)^2)) / (sw - 1)
}


##################
## dtruncnorm ####
##################
dtruncnorm <- function(x, a, b, mean, sd) {
  res <- numeric(length(x))

  for (i in seq_along(x)) {
    ca <- a
    cb <- b
    cx <- x[i %% length(x) + 1]

    if (ca <= cx && cx <= cb) {
      cmean <- mean[i %% length(mean) + 1]
      csd <- sd[i %% length(sd) + 1]

      c1 <- pnorm(ca, cmean, csd, lower.tail = TRUE, log.p = FALSE)
      c2 <- pnorm(cb, cmean, csd, lower.tail = TRUE, log.p = FALSE)
      c3 <- csd * (c2 - c1)
      c4 <- dnorm((cx - cmean) / csd, mean = 0, sd = 1, log = TRUE)

      if (!is.finite(log(c3))) {
        res[i] <- 1.0 / (cb - ca)
      } else {
        res[i] <- exp(c4 - log(c3))
      }
    } else {
      res[i] <- 0.0
    }
  }

  return(res)
}


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
