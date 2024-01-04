#' @title create_om_input
#'
#' @description function to create life history parameters for simulation model (without temporal or spatial components)
#'
#' Life history parameters
#' @param amax maximum age
#' @param linf von Bertalanffy asymptotic length
#' @param vbk von Bertalanffy growth coefficient
#' @param t0 default = -0.01; von Bertalanffy theoretical age at length = 0
#' @param M natural mortality
#' @param lwa length-weight a parameter
#' @param lwb length-weight b parameter
#' @param agefec default = NULL; fecundity-at-age vector
#' @param h steepness parameter
#' @param binwidth default = 1; width of length bins
#' @param start_ages default = 0; age to start (either 0 or 1)
#'
#' Fishery parameters
#' @param F_comm fishing mortality for commercial fishery
#' @param q catchability
#' @param fi_q catchability for FI fleet
#' @param ss_trip number of fish sampled from each trip
#' @param s1 first selectivity parameter for double normal
#' @param s2 second selectivity parameter for double normal
#' @param s3 third selectivity parameter for double normal
#' @param s4 fourth selectivity parameter for double normal
#' @param s5 fifth selectivity parameter for double normal
#' @param s6 sixth selectivity parameter for double normal
#' @param sizemin minimum size limit - for retention calculation
#' @param fi_50 age at 50% selectivity - fisheries independent fleet
#' @param fi_slope slope for selectivity - fisheries independent fleet
#'
#' Dynamic effort parameters
#' @param persis default = 1; stickyness parameter for effort dynamics (seen in van Poorten and Camp 2019)
#' @param sig1e shape parameter that describes how strongly fishing effort responds to changes in total utility
#' @param max_eff proportion of maximum effort given out
#' @param w_cost power weight for fisher cost function
#'
#' Utility parameters - for recreational fisheries
#' parameter order for "three_param": max PWU, steepness, inflection point
#' @param ut_cpue utility for CPUE
#' @param ut_hpue utility for harvest
#' @param ut_size utility for size
#' @param ut_dist dis-utility for distance
#' @param ut_crowd dis-utility for crowding
#'
#' Error parameters -
#' @param cv_len CV of the growth curve
#' @param cv_eff CV of total effort
#' @param sigma_c standard deviation - observation error of catch data
#' @param sigma_i standard deviation - observation error of index data
#' @param sigma_r standard deviation - process error for recruitment time series
#' @param sigma_fi standard deviation - process error for harvest time series
#' @param theta weighting parameter for dirichlet-multinomial
#' @param theta_fi weighting parameter for dirichlet-multinomial FI fleet
#'
#' Preference movement parameters
#' @param w_dep power weight for depth preference
#' @param w_hab power weight for habitat preference
#' @param w_dist power weight for fish movement distance
#'
#' AR parameters
#' @param ar_flag designed to be a flag for whether or not there artificial reefs are implemented
#' @param ar_prop_M percent change at artificial reef sites
#' @param ar_prop_q percent change of catchability at artificial reef sites
#' NR parameters
#' @param nr_flag designed to be a flag for whether or not there are natural reefs
#' @param nr_prop_M percent change at natural reef sites
#' @param nr_prop_q percent change of catchability at natural reef sites
#'
#' @param nseason specify number of sub-time periods per year; default=1
#'

#################################################
## Input list for create_om_input R function ####
#################################################
input_list <- list(
  amax = 20, # SEDAR 52
  linf = 85.6, # SEDAR 52
  vbk = 0.1919, # SEDAR 52
  t0 = -0.3925, # SEDAR 52
  M = c(2, 1.2, 0.19, 0.15, 0.129, 0.115, 0.106, 0.099, 0.095, 0.091, 0.088, 0.086, 0.085, 0.083, 0.082, 0.081, 0.081, 0.08, 0.08, 0.079, 0.078), # SEDAR 52
  # M = 0.08, # SEDAR 52
  lwa = 4.47e-04,
  lwb = 2.994,
  agefec = c(0, 0, 0.35E6, 2.62E6, 9.07E6, 20.3E6, 34.71E6, 49.95E6, 64.27E6, 76.76E6, 87.15E6, 95.53E6, 102.15E6, 107.3E6, 111.27E6, 114.3E6, 116.61E6, 118.36E6, 119.68E6, 120.67E6, 123.234591E6), # SEDAR 52
  h = 0.99, # SEDAR 52
  binwidth = 1,
  start_ages = 0,
  F_comm = 0.00082, # SEDAR 52 fleet 3 - LL (2016)
  q = 0.00033, #
  fi_q = 0.001,
  ss_trip = 30,
  s1 = 3.19257,
  s2 = -0.84110,
  s3 = -0.07316,
  s4 = -0.21389,
  s5 = -12.43890,
  s6 = -4.45326,
  fi_50 = 8.67398,
  fi_slope = 2.45164,
  persis = 0.5,
  sig1e = 1,
  sizemin = 40.64, # 16" regulation
  max_eff = 1069073.6, # calculated from data
  w_cost = 1,
  ut_cpue = c(1.5, 20, 0.2397095), # tuned
  ut_hpue = c(1.5, 25, 0.1726116), # tuned
  ut_size = c(1.5, 0.02, 50.45908), # tuned
  ut_dist = c(1.5, -0.006, 185.4731), # tuned
  ut_crowd = c(1.5, -0.0009, 745.7763), # tuned
  cv_len = 0.1,
  cv_eff = 0.25,
  sigma_c = 0.2,
  sigma_i = 0.2,
  sigma_r = 0.3,
  sigma_fi = 0.3,
  theta = 1,
  theta_fi = 1,
  w_dep = 1,
  w_hab = 1,
  w_dist = 1,
  ar_flag = 0,
  ar_prop_M = 0,
  ar_prop_q = 0,
  nr_flag = 1,
  nr_prop_M = -0.25,
  nr_prop_q = 0.25,
  nseason = 1
)



################################
## create_om_input function ####
################################
create_om_input <- function(input_list) {
  with(input_list, {
    ## Length at age
    # length bins
    mids <- seq((binwidth / 2), linf * 1.3, by = binwidth) # midlength bins
    highs <- mids + (binwidth / 2) # length bins high
    lows <- mids - (binwidth / 2) # length bins low

    # length and weight at age
    ages <- seq(start_ages, to = (amax + 1 - (1 / nseason)), by = (1 / nseason)) # all ages including age at 0 if start_ages = 0 and seasonality
    la <- linf * (1 - exp(-vbk * (ages - t0))) # length at age
    wa <- lwa * la^lwb # weight at age



    ## Selectivity
    s_fa <- rep(NA, length(ages)) # selectivity by age
    # double-normal selectivity
    # peak - endpoint where selectivity = 1.0
    peak <- s1 + binwidth + ((0.99 * amax - s1 - binwidth) / (1 + exp(-s2)))
    # joiner functions for asc and desc components
    jf1 <- (1 + exp(-20 * (ages - s1) / (1 + abs(ages - s1))))^-1
    jf2 <- (1 + exp(-20 * (ages - peak) / (1 + abs(ages - peak))))^-1
    # t1 and t2
    t1 <- exp(-(start_ages - s1)^2 / exp(s3))
    t2 <- exp(-(amax - peak)^2 / exp(s4))
    # asc and desc portions
    asc <- (1 + exp(-s5))^-1 + (1 - (1 + exp(-s5))^-1) * (exp(-(ages - s1)^2 / exp(s3)) - t1) / (1 - t1)
    desc <- 1 + (((1 + exp(-s6))^-1) - 1) * (exp(-(ages - peak)^2 / exp(s4)) - 1) / (t2 - 1)
    # selectivity by age
    s_fa <- asc * (1 - jf1) + jf1 * ((1 - jf2) + jf2 * desc)
    s_fa <- s_fa / max(s_fa)

    # FI selectivity
    s_fi <- 1 / (1 + exp(-(ages - fi_50) / fi_slope))
    s_fi <- s_fi / max(s_fi)

    # retention
    ret_sigma <- cv_len * sizemin # get sd from CV
    ret <- pnorm((la - sizemin) / ret_sigma) # normal distribution of size limit
    ret <- ret / max(ret)



    ###################
    ## output list ####
    ###################
    output <- NULL
    output <- list(
      nseason = nseason,
      amax = length(ages),
      ages = ages,
      linf = linf,
      vbk = vbk,
      t0 = t0,
      mids = mids,
      highs = highs,
      lows = lows,
      M = M,
      lwa = lwa,
      lwb = lwb,
      agefec = agefec,
      h = h,
      binwidth = binwidth,
      start_ages = start_ages,
      F_comm = F_comm,
      q = q,
      fi_q = fi_q,
      ss_trip = ss_trip,
      persis = persis,
      sig1e = sig1e,
      ret = ret,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      s4 = s4,
      s5 = s5,
      s6 = s6,
      fi_50 = fi_50,
      fi_slope = fi_slope,
      max_eff = max_eff,
      w_cost = w_cost,
      ut_cpue = ut_cpue,
      ut_hpue = ut_hpue,
      ut_size = ut_size,
      ut_dist = ut_dist,
      ut_crowd = ut_crowd,
      la = la,
      wa = wa,
      s_fa = s_fa,
      s_fi = s_fi,
      cv_len = cv_len,
      cv_eff = cv_eff,
      sigma_c = sigma_c,
      sigma_i = sigma_i,
      sigma_r = sigma_r,
      sigma_fi = sigma_fi,
      theta = theta,
      theta_fi = theta_fi,
      w_dep = w_dep,
      w_hab = w_hab,
      w_dist = w_dist,
      ar_flag = ar_flag,
      ar_prop_M = ar_prop_M,
      ar_prop_q = ar_prop_q,
      nr_flag = nr_flag,
      nr_prop_M = nr_prop_M,
      nr_prop_q = nr_prop_q
    ) # end of output list

    return(output)
  }) # end of with function
} # end of create_om_input function
