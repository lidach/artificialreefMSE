#' @title calc_ref
#'
#' @description calculate reference points
#'
#' @param Ages vector of ages
#' @param M_a natural mortality at age
#' @param FM fishing mortality
#' @param sel_list list of selectivity parameters for double normal
#' @param fec_a fecundity at age
#' @param R0 unfished recruitment
#' @param ref calculate reference point?
#'

calc_ref <- function(Ages,
                     M_a,
                     FM,
                     sel_a,
                     fec_a,
                     R0,
                     ref = FALSE) {
  N_af <- vector(length = length(Ages))
  N_af[1] <- R0
  N_a0 <- vector(length = length(Ages))
  N_a0[1] <- R0

  FM_a <- vector(length = length(Ages))
  FM_a <- sel_a * FM

  # equilibrium fished
  for (a in 2:length(Ages)) {
    if (a < length(Ages)) N_af[a] <- N_af[a - 1] * exp(-M_a[a] - FM_a[a - 1])
    if (a == length(Ages)) N_af[length(Ages)] <- (N_af[a - 1] * exp(-M_a[a] - FM_a[a - 1])) / (1 - exp(-M_a[a] - FM_a[a - 1]))
  }
  # equilibrium unfished
  for (a in 2:length(Ages)) {
    if (a < length(Ages)) N_a0[a] <- N_a0[a - 1] * exp(-M_a[a])
    if (a == length(Ages)) N_a0[length(Ages)] <- (N_a0[a - 1] * exp(-M_a[a])) / (1 - exp(-M_a[a]))
  }
  Nofish <- sum(N_a0 * fec_a)
  Fish <- sum(N_af * fec_a)

  ratio <- Fish / Nofish

  if (ref == FALSE) {
    return(ratio)
  }

  if (ref != FALSE) {
    diff <- ref - ratio
    return(diff)
  }
} # end of calc_ref function
