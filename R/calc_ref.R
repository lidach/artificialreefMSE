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
                    sel_list,
                    fec_a,
                    R0,
                    ref = FALSE)
{
  N_af <- vector(length = length(Ages))
    N_af[1] <- R0
  N_a0 <- vector(length = length(Ages))
    N_a0[1] <- R0

  # calculate selectivity
  sel_a <- rep(NA, length(Ages))
  Amax <- length(Ages)
  binwidth <- 1
  # peak - endpoint where selectivity = 1.0
    peak <- sel_list$p1+binwidth+((0.99*Amax-sel_list$p1-binwidth)/(1+exp(-sel_list$p2)))
  # joiner functions for asc and desc components
    jf1 <- (1+exp(-20*(Ages-sel_list$p1)/(1+abs(Ages-sel_list$p1))))^-1
    jf2 <- (1+exp(-20*(Ages-peak)/(1+abs(Ages-peak))))^-1
  # t1 and t2
    t1 <- exp(-(Ages[1]-sel_list$p1)^2/exp(sel_list$p3))
    t2 <- exp(-(Amax-peak)^2/exp(sel_list$p4))
  # asc and desc portions
    asc <- (1+exp(-sel_list$p5))^-1+(1-(1+exp(-sel_list$p5))^-1)*(exp(-(Ages-sel_list$p1)^2/exp(sel_list$p3))-t1)/(1-t1)
    desc <- 1+(((1+exp(-sel_list$p6))^-1)-1)*(exp(-(Ages-peak)^2/exp(sel_list$p4))-1)/(t2-1)
  # selectivity by age
    sel_a <- asc*(1-jf1)+jf1*((1-jf2)+jf2*desc)
    sel_a <- sel_a/max(sel_a)

  FM_a <- vector(length = length(Ages))
    FM_a <- sel_a*FM

  # equilibrium fished
  for(a in 2:length(Ages)){
    if(a < length(Ages)) N_af[a] <- N_af[a-1] * exp(-M_a[a]-FM_a[a-1])
    if(a == length(Ages)) N_af[length(Ages)] <- (N_af[a-1]*exp(-M_a[a]-FM_a[a-1]))/(1-exp(-M_a[a]-FM_a[a-1]))
  } 
  # equilibrium unfished
  for(a in 2:length(Ages)){
    if(a < length(Ages)) N_a0[a] <- N_a0[a-1] * exp(-M_a[a])
    if(a == length(Ages)) N_a0[length(Ages)] <- (N_a0[a-1]*exp(-M_a[a]))/(1-exp(-M_a[a]))
  } 
  Nofish <- sum(N_a0*fec_a)
  Fish <- sum(N_af*fec_a)

  ratio <- Fish/Nofish

  if(ref == FALSE) return(ratio)

  if(ref != FALSE){
    diff <- ref - ratio
    return(diff)
  }

} # end of calc_ref function

