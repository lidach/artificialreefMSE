#' @title TAC_calc
#'
#' @description calculates Total Allowable Catch (TAC)
#'
#' @param Fref fishing mortality reference points
#' @param Depletion spawning biomass over unfished spawning biomass (stock status from estimation model)
#' @param SSB0 unfished spawning biomass
#' @param SSB_curr spawning biomass from estimation model
#' @param year amount of years to calculate TAC
#' @param M_a maturity at age
#' @param W_a weight at age
#' @param sel_list selectivity list - parameters calculated from estimation model
#' @param alpha the levels of spawning biomass relative to the average unfished spawning biomass at which the catch limit is zero
#' @param beta the levels of spawning biomass relative to the average unfished spawning biomass at which the exploitation rate on which the catch limit is based is set to the proxy for Fmsy
#'

TAC_calc <-
  function(Fref,
           Depletion,
           SSB0,
           SSB_curr,
           VB,
           year,
           Ages,
           R0,
           M_a,
           W_a,
           sel_a,
           alpha,
           beta) {
    # Harvest control rule according to TAC
    TAC <- vector(length = length(year))
    TAC_flag <- NA
    N_a <- vector(length = length(Ages))
    SSB_curr <- SSB_curr[year]
    Depletion <- Depletion[year]
    VB_curr <- VB[year]

    Harv <- 1 - exp(-Fref)
    Z <- M_a + (Fref * sel_a)

    N_a[1] <- R0
    for (a in 2:length(Ages)) {
      N_a[a] <- N_a[a - 1] * exp(-Z[a])
      if (a == length(Ages)) N_a[a] <- (N_a[a - 1] * exp(Z[a])) / (1 - exp(-Z[a]))
    }

    for (i in 1:length(year)) {
      if (Depletion[i] <= alpha) {
        TAC[i] <- 0
        TAC_flag <- TRUE
      }
      if (Depletion[i] >= alpha & Depletion[i] < beta) {
        TAC[i] <- (Harv * VB_curr[i]) * ((SSB_curr[i] - alpha * SSB0) * ((beta * SSB0 / SSB_curr[i]) / (beta * SSB0 - alpha * SSB0)))
        TAC_flag <- TRUE
      }
      if (Depletion[i] > beta) {
        TAC[i] <- Harv * VB_curr[i]
        TAC_flag <- FALSE
      }
    }


    output <- NULL
    output$Ntemp <- N_a
    output$VBtemp <- VB
    output$TAC <- TAC
    output$TAC_flag <- TAC_flag

    return(output)
  } # end of TAC_rule
