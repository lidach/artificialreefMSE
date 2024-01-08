#' @title create_em_input
#'
#' @description function to create data and define parameters for the estimation model (ASAM in TMB)
#'
#' @param sim_data data from the operating model
#' @param lh_list operating model inputs
#' @param years_run subset n_year in operating model to run estimation model on
#'

create_em_input <-
  function(sim_data,
           lh_list,
           years_run) {
    ## Estimation model data
    # lh inputs
    n_year <- length(years_run) # total number of n_year
    n_age <- sim_data$amax # number of n_age

    # data from observation model
    harv <- sim_data$catch_obs_save[years_run] # observed total fishery harvest (catch)
    cpue <- sim_data$cpue_index_save[years_run] # observed abundance index (CPUE)
    pa <- sim_data$pa_save[years_run, ] # observed age compositon
    pa_fi <- sim_data$pa_fi_save[years_run, ] # observed age composition (fishery-independent - FI)
    ss <- rep(lh_list$amax * lh_list$ss_trip, n_year) # sample size (for age composition)
    ss_fi <- rep(lh_list$amax * lh_list$ss_trip, n_year) # sample size (for FI age composition

    ages<- lh_list$ages # vector of ages
    la <- lh_list$la # length-at-age
    wa <- lh_list$wa # weight-at-age
    M_vec <- lh_list$M # natural mortality-at-age
    fec_a <- lh_list$agefec # fecundity-at-age

    ## Estimation model parameters
    h <- lh_list$h # steepness
    log_ro <- sim_data$log_ro # unfished recruitment
    log_q <- log(lh_list$q) # catchability
    log_fi_q <- log(lh_list$fi_q) # catchability (for FI survey)

    # errors
    log_sigma_r <- log(lh_list$sigma_r) # recruitment variations
    log_sigma_c <- log(lh_list$sigma_c) # observation error in catch/harvest
    log_cv_i <- log(lh_list$sigma_i) # observation error in CPUE index
    log_theta <- log(lh_list$theta) # weighting parameter for dirichlet-multinomial
    log_theta_fi <- log(lh_list$theta_fi) # weighting parameter for dirichlet-multinomial (FI survey)

    log_recruit_devs <- log(rep(0.5, n_year)) # recruitment deviations
    log_Fint <- log(rep(0.1, n_year)) # fishing mortality vector

    # selectivity
    p1 <- lh_list$s1 # selectivity - peak: beginning size for plateau
    p2 <- lh_list$s2 # selectivity - top: width of plateau, as logistic between peak and max size
    p3 <- lh_list$s3 # selectivity - ascending width
    p4 <- lh_list$s4 # selectivity - descending width
    p5 <- lh_list$s5 # selectivity - initial: selectivity at first bin
    p6 <- lh_list$s6 # selectivity - final: selectivity at last bin
    pfi_50 <- lh_list$fi_50 # FI age at 50% selectivity
    pfi_slope <- lh_list$fi_slope # FI slope for selectivity


    ###################
    ## output list ####
    ###################
    output <- list( ## lh inputs and data
      data = list(
        n_year = n_year,
        n_age = n_age,
        harv = harv,
        obs_i = cpue,
        obs_pa = pa,
        obs_pa_fi = pa_fi,
        obs_ss = ss,
        obs_ss_fi = ss_fi,
        ages = ages,
        la = la,
        wa = wa,
        M_vec = M_vec,
        fec_a = fec_a,
        amin = min(n_age),
        amax = max(n_age)
      ),
      ## parameters
      pars = list(
        h = h,
        log_ro = log_ro,
        log_q = log_q,
        log_fi_q = log_fi_q,
        log_sigma_r = log_sigma_r,
        log_sigma_c = log_sigma_c,
        log_cv_i = log_cv_i,
        log_theta = log_theta,
        log_theta_fi = log_theta_fi,
        log_Fint = log_Fint,
        log_recruit_devs = log_recruit_devs,
        p1 = p1,
        p2 = p2,
        p3 = p3,
        p4 = p4,
        p5 = p5,
        p6 = p6,
        pfi_50 = pfi_50,
        pfi_slope = pfi_slope
      )
    ) # end of output list

    return(output)
  } # end of create_em_input function
