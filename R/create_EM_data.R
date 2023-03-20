#' @title create_EM_data
#'
#' @description function to create data for the estimation model (ASAM in TMB)
#' 
#' @param sim_data data from the operating model
#' @param lh_list operating model inputs
#' @param years_run subset years in operating model to run estimation model on
#'

create_EM_data <- 
function(sim_data = sim_data,
        lh_list = lh_list,
        years_run = years_run)
{
  ## Estimation model data
    # lh inputs
    Nyear <- length(years_run)                                                # total number of years
    Nage <- sim_data$settings$Amax                                            # number of ages

    # data from observation model
    harv_t <- sim_data$EM_data$Catch_obs_save[years_run]              # observed total fishery harvest (catch)
    I_t <- sim_data$EM_data$CPUE_index_save[years_run]                # observed abundance index (CPUE) 
    agecomp_ta <- sim_data$EM_data$Age_comp_save[years_run,]          # observed age compositon
    agecomp_FI_ta <- sim_data$EM_data$Age_comp_FI_save[years_run,]    # observed age composition (fishery-independent - FI)
    SS <- rep(lh_list$Amax*lh_list$SS_trip, Nyear)                            # sample size (for age composition)
    SS_FI <- rep(lh_list$Amax*lh_list$SS_trip, Nyear)                         # sample size (for FI age composition

    ages <- lh_list$Ages                                                      # vector of ages
    L_a <- lh_list$L_a                                                        # length-at-age
    W_a <- lh_list$W_a                                                        # weight-at-age
    M_vec <- lh_list$M                                                        # natural mortality-at-age
    Fec_a <- lh_list$agefec                                                   # fecundity-at-age

  ## Estimation model parameters
    h <-lh_list$h                                         # steepness
    log_R0 <- sim_data$recruit$log_R0                     # unfished recruitment
    log_q <- log(lh_list$q)                               # catchability
    log_FI_q <- log(lh_list$FI_q)                         # catchability (for FI survey)

    # errors
    log_SigR <- log(lh_list$errors$SigmaR)                # recruitment variations
    log_SigC <- log(lh_list$errors$SigmaC)                # observation error in catch/harvest
    log_CV_I <- log(lh_list$errors$SigmaI)                # observation error in CPUE index
    log_theta <- log(lh_list$errors$theta)                # weighting parameter for dirichlet-multinomial
    log_theta_FI <- log(lh_list$errors$theta_FI)          # weighting parameter for dirichlet-multinomial (FI survey)

    log_recruit_devs <- log(rep(0.5, Nyear))              # recruitment deviations
    log_Fint <- log(rep(0.1, Nyear))                      # fishing mortality vector

    # selectivity
    p1 <- lh_list$sel_param$S1                            # selectivity - peak: beginning size for plateau
    p2 <- lh_list$sel_param$S2                            # selectivity - top: width of plateau, as logistic between peak and max size
    p3 <- lh_list$sel_param$S3                            # selectivity - ascending width
    p4 <- lh_list$sel_param$S4                            # selectivity - descending width
    p5 <- lh_list$sel_param$S5                            # selectivity - initial: selectivity at first bin
    p6 <- lh_list$sel_param$S6                            # selectivity - final: selectivity at last bin
    pFI_50 <- lh_list$sel_param$FI_50                     # FI age at 50% selectivity
    pFI_slope <- lh_list$sel_param$FI_slope               # FI slope for selectivity


  ###################
  ## output list ####
  ###################
  output <- list(## lh inputs and data
                data = list(Nyear = Nyear,
                            Nage = Nage,
                            harv_t = harv_t,
                            I_t = I_t,
                            agecomp_ta = agecomp_ta,
                            agecomp_FI_ta = agecomp_FI_ta,
                            SS = SS,
                            SS_FI = SS_FI,
                            ages = ages,
                            L_a = L_a,
                            W_a = W_a,
                            M_vec = M_vec,
                            Fec_a = Fec_a,
                            Amin = min(ages),
                            Amax = max(ages)),
                ## parameters
                parameters = list(h = h,
                                  log_R0 = log_R0,
                                  log_q = log_q,
                                  log_FI_q = log_FI_q,
                                  log_SigR = log_SigR,
                                  log_SigC = log_SigC,
                                  log_CV_I = log_CV_I,
                                  log_theta = log_theta,
                                  log_theta_FI = log_theta_FI,
                                  log_Fint = log_Fint,
                                  log_recruit_devs = log_recruit_devs,
                                  p1 = p1,
                                  p2 = p2,
                                  p3 = p3,
                                  p4 = p4,
                                  p5 = p5,
                                  p6 = p6,
                                  pFI_50 = pFI_50,
                                  pFI_slope = pFI_slope)
  ) # end of output list

  return(output)
} # end of create_EM_data function