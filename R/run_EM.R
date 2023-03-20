#' @title run_EM
#'
#' @description function to run estimation model (stock assessment) - age structured assessment model (ASAM)
#' 
#' @param sim_data data from operating model
#' @param lh_list life history inputs (same as ones put into operating model)
#' @param fixed_list list of fixed parameters for ASAM
#' @param rand_list list of parameters for random effects
#' @param years_run subset years in operating model to run estimation model on
#'

run_EM <- 
function(sim_data = sim_data,
        lh_list = lh_list,
        fixed_list,
        rand_list,
        years_run)
{
  # create inputs for ASAM
  EM_inputs <- create_EM_data(sim_data = sim_data, lh_list = lh_list, years_run = years_run)

  # create objects to save best results
  obj_save <- NULL
  JNLL <- NULL
  opt_save <- NULL
  opt_save[["final_gradient"]] <- NA

  ## first run
  obj <- TMB::MakeADFun(EM_inputs$data, EM_inputs$parameters, DLL = "ASAM", map = fixed_list, random = rand_list, inner.control=list(maxit=1e3), silent = TRUE, hessian = FALSE)



  ################
  ## Settings ####
  ################
  lower <- obj$par
    lower[names(lower) == "h"] <- 0
    lower[names(lower) == "log_R0"] <- 10
    lower[names(lower) == "log_q"] <- -20
    lower[names(lower) == "log_FI_q"] <- -20
    lower[names(lower) == "log_SigR"] <- -5
    lower[names(lower) == "log_SigC"] <- -5
    lower[names(lower) == "log_CV_I"] <- -5
    lower[names(lower) == "log_theta"] <- -5
    lower[names(lower) == "log_theta_FI"] <- -5
    lower[names(lower) == "log_recruit_devs"] <- rep(-10,EM_inputs$data$Nyear)
    lower[names(lower) == "log_Fint"] <- rep(-20, EM_inputs$data$Nyear)
    lower[names(lower) == "p1"] <- -10
    lower[names(lower) == "p2"] <- -10
    lower[names(lower) == "p3"] <- -10
    lower[names(lower) == "p4"] <- -10
    lower[names(lower) == "p5"] <- -20
    lower[names(lower) == "p6"] <- -10
    lower[names(lower) == "pFI_50"] <- 1e-3
    lower[names(lower) == "pFI_slope"] <- 1e-3
  upper <- obj$par
    upper[names(upper) == "h"] <- 1
    upper[names(upper) == "log_R0"] <- 25 
    upper[names(upper) == "log_q"] <- 1
    upper[names(upper) == "log_FI_q"] <- 1
    upper[names(upper) == "log_SigR"] <- 2
    upper[names(upper) == "log_SigC"] <- 2
    upper[names(upper) == "log_CV_I"] <- 2
    upper[names(upper) == "log_theta"] <- 5
    upper[names(upper) == "log_theta_FI"] <- 5
    upper[names(upper) == "log_recruit_devs"] <- rep(0,EM_inputs$data$Nyear)
    upper[names(upper) == "log_Fint"] <- rep(0,EM_inputs$data$Nyear)
    upper[names(upper) == "p1"] <- 5
    upper[names(upper) == "p2"] <- 5
    upper[names(upper) == "p3"] <- 5
    upper[names(upper) == "p4"] <- 5
    upper[names(upper) == "p5"] <- 5
    upper[names(upper) == "p6"] <- 5
    upper[names(upper) == "pFI_50"] <- 21
    upper[names(upper) == "pFI_slope"] <- 10

  # new parameters if need to run again
  new_params <- list(h = jitter(EM_inputs$parameters$h, 10),
                    log_R0 = jitter(EM_inputs$parameters$log_R0, 10),
                    log_q = jitter(EM_inputs$parameters$log_q, 10),
                    log_FI_q = jitter(EM_inputs$parameters$log_FI_q, 10),
                    log_SigR = jitter(EM_inputs$parameters$log_SigR, 10),
                    log_SigC = EM_inputs$parameters$log_SigC,
                    log_CV_I = EM_inputs$parameters$log_CV_I,
                    log_theta = jitter(EM_inputs$parameters$log_theta, 10),
                    log_theta_FI = jitter(EM_inputs$parameters$log_theta_FI, 10),
                    log_Fint = jitter(EM_inputs$parameters$log_Fint, 10),
                    log_recruit_devs = jitter(EM_inputs$parameters$log_recruit_devs, 10),
                    p1 = jitter(EM_inputs$parameters$p1, 10),
                    p2 = jitter(EM_inputs$parameters$p2, 10),
                    p3 = jitter(EM_inputs$parameters$p3, 10),
                    p4 = jitter(EM_inputs$parameters$p4, 10),
                    p5 = jitter(EM_inputs$parameters$p5, 10),
                    p6 = jitter(EM_inputs$parameters$p6, 10),
                    pFI_50 = jitter(EM_inputs$parameters$pFI_50,10),
                    pFI_slope = jitter(EM_inputs$parameters$pFI_slope,10))



  #########################
  ## Run the optimizer ####
  #########################
  opt <- try(fit_tmb(obj = obj, startpar = obj$par, upper = upper, lower = lower, newtonstep = 1, quiet = TRUE))
  if(class(opt)[1] != "try-error"){
    JNLL <- obj$report()$JNLL
    opt[["final_gradient"]] <- obj$gr(opt$par)
    opt_save <- opt
    obj_save <- obj
    JNLL_save <- obj_save$report(obj$env$last.par.best)$JNLL
  }

  # loop try to get opt to run
  for(i in 1:5){
    if(class(opt)[1] == "try-error" | is.na(JNLL) | all(is.na(opt_save[["final_gradient"]]))){  
      obj <- TMB::MakeADFun(EM_inputs$data, EM_inputs$parameters, DLL = "ASAM", map = fixed_list, random = rand_list, inner.control=list(maxit=1e3), silent = TRUE, hessian = FALSE)
      opt <- try(fit_tmb(obj = obj, startpar = obj$par, newtonstep = 1, upper = upper, lower = lower, quiet = TRUE))
      JNLL <- obj$report(obj$env$last.par.best)$JNLL
    }
    if(all(is.na(opt)) == FALSE & is.na(JNLL) == FALSE){
      opt[["final_gradient"]] = obj$gr(opt$par)       
      opt_save <- opt
      obj_save <- obj
      JNLL_save <- JNLL
      break
    }
  }

  ## if opt ran
  if(all(is.na(opt_save)) == FALSE){
    for(i in 1:5){
      if(all(is.na(opt_save[["final_gradient"]])==FALSE)){
        if(max(abs(opt_save[["final_gradient"]]))>0.001){
          obj <- TMB::MakeADFun(EM_inputs$data, EM_inputs$parameters, DLL = "ASAM", map = fixed_list, random = rand_list, inner.control=list(maxit=1e3), silent = TRUE, hessian = FALSE)
          opt <- try(fit_tmb(obj = obj, startpar = obj$par, newtonstep = 1, upper = upper, lower = lower, quiet = TRUE))
          JNLL <- obj$report(obj$env$last.par.best)$JNLL
        }
      }
      if(all(is.na(opt))==FALSE & is.na(JNLL)==FALSE){
        if(is.null(JNLL_save)){
          opt[["final_gradient"]] = obj$gr(opt$par)       
          opt_save <- opt
          obj_save <- obj
          JNLL_save <- JNLL
        }
        if(is.null(JNLL_save)==FALSE){
          if(JNLL<=JNLL_save){
            opt[["final_gradient"]] = obj$gr(opt$par)       
            opt_save <- opt
            obj_save <- obj
            JNLL_save <- JNLL
          }
        }
        if(all(is.na(opt_save[["final_gradient"]]))==FALSE){
          if(max(abs(opt_save[["final_gradient"]]))<=0.001) break
        }
      }
    }
    if(all(is.na(opt_save)) == FALSE) df <- data.frame("gradient" = as.vector(opt_save$final_gradient), "parameter" = names(obj_save$par), "estimate" = opt_save$par, "transformed" = exp(opt_save$par))
  }



  #################
  ## To export ####
  #################
  # calculate new selectivity to find F26%
  sel_list <- list(p1 = df$estimate[which(df$parameter == "p1")],
                  p2 = df$estimate[which(df$parameter == "p2")],
                  p3 = df$estimate[which(df$parameter == "p3")],
                  p4 = df$estimate[which(df$parameter == "p4")],
                  p5 = df$estimate[which(df$parameter == "p5")],
                  p6 = df$estimate[which(df$parameter == "p6")])
  F26 <- tryCatch(uniroot(calc_ref, lower = 0, upper = 50, 
                          Ages = lh_list$Ages, M_a = lh_list$M, 
                          sel_list = sel_list, fec_a = lh_list$agefec, 
                          R0 = exp(df$estimate[which(df$parameter == "log_R0")]), 
                          ref = 0.26)$root,
                  error = function(e) NA)
  Report <- tryCatch(obj_save$report(obj_save$env$last.par.best), error = function(x) NA)
  sdReport <- tryCatch(sdreport(obj_save, bias.correct=TRUE), error = function(x) NA)
  Derived <- NULL
    Derived$F26 <- F26
    Derived$SSB0 <- Report$SSB0



  ###################
  ## output list ####
  ###################
  output <- NULL
  output <- list(EM_inputs = EM_inputs,
                param_df = df,
                Report = Report,
                sdReport = sdReport,
                Derived = Derived,
                obj = obj_save,
                opt = opt_save) # end of output list

  return(output)  
} # end of run_EM function