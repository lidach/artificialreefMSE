#' @title run_MSE
#'
#' @description function to run MSEs, parallel computing
#'
#' @param MSE_list list that has initial operating model, estimation model and MSE settings
#'

run_MSE <-
  function(MSE_list) {
    with(MSE_list, {
      # compile ASAM
      compile(paste0(EM_set$TMBdir, "/ASAM.cpp"))
      dyn.load(dynlib(paste0(EM_set$TMBdir, "/ASAM")))

      # initial OM
      OM_init <- OM_set$OM_run

      # initial EM
      EM_init <- try(run_EM(
        sim_data = OM_init,
        lh_list = lh_list,
        fixed_list = EM_set$fixed_list,
        rand_list = "log_recruit_devs",
        years_run = 1:OM_init$settings$nyears
      ), silent = TRUE)
      if (class(EM_init)[1] == "try-error") {
        for (i in 1:3) {
          EM_init <- try(run_EM(
            sim_data = OM_init,
            lh_list = lh_list,
            fixed_list = EM_set$fixed_list,
            rand_list = "log_recruit_devs",
            years_run = 1:OM_init$settings$nyears
          ), silent = TRUE)
          if (class(EM_init)[1] != "try-error") break
        }
      }
      # get SSB trend
      SSB <- EM_init$Report$SSB
      SSB <- SSB[(OM_init$settings$nyears - EM_set$EM_skip + 1):OM_init$settings$nyears]
      yr <- 1:EM_set$EM_skip
      SSB_fit <- lm(SSB ~ yr)$coefficients[2]
      if (SSB_fit > 0) {
        SSB_trend_init <- "increase"
      } else if (SSB_fit < 0) {
        SSB_trend_init <- "decrease"
      }
      # Depletion - last year only
      Depletion <- EM_init$Report$Depletion
      Depletion <- Depletion[length(Depletion)]
      if (Depletion < 0.26) {
        SSB_trend_init <- "decrease"
      }


      ################
      ## MSE loop ####
      ################
      EM_run_years <- seq(from = OM_init$settings$nyears, to = OM_init$settings$nyears + MSE_set$MSE_nyear_forward, by = EM_set$EM_skip) # run every 5 years like SEDAR
      OM_list <- OM_init
      EM_res <- EM_init
      SSB_trend <- NA
      SSB_trend[1] <- SSB_trend_init
      seed_vec <- sample(1:1e8, length(EM_run_years), replace = FALSE)

      ## loop over MSE years
      for (em in 2:length(EM_run_years)) {
        OM_list <- run_future_OM(
          seed = seed_vec[em],
          nyear_forward = EM_set$EM_skip,
          OM_list = OM_list,
          lh_list = lh_list,
          EM_res = EM_res,
          SSB_trend = SSB_trend[em - 1],
          mgmt_rule = MSE_set$mgmt_rule
        )
        EM_res <- try(run_EM(
          sim_data = OM_list,
          lh_list = lh_list,
          fixed_list = EM_set$fixed_list,
          rand_list = "log_recruit_devs",
          years_run = 1:OM_list$settings$nyears
        ), silent = TRUE)
        # get SSB trend
        SSB <- EM_res$Report$SSB
        SSB <- SSB[(OM_list$settings$nyears - EM_set$EM_skip + 1):OM_list$settings$nyears]
        yr <- 1:EM_set$EM_skip
        SSB_fit <- lm(SSB ~ yr)$coefficients[2]
        if (SSB_fit < 0) {
          if (SSB_trend[em - 1] == "decrease" | SSB_trend[em - 1] == "decrease_cont") {
            SSB_trend[em] <- "decrease_cont"
          } else {
            SSB_trend[em] <- "decrease"
          }
        } else if (SSB_fit > 0) {
          if (SSB_trend[em - 1] == "increase" | SSB_trend[em - 1] == "increase_cont") {
            SSB_trend[em] <- "increase_cont"
          } else {
            SSB_trend[em] <- "increase"
          }
        }
        # Depletion - only last year
        Depletion <- EM_res$Report$Depletion
        Depletion <- Depletion[length(Depletion)]
        if (Depletion < 0.26) {
          if (SSB_trend[em - 1] == "decrease" | SSB_trend[em - 1] == "decrease_cont") {
            SSB_trend[em] <- "decrease_cont"
          } else {
            SSB_trend[em] <- "decrease"
          }
        }
      } # end of MSE time dynamic



      ##############
      ## Output ####
      ##############
      ## Export values and sd from EM
      sdReport <- tryCatch(sdreport(EM_res$obj, bias.correct = TRUE), error = function(x) NA)
      sdReport <- summary(sdReport)
      rep_values <- rownames(sdReport)
      # SSB
      SSB_tab <- data.frame(value = sdReport[rep_values == "SSB", 1])
      SSB_tab$SE <- sdReport[rep_values == "SSB", 2]
      SSB_tab$min <- SSB_tab$value - 2 * SSB_tab$SE
      SSB_tab$max <- SSB_tab$value + 2 * SSB_tab$SE
      # Depletion
      Depletion_tab <- data.frame(value = sdReport[rep_values == "Depletion", 1])
      Depletion_tab$SE <- sdReport[rep_values == "Depletion", 2]
      Depletion_tab$min <- Depletion_tab$value - 2 * Depletion_tab$SE
      Depletion_tab$max <- Depletion_tab$value + 2 * Depletion_tab$SE
      # Catch
      Harv_tab <- data.frame(value = sdReport[rep_values == "Harv", 1])
      Harv_tab$SE <- sdReport[rep_values == "Harv", 2]
      Harv_tab$min <- Harv_tab$value - 2 * Harv_tab$SE
      Harv_tab$max <- Harv_tab$value + 2 * Harv_tab$SE
      # CPUE_index
      CPUE_index_tab <- data.frame(value = sdReport[rep_values == "CPUE_index", 1])
      CPUE_index_tab$SE <- sdReport[rep_values == "CPUE_index", 2]
      CPUE_index_tab$min <- CPUE_index_tab$value - 2 * CPUE_index_tab$SE
      CPUE_index_tab$max <- CPUE_index_tab$value + 2 * CPUE_index_tab$SE
      # Agecomp
      Agecomp_tab <- data.frame(value = sdReport[rep_values == "Agecomp", 1])
      Agecomp_tab$SE <- sdReport[rep_values == "Agecomp", 2]
      Agecomp_tab$min <- Agecomp_tab$value - 2 * Agecomp_tab$SE
      Agecomp_tab$max <- Agecomp_tab$value + 2 * Agecomp_tab$SE
      # Agecomp_FI
      Agecomp_FI_tab <- data.frame(value = sdReport[rep_values == "Agecomp_FI", 1])
      Agecomp_FI_tab$SE <- sdReport[rep_values == "Agecomp_FI", 2]
      Agecomp_FI_tab$min <- Agecomp_FI_tab$value - 2 * Agecomp_FI_tab$SE
      Agecomp_FI_tab$max <- Agecomp_FI_tab$value + 2 * Agecomp_FI_tab$SE

      ## output list
      output <- NULL
      output$OM <- list(
        Depletion = OM_list$time_dyn$SSB.SSB0_save,
        SSB = OM_list$time_dyn$SSB_save
      )
      output$OM$Data <- list(
        Catch = OM_list$EM_data$Catch_obs_save,
        Index = OM_list$EM_data$CPUE_index_save,
        Age_comp = OM_list$EM_data$Age_comp_save,
        Age_comp_FI = OM_list$EM_data$Age_comp_FI_save,
        Effort = OM_list$fishery$eff_space,
        Recruitment = OM_list$time_dyn$Recruit_save
      )
      output$OM$Spatial <- list(
        AR_sites = OM_list$spatial$AR_sites,
        NR_sites = OM_list$spatial$NR_sites,
        Space = OM_list$spatial$coord_file,
        AR_coords = OM_list$spatial$AR_coord_file
      )
      output$OM$Season_length <- OM_list$fishery$pre_season_save
      output$EM <- list(
        parameters = EM_res$param_df,
        Depletion = Depletion_tab,
        SSB = SSB_tab,
        Catch = Harv_tab,
        Index = CPUE_index_tab,
        Age_comp = Agecomp_tab,
        Age_comp_FI = Agecomp_FI_tab,
        SSB_trend = SSB_trend,
        Ftarget = EM_res$Derived$F26
      )

      return(output)
    }) # end of with function
  } # end of run_MSE function
