#' @title run_OM
#'
#' @description function to create operating model
#'
#' @param nyear_forward number of years to run operating model
#' @param OM_list list of outputs from initial operating model
#' @param lh_list data list for creating operating model, includes data from create_OM_input()
#' @param EM_res results from estimation model
#' @param SSB_trend spawning biomass trend obtained from estimation model
#' increase", "increase_cont", "decrease", "decrease_cont"
#' @param mgmt_rule rule on how to implement artificial reefs
#' AR_ref_1 - place artificial reefs anywhere (goes by TAC)
#' AR_ref_2 - place artificial reefs near historica ones (goes by TAC)
#' AR_strategy - use SSB trend from estimation model to place artificial reefs
#'    increase - place near historic artificial reefs
#'    increase_cont - place in areas of high effort
#'    decrease - place in areas of high recruitment
#'    decrease_cont - place in areas of high recruitment and no fishing at 10% of artificial reefs
#'

run_future_OM <- 
function(seed,
        nyear_forward = 5,
        OM_list = OM_list,
        lh_list = lh_list,
        EM_res = EM_res,
        SSB_trend,
        mgmt_rule)
{

  with(OM_list, {
    #####################
    ## Control Panel ####
    #####################
    set.seed(seed)

    nyears <- settings$nyears
    nsites <- settings$nsites
    Amax <- settings$Amax
    nlandings <- settings$nlandings

    nyears_beg <- nyears+1
    nyears_end <- nyears+nyear_forward
    nyears_future <- nyears_beg:nyears_end
    OM_years <- 1:nyears

    # AR sites
    AR_sites_old <- spatial$AR_sites
    NR_sites <- spatial$NR_sites
    AR_coord_file_old <- spatial$AR_coord_file
    # combining old and new
    AR_new <- list()
    AR_comp <- list()
    # filling in from previous years
    AR_num_new <- spatial$AR_num # compiling ARs each year (not total)
    AR_sites_new <- spatial$AR_sites # compiling ARs sites each year (not total)    
    # areas of high effort and recruitment
    eff_space_old <- colMeans(fishery$eff_space)
    recruit_space_old <- colMeans(time_dyn$Recruit_save)
    # AR fishery closures
    ARs_fishery_closure_old <- spatial$ARs_fishery_closure

    for(y in nyears_future){
      AR_new[[y]] <- AR_rule(coord_file = spatial$coord_file,
                            AR_coord_file = AR_coord_file_old,
                            ARs_fishery_closure = ARs_fishery_closure_old,
                            SSB_trend = "decrease_cont",
                            eff_space = eff_space_old,
                            recruit_space = recruit_space_old,
                            mgmt_rule = mgmt_rule)
      # fishery closure (if SSB trend keeps decreasing)
      ARs_fishery_closure_old <- AR_new[[y]]$ARs_fishery_closure
      if(mgmt_rule == "AR_strategy" & SSB_trend == "decrease_cont") AR_coord_file_old <- AR_coord_file_old[-ARs_fishery_closure_old$num,]
      AR_coord_file_old <- rbind(AR_coord_file_old, AR_new[[y]]$new_ARs_coords)
      AR_num_new <- convert_AR_coords(nsites = nsites, AR_coord_file = AR_coord_file_old, coord_file = spatial$coord_file)
      AR_sites_new <- which(AR_num_new>0)
      AR_comp[[y]] <- AR_num_new/sum(AR_num_new)
    }

    AR_coord_file_new <- AR_coord_file_old # compile all AR coordinates
    ARs_fishery_closure_new <- ARs_fishery_closure_old

    # NR sites
    NR <- spatial$coord_file$HB
    NR_comp <- NR/sum(NR)

    # recruitment
    lR0_FL <- recruit$log_R0

    ## errors
    IndexDev <- CatchDev <- SurveyDev <- vector(length = nyears_end)
    RecDev <- vector(length = nyears_end)
    RecDev[nyears_future] <- rnorm(nyear_forward, mean = -(lh_list$errors$SigmaR^2)/2, sd = lh_list$errors$SigmaR)
    IndexDev[nyears_future] <- rnorm(nyear_forward, mean = -(lh_list$errors$SigmaI^2)/2, sd = lh_list$errors$SigmaI)
    CatchDev[nyears_future] <- rnorm(nyear_forward, mean = -(lh_list$errors$SigmaC^2)/2, sd = lh_list$errors$SigmaC)
    SurveyDev[nyears_future] <- rnorm(nyear_forward, mean = 0, sd = lh_list$errors$SigmaFI)


    ## cost - for effort allocation
    # cost_null <- 50 + mvt$dist*0 # this represents low differences in costs between sites so HIGH movement of anglers (all sites have same cost)
    cost_mult <- 50 + mvt$dist*3 # intermediate # not yet exported
    # cost_exp <- 5 + mvt$dist^1.5 # this represents great differences in costs between sites, so LOW movement # not yet exported
    cost <- cost_mult


    ## Influence of ARs and NRs
    # assuming more ARs, better benefits
    # assuming ARs affect all ages equally
    AR_infl <- matrix(NA, nrow = nyears_end, ncol = nsites)
      for(y in nyears_future) AR_infl[y,] <- AR_comp[[y]]
    AR_infl_array <- array(NA, c(nyears_end, Amax, nsites))
      for(y in nyears_future) for(j in 1:Amax) AR_infl_array[y,j,] <- AR_comp[[y]]
    # natural mortality
    M_site <- array(NA, c(nyears_end, Amax, nsites))
      M_site[OM_years,,] <- time_dyn$M_site_save
    for(y in nyears_future) for(j in 1:nsites) M_site[y,,j] <- lh_list$M
    # catchability
    q_t <- matrix(0, nrow = nyears_end, ncol = nsites)
      q_t[OM_years,] <- fishery$qt_save
      q_t[nyears_future,] <- lh_list$q 

    if(lh_list$reef_flags$AR$AR_flag == 1){
      for(y in nyears_future){
        M_site[y,,AR_sites_new] <- M_site[y,,AR_sites_new]*(1+lh_list$reef_flags$AR$AR_prop_M)*(1+AR_infl_array[y,,AR_sites_new])
        q_t[y,AR_sites_new] <- q_t[y,AR_sites_new]*(1+lh_list$reef_flags$AR$AR_prop_q)*(1+AR_infl[y,AR_sites_new])
      }
    }

    if(lh_list$reef_flags$NR$NR_flag == 1){
      M_site[nyears_future,,NR_sites] <- M_site[nyears_future,,NR_sites]*(1+lh_list$reef_flags$NR$NR_prop_M)
      q_t[nyears_future,NR_sites] <- q_t[nyears_future,NR_sites]*(1+lh_list$reef_flags$NR$NR_prop_q)
    }



    #############################
    ## Recruitment processes ####
    #############################
    So_site <- exp(-M_site) # survival

    # survivorship
    S_site <- array(NA, c(nyears_end, Amax, nsites))
      S_site[,1,] <- 1
     for(a in 2:Amax) S_site[,a,] <- S_site[,a-1,]*So_site[,a-1,]
      S_site[,Amax,] <- S_site[,Amax,]/(1-exp(-M_site[,Amax,]))
    N0 <- exp(lR0_FL)*S_site
    SSB0 <- recruit$SSB0



    ##################################################################
    ## Creating empty vectors, matrices, arrays and initilization ####
    ##################################################################
    Nage <- array(0, dim = c(nyears_end, Amax, nsites)) # Actual population-at-age in each cell for each year before mortality
     Nage[OM_years,,] <- time_dyn$Nage_save
    Recruit <- matrix(NA, nrow = nyears_end, ncol = nsites)
      Recruit[OM_years,] <- time_dyn$Recruit_save
    SSB <- array(0, dim =nyears_end) # spawning biomass
      SSB[OM_years] <- time_dyn$SSB_save
    VB <- matrix(0, nrow = nyears_end, ncol = nsites) # vulnerable/exploitable biomass (for now used in gravity model)
      VB[OM_years,] <- time_dyn$VB_save
    Catch_bio <- matrix(0, nrow = nyears_end, ncol = nsites) # total catch (biomass)
      Catch_bio[OM_years,] <- fishery$Catch_bio_save
    Catch_num <- matrix(0, nrow = nyears_end, ncol = nsites) # catch at age
      Catch_num[OM_years,] <- fishery$Catch_num_save
    Catch_numage <- array(0, dim = c(nyears_end, Amax, nsites)) # total catch
      Catch_numage[OM_years,,] <- fishery$Catch_numage_save
    FM <- array(0, dim = c(nyears_end, Amax, nsites)) # fishing mortality
      FM[OM_years,,] <- fishery$FM_save
    pmove <- array(0, dim = c(nsites, nsites, Amax)) # movement
    pfish <- array(0, dim = c(nlandings, nsites, nyears_end)) # probability of fishing
    eff_space <- matrix(0, nrow = nyears_end, ncol = nsites) # effort expended
      eff_space[OM_years,] <- fishery$eff_space_save
    pre_season <- vector(length = nyears_end)
      pre_season[OM_years] <- fishery$pre_season_save
    exp_eff <- vector(length = nyears_end)
    eff_init <- vector(length = nsites)
    pmax_eff <- vector(length = nyears_end)
      pmax_eff[OM_years] <- fishery$pmax_eff_save

    # Utiltiy metrics for effort dynamic
    CPUE_Ut <- HPUE_Ut <- Size_Ut <- matrix(NA, nrow = nyears_end, ncol = nsites)
    U_cpue <- U_hpue <- U_size <- U_dist <- U_crd <- Tot_Ut <- matrix(NA, nrow = nyears_end, ncol = nsites) 
    Tot_Ut_sum <- NULL
    persis_pmax_eff <- NULL
    grav_wt <- matrix(NA, nrow = nyears_end, ncol = nsites)

    # effort dynamic
    tot_eff <- fishery$tot_eff_save

    # season length
    # derby effort
    int_eff <- as.numeric(spatial$RS_season_file[1])
    slope_eff <- as.numeric(spatial$RS_season_file[2])
    CF_eff <- as.numeric(spatial$RS_season_file[3])
    # theoretical daily effort
    eff_season <- (exp(int_eff)+exp(slope_eff)*1:365)*CF_eff 
    eff_cap <- pmin(eff_season, lh_list$eff_param$max_eff)

    # utility metrics
    Tot_Ut[OM_years,] <- fishery$Tot_Ut_save
    Tot_Ut_sum <- NA
      Tot_Ut_sum[nyears_beg-1] <- sum(Tot_Ut[nyears_beg-1,])
    Tot_init <- fishery$Tot_init_save

    ## Movement
    sub_mod <- list()
    hab_suit_per <- list()
    pred_data <- mvt$coord_dat
    for(y in nyears_future){
      pred_data$AR <- AR_comp[[y]]
      sub_mod[[y]] <- lm(RS_HS ~ AR+HB, data = pred_data)
      hab_suit_per[[y]] <- predict(sub_mod[[y]])
    }
    depth <- mvt$depth

    ## calculate new selectivity according to estimation model
    sel_list <- list(p1 = EM_res$param_df$estimate[which(EM_res$param_df$parameter == "p1")],
                    p2 = EM_res$param_df$estimate[which(EM_res$param_df$parameter == "p2")],
                    p3 = EM_res$param_df$estimate[which(EM_res$param_df$parameter == "p3")],
                    p4 = EM_res$param_df$estimate[which(EM_res$param_df$parameter == "p4")],
                    p5 = EM_res$param_df$estimate[which(EM_res$param_df$parameter == "p5")],
                    p6 = EM_res$param_df$estimate[which(EM_res$param_df$parameter == "p6")])
    # peak - endpoint where selectivity = 1.0
      peak <- sel_list$p1+lh_list$binwidth+((0.99*Amax-sel_list$p1-lh_list$binwidth)/(1+exp(-sel_list$p2)))
    # joiner functions for asc and desc components
      jf1 <- (1+exp(-20*(lh_list$Ages-sel_list$p1)/(1+abs(lh_list$Ages-sel_list$p1))))^-1
      jf2 <- (1+exp(-20*(lh_list$Ages-peak)/(1+abs(lh_list$Ages-peak))))^-1
    # t1 and t2
      t1 <- exp(-(lh_list$Ages[1]-sel_list$p1)^2/exp(sel_list$p3))
      t2 <- exp(-(Amax-peak)^2/exp(sel_list$p4))
    # asc and desc portions
     asc <- (1+exp(-sel_list$p5))^-1+(1-(1+exp(-sel_list$p5))^-1)*(exp(-(lh_list$Ages-sel_list$p1)^2/exp(sel_list$p3))-t1)/(1-t1)
     desc <- 1+(((1+exp(-sel_list$p6))^-1)-1)*(exp(-(lh_list$Ages-peak)^2/exp(sel_list$p4))-1)/(t2-1)
    # selectivity by age
      sel_a <- asc*(1-jf1)+jf1*((1-jf2)+jf2*desc)
       sel_a <- sel_a/max(sel_a)
    Sel <- matrix(sel_a, nrow = Amax, ncol = nsites)

    ## TAC mgmt rule - calculates FM
    TAC_catch <- TAC_calc(Fref = EM_res$Derived$F26,
                          Depletion = EM_res$Report$Depletion,
                          SSB0 = EM_res$Derived$SSB0,
                          SSB_curr = EM_res$Report$SSB,
                          VB = EM_res$Report$VB, 
                          year = (nyears-nyear_forward+1):nyears,
                          Ages = lh_list$Ages,
                          R0 = exp(EM_res$param_df$estimate[which(EM_res$param_df$parameter == "log_R0")]),
                          M_a = EM_res$EM_inputs$data$M_vec,
                          W_a = lh_list$W_a,
                          sel_a = sel_a,
                          alpha = 0.1, 
                          beta = 0.26)
    Ntemp <- TAC_catch$Ntemp
    VBtemp <- EM_res$Report$VB
    TAC <- TAC_catch$TAC
    TAC_flag <- TAC_catch$TAC_flag

    # calculate FM according to TAC
    temp <- TAC/(VBtemp+0.1*TAC)
    join <- 1/(1+exp(30*(temp-0.95)))
    temp2 <- join*temp+0.95*(1-join)
    Fnew <- -log(1-temp2[nyears_future-nyear_forward])

    # Hybrid method for estimating F from catch
    for(f in 1:4){
      for(y in 1:length(Fnew)){
        Z <- lh_list$M+Fnew[y]*sel_a
        alpha <- (1-exp(-Z))
        Ctemp <- sum((Fnew[y]/Z)*(Ntemp*lh_list$W_a*sel_a)*alpha)

        Zadj <- TAC[y]/(Ctemp+1e-4)

        Zprime <- lh_list$M+Zadj*(Z-lh_list$M)
        alpha <- (1-exp(-Zprime))/Zprime

        temp <- sum(Ntemp*lh_list$W_a*sel_a*alpha)
        Ftemp <- TAC[y]/(temp+1e-4)

        j2 <- 1/(1+exp(30*(Ftemp-0.95*1)))

        Fnew[y] <- j2*Ftemp+(1-j2)
      }
    }
    Fnew2 <- matrix(NA, nrow = nyears_end, ncol = Amax)
    Fnew2[nyears_future,] <- sapply(1:length(Fnew), function(x) Fnew[x])


    #####################
    ## Time Dynamics ####
    #####################   
    for(y in nyears_beg:nyears_end){
      ## Dynamic effort
      for(k in 1:nsites) grav_wt[y,k] <- (Tot_Ut[y-1,k]/cost[k])^lh_list$eff_param$w_cost
      pmax_eff[y] <- 1/(1+exp(-(Tot_Ut_sum[y-1]-Tot_init)/(lh_list$sig1e*Tot_init)))
      persis_pmax_eff[y] <- pmax_eff[y]*(1-lh_list$persis)+pmax_eff[y-1]*lh_list$persis # accounts for some "stickyness" in the effort time step to time step
      # season length
      for(k in 1:nsites) eff_init[k] <- (tot_eff*persis_pmax_eff[y]*grav_wt[y,k])/sum(grav_wt[y,])
      dens_eff <- dtruncnorm(eff_season, a = 0, mean = sum(eff_init), sd = sum(eff_init)*0.3)
      dens_eff <- dens_eff/sum(dens_eff) # needs to sum to 1 to prevent loss from density range
      exp_eff[y] <- sum(dens_eff*eff_cap) # expected "derby" effort
      pre_season[y] <- floor(uniroot(f = function(x){
                            (exp(int_eff)+exp(slope_eff)*x)*CF_eff - exp_eff[y]
                        }, lower = 0, upper = 365, extendInt = "yes")$root)
      # effort
      for(k in 1:nsites) eff_space[y,k] <- (tot_eff*persis_pmax_eff[y]*grav_wt[y,k])/sum(grav_wt[y,])

      # fishing mortality
      FM[y,,] <- q_t[y,]*t(t(Sel)*eff_space[y,])
        if(y == nyears_beg){
          FM[y,,] <- q_t[y,]*t(t(Sel)*eff_space[y,])
        } else {
          if(TAC_flag) FM[y,,] <- Fnew2[y,]*sel_a
        }

      ## Recruitment
      SSB[y] <- sum(rowSums(Nage[y-1,,]*lh_list$agefec))
      Recruit[y,] <- Nage[y,1,] <- (4*lh_list$h*exp(lR0_FL)*SSB[y])/(SSB0*(1-lh_list$h)+SSB[y]*(5*lh_list$h-1))*(mvt$pref_depth_std[round(depth),1]*hab_suit_per[[y]])/sum(mvt$pref_depth_std[round(depth),1]*hab_suit_per[[y]])*exp(RecDev[y])

      ## Mortality
      Nage[y,2:Amax,] <- Nage[y-1,1:(Amax-1),]*exp(-(M_site[y-1,1:(Amax-1),]+FM[y-1,1:(Amax-1),]))
      # plus group
      Nage[y,Amax,] <- Nage[y,Amax,] + Nage[y-1,Amax,]*exp(-(M_site[y-1,Amax,]+FM[y-1,Amax,]))

      ## Movement
      ind_vec <- ifelse(colSums(t(t(Nage[y,,])*lh_list$L_a^2)) > quantile(colSums(t(t(Nage[1,,])*lh_list$L_a^2)), 0.75), 1, 0)
      for(a in 1:Amax){
        pmove[,,a] <- t(t(exp(-mvt$lam*mvt$dist_mat)^lh_list$mvt_param$w_dist) * (mvt$pref_depth_std[round(depth),a]^lh_list$mvt_param$w_dep) * (hab_suit_per[[y]]^lh_list$mvt_param$w_hab) * ((1/colSums(t(t(Nage[y,,]*lh_list$L_a^2))/quantile(colSums(t(t(Nage[y,,])*lh_list$L_a^2)), 0.75))^0.5)^ind_vec))
        pmove[,,a] <- pmove[,,a]/rowSums(pmove[,,a])
      }
      for(a in 2:Amax) Nage[y,a,] <- colSums(t(sapply(1:nsites, function(x) rmultinom(n = 1, size = round(Nage[y,a,x]), prob = pmove[x,,a]))))

      ## Metrics
      VB[y,] <- colSums(Sel*Nage[y,,]*lh_list$W_a)
      Catch_bio[y,] <- colSums((FM[y,,]/(FM[y,,]+M_site[y,,]))*Nage[y,,]*(1-exp(-(M_site[y,,]+FM[y,,])))*lh_list$W_a)
      Catch_num[y,] <- colSums((FM[y,,]/(FM[y,,]+M_site[y,,]))*Nage[y,,]*(1-exp(-(M_site[y,,]+FM[y,,]))))
      Catch_numage[y,,] <- (FM[y,,]/(FM[y,,]+M_site[y,,]))*Nage[y,,]*(1-exp(-(M_site[y,,]+FM[y,,])))

      ## Utility calcs
      CPUE_Ut[y,] <- Catch_num[y,]/eff_space[y,]
        CPUE_Ut[is.na(CPUE_Ut)] <- 0
      HPUE_Ut[y,] <- colSums(Catch_numage[y,,]*lh_list$ret)/eff_space[y,]
        HPUE_Ut[is.na(HPUE_Ut)] <- 0
      Size_Ut[y,] <- colSums(Nage[y,,]*lh_list$L_a*sel_a)/colSums(Nage[y,,]*sel_a)
       Size_Ut[is.na(Size_Ut)] <- 0
      U_cpue[y,] <- lh_list$Ut_list$cpue[1]/(1+exp(lh_list$Ut_list$cpue[2]*(lh_list$Ut_list$cpue[3]-CPUE_Ut[y,])))
      U_hpue[y,] <- lh_list$Ut_list$hpue[1]/(1+exp(lh_list$Ut_list$hpue[2]*(lh_list$Ut_list$hpue[3]-HPUE_Ut[y,])))
      U_size[y,] <- lh_list$Ut_list$size[1]/(1+exp(lh_list$Ut_list$size[2]*(lh_list$Ut_list$size[3]-Size_Ut[y,])))
      U_dist[y,] <- lh_list$Ut_list$dist[1]/(1+exp(lh_list$Ut_list$dist[2]*(lh_list$Ut_list$dist[3]-mvt$dist))) 
      U_crd[y,] <- lh_list$Ut_list$crowd[1]/(1+exp(lh_list$Ut_list$crowd[2]*(lh_list$Ut_list$crowd[3]-eff_space[y,])))
      Tot_Ut[y,] <- U_cpue[y,]+U_hpue[y,]+U_size[y,]+U_dist[y,]+U_crd[y,]
      Tot_Ut_sum[y] <- sum(Tot_Ut[y,])
    } # end of y loop



    #################################
    ## Data for Estimation model ####
    #################################         
    ## catch data
    Catch_obs <- vector(length = nyears_end)
      Catch_obs[OM_years] <- EM_data$Catch_obs_save
      Catch_obs[nyears_future] <- rowSums(Catch_bio[nyears_future,])*exp(CatchDev[nyears_future])

    ## catch index
    CPUE_index <- vector(length = nyears_end) 
      CPUE_index[OM_years] <- EM_data$CPUE_index_save
      CPUE_index[nyears_future] <- apply(Catch_bio[nyears_future,]/eff_space[nyears_future,], 1, sum, na.rm=TRUE)*exp(IndexDev[nyears_future])

    ## observed effort
    Eff_obs <- rnorm(nyears_end, mean = rowSums(eff_space), sd = lh_list$errors$CVeff*rowSums(eff_space))

    ## age comps
    Age_comp_temp <- array(NA, c(nyears_end, Amax, nsites))
    for(y in nyears_future){
      for(k in 1:nsites){
        Age_comp_temp[y,,k] <- rmultinom(n = 1, size = lh_list$SS_trip, prob = Catch_numage[y,,k]/sum(Catch_numage[y,,k]))
      }
    }
    Age_comp_temp <- apply(X = Age_comp_temp, MARGIN = c(1,2), FUN = sum)
      Age_comp_temp <- Age_comp_temp/rowSums(Age_comp_temp)
    Age_comp <- matrix(NA, nrow = nyears_end, ncol = Amax)
      Age_comp[OM_years,] <- EM_data$Age_comp_save
      Age_comp[nyears_future,] <- Age_comp_temp[nyears_future,]

    # age comps - FI
    Age_comp_FI_temp <- array(NA, c(nyears_end, Amax, nsites))
    for(y in nyears_future){
      for(k in 1:nsites){
        Age_comp_FI_temp[y,,k] <- rmultinom(n = 1, size =  lh_list$SS_trip, prob = (lh_list$S_FI*lh_list$FI_q*Nage[y,,k])/sum(lh_list$S_FI*lh_list$FI_q*Nage[y,,k]))
      }
    }
    Age_comp_FI_temp <- apply(X = Age_comp_FI_temp, MARGIN = c(1,2), FUN = sum)
      Age_comp_FI_temp <- Age_comp_FI_temp/rowSums(Age_comp_FI_temp)
    Age_comp_FI <- matrix(NA, nrow = nyears_end, ncol = Amax)
      Age_comp_FI[OM_years,] <- EM_data$Age_comp_FI_save
      Age_comp_FI[nyears_future,] <- Age_comp_FI_temp[nyears_future,]        



    #################################
    ## Out of loop calculations ####
    #################################
    Tot_Ut_save <- Tot_Ut[1:nyears_end,]
    SSB.SSB0 <- SSB[1:nyears_end]/SSB0
    SSB_save <- SSB[1:nyears_end]
    VB_save <- VB[1:nyears_end,]
    Nage_save <- Nage[1:nyears_end,,]
    Recruit_save <- Recruit[1:nyears_end,]
    settings$nyears <- nyears_end


    ###################
    ## output list ####
    ###################
    output <- NULL
      output<- list(settings = settings,
                    mvt = mvt,
                    spatial = list(AR_sites = AR_sites_new,
                                  NR_sites = NR_sites,
                                  AR_num = AR_num_new,
                                  NR_num = NR,
                                  coord_file = spatial$coord_file,
                                  AR_coord_file = AR_coord_file_new,
                                  RS_season_file = spatial$RS_season_file,
                                  ARs_fishery_closure = ARs_fishery_closure_new),
                    recruit = recruit,
                    fishery = list(qt_save = q_t,
                                  Eff_obs_save = Eff_obs,
                                  eff_space_save = eff_space,
                                  tot_eff_save = tot_eff,
                                  FM_save = FM,
                                  Tot_Ut_save = Tot_Ut_save,
                                  Tot_init_save = Tot_init,
                                  Catch_num_save = Catch_num,
                                  Catch_bio_save = Catch_bio,
                                  Catch_numage_save = Catch_numage,
                                  pre_season_save = pre_season,
                                  pmax_eff_save = pmax_eff),
                    EM_data = list(Catch_obs_save = Catch_obs,
                                  CPUE_index_save = CPUE_index,
                                  Age_comp_save = Age_comp,
                                  Age_comp_FI_save = Age_comp_FI),
                    time_dyn = list(SSB.SSB0_save = SSB.SSB0,
                                    SSB_save = SSB_save,
                                    VB_save = VB_save,
                                    Nage_save = Nage_save,
                                    M_site_save = M_site,
                                    Recruit_save = Recruit_save)
      ) # end of output list
    return(output)
  }) # end of with function
} # end run_OM function