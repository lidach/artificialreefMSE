#' @title run_init_OM
#'
#' @description function to create initial operating model
#'
#' @param seed for random seed generation
#' @param nyears_forward run model for how many years (excluding burn in years)
#' @param nyears_init how many years to run burn in period
#' @param coord_file data frame with red snapper habitat suitability, NR, depth, longitude, and latitude coordinates for each cell/site
#' @param AR_coord_file data frame with AR coordinates (up to 2018)
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param RS_tag_file tagging data used for movement parameterization
#' @param RS_data_file numbers at age 0-10 for red snapper, used for movement parameterization
#' @param RS_season_file effort parameters to predict season length (based on data from MRIP)
#' @param lh_list data list for creating operating model, includes data from create_OM_input()
#'

run_init_OM <-
  function(seed = 24,
           nyears_forward = 50,
           nyears_init = 50,
           coord_file,
           AR_coord_file,
           landing_site_file,
           RS_tag_file,
           RS_data_file,
           RS_season_file,
           lh_list) {
    with(lh_list, {
      #####################
      ## Control Panel ####
      #####################
      set.seed(seed)

      ## Temporal parameters
      model_years <- (nyears_init * 2 + 1):(nyears_forward + nyears_init * 2)
      fish_strt <- (nyears_init + 1):model_years[length(model_years)]
      nyears <- length(fish_strt) + nyears_init

      ## Spatial inputs
      nsites <- nrow(coord_file) # number of cells/fishing sites
      nlandings <- nrow(landing_site_file) # number of landing sites
      long <- coord_file$long
      lat <- coord_file$lat
      depth <- coord_file$depth

      # distance for landing sites
      dist <- matrix(NA, nrow = nsites, ncol = nrow(landing_site_file))
      for (j in 1:ncol(dist)) {
        for (i in 1:nsites) {
          xone <- long[i]
          yone <- lat[i]
          dist[i, j] <- earth.dist(landing_site_file$long[j], landing_site_file$lat[j], xone, yone)
        }
      }
      # weighted mean of landing sites to each site
      dist.wt <- apply(dist, 1, function(x) weighted.mean(x, landing_site_file$pop))

      # AR and NR sites
      # convert amount of artificial reefs and their coordinates to spatial cells
      AR <- convert_AR_coords(nsites = nsites, AR_coord_file = AR_coord_file, coord_file = coord_file)
      NR <- coord_file$HB # not amount, habitat suitability measure for each "natural reef", assuming one cell = natural reef

      # which cell number has artificial and/or natural reef and composition % (# reef/# cells)
      AR_sites <- which(AR > 0)
      AR_comp <- AR / sum(AR)
      NR_sites <- which(NR > 0)
      NR_comp <- NR / sum(NR)

      ## R0
      Peast <- 0.23 # Proportion of recruits allocated to east of the Mississippi (per RS assessment)
      R0 <- 1.63e8 # unfished recruitment (SEDAR 52)
      FLsubset <- RS_data_file[RS_data_file$Lon > -85.5 & RS_data_file$Lon < -84.7 & RS_data_file$Lat > 28.8 & RS_data_file$Lat < 29.6 & RS_data_file$Depth < 70, ]
      # Subsetting dataset for east of the Mississippi river to apportion recruitment as the number of FL cells / number of cells east of mississippi
      EastMS <- RS_data_file[RS_data_file$Lon > -89 & RS_data_file$Depth < 70, ]
      # Proportion of recruitment that goes into Florida cells (rough approximation)
      PropR_FL <- dim(FLsubset)[1] / dim(EastMS)[1]
      R0_FL <- R0 * Peast * PropR_FL # unfished recruitment for Florida and spatial area
      lR0_FL <- log(R0_FL)

      ## Process and observation errors
      # recruitment
      RecDev <- rnorm(nyears + 1, mean = -(errors$SigmaR^2) / 2, sd = errors$SigmaR)
      # abundance index
      IndexDev <- c(rep(0, nyears_init), rnorm(length(fish_strt), mean = -(errors$SigmaI^2) / 2, sd = errors$SigmaI))
      # catch
      CatchDev <- c(rep(0, nyears_init), rnorm(length(fish_strt), mean = -(errors$SigmaC^2) / 2, sd = errors$SigmaC))

      ## Cost - for effort allocation
      cost_null <- 50 + dist.wt*0 # this represents low differences in costs between sites so HIGH movement of anglers (all sites have same cost)
      # cost_mult <- 50 + dist.wt * 3 # intermediate # not yet exported
      # cost_exp <- 5 + dist^1.5 # this represents great differences in costs between sites, so LOW movement # not yet exported
      cost <- cost_mult

      ## Influence of ARs and NRs
      # assuming more ARs, better benefits
      AR_infl <- matrix(rep(AR_comp, each = nyears), nrow = nyears, ncol = nsites)
      AR_infl_array <- array(NA, c(nyears, Amax, nsites))
      for (i in 1:nyears) for (j in 1:Amax) AR_infl_array[i, j, ] <- AR_comp
      # natural mortality
      M_site <- array(NA, c(nyears, Amax, nsites))
      for (i in 1:nyears) for (j in 1:nsites) M_site[i, , j] <- M
      # catchability
      q_t <- matrix(0, nrow = nyears, ncol = nsites)
      q_t[fish_strt, ] <- q # no fishing before fish_strt

      if (reef_flags$AR$AR_flag == 1) {
        M_site[, , AR_sites] <- M_site[, , AR_sites] * (1 + reef_flags$AR$AR_prop_M) * (1 + AR_infl_array[, , AR_sites])
        q_t[, AR_sites] <- q_t[, AR_sites] * (1 + reef_flags$AR$AR_prop_q) * (1 + AR_infl[, AR_sites])
      }

      if (reef_flags$NR$NR_flag == 1) {
        M_site[, , NR_sites] <- M_site[, , NR_sites] * (1 + reef_flags$NR$NR_prop_M)
        q_t[, NR_sites] <- q_t[, NR_sites] * (1 + reef_flags$NR$NR_prop_q)
      }



      #############################
      ## Recruitment processes ####
      #############################
      So <- c(1, cumprod(exp(-M))[1:Amax - 1]) # survival
      So[Amax] <- So[Amax] / (1 - exp(-M[Amax]))
      N0 <- exp(lR0_FL) * So
      # unfished spawning biomass
      SSB0 <- sum(N0 * agefec)



      ################
      ## Movement ####
      ################
      ## movement matrix - dist_mat
      dist_mat <- matrix(NA, nrow = nsites, ncol = nsites)
      for (i in 1:nsites) {
        for (j in 1:nsites) {
          dist_mat[i, j] <- earth.dist(long[i], lat[i], long[j], lat[j])
        }
      }

      ## Depth preference
      # aggregate, mean numbers caught at depth and sd
      agg_list <- list()
      for (i in 1:11) {
        agg_list[[i]] <- do.call(data.frame, aggregate(RS_data_file[, 50 + i], by = list(RS_data_file$Depth), FUN = function(x) c(Av = mean(x), sig = sd(x))))
      }

      # fitting mean depth at age and variance using VBF, then assuming those follow a normal
      means <- sapply(agg_list, FUN = function(x) weighted.mean(x[, 1], w = x[, 2]))
      mod_means <- nls(means ~ C * (1 - exp(-r * (0:10 - a))), start = list(C = 80, r = 1, a = 1))
      pred_means <- predict(mod_means)[1:11]
      vars <- sapply(agg_list, FUN = function(x) wtd.var(x[, 1], w = x[, 2]))
      mod_vars <- nls(vars ~ C * (1 - exp(-r * (0:10 - a))), start = list(C = 1000, r = 1, a = 1))
      pred_vars <- predict(mod_vars)[1:11]

      # calculations
      pref_depth <- matrix(NA, ncol = Amax, nrow = ceiling(max(depth)))
      pref_depth_cumd <- matrix(NA, ncol = Amax, nrow = ceiling(max(depth)))

      for (i in 1:length(pred_means)) {
        pref_depth[, i] <- dnorm(seq(1, ceiling(max(depth))), mean = pred_means[i], sd = sqrt(pred_vars[i]))
        for (j in 1:length(pred_means)) {
          pref_depth_cumd[j, i] <- pnorm(j + 0.5, mean = pred_means[i], sd = sqrt(pred_vars[i])) - pnorm(j - 0.5, mean = pred_means[i], sd = sqrt(pred_vars[i]))
        }
      }
      pref_depth[, 12:Amax] <- pref_depth[, 11]

      # every column has max at 1
      pref_depth_std <- t(t(pref_depth) / apply(pref_depth, 2, max))
      row.names(pref_depth_std) <- seq(1, ceiling(max(depth)))

      ## Substrate preference
      coord_dat <- data.frame(
        AR = AR / sum(AR),
        HB = coord_file$HB / sum(coord_file$HB),
        RS_HS = coord_file$RS_HS / sum(coord_file$RS_HS)
      )
      sub_mod <- lm(RS_HS ~ AR + HB, data = coord_dat)
      hab_suit_per <- predict(sub_mod)

      ## Distance preference
      RS_tag <- RS_tag_file[!is.na(RS_tag_file[, 1]), ]
      RS_tag$travdist <- earth.dist(RS_tag[, "Lon1"], RS_tag[, "Lat1"], RS_tag[, "Lon2"], RS_tag[, "Lat2"]) # Calculating the distance a red snapper traveled
      RS_tag$travdist <- RS_tag$travdist * (365 / RS_tag$days_at_large) # Calculating the distance they would have traveled in a full year based on their days at large
      move_RS <- function(theta) {
        lam_dist <- exp(theta[1])
        NLL <- -1 * sum(dexp(RS_tag$travdist[RS_tag$days_at_large > 200], lam_dist, log = TRUE))
        return(NLL)
      }
      move_distcost <- optim(log(0.05), fn = move_RS, method = "BFGS")
      lam <- exp(move_distcost$par)


      ## movement matrix without threshold preference (initial year)
      pmove_init <- array(0, dim = c(nsites, nsites, Amax))
      for (a in 1:Amax) {
        pmove_init[, , a] <- t(t(exp(-lam * dist_mat)^mvt_param$w_dist) * (pref_depth_std[round(depth), a]^mvt_param$w_dep) * (hab_suit_per^mvt_param$w_hab))
        pmove_init[, , a] <- pmove_init[, , a] / rowSums(pmove_init[, , a])
      }



      ################################################
      ## Creating empty vectors, matrices, arrays ####
      ################################################
      Nage <- array(0, dim = c(nyears, Amax, nsites)) # Actual population-at-age in each cell for each year before mortality
      Recruit <- matrix(0, nrow = nyears, ncol = nsites) # recruitment
      SSB <- array(0, dim = nyears) # spawning biomass
      VB <- matrix(0, nrow = nyears, ncol = nsites) # vulnerable/exploitable biomass (for now used in gravity model)
      Catch_bio <- matrix(0, nrow = nyears, ncol = nsites) # total catch (biomass)
      Catch_num <- matrix(0, nrow = nyears, ncol = nsites) # catch at age
      Catch_numage <- array(0, dim = c(nyears, Amax, nsites)) # total catch
      FM <- array(0, dim = c(nyears, Amax, nsites)) # fishing mortality
      pmove <- array(0, dim = c(nsites, nsites, Amax)) # movement
      eff_space <- matrix(c(rep(0, nsites * nyears_init), rep(1, nsites * length(fish_strt))), nrow = nyears, ncol = nsites, byrow = TRUE) # effort expended
      pre_season <- vector(length = nyears) # pre season prediction (based on derby effort)
      exp_eff <- vector(length = nyears) # expected "derby" effort
      eff_init <- vector(length = nsites) # initial effort estimate (recreational effort)
      Sel <- matrix(S_fa, nrow = Amax, ncol = nsites) # selectivity by space

      # Utiltiy metrics for effort dynamic
      CPUE_Ut <- HPUE_Ut <- Size_Ut <- matrix(NA, nrow = nyears, ncol = nsites)
      U_cpue <- U_hpue <- U_size <- U_dist <- U_crd <- Tot_Ut <- matrix(NA, nrow = nyears, ncol = nsites)
      Tot_Ut_sum <- NULL
      pmax_eff <- persis_pmax_eff <- NULL
      grav_wt <- matrix(NA, nrow = nyears, ncol = nsites)



      ##############
      ## Effort ####
      ##############
      # effort dynamic
      tot_eff <- eff_param$max_eff / 2
      pmax_eff[1] <- sum(eff_space[1, ]) / tot_eff
      persis_pmax_eff[1] <- pmax_eff[1]

      # derby effort
      int_eff <- RS_season_file[1]
      slope_eff <- RS_season_file[2]
      CF_eff <- RS_season_file[3]
      # theorectical daily effort
      eff_season <- (exp(int_eff) + exp(slope_eff) * 1:365) * CF_eff
      eff_cap <- pmin(eff_season, tot_eff)


      ############################
      ## Model Initialization ####
      ############################
      Recruit[1, ] <- Nage[1, 1, ] <- exp(lR0_FL) * ((pref_depth_std[round(depth), 1] * hab_suit_per) / sum(pref_depth_std[round(depth), 1] * hab_suit_per))
      Nage[1, , ] <- exp(lR0_FL) * ((pref_depth_std[round(depth), 1] * hab_suit_per) / sum(pref_depth_std[round(depth), 1] * hab_suit_per))

      # initial movement
      for (a in 2:Amax) Nage[1, a, ] <- Nage[1, a, ] %*% pmove_init[, , a]
      VB[1, ] <- colSums(S_fa * Nage[1, , ] * W_a)
      SSB[1] <- sum(rowSums(Nage[1, , ]) * agefec)


      # utility metrics
      CPUE_Ut[1, ] <- colSums(Nage[1, , ] * S_fa)
      HPUE_Ut[1, ] <- colSums(Nage[1, , ] * S_fa * ret)
      Size_Ut[1, ] <- colSums(Nage[1, , ] * S_fa * L_a) / colSums(Nage[1, , ] * S_fa)
      U_cpue[1, ] <- Ut_list$cpue[1] / (1 + exp(Ut_list$cpue[2] * (Ut_list$cpue[3] - CPUE_Ut[1, ])))
      U_hpue[1, ] <- Ut_list$hpue[1] / (1 + exp(Ut_list$hpue[2] * (Ut_list$hpue[3] - HPUE_Ut[1, ])))
      U_size[1, ] <- Ut_list$size[1] / (1 + exp(Ut_list$size[2] * (Ut_list$size[3] - Size_Ut[1, ])))
      U_dist[1, ] <- Ut_list$dist[1] / (1 + exp(Ut_list$dist[2] * (Ut_list$dist[3] - dist.wt)))
      U_crd[1, ] <- Ut_list$crowd[1] / (1 + exp(Ut_list$crowd[2] * (Ut_list$crowd[3] - eff_space[1, ])))
      Tot_Ut[1, ] <- U_cpue[1, ] + U_hpue[1, ] + U_size[1, ] + U_dist[1, ] + U_crd[1, ]
      Tot_Ut_sum[1] <- sum(Tot_Ut[1, ])
      Tot_init <- 1 * Tot_Ut_sum[1]



      #####################
      ## Time Dynamics ####
      #####################
      for (y in 2:nyears) {
        ## Effort dynamics
        for (k in 1:nsites) grav_wt[y, k] <- (Tot_Ut[y - 1, k] / cost[k])^eff_param$w_cost
        pmax_eff[y] <- 1 / (1 + exp(-(Tot_Ut_sum[y - 1] - Tot_init) / (sig1e * Tot_init)))
        persis_pmax_eff[y] <- pmax_eff[y] * (1 - persis) + pmax_eff[y - 1] * persis # accounts for some "stickyness" in the effort time step to time step
        eff_space[y, ] <- 0
        # season length implementation
        if (y >= (fish_strt[1])) {
          for (k in 1:nsites) eff_init[k] <- (tot_eff * persis_pmax_eff[y] * grav_wt[y, k]) / sum(grav_wt[y, ])
          dens_eff <- dtruncnorm(eff_cap, a = 0, mean = sum(eff_init), sd = sum(eff_init) * 0.3)
          dens_eff <- dens_eff / sum(dens_eff) # needs to sum to 1 to prevent loss from density range
          exp_eff[y] <- sum(dens_eff * eff_cap) # expected "derby" effort
          pre_season[y] <- floor(uniroot(f = function(x) {
            (exp(int_eff) + exp(slope_eff) * x) * CF_eff - exp_eff[y]
          }, lower = 0, upper = 365, extendInt = "yes")$root)
          # effort
          for (k in 1:nsites) eff_space[y, k] <- (tot_eff * persis_pmax_eff[y] * grav_wt[y, k]) / sum(grav_wt[y, ])
        }

        ## Fishing mortality
        FM[y, , ] <- q_t[y, ] * t(t(Sel) * eff_space[y, ])

        ## Recruitment
        SSB[y] <- sum(rowSums(Nage[y - 1, , ] * agefec))
        Recruit[y, ] <- Nage[y, 1, ] <- (4 * h * exp(lR0_FL) * SSB[y]) / (SSB0 * (1 - h) + SSB[y] * (5 * h - 1)) * (pref_depth_std[round(depth), 1] * hab_suit_per) / sum(pref_depth_std[round(depth), 1] * hab_suit_per) * exp(RecDev[y])

        ## Mortality
        Nage[y, 2:Amax, ] <- Nage[y - 1, 1:(Amax - 1), ] * exp(-(M_site[y - 1, 1:(Amax - 1), ] + FM[y - 1, 1:(Amax - 1), ]))
        # plus group
        Nage[y, Amax, ] <- Nage[y, Amax, ] + Nage[y - 1, Amax, ] * exp(-(M_site[y - 1, Amax, ] + FM[y - 1, Amax, ]))

        ## Movement
        ind_vec <- ifelse(colSums(t(t(Nage[y, , ]) * L_a^2)) > quantile(colSums(t(t(Nage[1, , ]) * L_a^2)), 0.75), 1, 0)
        for (a in 1:Amax) {
          pmove[, , a] <- t(t(exp(-lam * dist_mat)^mvt_param$w_dist) * (pref_depth_std[round(depth), a]^mvt_param$w_dep) * (hab_suit_per^mvt_param$w_hab) * ((1 / colSums(t(t(Nage[y, , ] * L_a^2)) / quantile(colSums(t(t(Nage[y, , ]) * L_a^2)), 0.75))^0.5)^ind_vec))
          pmove[, , a] <- pmove[, , a] / rowSums(pmove[, , a])
        }
        for (a in 2:Amax) Nage[y, a, ] <- colSums(t(sapply(1:nsites, function(x) rmultinom(n = 1, size = round(Nage[y, a, x]), prob = pmove[x, , a]))))

        ## Metrics
        VB[y, ] <- colSums(Sel * Nage[y, , ] * W_a)
        Catch_bio[y, ] <- colSums((FM[y, , ] / (FM[y, , ] + M_site[y, , ])) * Nage[y, , ] * (1 - exp(-(M_site[y, , ] + FM[y, , ]))) * W_a)
        Catch_num[y, ] <- colSums((FM[y, , ] / (FM[y, , ] + M_site[y, , ])) * Nage[y, , ] * (1 - exp(-(M_site[y, , ] + FM[y, , ]))))
        Catch_numage[y, , ] <- (FM[y, , ] / (FM[y, , ] + M_site[y, , ])) * Nage[y, , ] * (1 - exp(-(M_site[y, , ] + FM[y, , ])))

        ## Utility calcs
        CPUE_Ut[y, ] <- Catch_num[y, ] / eff_space[y, ]
        CPUE_Ut[is.na(CPUE_Ut)] <- 0
        HPUE_Ut[y, ] <- colSums(Catch_numage[y, , ] * ret) / eff_space[y, ]
        HPUE_Ut[is.na(HPUE_Ut)] <- 0
        Size_Ut[y, ] <- colSums(Nage[y, , ] * L_a * S_fa) / colSums(Nage[y, , ] * S_fa)
        Size_Ut[is.na(Size_Ut)] <- 0
        U_cpue[y, ] <- Ut_list$cpue[1] / (1 + exp(Ut_list$cpue[2] * (Ut_list$cpue[3] - CPUE_Ut[y, ])))
        U_hpue[y, ] <- Ut_list$hpue[1] / (1 + exp(Ut_list$hpue[2] * (Ut_list$hpue[3] - HPUE_Ut[y, ])))
        U_size[y, ] <- Ut_list$size[1] / (1 + exp(Ut_list$size[2] * (Ut_list$size[3] - Size_Ut[y, ])))
        U_dist[y, ] <- Ut_list$dist[1] / (1 + exp(Ut_list$dist[2] * (Ut_list$dist[3] - dist.wt)))
        U_crd[y, ] <- Ut_list$crowd[1] / (1 + exp(Ut_list$crowd[2] * (Ut_list$crowd[3] - eff_space[y, ])))
        Tot_Ut[y, ] <- U_cpue[y, ] + U_hpue[y, ] + U_size[y, ] + U_dist[y, ] + U_crd[y, ]
        Tot_Ut_sum[y] <- sum(Tot_Ut[y, ])
      } # end of y loop, time dynamics



      #################################
      ## Data for Estimation model ####
      #################################
      ## fishing mortality
      FM <- FM[model_years, , ]

      ## catch data
      Catch_obs <- rowSums(Catch_bio[model_years, ]) * exp(CatchDev[model_years])

      ## catch index
      CPUE_index <- apply(Catch_bio[model_years, ] / eff_space[model_years, ], 1, sum, na.rm = TRUE) * exp(IndexDev[model_years])

      ## observed effort
      Eff_obs <- rnorm(length(model_years), mean = rowSums(eff_space[model_years, ]), sd = errors$CVeff * rowSums(eff_space[model_years, ]))

      ## age comps
      Age_comp <- array(NA, c(nyears, Amax, nsites))
      for (y in model_years) {
        for (k in 1:nsites) {
          Age_comp[y, , k] <- rmultinom(n = 1, size = SS_trip, prob = Catch_numage[y, , k] / sum(Catch_numage[y, , k]))
        }
      }
      Age_comp <- apply(X = Age_comp[model_years, , ], MARGIN = c(1, 2), FUN = sum)
      Age_comp <- Age_comp / rowSums(Age_comp)

      ## age comps - FI
      Age_comp_FI <- array(NA, c(nyears, Amax, nsites))
      for (y in model_years) {
        for (k in 1:nsites) {
          Age_comp_FI[y, , k] <- rmultinom(n = 1, size = SS_trip, prob = (S_FI * FI_q * Nage[y, , k]) / sum(S_FI * FI_q * Nage[y, , k]))
        }
      }
      Age_comp_FI <- apply(X = Age_comp_FI[model_years, , ], MARGIN = c(1, 2), FUN = sum)
      Age_comp_FI <- Age_comp_FI / rowSums(Age_comp_FI)



      #################################
      ## Out of loop calculations ####
      #################################
      FM_save <- FM
      qt_save <- q_t[model_years, ]
      SSB.SSB0 <- SSB[model_years] / SSB0
      SSB_save <- SSB[model_years]
      VB_save <- VB[model_years, ]
      Nage_save <- Nage[model_years, , ]
      Recruit_save <- Recruit[model_years, ]
      Catch_bio_save <- Catch_bio[model_years, ]
      Catch_num_save <- Catch_num[model_years, ]
      Catch_numage_save <- Catch_numage[model_years, , ]
      RecDev_save <- RecDev[model_years]
      M_site_save <- M_site[model_years, , ]
      eff_space_save <- eff_space[model_years, ]
      Tot_Ut_save <- Tot_Ut[model_years, ]
      pre_season_save <- pre_season[model_years]
      pmax_eff_save <- pmax_eff[model_years]
      CPUE_Ut_save <- CPUE_Ut[model_years, ]
      HPUE_Ut_save <- HPUE_Ut[model_years, ]
      Size_Ut_save <- Size_Ut[model_years, ]

      # for ARs
      AR_sites_list <- list()
      AR_num_list <- list()
      for (y in 1:nyears) {
        AR_sites_list[[y]] <- AR_sites
        AR_num_list[[y]] <- AR
      }



      ###################
      ## output list ####
      ###################
      output <- NULL
      output <- list(
        settings = list(
          seed = seed,
          nyears = length(model_years),
          nsites = nsites,
          nlandings = nlandings,
          Amax = Amax
        ),
        mvt = list(
          dist = dist.wt,
          dist_mat = dist_mat,
          pref_depth_std = pref_depth_std,
          depth = depth,
          coord_dat = coord_dat,
          lam = lam
        ),
        spatial = list(
          AR_sites = AR_sites_list,
          AR_num = AR_num_list,
          NR_sites = NR_sites,
          coord_file = coord_file,
          AR_coord_file = AR_coord_file,
          RS_season_file = RS_season_file,
          ARs_fishery_closure = NULL
        ),
        recruit = list(
          log_R0 = lR0_FL,
          SSB0 = SSB0
        ),
        fishery = list(
          qt_save = qt_save,
          Eff_obs_save = Eff_obs,
          eff_space_save = eff_space_save,
          tot_eff_save = tot_eff,
          FM_save = FM_save,
          Tot_Ut_save = Tot_Ut_save,
          Tot_init_save = Tot_init,
          Catch_num_save = Catch_num_save,
          Catch_bio_save = Catch_bio_save,
          Catch_numage_save = Catch_numage_save,
          pre_season_save = pre_season_save,
          pmax_eff_save = pmax_eff_save
        ),
        EM_data = list(
          Catch_obs_save = Catch_obs,
          CPUE_index_save = CPUE_index,
          Age_comp_save = Age_comp,
          Age_comp_FI_save = Age_comp_FI,
          RecDev_save = RecDev_save
        ),
        time_dyn = list(
          SSB.SSB0_save = SSB.SSB0,
          SSB_save = SSB_save,
          VB_save = VB_save,
          Nage_save = Nage_save,
          M_site_save = M_site_save,
          Recruit_save = Recruit_save
        ),
        utility = list(
          CPUE_Ut_save = CPUE_Ut_save,
          HPUE_Ut_save = HPUE_Ut_save,
          Size_Ut_save = Size_Ut_save
        )
      ) # end of output list

      return(output)
    }) # end with function
  } # end of run_OM function
