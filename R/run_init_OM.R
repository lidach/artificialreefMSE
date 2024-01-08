#' @title run_init_om
#'
#' @description function to create initial operating model
#'
#' @param seed for random seed generation
#' @param nyears_forward run model for how many years (excluding burn in years)
#' @param nyears_init how many years to run burn in period
#' @param coord_file data frame with red snapper habitat suitability, NR, depth, longitude, and latitude coordinates for each cell/site
#' @param ar_coord_file data frame with AR coordinates (up to 2018)
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param rs_tag_file tagging data used for movement parameterization
#' @param rs_data_file numbers at age 0-10 for red snapper, used for movement parameterization
#' @param rs_season_file effort parameters to predict season length (based on data from MRIP)
#' @param lh_list data list for creating operating model, includes data from create_OM_input()
#'

run_init_om <-
  function(seed = 24,
           nyears_forward = 50,
           nyears_init = 50,
           coord_file,
           ar_coord_file,
           landing_site_file,
           rs_tag_file,
           rs_data_file,
           rs_season_file,
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
          dist[i, j] <- earth_dist(landing_site_file$long[j], landing_site_file$lat[j], xone, yone)
        }
      }
      # weighted mean of landing sites to each site
      dist_wt <- apply(dist, 1, function(x) weighted.mean(x, landing_site_file$pop))

      # AR and NR sites
      # convert amount of artificial reefs and their coordinates to spatial cells
      ar <- convert_ar_coords(nsites = nsites, ar_coord_file = ar_coord_file, coord_file = coord_file)
      nr <- coord_file$HB # not amount, habitat suitability measure for each "natural reef", assuming one cell = natural reef

      # which cell number has artificial and/or natural reef and composition % (# reef/# cells)
      ar_sites <- ar$ar_sites
      ar_comp <- ar$ar_comp / sum(ar$ar_comp)
      nr_sites <- which(nr > 0)
      nr_comp <- nr / sum(nr)

      ## ro
      # code from N. Fisch
      p_east <- 0.23 # Proportion of recruits allocated to east of the Mississippi (per RS assessment)
      ro <- 1.63e8   # unfished recruitment (SEDAR 52)
      # subsetting Florida within model area
      FLsubset <- rs_data_file[rs_data_file$Lon > -85.5 & 
                                rs_data_file$Lon < -84.7 & 
                                rs_data_file$Lat > 28.8 & 
                                rs_data_file$Lat < 29.6 & 
                                rs_data_file$Depth < 70, ]
      # Subsetting dataset for east of the Mississippi river to apportion recruitment as the number of FL cells / number of cells east of mississippi
      eastMS <- rs_data_file[rs_data_file$Lon > -89 & rs_data_file$Depth < 70, ]
      # Proportion of recruitment that goes into Florida cells (rough approximation)
      propr_FL <- dim(FLsubset)[1] / dim(eastMS)[1]
      ro_FL <- ro * p_east * propr_FL # unfished recruitment for Florida and spatial area
      lro_FL <- log(ro_FL)

      ## Process and observation errors
      # recruitment
      recdev <- rnorm(nyears + 1, mean = -(sigma_r^2) / 2, sd = sigma_r)
      # abundance index
      indexdev <- c(
        rep(0, nyears_init),
        rnorm(length(fish_strt), mean = -(sigma_i^2) / 2, sd = sigma_i)
      )
      # catch
      catchdev <- c(
        rep(0, nyears_init),
        rnorm(length(fish_strt), mean = -(sigma_c^2) / 2, sd = sigma_c)
      )

      ## Cost - for effort allocation
      cost_null <- 50 + dist_wt * 0   # this represents low differences in costs between sites so HIGH movement of anglers (all sites have same cost)
      # cost_mult <- 50 + dist_wt * 3 # intermediate                                                                                                 # not yet exported
      # cost_exp <- 5 + dist^1.5      # this represents great differences in costs between sites, so LOW movement                                    # not yet exported
      cost <- cost_null

      ## Influence of ARs and NRs
      # assuming more ARs, better benefits
      ar_infl <- matrix(rep(ar_comp, each = nyears), nrow = nyears, ncol = nsites)
      ar_infl_array <- array(NA, c(nyears, amax, nsites))
      for (i in 1:nyears) for (j in 1:amax) ar_infl_array[i, j, ] <- ar_comp
      # natural mortality
      M_site <- array(NA, c(nyears, amax, nsites))
      for (i in 1:nyears) for (j in 1:nsites) M_site[i, , j] <- M
      # catchability
      q_t <- matrix(0, nrow = nyears, ncol = nsites)
      q_t[fish_strt, ] <- q # no fishing before fish_strt
      if (ar_flag == 1) {
        M_site[, , ar_sites] <- M_site[, , ar_sites] * (1 + ar_prop_M) * 
                                                        (1 + ar_infl_array[, , ar_sites])
        q_t[, ar_sites] <- q_t[, ar_sites] * (1 + ar_prop_q) * 
                                              (1 + ar_infl[, ar_sites])
      }
      if (nr_flag == 1) {
        M_site[, , nr_sites] <- M_site[, , nr_sites] * (1 + nr_prop_M)
        q_t[, nr_sites] <- q_t[, nr_sites] * (1 + nr_prop_q)
      }


      #######################
      ## Unfished states ####
      #######################
      so <- c(1, cumprod(exp(-M))[1:amax - 1]) # survival
      so[amax] <- so[amax] / (1 - exp(-M[amax]))
      no <- exp(lro_FL) * so
      # unfished spawning biomass
      ssbo <- sum(no * agefec)


      ################
      ## Movement ####
      ################
      ## movement matrix based on distance from site to site 
      dist_mat <- matrix(NA, nrow = nsites, ncol = nsites)
      for (i in 1:nsites) {
        for (j in 1:nsites) {
          dist_mat[i, j] <- earth_dist(long[i], lat[i], long[j], lat[j])
        }
      }

      ## Depth preference
      # most of this from N. Fisch
      # aggregate, mean numbers caught at depth and sd
      agg_list <- list()
      for (i in 1:11) {
        agg_list[[i]] <- do.call(data.frame, aggregate(rs_data_file[, 50 + i], by = list(rs_data_file$Depth), FUN = function(x) c(Av = mean(x), sig = sd(x))))
      }
      # fitting mean depth at age and variance using von Bert function, then assuming those follow a normal
      means <- sapply(agg_list, FUN = function(x) weighted.mean(x[, 1], w = x[, 2]))
      mod_means <- nls(means ~ C * (1 - exp(-r * (0:10 - a))), start = list(C = 80, r = 1, a = 1))
      pred_means <- predict(mod_means)[1:11]
      vars <- sapply(agg_list, FUN = function(x) wtd_var(x[, 1], w = x[, 2]))
      mod_vars <- nls(vars ~ C * (1 - exp(-r * (0:10 - a))), start = list(C = 1000, r = 1, a = 1))
      pred_vars <- predict(mod_vars)[1:11]
      # calculations for depth preference using mean depth at age
      pref_depth <- matrix(NA, ncol = amax, nrow = ceiling(max(depth)))
      pref_depth_cumd <- matrix(NA, ncol = amax, nrow = ceiling(max(depth)))
      # depth preference following normal distribution
      for (i in 1:length(pred_means)) {
        pref_depth[, i] <- dnorm(seq(1, ceiling(max(depth))), mean = pred_means[i], sd = sqrt(pred_vars[i]))
        for (j in 1:length(pred_means)) {
          pref_depth_cumd[j, i] <- pnorm(j + 0.5, mean = pred_means[i], sd = sqrt(pred_vars[i])) - 
                                    pnorm(j - 0.5, mean = pred_means[i], sd = sqrt(pred_vars[i]))
        }
      }
      # ages 12 to max age have the same preference as age 11
      pref_depth[, 12:amax] <- pref_depth[, 11]
      # every column has max at 1
      pref_depth_std <- t(t(pref_depth) / apply(pref_depth, 2, max))
      row.names(pref_depth_std) <- seq(1, ceiling(max(depth)))

      ## Substrate preference
      # assuming that fish prefer spaces with artificial and natural reefs
      # preference based on # of reefs at site (more reefs -> higher preference)
      coord_dat <- data.frame(
        ar = ar_comp,
        HB = coord_file$HB / sum(coord_file$HB),
        rs_hs = coord_file$RS_HS / sum(coord_file$RS_HS)
      )
      sub_mod <- lm(rs_hs ~ ar + HB, data = coord_dat)
      hab_suit_per <- predict(sub_mod)

      ## Distance preference
      # movement parameter (lam) calculated from tagging data
      rs_tag <- rs_tag_file[!is.na(rs_tag_file[, 1]), ]
      rs_tag$travdist <- earth_dist(rs_tag[, "Lon1"], rs_tag[, "Lat1"], rs_tag[, "Lon2"], rs_tag[, "Lat2"]) # Calculating the distance a red snapper traveled
      rs_tag$travdist <- rs_tag$travdist * (365 / rs_tag$days_at_large) # Calculating the distance they would have traveled in a full year based on their days at large
      move_RS <- function(theta) {
        lam_dist <- exp(theta[1])
        NLL <- -1 * sum(dexp(rs_tag$travdist[rs_tag$days_at_large > 200], lam_dist, log = TRUE))
        return(NLL)
      }
      move_distcost <- optim(log(0.05), fn = move_RS, method = "BFGS")
      lam <- exp(move_distcost$par)
      
      ## Movement matrix without threshold preference (initial year)
      pmove_init <- array(0, dim = c(nsites, nsites, amax))
      for (a in 1:amax) {
        pmove_init[, , a] <- t(t(exp(-lam * dist_mat)^w_dist) * 
                                  (pref_depth_std[round(depth), a]^w_dep) * (hab_suit_per^w_hab))
        pmove_init[, , a] <- pmove_init[, , a] / rowSums(pmove_init[, , a])
      }


      ################################################
      ## Creating empty vectors, matrices, arrays ####
      ################################################
      nage <- array(0, dim = c(nyears, amax, nsites))         # Actual population-at-age in each cell for each year before mortality
      recruit <- matrix(0, nrow = nyears, ncol = nsites)      # recruitment
      ssb <- array(0, dim = nyears)                           # spawning biomass
      vb <- matrix(0, nrow = nyears, ncol = nsites)           # vulnerable/exploitable biomass (for now used in gravity model)
      catch_bio <- matrix(0, nrow = nyears, ncol = nsites)    # total catch (biomass)
      catch_num <- matrix(0, nrow = nyears, ncol = nsites)    # catch at age
      catch_numage <- array(0, dim = c(nyears, amax, nsites)) # total catch
      FM <- array(0, dim = c(nyears, amax, nsites))           # fishing mortality
      pmove <- array(0, dim = c(nsites, nsites, amax))        # movement
      eff_space <- matrix(c(rep(0, nsites * nyears_init), 
                              rep(1, nsites * length(fish_strt))), nrow = nyears, 
                              ncol = nsites, byrow = TRUE)    # effort expended
      pre_season <- vector(length = nyears)                   # pre season prediction (based on derby effort)
      exp_eff <- vector(length = nyears)                      # expected "derby" effort
      eff_init <- vector(length = nsites)                     # initial effort estimate (recreational effort)
      sel <- matrix(s_fa, nrow = amax, ncol = nsites)         # selectivity by space

      # Utiltiy metrics for effort dynamic
      cpue_ut <- hpue_ut <- size_ut <- matrix(NA, nrow = nyears, ncol = nsites)
      u_cpue <- u_hpue <- u_size <- u_dist <- u_crd <- tot_ut <- matrix(NA, nrow = nyears, ncol = nsites)
      tot_ut_sum <- NULL
      pmax_eff <- persis_pmax_eff <- NULL
      grav_wt <- matrix(NA, nrow = nyears, ncol = nsites)


      ##############
      ## Effort ####
      ##############
      # effort dynamic
      tot_eff <- max_eff
      pmax_eff[1] <- sum(eff_space[1, ]) / tot_eff
      persis_pmax_eff[1] <- pmax_eff[1]

      # derby effort
      int_eff <- rs_season_file[1]
      slope_eff <- rs_season_file[2]
      cf_eff <- rs_season_file[3]
      # theorectical daily effort
      eff_season <- (exp(int_eff) + exp(slope_eff) * 1:365) * cf_eff
      eff_cap <- pmin(eff_season, tot_eff)
      eff_cap[(which(match(unique(eff_cap), tot_eff) == 1)+1):length(eff_cap)] <- 0


      ############################
      ## Model Initialization ####
      ############################
      recruit[1, ] <- nage[1, 1, ] <- exp(lro_FL) * ((pref_depth_std[round(depth), 1] * hab_suit_per) 
                                                      / sum(pref_depth_std[round(depth), 1] * hab_suit_per))
      nage[1, , ] <- exp(lro_FL) * ((pref_depth_std[round(depth), 1] * hab_suit_per) / sum(pref_depth_std[round(depth), 1] * hab_suit_per))

      # initial movement
      for (a in 2:amax) nage[1, a, ] <- nage[1, a, ] %*% pmove_init[, , a]
      vb[1, ] <- colSums(s_fa * nage[1, , ] * wa)
      ssb[1] <- sum(rowSums(nage[1, , ]) * agefec)


      # utility metrics
      cpue_ut[1, ] <- colSums(nage[1, , ] * wa * q)
      hpue_ut[1, ] <- colSums(nage[1, , ] * wa * q * ret)
      size_ut[1, ] <- colSums(nage[1, , ] * s_fa * la) / colSums(nage[1, , ] * s_fa)
      u_cpue[1, ] <- ut_cpue[1] / (1 + exp(ut_cpue[2] * (ut_cpue[3] - cpue_ut[1, ])))
      u_hpue[1, ] <- ut_hpue[1] / (1 + exp(ut_hpue[2] * (ut_hpue[3] - hpue_ut[1, ])))
      u_size[1, ] <- ut_size[1] / (1 + exp(ut_size[2] * (ut_size[3] - size_ut[1, ])))
      u_dist[1, ] <- ut_dist[1] / (1 + exp(ut_dist[2] * (ut_dist[3] - dist_wt)))
      u_crd[1, ] <- ut_crowd[1] / (1 + exp(ut_crowd[2] * (ut_crowd[3] - eff_space[1, ])))
      tot_ut[1, ] <- u_cpue[1, ] + u_hpue[1, ] + u_size[1, ] + u_dist[1, ] + u_crd[1, ]
      tot_ut_sum[1] <- sum(tot_ut[1, ])
      tot_init <- 1 * tot_ut_sum[1]


      #####################
      ## Time Dynamics ####
      #####################
      for (y in 2:nyears) {
        ## Effort dynamics
        for (k in 1:nsites) grav_wt[y, k] <- (tot_ut[y - 1, k] / cost[k])^w_cost
        pmax_eff[y] <- 1 / (1 + exp(-(tot_ut_sum[y - 1] - tot_init) / (sig1e * tot_init)))
        persis_pmax_eff[y] <- pmax_eff[y] * (1 - persis) + pmax_eff[y - 1] * persis # accounts for some "stickyness" in the effort time step to time step
        eff_space[y, ] <- 0
        # season length implementation
        if (y >= (fish_strt[1])) {
          for (k in 1:nsites) eff_init[k] <- (tot_eff * persis_pmax_eff[y] * grav_wt[y, k]) / sum(grav_wt[y, ])
          dens_eff <- dtruncnorm(x = eff_cap, a = 0, b = Inf, mean = sum(eff_init), sd = sum(eff_init) * 0.3)
          dens_eff <- dens_eff / sum(dens_eff) # needs to sum to 1 to prevent loss from density range
          exp_eff[y] <- sum(dens_eff * eff_cap) # expected "derby" effort
          pre_season[y] <- floor(uniroot(f = function(x) {
            ((exp(int_eff) + exp(slope_eff) * x) * cf_eff) - exp_eff[y]
          }, lower = 0, upper = 365, extendInt = "yes")$root)
          # effort
          for (k in 1:nsites) eff_space[y, k] <- exp_eff[y] * ((grav_wt[y,k]) / sum(grav_wt[y, ]))        
        }

        ## Fishing mortality
        FM[y, , ] <- q_t[y, ] * t(t(sel) * eff_space[y, ])
        FM[y, , ] <- FM[y, , ] + F_comm / nsites # add commerical fishery fishing mortality

        ## recruitment
        ssb[y] <- sum(rowSums(nage[y - 1, , ] * agefec))
        recruit[y, ] <- (4 * h * exp(lro_FL) * ssb[y]) / (ssbo * (1 - h) + ssb[y] * (5 * h - 1)) *
          (pref_depth_std[round(depth), 1] * hab_suit_per) /
          sum(pref_depth_std[round(depth), 1] * hab_suit_per) * exp(recdev[y])
        nage[y, 1, ] <- recruit[y, ]

        ## Mortality
        nage[y, 2:amax, ] <- nage[y - 1, 1:(amax - 1), ] * exp(-(M_site[y - 1, 1:(amax - 1), ] + FM[y - 1, 1:(amax - 1), ]))
        # plus group
        nage[y, amax, ] <- nage[y, amax, ] + nage[y - 1, amax, ] * exp(-(M_site[y - 1, amax, ] + FM[y - 1, amax, ]))

        ## Movement
        ind_vec <- ifelse(colSums(t(t(nage[y, , ]) * la^2)) > quantile(colSums(t(t(nage[1, , ]) * la^2)), 0.75), 1, 0)
        for (a in 1:amax) {
          pmove[, , a] <- t(t(exp(-lam * dist_mat)^w_dist) * (pref_depth_std[round(depth), a]^w_dep) *
                                                              (hab_suit_per^w_hab) * ((1 / colSums(t(t(nage[y, , ] * la^2)) 
                                                              / quantile(colSums(t(t(nage[y, , ]) * la^2)), 0.75))^0.5)^ind_vec))
          pmove[, , a] <- pmove[, , a] / rowSums(pmove[, , a])
        }
        for (a in 2:amax) nage[y, a, ] <- colSums(t(sapply(1:nsites, function(x) rmultinom(n = 1, size = round(nage[y, a, x]), prob = pmove[x, , a]))))

        ## Metrics
        vb[y, ] <- colSums(sel * nage[y, , ] * wa)
        catch_bio[y, ] <- colSums((FM[y, , ] / (FM[y, , ] + M_site[y, , ])) * nage[y, , ] * (1 - exp(-(M_site[y, , ] + FM[y, , ]))) * wa)
        catch_num[y, ] <- colSums((FM[y, , ] / (FM[y, , ] + M_site[y, , ])) * nage[y, , ] * (1 - exp(-(M_site[y, , ] + FM[y, , ]))))
        catch_numage[y, , ] <- (FM[y, , ] / (FM[y, , ] + M_site[y, , ])) * nage[y, , ] * (1 - exp(-(M_site[y, , ] + FM[y, , ])))

        ## Utility calcs
        cpue_ut[y, ] <- catch_num[y, ] / eff_space[y, ]
        cpue_ut[is.na(cpue_ut)] <- 0
        hpue_ut[y, ] <- colSums(catch_numage[y, , ] * ret) / eff_space[y, ]
        hpue_ut[is.na(hpue_ut)] <- 0
        size_ut[y, ] <- colSums(nage[y, , ] * la * s_fa) / colSums(nage[y, , ] * s_fa)
        size_ut[is.na(size_ut)] <- 0
        u_cpue[y, ] <- ut_cpue[1] / (1 + exp(ut_cpue[2] * (ut_cpue[3] - cpue_ut[y, ])))
        u_hpue[y, ] <- ut_hpue[1] / (1 + exp(ut_hpue[2] * (ut_hpue[3] - hpue_ut[y, ])))
        u_size[y, ] <- ut_size[1] / (1 + exp(ut_size[2] * (ut_size[3] - size_ut[y, ])))
        u_dist[y, ] <- ut_dist[1] / (1 + exp(ut_dist[2] * (ut_dist[3] - dist_wt)))
        u_crd[y, ] <- ut_crowd[1] / (1 + exp(ut_crowd[2] * (ut_crowd[3] - eff_space[y, ])))
        tot_ut[y, ] <- u_cpue[y, ] + u_hpue[y, ] + u_size[y, ] + u_dist[y, ] + u_crd[y, ]
        tot_ut_sum[y] <- sum(tot_ut[y, ])
      } # end of y loop, time dynamics


      #################################
      ## Data for Estimation model ####
      #################################
      ## fishing mortality
      FM <- FM[model_years, , ]

      ## catch data
      catch_obs <- rowSums(catch_bio[model_years, ]) * exp(catchdev[model_years])

      ## catch index
      cpue_index <- apply(catch_bio[model_years, ] / eff_space[model_years, ], 1, sum, na.rm = TRUE) * exp(indexdev[model_years])

      ## observed effort
      eff_obs <- rnorm(length(model_years), mean = rowSums(eff_space[model_years, ]), sd = cv_eff * rowSums(eff_space[model_years, ]))

      ## age comps
      pa <- array(NA, c(nyears, amax, nsites))
      for (y in model_years) {
        for (k in 1:nsites) {
          pa[y, , k] <- rmultinom(n = 1, size = ss_trip, prob = catch_numage[y, , k] / sum(catch_numage[y, , k]))
        }
      }
      pa <- apply(X = pa[model_years, , ], MARGIN = c(1, 2), FUN = sum)
      pa <- pa / rowSums(pa)

      ## age comps - FI
      pa_fi <- array(NA, c(nyears, amax, nsites))
      for (y in model_years) {
        for (k in 1:nsites) {
          pa_fi[y, , k] <- rmultinom(n = 1, size = ss_trip, prob = (s_fi * fi_q * nage[y, , k]) / sum(s_fi * fi_q * nage[y, , k]))
        }
      }
      pa_fi <- apply(X = pa_fi[model_years, , ], MARGIN = c(1, 2), FUN = sum)
      pa_fi <- pa_fi / rowSums(pa_fi)


      #################################
      ## Out of loop calculations ####
      #################################
      FM_save <- FM
      qt_save <- q_t[model_years, ]
      ssb_ssbo <- ssb[model_years] / ssbo
      ssb_save <- ssb[model_years]
      vb_save <- vb[model_years, ]
      nage_save <- nage[model_years, , ]
      recruit_save <- recruit[model_years, ]
      catch_bio_save <- catch_bio[model_years, ]
      catch_num_save <- catch_num[model_years, ]
      catch_numage_save <- catch_numage[model_years, , ]
      recdev_save <- recdev[model_years]
      M_site_save <- M_site[model_years, , ]
      eff_space_save <- eff_space[model_years, ]
      tot_ut_save <- tot_ut[model_years, ]
      pre_season_save <- pre_season[model_years]
      pmax_eff_save <- pmax_eff[model_years]
      cpue_ut_save <- cpue_ut[model_years, ]
      hpue_ut_save <- hpue_ut[model_years, ]
      size_ut_save <- size_ut[model_years, ]

      # for ARs
      ar_sites_list <- list()
      ar_num_list <- list()
      for (y in 1:nyears) {
        ar_sites_list[[y]] <- ar_sites
        ar_num_list[[y]] <- ar
      }


      ###################
      ## output list ####
      ###################
      output <- NULL
      output <- list(
        seed = seed,
        nyears = length(model_years),
        nsites = nsites,
        nlandings = nlandings,
        amax = amax,
        dist = dist_wt,
        dist_mat = dist_mat,
        pref_depth_std = pref_depth_std,
        depth = depth,
        coord_dat = coord_dat,
        lam = lam,
        ar_sites = ar_sites_list,
        ar_num = ar_num_list,
        nr_sites = nr_sites,
        coord_file = coord_file,
        ar_coord_file = ar_coord_file,
        rs_season_file = rs_season_file,
        ars_fishery_closure = NULL,
        log_ro = lro_FL,
        ssbo = ssbo,
        qt_save = qt_save,
        eff_obs_save = eff_obs,
        eff_space_save = eff_space_save,
        tot_eff_save = tot_eff,
        FM_save = FM_save,
        tot_ut_save = tot_ut_save,
        tot_init_save = tot_init,
        catch_num_save = catch_num_save,
        catch_bio_save = catch_bio_save,
        catch_numage_save = catch_numage_save,
        pre_season_save = pre_season_save,
        pmax_eff_save = pmax_eff_save,
        catch_obs_save = catch_obs,
        cpue_index_save = cpue_index,
        pa_save = pa,
        pa_fi_save = pa_fi,
        recdev_save = recdev_save,
        ssb_ssbo_save = ssb_ssbo,
        ssb_save = ssb_save,
        vb_save = vb_save,
        nage_save = nage_save,
        M_site_save = M_site_save,
        recruit_save = recruit_save,
        cpue_ut_save = cpue_ut_save,
        hpue_ut_save = hpue_ut_save,
        size_ut_save = size_ut_save
      ) # end of output list

      return(output)
    }) # end with function
  } # end of run_OM function
