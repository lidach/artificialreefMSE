#' @title run_all_fn
#'
#' @description function to open all data files and code for MSE
#' doesn't run like function (sourced)
#'
#' @param datadir directory where all data are stored
#' @param Rdir directory where all the R code is stored
#' @param TMBdir directory where TMB code is stored

########################################
## Required packages for all models ####
########################################
require(sp)
require(rgdal)
require(raster)
require(rworldmap)
require(maps)
require(mapdata)
require(Hmisc)
require(TMB)
require(hutilscpp)
require(sf)
require(dplyr)
require(MASS)
require(truncnorm)
require(parallel)



#################
## Data sets ####
#################
# will be sourced externally
# spatial extent, HB sites
load(file.path(datadir, "GulfBay_dat.RData"))
# AR sites
load(file.path(datadir, "MSE_AR_sites.RData"))
# landing sites
load(file.path(datadir, "landingsites_dat.RData"))
# red snapper spatial catch at age data (bottom longline and vertical line observer database)
load(file.path(datadir, "RS_ssdf.Rdata"))
# red snapper tagging data
RS_tag_file <- read.csv(file.path(datadir, "Tag_snapper.csv"))
# red snapper effort (for predicting future effort)
load(file.path(datadir, "pred_eff_NW_FL.RData"))



# ###########################
# ## Run Operating model ####
# ###########################
# create parameters for operating models
source(file.path(Rdir, "create_OM_input.R"))
# additional functions
source(file.path(Rdir, "extra_fn.R"))
# # run initial operating model (set number of artificial reefs and no management strategy)
source(file.path(Rdir, "run_init_OM.R"))
# run operating model into future (create artificial reefs and use management strategies)
source(file.path(Rdir, "run_future_OM.R"))



# ############################
# ## Run Estimation model ####
# ############################
# # create parameters for estimation model
source(file.path(Rdir, "create_EM_data.R"))
# run age-structured assessment model (estimation model)
source(file.path(Rdir, "run_EM.R"))
# compile ASAM
compile(file.path(TMBdir, "ASAM.cpp"))
dyn.load(dynlib(paste0(TMBdir,"/ASAM")))



# #############################
# ## Run Management models ####
# #############################
# calculate reference point (for TAC)
source(file.path(Rdir, "calc_ref.R"))
# calculate total allowable catch (TAC) - for reference strategies
source(file.path(Rdir, "TAC_calc.R"))
# rules to place artificial reefs
source(file.path(Rdir, "AR_rule.R"))
# run overall MSE model
source(file.path(Rdir, "run_MSE.R"))

