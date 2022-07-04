#


# libraries required
# require(jagsUI)
# however jagsUI requires the installation of JAGS and other packages (e.g. rjags)
#




#*****************************************************************************
#
# Functions to run JAGS
#
#*****************************************************************************
#
#' Sets up standard parameter lists for jags and 
#' the standard model names for the jags models used by this package
#' u
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `get.params()`,  `get.jagsmodelfile()`
#' Helper functions to set up common parameter lists
#' and the jags model names used by this package
#' Multiple files and lists are used due to models
#' with and without detection and with varying numbers of covariates
#' 
#' #' @section Common parameters used by real data functions
#' @param param_type       String indicating parameter set type to be used.  
#'                         Parameter set types are related to the model type.
#'                         Current values inclde 
#'                         - "dep.det" model with detection
#'                         - "nocov.dep.det" model with detection but no covariates
#'                         - "dep.nodet" model without detection
#'                         - "nocov.dep" model without detection and no covariates
#'                         - "dep.det.ch" model with detection using cholesky appraoch
#'                         - "dep.nodet.ch" model without detection using cholesky approach
#'                         
#' @param model_type       String indicating model type to be used  
#'                         Current values inclde 
#'                         - "dep.det" model with detection
#'                         - "dep.det.Kcov" model with detection with K covariates, K > 1
#'                         - "dep.nodet" model without detection
#'                         - "dep.nodet.Kcov" model without detection with K covariates, K > 1
#'                         - "dep.det.ch" model with detection using cholesky appraoch
#'                         - "dep.nodet.ch" model without detection using cholesky approach
#'                         


# put commonly used model name and parameter lists here

####################################
# get.params
# gets standard input parameter values
# assumes that not requiring super population support
#
#' @section `get.params()` get jags parameter list
#' Returns: a list of parameters by name.  Must only contain parameters that exist in the model
#' #'
#' @export get.params
get.params <- function(param_type) {
  if ( param_type== "dep.det") {
    params <- c( "Tau",
                 "lpsi","lp","betalpsi","alphalp",
                 "mu.lpsi", "sd.lpsi",
                 "mu.betalpsi", "sd.betalpsi",
                 "mu.lp", "sd.lp",
                 "mu.alphalp", "sd.alphalp",
                 "mu.occ", "p","y","z","V" )
  } else  if ( param_type== "nocov.dep.det") {
    params <- c( "Tau",
                 "lpsi","lp",
                 "mu.lpsi", "sd.lpsi",
                 "mu.lp", "sd.lp",
                 "mu.occ", "p","y","z","V" )
  } else if ( param_type== "dep.nodet" ){
    params <-  c( "Tau",
                  "lpsi","betalpsi",
                  "mu.lpsi", "sd.lpsi",
                  "mu.betalpsi", "sd.betalpsi",
                  "mu.occ","y","V")
  } else if ( param_type== "nocov.dep" ){
    params <-  c( "Tau",
                  "lpsi",
                  "mu.lpsi", "sd.lpsi",
                  "mu.occ","y","V")
  } else if ( param_type== "dep.det.ch" ) {
    params <- c( "Tau","Sigma","Rho","A","B","Bj","Bk",
                 "lpsi","lp","betalpsi","alphalp",
                 "mu.lpsi", "sd.lpsi",
                 "mu.betalpsi", "sd.betalpsi",
                 "mu.lp", "sd.lp",
                 "mu.alphalp", "sd.alphalp",
                 "mu.occ", "p","y","z","V" )
  } else if ( param_type== "dep.nodet.ch" ) {
    params <-  c( "Tau","Sigma","Rho","A","B","Bj","Bk",
                  "lpsi","betalpsi",
                  "mu.lpsi", "sd.lpsi",
                  "mu.betalpsi", "sd.betalpsi",
                  "mu.occ","y","V")
  } else { params <- NULL}
  return( params )
}


# params.dep.det <- c( "Tau",
#                  "lpsi","lp","betalpsi","alphalp",
#                  "mu.lpsi", "sd.lpsi",
#                  "mu.betalpsi", "sd.betalpsi",
#                  "mu.lp", "sd.lp",
#                  "mu.alphalp", "sd.alphalp",
#                  "mu.occ", "p","y","z","V" )
# 
# params.dep.nodet <- c( "Tau",
#                  "lpsi","betalpsi",
#                  "mu.lpsi", "sd.lpsi",
#                  "mu.betalpsi", "sd.betalpsi",
#                  "mu.occ","y","V")
# 
# params.dep.det.ch <- c( "Tau","Sigma","Rho","A","B","Bj","Bk",
#                      "lpsi","lp","betalpsi","alphalp",
#                      "mu.lpsi", "sd.lpsi",
#                      "mu.betalpsi", "sd.betalpsi",
#                      "mu.lp", "sd.lp",
#                      "mu.alphalp", "sd.alphalp",
#                      "mu.occ", "p","y","z","V" )
# 
# params.dep.nodet.ch <- c( "Tau","Sigma","Rho","A","B","Bj","Bk",
#                        "lpsi","betalpsi",
#                        "mu.lpsi", "sd.lpsi",
#                        "mu.betalpsi", "sd.betalpsi",
#                        "mu.occ","y","V")



####################################
# get.jagsmodelfile
# gets name of the Jags model file
#
# note due odd jags reading these files only support one occupancy cov
# need to write second version for more than one
# an alternative is to provide the funciton that writes the file text
# but it will then need to be saved somewhere
#
#' @section `get.jagsmodelfile()` get jags model file name
#' Returns: a file name
#' #'
#' @export get.jagsmodelfile
get.jagsmodelfile <- function(model_type) {
  if ( model_type == "dep.det" ) {
    model_name  <-  "model11_10_sehII_dep_mvp_cov_v05b.txt"
  } else if ( model_type == "dep.det.Kcov" ) {
    model_name  <-  "model11_10_sehII_dep_mvp_Kcov_v05b.txt"
  } else if ( model_type == "dep.nodet" ) {
    model_name  <-  "model11_10_sehII_dep_nodet_mvp_cov_v05c.txt"
  } else if ( model_type == "dep.nodet.Kcov" ) {
    model_name  <-  "model11_10_sehII_dep_nodet_mvp_Kcov_v05c.txt"
  } else if ( model_type == "dep.det.ch" ) {
    model_name  <-  "model11_10_sehII_dep_mvp_cov_chv05b.txt"
  } else if ( model_type == "dep.nodet.ch" ) {
    model_name  <-  "model11_10_sehII_dep_nodet_mvp_cov_chv05c.txt"
  } else { model_name <- NULL }
  return( model_name )
}

# model_name.dep.det <- "model11_10_sehII_dep_mvp_cov_v05b.txt"
# model_name.dep.det.Kcov <- "model11_10_sehII_dep_mvp_Kcov_v05b.txt"
# model_name.dep.nodet <- "model11_10_sehII_dep_nodet_mvp_cov_v05c.txt"
# model_name.dep.nodet.Kcov <- "model11_10_sehII_dep_nodet_mvp_Kcov_v05c.txt"
# 
# model_name.dep.det.ch <- "model11_10_sehII_dep_mvp_cov_chv05b.txt"
# model_name.dep.nodet.ch <- "model11_10_sehII_dep_nodet_mvp_cov_chv05c.txt"

# set the following parameters and models
# with detection
# params = params.dep.det          # get.params("dep.det")
# model_name = model_name.dep.det  # get.jagsmodelfile("dep.det")
# no detection (used for collapsed data and single observation)  
# params.dep.nodet                 # get.params("dep.nodet")
# model_name.dep.nodet             # get.jagsmodelfile("dep.nodet")



# global variables
# model standard (without superpop varaibles)
# params.idep
# model_name.idep
# params.dep
# model_name.dep
#
# model without covariates (note: warnings in current windat data)
# params.nocov.idep
# model_name.nocov.idep
# params.nocov.dep
# model_name.nocov.dep
#
# model with superpop varaibles
# params.super.idep
# model_name.super.idep
# params.super.dep
# model_name.super.dep
#
#************************************************************************************
