#************************************************************************************
#
# Standardise covariate data and generate sequential test ranges
#
#
# Functions
# 1. std_occ_vars - standarised covariates for occurence (habitat)
# 2. std_det_vars - standarised covariates for detection (detect)
# 3. std_covars   - standarised covariates for any 2D array, matrix or dataframe
#
#************************************************************************************
#
# General Data prep
#
# 1. covariates for occurence (habitat) and detection (detect)
# 2. Observation data (including actual, collapsed and augmented forms)
#
#************************************************************************************
#'
#' General Data Utilities
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `std_occ_vars()`,  `std_det_vars()`, and `std_covars` are tools to standardise covariate
#' data, it returns the standardised data, mean, standard devation, simulated sequence covariate test data
#'
#' `std_occ_vars`    standardises occupancy covariates - assumes covaraties in 2D array
#' `std_det_vars`    standardises detection covariates - assumes covaraties in 3D array
#' `std_covars`      standardises covariates - generalisation of std_occ_vars - assumes covaraties in 2D array
#'
#' @section this is a group of functions that provide similar capabilities for standarising
#' covariate data.  `std_covars` and `std_occ_vars` provide the same functionality,
#' `std_occ_vars` uses the older "habitat" name for occupancy covariates.
#'
#' @param  habitat Occupancy covariates in 2D array format (Nsite by Ncov),
#'               *  Nsite is the number of sites,
#'               *  Ncov is the nuumber of covariates
#' @param  detect  Detection covariates in 3D array format - (Nsite by Nrep by Ncov) assume change over time
#'               *  Nrep is the number of survey replications (at different time periods)
#' @param  covars  General covariates in 2D array format  (Nsite by Ncov).
#'
#' @return list  Each function returns a list with the following components.
#'      * .range  Sequence of test data over covariate range (hab.range,det.range,cov.range).
#'      * .std    Standardised data for all covariates (hab.std,det.std,cov.std).
#'      * sd     Vector of standard deviations for each covariate.
#'      * mean   Vector of means for each covariate.


##################################
#
# std_occ_vars <- function(mydata) {
#
# data: uses occupancy variables as: habitat <- mydata$detect
# assumes that data is 2D matrix, array, or dataframe
#' @section `std_occ_vars` standardise occupancy data (habitat) in 2D array
#'   Returns list including:
#'   - hab.range simulated range data det.range,
#'   - hab.std   standardised data,
#'   - sd        standard deviation of given data,
#'   - mean      mean of given data.
#'
#' @export std_occ_vars
std_occ_vars <- function(habitat) {
  datarange <- seq(-2,2,0.1)
  HAB <- array(data=NA, dim(habitat))
  HAB.datarange  <- array(data=NA, c (length(datarange),dim(habitat)[2] ) )
  i=1
  for (i in 1:dim(habitat)[2] ) {
    hab.max = max(habitat[,i],na.rm=TRUE)
    hab.min = max(habitat[,i],na.rm=TRUE)
    hab.mean = mean(habitat[,i],na.rm=TRUE)
    hab.sd  =  sd(habitat[,i],na.rm=TRUE)
    HAB[,i] <- (habitat[,i] - hab.mean) / hab.sd
    HAB.datarange[,i] <- datarange
    #HAB.datarange[,i] <- (datarange- hab.mean) / hab.sd
  }
  return(list(hab.range = HAB.datarange, hab.std=HAB, sd=hab.sd, mean=hab.mean))
}


##################################
#
# std_det_vars <- function(mydata) {
#
# data: uses detection variables as: detect <- mydata$detect
# assumes that data is time varying
#' @section `std_det_vars` standardise detection data in 3D array
#'   Returns list including:
#'   - det.range simulated range data det.range,
#'   - det.std   standardised data,
#'   - sd        standard deviation of given data,
#'   - mean      mean of given data.
#'
#' @export std_det_vars

std_det_vars <- function(detect) {
  datarange <- seq(-2,2,0.1)
  DET <- array(data=NA, dim(detect))
  DET.datarange  <- array(data=NA, c (length(datarange),dim(detect)[2],dim(detect)[3] ) )
  dim(DET.datarange) ; dim(DET) ; dim(detect)
  #i=1; j=1
  for (i in 1:dim(detect)[3] ) { # no covariates
    det.max = max(detect[,,i],na.rm=TRUE)
    det.min = max(detect[,,i],na.rm=TRUE)
    det.mean = mean(detect[,,i],na.rm=TRUE)
    det.sd  =  sd(detect[,,i],na.rm=TRUE)
    DET[,,i] <- (detect[,,i] - det.mean) / det.sd
    for( j in 1:dim(detect)[2] ) {  # no repeats
      if ( j == 1 ) {
        #DET.datarange[,j,i] <- (datarange - det.mean) / det.sd
        DET.datarange[,j,i] <- datarange
      } else {
        DET.datarange[,j,i] <- sample(DET.datarange[,1,i],size=length(datarange), replace = TRUE)
      }
    }
  }
  return(list(det.range = DET.datarange, det.std=DET, sd=det.sd, mean=det.mean))
}# end function


#################################
# std_covars
#
# general standardisation, assumes covariates in data frame, matrix or 2d array
# if not use standard covars instead (same has hab.)
#' @section `std_covars` standardise covariate data in 2D array
#'   Returns list including:
#'   - cov.range simulated range data det.range,
#'   - cov.std   standardised data,
#'   - sd        standard deviation of given data,
#'   - mean      mean of given data.
#'
#' @export std_covars
std_covars <- function(covars) {
  datarange <- seq(-2,2,0.1)
  COV <- array(data=NA, dim(covars))
  COV.datarange  <- array(data=NA, c (length(datarange),dim(covars)[2] ) )
  i=1
  for (i in 1:dim(covars)[2] ) {
    cov.max = max(covars[,i],na.rm=TRUE)
    cov.min = max(covars[,i],na.rm=TRUE)
    cov.mean = mean(covars[,i],na.rm=TRUE)
    cov.sd  =  sd(covars[,i],na.rm=TRUE)
    COV[,i] <- (covars[,i] - cov.mean) / cov.sd
    COV.datarange[,i] <- datarange
  }
  return(list(cov.range = COV.datarange, cov.std=COV, sd=cov.sd, mean=cov.mean))
}


