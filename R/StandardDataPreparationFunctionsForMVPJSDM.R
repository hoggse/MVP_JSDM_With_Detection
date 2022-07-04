#************************************************************************************
#
# Formats general observation and covariate data for MVP JSDM
# assuming "mydata" simulation format
#
#
# Functions
# 1. set.win.dat  -  Set up data (with covariates) for dependent and independent models
# 2. add_cov_data -  add covariate occupancy and detection data for multiple model types and different numbers of replciations
# 3. add_cov_data_1type - Add covariate occupancy and detection data for single model type and a single give number of replications
#
#
#'
#' Standard (Simulation) Data Preparation Functions for MVP JSDM
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `set.win.cov()`,  `add_cov_data()`,
#' `add_cov_data_1type`, and `Windat.type` provide tools for formating data from simulated mydata format
#'  including both species observation and covariate data to be ready to use with the
#'  JAGs Multiviate Probit Joint Species Distribution Model (MVP JSDM).
#'
#' `set.win.cov()` set up species and coviarate data for required number of surey replications and
#'     and requested data types (replocated surveys (1), collapsed (2), single survey (random -3), or all (0))
#' `add_cov_data()`  adds formated covariate data to existing formated species data
#'      for multiple data types using same types as `set.win.cov()`
#' `add_cov_data_1type`  adds formated covariate data to existing formated species data
#'      for a single requested data type using same types as `set.win.cov()`
#'
#'
#' @param mydata  Base species observations and covariate data with fields.
#' *  J Number of survey replications where J=Nreps.
#' *  Y 3D array of observed data (Nsites by Nreps by Nspec).
#' *  Z 2D array of actual occupancy data (Nsites by Nspec) -
#'      valid for simulated data only.
#' *  habitat - List of occupancy covariates - assumed same for all survey replications.
#' *  detect - List of detection covariates - these may change with survey replication.
#' @param nz      Not zero species (extra species not observed, relec from kery community model).
#' @param df      Number of degrees of freedom for wishart distribution,
#'        df=1 => flat, uninformative distribution.
#' @param type    DAta type require set by integer values:
#'     1 - prepares data for nsite X nrep X nspec species (or nspec+nz if nz !=0);
#'     2 - prepares collapsed data where y* = max(y[1,3]) (i.e. across replications);
#'     3 - random selection of replciation (but for real data, always choses first  set replications);
#'     0 - All 3 of the above options are prepared for fitting.
#' @param dep    Osbserved species dependent data flag,
#'     TRUE by default.  Sets df = nspec + df, sets id to square diag matrix of size M (M=nspec)
#' @param preps  Array of the numbers of survey replications to use
#'     Can set multiple numbers, which must be less or equal to actual number of replications
#'     Default NULL - use the number of actual replications
#'     Used to investigate impact of replications on model performance
#' @param independent  Independent data flag for detection versus occupancy covariates
#'     TRUE by default, assumes occupancy and detection covariates are distinct for collapsed data
#' @param windat.sets  Formated species observations data for a single number of replciations, but mutliple data types
#' @param windat.type_set Formated species observations data for a single number of replciations and single data type
#' @param habitat list of occupancy covariates - assumed same for all survey replications
#' @param detect  list of detection covariates - these may change with survey replication
#' @param reps  Number of survey replications to be used for the data (esp detection covariates)
#'
#' @return
#'   * windats A list of formated datasets including observations and covariates
#'          for each requested windows datatype and each requested number of survey replciations
#'          first set of list - list by number of reps
#'          Inside rep list is
#'   * windat.sets A list of formated data sets
#'

#************************************************************************************


#************************************************************************************
#


########################################
#
# some standard variables
# Windat.type = c("windat.std","windat.clp","windat.rdm","windat.occ")
#
# Windat.type <- function () c("windat.std","windat.clp","windat.rdm","windat.occ")




#########################################################################
#
# set.win.cov
#
# set up species and coviarate data for required number of surey replications
# multiple survey replications can be requested using preps
#
# will use the aug_y_cases and aug_z_cases from BasicSimsJagsPrepV04.R
# will follow pattern of windep setups from BasicSimJagsPrep04
# which includes setting the matrices and degrees of freedom
# main difference is now have covariate data for occ and dec

# creating parameters to pass to JAGS/BUGS
# my data will now be the output from the AHM...Depv7c with some id extras.
# it uses windat.sets, aug_z_cases, aug_y_cases from BasicSimsJagsPrepV04
# also make call to add_cov_data (this file) for adding covariate approp
#
#'
#' @section `set.win.cov` creates list of datasets for requested datatypes and survey replication
#'     from data in (simulated) mydata format with covaraite and species observatin data.
#'
#' @export set.win.cov
#' @example
#' require(gtools)
#' mydata = list()
#' nrep = 3
#' mydata$J = nrep
#' nsite = 100
#' nspec = 2
#' poccupy =0.5
#' pdetect=0.6
#' ncov=2
#' mydata$habitat = matrix(data=rnorm(nsite*ncov),nrow=nsite)
#' mydata$detect = array(data=rnorm(nsite*nrep*ncov),dim=c(nsite,nrep,ncov))
#' mydata$Z = matrix(data=rbinom(nsite*nspec,1,prob=inv.logit(logit(poccupy)+rowSums(mydata$habitat))),nrow=nsite)
#' mydata$Y = array(dim=c(nsite,mydata$J,nspec))
#' for ( i in 1:mydata$J) mydata$Y[,i,] = mydata$Z * matrix(data=rbinom(nsite*nspec,1,prob=inv.logit(logit(pdetect)+rowSums(mydata$detect[,i,]))),nrow=nsite)
#' myFormatedData = set.win.cov(mydata)
#'

set.win.cov <- function(mydata, nz=0,df=1,type=0,dep=TRUE,preps=NULL, independent=TRUE) {

  # if requested alternative number of reps (e.g. 3 not 5, or two sets of reps 3 and 5)
  nreps = mydata$J
  multireps=NULL
  #nreps = ifelse ( (!is.null(preps) && preps>0 && preps<nreps ), preps, nreps)
  if ( !is.null(preps) ) {
    if ( length(preps) == 1 ) {
      nreps = ifelse ( (preps[1]>0 && preps[1]<nreps ), preps, nreps)
    } else {
      # now are they distinct from each other?
      if (max(preps)==min(preps)) {  # no then, just set as  before
        nreps = ifelse ( (preps[1]>0 && preps[1]<nreps ), preps, nreps)
      } else {
        multireps = rep(nreps,length(preps) )
        for (i in 1:length(preps)) { multireps[i] = ifelse(preps[i]>0 && preps[i]<nreps, preps[i],nreps) }
      }
    }# end multireps case
  }# end non null preps

  # get the base win data data by sending as follows:
  windats <- list()
  if (is.null(multireps)) {
    windat.sets <- set.windat(y=mydata$Y[,1:nreps,],z=mydata$Z,nz=nz,df=df,dep=dep,type=type) # get the windat and other data
    windats[[1]] <- add_cov_data(windat.sets=windat.sets,habitat=mydata$habitat,detect=mydata$detect,type=type,reps=nreps,independent=independent)
  } else {
    for ( i in 1:length(multireps) ) {
      windat.sets <- set.windat(y=mydata$Y[,1:multireps[i],],z=mydata$Z,nz=nz,df=df,dep=dep,type=type) # get the windat and other data
      #str(windat.sets)
      windats[[i]] <- add_cov_data(windat.sets=windat.sets,habitat=mydata$habitat,detect=mydata$detect,type=type,reps=multireps[i])
      #str( windats[[i]])
    }
  } # end get all windat data
  return(windats)
}


#################################
#
# add_cov_data
#
# sets up the occupancy and detection covariate data for a given data set
# can set up the data for a number of different model types including:
#      multiple surveys, collapsed, single replication
# type0 means a data set is setup for each model type
#
# Add covariate data for a single
# v07_02: updated so that variables are set correctly for the no detection models
#  for type 1: detection model with detection parameters
#  for type 2: take mean of detection parameters and add to habitat parameters (assumes different)
#  for type 3: takes matching detection for the randomly selected y and adds to habitat parameters
#  for type 4: occupancy: no need to allow for detection, so habitat vars only.
#'
#' @section `add_cov_data` Adds covariate data for a set of one or more datatypes to related dataset for each type.
#'   Returns a windata set for each datatype requested.
#'
#' @export
add_cov_data <- function(windat.sets,habitat,detect,type,reps, independent=TRUE) {
  # now add in the covariates
  Windat.types <- Windat.type()
  if (is.null(windat.sets)) {print("set.win.cov::add_cov_data: warning windat set not generated"); return(NULL) }
  if ( type!= 0 ) {  # single type set of value 1 to lenght(Windata.type)
    windat.sets[[Windat.types[type] ]] <- add_cov_data_1type(windat.type_set=windat.sets[[Windat.types[type] ]],
                                                            habitat=habitat,detect=detect,type=type,reps=reps,
                                                            independent=independent)

  } else {  # or if 0
    for (type_cnt in 1:length(Windat.types)) {
      windat.sets[[Windat.types[type_cnt] ]] <- add_cov_data_1type(windat.type_set=windat.sets[[Windat.types[type_cnt] ]],
                                                                  habitat=habitat,detect=detect,type=type_cnt,reps=reps,
                                                                  independent=independent)
    }# end type cnt loop
  }# end of type 0
  return(windat.sets)
}


#############################################
#
# add_cov_data_1type
#
# add covariate data (detection and occupancy) for a single model type

#'
#' @section `add_cov_data_1type` Adds covariate data for one datatype and its related datasets.
#'     Returns a windatset specified by a list with the following entries:
#'     * hab         Occupancy covariates
#'     * wind        Detection covariates (for type 1 and type 3 only)
#'     * novar.occ   Number of occupancy covariates
#'     * novar.det   Number of detection covariates (for type 1 and type 3 only)
#'     * Data set up previously by `Set.windat()` which includes
#'     *  * y    Observed occupancy data
#'     *  * z    Actual occupancy data (simulated data only)
#'     *  * nsite  Number of sites
#'     *  * nspec  Number of species
#'     *  * nrep   Number of survey replications
#'
#'
#' @export
add_cov_data_1type <- function(windat.type_set,habitat,detect,type,reps, independent=TRUE) {
  # assume single type: ie type = 1,2,3 or 4
  if (type == 1 ) {
    windat.type_set$hab= habitat
    windat.type_set$wind = detect[,1:reps,]
    windat.type_set$novar.occ = dim(habitat)[2]
    windat.type_set$novar.det = dim(detect)[3]
  } else if (isTRUE(independent)) {
    if ( type == 2 ) {
      windat.type_set$hab= cbind(habitat,apply( detect[,1:reps,],c(1,3),mean,na.rm=TRUE ))
      windat.type_set$novar.occ = dim(habitat)[2] + dim(detect)[3]
    } else if( type == 3 ) {
      rrep = ifelse( !is.null(windat.sets$other$rrep) && windat.sets$other$rrep >0 ,windat.sets$other$rrep ,1)
      rrep = ifelse( rrep > dim(detect)[2],dim(detect)[2],rrep)
      windat.type_set$hab = cbind(habitat,detect[,rrep,])
      windat.type_set$novar.occ = dim(habitat)[2] + dim(detect)[3]
    } else  { # type 4: occupancy - so no need to include detection vars at all
      windat.type_set$hab= habitat
      windat.type_set$novar.occ = dim(habitat)[2]
    } # end type 4 (occupancy) t
  } else {  # not dependent, so only use habitat vars
    windat.type_set$hab= habitat
    windat.type_set$novar.occ = dim(habitat)[2]
  }
  #str(windat.type_set)
  return(windat.type_set)
}


