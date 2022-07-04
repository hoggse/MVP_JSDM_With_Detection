#************************************************************************************
#
# Functions to run for Multivariate Probit JSDM for real data sets
#
# independent: hierarchial model, probit occurrent, logit detection
# dependent: P(Z = 1) = P(V>0), V ~ MVN(mu, Sigma)
#
# Functions
# 1. aug_y_cases.real -  processing of Observation data (from real data for collapsed and augmented forms)
# 2. set.windat.real  -  sets up data for a training or validation window on the data
# 3. set.win.cov.real -  sets up covariate data given real data
# 4. windats_id.list  -  gets list of identifying data for real data
# 5. Run_jags.clist.cov.rl  -  Runs  MVP Jags model analysis over a list of real data sets
# 6. saveJagsRes.cov.rl     -  saves MVP Jags model analysis results for real data
#
#'
#' Prepare Real Covariate and Species Data for MVP JSDM Jags model
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `set.win.cov.real()`,  `set.windat.real()`, `aug_y_cases.real()`,
#' `add_cov_data.real()`, provide tools for for reading in real species observations
#' and covariate data in csv format and putting into a format
#' ready to use with the JAGs MVP JSDM.
#'
#' `set.win.cov.real()` # setup real species observation data with occupancy and detection
#' covariates for JAGS/BUGS JSDM MGP model
#' `set.windat.real()` # creates window of species observation data for jags with training and
#' validation component for real data. removes z which is not available for real data
#' `aug_y_cases.real()` creates augmented species observation data  with training and
#' validation component for real data.  Substitutes first observation for random observation
#' `add_cov_data.real` adds covariate data to real observation data for single data type
#'
#' @section Common parameters used by real data functions
#' @param y       Raw species observation data (3 dimensional array nsite x nrep x nspec),
#'                either entire dataset or training dataset (if validation dataset present).
#' @param n.valid Length of the validation data (set to 0 or NULL if no training/validation)
#' @param nz      Non-zero species (extra species not observed, relec from kery community model)
#' @param y.valid Raw species observation data for validation data set
#' @param n.valid Length of the validation data (set to 0 if no training/validation)
#' @param nz      Not zero species (extra species not observed, relec from kery community model)
#' @param df      Number of degrees of freedom for wishart distribution, df=1 => flat, uninformative distribution
#' @param type    DAta type require set by integer values
#'     1 - prepares data for nsite X nrep X nspec species (or nspec+nz if nz !=0)
#'     2 - prepares collapsed data where y* = max(y[1,3]) (i.e. across replications)
#'     3 - random selection of replciation (but for real data, always choses first  set replications)
#'     0 - All 3 of the above options are prepared for fitting
#' @param dep    Dependent data flag, TRUE by default.  Sets df = nspec + df, sets id to square diag matrix of size M (M=nspec)
#'
#'




#**********************************************************************
# set.win.cov.real
#
# setup real observation data with occupancy/detection covariates for JAGS/BUGS
# using real data, so we have no mydata (simulated data output)
# myrealdat needs to include the following
# obs.tr, obs.va    Training and validation sets for y (observation) data
# hab.tr, hab.va    Habitation covariates training and validation data sets
# det.tr, det.va    Detection covariates training and validation data sets
# n.valid,n.train   lengths of the sets: just useful to have the set lengths here
#
# This function will set up data as follows:
#   Full covariate data: validation data appended to training data for covariates
#   Augmented y data: training data with NA appended for validation data
#   windats: data sets including both of these set up for jags for different modelling cases
#
#
#' @section `set.win.cov.real()` prepares species observation data depending on the requested data types.
#' @return A list of datasets for each requested datatypes plus an other entry.
#' For each data type, a data set in the form of a list which includes:
#' nsite, nspec, nrep, nz, df, M (augmented number of species nspec+nz)
#' The list also returns the y (observed data) and z (actual occupancy - simulated)
#' For each type the following item is returned
#' * Type 1: replicated survey data  y - 3D array
#' * Type 2: collapsed data y* - 2D array of summarised data where y*[i,k]=max(y[i,,k]) for species k at location i
#' * Type 3: random survey replication y* = [,rand_rep,]
#' * Type 4: use the actual observation data y* = z
#' The other data provides y, y* collapsed,
#'
#' @export set.win.cov.real

set.win.cov.real <- function(myrealdat, nz=0,df=1,type=0,dep=TRUE,preps=NULL, independent=TRUE) {

  # if requested alternative number of reps (e.g. 3 not 5, or two sets of reps 3 and 5)
  nreps = myrealdat$n.rep
  n.valid=myrealdat$n.valid
  multireps=NULL
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
    windat.sets <- set.windat.real(y=myrealdat$obs.tr[,1:nreps,],y.valid=myrealdat$obs.va[,1:nreps,],n.valid=n.valid,nz=nz,df=df,dep=dep,type=type) # get the windat and other data for training data
    windats[[1]] <- add_cov_data.real(windat.sets=windat.sets,covdat=myrealdat,type=type,reps=nreps,independent=independent)
  } else {
    for ( i in 1:length(multireps) ) {
      windat.sets <- set.windat.real(y=myrealdat$obs.tr[,1:multireps[i],],y.valid=myrealdat$obs.va[,1:nreps,],n.valid=n.valid,nz=nz,df=df,dep=dep,type=type) # get the windat and other data
      #str(windat.sets)
      windats[[i]] <- add_cov_data.real(windat.sets=windat.sets,covdat=myrealdat,type=type,reps=multireps[i])
      #str( windats[[i]])
    }
  } # end get all windat data

  return(windats)
}



#**********************************************************************
# set.windat.real
#
# creates window data for jags with training/validation component for real data
# removes z which is not available for real data
#
# y here is the training data
#' @section `set.windat.real` prepares real species observation data depending on the requested data types.
#' @return A list of datasets for each requested datatypes plus an other entry.
#' For each data type, a data set in the form of a list which includes:
#' nsite, nspec, nrep, nz, df, M (augmented number of species nspec+nz)
#' The list also returns the y (observed data) and z (actual occupancy - simulated)
#' For each type the following item is returned
#' * Type 1: replciated survey data  y - 3D array
#' * Type 2: collapsed data y' - 2D array of summarised data where y'[i,k]=max(y[i,,k]) for species k at location i
#' * Type 3: first survey replication y' = [,1,]
#' The other data provides y, y' collapsed, or y' first replication:
#' selecting items is not returned as y in the windat set.
#'
#'@export  set.windat.real
#' @example
#' require(gtools)
#' mydata = list()
#' nrep = 3
#' mydata$J = nrep
#' nsite = 100
#' nspec = 2
#' poccupy =0.5
#' pdetect=0.6
#' Z = matrix(data=rbinom(nsite*nspec,1,prob=poccupy),nrow=nsite)
#' Y = array(dim=c(nsite,mydata$J,nspec))
#' for ( i in 1:mydata$J) Y[,i,] = Z * matrix(data=rbinom(nsite*nspec,1,prob=pdetect),nrow=nsite)
#' myObsData = set.windat(y=Y,z=Z)
#'

set.windat.real <- function(y,y.valid=NULL, n.valid=0,nz=0,df=1,type=0,dep=TRUE) {
  #print(paste0("set windata dep",dep))
  caseData <- list()
  aug_y_dat <- aug_y_cases.real(y=y,n.valid=n.valid,nz=nz)
  if (is.null(y.valid)) {aug_y_valid=NULL} else {
    aug_y_valid <- aug_y_cases.real(y=y.valid,n.valid=0,nz=nz)  # get validation data into checking format
  }
  # case 0 = 1, 1 = standard y, 2 = clp, 3 =  first obs only)
  # note for our real data, only first obs will do for 3, as 2/3 obs may be missing
  if ( type ==0 || type ==1  ) {  # standard y with nrep reps
    caseData[[Windat.type[1] ]] <- list(y=aug_y_dat$y,nsite=aug_y_dat$nsite,
                                        #ntrain=aug_y_dat$ntrain,
                                        nspec=aug_y_dat$nspec,nrep=aug_y_dat$nrep,
                                        nz=nz,M=(aug_y_dat$nspec + nz) )
    if ( isTRUE(dep)) {
      #print("dependent add df and id")
      caseData[[Windat.type[1] ]]$df <- caseData[[Windat.type[1] ]]$M + df
      caseData[[Windat.type[1] ]]$id <- diag(caseData[[Windat.type[1] ]]$M)
    }
  }
  if ( type ==0 || type ==2  ) {  # collapsed y (max y over nrep reps)
    caseData[[Windat.type[2] ]] <- list(y=aug_y_dat$y_clp,nsite=aug_y_dat$nsite,
                                        #ntrain=aug_y_dat$ntrain,
                                        nspec=aug_y_dat$nspec,nrep=aug_y_dat$nrep,
                                        nz=nz,M=(aug_y_dat$nspec + nz) )
    if ( isTRUE(dep)) {
      #print("dependent add df and id")
      caseData[[Windat.type[2] ]]$df <- caseData[[Windat.type[2] ]]$M + df
      caseData[[Windat.type[2] ]]$id <- diag(caseData[[Windat.type[2] ]]$M)
    }
  }
  if ( type ==0 || type == 3  ) {  # random y (random 1 of nrep reps)
    caseData[[Windat.type[3] ]] <- list(y=aug_y_dat$y_rdm,nsite=aug_y_dat$nsite,
                                        #ntrain=aug_y_dat$ntrain,
                                        nspec=aug_y_dat$nspec,nrep=1,  # aug_y_dat$nrep,
                                        nz=nz,M=(aug_y_dat$nspec + nz) )
    if ( isTRUE(dep)) {
      #print("dependent add df and id")
      caseData[[Windat.type[3] ]]$df <- caseData[[Windat.type[3] ]]$M + df
      caseData[[Windat.type[3] ]]$id <- diag(caseData[[Windat.type[3] ]]$M)
    }
  }
  # no z data so have removed type==4, this type is not valid for real data

  caseData$other <- list(y_sum=aug_y_dat$y_sum,rrep=aug_y_dat$rrep,y_valid=aug_y_valid)
  if ( type == 1)       { caseData$other$y_clp= aug_y_dat$y_clp; caseData$other$y_rdm= aug_y_dat$y_rdm
  } else if (type == 2) {caseData$other$y= aug_y_dat$y; caseData$other$y_rdm= aug_y_dat$y_rdm
  } else if (type == 3) {caseData$other$y= aug_y_dat$y; caseData$other$y_clp= aug_y_dat$y_clp
  } # remove type 4 as no longer valid
  return(caseData)
}



#**********************************************************************
# aug_y_cases.real
#
# creates augmented data with training/validation component for real data
# substitutes first observation for random observation, as 1st is best for my data
#
#' @section `aug_y_cases.real` Processes 3D array of real species observation data.
#'                    Similar to aug_y_cases but always use first survey replication for rdm data.
#' @return  List of observation derived data:
#' *  y        Observed species data, optionally augmented with nz extra species
#' *  y_sum    Summed data, summarised by adding across all survey replciations for each species and site
#' *  y_clp    Collapsed data, summarised by taking maximum across all survey replciations for each species and site
#' *  y_rdm    Random survey replication data
#' *  nsite    The number of survey sites
#' *  nspec    The number of species in the survey (actually observed)
#' *  nrep     The number of survey replications
#' *  rrep     The random survey replication used for the y_rdm data
#'
#' @export aug_y_cases.real

aug_y_cases.real <- function(y,n.valid=0,nz=0) {
  # packaging Y values from getCorDataAllSigma or CreateOccAndObsAllSigs
  nsite = dim(y)[1]
  nspec = dim(y)[3]
  nrep = dim(y)[2]  # note this is the order used elsewhere!
  rrep = 1  # only the first observation is used, not random as for test data
  #rrep = sample.int(n=nrep,size=1,replace=FALSE)  # random repetition -
  # will no
  y_sum <- apply(y, c(1,3), sum, na.rm = T) # Collapse to detection frequency
  y_clp <- apply(y, c(1,3), max, na.rm = T) # Collapse to detected if detected anywhere, not otherwise
  if (nz > 0 ) {  # augments the number of species
    # augmented sum data (collapsed data to detection frequency)
    y_s.aug <- cbind(y_sum, array(0, dim=c(nsite, nz)))
    y_c.aug <- cbind(y_clp, array(0, dim=c(nsite, nz)))
    # augmented observation values (no collapsing of data)
    y_o.aug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
    y_o.aug[,,1:nspec] <- y # copy into it the observed data
    # Create same NA pattern in augmented species as in the observed species
    missings <- is.na(y_o.aug[,,1]) # e.g., third survey in high-elevation quads
    for(k in (nspec+1):(nspec+nz)){
      y_o.aug[,,k][missings] <- NA
    }
  } else {
    y_s.aug <- y_sum
    y_c.aug <- y_clp
    y_o.aug <- y
  }
  if ( n.valid > 0 ) {
    # no species in aug data
    n.aug.spec = nspec + nz
    y_s.aug <- rbind(y_s.aug,array(NA,dim=c(n.valid,n.aug.spec)))
    y_c.aug <- rbind(y_c.aug,array(NA,dim=c(n.valid,n.aug.spec)))
    y_o.aug1 <- array(NA, dim=c(nsite+n.valid, nrep, n.aug.spec))
    y_o.aug1[1:nsite,,1:n.aug.spec] <- y_o.aug
    y_o.aug <- y_o.aug1
  }
  y_r.aug <- y_o.aug[,rrep,]  # first observation y (previously random y)

  return(list(y=y_o.aug, y_sum=y_s.aug,y_clp=y_c.aug,y_rdm=y_r.aug,
              ntrain=nsite,nsite=nsite+n.valid,nspec=nspec,nrep=nrep,rrep=rrep) )
}


#**********************************************************************
# add_cov_data.real
#
# add covariates to the windat.sets data for real data
# need to amalgamate the training and validation set 1 after the other
# covdat = myrealdat
# nreps;type
# test.win.cov <- add_cov_data.real(windat.sets=windat.sets,covdat=myrealdat,type=type,reps=nreps)
# str(test.win.cov)
# covdat=myrealdat
#
#' @section `add_cov_data.real` Adds real covariate data for a set of one or more datatypes to related dataset for each type.
#'   Returns a windata set for each datatype requested.
#'
#' @export add_cov_data.real

add_cov_data.real <- function(windat.sets,covdat,type,reps, independent=TRUE) {
  # now add in the covariates
  if (is.null(windat.sets)) {print("set.win.cov::add_cov_data: warning windat set not generated"); return(NULL) }
  # regardless of what they do next, you need to set up the covariate data
  # habitat data: training data with validation data appended
  habitat <- rbind(covdat$hab.tr, covdat$hab.va)
  dim(habitat)
  # detection data: training data with validation data appended
  no_dec_vars_tr = ifelse ( (!is.null(covdat$det.tr) && !is.null(dim(covdat$det.tr))),ifelse( is.na(dim(covdat$det.tr)[3]),1,dim(covdat$det.tr)[3]),0 )
  reps_tr = ifelse ( (!is.null(covdat$det.tr) && !is.null(dim(covdat$det.tr))),ifelse( is.na(dim(covdat$det.tr)[2]),1,dim(covdat$det.tr)[2]),0 )
  # assume that this is the same for validation as well as training set
  detect <- array(NA, dim=c((covdat$n.train+covdat$n.valid),reps_tr,no_dec_vars_tr ) )
  detect[1:covdat$n.train,,] <- covdat$det.tr
  detect[(covdat$n.train+1):(covdat$n.train+covdat$n.valid),,] <- covdat$det.va
  s.reps <- apply(detect[,,1],1, function(x) {sum(!is.na(x))} )  # set site repetitions
  dim(detect)
  # note for our data: if a variable is missing for one repetition,
  #      the all detection covarites and observations are missing for that repetition
  #      observations are recorded in such a way as 1=1st made, 2 = 2nd made etc
  #              so if there is no 3rd repetition, 3 is NA, but 2nd is not impacted
  if (min(s.reps) == max(s.reps) ) s.reps=NULL;  # all obs there:not adding data

  if ( type!= 0 ) {  # single type set of value 1 to lenght(Windata.type)
    windat.sets[[Windat.type[type] ]] <- add_cov_data_1type(windat.type_set=windat.sets[[Windat.type[type] ]],
                                                            habitat=habitat,detect=detect,type=type,reps=reps,
                                                            independent=independent)
    if (type==1) { windat.sets[[Windat.type[type] ]]$s.reps =s.reps }  # set site reps for detection
    #str(windat.sets[[Windat.type[type] ]])
  } else {  # or if 0
    # type_cnt=1
    for (type_cnt in 1:(length(Windat.type)-1) ) {  #-1 as not including the occurrence case
      windat.sets[[Windat.type[type_cnt] ]] <- add_cov_data_1type(windat.type_set=windat.sets[[Windat.type[type_cnt] ]],
                                                                  habitat=habitat,detect=detect,type=type_cnt,reps=reps,
                                                                  independent=independent)
      if (type_cnt==1) { windat.sets[[Windat.type[type_cnt] ]]$s.reps =s.reps }  # set site reps for detection
    }# end type cnt loop
  }# end of type 0
  return(windat.sets)
}


