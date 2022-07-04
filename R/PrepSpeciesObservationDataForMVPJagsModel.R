###################################################################################
#
# Preparing species data for MVP JSDM
#
# This file contains functions that prep data from either real or simulated sources
# for analysis with JAGS

#************************************************************************************
#
# Functions
#
# from the augmented data, set up the rest of the windat for jags
# other data is put into another component of the return list incase needed later
# current types allowed are 0 to 3, where 0 runs all cases, 1 runs just standard, 2 runs just collapsed y, 3 = runs just random y
# global variable = required - note names can be changed if desired later.

#
#'
#' Prepare Species Observation Data  Functions for MVP JSDM
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `set.windat()`,  `aug_y_cases()`, `aug_z_cases`, and `Windat.type()` provide tools for
#' formating data from simulated mydata format that includes both species observation and
#' covariate data to be ready to use with the JAGs Multiviate Probit Joint Species
#' Distribution Model (MVP JSDM).
#'
#' `set.windat()` format species observation data for an "observation window" windat
#'     where y is the observed data over multiple observation and z is the actual occupancy data (simulations only)
#' `Windat.type`    Model data type required, set by string, with type values:
#'  *   "windat.std" - prepares data for nsite X nrep X nspec species (or nspec+nz if nz !=0)
#'  *   "windat.clp" - prepares collapsed data where y* = max(y[1,3]) (i.e. across replications)
#'  *   "windat.rdm" - random selection of replciation (but for real data, always choses first  set replications)
#' `aug_y_cases()` Processes observed occupancy data: required for real and simulation data
#'     Returns collapsed data, summed data, and allows for optional augmentation with unknown species
#' `aug_z_cases`   Processes latent occupancy data: required for simulation data
#'
#' @param y      Observed species data (3D array).
#' @param z      Actual species data (2D array).
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
#'
#'
#' @return Each function returns a list of data
#'     * windat
#'


#'
#' @section `set.windat` prepares species observation data depending on the requested data types.
#' @return A list of datasets for each requested datatypes plus an other entry.
#' For each data type, a data set in the form of a list which includes:
#' nsite, nspec, nrep, nz, df, M (augmented number of species nspec+nz)
#' The list also returns the y (observed data) and z (actual occupancy - simulated)
#' For each type the following item is returned
#' * Type 1: replciated survey data  y - 3D array
#' * Type 2: collapsed data y* - 2D array of summarised data where y*[i,k]=max(y[i,,k]) for species k at location i
#' * Type 3: random survey replication y* = [,rand_rep,]
#' * Type 4: use the actual observation data y* = z
#' The other data provides y, y* collapsed, y' random, or z (actual) data:
#' selecting items is not returned as y in the windat set.
#'


#'@export  set.windat
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

set.windat <- function(y,z,nz=0,df=1,type=0,dep=TRUE) {
  print(paste0("set windata dep",dep))
  caseData <- list()
  aug_y_dat <- aug_y_cases(y=y,nz=nz)
  aug_z_dat <- aug_z_cases(z=z,nz=nz)
  rrep=0
  Windat.types <- Windat.type()
  # case 0 = 1, 1 = standard y, 2 = clp, 3 = random)
  if ( type ==0 || type ==1  ) {  # standard y with nrep reps
    caseData[[Windat.types[1] ]] <- list(y=aug_y_dat$y,nsite=aug_y_dat$nsite,
                                         nspec=aug_y_dat$nspec,nrep=aug_y_dat$nrep,
                                         nz=nz,M=(aug_y_dat$nspec + nz) )
    if ( isTRUE(dep)) {
      #print("dependent add df and id")
      caseData[[Windat.types[1] ]]$df <- caseData[[Windat.types[1] ]]$M + df
      caseData[[Windat.types[1] ]]$id <- diag(caseData[[Windat.types[1] ]]$M)
    }
  }
  if ( type ==0 || type ==2  ) {  # collapsed y (max y over nrep reps)
    caseData[[Windat.types[2] ]] <- list(y=aug_y_dat$y_clp,nsite=aug_y_dat$nsite,
                                         nspec=aug_y_dat$nspec,nrep=aug_y_dat$nrep,
                                         nz=nz,M=(aug_y_dat$nspec + nz) )
    rrep= aug_y_dat$rrep
    if ( isTRUE(dep)) {
      #print("dependent add df and id")
      caseData[[Windat.types[2] ]]$df <- caseData[[Windat.types[2] ]]$M + df
      caseData[[Windat.types[2] ]]$id <- diag(caseData[[Windat.types[2] ]]$M)
    }
  }
  if ( type ==0 || type == 3  ) {  # random y (random 1 of nrep reps)
    caseData[[Windat.types[3] ]] <- list(y=aug_y_dat$y_rdm,nsite=aug_y_dat$nsite,
                                         nspec=aug_y_dat$nspec,nrep=aug_y_dat$nrep,
                                         nz=nz,M=(aug_y_dat$nspec + nz) )
    if ( isTRUE(dep)) {
      #print("dependent add df and id")
      caseData[[Windat.types[3] ]]$df <- caseData[[Windat.types[3] ]]$M + df
      caseData[[Windat.types[3] ]]$id <- diag(caseData[[Windat.types[3] ]]$M)
    }
  }
  if (type== 0 || type == 4 ) {  # put the z data into correct format for nodetection run
    caseData[[Windat.types[4] ]] <- list(y=aug_z_dat$y,nsite=aug_z_dat$nsite,
                                         nspec=aug_z_dat$nspec,nrep=aug_z_dat$nrep,
                                         nz=nz,M=(aug_z_dat$nspec + nz) )
    if ( isTRUE(dep)) {
      caseData[[Windat.types[4] ]]$df <- caseData[[Windat.types[4] ]]$M + df
      caseData[[Windat.types[4] ]]$id <- diag(caseData[[Windat.types[4] ]]$M)
    }
  }

  y_cors <- list(cor(apply(aug_y_dat$y,c(1,3),max)),cor(aug_y_dat$y_clp),cor(aug_y_dat$y_rdm) )
  y_tcors <- list(tetrachoric.try.rho(apply(aug_y_dat$y,c(1,3),max)),
                  tetrachoric.try.rho(aug_y_dat$y_clp),
                  tetrachoric.try.rho(aug_y_dat$y_rdm))
  caseData$other <- list(rrep=rrep,y_sum=aug_y_dat$y_sum,
                         y_cors=y_cors,y_tcors=y_tcors)
  if ( type == 1)       { caseData$other$y_clp= aug_y_dat$y_clp; caseData$other$y_rdm= aug_y_dat$y_rdm; caseData$other$y_occ = aug_z_dat$y
  } else if (type == 2) {caseData$other$y= aug_y_dat$y; caseData$other$y_rdm= aug_y_dat$y_rdm; caseData$other$y_occ = aug_z_dat$y
  } else if (type == 3) {caseData$other$y= aug_y_dat$y; caseData$other$y_clp= aug_y_dat$y_clp; caseData$other$y_occ = aug_z_dat$y
  } else if (type == 4) {caseData$other$y= aug_y_dat$y; caseData$other$y_clp= aug_y_dat$y_clp; caseData$other$y_rdm= aug_y_dat$y_rdm }
  return(caseData)
}


########################################
#
# some standard variables - now as a function
# Windat.type = c("windat.std","windat.clp","windat.rdm","windat.occ")
# Windat.type = c("windat.std","windat.clp","windat.rdm","windat.occ")
#
#' @section `Windat.type()`  function that returns the window data types
#'    which are used by 'set.windat' and other package functions
#'
#' @export

Windat.type <- function () c("windat.std","windat.clp","windat.rdm","windat.occ")






################################
#
# aug_y_cases
#
# Processing observed occupancy data: required for real and simulation data
# Create augmented and sum (collapsed data) in case required
#
#' @section `aug_y_cases()`  Processes species observation data for MVP JSDM.
#'       Returns 3D observed data, collapsed observation, and random survey replication data
#'       which may be augmented to allow for unknown species (following Kery).
#'       This data allows comparisions of replicated data, single survey data and collapsed data.
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
#'
#' @export
#'
aug_y_cases <- function(y,nz=0) {
  # packaging Y values from getCorDataAllSigma or CreateOccAndObsAllSigs
  nsite = dim(y)[1]
  nspec = dim(y)[3]
  nrep = dim(y)[2]  # note this is the order used elsewhere!
  rrep = sample.int(n=nrep,size=1,replace=FALSE)  # random repetition
  # will no
  y_sum <- apply(y, c(1,3), sum, na.rm = T) # Collapse to detection frequency
  y_clp <- apply(y, c(1,3), max, na.rm = T) # Collapse to detected if detected anywhere, not otherwise
  if (nz > 0 ) {
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
  y_r.aug <- y_o.aug[,rrep,]  # random y
  return(list(y=y_o.aug, y_sum=y_s.aug,y_clp=y_c.aug,y_rdm=y_r.aug,
              nsite=nsite,nspec=nspec,nrep=nrep,rrep=rrep) )
}

################################
#
# aug_z_cases
#
# Processing latent occupancy data: required for simulation data
#' @section `aug_z_cases()`  Processes actual species data for MVP JSDM (for simulations).
#'       Returns optionally augmented occupancy data (augmented to allow for extra unobserved species).
#'       This data allows comparisions of observed data with actual data for simulations.
#' @return  List of observation derived data:
#' *  y        Observed species data, optionally augmented with nz extra species
#' *  nsite    The number of survey sites
#' *  nspec    The number of species in the survey (actually observed)
#' *  nrep     The number of survey replications - set to 1
#'
#'
#' @export
aug_z_cases <- function(z,nz=0) {
  # packaging Z values from getCorDataAllSigma or CreateOccAndObsAllSigs
  # note no collapse so can get rid of other
  nsite = dim(z)[1]
  nspec = dim(z)[2]
  if (nz > 0 ) {
    # augmented occurence values values (no reps)
    z_o.aug <- cbind(z, array(0, dim=c(nsite, nz)))
  } else {
    z_o.aug <- z
  }
  return(list(y=z_o.aug, nsite=nsite,nspec=nspec,nrep=1) )  # note y set to z
  # id=list(mu,rho,sim,p,fid=rs)) ) - added on return
}

