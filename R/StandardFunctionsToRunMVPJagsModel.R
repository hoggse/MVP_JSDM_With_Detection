#************************************************************************************
#
#
# Functions to run dependent and independent JSDM based on multivariate probit model
#
# independent: hierarchial model, probit occurrence, logit detection
# dependent (Multivariate probit model): P(Z = 1) = P(V>0), V ~ MVN(mu, Sigma)
#
# Functions
# 1. runJags.case.cov - Run the MVP model: set up initial values, params, call jags
# 2. saveJagsRes.cov  - Save the result file from Jags analysis
# 3. Run_jags.clist.cov - Run the MVP model over a list of files
#
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `runJags.case.cov()`,  `saveJagsRes.cov()`, and `Run_jags.clist.cov` provide tools for
#'  for running the JAGs Multiviate Probit Joint Species Distribution Model (MVP JSDM)
#'  assuming data in a species observation and covariate "windat" format
#'
#' `runJags.case.cov`     Runs JAGS analysis on a single combination of covariates for given set of species data
#' `saveJagsRes.cov'      Saves the results of the JAGS analysis
#' `Run_jags.clist.cov`   Runs a list of cases for a given set of data - may set different covariates
#'                        different modelling approaches (e.g. with/without detection o
#'                        or do/don't usecholesky appraoch)
#'
#' @section Common parameters used by real data functions
#' @param obscase      Data for the covariate case including species observations and covariate data
#' @param params       List of parameters for the JAGS analysis to return.
#'                     Use get.params() function for standard sets.
#' @param model_name   Name of the Jags model file
#'                     Use get.jagsmodelfile() function for standard models.
#'                     This parameter sets type of modelling approach (e.g with/without detection)
#' @param choleskyFlg  Set true to use Cholesky approach (needs specific priors and model)
#' @param h            Parameter for Cholesky estiamtes
#' @param g            Parameter for Cholesky estiamtes
#' @param ni           Number of iterations for the jags pas
#' @param nt           Thining rate for jags
#' @param nb           Number of burn in iterations for jags
#' @param nc           Number of chains used by jags - usually 3 or 4
#' @param debug        Use JAGS debug



#****************************************************************************
# runJags.case.cov
#
# function to run jags on a single case from list of observation cases of 2 species
# assumes covariates for occurrence and detection
# passed with windat data ready to go
# creates the initialisation variables, depending on nz and detection or no detection.
# returns the jags output
# type determines type of model to use: 1 -> standard y wiht n reps, 2) collapsed y, 3) random singl rep of y
#              1 is used for models with detection, 2 and 3 are for models with no detection
# run the following setsps
#   step - 2.1 extract windat data (note that .df and .ID created with windata data now)
#   step - 2.2 create Tau and sigma (always dependendent model for mo), zst and wst
#   step-  2.3 create inits.dep
#   step-  2.4 run jags
#   step-  2.5 any preliminary checks - none so far
#
# note wnd.idx and type work exactly the same as for in BasicSimJagsPrepV04

# set the following parameters and models no cholesky
# detection
# params = params.dep.det
# model_name = model_name.dep.det
# no detection
# params = params.dep.nodet
# model_name = model_name.dep.nodet

# set the following parameters and models with cholesky
# detection
# params = params.dep.det.ch
# model_name = model_name.dep.det.ch
# no detection
# params = params.dep.nodet.ch
# model_name = model_name.dep.nodet.ch
#
#' @section `runJags.case.cov()` prepares species observation data depending on the requested data types.
#' @return Returns jagsui::jags() call output.
#' A list that includes the MCMC chainout puts for each parameter
#' and other data relating to the JAGS run
#'
#' @export runJags.case.cov


runJags.case.cov <- function(obscase, type=1, wnd.idx=1,
                         params= get.params("dep.det"),  #params.nocov.dep,
                         model_name=get.jagsmodelfile("dep.det"), # model_name.nocov.dep,
                         choleskyFlg=FALSE, h=NULL, g=NULL,
                         ni=1500, nt=2, nb=500, nc=3,na=500, debug=FALSE ) {
  myId = obscase$id
  # step 2.1 and 2.2 (already done on extraction steps)
  # so get the right windat.dep data (which will have df and id set already)
  windat.dep = obscase$windat[[wnd.idx]][[Windat.type[type] ]]
  if (type==1 && length(dim(windat.dep$wind) ) == 2 ) {
    # needs to be free.
    tmp <- array(dim=c(dim(windat.dep$wind),1) )
    tmp[,,1] = windat.dep$wind
    windat.dep$wind=tmp
  }
  # set up for cholesky
  if ( isTRUE(choleskyFlg) && !is.null(h) && ! is.null(g) ) {
    print(paste0("runJags.case.cov: ","Setting up for tuneable cholesky"))
    windat.dep$h = h; windat.dep$g=g
  } else { if(isTRUE(choleskyFlg)) { print(paste0("runJags.case.cov: ","Could not setup tuneable cholesky: need correct g and h"))} }
  Tau.dep = rwish(windat.dep$df, windat.dep$id); # step 2.3
  Sigma.dep = solve(Tau.dep)
  wst.dep <- rep(1, windat.dep$nspec+windat.dep$nz) # Simply set everybody as occurring
  zst.dep <- array(1, dim = c(windat.dep$nsite, windat.dep$M)) # ditto

  # pollock suggested initial variable values
  .X = windat.dep$hab # occ covariates
  .X <- cbind(rep(1, nrow(windat.dep$y)), .X)
  dim(.X)
  V <- abs(t(replicate(windat.dep$nsite, mvrnorm(1, rep(0, windat.dep$M), Sigma.dep)))) # males all - sima used here
  if( type != 1 ) {
    V <- ifelse(as.matrix(windat.dep$y), V, -1 * V)  # setting postive or negative correlation
  } else {
    Z = apply(windat.dep$y,c(1,3),max,na.rm=TRUE)
     V <- ifelse(as.matrix(Z), V, -1 * V)  # setting postive or negative correlation
  }
  if( !is.null(obscase$windat[[wnd.idx]]$other$y_valid) ) {
    nvalid <- obscase$windat[[wnd.idx]]$other$y_valid$nsite
  } else {nvalid <- 0}
  if ( nvalid > 0 ) {
    for ( i in (windat.dep$nsite-nvalid+1):windat.dep$nsite) {
      if( type == 1 ) { Z[i,] = NA;}
      V[i,]=NA
    }
  }
  Sigma <- mvnXXX(V[1:(windat.dep$nsite-nvalid),])$parameters$variance$sigma[, , 1]
  if( type != 1 ) {
    Beta <- t(sapply(seq_len(ncol(windat.dep$y)),
                     function(x) {unname(coef(glm(windat.dep$y[, x] ~ 1,
                                                  family=stats::binomial(link=probit))))}))
    lpsi <- Beta * sqrt(diag(Sigma))
    lpsi <- as.vector(lpsi)
    mu.lpsi = mean(Beta)
    sd.lpsi = ifelse(sd(Beta) < 0.1 , 0.1, sd(Beta))
  } else {
    # no covariates - just intercept
    Beta <- t(sapply(seq_len(ncol(Z)),
                     function(x) {unname(coef(glm(Z[, x] ~ 1,
                                                  family=stats::binomial(link=probit))))}))
    lpsi <- Beta * sqrt(diag(Sigma))
    lpsi <- as.vector(lpsi)
    mu.lpsi = mean(Beta)
    sd.lpsi = ifelse(sd(Beta) <0.1 ,0.1, sd(Beta))
    lp= apply((apply(windat.dep$y,c(1,3),sum,na.rm=TRUE)/dim(windat.dep$y)[2]),2,sum,na.rm=TRUE)/apply(Z,2,sum,na.rm=TRUE)
    mu.lp = mean(lp)
    sd.lp = ifelse(sd(lp) < 0.1 , 0.1, sd(lp))
  }
  Tau=solve(Sigma)
  # create inits - step 2.4
  if ( windat.dep$nz > 0 ) {  # assume you need super population data if nz>0
    if (type == 1 ) {
      if ( isTRUE(choleskyFlg)){ # nz>0, type=1, cholesky
        inits.dep <- function() { list(
          w=wst.dep,
          Z = Z,
          V = V,
          lpsi = lpsi, mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
          lp = lp, mu.lp = mu.lp, sd.lp = sd.lp,
          betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                            nrow=windat.dep$M),
          alphalp = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.det),
                           nrow=windat.dep$M)
        ) }
      } else {# nz>0, type=1, not cholesky
        inits.dep <- function() { list(
          w=wst.dep,
          Z = Z,
          Tau = Tau,
          V = V,
          lpsi = lpsi, mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
          lp = lp, mu.lp = mu.lp, sd.lp = sd.lp,
          betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                            nrow=windat.dep$M),
          alphalp = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.det),
                           nrow=windat.dep$M)
        ) }
        }# end nz>0, type=1, not choleksy
    } else {  # type!=1: no y, so y collapses to z
      if ( isTRUE(choleskyFlg)){ # nz>0, type!=1, cholesky
        inits.dep <- function() { list(
          w=wst.dep,
          V = V,
          Tau = Tau,
          lpsi = lpsi,
          mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
          betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                            nrow=windat.dep$M)
        ) }
      } else { # nz>0, type!=1, not cholesky
        inits.dep <- function() { list(
          w=wst.dep,
          V = V,
          lpsi = lpsi,
          mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
          betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                            nrow=windat.dep$M)
        ) }
      }# end nz>0, type!=1, not cholesky
    } # end type !=1, nz >0
  } else {  # but no super population otherwise:  nz=0
    if (type == 1 ) {
      if ( isTRUE(choleskyFlg)){ # nz=0, type=1, cholesky
        inits.dep <- function() { list(
          Z = Z,
          V = V,
          lpsi = lpsi, mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
          lp = lp, mu.lp = mu.lp, sd.lp = sd.lp,
          betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                            nrow=windat.dep$M),
          alphalp = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.det),
                           nrow=windat.dep$M)
        ) }
      } else { # nz=0, type=1, not cholesky
        inits.dep <- function() { list(
          Z = Z,
          Tau = Tau,
          V = V,
          lpsi = lpsi, mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
          lp = lp, mu.lp = mu.lp, sd.lp = sd.lp,
          betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                            nrow=windat.dep$M),
          alphalp = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.det),
                           nrow=windat.dep$M)
        ) }
      } # end nz=0, type=1, not cholesky
      } else { #  start nz=0, type!=1  (end type==1)
        if ( isTRUE(choleskyFlg)){ # nz=0, type!=1, cholesky
          inits.dep <- function() { list(#z = zst.dep,
            V = V,
            lpsi = lpsi,
            mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
            betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                              nrow=windat.dep$M)
          ) }
        } else { # nz=0, type!=1, not cholesky
          inits.dep <- function() { list(#z = zst.dep,
            Tau = Tau,
            V = V,
            lpsi = lpsi,
            mu.lpsi = mu.lpsi, sd.lpsi = sd.lpsi,
            betalpsi = matrix(data=rnorm(n = windat.dep$M*windat.dep$novar.occ),
                              nrow=windat.dep$M)
          ) }
        }# end nz=0, type!=1, not cholesky
    }# end type!=1, nz=0
  }# end nz=0

  # now run JAGS model
  #windat.dep$wind[1,1,3] - errors with NAs, need to work out how to deal
  occ_dat.jags <- jagsUI::jags(windat.dep, inits.dep, params, model_name,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.adapt=na,
                       n.burnin = nb, parallel = !debug)
  return(occ_dat.jags )
}


############################################
#
# saveJagsRes.cov
# saving the results
#
#' @section `saveJagsRes.cov()` Saves analysis output from JAGs run to file, assuming part of simulation run.
#' @return Lists including.
#' - resFile - Name of file where results are stored
#' - resID   - Id data for the results file (based on the simulated data info)
#' and other data relating to the JAGS run
#'
#' @export runJags.case.cov

saveJagsRes.cov <- function(jagsRes,myId,mu_case,p_case,nreps,type=1,fid="",choleskyFlg=FALSE,savedir="") {
  #print(paste0("my id ",myId))
  ch_str = ifelse(isTRUE(choleskyFlg),"_ch_","")
  Filename = paste0("JagsRes_","CovCs",myId$cov_case,"Rho",myId$rho,"MuCs",mu_case,
                    "DetpCs",p_case,"Sim",myId$sim_no,"_",myId$dsim_no,"_tp",type,
                    ch_str,"rep",nreps,fid,".RDS",sep="")
  Fullfilename = paste0(savedir,Filename,sep="")
  print(paste0("Full filename ",Fullfilename,sep=""))

  # may want to extract summary of results to return, e.g.
  # tau, Mu (lpsi), p (lp), get median, mean, 95% and 5% results (look at others)
  myId$my_case = mu_case
  myId$p_case = p_case
  myId$nreps = nreps
  myId$type =1
  myId$fid = fid

  saveRDS(list(jagsRes=jagsRes,resId=myId), file=Fullfilename )
  #Sys.sleep(1)
  return(list(resFile=Filename,resID=myId))
  #print(paste0("Case ",case,"saved as ",Filename,sep=" "))
}





##########################################
#
# Run_jags.clist.cov
#
# running the MVP jags model over multiple datasets in a list

#' @section `Run_jags.clist.cov()`
#' Runs JAGS analysis on a list of different data sets generated by simulation.
#' Results are stored to file, with the names of the results file being returned to this function
#' This file then saves the list of filenames in a .RDS format
#'
#' @return Lists of case names - this is the list of file names and identifiers from earlier storage
#' - resFile - Name of file where results are stored
#'
#' @export Run_jags.clist.cov
Run_jags.clist.cov <- function(caseList,  type=1, wnd.idx=1,
                           params=params.dep.det,
                           model_name=model_name.dep.det,
                           choleskyFlg=FALSE, h=NULL, glist=NULL,
                           ni=1500,nt=2,nb=500,nc=3,na=1000, debug=FALSE, # note testing values for ni, nb and nc
                           fid="",savedir="", store_filename_cases=FALSE ) {
  if ( is.null(type) || type <1 || type > length(Windat.type)) {
    print(paste0("Warning: type is wrong  must be an integer between 1 and ",length(Windat.type)))
    return(NULL)
  } else if ( is.null(wnd.idx) || wnd.idx < 1 || wnd.idx > 15) {
    print(paste0("Warning: win.idx is wrong  must be an integer between 1 and length of preps (1 if preps=NULL)"))
    return(NULL)
  }
  caseNames <- list()  # list of file names for saved jags output files
  for (case in 1:length(caseList)) {
    print(paste("Case",case,"of",length(caseList),sep=" "))
    caseData <- caseList[[case]]  # 2.1
    myCaseId <- caseData$id
    g=NULL
    if (isTRUE(choleskyFlg) ) {
      if ( length(glist)== 1 ) {g=glist[[1]]} else {g=glist[[case]]}
    }
    # create inits, and run jags with given model, params, and windata (caseData)
    occ_dat.jags <- runJags.case.cov(caseData,type=type,wnd.idx=wnd.idx,
                                 params=params,model_name=model_name,
                                 choleskyFlg=choleskyFlg,h=h,g=g,
                                 ni=ni,nt=nt,nb=nb,nc=nc,na=na, debug=debug)
    # now save it
    # get filename and add to the list
    caseNames[[case]] = saveJagsRes.cov(occ_dat.jags,nreps=caseData$windat[[wnd.idx]][[Windat.type[type] ]]$nrep,
                                        type=type, myId=myCaseId,mu_case=caseData$mu_case,
                                        p_case=caseData$p_case,fid=fid,choleskyFlg=choleskyFlg,savedir=savedir)
    print(paste0("Case ",case," saved as ",caseNames[[case]]$resFile,sep=" "))
    if ( store_filename_cases==TRUE) {
      if (is.null(caseData$windat[[wnd.idx]]$filename)) {caseData$windat[[wnd.idx]]$filename <- list()}
      caseData$windat[[wnd.idx]]$filename[[Windat.type[type]]] <- caseNames[[case]]$resFile
    }
  }# end of case list loop
  if ( store_filename_cases==TRUE ) {
    return(list(resNames = caseNames,resCases=caseList))
  }
  return(caseNames)
}

#################################################



