#************************************************************************************
#
# Functions to run for Multivariate Probit JSDM for real data sets
#
# independent: hierarchial model, probit occurrent, logit detection
# dependent: P(Z = 1) = P(V>0), V ~ MVN(mu, Sigma)
#
# Functions
# 1. windats_id.list  -  gets list of identifying data for real data
# 2. Run_jags.clist.cov.rl  -  Runs  MVP Jags model analysis over a list of real data sets
# 3. saveJagsRes.cov.rl     -  saves MVP Jags model analysis results for real data
#
#'
#' Run the MVP JSDM Jags model with Real Data
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `Run_jags.clist.cov.rl`,  `windats_id.list`, and `saveJagsRes.cov.rl` tools for
#' for fitting real data in the windats format to the Jags Multivariate Probit
#' Joint Species Distribution Model (MVP JSDM).
#'
#'
#' `Run_jags.clist.cov.rl()` runs jags JSDM MVP model for a list of covariate cases assuming real data.
#' `windats_id.list`   gets identifiers from real data to identify later results files
#' `saveJagsRes.cov.rl()` saves the results of real data: no simulation cases.
#'
#'
#' @section `Run_jags.clist.cov.rl()` Runs the MVP JSDM model for a list of dataset drawn from real data.
#' @param caseList        List of datasets in windat format
#' @param type            DAta type require set by integer values
#'     *     1 - prepares data for nsite X nrep X nspec species (or nspec+nz if nz !=0)
#'     *     2 - prepares collapsed data where y' = max(y[1,3]) (i.e. across replications)
#'     *     3 - random selection of replciation (but for real data, always choses first  set replications)
#' @param wnd.idx         Index to get dataset with specified number of survey replications
#' @param model_name      Length of the validation data (set to 0 if no training/validation)
#' @param choleskyFlg     If TRUE using Cholesky decomposition and uniform priors for covariance matrix
#'                        If FALSE assumes the model uses the Wishart prior for covariance matrix
#'                        Default is FALSE (assume Wishart or similar prior).  TRUE requires g and h to be set.
#' @param h               Number of degrees of freedom for wishart distribution, df=1 => flat, uninformative distribution
#' @param g               Number of degrees of freedom for wishart distribution, df=1 => flat, uninformative distribution
#' @param glist
#' @param ni              Number of interations for Bayesian fitting with Jags (see Jags and BUGs software)
#' @param nb              Number of interations of burn in for Bayesian fitting with Jags (see Jags and BUGs software)
#' @param nc              Number of chains for Bayesian fitting with Jags (see Jags and BUGs software)
#' @param na              Number of interations for adaption for Bayesian fitting with Jags (see Jags and BUGs software)
#' @param debug           Debug flag for Jags (default set to FALSE), if true runs in debug mode
#' @param fid             User provided unique Id string ensure any filenames to store data are unique
#' @param savedir         Directory to save results in (results are very large and best stored as files)
#' @param store_filename_cases  If true, returns stored case file names with input cases
#'                              If FALSE, just reutrns the stored case file name (same order as )
#'
# # running JAGS
# require(R2WinBUGS)
# require(rjags)
# require(jagsUI)
# # functions to support bayesian analysis
# require(arm)
# require(MCMCpack)
# # multivariate normal related functions and distributions
# require(MASS)
# require(mvtnorm)
# require(monomvn)
# # other regression requirements
# require(mclust)





#*******************************************************************
#
# Run_jags.clist.cov.rl
#
# Runs jags JSDM MVP model for a list of covariate cases assuming real data
# from Run_jags.clist.cov
# the main change is the saving the data:
# id is now changed it tells us about species and is currently simple string
# may add other post processing as confident

#'
#' @export Run_jags.clist.cov.rl

Run_jags.clist.cov.rl <- function(caseList,  type=1, wnd.idx=1,
                               params=params.dep.det,
                               model_name=model_name.dep.det,
                               choleskyFlg=FALSE, h=NULL, glist=NULL,
                               ni=1500,nt=2,nb=500,nc=3,na=1000, debug=FALSE, # note testing values for ni, nb and nc
                               fid="",savedir="", store_filename_cases=FALSE ) {
  # runJags.case <- function(obscase, type=1, win.idx=1,
  #                          params=params.nocov.dep,
  #                          model_name=model_name.nocov.dep,
  #                          ni=1500, nt=10, nb=500, nc=3, debug=FALSE ) {
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
    case_Data <- caseList[[case]]  # 2.1
    myCaseId <- case_Data$id#case_Data[[wnd.idx]][[type]]$id
    myCaseId <- ifelse ( is.null(myCaseId), paste0("real_case_",case),myCaseId )
    g=NULL
    if (isTRUE(choleskyFlg) ) {
      if ( length(glist)== 1 ) {g=glist[[1]]} else {g=glist[[case]]}
    }
    # create inits, and run jags with given model, params, and windata (caseData)
    # for study 2: must be creating data differently
    Case_Data_jags <- list ( id=fid, windat=case_Data)
    occ_dat.jags <- runJags.case.cov(Case_Data_jags,type=type,wnd.idx=wnd.idx,
                                     params=params,model_name=model_name,
                                     choleskyFlg=choleskyFlg,h=h,g=g,
                                     ni=ni,nt=nt,nb=nb,nc=nc,na=na, debug=debug)
    # now save it
    # get filename and add to the list
    nvalid= ifelse( is.null(case_Data$windat[[wnd.idx]]$other$y_valid),0,
                    case_Data$windat[[wnd.idx]]$other$y_valid$nsite )

    caseNames[[case]] = saveJagsRes.cov.rl(occ_dat.jags,myId=myCaseId,
                                           windat=case_Data$windat[[wnd.idx]][[Windat.type[type] ]],
                                           nvalid=nvalid,type=type,
                                           fid=fid,choleskyFlg=choleskyFlg,savedir=savedir)
    print(paste0("Case ",case," saved as ",caseNames[[case]]$resFile,sep=" "))
    if ( store_filename_cases==TRUE) {
      if (is.null(case_Data$windat[[wnd.idx]]$filename)) {case_Data$windat[[wnd.idx]]$filename <- list()}
      case_Data$windat[[wnd.idx]]$filename[[Windat.type[type]]] <- caseNames[[case]]$resFile
    }
  }# end of case list loop
  if ( store_filename_cases==TRUE ) {
    return(list(resNames = caseNames,resCases=caseList))
  }
  return(caseNames)
}





#*******************************************************************
#
# windats_id.list
#
# just trying to get data into right format for runJags.case.cov.rl
#' @section `windats_id.list()`
#' Gets windats list into correct format for Runjags by adding id data
#' Returns:
#'    - case.list - list of windats data, each entry containing id data
#' @export windats_id.list
windats_id.list <- function(windats.list,id.list) {
  case.list <- list();
  for ( i in 1:length(windats.list) ) {
    case.list[[i]] <- list(id=id.list[[i]],windat=windats.list[[i]]);
  }
  return(case.list)
}


#*******************************************************************
#
# saveJagsRes.cov.rl
#
# saving the results of real data: no simulation cases
#
#' @section `saveJagsRes.cov.rl()`
#' Saves results from real data when there are no simulation cases to file in RDS format.
#' File name is determined from identification info.
#' Returns:
#'    - resFile Filename of saved file,
#'    - resID   Identifying information
#' @export saveJagsRes.cov.rl
saveJagsRes.cov.rl <- function(jagsRes,myId,windat,nvalid=0, type=1,fid="",choleskyFlg=FALSE,savedir="") {
  #print(paste0("my id ",myId))
  ch_str = ifelse(isTRUE(choleskyFlg),"ch","")
  nspec=windat$nspec
  habvar=windat$novar.occ
  detvar = ifelse (is.null(windat$novar.det),0,windat$novar.det )
  nreps=windat$nrep
  Filename = paste("JagsRes","ID",myId,"tp",type,
                    "nspec",nspec,"nhabvar",habvar,"ndetvar",detvar,
                    ch_str,"rep",nreps,fid,sep="_")
  Filename_ext = paste0(Filename,".RDS")
  Fullfilename = paste0(savedir,Filename_ext,sep="")
  print(paste0("Full filename ",Fullfilename))


  # may want to extract summary of results to return, e.g.
  # tau, Mu (lpsi), p (lp), get median, mean, 95% and 5% results (look at others)
  id.dat = list()
  id.dat$idStr = myId
  id.dat$nspec = nspec
  id.dat$nsite = windat$nsite
  id.dat$nvalid = nvalid
  id.dat$ntrain = windat$nsite - nvalid
  id.dat$df = windat$df   # degrees of freedom used when estimating dwishart
  id.dat$nz = windat$nz   # missing species that can be added in
  id.dat$nreps = nreps
  id.dat$type = type
  id.dat$ch_flag = choleskyFlg
  id.dat$fid = fid  # series id

  saveRDS(list(jagsRes=jagsRes,resId=id.dat), file=Fullfilename )
  return(list(resFile=Filename,resID=myId))
}
