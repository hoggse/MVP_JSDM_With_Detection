#************************************************************************************
#
# Functions to process JSDM based on probit model for plotting and printing
#
#
# 1.  test_jagoutput         Checks at rhat and neff (autocorrelation values)
# 2.  getCorTauList          Calculate correlation for a list of Tau values
# 3.  getCor                 Calculate correlation for a single set of Tau values
# 4.  jagrespsforplot        Gets correlation estimates for plotting, and simple uncorrected beta/alpha intercept summary estimates
# 5.  jagcovrespsforplot     Gets correlation estimates for plotting, and simple uncorrected beta/alpha summary estimates for multiple covariates
# 6.  jagcovresps_beta_est   Gets corrected beta estimates for plotting (intercepts and slope parameters)
# 7.  printouts_rl1_det      summary data raw parameters (tau, beta, gamma, hyperparameter)
# 8.  printouts_rl2          summary data raw parameters (tau, beta, gamma, hyperparameter) and training sets for printing and to calculate ROC
# 9.  getOccDecProb          calc the ROC for each of a set of jagfiles
# 10. combineAlphaBetaJagsParamEstLists     combine corrected beta and alpha mean, median, q 2.5, q 97.5 estimates for lists of returns from jagcovrespsforplot and jagcovresps_beta_est
# 11. combineAlphaBetaJagsParamEsts         combine corrected beta and alpha mean, median, q 2.5, q 97.5 estimates for returns from jagcovrespsforplot and jagcovresps_beta_est
#
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `test_jagoutput()`,  `getCorTauList()`, `getCorTauList.multspec()`,`getCor()`,
#' `jagrespsforplot()`, `jagcovrespsforplot()`,
#' `jagcovresps_beta_est()`, `combineAlphaBetaJagsParamEstLists()`,`combineAlphaBetaJagsParamEsts()`
#'  `printouts_rl1_det()`,`printouts_rl2()`,
#' `getOccDecProb()`
#' provide tools for processing the reults of the JAGS analysis for the multivariate probit model
#'
#' `test_jagoutput()` Calculates rhat and effective n for a single jags response
#'
#' @section Common parameters used by jags data processing functions
#' @param jag_resp        Raw results data from jags
#' @param base_rhat_val   Base value for comparison
#' @param base_neff_value Base value for comparison
#' @param taulist         List of tau precision matrices
#' @param taures          Single tau precision matrix
#' @param list_jagfiles   List of jags results files
#' @param traceplot       Provide a traceplot if true, restricted to set of paramters of interst
#' @param savedircov      Directory to save analysis results
#' @param testFlag        Carry out test or not
#' @param test_params     Set of parameters to consider, given by indices
#' @param myBetaData      List of beta parameter estimates (occupancy)
#' @param myAlphaData     Lits of alpha parameter estiamtes (detection)
# setup real species observation data with occupancy and detection
# covariates for JAGS/BUGS JSDM MGP model
#

# packages required
# require(mvtnorm)
# require(ROCR)
# NB may need to download ROCR from http://ipa-tys.github.io/ROCR/ if you want to calculate AUC


#######################################
#
# test_jagoutput
#
# TESt output of jags run, looking for convergence and correlation found
# returns table flagging trouble some results
#' @section `test_jagoutput()`
#' Calculates rhat and effective n for a single jags response.
#' These are used to determine if the model is deemed to have converged and perform ok.
#' @return List of rhat and n.eff for each paramter in the model
#' @export test_jagoutput

test_jagoutput <- function (jag_resp,base_rhat_val=1.1, base_neff_value=1000, test_params=1:50) {
  converge=TRUE;
  CorrelationOk=TRUE;
  rhat.max=NULL;n.eff.min= NULL;
  if ( !is.null(test_params) ) { rhat_len = length(test_params) } else { rhat_len = length(jag_resp$summary[,"Rhat"])}
  for (i in 1:rhat_len ) {
    if ( is.null(jag_resp$summary[i,"Rhat"]) || is.na(jag_resp$summary[i,"Rhat"]) ){
      converge=FALSE ; print(paste0("Rhat ",i," is NULL or NA"))
    } else {
      if ( jag_resp$summary[i,"Rhat"] > base_rhat_val ) {converge=FALSE}  # can set this to preferred value
      if ( is.null(rhat.max) ) {rhat.max= jag_resp$summary[i,"Rhat"]} else {
        rhat.max=max(rhat.max,jag_resp$summary[i,"Rhat"])}
    }
    if ( is.null(jag_resp$summary[i,"n.eff"]) || is.na(jag_resp$summary[i,"n.eff"]) ){
      CorrelationOk=FALSE; print(paste0("n.eff ",i," is NULL or NA"))
    } else {
      if ( jag_resp$summary[i,"n.eff"] < base_neff_value ) {CorrelationOk=FALSE}
      if ( is.null(n.eff.min) ) {n.eff.min= jag_resp$summary[i,"n.eff"]} else {
        n.eff.min=min(n.eff.min,jag_resp$summary[i,"n.eff"])}
    }
  }
  print(paste("Converges",converge,"No Auto correlation",CorrelationOk,"Rhat max",round(rhat.max,4),"n.eff min",round(n.eff.min,4),sep=" "))
  return(list(convergeOk=converge, noAutoCor=CorrelationOk, rhat.max=rhat.max,n.eff.min=n.eff.min))
}


#######################################
#
# getCorTauList
#
# Calculate correlation for a list of Tau values for 2 species
# a single correlation value for the pair
# as returned by JAGS MVP model

# get correlation for list of Tau
# as in lists returned in jags_res$sims.list$Tau
#' @section `getCorTauList()`
#' Calcualtes the correlation matrices for a a list of Tau presision matrices000000
#' assuming only 2 species
#' @return List of correlation matrices
#' @export getCorTauList

getCorTauList <- function(taulist) {
  sigmalist <- apply(taulist,1,solve)
  dim(sigmalist) <- rev(dim(taulist))
  sigmalist <- aperm(sigmalist,c(3,2,1))
  estcorlist <- apply(sigmalist,1,cov2cor)
  dim(estcorlist) <- rev(dim(sigmalist))
  estcorlist <- aperm(estcorlist,c(3,2,1))
  nonmatch=0
  for ( i in dim(estcorlist)[1] ) {nonmatch = nonmatch + (estcorlist[i,1,2]==estcorlist[i,2,1])}
  print( paste("min cor",round(min(estcorlist[,2,1]),4),"max cor",round(max(estcorlist[,2,1]),4),
               "mean cor",round(mean(estcorlist[,2,1]),4),"median cor ",round(median(estcorlist[,2,1]),4),
               "q2.5",round(quantile(estcorlist[,2,1],0.025),4),"q97.5",round(quantile(estcorlist[,2,1],0.975),4),"cor sd",round(sd(estcorlist[,1,2]),4),
               "# nulls",sum(is.null(estcorlist[,2,1])),"# NAs",sum(is.na(estcorlist[,2,1])),"no match cor",nonmatch,sep=" " ) )
  return(list(est_rho=estcorlist[,2,1],
              est_mean_rho=mean(estcorlist[,2,1]),
              est_median_rho=median(estcorlist[,2,1]),
              est_q2.5_rho=quantile(estcorlist[,2,1],0.025),
              est_q97.5_rho=quantile(estcorlist[,2,1],0.975),
              est_sd_rho=sd(estcorlist[,1,2])
  ) )
}


#######################################
#
# getCorTauList.multspec
#
# Calculate correlation for a list of Tau values for >2 species
# multiple correlation values for each species pair
# as returned by JAGS MVP model

# get correlation for list of Tau
# as in lists returned in jags_res$sims.list$Tau
#' @section `getCorTauList.multspec()`
#' Calcualtes the correlation matrices for a a list of Tau presision matrices
#' for multiple species (> 2)
#' @return List of correlation matrices
#' @export getCorTauList.multspec

getCorTauList.multspec <- function(taulist) {
  nspec = dim(taulist)[2]
  sigmalist <- apply(taulist,1,solve)
  dim(sigmalist) <- rev(dim(taulist))
  sigmalist <- aperm(sigmalist,c(3,2,1))
  estcorlist <- apply(sigmalist,1,cov2cor)
  dim(estcorlist) <- rev(dim(sigmalist))
  estcorlist <- aperm(estcorlist,c(3,2,1))
  var.est <- array(dim=c(nspec,6))
  covar.est <- array(dim=c(nspec,nspec,6))
  corr.est <- array(dim=c(nspec,nspec,6))
  names <- c("mean","median","q2.5","q97.5","min","max")
  colnames(var.est) <- names
  j=1;k=2
  for ( j in 1:nspec ) {
    # diagnals are the variance
    var.est[j,] = cbind(mean(sigmalist[,j,j]), median(sigmalist[,j,j]),quantile(sigmalist[,j,j],0.025),quantile(sigmalist[,j,j],0.975),min(sigmalist[,j,j]),max(sigmalist[,j,j]))
    covar.est[j,j,] = var.est[j,]
    corr.est[j,j,] = 1
    for ( k in (j+1):nspec) {
      if ( j< nspec ) {
        covar.est[j,k,] = cbind(mean(sigmalist[,j,k]), median(sigmalist[,j,k]),quantile(sigmalist[,j,k],0.025),quantile(sigmalist[,j,k],0.975),min(sigmalist[,j,k]),max(sigmalist[,j,k]) )
        corr.est[j,k,] = cbind(mean(estcorlist[,j,k]), median(estcorlist[,j,k]),quantile(estcorlist[,j,k],0.025),quantile(estcorlist[,j,k],0.975),min(estcorlist[,j,k]),max(estcorlist[,j,k]) )
        covar.est[k,j,] = covar.est[j,k,]
        corr.est[k,j,] = corr.est[j,k,]
      }
    }
  }
  return(list(var.est=var.est,
              covar.est=covar.est,
              rho.est=corr.est,
              estcorlist=estcorlist,
              estsigmalist=sigmalist )
  )
}



#######################################
#
# getCor
#
# for correlation matrix for single values of Tau
# e.g. jags_res$mean$Tau
#' @section `getCor()`
#' Calcualtes the correlation matrix for a single Tau presision matrix
#' @return correlation matrix
#' @export getCor
getCor <- function(tau_res) {
  sigma_res <- solve(tau_res)
  cor_res <- cov2cor(sigma_res)
  return(cor_res)
}


########################################
#
# jagrespsforplot
#
# gets correlation data and intercepts:
# note beta intercepts not standardised and should not be used for simulation comparisons
#' @section `jagrespsforplot()`
#' Extracts parameter estimates from jags files including correlation, and alpha parameters
#' Beta parameters are not standarised.
#' @return List of parameter estimates for each jags files
#' list includes: mean, 97.5% percentile, 2.5$ percentile, 50% percentile for each paramter
#' plus list of all calculated correlation estimates
#' @export jagrespsforplot

jagrespsforplot <- function(list_jagfiles,traceplot=TRUE,savedircov="",testFlag=TRUE, test_params=1:50) {
  responses <- list()
  for ( i in 1:length(list_jagfiles)) {
    file_base = list_jagfiles[[i]][[1]]
    jags_file <- readRDS(paste0(savedircov,file_base) )
    if (testFlag ) {
      test_res <- test_jagoutput(jags_file[[1]],test_params=test_params)
    } else { test_res=NULL}
    rho_mean.mx <- getCor(jags_file[[1]]$mean$Tau)
    rho_q97.5.mx <- getCor(jags_file[[1]]$q97.5$Tau)
    rho_q02.5.mx <- getCor(jags_file[[1]]$q2.5$Tau)
    given_rho = jags_file$resId$rho
    if ( dim(jags_file[[1]]$sims.list$Tau)[2] == 2 ) { # single correlation
      rho_list <- getCorTauList(jags_file[[1]]$sims.list$Tau)
    } else {  # multiple correlations
      rho_list <- getCorTauList.multspec(jags_file[[1]]$sims.list$Tau)
    }
    given_mu = jags_file$resId$mu
    mu_mean = jags_file[[1]]$mean$lpsi
    mu_q97.5 = jags_file[[1]]$q97.5$lpsi
    mu_q02.5 = jags_file[[1]]$q2.5$lpsi
    given_detp = jags_file$resId$detp
    detp_mean = inv.logit(jags_file[[1]]$mean$lp)
    detp_q97.5 = inv.logit(jags_file[[1]]$q97.5$lp)
    detp_q02.5 = inv.logit(jags_file[[1]]$q2.5$lp)
    # Perform trace plots if requested
    # return response data in in list of easy to use matrices
    file_png <- gsub(".RDS","traces.PNG",file_base)
    paste0(savedircov,file_png)
    if ("lp"  %in% jags_file[[1]]$parameters ) {
      if( isTRUE(traceplot) ) {
        png(file_png)
        par(mfrow=c(4,3));par(ask=F) # works better for file
        jagsUI::traceplot(jags_file$jagsRes ,c("Tau","lpsi","mu.lpsi","sd.lpsi","lp","mu.lp","sd.lp") )
        dev.off()
      }
      responses[[i]] = list(resid=jags_file$resId,test_res=test_res,rho_mean=rho_mean.mx,rho_list=rho_list,
                            given_rho=given_rho, est_mean_rho=rho_list$est_mean_rho, est_q2.5_rho=rho_list$est_q2.5_rho, est_q97.5_rho=rho_list$est_q97.5_rho,
                            given_mu=given_mu, mu_mean=mu_mean, mu_q97.5=mu_q97.5, mu_q02.5=mu_q02.5,
                            given_detp=given_detp, detp_mean=detp_mean, detp_q97.5=detp_q97.5, detp_q02.5=detp_q02.5,
                            trace_file=file_png )
    } else {
      if(isTRUE(traceplot)) {
        png(file_png)
        par(mfrow=c(3,3));par(ask=F)
        jagsUI::traceplot(jags_file$jagsRes ,c( "Tau","lpsi","mu.lpsi","sd.lpsi") )
        dev.off()
      }
      responses[[i]] = list(resid=jags_file$resId,test_res=test_res,rho_mean=rho_mean.mx,rho_list,
                            given_rho=given_rho, est_mean_rho=rho_list$est_mean_rho, est_q2.5_rho=rho_list$est_q2.5_rho, est_q97.5_rho=rho_list$est_q97.5_rho,
                            given_mu=given_mu, mu_mean=mu_mean, mu_q97.5=mu_q97.5, mu_q02.5=mu_q02.5,
                            trace_file=file_png )
    }
  }
  return(responses)
}


######################################
#
# jagcovrespsforplot
#
# returns correlation with raw beta/alpha estimates
# note beta estimate have not been standardised
# so should not be used for comparison purposes for simulation parameters
# also returns site esimates of mean occupance and probability of detection
#' @section `jagcovrespsforplot()`
#' Extracts parameter estiamtes from jags files including Covariance, correlation, beta and alpha parameters
#' Beta parameters are not standarised.
#' @return List of parameter estimates for each jags files
#' list includes: mean, 97.5% percentile, 2.5$ percentile, 50% percentile for each paramter
#' @export jagcovrespsforplot

jagcovrespsforplot <- function(list_jagfiles,traceplot=FALSE,savedircov,testFlag=FALSE, test_params=1:50) {
  responses <- list()
  for ( i in 1:length(list_jagfiles)) {
    file_base = list_jagfiles[[i]][[1]]
    jags_file <- readRDS(paste0(savedircov,file_base) )
    # get errors
    if (testFlag )test_res <- test_jagoutput(jags_file[[1]],test_params=test_params) else test_res=NULL
    # correlation parameters
    rho_mean.mx <- getCor(jags_file[[1]]$mean$Tau)
    rho_med.mx <- getCor(jags_file[[1]]$q50$Tau)
    rho_q97.5.mx <- getCor(jags_file[[1]]$q97.5$Tau)
    rho_q02.5.mx <- getCor(jags_file[[1]]$q2.5$Tau)
    sigma_mean <- solve(jags_file[[1]]$mean$Tau)
    sigma_med <- solve(jags_file[[1]]$q50$Tau)
    sigma_q97.5 <- solve(jags_file[[1]]$q97.5$Tau)
    sigma_q02.5 <- solve(jags_file[[1]]$q2.5$Tau)
    # given parameter: resolve to NULL if real data or not provided
    given_rho = jags_file$resId$rho
    calc_vcor = jags_file$resId$V_cor
    calc_zcor = jags_file$resId$Z_cor
    if ( dim(jags_file[[1]]$sims.list$Tau)[2] == 2 ) {
      rho_list <- getCorTauList(jags_file[[1]]$sims.list$Tau)
    } else {
      rho_list <- getCorTauList.multspec(jags_file[[1]]$sims.list$Tau)
    }

    # beta data (raw, no correction)
    nspec = length(jags_file[[1]]$mean$lpsi)
    if ( !is.null(jags_file[[1]]$mean$betalpsi) ) {
      nparam = ifelse(is.matrix(jags_file[[1]]$mean$betalpsi),dim(jags_file[[1]]$mean$betalpsi)[2],ifelse (nspec==1,length(jags_file[[1]]$mean$betalpsi),1) )
    } else {nparam=0}
    # if there is given data
    if ( !is.null(jags_file$resId$beta1) ){
      if ( is.matrix(jags_file$resId$beta1)) {
        for ( p_cnt in 1:nparam) {
          if ( p_cnt == 1) {
            beta_given <- rbind(jags_file$resId$beta0,jags_file$resId$beta1[,p_cnt])
          } else {
            beta_given <- rbind(beta_given,jags_file[[1]]$resId$beta1[,p_cnt])
          }
        }
      } else {
        if (nspec == 1) {
          beta_given <-matrix(c(jags_file$resId$beta0,jags_file$resId$beta1),ncol=1)
        } else {
          beta_given <- rbind(jags_file$resId$beta0,jags_file$resId$beta1)
        }
      }  #otherwise, sets betagive to intercepts only or NULL if no given data
    } else if  ( !is.null(jags_file$resId$beta0) ) {
      beta_given <- matrix(jags_file$resId$beta0,nrow=1)
    } else {beta_given = NULL}

    # now set up the estimated beta data (raw, not corrected)
    # if use raw data for predictions need to use with covariance data
    if ( !is.null(jags_file[[1]]$mean$betalpsi) ) { #if there are covariates
      # matrix of covariates
      if ( is.matrix(jags_file[[1]]$mean$betalpsi) ) { # matrix of data
        for ( p_cnt in 1:nparam) {
          if ( p_cnt == 1) {
            beta_mean <- rbind(jags_file[[1]]$mean$lpsi,jags_file[[1]]$mean$betalpsi[,p_cnt])
            beta_med <- rbind(jags_file[[1]]$q50$lpsi,jags_file[[1]]$q50$betalpsi[,p_cnt])
            beta_q97.5 <- rbind(jags_file[[1]]$q97.5$lpsi,jags_file[[1]]$q97.5$betalpsi[,p_cnt])
            beta_q02.5 <- rbind(jags_file[[1]]$q2.5$lpsi,jags_file[[1]]$q2.5$betalpsi[,p_cnt])
            beta_sd    <- rbind(jags_file[[1]]$sd$lpsi,jags_file[[1]]$sd$betalpsi[,p_cnt])
          } else {
            beta_mean <- rbind(beta_mean,jags_file[[1]]$mean$betalpsi[,p_cnt])
            beta_med  <- rbind(beta_med,jags_file[[1]]$q50$betalpsi[,p_cnt])
            beta_q97.5 <- rbind(beta_q97.5,jags_file[[1]]$q97.5$betalpsi[,p_cnt])
            beta_q02.5 <- rbind(beta_q02.5,jags_file[[1]]$q2.5$betalpsi[,p_cnt])
            beta_sd    <- rbind(beta_sd,jags_file[[1]]$sd$betalpsi[,p_cnt])
          }
        }# for each covariate
      } else {  # otherwise must be a vector of covariate estimates
        if ( nspec == 1) { # multiple covariates and a single specvies
          beta_mean = matrix(c(jags_file[[1]]$mean$lpsi,jags_file[[1]]$mean$betalpsi),ncol=1)
          beta_med  = matrix(c(jags_file[[1]]$q50$lpsi,jags_file[[1]]$q50$betalpsi),ncol=1)
          beta_q97.5 = matrix(c(jags_file[[1]]$q97.5$lpsi,jags_file[[1]]$q97.5$betalpsi),ncol=1)
          beta_q02.5 = matrix(c(jags_file[[1]]$q2.5$lpsi,jags_file[[1]]$q2.5$betalpsi),ncol=1)
          beta_sd    = matrix(c(jags_file[[1]]$sd$lpsi,jags_file[[1]]$sd$betalpsi),ncol=1)
        } else {  # otherwise multiple species but a single parameter
          beta_mean <- rbind(jags_file[[1]]$mean$lpsi,jags_file[[1]]$mean$betalpsi)
          beta_med <-  rbind(jags_file[[1]]$q50$lpsi,jags_file[[1]]$q50$betalpsi)
          beta_q97.5 <- rbind(jags_file[[1]]$q97.5$lpsi,jags_file[[1]]$q97.5$betalpsi)
          beta_q02.5 <- rbind(jags_file[[1]]$q2.5$lpsi,jags_file[[1]]$q2.5$betalpsi)
          beta_sd <- rbind(jags_file[[1]]$sd$lpsi,jags_file[[1]]$sd$betalpsi)
        }
      }# end matrix or vector for lpsi
    } else { # no covariates: set beta estimates to intercepts only
      beta_mean <- matrix(jags_file[[1]]$mean$lpsi,nrow=1)
      beta_med  = matrix(jags_file[[1]]$q50$lpsi,nrow=1)
      beta_q97.5 = matrix(jags_file[[1]]$q97.5$lpsi,nrow=1)
      beta_q02.5 = matrix(jags_file[[1]]$q2.5$lpsi,nrow=1)
      beta_sd    = matrix(jags_file[[1]]$sd$lpsi,nrow=1)
    }# end setting up beta estimates

    # alpha parameters
    # given parameter values (if relevant)
    if ( !is.null(jags_file[[1]]$mean$alphalp) ) {
      nparam.det = ifelse(is.matrix(jags_file[[1]]$mean$alphalp),dim(jags_file[[1]]$mean$alphalp)[2],ifelse (nspec==1,length(jags_file[[1]]$mean$alphalp),1) )
    } else {nparam.det=0}
    # if there is given data
    if ( !is.null(jags_file$resId$alpha1) ){
      if ( is.matrix(jags_file$resId$alpha1)) {
        for ( p_cnt in 1:nparam.det) {
          if ( p_cnt == 1) {
            alpha_given <- rbind(jags_file$resId$alpha0,jags_file$resId$alpha1[,p_cnt])
          } else {
            alpha_given <- rbind(alpha_given,jags_file[[1]]$resId$alpha1[,p_cnt])
          }
        }
      } else {
        if (nspec == 1) {
          alpha_given <-matrix(c(jags_file$resId$alpha0,jags_file$resId$alpha1),ncol=1)
        } else {
          alpha_given <- rbind(jags_file$resId$alpha0,jags_file$resId$alpha1)
        }
      }  #otherwise, sets betagive to intercepts only or NULL if no given data
    } else if  ( !is.null(jags_file$resId$alpha0) ) {
      alpha_given <- matrix(jags_file$resId$alpha0,nrow=1)
    } else {alpha_given = NULL}

    if ("lp"  %in% jags_file[[1]]$parameters ) {
      if ( !is.null(jags_file[[1]]$mean$alphalp) ) { #if there are covariates
        # matrix of covariates
        if ( is.matrix(jags_file[[1]]$mean$alphalp) ) { # matrix of data
          for ( p_cnt in 1:nparam.det) {
            if ( p_cnt == 1) {
              alpha_mean <- rbind(jags_file[[1]]$mean$lp,jags_file[[1]]$mean$alphalp[,p_cnt])
              alpha_med <- rbind(jags_file[[1]]$q50$lp,jags_file[[1]]$q50$alphalp[,p_cnt])
              alpha_q97.5 <- rbind(jags_file[[1]]$q97.5$lp,jags_file[[1]]$q97.5$alphalp[,p_cnt])
              alpha_q02.5 <- rbind(jags_file[[1]]$q2.5$lp,jags_file[[1]]$q2.5$alphalp[,p_cnt])
              alpha_sd <- rbind(jags_file[[1]]$sd$lp,jags_file[[1]]$sd$alphalp[,p_cnt])
            } else {
              alpha_mean <- rbind(alpha_mean,jags_file[[1]]$mean$alphalp[,p_cnt])
              beta_med  <- rbind(alpha_med,jags_file[[1]]$q50$alphalp[,p_cnt])
              alpha_q97.5 <- rbind(alpha_q97.5,jags_file[[1]]$q97.5$alphalp[,p_cnt])
              alpha_q02.5 <- rbind(alpha_q02.5,jags_file[[1]]$q2.5$alphalp[,p_cnt])
              alpha_sd <- rbind(alpha_sd,jags_file[[1]]$sd$alphalp[,p_cnt])
            }
          }# for each covariate
        } else {  # otherwise must be a vector of covariate estimates
          if ( nspec == 1) { # multiple covariates and a single specvies
            alpha_mean = matrix(c(jags_file[[1]]$mean$lp,jags_file[[1]]$mean$alphalp),ncol=1)
            alpha_med  = matrix(c(jags_file[[1]]$q50$lp,jags_file[[1]]$q50$alphalp),ncol=1)
            alpha_q97.5 = matrix(c(jags_file[[1]]$q97.5$lp,jags_file[[1]]$q97.5$alphalp),ncol=1)
            alpha_q02.5 = matrix(c(jags_file[[1]]$q2.5$lp,jags_file[[1]]$q2.5$alphalp),ncol=1)
            alpha_sd    = matrix(c(jags_file[[1]]$sd$lp,jags_file[[1]]$sd$alphalp),ncol=1)
          } else {  # otherwise multiple species but a single parameter
            alpha_mean <- rbind(jags_file[[1]]$mean$lp,jags_file[[1]]$mean$betalpsi)
            alpha_med <-  rbind(jags_file[[1]]$q50$lp,jags_file[[1]]$q50$betalpsi)
            alpha_q97.5 <- rbind(jags_file[[1]]$q97.5$lp,jags_file[[1]]$q97.5$betalpsi)
            alpha_q02.5 <- rbind(jags_file[[1]]$q2.5$lp,jags_file[[1]]$q2.5$betalpsi)
            alpha_sd    <- rbind(jags_file[[1]]$sd$lp,jags_file[[1]]$sd$betalpsi)
          }# end if nspec > 1
        }# if matrix/vector for alphlp
      } else { # no covariates: set beta estimates to intercepts only
        alpha_mean <- matrix(jags_file[[1]]$mean$lp,nrow=1)
        alpha_med  = matrix(jags_file[[1]]$q50$lp,nrow=1)
        alpha_q97.5 = matrix(jags_file[[1]]$q97.5$lp,nrow=1)
        alpha_q02.5 = matrix(jags_file[[1]]$q2.5$lp,nrow=1)
        alpha_sd    = matrix(jags_file[[1]]$sd$lp,nrow=1)
      }# end setting up alpha estimates
    }# end of if lp estimates available (not available if no explicit detection)

    # not getting the hyper params here
    # get the recorded values
    # first the mean occupancy estimates per site/species
    mu_mean = jags_file[[1]]$mean$mu.occ  # average for each site for each species
    mu_med = jags_file[[1]]$q50$mu.occ  # average for each site for each species
    mu_q97.5 = jags_file[[1]]$q97.5$mu.occ
    mu_q02.5 = jags_file[[1]]$q2.5$mu.occ
    # first the detection probability estimates per site/species
    if ("p"  %in% jags_file[[1]]$parameters ) {
      detp_mean = inv.logit(jags_file[[1]]$mean$p)  # average for each site for each species
      detp_med = inv.logit(jags_file[[1]]$q50$p)  # average for each site for each species
      detp_q97.5 = inv.logit(jags_file[[1]]$q97.5$p)
      detp_q02.5 = inv.logit(jags_file[[1]]$q2.5$p)
    }

    # Perform trace plots if requested
    # return response data in in list of easy to use matrices
    file_png <- gsub(".RDS","traces.PNG",file_base)
    paste0(savedircov,file_png)
    # get the probility of detection per site/species
    if ("lp"  %in% jags_file[[1]]$parameters ) {
      if( isTRUE(traceplot) ) {
        png(file_png)
        par(mfrow=c(4,3));par(ask=F) # works better for file
        jagsUI::traceplot(jags_file$jagsRes ,c("Tau","lpsi","mu.lpsi","sd.lpsi","lp","mu.lp","sd.lp") )
        dev.off()
      }
      if ("p"  %in% jags_file[[1]]$parameters ) {  # if detp, p etc are returned
        responses[[i]] = list(resid=jags_file$resId,test_res=test_res,rho_list=rho_list,
                              given_rho=given_rho, est_mean_rho=rho_list$est_mean_rho, est_q2.5_rho=rho_list$est_q2.5_rho, est_q97.5_rho=rho_list$est_q97.5_rho,
                              rho_mean=rho_mean.mx,rho_med.mx=rho_med.mx,rho_q97.5.mx=rho_q97.5.mx,rho_q02.5.mx=rho_q02.5.mx,
                              sigma_mean=sigma_mean,sigma_med=sigma_med,sigma_q97.5=sigma_q97.5,sigma_q02.5=sigma_q02.5,
                              beta_given=beta_given,beta_mean=beta_mean,beta_med=beta_med,beta_q97.5=beta_q97.5,beta_q02.5=beta_q02.5,beta_sd=beta_sd,
                              alpha_given=alpha_given,alpha_mean=alpha_mean,alpha_med=alpha_med,alpha_q97.5=alpha_q97.5,alpha_q02.5=alpha_q02.5,alpha_sd=alpha_sd,
                              mu_mean=mu_mean, mu_med=mu_med, mu_q97.5=mu_q97.5, mu_q02.5=mu_q02.5,
                              detp_mean=detp_mean, detp_med=detp_med,detp_q97.5=detp_q97.5, detp_q02.5=detp_q02.5,
                              trace_file=file_png )
      } else { # if p and detp not returned, cut down set returned
        responses[[i]] = list(resid=jags_file$resId,test_res=test_res,rho_list=rho_list,
                              given_rho=given_rho, est_mean_rho=rho_list$est_mean_rho, est_q2.5_rho=rho_list$est_q2.5_rho, est_q97.5_rho=rho_list$est_q97.5_rho,
                              rho_mean=rho_mean.mx,rho_med.mx=rho_med.mx,rho_q97.5.mx=rho_q97.5.mx,rho_q02.5.mx=rho_q02.5.mx,
                              sigma_mean=sigma_mean,sigma_med=sigma_med,sigma_q97.5=sigma_q97.5,sigma_q02.5=sigma_q02.5,
                              beta_given=beta_given,beta_mean=beta_mean,beta_med=beta_med,beta_q97.5=beta_q97.5,beta_q02.5=beta_q02.5,beta_sd=beta_sd,
                              alpha_given=alpha_given,alpha_mean=alpha_mean,alpha_med=alpha_med,alpha_q97.5=alpha_q97.5,alpha_q02.5=alpha_q02.5,alpha_sd=alpha_sd,
                              mu_mean=mu_mean, mu_med=mu_med, mu_q97.5=mu_q97.5, mu_q02.5=mu_q02.5,
                              #detp_mean=detp_mean, detp_med=detp_med,detp_q97.5=detp_q97.5, detp_q02.5=detp_q02.5,
                              trace_file=file_png )
      }
    } else {
      if(isTRUE(traceplot)) {
        png(file_png)
        par(mfrow=c(3,3));par(ask=F)
        jagsUI::traceplot(jags_file$jagsRes ,c( "Tau","lpsi","mu.lpsi","sd.lpsi") )
        dev.off()
      }
      responses[[i]] = list(resid=jags_file$resId,test_res=test_res,rho_list=rho_list,
                            given_rho=given_rho, est_mean_rho=rho_list$est_mean_rho, est_q2.5_rho=rho_list$est_q2.5_rho, est_q97.5_rho=rho_list$est_q97.5_rho,
                            rho_mean=rho_mean.mx,rho_med.mx=rho_med.mx,rho_q97.5.mx=rho_q97.5.mx,rho_q02.5.mx=rho_q02.5.mx,
                            sigma_mean=sigma_mean,sigma_med=sigma_med,sigma_q97.5=sigma_q97.5,sigma_q02.5=sigma_q02.5,
                            beta_given=beta_given,beta_mean=beta_mean,beta_med=beta_med,beta_q97.5=beta_q97.5,beta_q02.5=beta_q02.5,beta_sd=beta_sd,
                            mu_mean=mu_mean, mu_med=mu_med, mu_q97.5=mu_q97.5, mu_q02.5=mu_q02.5,
                            trace_file=file_png )
    }
  }
  return(responses)
}


######################################
#
# jagcovresps_beta_est
#
# returns standardised beta parameter estimates
# This data is required for meaningful comparisons
#' @section `jagcovresps_beta_est()`
#' Extracts beta parameter estimates for a list of jag results files
#' Beta estimates are standardised using the sd from estimated covariance matrix
#' @return List of beta parameter estimates for each jags file
#'
#' @export jagcovresps_beta_est

jagcovresps_beta_est <- function(list_jagfiles, savedircov) {
  if (is.null(list_jagfiles) ) {return(NULL)}
  beta_responses <- list()
  for ( i in 1:length(list_jagfiles)) {
    file_base = list_jagfiles[[i]][[1]]
    jags_file <- readRDS(paste0(savedircov,file_base) )
    # get sigma_list (need the diagonal values)
    Sigma2 <- apply(jags_file[[1]]$sims.list$Tau, 1, solve)
    dim(Sigma2) <- rev(dim(jags_file[[1]]$sims.list$Tau))
    Sigma2 <- aperm(Sigma2, c(3, 2, 1))
    Diag_Sigma_2 <- apply(Sigma2,1,diag)

    # get standardised beta estimates from raw beta estimates
    beta_raw_list <- jags_file[[1]]$sims.list$betalpsi #  sims.list$Beta
    lpsi_raw_list <- jags_file[[1]]$sims.list$lpsi
    beta_list <- apply(beta_raw_list, 3, function(x) x / t(sqrt(Diag_Sigma_2)))
    dim(beta_list) <- dim(beta_raw_list)
    #dim(beta_list);length(beta_list)
    lpsi_list <- lpsi_raw_list / t(sqrt(Diag_Sigma_2))
    dim(lpsi_list) <- dim(lpsi_raw_list)
    #length(lpsi_list);dim(lpsi_list)

    # find statistics for standardised slope parameters
    est_mean_beta = apply(beta_list,c(2,3),mean)
    est_median_beta = apply(beta_list,c(2,3),median)
    est_q2.5_beta = apply(beta_list,c(2,3),quantile,0.025)
    est_q97.5_beta = apply(beta_list,c(2,3),quantile,0.975)
    est_sd_beta=apply(beta_list,c(2,3),sd)
    est_mean = apply(beta_list,c(2),mean)
    est_mean = apply(beta_list,c(1),mean)
    # dim(est_mean)

    # find statistics for standardised intercept parameters
    est_mean_beta0 = apply(lpsi_list,2,mean)
    est_median_beta0 = apply(lpsi_list,2,median)
    est_q2.5_beta0 = apply(lpsi_list,2,quantile,0.025)
    est_q97.5_beta0 = apply(lpsi_list,2,quantile,0.975)
    est_sd_beta0=apply(lpsi_list,2,sd)


    # put summary data into matrix format
    nspec = length(est_mean_beta0)
    nparam = ifelse(is.matrix(est_mean_beta),dim(est_mean_beta)[2],ifelse (nspec==1,length(est_mean_beta),1) )
    for ( p_cnt in 1:nparam) {
      if (p_cnt == 1 ) { # initial spec_cnt - setup betaAll
        if ( is.matrix(est_mean_beta)) {  # beta is a matrix
          est_mean_betaAll=rbind(est_mean_beta0,est_mean_beta[,p_cnt])
          est_median_betaAll=rbind(est_median_beta0,est_median_beta[,p_cnt])
          est_q2.5_betaAll=rbind(est_q2.5_beta0,est_q2.5_beta[,p_cnt])
          est_q97.5_betaAll=rbind(est_q97.5_beta0,est_q97.5_beta[,p_cnt])
          est_sd_betaAll=rbind(est_sd_beta0,est_sd_beta[,p_cnt])
        } else { # o/w beta is a vector
          est_mean_betaAll=rbind(est_mean_beta0,est_mean_beta[p_cnt])
          est_median_betaAll=rbind(est_median_beta0,est_median_beta[p_cnt])
          est_q2.5_betaAll=rbind(est_q2.5_beta0,est_q2.5_beta[p_cnt])
          est_q97.5_betaAll=rbind(est_q97.5_beta0,est_q97.5_beta[p_cnt])
          est_sd_betaAll=rbind(est_sd_beta0,est_sd_beta[p_cnt])
        }# end not matrix
      } else { # o/w spec_cnt > 1
        if ( is.matrix(est_mean_beta)) {  # beta is a matrix?  (K>1)
          est_mean_betaAll=rbind(est_mean_betaAll,est_mean_beta[,p_cnt])
          est_median_betaAll=rbind(est_median_betaAll,est_median_beta[,p_cnt])
          est_q2.5_betaAll=rbind(est_q2.5_betaAll,est_q2.5_beta[,p_cnt])
          est_q97.5_betaAll=rbind(est_q97.5_betaAll,est_q97.5_beta[,p_cnt])
          est_sd_betaAll=rbind(est_sd_betaAll,est_sd_beta[,p_cnt])
        } else { # o/w is a vector
          est_mean_betaAll=rbind(est_mean_betaAll,est_mean_beta[p_cnt])
          est_median_betaAll=rbind(est_median_betaAll,est_median_beta[p_cnt])
          est_q2.5_betaAll=rbind(est_q2.5_betaAll,est_q2.5_beta[p_cnt])
          est_q97.5_betaAll=rbind(est_q97.5_betaAll,est_q97.5_beta[p_cnt])
          est_sd_betaAll=rbind(est_sd_betaAll,est_sd_beta[p_cnt])
        }# end is a vector
      }# end spec_cnt != 1
    }# end for loop
    beta_responses[[i]] <- list(est_mean_beta=est_mean_betaAll,
                                est_median_beta=est_q2.5_betaAll,
                                est_q2.5_beta=est_q2.5_betaAll,
                                est_q97.5_beta=est_q97.5_betaAll,
                                est_sd_beta=est_sd_betaAll
    )
  }
  return(beta_responses)
}


######################################
#
# combineAlphaBetaJagsParamEstLists
#
# combine the corrected beta parameter estiamtes with alpha estimates for plotting
# this funciton will run the extraction over lists of returns from:
# jagcovresps_beta_est - beta data
# jagcovrespsforplot - alpha data
# Allows for comparing different lists of responses
#
#' @section `combineAlphaBetaJagsParamEstLists()`
#' Combines previously extracted lists of beta and alpha paramters
#' @return List of combined set of lists
#'
#' @export combineAlphaBetaJagsParamEstLists

combineAlphaBetaJagsParamEstLists <- function(myBetaData, myAlphaData) {
  if ( is.null(myBetaData) || is.null(myAlphaData) ) {
    print(paste0("combineAlphaBetaJagsParamEsts: ", "Missing beta or alpha input data"))
    return(NULL)
  }
  if ( length(myBetaData) != length(myAlphaData) ) {
    print(paste0("combineAlphaBetaJagsParamEsts: ", "Beta alpha list mismatch"))
    return(NULL)
  }
  myCombData <- list()
  for ( idx in 1:length(myBetaData) ) {
    myCombData[[idx]] <- combineAlphaBetaJagsParamEsts(myBetaData[[idx]],myAlphaData[[idx]])
  }
  return(myCombData)
}


######################################
#
# combineAlphaBetaJagsParamEsts
#
# combine the corrected beta parameter estimates with alpha estimates for plotting
# this funciton will run the extraction over returns from:
# jagcovresps_beta_est - beta data
# jagcovrespsforplot - alpha data
#
#' @section `combineAlphaBetaJagsParamEsts()`
#' Combines previously extracted lists of beta and alpha paramters
#' @return combined set of lists
#'
#' @export combineAlphaBetaJagsParamEsts

combineAlphaBetaJagsParamEsts <- function(myBetaData, myAlphaData) {
  if ( is.null(myBetaData) || is.null(myAlphaData) ) {
    print(paste0("combineAlphaBetaJagsParamEsts: ", "Missing beta or alpha input data"))
    return(NULL)
  }
  if ( length(myBetaData) != length(myAlphaData) ) {
    print(paste0("combineAlphaBetaJagsParamEsts: ", "Beta alpha list mismatch"))
    return(NULL)
  }
  myCombData <- myBetaData
  for ( idx in 1:length(myCombData) ) {
    myCombData[[idx]]$alpha_given      <- myAlphaData[[idx]]$alpha_given
    myCombData[[idx]]$est_mean_alpha    <- myAlphaData[[idx]]$alpha_mean
    myCombData[[idx]]$est_med_alpha        <- myAlphaData[[idx]]$alpha_med
    myCombData[[idx]]$est_q97.5_alpha        <- myAlphaData[[idx]]$alpha_q97.5
    myCombData[[idx]]$est_q02.5_alpha        <- myAlphaData[[idx]]$alpha_q02.5
    myCombData[[idx]]$est_sd_alpha        <- myAlphaData[[idx]]$alpha_sd
  }
  return(myCombData)
}



###########################################
#
# printouts_rl1_det
#
# summary data raw parameters (tau, beta, gamma, hyperparameter) for printing
#
#' @section `printouts_rl1_det()`
#' Prints out raw parameter estimates for each entry in a list of jags result files.
#' @return A list of estimates for each file including
#' - DIC
#' - Beta parameters
#' - Alpha parameters
#' - Tau parameters (inverse of covariance matrix, Tau is estiamted by JAGS)
#'
#' @export printouts_rl1_det

printouts_rl1_det <- function (list_jagfiles,windat.list,wnd.idx=1,savedir="") {
  reslen <- length(list_jagfiles)
  param_sums <- list()
  for ( res in 1:reslen) {
    file_base = list_jagfiles[[res]][[1]]
    jags_res <- readRDS(paste0(savedir,file_base) )
    resid <- jags_res$resId
    nspec <- resid$nspec
    type <- resid$type
    windat <- windat.list[[1]]$windat[[wnd.idx]][[Windat.type[type]]]

    noccvar <- windat$novar.occ
    ndetvar <- ifelse( is.null(windat$novar.det), 0, windat$novar.det )
    resid$noccvar <- noccvar
    resid$ndetvar <- ndetvar

    # get precision
    tau_sum <- list()
    tau_sum[["mean"]] = jags_res$jagsRes$mean[["Tau"]]
    tau_sum[["median"]] = jags_res$jagsRes$q50[["Tau"]]
    tau_sum[["q2.5"]] = jags_res$jagsRes$q2.5[["Tau"]]
    tau_sum[["q97.5"]] = jags_res$jagsRes$q97.5[["Tau"]]
    tau_sum[["rhat"]] = jags_res$jagsRes$Rhat[["Tau"]]

    tau_list <- jags_res$jagsRes$sims.list$Tau
    sigrhoests <- getCorTauList.multspec(tau_list)

    # get the occurrence parameters
    names = c("mean","median","q2.5","q97.5","rhat")
    varrownames <- c("beta0")
    for ( j in 1:noccvar) { varrownames=c(varrownames,paste0("beta",j)) }
    beta_sum <- list()
    i=1;j=1
    for ( i in 1:nspec ) {
      beta_sum[[i]] <- array(NA,dim=c(noccvar+1,5))
      colnames(beta_sum[[i]]) <- names
      beta_sum[[i]][1,] = cbind(jags_res$jagsRes$mean[["lpsi"]][i],jags_res$jagsRes$q50[["lpsi"]][i],jags_res$jagsRes$q2.5[["lpsi"]][i],jags_res$jagsRes$q97.5[["lpsi"]][i],jags_res$jagsRes$Rhat[["lpsi"]][i])
      for ( j in 1:noccvar) {
        beta_sum[[i]][j+1,] = cbind(jags_res$jagsRes$mean[["betalpsi"]][i,j],jags_res$jagsRes$q50[["betalpsi"]][i,j],jags_res$jagsRes$q2.5[["betalpsi"]][i,j],jags_res$jagsRes$q97.5[["betalpsi"]][i,j],jags_res$jagsRes$Rhat[["betalpsi"]][i,j])
      }
      rownames(beta_sum[[i]]) <- varrownames
    }

    # hyper occurrence parameters
    dnames=list(names,varrownames,c("mu","sd"))
    dname = list("Covars","Estimators","HyperParams")
    hyper_beta_sum <- array(NA,dim=c(noccvar+1,5,2) )
    hyper_beta_sum[1,,1] = cbind(jags_res$jagsRes$mean[["mu.lpsi"]],jags_res$jagsRes$q50[["mu.lpsi"]],jags_res$jagsRes$q2.5[["mu.lpsi"]],jags_res$jagsRes$q97.5[["mu.lpsi"]],jags_res$jagsRes$Rhat[["mu.lpsi"]])
    hyper_beta_sum[1,,2] = cbind(jags_res$jagsRes$mean[["sd.lpsi"]],jags_res$jagsRes$q50[["sd.lpsi"]],jags_res$jagsRes$q2.5[["sd.lpsi"]],jags_res$jagsRes$q97.5[["sd.lpsi"]],jags_res$jagsRes$Rhat[["sd.lpsi"]])
    for ( j in 1:noccvar) {
      # print(j)
      hyper_beta_sum[j+1,,1] = cbind(jags_res$jagsRes$mean[["mu.betalpsi"]][j],jags_res$jagsRes$q50[["mu.betalpsi"]][j],jags_res$jagsRes$q2.5[["mu.betalpsi"]][j],jags_res$jagsRes$q97.5[["mu.betalpsi"]][j],jags_res$jagsRes$Rhat[["mu.betalpsi"]][j])
      hyper_beta_sum[j+1,,2] = cbind(jags_res$jagsRes$mean[["sd.betalpsi"]][j],jags_res$jagsRes$q50[["sd.betalpsi"]][j],jags_res$jagsRes$q2.5[["sd.betalpsi"]][j],jags_res$jagsRes$q97.5[["sd.betalpsi"]][j],jags_res$jagsRes$Rhat[["sd.betalpsi"]][j])
    }

    # get detection variables
    alpha_sum <- NULL # return val if no detection
    hyper_alpha_sum <- NULL
    if( type == 1) {  # if type 1 then there is detection
      # detection parameters
      varrownames <- c("alpha0")
      for ( j in 1:ndetvar) { varrownames=c(varrownames,paste0("alpha",j)) }
      alpha_sum <- list()
      j=1;i=1
      for ( i in 1:nspec ) {
        alpha_sum[[i]] <- array(NA,dim=c(ndetvar+1,5))
        colnames(alpha_sum[[i]]) <- names
        alpha_sum[[i]][1,] = cbind(jags_res$jagsRes$mean[["lp"]][i],jags_res$jagsRes$q50[["lp"]][i],jags_res$jagsRes$q2.5[["lp"]][i],jags_res$jagsRes$q97.5[["lp"]][i],jags_res$jagsRes$Rhat[["lp"]][i])
        for ( j in 1:ndetvar) {
          alpha_sum[[i]][j+1,] = cbind(jags_res$jagsRes$mean[["alphalp"]][i,j],jags_res$jagsRes$q50[["alphalp"]][i,j],jags_res$jagsRes$q2.5[["alphalp"]][i,j],jags_res$jagsRes$q97.5[["alphalp"]][i,j],jags_res$jagsRes$Rhat[["alphalp"]][i,j])
        }
        rownames(alpha_sum[[i]]) <- varrownames
      }

      # hyper detection parameters
      hyper_alpha_sum <- array(NA,dim=c(ndetvar+1,5,2) )
      hyper_alpha_sum[1,,1] = cbind(jags_res$jagsRes$mean[["mu.lp"]],jags_res$jagsRes$q50[["mu.lp"]],jags_res$jagsRes$q2.5[["mu.lp"]],jags_res$jagsRes$q97.5[["mu.lp"]],jags_res$jagsRes$Rhat[["mu.lp"]])
      hyper_alpha_sum[1,,2] = cbind(jags_res$jagsRes$mean[["sd.lp"]],jags_res$jagsRes$q50[["sd.lp"]],jags_res$jagsRes$q2.5[["sd.lp"]],jags_res$jagsRes$q97.5[["sd.lp"]],jags_res$jagsRes$Rhat[["sd.lp"]])
      for ( j in 1:ndetvar) {
        hyper_alpha_sum[j+1,,1] = cbind(jags_res$jagsRes$mean[["mu.alphalp"]][j],jags_res$jagsRes$q50[["mu.alphalp"]][j],jags_res$jagsRes$q2.5[["mu.alphalp"]][j],jags_res$jagsRes$q97.5[["mu.alphalp"]][j],jags_res$jagsRes$Rhat[["mu.alphalp"]][j])
        hyper_alpha_sum[j+1,,2] = cbind(jags_res$jagsRes$mean[["sd.alphalp"]][j],jags_res$jagsRes$q50[["sd.alphalp"]][j],jags_res$jagsRes$q2.5[["sd.alphalp"]][j],jags_res$jagsRes$q97.5[["sd.alphalp"]][j],jags_res$jagsRes$Rhat[["sd.alphalp"]][j])
      }
    }# end if type==1 (i.e. det model)
    param_sums[[res]]=list(resid=resid,tau_sum=tau_sum,
                           var.est=sigrhoests$var.est,covar.est=sigrhoests$covar.est,rho.est=sigrhoests$rho.est,
                           beta_sum=beta_sum,hyper_beta_sum=hyper_beta_sum,
                           alpha_sum=alpha_sum,hyper_alpha_sum=hyper_alpha_sum,
                           sig.list=sigrhoests$estsigmalist,rho.list=sigrhoests$estcorlist
    )
  } # end for res in reslen
  return(param_sums)

}


###########################################
#
# printouts_rl2
#
# summary data raw parameters (tau, beta, gamma, hyperparameter) and training sets for printing
# data in this format used to calculate ROC AUC in getOccDecProb
#
#' @section `printouts_rl2()`
#' Prints out raw parameter estimates for each entry in a list of jags result files.
#' @return A list of estimates for each file including
#' - DIC
#' - Beta parameters
#' - Alpha parameters
#' - Tau parameters (inverse of covariance matrix, Tau is estiamted by JAGS)
#'
#' @export printouts_rl2

printouts_rl2 <- function (list_jagfiles,windat.list,n.train,wnd.idx=1,savedir="") {
  reslen <- length(list_jagfiles)
  param_sums <- list()
  for ( res in 1:reslen) {
    reSave=FALSE
    file_base = list_jagfiles[[res]][[1]]
    jags_res <- readRDS(paste0(savedir,file_base) )
    resid <- jags_res$resId
    nspec <- resid$nspec
    type <- resid$type
    ntrain <- resid$ntrain
    windat <- windat.list[[1]][[wnd.idx]][[Windat.type[type]]]
    if (is.null(nspec)) {nspec <- windat$nspec ;jags_res$resId$nspec=nspec ; reSave=TRUE }
    if (is.null(ntrain) || identical(ntrain,numeric(0)) ) {
      if ( !is.null(n.train) ) {ntrain=n.train} else {ntrain=windat$nsite}
      jags_res$resId$ntrain=ntrain ;
      jags_res$resId$nvalid = windat$nsite-ntrain
      reSave=TRUE
    }
    noccvar <- windat$novar.occ
    ndetvar <- ifelse( is.null(windat$novar.det), 0, windat$novar.det )
    resid$noccvar <- noccvar
    resid$ndetvar <- ndetvar
    if (reSave==TRUE ) {  # update the jags file with the training/validation/nspec data
      saveRDS(jags_res,paste0(savedir,file_base) )
    }
    DIC=jags_res$jagsRes$DIC
    # get precision
    tau_sum <- list()
    tau_sum[["mean"]] = jags_res$jagsRes$mean[["Tau"]]
    tau_sum[["median"]] = jags_res$jagsRes$q50[["Tau"]]
    tau_sum[["q2.5"]] = jags_res$jagsRes$q2.5[["Tau"]]
    tau_sum[["q97.5"]] = jags_res$jagsRes$q97.5[["Tau"]]
    tau_sum[["rhat"]] = jags_res$jagsRes$Rhat[["Tau"]]

    tau_list <- jags_res$jagsRes$sims.list$Tau
    sigrhoests <- getCorTauList.multspec(tau_list)

    names = c("mean","median","q2.5","q97.5","rhat")
    varrownames <- c("beta0")
    for ( j in 1:noccvar) { varrownames=c(varrownames,paste0("beta",j)) }
    beta_sum <- list()
    i=1;j=1
    for ( i in 1:nspec ) {
      beta_sum[[i]] <- array(NA,dim=c(noccvar+1,5))
      colnames(beta_sum[[i]]) <- names
      beta_sum[[i]][1,] = cbind(jags_res$jagsRes$mean[["lpsi"]][i],jags_res$jagsRes$q50[["lpsi"]][i],jags_res$jagsRes$q2.5[["lpsi"]][i],jags_res$jagsRes$q97.5[["lpsi"]][i],jags_res$jagsRes$Rhat[["lpsi"]][i])
      for ( j in 1:noccvar) {
        beta_sum[[i]][j+1,] = cbind(jags_res$jagsRes$mean[["betalpsi"]][i,j],jags_res$jagsRes$q50[["betalpsi"]][i,j],jags_res$jagsRes$q2.5[["betalpsi"]][i,j],jags_res$jagsRes$q97.5[["betalpsi"]][i,j],jags_res$jagsRes$Rhat[["betalpsi"]][i,j])
      }
      rownames(beta_sum[[i]]) <- varrownames
    }

    # hyper occurrence parameters
    dnames=list(names,varrownames,c("mu","sd"))
    dname = list("Covars","Estimators","HyperParams")
    hyper_beta_sum <- array(NA,dim=c(noccvar+1,5,2) )
    hyper_beta_sum[1,,1] = cbind(jags_res$jagsRes$mean[["mu.lpsi"]],jags_res$jagsRes$q50[["mu.lpsi"]],jags_res$jagsRes$q2.5[["mu.lpsi"]],jags_res$jagsRes$q97.5[["mu.lpsi"]],jags_res$jagsRes$Rhat[["mu.lpsi"]])
    hyper_beta_sum[1,,2] = cbind(jags_res$jagsRes$mean[["sd.lpsi"]],jags_res$jagsRes$q50[["sd.lpsi"]],jags_res$jagsRes$q2.5[["sd.lpsi"]],jags_res$jagsRes$q97.5[["sd.lpsi"]],jags_res$jagsRes$Rhat[["sd.lpsi"]])
    for ( j in 1:noccvar) {
      hyper_beta_sum[j+1,,1] = cbind(jags_res$jagsRes$mean[["mu.betalpsi"]][j],jags_res$jagsRes$q50[["mu.betalpsi"]][j],jags_res$jagsRes$q2.5[["mu.betalpsi"]][j],jags_res$jagsRes$q97.5[["mu.betalpsi"]][j],jags_res$jagsRes$Rhat[["mu.betalpsi"]][j])
      hyper_beta_sum[j+1,,2] = cbind(jags_res$jagsRes$mean[["sd.betalpsi"]][j],jags_res$jagsRes$q50[["sd.betalpsi"]][j],jags_res$jagsRes$q2.5[["sd.betalpsi"]][j],jags_res$jagsRes$q97.5[["sd.betalpsi"]][j],jags_res$jagsRes$Rhat[["sd.betalpsi"]][j])
    }

    # get detection variables
    alpha_sum <- NULL # return val if no detection
    hyper_alpha_sum <- NULL
    if( type == 1) {  # if type 1 then there is detection
      # detection parameters
      varrownames <- c("alpha0")
      for ( j in 1:ndetvar) { varrownames=c(varrownames,paste0("alpha",j)) }
      alpha_sum <- list()
      j=1;i=1
      for ( i in 1:nspec ) {
        alpha_sum[[i]] <- array(NA,dim=c(ndetvar+1,5))
        colnames(alpha_sum[[i]]) <- names
        alpha_sum[[i]][1,] = cbind(jags_res$jagsRes$mean[["lp"]][i],jags_res$jagsRes$q50[["lp"]][i],jags_res$jagsRes$q2.5[["lp"]][i],jags_res$jagsRes$q97.5[["lp"]][i],jags_res$jagsRes$Rhat[["lp"]][i])
        for ( j in 1:ndetvar) {
          alpha_sum[[i]][j+1,] = cbind(jags_res$jagsRes$mean[["alphalp"]][i,j],jags_res$jagsRes$q50[["alphalp"]][i,j],jags_res$jagsRes$q2.5[["alphalp"]][i,j],jags_res$jagsRes$q97.5[["alphalp"]][i,j],jags_res$jagsRes$Rhat[["alphalp"]][i,j])
        }
        rownames(alpha_sum[[i]]) <- varrownames
      }

      # hyper detection parameters
      hyper_alpha_sum <- array(NA,dim=c(ndetvar+1,5,2) )
      hyper_alpha_sum[1,,1] = cbind(jags_res$jagsRes$mean[["mu.lp"]],jags_res$jagsRes$q50[["mu.lp"]],jags_res$jagsRes$q2.5[["mu.lp"]],jags_res$jagsRes$q97.5[["mu.lp"]],jags_res$jagsRes$Rhat[["mu.lp"]])
      hyper_alpha_sum[1,,2] = cbind(jags_res$jagsRes$mean[["sd.lp"]],jags_res$jagsRes$q50[["sd.lp"]],jags_res$jagsRes$q2.5[["sd.lp"]],jags_res$jagsRes$q97.5[["sd.lp"]],jags_res$jagsRes$Rhat[["sd.lp"]])
      for ( j in 1:ndetvar) {
        hyper_alpha_sum[j+1,,1] = cbind(jags_res$jagsRes$mean[["mu.alphalp"]][j],jags_res$jagsRes$q50[["mu.alphalp"]][j],jags_res$jagsRes$q2.5[["mu.alphalp"]][j],jags_res$jagsRes$q97.5[["mu.alphalp"]][j],jags_res$jagsRes$Rhat[["mu.alphalp"]][j])
        hyper_alpha_sum[j+1,,2] = cbind(jags_res$jagsRes$mean[["sd.alphalp"]][j],jags_res$jagsRes$q50[["sd.alphalp"]][j],jags_res$jagsRes$q2.5[["sd.alphalp"]][j],jags_res$jagsRes$q97.5[["sd.alphalp"]][j],jags_res$jagsRes$Rhat[["sd.alphalp"]][j])
      }
    }# end if type==1 (i.e. det model)
    param_sums[[res]]=list(resid=resid,DIC=DIC,tau_sum=tau_sum,
                           var.est=sigrhoests$var.est,covar.est=sigrhoests$covar.est,rho.est=sigrhoests$rho.est,
                           beta_sum=beta_sum,hyper_beta_sum=hyper_beta_sum,
                           alpha_sum=alpha_sum,hyper_alpha_sum=hyper_alpha_sum,
                           sig.list=sigrhoests$estsigmalist,rho.list=sigrhoests$estcorlist
    )

  } # end for res in reslen
  return(param_sums)
}



###########################################
#
# getOccDecProb
#
# calc the ROC for each of a set of jagfiles
#
#' @section `getOccDecProb()` Determines AUC ROC data for list of jag files
#' Plots ROC data in series of plots
#' Prints ROC outputs
#' @return A list of ROC results data for each jag result file in the input list
#' For each jagsfile, the following data is included (in list form)
#' - covariance and correlation matrix estimates
#' - ROC data for the validation dataset (average and separate for each survey replication for each species)
#' - ROC data for the training dataset (average and separate for each survey replication for each species)
#'
#' @export getOccDecProb


getOccDecProb <- function (list_jagfiles,list_prresp,Valid_set,savedir="",plot_muocc_flag=FALSE,view_pred_flag=FALSE) {
  reslen <- length(list_jagfiles)
  param_sums <- list()
  for ( res in 1:reslen) {
    param_sums[[res]] <- list()
    file_base = list_jagfiles[[res]][[1]]
    jags_res <- readRDS(paste0(savedir,file_base) )
    resid <- jags_res$resId
    nspec <- resid$nspec
    type <- resid$type
    nsite <- resid$nsite
    ntrain <- resid$ntrain
    nvalid <- resid$nvalid
    if (is.null(nsite)) {nsite=ntrain+nvalid}
    if ( plot_muocc_flag ) {
      #mean.mu.occ <- jags_res$jagsRes$mean[["mu.occ"]]
      jset <- sample.int(ntrain, size=6)
      par(mfrow=c(3,nspec))
      j=jset[1];i=1;k=2
      for ( j in jset ) {
        for (i in 1:nspec ) {
          title1=paste("Histogram mu occupancy training, species",i,"site",j,sep=" ")
          lab1 = paste("Mu species",i,"site",j,sep=" ")
          hist(jags_res$jagsRes$sims.list$mu.occ[,j,i],main=list(title1,cex=0.9),xlab=list(lab1,cex=0.8) )
          dim(jags_res$jagsRes$mean[["mu.occ"]])
          abline(v=jags_res$jagsRes$mean[["mu.occ"]][j,i],col="red",lty="solid")
          abline(v=jags_res$jagsRes$q50[["mu.occ"]][j,i],col="blue",lty="dashed")
          abline(v=jags_res$jagsRes$q2.5[["mu.occ"]][j,i],col="grey",lty="dotted")
          abline(v=jags_res$jagsRes$q97.5[["mu.occ"]][j,i],col="grey",lty="dotted")
        }# end for
      }# end for
      jset <- sample((ntrain+1):nsite, size=6)
      par(mfrow=c(3,nspec))
      for ( j in jset ) {
        for (i in 1:nspec ) {
          title1=paste("Histogram mu occupancy validation, species",i,"site",j,sep=" ")
          lab1 = paste("Mu species",i,"site",j,sep=" ")
          hist(jags_res$jagsRes$sims.list$mu.occ[,j,i],main=list(title1,cex=0.9),xlab=list(lab1,cex=0.8) )
          abline(v=jags_res$jagsRes$mean[["mu.occ"]][j,i],col="red",lty="solid")
          abline(v=jags_res$jagsRes$q50[["mu.occ"]][j,i],col="blue",lty="dashed")
          abline(v=jags_res$jagsRes$q2.5[["mu.occ"]][j,i],col="grey",lty="dotted")
          abline(v=jags_res$jagsRes$q97.5[["mu.occ"]][j,i],col="grey",lty="dotted")
        }# end for
      }# end for
    }# end plot mu_occ

    # get p occ values
    # get a correlation or covariance mx
    corr_mx <- array(dim=c(nspec,nspec))
    cov_mx <- array(dim=c(nspec,nspec))
    for ( i in 1:nspec ) {
      corr_mx[i,i] = 1
      cov_mx[i,i] = list_prresp[[res]]$var.est[i,1]
      if ( i < nspec ) {
        for ( k in (i+1):nspec ){
          corr_mx[i,k] <- list_prresp[[res]]$rho.est[i,k,1]
          corr_mx[k,i] <- list_prresp[[res]]$rho.est[i,k,1]
          cov_mx[i,k] <- list_prresp[[res]]$covar.est[i,k,1]
          cov_mx[k,i] <- list_prresp[[res]]$covar.est[i,k,1]
        }
      }
    }

    # we have 4 species, want marginal prob for each
    #species 1
    #lowerbnd = rbind(c(0,-Inf,-Inf,-Inf),c(-Inf,0,-Inf,-Inf),c(-Inf,-Inf,0,-Inf),c(-Inf,-Inf,-Inf,0))
    lowerbnd=array(dim=c(nspec,nspec))
    for ( i in 1:nspec) {
      lowerbnd[i,i]=0
      for ( k in i+1:nspec) {
        if ( k <= nspec) {lowerbnd[i,k]=lowerbnd[k,i]=-Inf}
      }
    }
    upperbnd = rep(Inf,nspec) #c(Inf,Inf,Inf,Inf)
    psi_occ_mean <- array(dim=c(nsite,nspec))
    #spec=1;site=1 ; site=2
    for (site in 1:nsite ) {
      for ( spec in 1:nspec ) {
        mean_site = jags_res$jagsRes$mean[["mu.occ"]][site,]
        prob_occ_mean <- mvtnorm::pmvnorm(lowerbnd[spec,],upperbnd,mean_site,sigma=cov_mx)
        psi_occ_mean[site,spec] <- prob_occ_mean[[1]]
      }
    }
    param_sums[[res]]$cor_mx.mn <- corr_mx
    param_sums[[res]]$cov_mx.mn <- cov_mx
    param_sums[[res]]$psi_occ_mean.mn <- psi_occ_mean

    # get predictions
    # get estimated V (latent V) values
    V_mean <- jags_res$jagsRes$mean[["V"]]
    # get Y values: these are the given Y values for 1:ntrain
    Y_mean <- jags_res$jagsRes$mean[["y"]]
    if (type == 1 ) {
      # get the estimated Z values
      Z_mean <- jags_res$jagsRes$mean[["z"]]
      # get the estimated detection probability values
      mean.p.det <- jags_res$jagsRes$mean[["p"]]
      #dim(mean.p.det)
      # get the probability of observed occupancy
      predictor_mean <- array(dim=dim(mean.p.det))
      predictor_mean.av <- array(dim=c(dim(mean.p.det)[1],dim(mean.p.det)[3]))
      for ( i in 1:dim(mean.p.det)[2]) {
        predictor_mean[,i,] <- mean.p.det[,i,]*psi_occ_mean
      }
      param_sums[[res]]$prob_occ_det.mn <- predictor_mean
      param_sums[[res]]$actual_tr = Y_mean[1:ntrain,,]
      param_sums[[res]]$actual_va=  Valid_set$y

      # average predictions: get average observed occupancy
      predictor_mean.av <- apply(predictor_mean,c(1,3),mean,na.rm=TRUE)
      # get the collapsed y actual values for training set (mean_y)
      actual.clp.tr <- apply(Y_mean[1:ntrain,,],c(1,3),max,na.rm=TRUE)
      # get the collapsed y actual values for validation set (windat$other$valid)
      actual.clp.va <- Valid_set$y_clp
      param_sums[[res]]$prob_occ_det.mn.av <- predictor_mean.av
      param_sums[[res]]$actual_tr.av <- actual.clp.tr
      param_sums[[res]]$actual_va.av <- actual.clp.va

      # long predictions: rbind the all the observations sets
      # for both training and validation
      # get long probabilities of observed occupancy for training set
      predictor_mean.long.tr = rbind(predictor_mean[1:ntrain,1,],predictor_mean[1:ntrain,2,],predictor_mean[1:ntrain,3,])
      # get the long  y actual values for training set (mean_y)
      actual.long.tr = rbind(Y_mean[1:ntrain,1,],Y_mean[1:ntrain,2,],Y_mean[1:ntrain,3,])
      # get long probabilities of observed occupancy for validation set
      predictor_mean.long.va <- rbind(predictor_mean[(ntrain+1):nsite,1,], predictor_mean[(ntrain+1):nsite,2,],predictor_mean[(ntrain+1):nsite,3,])
      # get the long  y actual values for validation set  (windat$other$valid)
      actual.long.va <- rbind(Valid_set$y[,1,], Valid_set$y[,2,],Valid_set$y[,3,])
      tr.complete <- complete.cases(actual.long.tr)
      va.complete <- complete.cases(actual.long.va)
      #dim(predictor_mean.long.tr)
      predictor_mean.long.tr <- predictor_mean.long.tr[tr.complete,]
      actual.long.tr <- actual.long.tr[tr.complete,]
      predictor_mean.long.va <- predictor_mean.long.va[va.complete,]
      actual.long.va <- actual.long.va[va.complete,]
      param_sums[[res]]$prob_occ_det.mn.lg.tr <- predictor_mean.long.tr
      param_sums[[res]]$prob_occ_det.mn.lg.va <- predictor_mean.long.va
      param_sums[[res]]$actual_tr.lg <- actual.long.tr
      param_sums[[res]]$actual_va.lg <- actual.long.va
      if ( view_pred_flag ) {
        View(predictor_mean.long.tr)
        View(predictor_mean.long.va)
      }


      pred.tr.list <- list()
      perf.tr.list <- list()
      pred.va.list <- list()
      perf.va.list <- list()
      par(mfrow=c(4,nspec))
      # i=1
      for ( i in 1:nspec ) {
        tr.complete.nreps = complete.cases(actual.spec.tr)
        va.complete.nreps = complete.cases(actual.spec.va)
        predictor_mean.spec.tr = predictor_mean[1:ntrain,1:2,i]
        actual.spec.tr = Y_mean[1:ntrain,1:2,i]
        predictor_mean.spec.va = predictor_mean[(ntrain+1):nsite,1:2,i]
        actual.spec.va= Valid_set$y[,1:2,i]
        dim(predictor_mean);dim(Y_mean[1:ntrain,,]);dim(Valid_set$y)
        pred.tr.list[[i]] <- ROCR::prediction(predictor_mean.spec.tr[tr.complete.nreps,],actual.spec.tr[tr.complete.nreps,] )
        perf.tr.list[[i]] <- list()
        perf.tr.list[[i]]$roc <- ROCR::performance(pred.tr.list[[i]],"tpr","fpr")
        perf.tr.list[[i]]$auc <- ROCR::performance(pred.tr.list[[i]],"auc")

        pred.va.list[[i]] <- ROCR::prediction(predictor_mean.spec.va[va.complete.nreps,],actual.spec.va[va.complete.nreps,] )
        perf.va.list[[i]] <- list()
        perf.va.list[[i]]$roc <- ROCR::performance(pred.va.list[[i]],"tpr","fpr")
        perf.va.list[[i]]$auc <- ROCR::performance(pred.va.list[[i]],"auc")

        title1.tr <- paste("ROC species",i,"training, 1 per obs",sep=" ")
        title2.tr.av <- paste("ROC species",i,"training average",sep=" ")
        plot(perf.tr.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.tr,cex=0.9))
        AUC1 <- as.numeric(perf.tr.list[[i]]$auc@y.values[1])
        AUC2 <- as.numeric(perf.tr.list[[i]]$auc@y.values[2])
        text(0.6,0.2,paste("AUC tr obs1:",round(AUC1,4),"obs2:",round(AUC2,4) ) )
        plot(perf.tr.list[[i]]$roc,colorize=TRUE,avg="threshold",main=list(title2.tr.av,cex=0.9))
        text(0.6,0.2,paste("AUC tr mean obs1 and 2:",round(mean(AUC1,AUC2),4) ) )
        title1.va<- paste("ROC species",i,"validation, 1 per obs",sep=" ")
        title2.va.av <- paste("ROC species",i,"validation average",sep=" ")
        plot(perf.va.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.va,cex=0.9))
        AUC1 <- as.numeric(perf.va.list[[i]]$auc@y.values[1])
        AUC2 <- as.numeric(perf.va.list[[i]]$auc@y.values[2])
        text(0.6,0.1,paste("AUC val obs1:",round(AUC1,4),"obs2:",round(AUC2,4) ) )
        plot(perf.va.list[[i]]$roc,colorize=TRUE,avg="threshold",main=list(title2.va.av,cex=0.9))
        text(0.6,0.1,paste("AUC val mean obs1 and 2:",round(mean(AUC1,AUC2),4) ) )
      }
      param_sums[[res]]$pred.tr.list=  pred.tr.list
      param_sums[[res]]$perf.tr.list=  perf.tr.list
      param_sums[[res]]$pred.va.list=  pred.va.list
      param_sums[[res]]$perf.va.list=  perf.va.list

      for ( i in 1:nspec ) {
        AUC.tr.ob1 <- as.numeric(perf.tr.list[[i]]$auc@y.values[1])
        AUC.tr.ob2 <- as.numeric(perf.tr.list[[i]]$auc@y.values[2])
        AUC.va.ob1 <- as.numeric(perf.va.list[[i]]$auc@y.values[1])
        AUC.va.ob2 <- as.numeric(perf.va.list[[i]]$auc@y.values[2])
        # just using ob1 and 2 due to small number of 3rd values
        # which could not always be computed.
        print(paste("AUC (ROC) for species",i,"obs1","training",round(AUC.tr.ob1,4),"validation",round(AUC.va.ob1,4),sep=" "))
        print(paste("AUC (ROC) for species",i,"obs2","training",round(AUC.tr.ob2,4),"validation",round(AUC.va.ob2,4),sep=" "))
        if( length(perf.tr.list[[i]]$auc@y.values)>=3 && length(perf.va.list[[i]]$auc@y.values)>=3 ) {
          AUC.tr.ob3 <- as.numeric(perf.tr.list[[i]]$auc@y.values[3])
          AUC.va.ob3 <- as.numeric(perf.va.list[[i]]$auc@y.values[3])
          print(paste("AUC (ROC) for species",i,"obs3","training",round(AUC.tr.ob3,4),"validation",round(AUC.va.ob3,4),sep=" "))

        }
      }

      # plots and reports for performance data (av and combined)
      # performance and prediction lists for long and average
      # for species separately
      pred.tr.av.list <- list()
      perf.tr.av.list <- list()
      pred.va.av.list <- list()
      perf.va.av.list <- list()
      pred.tr.long.list <- list()
      perf.tr.long.list <- list()
      pred.va.long.list <- list()
      perf.va.long.list <- list()
      # i=1
      par(mfrow=c(4,nspec))
      # individual species ROC cacls for long and average lists
      for ( i in 1:nspec ) {
        # long: combine observations using rbind
        pred.tr.long.list[[i]] <- ROCR::prediction(predictor_mean.long.tr[,i],actual.long.tr[,i] )
        perf.tr.long.list[[i]] <- list()
        perf.tr.long.list[[i]]$roc <- ROCR::performance(pred.tr.long.list[[i]],"tpr","fpr")
        if ( nspec != 3 && nspec != 4 ) perf.tr.long.list[[i]]$cal<- ROCR::performance(pred.tr.long.list[[i]],"cal")
        perf.tr.long.list[[i]]$auc <- ROCR::performance(pred.tr.long.list[[i]],"auc")

        pred.va.av.list[[i]] <- ROCR::prediction(predictor_mean.long.va[,i],actual.long.va[,i] )
        perf.va.long.list[[i]] <- list()
        perf.va.long.list[[i]]$roc <- ROCR::performance(pred.va.av.list[[i]],"tpr","fpr")
        if ( nspec != 3 && nspec != 4 ) perf.va.long.list[[i]]$cal<- ROCR::performance(pred.va.av.list[[i]],"cal")
        perf.va.long.list[[i]]$auc <- ROCR::performance(pred.va.av.list[[i]],"auc")

        # av: collapsed y, and average of p.psi values
        pred.tr.av.list[[i]] <- ROCR::prediction(predictor_mean.av[1:ntrain,i],actual.clp.tr[,i] )
        perf.tr.av.list[[i]] <- list()
        perf.tr.av.list[[i]]$roc <- ROCR::performance(pred.tr.av.list[[i]],"tpr","fpr")
        if ( nspec != 3 && nspec != 4 ) perf.tr.av.list[[i]]$cal<- ROCR::performance(pred.tr.av.list[[i]],"cal")
        perf.tr.av.list[[i]]$auc <- ROCR::performance(pred.tr.av.list[[i]],"auc")

        pred.va.av.list[[i]] <- ROCR::prediction(predictor_mean.av[(ntrain+1):nsite,i],actual.clp.va[,i] )
        perf.va.av.list[[i]] <- list()
        perf.va.av.list[[i]]$roc <- ROCR::performance(pred.va.av.list[[i]],"tpr","fpr")
        if ( nspec != 3 && nspec != 4 ) perf.va.av.list[[i]]$cal<- ROCR::performance(pred.va.av.list[[i]],"cal")
        perf.va.av.list[[i]]$auc <- ROCR::performance(pred.va.av.list[[i]],"auc")

        # plots
        title1.tr.l <- paste("ROC species",i,"training, combined obs",sep=" ")
        title1.tr.a <- paste("ROC species",i,"training averaged obs",sep=" ")
        plot(perf.tr.long.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.tr.l,cex=0.9))
        AUC <- as.numeric(perf.tr.long.list[[i]]$auc@y.values)
        text(0.8,0.2,paste("AUC = ",round(AUC,4)))
        plot(perf.tr.av.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.tr.a,cex=0.9))
        AUC <- as.numeric(perf.tr.av.list[[i]]$auc@y.values)
        text(0.8,0.2,paste("AUC = ",round(AUC,4)))
        title1.va.l <- paste("ROC species",i,"validation, combined obs",sep=" ")
        title1.va.a <- paste("ROC species",i,"validation averaged obs",sep=" ")
        plot(perf.va.long.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.va.l,cex=0.9))
        AUC <- as.numeric(perf.va.long.list[[i]]$auc@y.values)
        text(0.8,0.2,paste("AUC = ",round(AUC,4)))
        plot(perf.va.av.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.va.a,cex=0.9))
        AUC <- as.numeric(perf.va.av.list[[i]]$auc@y.values)
        text(0.8,0.2,paste("AUC = ",round(AUC,4)))
      } # averaged/combined plots
      # export ROC_resp3_all_c_av_det_16panel

      param_sums[[res]]$pred.tr.long.list=  pred.tr.long.list
      param_sums[[res]]$perf.tr.long.list=  perf.tr.long.list
      param_sums[[res]]$pred.va.long.list=  pred.va.long.list
      param_sums[[res]]$perf.va.long.list=  perf.va.long.list
      param_sums[[res]]$pred.tr.av.list=  pred.tr.av.list
      param_sums[[res]]$perf.tr.av.list=  perf.tr.av.list
      param_sums[[res]]$pred.va.av.list=  pred.va.av.list
      param_sums[[res]]$perf.va.av.list=  perf.va.av.list

      # export: ROC_comb_av_allspec_det_16panel
      for ( i in 1:nspec ) {
        AUC.tr.lg <- as.numeric(perf.tr.long.list[[i]]$auc@y.values)
        AUC.tr.av <- as.numeric(perf.tr.av.list[[i]]$auc@y.values)
        AUC.va.lg <- as.numeric(perf.va.long.list[[i]]$auc@y.values)
        AUC.va.av <- as.numeric(perf.va.av.list[[i]]$auc@y.values)
        print(paste("AUC (ROC) for species",i,"combined","training",round(AUC.tr.lg,4),"validation",round(AUC.va.lg,4),sep=" "))
        print(paste("AUC (ROC) for species",i,"averaged","training",round(AUC.tr.av,4),"validation",round(AUC.va.av,4),sep=" "))
      }# averaged/combined printout

      # end type1 pred/perf
    } else { # else type != 1 -> no detection pred/perf
      predictor_mean <- psi_occ_mean
      predictor_mean.av <- psi_occ_mean
      actual.clp.tr <- Y_mean[1:ntrain,] #apply(Y_mean[1:ntrain,,],c(1,3),max,na.rm=TRUE)
      # get the collapsed y actual values for validation set (windat$other$valid)
      if( type == 2 ) {
        actual.clp.va <- Valid_set$y_clp
      } else {actual.clp.va <- Valid_set$y_rdm}
      param_sums[[res]]$prob_occ_det.mn <- predictor_mean
      param_sums[[res]]$prob_occ_det.mn.av <- predictor_mean.av
      param_sums[[res]]$actual_tr.av <- actual.clp.tr
      param_sums[[res]]$actual_va.av <- actual.clp.va

      pred.tr.av.list <- list()
      perf.tr.av.list <- list()
      pred.va.av.list <- list()
      perf.va.av.list <- list()
      i=1
      i=2
      for ( i in 1:nspec ) {
        # av: collapsed y, and average of p.psi values
        tr.complete.clp = complete.cases(actual.clp.tr[,i])
        va.complete.clp = complete.cases(actual.clp.va[,i])
        tr.predictor_mean.av.sp <- predictor_mean.av[1:ntrain,i]
        va.predictor_mean.av.sp <- predictor_mean.av[(ntrain+1):nsite,i]

        pred.tr.av.list[[i]] <- ROCR::prediction(tr.predictor_mean.av.sp[tr.complete.clp],actual.clp.tr[tr.complete.clp,i] )
        #pred.tr.av.list[[i]] <- prediction(predictor_mean.av[1:ntrain,i],actual.clp.tr[,i] )
        perf.tr.av.list[[i]] <- list()
        perf.tr.av.list[[i]]$roc <- ROCR::performance(pred.tr.av.list[[i]],"tpr","fpr")
        perf.tr.av.list[[i]]$cal<- ROCR::performance(pred.tr.av.list[[i]],"cal")
        perf.tr.av.list[[i]]$auc <- ROCR::performance(pred.tr.av.list[[i]],"auc")

        pred.va.av.list[[i]] <- ROCR::prediction(va.predictor_mean.av.sp[va.complete.clp],actual.clp.va[va.complete.clp,i] )
        #pred.va.av.list[[i]] <- prediction(predictor_mean.av[(ntrain+1):nsite,i],actual.clp.va[,i] )
        perf.va.av.list[[i]] <- list()
        perf.va.av.list[[i]]$roc <- ROCR::performance(pred.va.av.list[[i]],"tpr","fpr")
        #perf.va.av.list[[i]]$cal<- performance(pred.va.av.list[[i]],"cal")
        perf.va.av.list[[i]]$auc <- ROCR::performance(pred.va.av.list[[i]],"auc")

        # plots
        title1.tr.a <- paste("ROC species",i,"training averaged obs",sep=" ")
        plot(perf.tr.av.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.tr.a,cex=0.9))
        AUC <- as.numeric(perf.tr.av.list[[i]]$auc@y.values)
        text(0.8,0.2,paste("AUC = ",round(AUC,4)))
        title1.va.a <- paste("ROC species",i,"validation averaged obs",sep=" ")
        plot(perf.va.av.list[[i]]$roc,colorize=TRUE,avg="none",main=list(title1.va.a,cex=0.9))
        AUC <- as.numeric(perf.va.av.list[[i]]$auc@y.values)
        text(0.8,0.2,paste("AUC = ",round(AUC,4)))
      }# averaged plots
      param_sums[[res]]$pred.tr.av.list=  pred.tr.av.list
      param_sums[[res]]$perf.tr.av.list=  perf.tr.av.list
      param_sums[[res]]$pred.va.av.list=  pred.va.av.list
      param_sums[[res]]$perf.va.av.list=  perf.va.av.list
      for ( i in 1:nspec ) {
        AUC.tr.av <- as.numeric(perf.tr.av.list[[i]]$auc@y.values)
        AUC.va.av <- as.numeric(perf.va.av.list[[i]]$auc@y.values)
        print(paste("AUC (ROC) for species",i,"averaged","training",round(AUC.tr.av,4),"validation",round(AUC.va.av,4),sep=" "))
      }# averaged printout
      # export ROC_resp3_all_c_av_clp_12panel
      # export ROC_resp3_all_c_av_fst_12panel
    }# end not type1 pred/perf

  }# end for res
  return(param_sums)
}

