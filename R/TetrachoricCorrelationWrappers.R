###################################################################################
#
# Tetrachoric functions for MGP JSDM
#
# This file contains functions that prep data from either real or simulated sources
# for analysis with JAGS

#************************************************************************************
#
# Packages
#
# require(psych)
#

# Try tetrachoric functions for  MVP JSDM
#
#' Try tetrachoric functions for  MVP JSDM
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `tetrachoric.try()`,  `tetrachoric.try.rho()`, `tetrachoric.marg.try()` and
#' `tetrachoric.marg.try.rho` are wrapper functions for the psych tetrachoric
#'
#'
#' `tetrachoric.try()` calls the tetrachoric function using 4 joint probabilities
#'                   and captures any errors or warnings
#' `tetrachoric.try.rho()` uses tetrachoric.try with 4 joint probabilities
#'                   and returns just rho, the correlation value.
#' `tetrachoric.marg.try()` calls the tetrachoric function using 2 marginal probabilities
#'                   and one joint probability of co-ocurrence.
#'                   It captures any errors or warnings.
#' `tetrachoric.marg.try.rho()` uses tetrachoric.marg.try with 2 marginal probabilities
#'                   and one joint probability of co-ocurrence.
#'                   It returns just rho, the correlation value.
#'
#'
#'
#' @section this is a group of functions that provide wrappers for tetrachoric functions
#'   which are used by the MVP JSDM modelling and simulation functions.
#'
#' @param  x A vector of data acceptable to the tetrachoric, usually
#'          is a vector of 4 joint probabilities of occupancy which are
#'               *  probability both species co-occur,
#'               *  probability species 1 occurs only
#'               *  probability species 2 occurs only
#'               *  probability neither species occurs
#' @param m1  Marginal probability of occupancy for species 1
#' @param m2  Marginal probability of occupancy for species 2
#' @param j12  Joint probability of occupancy both species 1 and 2 occur
#'
#' @return
#'   * `tetrachoric.try()` and `tetrachoric.marg.try()` return psych tetrachoric strucutre
#'   * NA is returned if an error occurs
#'   * `tetrachoric.try.rho()` and `tetrachoric.marg.try.rho()` return rho - a single numeric value of correaltion
#'
#'


#************************************************************************************
#
# Try tetrachoric correlation
# function catches errors to prevent tetrachor crashing the simulation runs
#
#' @export
tetrachoric.try <- function(x) {
  res <- tryCatch(
    { # try componnent
      psych::tetrachoric(x)
    },
    warning = function(cond) {
      print("warning")
      return(psych::tetrachoric(x) )
      #return(cond)
    },
    error = function(cond) {
      print("error")
      return(NA)
    }
  )
  return(res)
}
#tetrachoric.try(x)
#tetrachoric(x)$rho

#' @export
tetrachoric.try.rho <- function(x) {
  res <- tetrachoric.try(x)
  rho <- ifelse(is.na(res) || is.null(res),NA,res$rho)
}


#' @export
tetrachoric.marg.try <- function(m1,m2,j12) {
  res <- tryCatch(
    { # try componnent
      psych::tetrachoric(c(m1,m2,j12))
    },
    warning = function(cond) {
      print("warning")
      return(tetrachoric(c(m1,m2,j12)) )
      #return(cond)
    },
    error = function(cond) {
      print("error")
      return(NA)
    }
  )
  return(res)
}

#' @export
tetrachoric.marg.try.rho <- function(m1,m2,j12) {
  res <- tetrachoric.try(m1,m2,j12)
  rho <- ifelse(is.na(res) || is.null(res),NA,res$rho)
}

