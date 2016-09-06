# Project: ssimf

#'   Significance
#'
#' that is used to obtain the "percenta" line based on Wu and Huang (2004).
#'
#' @details For this program to work well, the minimum data length is 36
#'
#' @param percenta a parameter having a value between 0.0 ~ 1.0, e.g., 0.05 
#'                 represents 95#' confidence level (upper bound); and 0.95 
#'                 represents 5#' confidence level (lower bound)
#' @param imfs the true IMFs from running EMD code. The first IMF must
#'                 be included for it is used to obtain the relative mean
#'                 energy for other IMFs. The trend is not included.
#'
#' @return \item{sigline}{a two column matrix, with the first column the natural
#'                 logarithm of mean period, and the second column the
#'                 natural logarithm of mean energy for significance line}
#'          \item{logep}{a two colum matrix, with the first column the natural
#'                 logarithm of mean period, and the second column the
#'                 natural logarithm of mean energy for all IMFs}
#'
#' @references Wu, Z., and N. E Huang (2008). Ensemble Empirical Mode Decomposition: a noise-assisted data analysis method.
#' Advances in Adaptive Data Analysis. Vol.1, No.1. 1-41.  
#' 
#' @import RcppOctave
#' 
#' @export
#'

significance <- function(imfs, percenta){
  o_clear(all=TRUE)
  o_source("inst/octave.m")
  o_assign(warning('off', 'Octave:possible-matlab-short-circuit-operator'))
  o_assign(imfs=imfs)
  o_assign(percenta=percenta)
  signi <- o_eval("[sigline, logep] = significance(imfs, percenta)")
  return(signi)
}

