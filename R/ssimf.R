# Project: ssimf

#' Ensemble Empirical Mode Decomposition
#' 
#' This is an EMD/EEMD program
#'
#' @param Y Inputted data;1-d data only
#' @param Nstd ratio of the standard deviation of the added noise and that of Y;
#' @param NE Ensemble number for the EEMD
#' 
#' @return A matrix of N*(m+1) matrix, where N is the length of the input
#'         data Y, and m=fix(log2(N))-1. Column 1 is the original data, columns 2, 3, ...
#'         m are the IMFs from high to low frequency, and comlumn (m+1) is the
#'         residual (over all trend).
#'
#' @details It should be noted that when Nstd is set to zero and NE is set to 1, the
#' program degenerates to a EMD program.(for EMD Nstd=0,NE=1)
#' This code limited sift number=10 ,the stoppage criteria can't change. 
#'
#' @references Wu, Z., and N. E Huang (2008). Ensemble Empirical Mode Decomposition: a noise-assisted data analysis method.
#' Advances in Adaptive Data Analysis. Vol.1, No.1. 1-41.  
#' 
#' @import  RcppOctave
#' 
#' @export

eemd <- function(Y, Nstd, NE){
  o_clear(all=TRUE)
  o_source("octave.m")
  o_assign(Y=Y)
  o_assign(Nstd=Nstd)
  o_assign(NE=NE)
  ans <- o_eval("eemd(Y,Nstd,NE)")
  return(ans)
}


