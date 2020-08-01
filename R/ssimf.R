# Project: ssimf

#' Instantanesous Frequency
#'
#' this is a function to calculate instantaneous based on EMD method
#'
#' @details this is a function to calculate instantaneous based on EMD method--
#' normalize the absolute values and find maximum envelope for 5 times
#' then calculate the Quadrature ,Phase angle,then take difference to them
#' finally,the instantaneous frequency values of an IMF is found. 
#'
#' @param vimf an IMF
#' @param dt time interval of the imputted data
#' 
#' @return omega: instantanesous frequency, which is 2*PI/T, where T
#' is the period of an ascillation
#'
#' @references Wu, Z., & Huang, N. (2004). A study of the characteristics of white noise using the 
#' empirical mode decomposition method. Proceedings of the Royal Society of London. Series A: Mathematical, 
#' Physical and Engineering Sciences, 460(2046), 1597–1611. https://doi.org/10.1098/rspa.2003.1221 
#' 
#' @export
#'

ifndq <- function(vimf, dt){
  #
  #1.set initial parameters
  Nnormal <- 5#number of spline envelope normalization for AM
  rangetop <- 0.90#the threshold of outliner remove for instantaneous frequency values
  vlength <- max(matlab::size(vimf))
  vlength_1 <- vlength -1
  
  #2.find absolute values
  abs_vimf <- vector(length = length(vimf))
  for (i in 1:vlength){
    abs_vimf[i] <- vimf[i]
    if (abs_vimf[i] < 0){
      abs_vimf[i] <- -vimf[i]
    }
  }
  
  #3.find the spline envelope for AM,take those out-loop start
  for (jj in 1:Nnormal){
    ext <-  extrema(abs_vimf)
    spmax <- ext$spmax
    spmin <- ext$spmin
    dd <- 1:vlength
    upper <-  pracma::pchip(spmax[,1],spmax[,2],dd)
    
    #4.Normalize the envelope out
    for (i in 1:vlength){
      abs_vimf[i] <- abs_vimf[i]/upper[i]
    }
  }
  #3.find the spline envelope for AM,take those out-loop end
  
  #5.flip back those negative values after AM been removed
  nvimf <- vector()
  for (i in 1:vlength){
    nvimf[i] <- abs_vimf[i]
    if (vimf[i] < 0){
      nvimf[i] <- -abs_vimf[i]
    }
  }
  
  #6.Calculate the quadrature values
  dq <- vector()
  for (i in 1:vlength){
    dq[i] <- sqrt(1-(nvimf[i]*nvimf[i]))
  }
  
  #7.Calculate the differece of the phase angle and give +/- sign
  devi <- vector()
  for (i in 2:vlength_1){
    devi[i] <- nvimf[i+1]-nvimf[i-1]
    if (devi[i]>0 & nvimf[i]<1){
      dq[i] <- -dq[i]
    }
  }
  
  #8.create a algorithm to remove those outliner
  omgcos <- vector()
  rangebot <- -rangetop
  for (i in 2:(vlength-1)){
    if (nvimf[i]>rangebot & nvimf[i] < rangetop){
      #good original value,direct calculate instantaneous frequency
      omgcos[i] <- abs(nvimf[i+1]-nvimf[i-1])*0.5/sqrt(1-nvimf[i]*nvimf[i])
    } else {
      #bad original value,direct set -9999,mark them
      omgcos[i] <- -9999
    }
  }
  omgcos[1] <- -9999
  omgcos[vlength] <- -9999
  
  #9.remove those marked outliner
  ddd <- temp <- vector()
  jj <- 1
  for (i in 1:vlength){
    if (omgcos[i]>-1000){
      ddd[jj] <- i
      temp[jj] <- omgcos[i]
      jj <- jj+1
    }
  }
  
  #10.use cubic spline to smooth the instantaneous frequency values
  temp2 <- pracma::pchip(ddd,temp,dd)
  omgcos <- temp2
  
  #11.return the values back
  omega <- vector()
  for (i in 1:vlength){
    omega[i] <- omgcos[i]
  }
  pi2 <- pi*2
  omega <- omega/dt
  return(omega)
}


#'  Extrema
#'  
#' This is a utility program for cubic spline envelope,
#' the code is to  find out max values and max positions
#' min values and min positions 
#' 
#' @param in_data Inputted data, a time series to be sifted
#' 
#' @return 
#' \item{spmax}{The locations (col 1) of the maxima and its corresponding
#' values (col 2)}
#' \item{spmin}{The locations (col 1) of the minima and its corresponding 
#' values (col 2)}
#'
#' @details EMD uses Cubic Spline to be the Maximun and Minimum Envelope for
#' the data.Besides finding spline,end points should be noticed.
#'
#' @references Wu, Z., & Huang, N. (2004). A study of the characteristics of white noise using the 
#' empirical mode decomposition method. Proceedings of the Royal Society of London. Series A: Mathematical, 
#' Physical and Engineering Sciences, 460(2046), 1597–1611. https://doi.org/10.1098/rspa.2003.1221 
#'  
#' @export
#' 


extrema <- function(in_data){
  ext <-  Rlibeemd::extrema(in_data)
  spmax <- ext$maxima
  spmax[,1] <- spmax[,1] + 1
  spmin <- ext$minima
  spmin[,1] <- spmin[,1] + 1
  return(list(spmax=spmax, spmin=spmin))
}



#' Statistic
#'
#' that is used to obtain the "percenta" line based on Wu and Huang (2004).
#'
#' @details For this program to work well, the minimum data length is 36
#'
#' @param imfs the true IMFs from running EMD code. The first IMF must
#'                 be included for it is used to obtain the relative mean
#'                 energy for other IMFs. The trend is not included.
#'
#' @return \item{logep}{a two colum matrix, with the first column the natural
#'                 logarithm of mean period, and the second column the
#'                 natural logarithm of mean energy for all IMFs}
#'
#' @references Wu, Z., & Huang, N. (2004). A study of the characteristics of white noise using the 
#' empirical mode decomposition method. Proceedings of the Royal Society of London. Series A: Mathematical, 
#' Physical and Engineering Sciences, 460(2046), 1597–1611. https://doi.org/10.1098/rspa.2003.1221 
#' 
#' 
#' 
#' @export
#'

statistic <- function(imfs){
  
  nDof <- length(imfs[,1])
  columns <- length(imfs[1,])
  logep <- matrix(NA, nrow = columns, ncol = 2)
  for (i in 1:columns){
    logep[i,2] <- 0
    logep[i,1] <- 0
    for (j in 1:nDof){
      logep[i,2] <- logep[i,2]+imfs[j,i]*imfs[j,i]
    }
    logep[i,2] <- logep[i,2]/nDof
  }
  sfactor <- logep[1,2]
  for (i in 1:columns){
    logep[i,2] <- 0.5636*logep[i,2]/sfactor# 0.6441
  }
  for (i in 1:3){
    ext <-  Rlibeemd::extrema(imfs[,i])
    spmax <- ext$maxima
    spmin <- ext$minima
    temp <- length(spmax[,1])-1
    logep[i,1] <- nDof/temp
  }
  for (i in 4:columns){
    omega <- ifndq(imfs[,i],1)
    sumomega <- 0
    for (j in 1:nDof){
      sumomega <- sumomega+omega[j]
    }
    logep[i,1] <- nDof*2*pi/sumomega
  }
  logep <- 1.4427*log(logep)
  return(logep)
}

#' Critical Value
#'
#' that is used to obtain the "percenta" line based on Wu and
#' Huang (2004).
#'
#' @details For this program to work well, the minimum data length is 36
#'
#'  
#' @param percenta a parameter having a value between 0.0 ~ 1.0, e.g., 0.05 
#' represents 95#' confidence level (upper bound); and 0.95 
#' represents 5#' confidence level (lower bound) 
#' @param imfs the true IMFs from running EMD code. The first IMF must
#' be included. The trend is not included.
#'   
#' @return sigline: a two column matrix, with the first column the natural
#' logarithm of mean period, and the second column the
#' natural logarithm of mean energy for significance line
#' 
#' @references Wu, Z., & Huang, N. (2004). A study of the characteristics of white noise using the 
#' empirical mode decomposition method. Proceedings of the Royal Society of London. Series A: Mathematical, 
#' Physical and Engineering Sciences, 460(2046), 1597–1611. https://doi.org/10.1098/rspa.2003.1221 
#'
#' @references 
#' 
#' @export
#'

criticalvalue <- function(imfs, percenta){
  
  nDof <- length(imfs[,1])
  pdMax <- matlab::fix(log(nDof))+1
  
  pdIntv <- matlab::linspace(1,pdMax,100)
  yBar <- -pdIntv
  
  yUpper <- rep(0, 100)
  yLower <-   -3-pdIntv*pdIntv 
  
  
  sigline <- matrix(NA, nrow = 100, ncol = 2)
  for (i in 1:100){
    sigline[i,1] <- pdIntv[i]
    
    yPos <- matlab::linspace(yUpper[i],yLower[i],5000)
    dyPos <- yPos[1]-yPos[2]
    yPDF <- dist_value(yPos,yBar[i],nDof)
    sum <- sum(yPDF)
    
    jj1 <- 0
    jj2 <- 1
    psum1 <- 0.0
    psum2 <- yPDF[1]
    pratio1 <- psum1/sum
    pratio2 <- psum2/sum
    
    while (pratio2 < percenta){
      jj1 <- jj1+1
      jj2 <- jj2+1
      psum1 <- psum1+yPDF[jj1]
      psum2 <- psum2+yPDF[jj2]
      pratio1 <- psum1/sum
      pratio2 <- psum2/sum
      yref <- yPos[jj1]
    }
    sigline[i,2] <- yref + dyPos*(pratio2-percenta)/(pratio2-pratio1)
    sigline[i,2] <- sigline[i,2] + 0.066*pdIntv[i] + 0.12
  }
  sigline <- 1.4427*sigline
  return(sigline)
}



#' Distance Value
#' 
#' calculate the PDF value  for 'sigline' in function critical value
#' 
#' @param yPos An input array at which PDF values are calculated---y value
#' yPos is an 1D 5000 pt matrix for y value in interval - [yUpper,yLower]
#' @param yBar The mean value of y ---exp(yBar)=E-bar
#' @param nDof The number of degree of freedom---Ndof=Npt
#' 
#' @return PDF: a normalized output array -about the PDF distribution
#' 
#' @details This is a utility program being called by "confidence".
#' this code calculate PDF value under yPos range
#' within the mean value of y(y-bar) of a chi-square distribution
#' here main job is calculating equation(3.4)--PDF formula from the reference paper
#' 
#' @references Wu, Z., & Huang, N. (2004). A study of the characteristics of white noise using the 
#' empirical mode decomposition method. Proceedings of the Royal Society of London. Series A: Mathematical, 
#' Physical and Engineering Sciences, 460(2046), 1597–1611. https://doi.org/10.1098/rspa.2003.1221 
#' 
#' 
#' 
#' @export

dist_value <- function(yPos, yBar, nDof){
  
  ylen <- length(yPos)
  
  #1.start to form equation(3.4)--PDF formula
  eBar <- exp(yBar)#E-bar
  evalue <- exp(yPos)#E
  
  #2.calculate PDF value for every y value
  tmp3 <- vector(length = ylen)
  for (i in 1:ylen){
    tmp1 <- evalue[i]/eBar-yPos[i]#calculate---tmp1=(E/E-bar)-y
    tmp2 <- -tmp1*nDof*eBar/2#calculate--------tmp2=(-1/2)*E-bar*Ndof*tmp1
    tmp3[i] <- 0.5*nDof*eBar*log(nDof) + tmp2
  }
  
  #3.to ensure the converge of the calculation,divide by rscale
  #   because the PDF is a relative value,not a absolute value
  rscale <- max(tmp3)
  
  tmp4 <- tmp3 - rscale#minus means divide,after we take expontial
  PDF <-  exp(tmp4)
  return(PDF)
}


#' statistic significance test for imf
#' 
#' @param demd Time series object of class "mts" where 
#' series corresponds to IMFs of the input signal, with 
#' the last series being the final residual
#' @param alpha significance level, default 0.05 
#'   
#' @return A list with class attribute 'ssimf' holding the following elements:
#' \itemize{
#' \item{\code{signal}}{The input signal.}
#' \item{\code{noise}}{The extracted noise.}
#' \item{\code{denoise}}{The input signal without the extracted noise.}
#' \item{\code{sele}}{Non-significant IMFs.}
#' \item{\code{id}}{Index of all MFIs.}
#' \item{\code{stat}}{a two colum matrix, with the first column the natural
#'                 logarithm of mean period, and the second column the
#'                 natural logarithm of mean energy for all IMFs}
#' \item{\code{cv.lower}}{critical value (lower bound)}
#' \item{\code{cv.upper}}{critial value (upper bound)}
#' \item{\code{demd}}{IMFs of the input signal, with 
#' the last series being the final residual}
#' }
#' 
#' 
#' @importFrom utils tail
#' 
#' @export


significance <- function(demd, alpha = 0.05){
  imfs <- demd[, -dim(demd)[2]]
  id <- colnames(imfs)
  stat <- statistic(imfs)
  cv.lower <- criticalvalue(imfs, alpha/2)
  cv.upper <- criticalvalue(imfs, (1- (alpha/2)))
  # seleciona as imfs com base na estatistica 
  sele <- id[1]
  for(i in 2:length(id)){
    #imfs dentro do zona ruido
    logT <- tail(which(cv.upper[,1] <= stat[i,1]), 1)
    if(sum(stat[i,2] >= cv.upper[logT,2] & stat[i,2] <= cv.lower[logT,2])!=0){
      sele <- cbind(sele, id[i])
    }
  }
  signal <- apply(demd, 1, sum)
  noise <- apply(as.matrix(demd[,sele]), 1, sum)
  attributes(signal) <- attributes(noise) <- attributes(demd[,1])
  denoise <- signal - noise
  return(structure(list(signal=signal, noise=noise, denoise=denoise, 
                        sele=sele, id=id, stat=stat, cv.lower=cv.lower, 
                        cv.upper=cv.upper, demd=demd), 
                   class = "ssimf"))
}


#' Plot object class ssimf
#' 
#' @param x object with class attribute 'ssimf'
#' @param type character string giving the type of plot to be computed. 
#' Allowed values are "decompositon" (the default), "test" for signficance test
#'
#' @return plot significance graph. Shaded area isn't significance.
#' 
#' @import ggplot2 
#' 
#' @importFrom ggplot2 autoplot
#' 
#' @export


autoplot.ssimf <- function(x,  type = c("decomposition", "test")){
  type <- match.arg(type)
  # significance
  if(type == "test"){
    dados <- data.frame(id=x$id, logT=x$stat[,1], logE=x$stat[,2]) 
    dados2 <- data.frame(logT=x$cv.lower[,1], cv.lower=x$cv.lower[,2], cv.upper=x$cv.upper[,2]) 
    g <- ggplot2::ggplot()
    g <- g + ggplot2::geom_point(ggplot2::aes(x=logT, y=logE), data = dados)
    g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = cv.lower, ymax = cv.upper, x = logT), data = dados2, fill = "grey70", alpha = 0.5)
    g <- g + ggplot2::geom_text(ggplot2::aes(x=logT-0.3, y=logE+1.2, label=id), data = dados, size=4, hjust=0) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position="none")
  }
  if(type == "decomposition"){
    noise <- x$noise
    denoise <- x$denoise
    signal <- x$signal
    suppressWarnings(g <- tsplot(cbind(denoise, noise), cbind(signal,rep(NA, length(signal)))))
  }
  return(g)
}

