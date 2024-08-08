## Script to perform automated peak detection ## 

###############
## LIBRARIES ##
###############

library(ggplot2)
library(reshape2)
library(pracma)
library(Matrix)

###############
## FUNCTIONS ##
###############

#' Read arw file as data.frame
#' 
#' @param path character, path to file
#' @return data.frame 
#' 
#' @author Emery Bosten
#' @export
readFile = function(path){
  
  data = read.table(path, skip = 2, sep = "")
  
  return(data)  
}


#' Combine multiple samples into a single chromatogram
#' 
#' @param paths array of characters, paths to files
#' @param filenames array of characters, names of files
#' @return chrom chromatogam (data.frame)
#' 
#' @author Emery Bosten, smathieu
#' @export
combineChroms = function(path,filenames){
  ind_blank <- grepl("blank", tolower(filenames))
  if(any(ind_blank)){
    filenames <- c(filenames[!ind_blank],  filenames[ind_blank])
  }
  combinedChrom <- readFile(paste0(path, "/", filenames[1]))
  t = combinedChrom[,1]
  
  if(any(ind_blank)){
    #with blank file
    
    #if more than 1 file
    if (length(filenames)>1){
      for (i in 2:(length(filenames) - 1)) {
        chrom <- readFile(paste0(path, "/", filenames[i]))
        if (nrow(combinedChrom) == nrow(chrom)) {
          combinedChrom <- combinedChrom + chrom
        } else if (nrow(combinedChrom) > nrow(chrom)) {
          combinedChrom <- combinedChrom[1:nrow(chrom), ] + chrom
        } else {
          combinedChrom <- combinedChrom + chrom[1:nrow(combinedChrom), ]
        }
      }
    }
    #blank
    blank <- readFile(paste0(path, "/", filenames[length(filenames)]))
    
    if (nrow(combinedChrom) == nrow(blank)) {
      combinedChrom <- combinedChrom - (length(filenames) - 1) * blank
    } else if (nrow(combinedChrom) > nrow(blank)) {
      combinedChrom <- combinedChrom[1:nrow(blank), ] - (length(filenames) - 1) * blank
    } else {
      combinedChrom <- combinedChrom - (length(filenames) - 1) * blank[1:nrow(combinedChrom), ]
    }
  } else {
    #no blank file
    if(length(filenames) > 1){
      for (i in 2:(length(filenames))) {
        chrom <- readFile(paste0(path, "/", filenames[i]))
        
        if (nrow(combinedChrom) == nrow(chrom)) {
          combinedChrom <- combinedChrom + chrom
        } else if (nrow(combinedChrom) > nrow(chrom)) {
          combinedChrom <- combinedChrom[1:nrow(chrom), ] + chrom
        } else {
          combinedChrom <- combinedChrom + chrom[1:nrow(combinedChrom), ]
        }
      }
    }
  }
  
  
  chrom = data.frame(t,combinedChrom[,2])
  colnames(chrom) <- c("t","abs")
  
  
  return(chrom)
}


#' Perform Whittaker smoother (sparse implementation)
#' 
#' @param y noisy signal
#' @param lambda regularization parameter
#' @param d degree parameter
#' @return smooth signal
#' 
#' @author Emery Bosten
#' @export
WSsparse <- function (y, lambda, d) {
  E <- Diagonal(length(y))
  as.numeric(solve(E + (lambda * crossprod(diff(E, 1, d))), t(y)))
}


#' Compute L-Curve terms 
#' 
#' @param z smooth signal
#' @param y noisy signal
#' @param d degree parameter
#' @return penalty and error term
#' 
#' @author Emery Bosten
#' @export
Lcurve <- function(z, y, d) {
  
  # Compute the error term
  E <- sum((y - z)^2)
  
  # Compute the smoothness term
  D <- sum(diff(z, differences = d)^2)
  
  return(list(E = E, D = D))
}


#' Automatically denoise signal using WS-VCurve
#' 
#' @param y noisy signal
#' @param lambdas set of regularization parameters
#' @return z smooth signal
#' 
#' @author Emery Bosten
#' @export
VcurveWS <- function(y, lambdas) {
  
  d <- 2  # d : filter order parameter (d = 1, 2, or 3)
  
  Ps <- numeric(length(lambdas))  # store penalties for each lambda
  Rs <- numeric(length(lambdas))  # store residuals for each lambda
  xs <- matrix(0, nrow = length(y), ncol = length(lambdas))  # store estimates for each lambda
  
  # L-curve
  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]
    z <- WSsparse(t(y), lam, d)
    result <- Lcurve(z, t(y), d)
    Ps[i] <- result$D
    Rs[i] <- result$E
    xs[, i] <- z
  }
  
  # V-Curve
  # Implementation follows "L- and V-curves for optimal smoothing"
  # by Gianluca Frasso and Paul HC Eilers
  d3Mat <- cbind(log10(Rs), log10(Ps), lambdas)
  d3MatB <- d3Mat[order(d3Mat[, 1]), ]
  
  # Compute V curve 
  distance <- numeric(length(lambdas) - 1)
  lam <- numeric(length(lambdas) - 1)
  for (i in 1:(length(lambdas) - 1)) {
    distance[i] <- sqrt((d3MatB[i + 1, 2] - d3MatB[i, 2])^2 + (d3MatB[i + 1, 1] - d3MatB[i, 1])^2)
    lam[i] <- (d3MatB[i + 1, 3] + d3MatB[i, 3]) / 2
  }
  
  # Find minimum distance
  I <- which.min(distance)
  
  # run WS with optimal lambda
  z <- WSsparse(t(y), lam[I], d)
  
  return(z)
}

#' Perform weighted WS
#' 
#' @param x signal
#' @param w weights 
#' @param lambda regularization parameter
#' @param differences degree parameter
#' @return background
#' 
#' @author Emery Bosten
#' @export
WhittakerSmooth <- function(x,w,lambda,differences=2) {
  x=matrix(x,nrow = 1, ncol=length(x))
  L=length(x)
  E=Matrix::spMatrix(L,L,i=seq(1,L),j=seq(1,L),rep(1,L))
  D=as(Matrix::diff(E,1,differences),"CsparseMatrix")
  W=as(Matrix::spMatrix(L,L,i=seq(1,L),j=seq(1,L),w),"CsparseMatrix")
  background=Matrix::solve((W+lambda*Matrix::t(D)%*%D),Matrix::t((w*x)));
  return(as.vector(background))
}

#' Perform adaptive iteratively reweighted Penalized Least Squares  
#' 
#' @param x signal
#' @param lambda regularization parameter
#' @param difference degree parameter
#' @param itermax number of iterations
#' @return z baseline
#' 
#' @author Emery Bosten
#' @export
airPLS <- function(x,lambda,differences,itermax){
  
  x = as.vector(x)
  m = length(x)
  w = rep(1,m)
  control = 1
  i = 1
  while(control==1){
    z = WhittakerSmooth(x,w,lambda,differences)
    d = x-z
    sum_smaller = abs(sum(d[d<0])) 
    if(sum_smaller<0.001*sum(abs(x))||i==itermax)
    {
      control = 0
    }
    w[d>=0] = 0
    w[d<0] = exp(i*abs(d[d<0])/sum_smaller)
    w[1] = exp(i*max(d[d<0])/sum_smaller)
    w[m] = exp(i*max(d[d<0])/sum_smaller)
    i=i+1
  }
  return(z) 
}

#' Automated airPLS baseline correction function derived from paper 
#' "An Automatic Baseline Correction Method Based on the Penalized Least Squares Method"  
#' 
#' @param y signal
#' @param lambda regularization parameter
#' @param difference degree parameter
#' @param itermax number of iterations
#' @return err fitting errors 
#' 
#' @author Emery Bosten
#' @export
erPLSrmse <- function(x,lambda,differences, itermax) {
  # Automated erPLS baseline correction function derived from paper "An Automatic Baseline Correction Method Based on
  # the Penalized Least Squares Method"
  
  N <- length(y)
  om <- round(N / 50)
  W <- round(N / 25)
  H <- max(y)
  
  p <- lm(y[(N-om+1):N] ~ seq(1, om))$coefficients
  ye <- p[1] + p[2] * seq(1, om + W)
  yext <- ye[(om + 1):(om + W)]
  
  gauss <- function(x, mean, sd) {
    exp(-((x - mean)^2) / (2 * sd^2))
  }
  
  yg <- (H - yext[W]) * gauss(seq(1, W), W / 2, W / 20)
  
  ynew <- yext + yg
  yadd <- c(y, ynew)
  plot(yadd,type = 'l')
  z <- airPLS(yadd, lambda, differences, itermax)
  
  
  rmse <- function(actual, predicted) {
    sqrt(mean((actual - predicted)^2))
  }
  
  err <- rmse(yext, z[(length(z) - W + 1):length(z)])
  
  return(err)
}

#' Locate minima in second derivative 
#' 
#' @param deriv2 array, second derivative of signal
#' @param y array, original signal
#' @param noise_thresh numeric, peak needs to be at least 'noise_thresh' times higher than estimated noise
#' @param height_ratio numeric, peak needs to be at least 1/'height_ratio' from highest peak
#' @param minimum_depht numeric, minima needs a y value below '-minimum_depht'
#' @return positions of local minima
#' 
#' @author Emery Bosten
findMin <- function(deriv2, y, noise_thresh, height_ratio, minimum_depht) {
  
  # normalization (between -1 and 1)
  max_deriv2 <- max(abs(deriv2)) 
  deriv2 <- deriv2 / max_deriv2 
  
  max_y_noise <- median(diff(y,1))^2 # estimate of noise level in original signal
  max_y <- max(y) # highest peak in original signal
  
  apex_loc <- c()
  for (i in 4:(length(deriv2)-3)) {
    if (deriv2[i-1] > deriv2[i] && deriv2[i-2] > deriv2[i] && deriv2[i-3] > deriv2[i] && 
        deriv2[i+1] > deriv2[i] && deriv2[i+2] > deriv2[i] && deriv2[i+3] > deriv2[i] && 
        deriv2[i] < -minimum_depht && y[i] > max(noise_thresh * max_y_noise, max_y / height_ratio)) {
      
      for (j in 1:20) { # why 20 ?
        if (y[i+j] > y[i + j + 1] && y[i+j] > y[i + j - 1]) {
          i <- i + j
        }
        if (y[i-j] > y[i - j + 1] && y[i-j] > y[i - j - 1]) {
          i <- i - j
        }
      }
      apex_loc <- c(apex_loc, i)
    }
  }
  
  apex_loc <-  unique(apex_loc)
  
  return(apex_loc)
}

#' Automated peak detection in chromotograms 
#' 
#' @param chrom_PDA data.frame, chromatogam in PDA format
#' @param duration numeric, time of gradient (?)
#' @param noise_thresh numeric, peak needs to be at least 'noise_thresh' times higher than estimated noise
#' @param height_ratio numeric, peak needs to be at least 1/'height_ratio' from highest peak
#' @param minimum_depht numeric, minima needs a y value below '-minimum_depht'
#' @return list with the chromatogram, number of peaks, locations, 
#' time axis and plot
#' 
#' @author Emery Bosten
#' @export
peakDetect <-  function(y, x, noise_thresh,
                        height_ratio, minimum_depht){
  
  # First derivative
  dz1 <- diff(y)
  
  # Second derivative
  dz2 <- diff(dz1)
  
  # Locate peaks
  apex_loc <- findMin(dz2, y, noise_thresh, height_ratio, minimum_depht)
  numPeaks <- length(apex_loc)
  
  # Plot chromatogram with detected peaks
  duration = x[length(x)]
  chrom <-  data.frame(x,y)
  peak_pos_x <-  x[apex_loc]
  peak_pos_y <-  y[apex_loc]
  peaks <-  data.frame(peak_pos_x, peak_pos_y)
  plot <-  ggplot(chrom, aes(x = x, y = y)) +
    geom_line(color = "blue", linewidth = 0.75) + labs(x = "Time", y = "Absorbance") +
    geom_point(data = peaks, aes(x = peak_pos_x, y = peak_pos_y), colour = "red")
  
  
  return(list("chrom" = y, 
              "peak_locations" = apex_loc, 
              "number_peaks" = numPeaks, 
              "time_axis" = t, 
              "plot" = plot))
}

#' Generates gaussian peak
#' 
#' @param h numeric, height of peak
#' @param x array of numeric, x-axis (time)
#' @param pos numeric, position of peak
#' @param wid numeric, half-width of peak
#' @return gaussian peak 
#' 
#' @author Emery Bosten
gaussianPeak= function(h,x,pos,wid){
  
  g = h*(exp(-((x-pos)/(0.60056120439323*wid))^2))
  return(g)
}

#' Initializes and show initial conditions for peak fitting
#' 
#' @param h array of numeric, height of peaks
#' @param t array of numeric, x-axis (time)
#' @param pos array of numeric, position of peaks
#' @param wid array of numeric, initial widths of peaks
#' @param numPeaks numeric, number of peaks
#' @return ggplot
#' 
#' @author Emery Bosten
#' @export
showInitCond <-  function(h, t, pos, wid, numPeaks){
  
  peaks <-  matrix(nrow = length(t), ncol = numPeaks, 0) # fit gaussian peaks
  
  # generate initial peaks
  for (i in 1:numPeaks){
    peaks[,i] <-  gaussianPeak(h = h[i], x = t, pos = pos[i], wid = wid[i])
  }
  
  ## plot peaks on top of each other
  peaks <-  as.data.frame(peaks)
  peaks <-  cbind(peaks, t)
  #peaks <-  cbind(peaks, y)
  df <-  melt(peaks, id = "t")
  
  plot <-  ggplot() + geom_line(data = df, aes(x = t, y = value, color = variable)) + 
    theme(legend.position = "none")
  
  return(plot)
}
#' Fits the gaussian peaks to the chrom data
#' 
#' @param t array of numeric, x-axis (time)
#' @param y array of numeric, chromotagram
#' @param numPeaks numeric, number of Peaks
#' @param h array of numeric, height of peaks
#' @param pos array of numeric, position of peaks
#' @param wid numeric, initial half-width of peaks
#' @return ggplot with fitted peaks and coefficients 
#' 
#' @author Emery Bosten
#' @export
fitPeaks <-  function(t, y, numPeaks, h , pos, wid){
  
  # fit peaks to data
  df <- data.frame(t, y) 
  
  if(numPeaks > 0){
    
    # generate formula for model
    hnam <- paste0("h", 1:numPeaks) # height parameter 
    posnam <- paste0("pos", 1:numPeaks) # position parameter
    wnam <- paste0("wid", 1:numPeaks) # width parameter
    
    # paste into single terms
    terms = paste("gaussianPeak","(", hnam,",", "t",",", posnam,",", wnam,")")
    # combine terms into formula
    formula <- as.formula(paste("y ~ ", paste(terms, collapse= "+")))
    
    # generate list of starting conditions
    start_conditions <-  list()
    names <-  list()
    for (i in 1:numPeaks){
      hs = hnam[i]
      ps = posnam[i]
      ws = wnam[i]
      start_conditions = c(start_conditions, list(hs = h[i], ps = pos[i], ws = wid[i]))
      names = c(names, c(hs, ps, ws))
    }
    names(start_conditions) <- names
    
    # model fitting
    fit <- tryCatch(nls(formula, data = df,
                        start = start_conditions, algorithm = "port"), error = function(e) e)
    
    if(inherits(fit, "error")){
      
      coeffs <- NULL
      plot <- ggplot(data = df, aes(x = t, y = y)) +
        geom_line(color = "blue", linewidth = 0.75) + labs(x = "Time", y = "Absorbance") 
      
      #showNotification(sprintf("Failed to fit peaks: %s", fit), type = "error")
      
      
    } else {
      # fitted heights, positions, and widths
      coeffs <-  coef(fit)
      coeffs <-  matrix(coeffs, 3, numPeaks, byrow = F)
      
      # Predict the fitted model
      df_fit <- data.frame(t)
      df_fit$y <- predict(fit, newdata = df_fit)
      
      # Plot the data with the model superimposed
      
      #      plot <- ggplot(data = df, aes(x = t, y = y)) +
      #          geom_point(size = 0.5) +
      #          geom_smooth(data = df_fit, stat = "identity", color = "red",  linewidth = 0.5)
      
      plot <- ggplot(data = df, aes(x = t, y = y)) +
        geom_point(size = 0.5) +
        geom_line(data = df_fit, color = "red",  linewidth = 0.5)
    }
    
  } else {
    
    coeffs <- NULL
    plot <- ggplot(data = df, aes(x = t, y = y)) +
      geom_line(color = "blue", linewidth = 0.75) + labs(x = "Time", y = "Absorbance") 
    
  }
  
  return(list(
    "coeffs" = coeffs, 
    "plot" = plot))
}


#' Computes resolution between two 'gaussian' peaks 
#' 
#' @param rt1 numeric, retention time compound 1
#' @param rt2 numeric, retention time compound 2
#' @param w1 numeric, half-width compound 1
#' @param w2 numeric, half-width compound 2
#' @return resolution 
#' 
#' @author Emery Bosten
#' @export
resolution <-  function(rt1, rt2, w1, w2){
  
  # wiki :  resolution is a measure of the separation of two peaks 
  # of different retention time t in a chromatogram
  resolution <-  abs(1.18 * (rt2 - rt1) / (w2 + w1))
  
  return(resolution)
}

#' Computes the resolution of a critical pair of a chromatogram and
#' the median resolution
#' 
#' @param coeffs array of numeric, coefficients of fitted peaks
#' @return critical resolution and median resolution
#' 
#' @author Emery Bosten
#' @export
criticalResolution <-  function(coeffs){
  
  crit_res <- NA
  med_res <-  NA
  
  pairs <-  ncol(coeffs) - 1
  if (length(pairs) == 0){
    pairs = 0
  }
  if(pairs > 0 ){
    resolutions <-  rep(0, pairs)
    # compute resolution for each peak pair
    for (i in 1:pairs){
      resolutions[i] <-  resolution(coeffs[2, i], coeffs[2, i+1], coeffs[3, i], coeffs[3, i+1])
    }
    
    crit_res <-  min(resolutions)  # compute critical resolution
    med_res <-  median(resolutions) # compute median resolutions 
  }
  
  return(list(
    "crit_res" = crit_res, 
    "med_res" = med_res))
}

#' Process chromatogram
#' 
#' @param x time array
#' @param y chrom
#' @param baselineCorrect True or False for baseline correction
#' @param height_ratio peak detection parameter
#' @param minimum_depht peak detection parameter
#' @param width peak fitting parameter 
#' @return number of Peaks, critical resolution and median resolution
#' 
#' @author Emery Bosten
#' @export
peakProcess = function(x,y,baselineCorrect=F,height_ratio=150,minimum_depht=0.005, width = 0.05){
  ########################
  ## DATA PREPROCESSING ##
  ########################
  
  ## optional ##
  if(baselineCorrect){
    # Baseline correction
    base_lambdas <- logspace(1, 12, 10)
    d <- 2
    n <- 20
    err <- numeric(length(base_lambdas))
    
    # Compute errors for each lambda
    for (i in seq_along(base_lambdas)) {
      lambda <- base_lambdas[i]
      err[i] <- erPLSrmse(y, lambda, d, n)
    }
    
    # Find the minimum error and corresponding index
    m <- min(err)
    ind <- which.min(err)
    
    # Perform baseline correction using the optimal lambda
    z <- airPLS(y, base_lambdas[ind], d, n)
    
    y = y - z
    plot(x,y,type="l")
    
    # Noise reduction
    WS_lambdas <- logspace(-3, 5, 25)
    ys <- VcurveWS(y, WS_lambdas) # smooth signal
    
    # Plot preprocessing the results
    plot(x, ys, type = "l", col = "blue", main = "Baseline Correction", xlab = "x", ylab = "y")
    y = ys
  }
  #####################
  ## DATA PROCESSING ##
  #####################
  
  ## automated peak detection ##
  noise_thresh = 10 # noise threshold 
  return_list = peakDetect(y,x, noise_thresh, height_ratio, minimum_depht)
  
  apex_loc = return_list$peak_location # locations of peak apex
  numPeaks = return_list$number_peaks # number of Peaks
  
  ## fitting peak model ##
  
  # initialize and visualize starting conditions
  w = rep(width,numPeaks) # peak widths  
  h = y[apex_loc] # peak heights
  pos = x[apex_loc] # peak positions 
  
  # fit peaks to data
  return_list = fitPeaks(x, y, numPeaks, h , pos, w)
  coeffs = return_list$coeffs
  
  # critical resolution and median resolution 
  result_list = criticalResolution(coeffs)
  crit_res = result_list$crit_res
  med_res = result_list$med_res
  
  return(list("numPeaks" = numPeaks, "critRes" = crit_res, "medianRes" = med_res))
}


###################
## SET DIRECTORY ##
###################
setwd("C:/Users/emery/OneDrive - KU Leuven/R scripts/Peak detection/data/IL-17/MeOH/")


##########
## DATA ##
##########

## load data ##
pathnames = c("pH 3.07/C18", "pH 3.07/CSH 18", "pH 3.07/Fluoro-Phenyl", "pH 3.07/Phenyl", "pH 4.8/C18", "pH 4.8/CSH 18", "pH 4.8/Fluoro-Phenyl", "pH 4.8/Phenyl")
nconds = length(pathnames)
filenames = c("APIlight.arw", "Dia.arw", "SM.arw", "Blank2.arw")
chroms = list()
for (i in 1:nconds){
  chrom = combineChroms(paste0("./",pathnames[i]), filenames)
  y = chrom$abs # noisy signal
  x = chrom$t
  # Remove negative values
  y = replace(y, y<0, 0)
  plot(x,y,type = "l")
  chrom = data.frame(x,y)
  chroms=c(chroms,chrom)
} 

##############################
## PARAMETER INITIALIZATION ##
##############################

## Reference chromatogram ## 
x = chroms[13]$x
y = chroms[14]$y

########################
## DATA PREPROCESSING ##
########################

## optional ##
baselineCorrect = F
if (baselineCorrect){
  # Baseline correction
  base_lambdas <- logspace(1, 12, 10)
  d <- 2
  n <- 20
  err <- numeric(length(base_lambdas))
  
  # Compute errors for each lambda
  for (i in seq_along(base_lambdas)) {
    lambda <- base_lambdas[i]
    err[i] <- erPLSrmse(y, lambda, d, n)
  }
  
  # Find the minimum error and corresponding index
  m <- min(err)
  ind <- which.min(err)
  
  # Perform baseline correction using the optimal lambda
  z <- airPLS(y, base_lambdas[ind], d, n)
  
  y = y - z
  plot(x,y,type="l")
  
  # Noise reduction
  WS_lambdas <- logspace(-3, 5, 25)
  ys <- VcurveWS(y, WS_lambdas) # smooth signal
  
  # Plot preprocessing the results
  plot(x, ys, type = "l", col = "blue", main = "Baseline Correction", xlab = "x", ylab = "y")
  y = ys
}

#####################
## DATA PROCESSING ##
#####################

## automated peak detection ##
noise_thresh = 10 # noise threshold 
height_ratio = 150 # ratio with highest peak
minimum_depht = 0.005 # threshold depth of minimum

return_list = peakDetect(y,x, noise_thresh, height_ratio, minimum_depht)

apex_loc = return_list$peak_location # locations of peak apex
numPeaks = return_list$number_peaks # number of Peaks
g = return_list$plot # plot showing chromatogram with peak locations
#plot chrom
g

## fitting peak model ##

# initialize and visualize starting conditions
width = 0.05 # initial peak width 
w = rep(width,numPeaks) # peak widths  
h = y[apex_loc] # peak heights
pos = x[apex_loc] # peak positions 

g = showInitCond(h, x, pos, w, numPeaks)
g
# fit peaks to data
return_list = fitPeaks(x, y, numPeaks, h , pos, w)
g = return_list$plot
g
coeffs = return_list$coeffs

# Compute chromatogram metrics 
# number of peaks  
numPeaks

# critical resolution and median resolution 
result_list = criticalResolution(coeffs)
crit_res = result_list$crit_res
med_res = result_list$med_res

# critical resolution
crit_res
#median resolution
med_res

###################################################################
## Application on entire chromatogram set with preset parameters ##
###################################################################
Metrics = matrix(0,3,nconds)
for (i in 1:nconds){
  x = chroms[i*2-1]$x
  y = chroms[i*2]$y
  
  result = peakProcess(x,y,baselineCorrect,height_ratio,minimum_depht,width)
  Metrics[1,i] =  result$numPeaks
  Metrics[2,i] =  result$critRes
  Metrics[3,i] =  result$medianRes
  print(i/nconds)
}
Metrics

