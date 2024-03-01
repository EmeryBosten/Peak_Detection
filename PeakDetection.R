## Script to perform automated peak detection ## 

###############
## LIBRARIES ##
###############

library(ggplot2)
library(reshape2)

###############
## FUNCTIONS ##
###############

## openEmpowerFile ##
# function to open Empower PDA file
# filename = name of file
openEmpowerFile = function(filename){
data = read.table(filename, skip = 4, sep = "")
return(data)  
}
  
## combineChroms ##
# function to combine multiple samples into single chromatogram
# Blank sample must be referenced last
# path = path to files
# filenames = name of files
combineChroms = function(path,filenames){

combinedChrom <- openEmpowerFile(paste0(path, "/", filenames[1]))

for (i in 2:(length(filenames) - 1)) {
  chrom <- openEmpowerFile(paste0(path, "/", filenames[i]))
  if (nrow(combinedChrom) == nrow(chrom)) {
    combinedChrom <- combinedChrom + chrom
  } else if (nrow(combinedChrom) > nrow(chrom)) {
    combinedChrom <- combinedChrom[1:nrow(chrom), ] + chrom
  } else {
    combinedChrom <- combinedChrom + chrom[1:nrow(combinedChrom), ]
  }
}

blank <- openEmpowerFile(paste0(path, "/", filenames[length(filenames)]))

if (nrow(combinedChrom) == nrow(blank)) {
  combinedChrom <- combinedChrom - (length(filenames) - 1) * blank
} else if (nrow(combinedChrom) > nrow(blank)) {
  combinedChrom <- combinedChrom[1:nrow(blank), ] - (length(filenames) - 1) * blank
} else {
  combinedChrom <- combinedChrom - (length(filenames) - 1) * blank[1:nrow(combinedChrom), ]
}

return(combinedChrom)
}

## findMin ## 
# function to locate minima in second derivative 
findMin <- function(deriv2, y, thresh, thresh2, thresh3) {
  # function to find peak apex based on minima in second derivative
  # deriv2 = second derivative of signal
  # y = original signal
  # thresh = peak needs to be at least 'thresh' times higher than estimated noise
  # thresh2 = peak needs to be at least 1/'thresh2' from highest peak
  # thresh3 = minima needs a y value below '-thresh3'
  max_deriv2 <- max(abs(deriv2)) # absolute max value in second deric signal
  deriv2 <- deriv2 / max_deriv2 # normalization (between -1 and 1)
  max_y_noise <- max(c(y[1:25], y[(length(y)-25):length(y)])) # estimate of noise level in original signal
  max_y <- max(y) # highest peak in original signal
  ind <- 0
  apex_loc <- c()
  for (i in 4:(length(deriv2)-3)) {
    if (deriv2[i-1] > deriv2[i] && deriv2[i-2] > deriv2[i] && deriv2[i-3] > deriv2[i] && 
        deriv2[i+1] > deriv2[i] && deriv2[i+2] > deriv2[i] && deriv2[i+3] > deriv2[i] && 
        deriv2[i] < -thresh3 && y[i] > max(thresh * max_y_noise, max_y / thresh2)) {
      ind <- ind + 1
      for (j in 1:20) {
        if (y[i+j] > y[i + j + 1] && y[i+j] > y[i + j - 1]) {
          i <- i + j
        }
        if (y[i-j] > y[i - j + 1] && y[i-j] > y[i - j - 1]) {
          i <- i - j
        }
      }
      apex_loc <- c(apex_loc, i)
      apex_loc = unique(apex_loc)
    }
  }
  return(apex_loc)
}

## peakDetect ## 
# function to perform automated peak detection
# chrom_PDA = chromatogam in PDA format
# duration = time of gradient
# noise_thresh = peak needs to be at least 'thresh' times higher than estimated noise
# height_ratio = peak needs to be at least 1/'thresh2' from highest peak
# minimum_depht = minima needs a y value below '-thresh3'
max_deriv2 <- max(abs(deriv2)) # absolute max value in second deric signal
peakDetect = function(chrom_PDA, duration, noise_thresh, height_ratio, minimum_depht){

  # Peak detection
  y = rowMeans(chrom_PDA)
  x <- seq(0, length(y), length.out = length(y))
  
  # First derivative
  dz1 <- diff(y)
  
  # Second derivative
  dz2 <- diff(dz1)
  
  # Locate peaks
  apex_loc <- findMin(dz2, y, noise_thresh, height_ratio, minimum_depht)
  numPeaks <- length(apex_loc)
  
  # Plot chromatogram with detected peaks
  t = x/length(x) * duration
  chrom = data.frame(t,y)
  peak_pos_x = x[apex_loc]/length(x) * duration
  peak_pos_y = y[apex_loc]
  peaks = data.frame(peak_pos_x, peak_pos_y)
  g = ggplot(chrom, aes(x=t, y=y)) +
  geom_line(color = "blue", linewidth = 0.75) + labs(x = "Time", y = "Absorbance") +
  geom_point(data = peaks, aes(x = peak_pos_x, y = peak_pos_y), colour = "red")
  g
  return_list = list("chrom" = y, "peak_locations" = apex_loc, "number_peaks" = numPeaks, "time_axis" = t, "plot" = g)
  return(return_list)
}

## gaussian ##
# function to generate gaussian peak 
# h = height
# x = x-axis
# pos = position
# wid = half-width
gaussian = function(h,x,pos,wid){
  g = h*(exp(-((x-pos)/(0.60056120439323*wid))**2))
  return(g)
}

## visuInitCond ##
#initialize and visualize initial conditions for peak fitting
# h = height
# x = x-axis
# pos = position
# wid = half-width
# numPeaks = number of Peaks
visuInitCond = function(w, h, pos,t,numPeaks){
  peaks = matrix(nrow = length(t),ncol = numPeaks,0) # fit gaussian peaks
  # generate initial peaks
  for (i in 1:numPeaks){
    peaks[,i] = gaussian(h[i],t,pos[i],w[i])
  }
  ## plot peaks on top of each other
  peaks = as.data.frame(peaks)
  peaks = cbind(peaks,t)
  peaks = cbind(peaks,y)
  df = melt(peaks,id="t")
  
  g = ggplot() + geom_line(data = df, aes(x = t, y = value, color = variable)) + 
    theme(legend.position = "none")
  return(g)
}

## fitPeaks ##
# function to fit the gaussian peaks to the chrom data
# t = time
# y = signal
# numPeaks = number of Peaks
fitPeaks = function(t,y,numPeaks){
  # fit peaks to data
  df <- data.frame(t, y) 
  
  # generate formula for model
  hnam <- paste0("h", 1:numPeaks) # height parameter 
  posnam <- paste0("pos", 1:numPeaks) #position parameter
  wnam <- paste0("w", 1:numPeaks) #width parameter
  
  # paste into single terms
  terms = paste("gaussian","(",hnam,",", "t",",", posnam,",", wnam,")")
  # combine terms inti formula
  fmla <- as.formula(paste("y ~ ", paste(terms, collapse= "+")))
  
  # generate list of starting conditions
  start_conditions = list()
  namelist = list()
  for (i in 1:numPeaks){
    hs = hnam[i]
    ps = posnam[i]
    ws = wnam[i]
    start_conditions = c(start_conditions,list(hs=h[i], ps=pos[i], ws=w[i]))
    namelist = c(namelist,c(hs,ps,ws))
  }
  names(start_conditions) <- namelist
  start_conditions
  
  # model fitting
  fit <- nls(fmla, data=df,
             start= start_conditions, algorithm="port")  
  
  # fitted heights, positions, and widths
  coeffs = coef(fit)
  coeffs = matrix(coeffs, 3, numPeaks, byrow = F)
  
  # Predict the fitted model
  dffit <- data.frame(t)
  dffit$y <- predict(fit, newdata=dffit)
  
  # Plot the data with the model superimposed
  g = ggplot(df, aes(x=t, y=y)) + geom_point(size = 0.5) +
    geom_smooth(data=dffit, stat="identity", color="red", linewidth=0.5)
  return_list = list("coeffs"= coeffs, "plot"= g)
  return(return_list)
}

## resolution ## 
# function that computes resolution between two 'gaussian' peaks
# with half-widths
# rt1 = retention time compound 1
# rt2 = retention time compound 2
# w1 = half-width compound 1
# wé = half-width compound 2
res = function(rt1, rt2, w1, w2){
  resolution = 1.18*(rt2 - rt1)/(w2 + w1)
  return(resolution)
}

## Res ## 
# function that computes the resolution of the critical pair of a chrom
# and the median resolution
# coeffs = coefficients of fitting object
ResMet = function(coeffs){
  pairs = ncol(coeffs) - 1
  resolutions = rep(0,pairs)
  # compute resolution for each peak pair
  for (i in 1:pairs){
    r = res(coeffs[2,i], coeffs[2, i+1], coeffs[3,i], coeffs[3,i+1])
    resolutions[i] = r 
  }
  crit_res = min(resolutions) # compute critical resolution
  med_res = median(resolutions) # compute median resolutions 
  result_list = list("crit_res" = crit_res, "med_res" = med_res)
  return(result_list)
}

###################
## SET DIRECTORY ##
###################
setwd("C:/Users/emery/OneDrive - KU Leuven/R scripts/Peak detection")


##########
## DATA ##
##########

## load data ##
filenames = c("DS.arw", "Mix.arw", "Blank.arw")
chrom_PDA = combineChroms("./data", filenames)
#plot chrom
plot(rowMeans(chrom_PDA), type = 'l') 

#####################
## DATA PROCESSING ##
#####################

## automated peak detection ##

duration = 50 # time of gradient
noise_thresh = 2 # noise threshold 
height_ratio = 100 # ratio with highest peak
minimum_depht = 0.0005 # threshold depth of minimum

return_list = peakDetect(chrom_PDA,duration, noise_thresh, height_ratio, minimum_depht)

y = return_list$chrom #chromatogram
apex_loc = return_list$peak_location # locations of peak apex
numPeaks = return_list$number_peaks # number of Peaks
t = return_list$time_axis # time axis
g = return_list$plot # plot showing chromatogram with peak locations
#plot chrom
g

## fitting peak model ##

# initialize and visualize starting conditions
width = 0.15 # initial peak width 
w = rep(width,numPeaks) # peak widths  
h = y[apex_loc] # peak heights
pos = t[apex_loc] # peak positions 

g = visuInitCond(w, h, pos,t,numPeaks)
g

# fit peaks to data
return_list = fitPeaks(t,y,numPeaks)
g = return_list$plot
g
coeffs = return_list$coeffs

# Compute chromatogram metrics 
# number of peaks  
numPeaks
# critical resolution and median resolution 
result_list = ResMet(coeffs)
crit_res = result_list$crit_res
med_res = result_list$med_res
# critical resolution
crit_res
#median resolution
med_res

