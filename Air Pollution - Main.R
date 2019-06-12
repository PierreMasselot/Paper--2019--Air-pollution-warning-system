##########################################################################
#
#                   Script of the paper:
#        Toward an air pollution warning system in 
#                 the province of Quebec
#
#            Authors: Pierre Masselot, Charlotte Billy
#
#                      June 2019
#
##########################################################################

# Necessary packages to run the analysis
library(dlnm) # DLNM
library(splines) # For the 'ns' function
library(devtools) # For installing the hhws package
library(quantmod) # For the Lag function
library(mgcv) # For the DLNM
library(hhws) # To install: devtools::install_github("PierreMasselot/hhws")

# Working directory
# setwd("C:/...")

# Sourcing some plotting functions
source("Other functions.R")

#------------------------------------------------
# Part 0: DATA LOADING
#------------------------------------------------

#---- PARAMETERS
# Region of interest
region <- "CMM"
# Seasons definition (months) 
month <- list(Summer = c(5:9), Winter = c(10:12, 1:4))
# Number of seasons
ns <- length(month)

#---- LOADING DATA
dataread <- read.table(sprintf("%sdata.csv", region), header = T, sep = ";")

#---- DATA PREPROCESSING
# mean of Temperature over 3 days
dataread$Tmc <- filter(dataread$Tmean, filter = rep(1/3, 3), sides = 1)
# mean of Relative Humidity over 3 days 
dataread$RHc <- filter(dataread$RH, filter = rep(1/3, 3), sides = 1)
# Creation of a Date variable 
dataread$Date <- as.POSIXlt(with(dataread, 
  sprintf("%i-%02.0f-%02.0f", Year, Month, Day)))
# Day of week variable 
dataread$dow <- weekdays(as.POSIXlt(dataread$Date))

#---- SPLIT DATA TABLE FOR EACH SEASON
# Create list to store results
datatab <- vector("list", ns); names(datatab) <- names(month)
for (s in 1:ns){ # For each season
  # Select observations
  datatab[[s]] <- dataread[dataread$Month %in% month[[s]],]
  # Create a day of season variable
  datatab[[s]]$dos <- sequence(tapply(datatab[[s]]$Date, 
    datatab[[s]]$Year,length))
}

# Number of observations in each season
n <- sapply(datatab, nrow)
# Number of years of data  
ny <- length(unique(dataread$Year)) 

#------------------------------------------------
# Part 1: CHOOSE MAXIMUM LAG FOR AP INDICATORS
#------------------------------------------------
#---- PARAMETERS
# maximum lag to consider
maxlag <- 5 
# DLNM: number of DF for seasonality spline
dfseason <- 4
# DLNM: number of DF for trend spline
dftrend <- round(ny / 10)
# DLNM: number of DF for temperature spline
dftemp <- 3
# Slices at each to display the lag-response relationship (PM2.5, Ox)  
slices <- c(25, 50)  

#---- DISTRIBUTED LAG NONLINEAR MODEL (DLNM)
# Create lists to store results
surfPM <- surfOX <- vector("list", ns) 
names(surfPM) <- names(surfOX) <- names(month)
for (s in 1:ns){ # For each season
  # Select seasonal data
  datas <- datatab[[s]]
  
  ### PM2.5
  # Define crossbasis
  cbPM <- crossbasis(datas$PM25max, lag = maxlag, 
    argvar = list(fun = "ps", df = 9), arglag = list(fun = "ns", knots = 1:2))
  # Set penalization
  penPM <- cbPen(cbPM)
  # Model fitting
  fitPM <- gam(Death ~ cbPM + ns(Date$yday, dfseason) + ns(Year, dftrend) + dow
    + ns(Tmc, dftemp) + RHc, family = quasipoisson(), data = datas,
    paraPen = list(cbPM = penPM), method = "REML")
  # Creation of the final surface
  surfPM[[s]] <- crosspred(cbPM, fitPM, cen = 0, by = 1)
  
  ### Ox
  # Define crossbasis
  cbOX <- crossbasis(datas$OXmax, lag = maxlag, 
    argvar = list(fun = "ps", df = 9), arglag = list(fun = "ns", knots = 1:2))
  # Set penalization
  penOX <- cbPen(cbOX)
  # Model fitting
  fitOX <- gam(Death ~ cbOX + ns(Date$yday, dfseason) + ns(Year, dftrend) + dow
    + ns(Tmc, dftemp) + RHc, family = quasipoisson(), data = datas,
    paraPen = list(cbOX = penOX), method = "REML")
  # Creation of the final surface
  surfOX[[s]] <- crosspred(cbOX, fitOX, cen = 0, by = 1)
}

#---- PLOT SURFACE
x11(title = "FigureS4", width = 10, height = 10)
par(mfrow = c(2,ns), cex.lab = 1.3, cex.main = 1.3)
for (s in 1:ns){ # For each season
  # Plot the PM2.5 surface
  plot(surfPM[[s]], 
    main = substitute(expression(paste(a, ") ", b, ", ", PM[2.5])), 
      list(a = letters[1+(s-1)*ns], b = names(month)[s])),
    xlab = "PM2.5 (µg/m3)", zlab = "RR", ylab = "Lag (days)")
  # Plot the Ox Surface
  plot(surfOX[[s]], 
    main = substitute(expression(paste(a, ") ", b, ", ", O[x])), 
      list(a = letters[2+(s-1)*ns], b = names(month)[s])), 
    xlab = "Ox (ppb)", zlab = "RR", ylab = "Lag (days)")
}
### Saving the plot
dev.print(png, filename = "FigureS4.png", res = 600, 
  width = dev.size()[1], height = dev.size()[2], units = "in") 

# Slices
x11(title = "Figure1", width = 10, height = 8)
par(mfrow = c(2,ns), cex.lab = 1.3, cex.main = 1.5)
for (s in 1:ns){ # For each season
  # Slice of PM2.5 surface
  plot(surfPM[[s]], var = slices[1], ylab = "RR", cex = 1.5, xlab = "Lag (days)", 
    main = substitute(expression(paste(a, ") ", b, ", ", PM[2.5])), 
      list(a = letters[1+(s-1)*ns], b = names(month)[s])), 
    lwd = 3, ci = "bars", ci.arg = list(col = gray(.5), lwd = 2), type = "p", 
    pch = 16, col = "red")
  # Slice of Ox
  plot(surfOX[[s]], var = slices[2], ylab = "RR", cex = 1.5, xlab = "Lag (days)", 
    main = substitute(expression(paste(a, ") ", b, ", ", O[x])), 
      list(a = letters[2+(s-1)*ns], b = names(month)[s])), 
    lwd = 3, col = "red", ci = "bars", ci.arg = list(col = gray(.5), lwd = 2), 
    type = "p", pch = 16)
}
dev.print(png, filename = "Figure1.png", res = 600, 
  width = dev.size()[1], height = dev.size()[2], units = "in")

#------------------------------------------------
# Part 2: COMPUTE EXCESS MORTALITY FROM DATA
#------------------------------------------------
#---- PARAMETERS
# degrees of freedom (/ year) for the baseline spline
df.sp <- 8 

#---- BASELINE (EM)
em <- baseline(dataread$Death, dates = dataread$Date, smoothing.fun = "spline",
  df = df.sp * ny, nyear = 1)

#---- EXCESS MORTALITY FROM BASELINE
# For the whole dataset
OMall <- excess(dataread$Death, em)
# Seasonal separation
OM <- vector("list", ns); names(OM) <- names(month)
for (s in 1:ns){
  OM[[s]] <- OMall[dataread$Month %in% month[[s]]]
}

#------------------------------------------------
# Part 3: DETERMINE OVER-MORTALITY EPISODES
#------------------------------------------------
#---- PARAMETERS
# Maximum threshold sOM to test (summer, winter)
maxOMt <- c(80, 55)
# Preliminary threshold (constraint) on PM2.5
prelim.thresh <- 25
# Number of days separating two episodes to be considered distinct 
r <- 3
# Number of days before and after extreme OM days to be considered as an episode
l <- 3

#---- NUMBER OF EPISODES FOUND FOR DIFFERENT THRESHOLDS ON OM
Nepis <- vector("list", ns); names(Nepis) <- names(month) 
for (s in 1:ns){ # For each season
  # Sequence of tested thresholds (increment of 5)
  OMt.vec <- seq(30, maxOMt[s], by = 5)
  # Number of episode for each element in OMt.vec
  Nepis[[s]] <- sapply(OMt.vec, function(t){
    # With the constraint on PM2.5
    extCov <- episodes(OM[[s]], t, r = r, 
      covariates = cbind(datatab[[s]]$PM25max),#, -datatab$Tmax), 
      uc = prelim.thresh)
    # Without constraint
    ext <- episodes(OM[[s]], t, r = r)
    ret <- c(max(ext$episode), max(extCov$episode))
    names(ret) <- c("Without PM", "With PM")
    return(ret)
  })
  colnames(Nepis[[s]]) <- OMt.vec
}

#---- CHOSEN VALUE
OMt <- c(50, 40)

#---- PLOTTING THE NUMBER OF EPISODES FOR EACH TESTED THRESHOLD
x11(title = "Figure2", width = 10, height = 5)
par(mfrow = c(1,ns))
for (s in 1:ns){ # For each season
  # The curves
  matplot(seq(30, maxOMt[s], by = 5), t(Nepis[[s]]), type = "b", lwd = 2, 
    pch = 15:16, col = c("black", "cornflowerblue"), cex.axis = 1.2, 
    cex.lab = 1.3, cex.main = 1.3, ylab = "Episode number", lty = c(1,3), 
    xlab = expression(paste(s[OM], " (%)")), 
    main = sprintf("%s) %s", letters[s], names(month)[s]))
  # Vertical line indicating the chosen threshod
  abline(v = OMt[s], lwd = 3, lty = 2, col = "indianred")
  text(OMt[s], par("usr")[4], 
    substitute(paste(s[OM], " = ", i, "%"), list(i = OMt[s])), 
    cex = 1.3, col = "indianred", adj = c(1.1, 1.4))
}
legend("topright", expression("Without constraint", 
  paste("PM2.5 > 25", mu, "g/", m^3)), bty = "n", lty = c(1,3), pch = 15:16, 
  cex = 0.8, col = c("black", "cornflowerblue"))

dev.print(png, filename = "Figure2.png", res = 600, 
  width = dev.size()[1], height = dev.size()[2], units = "in")

#---- FINAL EPISODE EXTRACTION
OM.episodes <- vector("list", ns)
for (s in 1:ns){ # For each season
  OM.episodes[[s]] <- episodes(OM[[s]], OMt[s], r = r, l = l, 
    covariates = cbind(datatab[[s]]$PM25max), #, -datatab$Tmax), 
    uc = prelim.thresh)
}

#---- PLOT: EXCESS MORTALITY SERIES WITH IDENTIFIED EPISODES
cols <- c("red", "blue")
x11("Figure3", width = 10, height = 5)
# Initializing the plot
plot(dataread$Date, OMall, type = "n", xlab = "Time", ylab = "OM (%)", 
  cex.lab = 1.3, cex.axis = 1.2, xaxt = "n", ylim = range(OMall) * c(1, 1.1))
season.lims <- dataread$Date[which(dataread$Day == 1 & 
  dataread$Month %in% sapply(month, "[", 1))]
# Background color (custom function)
bg_color(xbreaks = as.numeric(season.lims), 
  col = c("cornflowerblue", "indianred"), transp = .7) # Background color
# Adding the series
lines(dataread$Date, OMall)  # Add OM series
# Adding the OM thresholds
segments(c(par("usr")[1], as.numeric(season.lims)), rev(OMt), 
  c(as.numeric(season.lims), par("usr")[2]), lwd = 2, col = rev(cols)) # Add thresholds
# Adding episodes
for (s in 1:ns){ 
  episodes <- OM.episodes[[s]][OM.episodes[[s]]$extreme,]
  episx <- datatab[[s]]$Date[episodes[,"t"],]
  episy <- episodes[,"value"]
  points(episx, episy, pch = 4, lwd = 2, col = cols[s])
  episMax <- aggregate(cbind(as.numeric(episx), episy), 
    by = list(episode = episodes[,"episode"]), max)
  text(episMax[,2], episMax[,3], as.character(episMax[,1]), 
    pos = 3, xpd = T, col = cols[s], cex = 1)
}
# X-axis
year.lims <- dataread$Date[which(dataread$Day == 1 & dataread$Month == 1)]
axis.intervals(side = 1, ticks = as.numeric(c(year.lims, max(dataread$Date))),
  labels = unique(dataread$Year))
# Legend
outerLegend("topcenter", c("Winter", "Summer"), fill = c("blue", "red"), 
  bty = "n", ncol = 2)  # Legend

dev.print(png, filename = "Figure3.png", res = 600, 
  width = dev.size()[1], height = dev.size()[2], units = "in")

#------------------------------------------------
# Part 4: CHOOSE THE BEST INDICATORS/THRESHOLDS
#------------------------------------------------
#---- PARAMETERS
# Lags chosen on step 1 above, respectively for: PM2.5, Ox
lags <- list(Summer = c(1, 1), Winter = c(1, 1))
# Are the sensitivity and FA computed on days or episodes ?
thinning <- "episodes" 

#---- PREPARE INDICATOR VARIABLES
X <- dataread[,c("PM25max", "OXmax")]
# Create object containing indicators
indicators <- vector("list", ns); names(indicators) <- names(month)
for (s in 1:ns){ # For each season
  indicators[[s]] <- vector("list", 2)
  names(indicators[[s]]) <- c("PM25", "Ox") 
  # Lag pollutant variables
  for (i in 1:2){
    indicators[[s]][[i]] <- sapply(0:lags[[s]][i], Lag, x = X[,i])  # Lag the series
  }
  # Select season
  indicators[[s]] <- lapply(indicators[[s]], "[", dataread$Month %in% month[[s]], ) 
}

#---- TEST THRESHOLDS AND WEIGHTINGS
tested <- vector("list", ns); names(tested) <- names(month) # tested thresholds
for (s in 1:ns){
  # Function computing different criteria for a large array of thresholds
  #   and weightings
  tested[[s]] <- find.threshold(indicators[[s]], episodes = OM.episodes[[s]], 
    u.grid = list(20:60, 20:60), order.result = "Episodes_found", 
    thinning = thinning, same.alphas = FALSE, trim = 30)
}

#---- FINAL THRESHOLDS AND WEIGHTINGS
# Lines in tested that correspond to the chosen APHWS
chosen.lines <- c(2, 6)
# Select the lines 
final <- vector("list", ns); names(final) <- names(month)
for (s in 1:ns) final[[s]] <- tested[[s]][chosen.lines[s],]

#---- OBTAIN FINAL ALARMS
alarms <- vector("list", ns); names(alarms) <- names(month)
for (s in 1:ns){
  # Thresholds values
  thresh <- final[[s]][grep("threshold_", names(final[[s]]))]
  # Weightings values
  alpha <- split(unlist(final[[s]][grep("alpha", names(final[[s]]))]), 
    rep(1:2, lags[[s]] + 1))
  # Alarms
  alarms[[s]] <- predict_alarms(indicators[[s]], alpha = alpha, s = thresh, 
    y = OM[[s]], r = 3)
  # Alarm dates
  alarms[[s]]$Date <- datatab[[s]]$Date[alarms[[s]][,"t"],] 
}

#------------------------------------------------
# SENSITIVITY ANALYSIS
#------------------------------------------------
# Prepare objects stocking results
onlyPM25 <- onlyOx <- both <- onlyPM25.l3 <- onlyOx.l3 <- both.l3 <- 
  vector("list", ns)

for (s in 1:ns){ # For each season
  #---- EVALUATE SYSTEM WITH SINGLE POLLUTANTS
  onlyPM25[[s]] <- find.threshold(indicators[[s]][[1]], episodes = OM.episodes[[s]], 
    u.grid = list(seq(20, max(datatab[[s]]$PM25max, na.rm = T), length.out = 50)), 
    order.result = "False_episodes", thinning = "episodes")
  onlyOx[[s]] <- find.threshold(indicators[[s]][[2]], episodes = OM.episodes[[s]], 
    u.grid = list(seq(20, max(datatab[[s]]$OXmax, na.rm = T), length.out = 50)), 
    order.result = "False_episodes", thinning = "episodes")
  both[[s]] <- find.threshold(indicators[[s]], episodes = OM.episodes[[s]], 
    u.grid = lapply(X, function(x) seq(20, max(x, na.rm = T), length.out = 50)), 
    order.result = "False_episodes", thinning = "episodes", same.alphas = FALSE)  
    # A bit long
      
  #---- EVALUATE SYSTEM WITH L = 2
  # Prepare indicators
  indic3 <- vector("list", length(lags))
  for (i in 1:length(lags)){
    indic3[[i]] <- sapply(0:3, Lag, x = X[,i])  # Lag the series
  }
  indic3 <- lapply(indic3, "[", dataread$Month %in% month[[s]], ) # Keep summer
  names(indic3) <- c("PM25", "Ox")
  
  # Compute values 
  onlyPM25.l3[[s]] <- find.threshold(indic3[[1]], episodes = OM.episodes[[s]], 
    u.grid = list(seq(20, max(datatab[[s]]$PM25max, na.rm = T), length.out = 50)), 
    order.result = "False_episodes", thinning = "episodes")
  onlyOx.l3[[s]] <- find.threshold(indic3[[2]], episodes = OM.episodes[[s]], 
    u.grid = list(seq(20, max(datatab[[s]]$OXmax, na.rm = T), length.out = 50)), 
    order.result = "False_episodes", thinning = "episodes")
  both.l3[[s]] <- find.threshold(indic3, episodes = OM.episodes[[s]], 
    u.grid = lapply(X, function(x) seq(20, max(x, na.rm = T), length.out = 20)), 
    order.result = "False_episodes", thinning = "episodes", 
    same.alphas = FALSE)
    # A bit long
}

#---- PLOT THE ROC CURVE FOR EACH TYPE OF SYSTEM
xlims <- c(100, 200)

x11(title = "Figure4", width = 14)
par(mfrow = c(1,2))
for (s in 1:ns){ # For each season
  # Two pollutant and L = 1
  plot(both[[s]][,"False_episodes"], both[[s]][,"Episodes_sensitivity"], 
    ylab = "Sensitivity (Episodes)", xlab = "Number of false episodes", 
    pch = 16, cex.lab = 1.5, cex.axis = 1.3, lwd = 1.5, xlim = c(0, xlims[s]), 
    ylim = c(0, 1), type = "b", cex.main = 1.5, 
    main = sprintf("%s) %s", letters[s], names(month)[s]))
  # Single pollutant systems and L = 1
  points(onlyPM25[[s]][,"False_episodes"], onlyPM25[[s]][,"Episodes_sensitivity"], 
    pch = 15, col = "forestgreen", type = "b", lwd = 1.5)
  points(onlyOx[[s]][,"False_episodes"], onlyOx[[s]][,"Episodes_sensitivity"], 
    pch = 17, col = "cornflowerblue", type = "b", lwd = 1.5)
  # L = 2 systems    
  points(both.l3[[s]][,"False_episodes"], both.l3[[s]][,"Episodes_sensitivity"], 
    pch = 1, lwd = 1.5, col = "black", type = "b", lty = 2)
  points(onlyPM25.l3[[s]][,"False_episodes"], onlyPM25.l3[[s]][,"Episodes_sensitivity"],
    pch = 0, lwd = 1.5, col = "forestgreen", type = "b", lty = 2)
  points(onlyOx.l3[[s]][,"False_episodes"], onlyOx.l3[[s]][,"Episodes_sensitivity"], 
    pch = 2, lwd = 1.5, col = "cornflowerblue", type = "b", lty = 2)
  # Final system   
  points(final[[s]][,"False_episodes"], final[[s]][,"Episodes_sensitivity"], 
    col = "red", pch = 4, lwd = 2, cex = 2)  
  # Legend
  legend("bottomright", c("PM2.5 only", "Ox only", "PM2.5 and Ox", 
    "Chosen system", "PM2.5: 3 days", "Ox: 3 days", "PM2.5, Ox: 3 days"), 
    pch = c(15, 17, 16, 4, 0, 2, 1), col = c("forestgreen", 
    "cornflowerblue", "black", "red"), bty = "n", ncol = 2, bg = "white") 
}

dev.print(png, filename = "Figure4.png", res = 600, 
  width = dev.size()[1], height = dev.size()[2], units = "in")

