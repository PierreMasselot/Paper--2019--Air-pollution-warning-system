#' Specify complex divisions of background color
bg_color <- function(xbreaks = NULL, ybreaks = NULL, col = NA, transp = 0, 
  border = NA, ...){
  # Add plot limits to breaks
  plot_lims <- par("usr")
  xbreaks <- c(plot_lims[1], xbreaks, plot_lims[2])
  nx <- length(xbreaks) - 1
  ybreaks <- c(plot_lims[3], ybreaks, plot_lims[4])
  ny <- length(ybreaks) - 1
  # Potentially add transparency to colors
  transp <- rep_len(transp, nx * ny)
  col <- rep_len(col, nx * ny)
  col <- mapply(transparency, col, transp)
  # Rectangles coordinates
  gr <- expand.grid(1:nx, 1:ny)
  rect(xbreaks[gr[,1]], ybreaks[gr[,2]], xbreaks[gr[,1] + 1], 
    ybreaks[gr[,2] + 1], col = col, border = border, ...)
  box()
}

#' Make a color transparent
transparency <- function(col, alpha = 0){
  rgb.val <- col2rgb(col)
  transp.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, 
    alpha = (1 - alpha) * 255)
  return(transp.col)
}

#' Add axis as intervals
axis.intervals <- function(side=1, ticks = axTicks(side), atLabels = NULL, labels = 1:length(atLabels), ...)
{
  stopifnot((side <- as.integer(side)) %in% 1:4)
  if (length(ticks) == 1){
     is.x <- side%%2 == 1
     usr <- par("usr")[1:2 + 2*!is.x]
     XY <- function(ch) paste0(if (is.x) "x" else "y", ch)
     axs <- par(XY("axs"))
     if (axs == "r"){
        per4 <- 4*diff(usr)/108
        usr[1] <- usr[1] + per4
        usr[2] <- usr[2] - per4
     }
     nIntervals <- floor(ticks)
     ticks <- seq(usr[1],usr[2],length.out = nIntervals + 1)
  } else {
     nIntervals <- length(ticks) - 1
  }    
  axis(side, at = ticks, labels = FALSE, ...)
  if (is.null(atLabels)){
     atLabels <- (ticks[-1] + ticks[1:nIntervals]) / 2
  } else {
     if (length(atLabels) != nIntervals) atLabels <- rep_len(atLabels,nIntervals)
  }
  if (length(labels) != nIntervals) labels <- rep_len(labels,nIntervals)
  axis(side, at = atLabels, labels = labels, tick = FALSE, ...)
}

#' Addition of a legend in the margin of an existing plot 
outerLegend <- function(location,...)
{
   if (length(dev.list())==0) stop("Aucun graphique ouvert") 
   olparams <- switch(tolower(location),
       topcenter = list(originalLocation = "top", direction = c(0,1)),
       topright = list(originalLocation = "topright", direction = c(0,1)),
       topleft = list(originalLocation = "topleft", direction = c(0,1)),
       rightcenter = list(originalLocation = "right", direction = c(1,0)),
       righttop = list(originalLocation = "topright", direction = c(1,0)),
       rightbottom = list(originalLocation = "bottomright", direction = c(1,0)),
       leftcenter = list(originalLocation = "left", direction = c(1,0)),
       lefttop = list(originalLocation = "topleft", direction = c(1,0)),
       leftbottom = list(originalLocation = "bottomleft", direction = c(1,0)),
       bottomcenter = list(originalLocation = "bottom", direction = c(0,1)),
       bottomright = list(originalLocation = "bottomright", direction = c(0,1)),
       bottomleft = list(originalLocation = "bottomleft", direction = c(0,1)),
       stop("Paramètre 'location' non reconnu")
   )
   leg <- legend(olparams$originalLocation,plot=F,...) 
   insetx <- leg$rect$w/(par("usr")[2]-par("usr")[1])
   insety <- leg$rect$h/(par("usr")[4]-par("usr")[3])
   insetpar <- -c(insetx,insety) * olparams$direction
   legend(olparams$originalLocation,xpd=NA,inset=insetpar,...)
}