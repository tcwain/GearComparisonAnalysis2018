#!/usr/bin/R

#Plume Project -- Length-Frequency Comparison of MMED trawl
#  comparisons Combined Data
# Sampling described in RC Dotson etal 2010 NOAA-TM-NMFS-SWFSC-455,
#  and NWFSC Plume Project cruise reports

# data.dir <- '.'  #Directory where MMED data resides

## ---- GetLenFreqData
# Restrict to species selected above, but no age groups for salmon
len.spec <- c("CHINOOK SALMON", "CHUM SALMON", "COHO SALMON",
             "MARKET SQUID", "NORTHERN ANCHOVY", "PACIFIC HERRING",
             "SEA NETTLE", "WATER JELLY")
lenData <- MMEDdata[ , c("Cruise","MMED","Species","Length","Number")]
lenData <- lenData[lenData$Species %in% len.spec, ]
lenData <- lenData[!is.na(lenData$Length), ]
lenData$Len5mm <- 5 * round(lenData$Length/5)
## print(summary(lenData)) ### NOT RUN

## ---- LenFreqAnal

for (excl in c("Up","Down")) {
  cat('**************** Excluder: ', excl, ' *****************\n')
  lD <- if(excl %in% "Up") {
    lenData[lenData$Cruise %in% c(41,43,50), ]
  } else {
    lenData[lenData$Cruise %in% 53, ]
  }
  lD$Species <- factor(as.character(lD$Species))
  .tab <- with(lD, tapply(Number, list(Species), sum, na.rm=T))
  print(.tab)
  lf.sel.spec <- names(.tab)[.tab >= 100]
  lenFreq <- with(lD, tapply(Number, list(Len5mm, MMED, Species), sum))
  .mfrow <- if(lndscp) c(3,3) else c(3,2)
  par(mfrow=.mfrow, omi=c(0.5,0.5,0,0), mar=c(3,3,2,1))
  for (sp in lf.sel.spec) {
    cat('****************', sp, '*****************\n')
    if (sp %in% dimnames(lenFreq)[[3]]) {
      .dat <- lenFreq[ , , sp]
      .maxN <- max(.dat,na.rm=T)
      .len <- lD[lD$Species %in% sp, ]
      .maxL <- max(.len$Length, na.rm=T)
      .minL <- min(.len$Length, na.rm=T)
      plot(as.numeric(rownames(.dat)), .dat[, 'None'], type='h', col=BLACK, lwd=1,
           xlim=c(.minL,.maxL), ylim=c(-.maxN,.maxN), axes=F,
           xlab='', ylab='')
      axis(side=1, lwd=0, lwd.ticks=1)
      axis(side=2, at=pretty(c(-max(.dat,na.rm=T),max(.dat,na.rm=T))),
           labels=abs(pretty(c(-max(.dat,na.rm=T),max(.dat,na.rm=T)))))
      points(as.numeric(rownames(.dat)), -.dat[, excl], type='h', 
             col={if(excl %in% "Up") RED else BLUE}, lwd=1)
      mtext(sp, side=3, cex=0.75, line=0)
      abline(h=0, lwd=1)
      .len.std <- .len[.len$MMED=='None', ]
      .len.mmed <- .len[.len$MMED==excl, ]
      .x <- rep(.len.std$Length, .len.std$Number)
      .y <- rep(.len.mmed$Length, .len.mmed$Number)
      boxplot(list(.x,.y), add=TRUE, horizontal=TRUE, boxwex=0.2*.maxN, 
              at=0.25*c(.maxN,-.maxN), lwd=2, col=rgb(1,1,1,0.5), 
              outline=FALSE, axes=FALSE)
      if ((length(.x) > 10) & (length(.y) > 10)) {
        print(wilcox.test(.x, .y, alt="two.sided"))
        print(ks.test(.x, .y))
      } else {
        cat('\n Insufficient data for KS test \n')
      } # if (length...)
    } else {
      cat('\n\tNO LENGTH DATA FOR ', sp, '\n')
    } # if (sp %in% ...)
  } # for(sp)
  mtext('Size (mm)', side=1, outer=TRUE, cex=1.5)
  mtext('Number Caught', side=2, outer=TRUE, cex=1.5)
} # for (excl)
