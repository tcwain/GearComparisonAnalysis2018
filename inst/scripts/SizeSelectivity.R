#!/usr/bin/R

#Plume Project -- Size-selectivity of MMED trawl
#  comparisons Combined Data
# Sampling described in RC Dotson etal 2010 NOAA-TM-NMFS-SWFSC-455,
#  and NWFSC Plume Project cruise reports

# data.dir <- '.'  #Directory where MMED data resides

## ---- GetSizeData
# Restrict to species selected above, but no age groups for salmon
len.spec <- c("CHINOOK SALMON", "CHUM SALMON", "COHO SALMON",
             "MARKET SQUID", "NORTHERN ANCHOVY", "PACIFIC HERRING",
             "SEA NETTLE", "WATER JELLY")
lenData <- MMEDdata[ , c("Cruise","MMED","Species","Length",
                         "Number","Distance","Haul")]
# Add subsample ratio:
# Total number by Haul (rows) and Species
cnt <- with (lenData, tapply(Number, list(Haul, Species), FUN=sum, simplify=T))
cnt[is.na(cnt)] <- 0
cnt <- as.data.frame(cnt)
## print(summary(cnt)) ### DEBUG ###
# Total measured by Haul & species
meas <- with(lenData, tapply(!is.na(Length), list(Haul, Species), FUN=sum, simplify=T))
meas[is.na(meas)] <- 0
meas <- as.data.frame(meas)
## print(summary(meas)) ### DEBUG ###
# Subsampling ratio by Haul & Species:
ssr <- meas / cnt
##print(summary(ssr)) ### DEBUG ###
# Adjusted Numbers (expanded by ssr)
lenData <- lenData[!is.na(lenData$Length), ] #remove non-measured counts
lenData$AdjNum <- lenData$Number / unlist(apply(lenData[c("Haul","Species")], 1,
                  function(x){ssr[x["Haul"], x["Species"]]}))
lenData <- lenData[lenData$Species %in% len.spec, ]
lenData <- lenData[!is.na(lenData$Length), ] #Already accounted for in AdjNum
# Length bin size (mm)
binsize <- 5
lenData$LenBin <- binsize * round(lenData$Length/binsize)
print(summary(lenData)) ### DEBUG ###

## ---- SizeSummary

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
  lenFreq <- with(lD, tapply(Number, list(LenBin, MMED, Species), sum))
  .mfrow <- if(lndscp) c(3,3) else c(3,2)
  par(mfrow=.mfrow, omi=c(0.5,0.5,0,0), mar=c(3,3,2,1))
  for (sp in lf.sel.spec) {
    cat('****************', sp, '*****************\n')
    if (sp %in% dimnames(lenFreq)[[3]]) {
      .dat <- lenFreq[ , , sp]
      .maxN <- max(.dat,na.rm=T)
      .len <- lD[lD$Species %in% sp, ]
      .maxL <- max(.len$LenBin, na.rm=T)
      .minL <- min(.len$LenBin, na.rm=T)
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

## ---- SizeSelMethods

# First, a function for each of three methods:
# (klwg = "Kotwicki, Lauth, Williams, Goodman")
# klwg_GLMM <- function(sfdat) {
#   print("GLMM Not Yet Implemented")
# }

klwg_SCMM <- function(sfdat) {
  NumTotL <- with(sfdat, tapply(AdjNum, list(LenBin, MMED), sum, na.rm=TRUE))
##  print(dim(NumTotL))   ### DEBUG ###
  EffTotL <- with(sfdat, tapply(Distance, list(LenBin, MMED),sum, na.rm=TRUE))
##  print(dim(EffTotL))   ### DEBUG ###
  cpue <- NumTotL/EffTotL
  cpue[is.na(cpue)] <- 0
  ##print(summary(cpue))   ### DEBUG ###
  std <- cpue[ , 1]
  tst <- cpue[ , 2]
  # Binomial weights based on number measured in both gears:
  Nmeas <- with(sfdat, tapply(Number, list(LenBin, MMED), sum, na.rm=TRUE))
  wts <- Nmeas[ , 1] + Nmeas[ , 2]
  ##print(summary(wts))    ### DEBUG ###
  p.L12 <- std / (std + tst)
  ##p.L12 <- p.L12[!is.na(p.L12)]
  ##print(summary(p.L12))   ### DEBUG ###
  ##wts <- wts[!is.na(wts)]
  L <- as.numeric(names(p.L12))
  old.opt <- options(warn = -1) # suppress warnings about non-integer weights
  res <- mgcv::gam(p.L12 ~ s(L, bs="cr", k=5), family=binomial, weights=wts)
  options(old.opt)
  return(res)
}

# klwg_BetaR <- function(sfdat) {
#   print("BetaR Not Yet Implemented")
# }

## ---- SizeSelAnal
for (excl in c("Up","Down")) {
  cat('**************** Excluder: ', excl, ' *****************\n')
  lD <- if(excl %in% "Up") {
    lenData[lenData$Cruise %in% c(41,43,50), ]
  } else {
    lenData[lenData$Cruise %in% 53, ]
  }
  lD$Species <- factor(as.character(lD$Species))
  .tab <- with(lD, tapply(AdjNum, list(Species), sum, na.rm=T))
  print(.tab)
  lf.sel.spec <- names(.tab)[.tab >= 100]
  lenFreq <- with(lD, tapply(AdjNum, list(LenBin, MMED, Species), sum, na.rm=TRUE))
  lenFreq[is.na(lenFreq)] <- 0
  ##print(summary(lenFreq))    ### DEBUG ###
  .mfrow <- if(lndscp) c(3,3) else c(3,2)
  par(mfrow=.mfrow, omi=c(0.5,0.5,0,0.5), mar=c(3,3,2,3))
  ##for(sp in c("COHO SALMON")) {       ### TESTING ###
  for (sp in lf.sel.spec) {
    cat('****************', sp, '*****************\n')
    if (sp %in% dimnames(lenFreq)[[3]]) {
      .dat <- lenFreq[ , , sp]
      .maxN <- max(.dat,na.rm=T)
      ##print(.dat)     ### DEBUG ###
      .len <- lD[lD$Species %in% sp, ]
      .maxL <- max(.len$Length, na.rm=T)
      .minL <- min(.len$Length, na.rm=T)
      .len.std <- .len[.len$MMED=='None', ]
      ##print(.len.std)     ### DEBUG ###
      .len.mmed <- .len[.len$MMED==excl, ]
      .x <- sum(.len.std$Number)
      ##print(.x)     ### DEBUG ###
      .y <- sum(.len.mmed$Number)
      ##print(.y)     ### DEBUG ###
      if ((.x > 40) & (.y > 40)) {
      ##if ((length(!is.na(.x)) > 40) & (length(!is.na(.y)) > 40)) {
        ##glmm <- klwg_GLMM(sfdat=.len)
        scmm <- klwg_SCMM(sfdat=.len)
        print(summary(scmm))
        L.pred <- seq(.minL, .maxL, binsize)
        pred <- predict(scmm, newdata=data.frame(L=L.pred, wts=1.0),
                        type="response", se.fit=TRUE)
        p.pred <- pred$fit
        p.UL <- p.pred + 1.96*pred$se.fit
        p.LL <- p.pred - 1.96*pred$se.fit
        CR.obs <- 1/scmm$model$p.L12 - 1
        CR.pred <- 1/p.pred - 1
        CR.UL <- 1/p.LL - 1
        CR.UL[CR.UL > 1000] <- 1000 #recode infinite values
        CR.LL <- 1/p.UL - 1
        CR.LL[CR.LL < 0.001] <- 0.001 # recode zeros for log plot
        plot(as.numeric(rownames(.dat)), -.dat[, 'None'], type='h', col=BLACK, lwd=1,
             xlim=c(.minL,.maxL), ylim=c(-.maxN,.maxN), axes=F,
             xlab='', ylab='')
        box()
        axis(side=1, lwd=0, lwd.ticks=1)
        axis(side=2, at=pretty(c(-max(.dat,na.rm=T),max(.dat,na.rm=T))),
             labels=abs(pretty(c(-max(.dat,na.rm=T),max(.dat,na.rm=T)))))
        points(as.numeric(rownames(.dat)), .dat[, excl], type='h',
               col={if(excl %in% "Up") RED else BLUE}, lwd=1)
        mtext(sp, side=3, cex=0.75, line=0)
        abline(h=0, lwd=1)
        # add Catch Ratio plot on right axis
        par(new=TRUE)
        plot(L.pred, CR.pred, log="y", type='l', lwd=2, axes=FALSE, bty="n",
             ylim=c(1/50, 50), xlab="", ylab="")
        axis(side=4, at=c(0.02, 0.2, 1, 5, 50),
             labels=c("0.02", "0.2", "1", "5", "50"))
        # abline(h=1.0, col="red")
        lines(L.pred, CR.UL, lty=2, lwd=2)
        lines(L.pred, CR.LL, lty=2, lwd=2)
        ##beta <- klwg_BetaR(sfdat=.len)
      } else {
        cat('\n Insufficient data \n')
      } # if (length...)
    } else {
      cat('\n\tNO LENGTH DATA FOR ', sp, '\n')
    } # if (sp %in% ...)
  } # for(sp)
  mtext('Size (mm)', side=1, outer=TRUE, at=c(0.5,0.5), cex=1.5)
  mtext('Observed Catch per km', side=2, outer=TRUE, at=c(0.5,0.5), cex=1.5)
  mtext('Estimated Catch Ratio', side=4, outer=TRUE, at=c(0.5,0.5), cex=1.5)
} # for (excl)

# Implementation of "double-bootstrap"
