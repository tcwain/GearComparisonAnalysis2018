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
lenData <- lenData[lenData$Species %in% len.spec, ]
# Add subsample ratio:
# Total number by Haul (rows) and Species (cols)
cnt <- with (lenData, tapply(Number, list(Haul, Species), FUN=sum, simplify=T))
cnt[is.na(cnt)] <- 0
cnt <- as.data.frame(cnt)
# Total measured by Haul & species
meas <- with(lenData, tapply(!is.na(Length), list(Haul, Species), FUN=sum, simplify=T))
meas[is.na(meas)] <- 0
meas <- as.data.frame(meas)
# Subsampling ratio by Haul & Species:
ssr <- meas / cnt
cat('\nSubsampling Ratios:\n')
print(summary(ssr))
# Adjusted Numbers (expanded by ssr)
lenData <- lenData[!is.na(lenData$Length), ] #remove non-measured counts
lenData$AdjNum <- lenData$Number / unlist(apply(lenData[c("Haul","Species")], 1,
                  function(x){ssr[x["Haul"], x["Species"]]}))
# Length bin size (mm)
binsize <- 5
lenData$LenBin <- binsize * round(lenData$Length/binsize)
# Remove size outliers for anchovy & water jelly (likely data errors)
lenData <- lenData[!((lenData$Species == "NORTHERN ANCHOVY") &
                        (lenData$Length < 100)), ]
lenData <- lenData[!((lenData$Species =="WATER JELLY") &
                        (lenData$Length > 150)), ]
cat('\nSummary of Length Data:\n')
print(summary(lenData))

## ---- SizeSelMethods

boot_GLM3P <- function(sdat, nrep=10, binsz=5, L.pr=NULL) {
  fit.model <- function(sdat) {
    NumTotL <- with(sdat, tapply(AdjNum, list(LenBin, MMED), sum,
                                 na.rm=TRUE, default=0))
    EffTotL <- with(sdat, tapply(Distance, list(LenBin, MMED),sum,
                                 na.rm=TRUE, default=0))
    cpue <- NumTotL/EffTotL
    cpue[is.na(cpue)] <- 0
    STD <- match("None", colnames(cpue))
    TST <- match("Up", colnames(cpue))
    if (is.na(TST)) TST <- match("Down", colnames(cpue))
    std <- cpue[ , STD]
    tst <- cpue[ , TST]
    p.L12 <- std / (std + tst)
    # Binomial weights based on number measured in both gears:
    Nmeas <- with(sdat, tapply(Number, list(LenBin, MMED), sum,
                               na.rm=TRUE, default=0))
    wts <- Nmeas[ , STD] + Nmeas[ , TST]
    L <- as.numeric(names(p.L12))
    old.opt <- options(warn = -1) # suppress warnings about non-integer values
    fit.glm <- glm(p.L12 ~ L + I(L^2) + I(L^3), family=binomial, weights=wts)
    options(old.opt)
    return(fit.glm)
  } # fit.model()

  # Fit the model to the original (full) dataset:
  fit.full <- fit.model(sdat)
  # Predictions of full model, with rough SE's
  if (is.null(L.pr)) L.pr <- seq(min(sdat$LenBin), max(sdat$LenBin), 5)
  pred.full <- predict(fit.full, newdata=data.frame(L=L.pr, wts=1.0),
                       type="response")
  names(pred.full) <- L.pr
  # Bootstrap predictions:
  bs <- matrix(NA, nrow=length(L.pr), ncol=nrep, dimnames=list(L.pr, NULL))
  rep <- 0
  while (rep < nrep) {
    hauls <- unique(sdat$Haul)
    hauls.samp <- sample(hauls, length(hauls), replace=TRUE)
    .data <- data.frame()
    for (h in hauls.samp) {
      .hdata <- sdat[sdat$Haul == h, ]
      ssr <- with(.hdata, sum(Number)/sum(AdjNum)) # subsample rate
      L.ex <- with(.hdata, rep(LenBin, Number))  #expand Number to indiv. lengths
      if(length(L.ex) > 1) {
        L.smp <- sample(L.ex, length(L.ex), replace=TRUE)  #resample lengths
      } else {
        L.smp <- L.ex  # sample() doesn't work for length 1 vector
      }
      new.freq <- as.data.frame(table(L.smp))
      .ndata <- data.frame(MMED=unique(.hdata$MMED),
                           Haul=unique(.hdata$Haul),
                           Distance=mean(.hdata$Distance),
                           LenBin=as.numeric(levels(new.freq$L.smp)),
                           Number=new.freq$Freq,
                           AdjNum=new.freq$Freq/ssr)
      .data <- rbind(.data, .ndata)
    } # for (h)
    names(.data) <- names(sdat)
    fit.rep <- fit.model(.data)
    rep <- rep+1
    bs[ , rep] <- predict(fit.rep, newdata=data.frame(L=L.pr, wts=1.0),
                          type="response")
  } # for (rep)
  rownames(bs) <- L.pr
  bs.mn <- apply(bs, 1, mean, na.rm=FALSE)
  bs.q <- t(apply(bs, 1, quantile, probs=c(0,0.05,0.25,0.50,0.75,0.95,1),
                  na.rm=FALSE))
  return(list(gam=fit.full, pred=pred.full, boot=bs,
              boot.sum=data.frame(mean=bs.mn, q=bs.q)))
} # boot_GLM3P()

## ---- SizeSelAnal

# Set number of bootstrap replicates:
nbsr <- 50   ### TESTING ###
## nbsr <- 1000 ### PRODUCTION ###

for (excl in c("Up","Down")) {
  cat('\n**************** Excluder: ', excl, ' *****************\n')
  lD <- if(excl %in% "Up") {
    lenData[lenData$Cruise %in% c(41,43,50), ]
  } else {
    lenData[lenData$Cruise %in% 53, ]
  }
  lD$Species <- factor(as.character(lD$Species))
  .tab <- with(lD, tapply(AdjNum, list(Species), sum, na.rm=T))
  cat('Total Adjusted Catch:\n')
  print(.tab)
  lf.sel.spec <- names(.tab)[.tab >= 100]
  lenFreq <- with(lD, tapply(AdjNum, list(LenBin, MMED, Species), sum, na.rm=TRUE))
  lenFreq[is.na(lenFreq)] <- 0
  .mfrow <- if(lndscp) c(3,3) else c(3,2)
  par(mfrow=.mfrow, omi=c(0.5,0.5,0,0.5), mar=c(3,3,2,3))
  for (sp in lf.sel.spec) {
    cat('\n****************', sp, '*****************\n')
    if (sp %in% dimnames(lenFreq)[[3]]) {
      .dat <- lenFreq[ , , sp]
      .maxN <- max(.dat,na.rm=T)
      .len <- lD[lD$Species %in% sp,
                 c("MMED", "Haul", "Distance", "LenBin", "Number", "AdjNum")]
      .maxL <- max(.len$LenBin, na.rm=T)
      .minL <- min(.len$LenBin, na.rm=T)
      .len.std <- .len[.len$MMED=='None', ]
      .len.mmed <- .len[.len$MMED==excl, ]
      .x <- rep(.len.std$LenBin, .len.std$Number)
      .y <- rep(.len.mmed$LenBin, .len.mmed$Number)
      # Run analysis only if > 40 measurements in each gear:
      if ((length(.x) > 40) & (length(.y) > 40)) {
        # Wilcox & KS test for overall difference in size-frequencies
        print(wilcox.test(.x, .y, alt="two.sided"))
        print(ks.test(.x, .y))
        # GLM fit of Catch Ratio to size:
        mod.fit <- boot_GLM3P(sdat=.len, nrep=nbsr, binsz=binsize)
        cat("\tSummary of GAM fit: \n")
        print(summary(mod.fit$gam))
        print(anova(mod.fit$gam, test="Chisq"))
        cat("\tSummary of bootstrap fits: \n")
        print(summary(mod.fit$boot.sum))
        p.pred <- mod.fit$pred
        L.pred <- as.numeric(names(p.pred))
        # Convert probability to Catch Ratio:
        CR.obs <- 1/mod.fit$gam$model$p.L12 - 1
        CR.pred <- 1/p.pred - 1
        CR.boot <- 1/mod.fit$boot.sum - 1
        CR.boot[CR.boot > 1000] <- 1000 #recode infinite values
        CR.boot[CR.boot < 1/1000] <- 1/1000 #recode zero values
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
        plot(L.pred, CR.pred, log="y", type='l', lwd=3, axes=FALSE, bty="n",
             ylim=c(1/50, 50), xlab="", ylab="")
        axis(side=4, at=c(0.02, 0.2, 1, 5, 50),
             labels=c("0.02", "0.2", "1", "5", "50"))
        lines(L.pred, CR.boot[ , "q.50."], lty=1, col='gray50', lwd=2) # bs median
        lines(L.pred, CR.boot[ , "q.5."], lty=2, col='gray50', lwd=2) # bs lower 5%
        lines(L.pred, CR.boot[ , "q.95."], lty=2, col='gray50', lwd=2) # bs upper 5%
      } else {
        cat('\n Insufficient data \n')
      } # if (length...)
    } else {
      cat('\n\tNO LENGTH DATA FOR ', sp, '\n')
    } # if (sp %in% ...)
  } # for(sp)
  mtext('Size (mm)', side=1, outer=TRUE, at=c(0.5,0.5), cex=1.5)
  mtext('Observed Catch Per Km', side=2, outer=TRUE, at=c(0.5,0.5), cex=1.5)
  mtext('Estimated Catch Ratio', side=4, outer=TRUE, at=c(0.5,0.5), cex=1.5)
} # for (excl)
