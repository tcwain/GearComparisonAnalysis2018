#!/usr/bin/R

#Plume Project -- Analysis of MMED trawl data
#  Compute catch-ratio statistics by various methods

## ---- StatsSetup

assign('gears', c('None', 'Up', 'Down'), 1) 
options(scipen=3) # Favor non-scientific notation

## ---- NegBinAnoDevFnc

GLMAnoDevEst <- function (dat, qprobs=0.5, nb=TRUE, init.theta=1, 
                          diag.plt=FALSE, plt.lab='') {
  rslt <- matrix(NA, ncol=length(gears)-1, nrow=1+length(qprobs), 
                 dimnames=list(c("Mean",as.character(qprobs)), gears[2:length(gears)]))
  library(MASS)
  fit.data <- data.frame(Count=as.vector(dat[,1]), 
                         Date=as.factor(dat$Date),
                         Station=as.factor(dat$Station),
                         Block=as.factor(dat$Block), 
                         Gear=as.factor(dat$Gear),
                         Offset=as.vector(dat$Effort))
  if (nb) { # Negative Binomial fit
    cat('\n\tInitial fit to estimate theta\n')
    fit.init <- try(glm.nb(Count ~ Block + Gear + offset(log(Offset)), 
                           data=fit.data, init.theta=init.theta,
                           control=list(epsilon=1e-03, maxit=500, trace=0)))
    if (inherits(fit.init, 'try-error')) {
      print(fit.init)
      warning("glm.nb failed to estimate theta; using default value")
      theta.init=init.theta  # initial estimate for herring from theta.ml
    } else {
      cat('Estimated theta: ', fit.init$theta, ', SE: ', fit.init$SE.theta, '\n')
      if(is.finite(fit.init$SE.theta)) {
        theta.init <- fit.init$theta
      } else { 
        warning("glm.nb failed to estimate theta; using default value")
        theta.init <- init.theta
      }
    } # if 'try-error'
    cat('\n\tFinal fit with theta = ', theta.init, '\n')
    fit.fin <- glm(Count ~ Block + Gear + offset(log(Offset)),
                   data=fit.data, family=negative.binomial(theta.init),
                   control=list(epsilon=1e-08, maxit=500, trace=FALSE))
  } else { # Poisson fit
    fit.fin <- glm(Count ~ Block + Gear + offset(log(Offset)),
                   data=fit.data, family=poisson,
                   control=list(epsilon=1e-08, maxit=500, trace=FALSE))
  } # if (nb)
  cat("\nFIT STATISTICS:")
  # print(summary(fit.fin))
  fit.anova <- anova(fit.fin, test="Chisq")
  print(fit.anova)
  geareffects <- paste('Gear',gears[2:length(gears)], sep='')
  .lmn <- summary(fit.fin)$coefficients[geareffects, "Estimate"]
  .lsd <- summary(fit.fin)$coefficients[geareffects, "Std. Error"]
  .df <- fit.fin$df.residual
  for (g in 1:length(.lmn)) {
    if (.lsd[g] > 1000) {  # Estimate blew up, just use the mean value
      .mn <- exp(.lmn[g])
      .qnt <- rep(NA, length(qprobs))
    } else {
      .mn <- exp(.lmn[g] + .lsd[g]^2 / 2)
      .qnt <- exp(qt(qprobs, .df)*.lsd[g]+.lmn[g])
    } # if (.lsd[g])
    rslt[ , g] <- c(.mn, .qnt)
  } # for (g)
  predCatch <- predict(fit.fin, type="response", se.fit=TRUE)
  if (diag.plt) {
    rs <- resid(fit.fin, type="deviance")
    op <- par(omi=c(0,0,0.25,0), mfrow=c(1,2), mar=c(4,4,1,1))
    plot(predCatch$fit, rs, xlab="Prediction", ylab="Deviance Resids")
    qqnorm(rs, ylab="Deviance Resids")
    qqline(rs)
    mtext(paste(plt.lab, ifelse(nb, "Neg. Binomial", "Poisson"), sep=' - '), 
          side=3, outer=TRUE)
    par(op)
  } # if (diat.plt)
  return(list(Smry=rslt, Pred=predCatch))
} # GLMAnoDevEst()

## ---- RunStats

stat.sum <- list()  # structure for storing summary results
GLM.pred <- list() # structure for storing GLM predictions
for (sp in sel.spec) {
    cat('\n*******************  ', sp, '  ********************\n')
  .sumtbl <- array(NA, dim=c(2, 6, length(gears)-1),
                   dimnames=list(c('GLM.Po', 'GLM.nb'), 
                                 c('Mean', 'q0.05', 'q0.25', 'Median', 'q0.75', 'q0.95'), 
                                 gears[2:length(gears)]))
  .qprobs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  cat("\n*** METHOD 1: GLM ANODEV, Poisson ***\n")
  .est <- GLMAnoDevEst(cbind(MMEDcnt[[sp]], MMEDhauls), .qprobs, 
                       nb=FALSE, diag.plt=TRUE, plt.lab=sp)
  print(.est$Smry)
  .sumtbl['GLM.Po', , ] <- .est$Smry
  
  cat("\n*** METHOD 2: GLM ANODEV, negative binomial ***\n")
  init.theta <- theta.ml(y=MMEDcnt[[sp]], mu=.est$Pred$fit, 
                         n=length(MMEDcnt[[sp]]), limit=100, trace=FALSE)
  cat('Initial Theta: ', init.theta, '\n')
  .est <- GLMAnoDevEst(cbind(MMEDcnt[[sp]], MMEDhauls), .qprobs, 
                       nb=TRUE, init.theta=init.theta, diag.plt=TRUE, plt.lab=sp)
  print(.est$Smry)
  .sumtbl['GLM.nb', , ] <- .est$Smry
  
  # Add species to summary lists
  GLM.pred[[sp]] <- .est$Pred
  .sumtbl <- round(.sumtbl, 3)  # round statistical results
  stat.sum[[sp]] <- .sumtbl
} # for(sp)

## ---- CRSummaryTable

print(stat.sum)

## ---- CRSummaryFig

par(mfrow=c(ceiling(length(stat.sum)/plcol), plcol), omi=c(0.5,0.5,0,0), mar=c(3,4,2,1))
for (sp in names(stat.sum)) {
  .sumtbl <- stat.sum[[sp]]
  .sumtbl[.sumtbl==Inf] <- 99 # recode infinite values as +99
  .minx <- 5e-3
  .sumtbl[.sumtbl<.minx] <- .minx # recode zeros as small pos. value (for log scale plots)
  .nstats <- dim(.sumtbl)[1]
  .ngears <- dim(.sumtbl)[3]
  .ny <- .ngears*.nstats # number of elements along y axis.
  .gclr <- rep(c('red','blue','green3')[1:.ngears], .nstats) # gear color codes
  .xmax <- ceiling(max(.sumtbl[ , 'q0.75', ], na.rm=TRUE)) # make sure the quartiles are covered
  .xmax <- max(.xmax, 2) #make sure upper bound is above 1
  .xmax <- min(.xmax, 10) #truncate high values so plot is readable
  .xlim <- c(.minx, 1/.minx)
  plot(t(.sumtbl[ ,'Mean', ]), log='x', 1:.ny, col=.gclr, pch=3, cex=0.75, axes=F,
       xlim=.xlim, ylim=c(0.5, .ny+0.5),
       main=sp, xlab='', ylab='')
  box()
  abline(v=1, lwd=2, col='black')
  axis(side=1, cex.axis=1.2)
  axis(side=2, at=seq(1,.ny,.ngears)+1/.ngears, labels=dimnames(.sumtbl)[[1]], las=2, cex.axis=1.2)
  points(t(.sumtbl[ ,'Median', ]), 1:.ny, col=.gclr, pch=5, cex=0.75)
  rect(t(.sumtbl[ ,'q0.25', ]), (1:.ny)-0.35, t(.sumtbl[ ,'q0.75', ]), (1:.ny)+0.35, border=.gclr)
  segments(t(.sumtbl[ ,'q0.05', ]), 1:.ny, t(.sumtbl[ ,'q0.25', ]), 1:.ny, col=.gclr, lwd=1)
  segments(t(.sumtbl[ ,'q0.75', ]), 1:.ny, t(.sumtbl[ ,'q0.95', ]), 1:.ny, col=.gclr, lwd=1)
} # for (sp)
mtext('Catch Ratio', outer=T, side=1, line=1, cex=1.5)
mtext('Method', outer=T, side=2, line=1, cex=1.5)

## ---- GLMSummaryFig

par(mfrow=c(1,1), omi=c(0.5,0.5,0,0), mar=c(4,4,1,1))
.sumtbl <- simplify2array(stat.sum)[ "GLM.nb", , , ] # array: probs x gear x species
.sumtbl[.sumtbl==Inf] <- 99 # recode infinite values as +99
.miny <- 1e-2
.sumtbl[.sumtbl<.minx] <- .miny # recode zeros as small pos. value (for log scale plots)
.ngears <- dim(.sumtbl)[2]
.nspecs <- dim(.sumtbl)[3]
.nx <- .ngears*.nspecs      # number of elements along y axis.
.gbox <- rep(c('red','blue','green3')[1:.ngears], .nspecs) # box colors
.gpnt <- rep(c('red','blue','grey80')[1:.ngears], .nspecs) # point colors
.gfill <- rep(c(NA,'grey80','green3')[1:.ngears], .nspecs) # fill colors
.maxy <- ceiling(max(.sumtbl[ 'q0.95', , ], na.rm=TRUE))
.ylim <- c(.miny, .maxy)
plot(1:.nx, .sumtbl['Mean', , ], type='n', log='y', xaxs='i', axes=F,
     ylim=.ylim, xlim=c(0.5, .nx+0.5),
     xlab='', ylab='Catch Ratio')
box()
abline(h=1, lwd=2, col='black')
# label "fake zero" as zero:
axis(side=2, at=c(.miny,0.10,1,10,100), labels=c('0','0.1','1','10','100'), las=2)
.labs <- dimnames(.sumtbl)[[3]]
.labs <- sub(' ', '\n', .labs)
axis(side=1, at=seq(1, .nx, .ngears) + 1/.ngears, tick=FALSE, labels=.labs, 
     las=2, cex.axis=0.8)
abline(v=seq(.ngears+1, .nx, .ngears) - 1/.ngears, col="grey")
rect((1:.nx)-0.35, .sumtbl['q0.25', , ], (1:.nx)+0.35, .sumtbl['q0.75', , ], 
     border=.gbox, col=.gfill, lwd=2)
segments(1:.nx, .sumtbl['q0.05', , ], 1:.nx, .sumtbl['q0.25', , ], col=.gbox, lwd=2)
segments(1:.nx, .sumtbl['q0.75', , ], 1:.nx, .sumtbl['q0.95', , ], col=.gbox, lwd=2)
points(1:.nx, .sumtbl['Median', , ], col=.gpnt, pch=5, cex=0.75)
points(1:.nx, .sumtbl['Mean', , ], col=.gpnt, pch=3, cex=0.75)
text(1, .miny, "Upward", col="red", adj=c(0,0.5), srt=90)
text(2, .miny, "Downward", col="blue", adj=c(0,0.5), srt=90)
