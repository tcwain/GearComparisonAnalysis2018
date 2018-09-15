#!/usr/bin/R

#Plume Project -- Analysis of MMED trawl data
#  Data summary tables

# !!! Data has been read by DataSummaryTables.R !!!

## ---- TabCountsByHaulSpec

# Total number by Haul (rows) and Species
MMEDcnt <- with (MMEDdata, tapply(Number, list(Haul, SpecAge), FUN=sum, simplify=T))
MMEDcnt[is.na(MMEDcnt)] <- 0
MMEDcnt <- as.data.frame(MMEDcnt)
## print(summary(MMEDcnt)) ### DEBUG ###

## ---- TabHaulData

MMEDhauls <- data.frame(HaulID = rownames(MMEDcnt))
MMEDhauls$Station <- with(MMEDdata, tapply(as.character(Station), list(Haul), FUN=unique))
MMEDhauls$Date <- with(MMEDdata,  tapply(as.character(Date), list(Haul), FUN=unique))
MMEDhauls$Effort <- with(MMEDdata, tapply(Distance, list(Haul), FUN=unique))
MMEDhauls$Gear <- factor(with(MMEDdata, tapply(as.character(MMED), list(Haul), FUN=unique)),
                            levels=c('None','Up','Down'))
MMEDhauls$Block <- with(MMEDdata, tapply(as.character(Block), list(Haul), FUN=unique))

## ---- CPUECalcs

MMEDhauls$PlotTime <- match(MMEDhauls$Block, LETTERS[1:13]) - 1 +
  c((1:10)/11, (1:10)/11, (1:4)/5, (1:4)/5, (1:10)/11, (1:4)/5, (1:4)/5, (1:8)/9,
    (1:8)/9, (1:8)/9, (1:4)/5, (1:4)/5, (1:8)/9)
## print(summary(MMEDhauls)) ### DEBUG ###
MMEDcpue=sweep(MMEDcnt, 1, MMEDhauls$Effort, '/')
## print(summary(MMEDcpue)) ### DEBUG ###

## ---- PlotSetup

bw <- FALSE          #Flag for black-and-white figures
BLACK <- 'black'
BLUE <- if(bw) 'black' else 'blue'   #color code for blue
RED <-  if(bw) 'black' else 'red'    #color code for red
lndscp <- FALSE                      #flag for landscape figures
plcol <- if(lndscp) 3 else 2         #number of columns for multi-plots

## ---- CPUEPlotFnc

cpue.plot <- function(t, y, dot.col=1, log.zero=FALSE, ...) {
  minpos <- 0
  if (log.zero) {
    minpos <- min(y[y>0]) #minimum positive value
    if(any(y<=0)) {
      y[y<=0] <- minpos/2 # recode
    } # if(any...
  } # if(log.zero)
  plot(t, y, log=ifelse(log.zero,"y",""), type='p',
       axes=F, ...)
  # label blocks at midpoint:
  blocks <- seq(round(min(t)), round(max(t)))
  axis(side=1, at=blocks, labels=NA) # Ticks at day boundaries
  abline(v=blocks, col='blue')
  axis(side=1, at=blocks[-1]-0.5,
       labels=sort(unique(as.character(MMEDhauls$Block))),
       tick=FALSE, cex.axis=0.8)
  if (log.zero && any(y<minpos)) {
    tck<- axisTicks(range(y),log=FALSE) #default tick locations
    axis(side=2, at=c(minpos/2,tck), labels=c(0,tck))
  } else {
    axis(side=2)
  } # if (log.zero && ...
  box()
} # cpue.plot()

## ---- DoCPUEPlots

par(mfrow=c(ceiling(length(sel.spec)/plcol), plcol), omi=c(0.5,0.5,0,0), mar=c(3,2,2,1))
for (sp in sel.spec) {
  cpue.plot(MMEDhauls$PlotTime, MMEDcpue[[sp]], log.zero=TRUE,
            col=c(BLACK,RED,BLUE)[as.numeric(MMEDhauls$Gear)],
            pch=1+4*(MMEDhauls$Gear!='None'), xlab="", ylab="", main=sp)
} # for (sp)
mtext('Survey Block', side=1, outer=TRUE, cex=1.5)
mtext('CPUE', side=2, outer=TRUE, cex=1.5)
