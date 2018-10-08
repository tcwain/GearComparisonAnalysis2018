#!/usr/bin/R

#Plume Project -- Analysis of MMED trawl data
#  Data summary tables

## ---- ReadData

#data.dir <- '.'  #Directory where the data files resides
MMEDdata <- read.csv(system.file('extdata','AllSppCounts&Lengths.csv',
                                 package = "GearComparisonAnalysis2018"))
# print(summary(MMEDdata))

## ---- FixData

# Haul ID is last three characters of station code:
MMEDdata$Station <- as.character(MMEDdata$Station) #remove factor levels
MMEDdata$Haul <- with(MMEDdata, substr(Station, nchar(Station)-2, nchar(Station)))
# Excluder code (Y or N) is last character of station code:
MMEDdata$Excluder <- with(MMEDdata, substr(Station, nchar(Station), nchar(Station)))
# Filter out "experimental" hauls:
MMEDdata <- MMEDdata[-grep('*X$', MMEDdata$Haul), ]
# Recode MMED types to Standard short labels
### NOTE: work around database error where some records with Excluder code "N"
###   have MMED code "Yes; Up"
MMEDdata$MMED <- as.character(MMEDdata$MMED) #remove factor levels
MMEDdata$MMED[MMEDdata$Excluder == "N"] <- "None"
MMEDdata$MMED[MMEDdata$Excluder == "Y" & grepl("Up", MMEDdata$MMED)] <- "Up"
MMEDdata$MMED[MMEDdata$Excluder == "Y" & grepl("Down", MMEDdata$MMED)] <- "Down"
# Fix a few species names, add age-classes for Chinook & coho
MMEDdata$Species <- toupper(as.character(MMEDdata$Species))    #all upper case
MMEDdata$Species[MMEDdata$Species %in% 'CALIFORNIA MARKET SQUID'] <- 'MARKET SQUID'
MMEDdata$Species[grepl('SMELT',MMEDdata$Species)] <- "SMELT SPP."
MMEDdata$SpecAge <- as.character(MMEDdata$Species)
.index <- MMEDdata$SpecAge %in% 'CHINOOK SALMON'
.newnames <- paste('CHINOOK', MMEDdata$AgeGp[.index])
MMEDdata$SpecAge[.index] <- .newnames
.index <- MMEDdata$SpecAge %in% 'COHO SALMON'
.newnames <- paste('COHO', MMEDdata$AgeGp[.index])
MMEDdata$SpecAge[.index] <- .newnames
MMEDdata$SpecAge[MMEDdata$SpecAge %in%
    c('CHINOOK subadult/adult','CHINOOK mixed age juvenile')] <- 'CHINOOK subadult'
MMEDdata$SpecAge[MMEDdata$SpecAge %in% 'COHO subadult/adult'] <- 'COHO subadult'
MMEDdata$SpecAge <- factor(as.character(MMEDdata$SpecAge))
# print(summary(MMEDdata))

## ---- DefineBlocks

# sequential along the sorted haul numbers,
blocks <- as.factor(c(rep('A',10), rep('B',10), rep('C', 4), rep('D', 4),
                      rep('E',10), rep('F', 4), rep('G', 4), rep('H', 8),
                      rep('I', 8), rep('J', 8), rep('K', 8), rep('L', 8)))
sortHauls <- sort(unique(MMEDdata$Haul))
MMEDdata$Block <- blocks[match(MMEDdata$Haul, sortHauls)]
# print(with(MMEDdata, t(tapply(as.character(Block), list(Haul,MMED), unique))))

## ---- TotCatchBySpecies

cat('\n*** Total Catch By Species and Gear Type ***\n')
tab1 <- with(MMEDdata, tapply(Number, list(SpecAge, MMED), sum, na.rm=T))
tab1[is.na(tab1)] <- 0  #Missing values are actually zero counts
tab1 <- cbind(tab1, Total=apply(tab1, 1, sum))
print(tab1[ , c('Down', 'Up', 'None', 'Total')])

## ---- FrequencyBySpecies

cat('\n*** Frequency of Catch By Species and Gear Type ***\n')
.tmp.all <- with(MMEDdata, table(SpecAge, Haul))
.tmp.std <- with(MMEDdata[MMEDdata$MMED=='None', ], table(SpecAge, Haul))
.tmp.up <- with(MMEDdata[MMEDdata$MMED=='Up', ], table(SpecAge, Haul))
.tmp.dwn <- with(MMEDdata[MMEDdata$MMED=='Down', ], table(SpecAge, Haul))
tab2 <- cbind(apply(.tmp.dwn>0, 1, sum),    # Num. occurrences in STD
              apply(.tmp.up>0, 1, sum),     # . . .
              apply(.tmp.std>0, 1, sum),     # . . .
              apply(.tmp.all>0, 1, sum))    # Total Num. occurrences
colnames(tab2) <- c('Down','Up','None','Total')
print(tab2)

## ---- SpeciesByCruise

cat('\n*** Total Catch By Species and Cruise ***\n')
tab3 <- with(MMEDdata, tapply(Number, list(SpecAge, Cruise), sum, na.rm=T))
tab3[is.na(tab3)] <- 0  #Missing values are actually zero counts
ngt0 <- apply(tab3>0, 1, sum)
ngt1 <- apply(tab3>1, 1, sum)
tab3 <- cbind(tab3, Ngt0=ngt0, Ngt1=ngt1)
print(tab3)

## ---- NumMeasBySpecies

cat('\n*** Total Number Measured By Species and Gear Type ***\n')
tab4 <- with(MMEDdata[!is.na(MMEDdata$Length), ],
             tapply(Number, list(SpecAge, MMED), sum, na.rm=T))
tab4[is.na(tab4)] <- 0  #Missing values are actually zero counts
tab4 <- cbind(tab4, Total=apply(tab4, 1, sum))
print(tab4[ , c('Down', 'Up', 'None', 'Total')])

## ---- SSRBySpecies

cat('\n*** Average Subsampling Rate By Species and Gear Type ***\n')
tab5 <- round(tab4/tab1,2)
print(tab5[ , c('Down', 'Up', 'None', 'Total')])

## ---- SelSpecies

cat('\n*** Species Selected for Analysis ***\n')
sel.spec <- rownames(tab1)[tab1[ ,"Total"] >= 100]
sel.spec <- sel.spec[sel.spec %in% rownames(tab3[tab3[,"Ngt1"]>=3, ])]
print(sel.spec)
