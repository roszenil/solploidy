### Compare breeding system data from 2014 compilation with the current file.
### (Don't worry about uncertain/polymorphic species for now.)

library(ape)
phy <- read.nexus("sarkinen1.nex")  # just the first of the 100 trees

#--------------------------------------------------
# Data files from 2014
#--------------------------------------------------

dat1 <- read.csv("statesICDP.csv", as.is=T)
dat1 <- subset(dat1, Species %in% phy$tip.label)

dat1 <- subset(dat1, !(ID==CD & CD==CP)) # remove 0.33 0.33 0.33
dat1$SI <- dat1$ID
dat1$SC <- dat1$CD + dat1$CP
dat1 <- subset(dat1, SI != SC) # remove 0.5 0.5
dat1[which(!(dat1$SI %in% c(0, 1))),] # Solanum_brevicaule Solanum_candolleanum Solanum_chacoense
dat1 <- dat1[-which(!(dat1$SI %in% c(0, 1))),] # drop further uncertainty
table(dat1$SI)
table(dat1$SC)
# 331 SC, 143 SI
subset(dat1, Species %in% c("Capsicum_schottianum", "Capsicum_scolnikianum", "Lycium_californicum")) # for comparison with dat2, below

# fully cleaned
dat1 <- data.frame(name=dat1$Species, bs1=c("SC", "SI")[dat1$SI+1])

#--------------------------------------------------
# Current data file
#--------------------------------------------------

dat2 <- read.csv("solanaceaerecords.csv", as.is=T)
dat2 <- subset(dat2, fullname %in% phy$tip.label)

dat2 <- unique(dat2[,c("fullname", "selfincomp")])
dat2 <- subset(dat2, selfincomp != "")
dat2[which(duplicated(dat2$fullname)),]
subset(dat2, fullname %in% c("Capsicum_schottianum", "Lycium_californicum")) # the Capsicum was coded as SC only in dat1; the Lycium is truly polymorphic

# fully cleaned
dat2 <- subset(dat2, !(fullname %in% c("Capsicum_schottianum", "Capsicum_scolnikianum", "Lycium_californicum")))
names(dat2) <- c("name", "bs2")

#--------------------------------------------------
# Compare
#--------------------------------------------------

dat <- merge(dat1, dat2, all=T)

# species that disagree between the two data sets:
subset(dat, bs1 != bs2)
#                      name bs1 bs2
#          Capsicum_eximium  SI  SC
#        Capsicum_pubescens  SI  SC
# Nicotiana_plumbaginifolia  SI  SC
#         Petunia_axillaris  SI  SC
#   Solanum_agrimoniifolium  SC  SI
#          Solanum_chilense  SI  SC
#          Solanum_obliquum  SI  SC
#          Solanum_palustre  SI  SC
#      Solanum_stipuloideum  SI  SC

# species that had a known breeding system before but don't now:
subset(dat, is.na(bs2))
# 108 rows
