setwd('~/Dropbox/solploidy/basicdata')
library("ape")
library("phytools")
library("treeplyr")

## Time-tree from Sarkinen et al. 2013 solaneaceae done with Beast
solanum.sarkinen<-read.tree(file="solanumbeast.tre")
## Names for the solaneaceae tree, the tree above was listed with numbers which was less than ideal
solanum.names<-read.csv("speciesnames2.csv",header=FALSE,sep=",",stringsAsFactors=FALSE)
## Putting the names on the tip labels
fixingtiplabels<-as.numeric(solanum.sarkinen$tip.label)
solanum.sarkinen$tip.label<-solanum.names[fixingtiplabels,2] 
# I double checked with their publication turns out to be exact tree shown in figure 12862_2013_2443_MOESM2_ESM.tiff

#Tips to drop/ repeated species list from EEG
remove.names<-read.csv("tipstodrop.csv",header=FALSE, sep=",",stringsAsFactors=FALSE)
length(remove.names[,1])#55 to drop
solanaceae.correcttips<-drop.tip(solanum.sarkinen,remove.names[,1])
# This corrected tree has 1023 tips, it dropped 53 out 55
write.tree(solanaceae.correcttips, file="solcorrecttip.tre")

##Plot
plotTree(solanaceae.correcttips,part=0.98,fsize=0.05, lwd=0.9)
# The tree is approximately 30my
h<-max(nodeHeights(solanaceae.correcttips))
#axis(1,labels=c(136.9,68.4,0), pos=-0.03*h,at=seq(0,h,by=h/2),lwd=0.2,cex=0.2)
