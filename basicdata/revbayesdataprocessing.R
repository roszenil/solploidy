# Processing for Revbayes files
library("ape")
library("phytools")
library("treeplyr")
solanaceae.tree<-read.tree("~/Dropbox/solploidy/basicdata/solcorrecttip.tre")
threestate.matrix <- read.csv("~/Dropbox/solploidy/basicdata/newstatesICDP.csv", as.is=T)
threestate.matrix<-threestate.matrix[,-1]

matched<-make.treedata(solanaceae.tree,threestate.matrix) #638 species!!!!
write.nexus(matched$phy,file="~/Dropbox/solploidy/basicdata/fullmatchtree.nex")

diploidsi<-which((matched$dat[,1]==1 & matched$dat[,2]==0 & matched$dat[,3]==0)==TRUE)# 144
diploidsc<-which((matched$dat[,1]==0 & matched$dat[,2]==1 & matched$dat[,3]==0)==TRUE)#157
polyploids<-which((matched$dat[,1]==0 & matched$dat[,2]==0 & matched$dat[,3]==1)==TRUE) #108

diploidsna<-which((matched$dat[,1]==1 & matched$dat[,2]==1 & matched$dat[,3]==0)==TRUE) #170
sc<-which((matched$dat[,1]==0 & matched$dat[,2]==1 & matched$dat[,3]==1)==TRUE)#70
idcp<-which(matched$dat[,1]==1 & matched$dat[,2]==0 & matched$dat[,3]==1)#2

## Data for threestate analyses
threestate<- rep(0,dim(matched$dat)[1])
threestate[diploidsc]="0"
threestate[polyploids]="1"
threestate[diploidsi]="2"
threestate[sc]="(0 1)"
threestate[diploidsna]<-"(0 2)"
threestate[idcp]<-"(1 2)"
threestate<-as.data.frame(threestate)
row.names(threestate)<-matched$phy$tip.label

write.table(threestate,file="~/Dropbox/solploidy/basicdata/threestate2.tsv",sep="\t")

##binary ploidy
binary.ploidy<-rep(0,dim(matched$dat)[1])
binary.ploidy[diploidsc]<-"0"
binary.ploidy[polyploids]<- "1"
binary.ploidy[diploidsi]<-"0"
binary.ploidy[sc]<-"?"
binary.ploidy[diploidsna]<-"0"
binary.ploidy[idcp]<-"?"
binary.ploidy<-as.data.frame(binary.ploidy)
row.names(binary.ploidy)<-matched$phy$tip.label

write.table(binary.ploidy,file="~/Dropbox/solploidy/basicdata/binaryploidy2.tsv",sep="\t")

##binary breeding system
binary.si<-rep(0,dim(matched$dat)[1])
binary.si[diploidsc]<-"0"
binary.si[diploidsi]<-"1"
binary.si[polyploids]<-"0"
binary.si[sc]<-"0"
binary.si[diploidsna]<-"?"
binary.si[idcp]<-"?"
binary.si<-as.data.frame(binary.si)
row.names(binary.si)<-matched$phy$tip.label

write.table(binary.si,file="~/Dropbox/solploidy/basicdata/binarysi2.tsv",sep="\t")

