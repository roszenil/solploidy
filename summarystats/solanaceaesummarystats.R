setwd('~/Dropbox/solploidy/summarystats')
library("ape")
library("phytools")
library("treeplyr")

## Solanaceae time-tree with corrected names Sarkinen et al 2013
solanaceae.tree<-read.tree("~/Dropbox/solploidy/basicdata/solcorrecttip.tre")

# Dataset created by RZF from CCDB and independent records by EEG
#idnumber, source,fullname, genus, species,cultivarhybrid,csomenumber,lifehistoryclean, ploidy selfincomp, ploidybyZF, namechangedfrom
solanaceae.allrecords<-read.csv("~/Dropbox/solploidy/basicdata/solanaceaerecords.csv",header=TRUE,sep=",", stringsAsFactors=FALSE, na.strings="")

cultivars<-which(solanaceae.allrecords$cultivarhybrid==1)
#How many species, how many genera
species<-unique(solanaceae.allrecords$fullname)
total.species<-length(species)# 2270 total species
genera<-unique(solanaceae.allrecords$genus)
total.genera<-length(genera)#99 genera


solanaceae.allrecords<-solanaceae.allrecords[-cultivars,]

#Calculating the mode(s)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux,incomparables=NA)))]
}

#Creating a very simple database with summary stats
mode.csome<-rep(0,total.species)
min.csome<-rep(0,total.species)
mode.ploidy<-rep(0,total.species)
min.ploidy<-rep(0,total.species)
self.incomp<-rep(0,total.species)
life.hist<-rep(0,total.species)

for (i in 1:total.species){
  index<-which(solanaceae.allrecords$fullname==species[i])
  aux1<- solanaceae.allrecords$csomenumber[index]
  if(all(is.na(aux1))){
  	 mode.csome[i]<-NA
  	 min.csome[i]<-NA  	
  }else{
  mode.csome[i]<-min(Mode(aux1))
  min.csome[i]<-min(aux1,na.rm=TRUE)
  }
  aux2<- solanaceae.allrecords$ploidy[index]
  if(all(is.na(aux2))){
  	 mode.ploidy[i]<-NA
  	 min.ploidy[i]<-NA  	
  }else{
  mode.ploidy[i]<-min(Mode(aux2))
  min.ploidy[i]<-min(aux2,na.rm=TRUE)
  }	
  self.incomp[i]<-solanaceae.allrecords$selfincomp[index[1]]
  life.hist[i]<-solanaceae.allrecords$lifehistoryclean[index[1]]
}

#Simplified contains summary stats for each spces
simplified.solanaceae<-data.frame(species, min.csome, mode.csome,min.ploidy,mode.ploidy,self.incomp, life.hist)


nor.ploidysi<-which(is.na(simplified.solanaceae[,5])==TRUE & is.na(simplified.solanaceae[,6])==TRUE)

norploidysi.solanaceae<-simplified.solanaceae[-c(nor.ploidysi),]

###########################################
### Matching tree with matching dataset if they have either ploidy or self-incomp information that is 


##Just checking here that all data has either ploidy or breeding system info
##nor.ploidysi2<-which(is.na(norploidysi.solanaceae[,5])==TRUE & is.na(norploidysi.solanaceae[,6])==TRUE)

matched<-make.treedata(solanaceae.tree,norploidysi.solanaceae) #635 species!!!!
write.nexus(matched$phy,file="~/Dropbox/solploidy/basicdata/fullmatchtree.nex")
#matched$dat=as.data.frame(matched$dat,row.names=matched$phy$tip.label)
#write.csv(matched$dat, file="~/Dropbox/solploidy/simplifiedtable.csv")

diploidsi<-which(matched$dat[,4]==2 & matched$dat[,5]=="SI")# 120
diploidsc<-which(matched$dat[,4]==2 & matched$dat[,5]=="SC")#159
diploidsna<-which(matched$dat[,4]==2 & is.na(matched$dat[,5])==TRUE) #192
polyploids<-which(matched$dat[,4]>2) #102
polyploidsc<-which(matched$dat[,4]>2 & matched$dat[,5]=="SC")#82
polyploidsna<-which(matched$dat[,4]>2 & is.na(matched$dat[,5])==TRUE)#19
si<-which(is.na(matched$dat[,4])==TRUE & matched$dat[,5]=="SI")#12
sc<-which(is.na(matched$dat[,4])==TRUE & matched$dat[,5]=="SC")#50
threestate<- rep(0,dim(matched$dat)[1])
threestate[diploidsc]="0"
threestate[polyploids]="1"
threestate[diploidsi]="2"
threestate[si]="2"
threestate[sc]="(0 1)"
threestate[diploidsna]<-"(0 2)"
threestate<-as.data.frame(threestate)
row.names(threestate)<-matched$phy$tip.label

write.table(threestate,file="~/Dropbox/solploidy/basicdata/threestate2.tsv",sep="\t")

##binary ploidy
binary.ploidy<-rep(0,dim(matched$dat)[1])
binary.ploidy[diploidsc]<-"0"
binary.ploidy[diploidsi]<-"0"
binary.ploidy[polyploids]<-"1"
binary.ploidy[si]<-"0"
binary.ploidy[sc]<-"?"
binary.ploidy[diploidsna]<-"0"
binary.ploidy<-as.data.frame(binary.ploidy)
row.names(binary.ploidy)<-matched$phy$tip.label

write.table(binary.ploidy,file="~/Dropbox/solploidy/basicdata/binaryploidy2.tsv",sep="\t")

##binary breeding system
binary.si<-rep(0,dim(matched$dat)[1])
binary.si[diploidsc]<-"0"
binary.si[diploidsi]<-"1"
binary.si[polyploids]<-"0"
binary.si[si]<-"1"
binary.si[sc]<-"0"
binary.si[diploidsna]<-"?"
binary.si<-as.data.frame(binary.si)
row.names(binary.si)<-matched$phy$tip.label

write.table(binary.si,file="~/Dropbox/solploidy/basicdata/binarysi2.tsv",sep="\t")


#########################
### Here are previous matches using full info but this make models non-comparable since the input is different sample sizes


no.ploidy=which(is.na(simplified.solanaceae$mode.ploidy)==TRUE) #1366
no.csome=which(is.na(simplified.solanaceae$mode.csome)==TRUE) #1248
no.si<-which(is.na(simplified.solanaceae$self.incomp)==TRUE)#1749

####Creating dataset with with ploidy and selfincompatibility info that matches the tree results in 220 species
noploidyandsi<-union(no.ploidy, no.si)#1960 without ploidy or selfincompatibility
#Dataset for all solanaceae that have ploidy and self incompatibility308 species
ploidysi.solanaceae<-simplified.solanaceae[-c(noploidyandsi),]

#Matching for everything
matched<-make.treedata(solanaceae.tree,ploidysi.solanaceae, name_column="species") #220 out of 308 that match the tree

ploidysi.dataset<-data.frame(matched$dat)
row.names(ploidysi.dataset)<-matched$phy$tip.label
ploidysi.tree<-matched$phy
plotTree(ploidysi.tree,type="fan",fsize=0.4,lwd=1) #220 tips
a<-hist(ploidysi.dataset$mode.ploidy)
barplot(ploidysi.dataset$self.incomp)

### Generating tree and values for ploidy and self-incomp analyses
write.nexus(ploidysi.tree,file="~/Dropbox/solploidy/basicdata/ploidysi.nex")
aux1<-which(ploidysi.dataset$mode.ploidy>2)# 33 polyploids
aux2<-which(ploidysi.dataset$mode.ploidy==2 & ploidysi.dataset$self.incomp=="SI") # 74 Self-Incompatible 
threestate.var<-rep(0,220) 
threestate.var[aux1]=1
threestate.var[aux2]=2
threestate.var<-data.frame(threestate.var)
threestate.table<-cbind(ploidysi.tree$tip.label, threestate.var)

write.csv(threestate.table,file="~/Dropbox/solploidy/basicdata/ploidysi.csv",row.names=FALSE)
#######################
####Creating dataset with with ploidy only that matches the tree results

ploidy.solanaceae<-simplified.solanaceae[-c(no.ploidy),]

#Matching for everything
matched<-make.treedata(solanaceae.tree,ploidy.solanaceae, name_column="species") #423

ploidy.dataset<-data.frame(matched$dat)
row.names(ploidy.dataset)<-matched$phy$tip.label
ploidy.tree<-matched$phy
plotTree(ploidy.tree,type="fan",fsize=0.2,lwd=1)
a<-hist(ploidy.dataset$mode.ploidy)

write.tree(ploidy.tree,file="~/Dropbox/solploidy/basicdata/ploidy.tre")
write.nexus(ploidy.tree, file="~/Dropbox/solploidy/basicdata/ploidy.nex")
tip.names<-ploidy.tree$tip.label
write.table(ploidy.dataset$mode.ploidy,file="~/Dropbox/solploidy/basicdata/ploidy.txt",sep=',',row.names=FALSE, col.names=TRUE)
aux1<-which(ploidy.dataset$mode.ploidy==2)#360 diploids
aux<-which(ploidy.dataset$mode.ploidy!=2)#63 polyploids 14.8% of the sample
binary.var<-rep(1,423) 
binary.var[aux1]=0
binary.var<-data.frame(binary.var)
binary.table<-cbind(tip.names, binary.var)
write.table(binary.table,file="~/Dropbox/solploidy/basicdata/binaryploidy.tsv",sep="\t",row.names=FALSE, col.names=FALSE)

########### Data set with only self incompatibility

si.solanaceae<-simplified.solanaceae[-c(no.si),]
matched<-make.treedata(solanaceae.tree,si.solanaceae, name_column="species") #347
si.dataset<-data.frame(matched$dat)
row.names(si.dataset)<-matched$phy$tip.label
si.tree<-matched$phy
tip.names<-si.tree$tip.label
plotTree(si.tree,type="fan",fsize=0.2,lwd=1)
write.tree(si.tree,file="~/Dropbox/solploidy/basicdata/si.tre")
write.nexus(si.tree, file="~/Dropbox/solploidy/basicdata/si.nex")
aux1<-which(si.dataset$self.incomp=="SC")#231 self compatible taxa
aux<-which(si.dataset$self.incomp=="SI")#116 self incompatible txa
binary.var<-rep(1,347) 
binary.var[aux1]=0
binary.var<-data.frame(binary.var)
binary.table<-cbind(tip.names, binary.var)
write.csv(binary.table,file="~/Dropbox/solploidy/basicdata/binarysi.csv",row.names=FALSE)


#Building a dataset with species and all of their csome numbers not only a single summary
datamat<-matrix(rep(NA,32*total.species),ncol=32)
for (i in 1:total.species){
  index<-which(solanaceae.allrecords$fullname==species[i])
  aux1<- table(solanaceae.allrecords$csomenumber[index])
  csomes<-as.numeric(names(aux1))
  max.chrom<-length(csomes)
  datamat[i, 31]<-max.chrom
  if(max.chrom>1){
  	datamat[i,32]<-1
  }
  if(max.chrom>0){
  datamat[i,1:length(csomes)]<-csomes
  freq<-as.numeric(aux1)
  datamat[i,16:(15+length(freq))]<-freq
  }
}
# I found one species Solanum nigrum with 15 different chromosome numbers (wow!)
csomefreq<-data.frame(species,datamat,stringsAsFactors=FALSE)
names(csomefreq)<-c("taxa","csome1","csome2","csome3","csome4","csome5","csome6","csome7","csome8","csome9","csome10","csome11","csome12","csome13","csome14","csome15","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freq11","freq12","freq13","freq14","freq15","levels","chromseries")
#179 chromosome series taxa


no.csome=which(is.na(csomefreq$csome1)==TRUE)
csomefreq2<-csomefreq[-no.csome,]
#matching every species on the database with the tree it's 846 taxa
matched1<-make.treedata(solanaceae.tree,csomefreq,name_column="taxa")

csomeseriesdat1<-data.frame(matched1$dat)
row.names(csomeseriesdat1)<-matched1$phy$tip.label

write.tree(matched1$phy,file="~/Dropbox/solploidy/csomeseries/soltreealltaxa.tre")
write.csv(csomeseriesdat1,file="~/Dropbox/solploidy/csomeseries/soldataalltaxa.csv")

#matching every species that actually has chromsome numbers with the tree results in 593. That is  a 253 species difference
matched2<-make.treedata(solanaceae.tree, csomefreq2, name_column="taxa")

csomeseriesdat2<-data.frame(matched2$dat)
row.names(csomeseriesdat2)<-matched2$phy$tip.label

write.tree(matched2$phy,file="~/Dropbox/solploidy/csomeseries/soltreecsometaxa.tre")
write.csv(csomeseriesdat2,file="~/Dropbox/solploidy/csomeseries/soldatacsometaxa.csv")


##TRYING TO PLOT THE DATA BUT NOT READY YET
library("wesanderson")
library("RColorBrewer")
tree<-reorder(solanum.tree[[1]],"cladewise")
X<-solanum.dataset[tree$tip.label,]
plotTree(solanum.tree[[1]],plot=FALSE)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
plotTree(solanum.tree[[1]],lwd=1,ylim=c(0,obj$y.lim[2]*1.05),xlim=c(0,obj$x.lim[2]*1.3),
    ftype="off")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
h<-max(obj$xx)
fsize<-0.3
for(i in 1:Ntip(tree)){
    lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted")
    text(h,obj$yy[i],tree$tip.label[i],cex=fsize,pos=4,font=3,offset=0.1)
}
s<-max(fsize*strwidth(tree$tip.label))
start.x<-1.05*h+s

cols<-list()
cols[[1]]<-c(wes_palette("Darjeeling1",2))
cols[[2]]<-c(wes_palette("Darjeeling1",5))

for(i in 1:ncol(X)){
    text(start.x,max(obj$yy)+1,paste("trait",colnames(X)[i]),pos=4,srt=60,
        cex=0.8,offset=0)
        
       cols[[i]]<-setNames(length(levels(X[[i]]))),levels(X[[i]]))
    for(j in 1:nrow(X)){
        xy<-c(start.x,obj$yy[j])
        y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
        asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
            par()$pin[2]/par()$pin[1]
        x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
        polygon(x,y,col=cols[[i]][as.character(X[[i]][j])])
    }
    start.x<-start.x+2*asp
}
start.y<-max(obj$yy)
for(i in 1:ncol(X)){
    text(start.x,start.y,paste("trait",colnames(X)[i]),pos=4,cex=0.9,
        offset=0)
    add.simmap.legend(colors=cols[[i]],shape="square",prompt=FALSE,
        x=start.x,y=start.y-2*strheight("W")*0.9,fsize=0.9)
    start.y<-start.y-1.5*0.9*strheight("W")*(length(cols[[i]])-1)-6
}
