setwd('~/Dropbox/solploidy/basicdata')
library("ape")
library("phytools")
library("treeplyr")
solanaceae.tree<-read.tree("solanaceae2013.tre")


#idnumber, source,fullname, genus, species,cultivarhybrid,csomenumber,lifehistoryclean, ploidy selfincomp, ploidybyZF
solanaceae.allrecords<-read.csv("solanaceaerecords.csv",header=TRUE,sep=",", stringsAsFactors=FALSE, na.strings="")

species<-unique(solanaceae.allrecords$fullname)
total.species<-length(species)# 2278 total species
genera<-unique(solanaceae.allrecords$genus)
total.genera<-length(genera)#97 genera

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

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


simplified.solanaceae<-data.frame(species, min.csome, mode.csome,min.ploidy,mode.ploidy,self.incomp, life.hist)

plotTree(solanaceae.tree,plot=FALSE,type="fan",fsize=0.05,lwd=1)
matched<-make.treedata(solanaceae.tree,simplified.solanaceae, name_column="species"). #864

no.ploidy=which(is.na(matched$dat[,3])==TRUE) #283
no.csome=which(is.na(matched$dat[,1])==TRUE) #257
no.si<-which(is.na(matched$dat[,6])==TRUE)#351
plotTree(matched$phy,type="fan",fsize=0.1,lwd=1)
plotTree(matched$phy,type="fan",fsize=0.1,lwd=1)
###########

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
