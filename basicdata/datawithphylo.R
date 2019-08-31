library("ggplot2")
library("ggtree")
library("treeio")
library("ggnewscale")
#library("RevGadgets")
library("wesanderson")
library("dyplr")

cols<-c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2"))
cols2<- c("#7b3294","#c2a5cf","#a6dba0","#008837","#ffffbf")
cols3<-c("#ece2f0", "#67a9cf","#e31a1c","#fd8d3c","#02818a","#014636")


soltree<-read.nexus("~/Dropbox/solploidy/basicdata/fullmatchtree.nex")
soldata<- read.csv("~/Dropbox/solploidy/basicdata/discretematrix.csv")
dp<-data.frame(soldata[,2:4])
rownames(dp)<-soltree$tip.label

circ <- ggtree(soltree, layout = "circular")

p1 <- gheatmap(circ, dp[, "Ploidy", drop=F], offset=.8, width=.1,
               colnames_angle=90, colnames_offset_y = .25)
p2 <- gheatmap(p1, dp[, "BS", drop=F], offset=5, width=.1,
               colnames_angle=90, colnames_offset_y = .25)
p3 <- gheatmap(p2, dp[, "Ploidybs", drop=F], offset=10, width=.1,
               colnames_angle=90, colnames_offset_y = .25)

col <- c(cols2[2],cols[7],cols2[1],cols[1],cols[9],cols[9],cols2[3],cols2[4],"yellow",cols2[3],cols[6],cols[6],"grey")
names(col) = c("C","CD","Compatible","CP","D","Diploid","I","ID","IDCP","Incompatible","P","Polyploid","Unknown")
#D = cols[9]
#%P=cols[6]

#CD cols[7]
#CP cols[1]

#Compatible=cols2[1],
#C= cols2[2]
#Incompatible= cols2[3]
#I cols2[3])
#ID= cols2[4]
#IDCP= "yellow"
#Diploid= cols[9]
#Polyploid= cols[6]

pp <- p3 + scale_fill_manual(values=col)
require(phytools)
plotTree(soltree,type="fan", cex=0.18,fsize=0.19,label.offset= 2)
