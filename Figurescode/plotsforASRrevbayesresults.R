#biocLite("ggtree")
#biocLite("treeio")
library("ggplot2")
library("RevGadgets")
library("wesanderson")

### Colors
cols<-c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2"))
cols2<- c("#7b3294","#c2a5cf","#a6dba0","#008837","#ffffbf")
cols3<-c("#ece2f0", "#67a9cf","#e31a1c","#fd8d3c","#02818a","#014636")
#sampling<-seq(1,250000,100)
### Plots for diversification rates Bisse DP with diploidization
#setwd("~/Dropbox/solploidypersonal/Figures")
source("multiplot.R")
#source("plot_ancestral_states_2.R")
output.sse<-read.table("~/Dropbox/solploidypersonal/bisse250K/output/BiSSE_polydip250K.log", header=TRUE)
sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$rate_12,output.sse$rate_21),Type=rep(c("Polyploidization","Diploidization"),each=length(output.sse$rate_1)))


p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[1]))

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols[7],cols[1]))

p3.1<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[1]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[1]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[2],cols[4]))
multiplot(p1,p3.1,p5,p2,p4, cols=2)

############################################################################################
### Plots for diversification rates Bisse DP without diploidization~
output.sse<-read.table("~/Dropbox/solploidypersonal/bissenodip250K/output/BiSSE_polynodip250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2),Type=rep(c("Diploid","Polyploid"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=output.sse$rate_12,Type=rep("Polyploidization",length(output.sse$rate_1)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[1]))


p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols[7],cols[1]))

p3.1nodip<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[1]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[1]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols[4]))
multiplot(p1,p3.1nodip,p5,p2,p4, cols=2)


##########


############################################################################################
### Plots for diversification rates Hisse DP with diploidization~/Dropbox/solploidypersonal/hisse250K/output/HiSSE_polydip250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/hisse250K/output/HiSSE_polydip250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3,output.sse$extinction.4),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3,output.sse$speciation.4),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3,output.sse$speciation.4-output.sse$extinction.4),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3,output.sse$extinction.3/output.sse$speciation.3),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$Q.1,output.sse$Q.2) ,Type=rep(c("Diploidization", "Polyploidization"),each=length(output.sse$Q.1)))
hidden.rate<-data.frame(dens=output.sse$R.1 ,Type=rep("AB=BA",length(output.sse$rate_1)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[9],cols[1],cols[6]))

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols[7],cols[9],cols[1],cols[6]))

p3.2<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[9],cols[1],cols[6]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[9],cols[1],cols[6]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[4],cols[2]))

p6<-ggplot(hidden.rate, aes(x=dens, fill=Type))+labs(title="Hidden trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols[10]))
multiplot(p1,p3.2,p5,p2,p4,p6, cols=2)
##########
############################################################################################
### Plots for diversification rates Hisse DP without diploidization~/Dropbox/solploidypersonal/hissenodip250K/output/HiSSE_polynodip250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/hissenodip250K/output/HiSSE_polynodip250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3,output.sse$extinction.4),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3,output.sse$speciation.4),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3,output.sse$speciation.4-output.sse$extinction.4),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3,output.sse$extinction.3/output.sse$speciation.3),Type=rep(c("Diploid A","Polyploid_A", "Diploid_B", "Polyploid_B"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=(output.sse$Q.1) ,Type=rep("Polyploidization",length(output.sse$Q.1)))
hidden.rate<-data.frame(dens=output.sse$R.1 ,Type=rep("AB=BA",length(output.sse$rate_1)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[9],cols[1],cols[6]))

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols[7],cols[9],cols[1],cols[6]))+xlim(0,2)

p3.2nodip<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[9],cols[1],cols[6]))+xlim(-0.5,1.5)

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[7],cols[9],cols[1],cols[6]))+xlim(0,2)


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[4]))+xlim(0,0.15)


p6<-ggplot(hidden.rate, aes(x=dens, fill=Type))+labs(title="Hidden trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols[10]))+xlim(0,0.15)
multiplot(p1,p3.2nodip,p5,p2,p4,p6, cols=2)
##########


############################################################################################
### Plots for diversification rates Bisse Breeding systems ~/Dropbox/solploidypersonal/bisseSInoreturn250K/output/BiSSE_sinoreturn250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/bisseSInoreturn250K/output/BiSSE_sinoreturn250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2),Type=rep(c("SC","SI"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2),Type=rep(c("SC","SI"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2),Type=rep(c("C","I"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2),Type=rep(c("SC","SI"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=output.sse$rate_21,Type=rep("SI to SC",length(output.sse$rate_21)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[3]))


p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols2[1],cols2[3]))

p3.3<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[3]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[3]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols2[5]))
multiplot(p1,p3.3,p5,p2,p4, cols=2)
##########

############################################################################################
### Plots for diversification rates Bisse Breeding systems with return ~/Dropbox/solploidypersonal/bisseSI250K/output/BiSSE_selfincomp250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/bisseSI250K/output/BiSSE_selfincomp250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2),Type=rep(c("SC","SI"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2),Type=rep(c("SC","SI"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2),Type=rep(c("C","I"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2),Type=rep(c("SC","SI"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$rate_12,output.sse$rate_21),Type=rep(c("SC to SI","SI to SC"),each=length(output.sse$rate_21)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[3]))


p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols2[1],cols2[3]))

p3.31<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[3]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[3]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[4],cols2[5]))
multiplot(p1,p3.31,p5,p2,p4, cols=2)
##########


############################################################################################
### Plots for diversification rates Hisse for breeding systems ~/Dropbox/solploidypersonal/hisseSInoreturn250K/output/HiSSE_sinoret250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/hisseSInoreturn250K/output/HiSSE_sinoret250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3,output.sse$extinction.4),Type=rep(c("SC A","SI A", "SC B", "SI B"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3,output.sse$speciation.4),Type=rep(c("SC A","SI A", "SC B", "SI B"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3,output.sse$speciation.4-output.sse$extinction.4),Type=rep(c("SC A","SI A", "SC B", "SI B"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3,output.sse$extinction.3/output.sse$speciation.3),Type=rep(c("SC A","SI A", "SC B", "SI B"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=(output.sse$Q.2) ,Type=rep("SI to SC",length(output.sse$Q.2)))
hidden.rate<-data.frame(dens=output.sse$R.1 ,Type=rep("AB=BA",length(output.sse$R.1)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[2],cols2[4],cols2[3]))

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[2],cols2[4],cols2[3]))

p3.4<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =  c(cols2[1],cols2[2],cols2[4],cols2[3]))



p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[1],cols2[2],cols2[4],cols2[3]))

p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols2[5]))


p6<-ggplot(hidden.rate, aes(x=dens, fill=Type))+labs(title="Hidden trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols[10]))
multiplot(p1,p3.4,p5,p2,p4,p6, cols=2)
##########

############################################################################################~
### Plots for diversification rates Musse DP and breeding systems~/Dropbox/solploidypersonal/mussefull250k/output/MuSSE_ploidysi250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/mussefull250k/output/MuSSE_ploidysi250K.log", header=TRUE)


sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$rate_12,output.sse$rate_21,output.sse$rate_31,output.sse$rate_32) ,Type=rep(c("Polyploidization from SC","Diploidization","SI to SC","Polyploidization from SI"),each=length(output.sse$rate_12)))


p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))

p3.5<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =c(cols[2],cols[4],"pink",cols2[5]))

multiplot(p1,p3.5,p5,p2,p4, cols=2) 
#################################################################################################
### Plots for diversification rates Musse DP and breeding systems~/Dropbox/solploidypersonal/mussenodip250k/output/MuSSE_ploidysi250K.log

output.sse<-read.table("~/Dropbox/solploidypersonal/mussenodip250k/output/MuSSE_ploidysinodip250K.log", header=TRUE)

sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3),Type=rep(c("SC-Diploid","Polyploid", "SI-Diploid"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$rate_12,output.sse$rate_31,output.sse$rate_32) ,Type=rep(c("Polyploidization from SC","SI to SC","Polyploidization from SI"),each=length(output.sse$rate_12)))


p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))

p3.5nodip<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols3[3],cols3[2],cols2[4]))


p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[4],"pink",cols2[5]))

multiplot(p1,p3.5nodip,p5,p2,p4, cols=2)
#########


############################################################################################~
### Plots for diversification rates MuHiSSE DP and breeding systems~/Dropbox/solploidypersonal/muhisse250k/output/MuHiSSE_ploidysi250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/muhisse250k/output/MuHiSSE_ploidysi250K.log", header=TRUE)


sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3, output.sse$extinction.4,output.sse$extinction.5,output.sse$extinction.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3,output.sse$speciation.4,output.sse$speciation.5,output.sse$speciation.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$speciation.1)))

#sse.speciation1<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2),Type=rep(c("SC-D A","Polyploid A"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3,output.sse$speciation.4-output.sse$extinction.4,output.sse$speciation.5-output.sse$extinction.5, output.sse$speciation.6-output.sse$extinction.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3,output.sse$extinction.4/output.sse$speciation.4,output.sse$extinction.5/output.sse$speciation.5,output.sse$extinction.6/output.sse$speciation.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$rate_12,output.sse$rate_21,output.sse$rate_31,output.sse$rate_32) ,Type=rep(c("Polyploidization from SC","Diploidization","SI to SC","Polyploidization from SI"),each=length(output.sse$rate_12)))

hidden.rate<-data.frame(dens=output.sse$R.1 ,Type=rep("AB=BA",length(output.sse$R.1)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(0,2)

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(0,2)

p3.6<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(-2,2)

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(0,3)

p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[2],cols[4],"pink",cols2[5]))

p6<-ggplot(hidden.rate, aes(x=dens, fill=Type))+labs(title="Hidden trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols[10]))

multiplot(p1,p3.6,p5,p2,p4,p6, cols=2)

########

############################################################################################~
### Plots for diversification rates MuHiSSE DP and breeding systems no diploidization ~/Dropbox/solploidypersonal/muhisse250k/output/MuHiSSE_ploidysi250K.log
output.sse<-read.table("~/Dropbox/solploidypersonal/muhissenodip250k/output/MuHiSSEnodip_ploidysi250K.log", header=TRUE)


sse.extinction<-data.frame(dens=c(output.sse$extinction.1,output.sse$extinction.2,output.sse$extinction.3, output.sse$extinction.4,output.sse$extinction.5,output.sse$extinction.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$extinction.1)))

sse.speciation<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2,output.sse$speciation.3,output.sse$speciation.4,output.sse$speciation.5,output.sse$speciation.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$speciation.1)))

#sse.speciation1<-data.frame(dens=c(output.sse$speciation.1,output.sse$speciation.2),Type=rep(c("SC-D A","Polyploid A"),each=length(output.sse$speciation.1)))

sse.netdiv<-data.frame(dens=c(output.sse$speciation.1-output.sse$extinction.1,output.sse$speciation.2-output.sse$extinction.2, output.sse$speciation.3-output.sse$extinction.3,output.sse$speciation.4-output.sse$extinction.4,output.sse$speciation.5-output.sse$extinction.5, output.sse$speciation.6-output.sse$extinction.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$speciation.1)))

sse.reldiv<-data.frame(dens=c(output.sse$extinction.1/output.sse$speciation.1,output.sse$extinction.2/output.sse$speciation.2,output.sse$extinction.3/output.sse$speciation.3,output.sse$extinction.4/output.sse$speciation.4,output.sse$extinction.5/output.sse$speciation.5,output.sse$extinction.6/output.sse$speciation.6),Type=rep(c("SC-D A","Polyploid A", "SI-Diploid A","SC-D B","Polyploid B", "SI-Diploid B"),each=length(output.sse$speciation.1)))

trait.rates<-data.frame(dens=c(output.sse$rate_12,output.sse$rate_31,output.sse$rate_32) ,Type=rep(c("Polyploidization from SC","SI to SC","Polyploidization from SI"),each=length(output.sse$rate_12)))

hidden.rate<-data.frame(dens=output.sse$R.1 ,Type=rep("AB=BA",length(output.sse$R.1)))

p1<-ggplot(sse.speciation, aes(x=dens, fill=Type))+labs(title="Speciation",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(0,2)

p2<-ggplot(sse.extinction, aes(x=dens, fill=Type))+labs(title="Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(0,2)

p3.6nodip<-ggplot(sse.netdiv, aes(x=dens, fill=Type))+labs(title="Net Diversification",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))

p4<-ggplot(sse.reldiv, aes(x=dens, fill=Type))+labs(title="Relative Extinction",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[1],cols[6],cols[7],cols[9],cols2[4],cols2[3]))+xlim(0,2)

p5<-ggplot(trait.rates, aes(x=dens, fill=Type))+labs(title="Trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = c(cols[4],"pink",cols2[5]))

p6<-ggplot(hidden.rate, aes(x=dens, fill=Type))+labs(title="Hidden trait change",x="Rate", y="Posterior Density")+geom_density(alpha=0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = (cols[10]))

multiplot(p1,p3.6nodip,p5,p2,p4,p6, cols=2)

########
multiplot(p3.1nodip, p3.5nodip, p3.2nodip, p3.6nodip, cols=2)






## Net diversification plots only
multiplot(p3.1,p3.3,p3.5,p3.2,p3.4,p3.6, cols=2)







# Plots of ancestral state reconstruction Solanaceae

#1. BiSSE no diploidization 250K runs of MCMC
p = plot_ancestral_states("~/Dropbox/solploidypersonal/bissenodip250K/output/anc_states_summaryBiSSEtreepolynodip130K.tree", 
                          size=0.4, include_start_states=FALSE, 
                          summary_statistic="MAP",tip_label_size=0.45, 
                          tip_label_offset=0,
                          #tip_label_offset=0.4,
                          tip_label_italics=TRUE,
                          node_label_size=FALSE,node_label_nudge_x=0.1,
                          shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6),
                          tree_layout="circular")
print (p)


#2. Full BiSSE with 250K runs of MCMC

p = plot_ancestral_states("~/Dropbox/solploidypersonal/bisse250K/output/anc_states_summaryBiSSEtree250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6),state_colors=c(cols[7],cols[1]))
print(p)

#3. Full HiSSE no diploidization with 250K runs of MCMC

p = plot_ancestral_states("~/Dropbox/solploidypersonal/hissenodip250K/output/anc_states_summaryHiSSEnodiptree250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6),state_colors=c(cols[7],cols[1], cols[9],cols[6]))
print(p)

#4. Full Hisse wwith 250K runs of MCMC
p = plot_ancestral_states("~/Dropbox/solploidypersonal/hisse250K/output/anc_states_summaryHiSSEtree250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6),state_colors=c(cols[7],cols[1], cols[9],cols[6]))
print(p)

#5. Full BiSSE for self incompatibility 0=SC 1=SI with 250K runs of MCMC

p = plot_ancestral_states("~/Dropbox/solploidypersonal/bisseSInoreturn250K/output/anc_states_summaryBiSSEtreesinoret250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6),state_colors=c(cols2[2],cols2[3]))
print(p)

cols2[1],cols2[4],cols2[2],cols2[3]


#6. HiSSE SI
p = plot_ancestral_states("~/Dropbox/solploidypersonal/hisseSInoreturn250K/output/anc_states_summaryHiSSEsinorettree250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6),state_colors=c(cols2[2],cols2[3],cols2[1],cols2[4]))
print(p)


#7.  MuSSE for self incompatibility without SI evolving into SC 0=SC 1=SI and with diploidization.  250K runs of MCMC

p = plot_ancestral_states("~/Dropbox/solploidypersonal/mussefull250k/output/anc_states_summaryMuSSEploidysi250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6), state_colors=c(cols3[2],cols3[3],cols2[4]))

#8.  MuSSE for self incompatibility without SI evolving into SC 0=SC 1=SI and no diploidization.  250K runs of MCMC
p = plot_ancestral_states("~/Dropbox/solploidypersonal/mussenodip250K/output/anc_states_summaryMuSSEploidysinodip250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.3,node_size_range=c(2, 6), state_colors=c(cols3[2],cols3[3],cols2[4]))
print(p)

#9. MuHiSSE for self incompatibility without SI evolving into SC 0=SC 1=SI and with diploidization.  250K runs of MCMC

p = plot_ancestral_states("~/Dropbox/solploidypersonal/muhisse250K/output/anc_states_summaryMuHiSSEploidysi250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.5,node_size_range=c(2, 6),state_colors=c(cols[7],cols[1], cols2[4],cols[9],cols[6],cols2[3]))

print(p)

#9. MuHiSSE for self incompatibility without SI evolving into SC 0=SC 1=SI and no diploidization.  250K runs of MCMC
p = plot_ancestral_states("~/Dropbox/solploidypersonal/muhissenodip250K/output/anc_states_summaryMuHiSSEploidysinodip250K.tree",summary_statistic="MAP", include_start_states=FALSE, tip_label_size=0.3, tip_label_offset=0.5,node_label_size=0.1,node_label_nudge_x=0.1,shoulder_label_size=2,alpha=.5,node_size_range=c(2, 6),state_colors=c(cols[7],cols[1], cols2[4],cols[9],cols[6],cols2[3]))

print(p)
