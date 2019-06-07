library("ape")
library("phytools")
library("treeplyr")

## Solanaceae time-tree with corrected names Sarkinen et al 2013
solanaceae.tree<-read.tree("~/Dropbox/solploidy/basicdata/solcorrecttip.tre")

cultivars<-which(solanaceae.allrecords$cultivarhybrid==1)
#How many species, how many genera
species<-unique(solanaceae.allrecords$fullname)
total.species<-length(species)# 2270 total species
genera<-unique(solanaceae.allrecords$genus)
total.genera<-length(genera)#99 genera

# Dataset created by RZF from CCDB and independent records by EEG
#idnumber, source,fullname, genus, species,cultivarhybrid,csomenumber,lifehistoryclean, ploidy selfincomp, ploidybyZF, namechangedfrom
solanaceae.allrecords<-read.csv("~/Dropbox/solploidy/basicdata/solanaceaerecords.csv",header=TRUE,sep=",", stringsAsFactors=FALSE, na.strings="")

cultivars<-which(solanaceae.allrecords$cultivarhybrid==1)
solanaceae.allrecords<-solanaceae.allrecords[-cultivars,]

uncertain.species<-c("Nierembergia_rigida","Anisodus_luridus","Lycium_europaeum","Lycium_chilense","Lycium_ciliatum","Chamaesaracha_coronopus","Chamaesaracha_sordida","Physalis_hederifolia","Solanum_immite","Solanum_acroscopicum","Solanum_multiinterruptum","Solanum_brevicaule","Solanum_bahamense","Solanum_juvenale","Solanum_hieronymi","Solanum_campylacanthum","Solanum_discolor","Mandragora_officinarum")

long1<-length(uncertain.species)
for (i in 1:long1){
	aux1<-which(solanaceae.allrecords$fullname==uncertain.species[i])
	solanaceae.allrecords<-solanaceae.allrecords[-aux1,]
}

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
simplified.solanaceae<-data.frame(species, min.csome, mode.csome,min.ploidy,mode.ploidy,self.incomp, life.hist)


nor.ploidysi<-which(is.na(simplified.solanaceae[,5])==TRUE & is.na(simplified.solanaceae[,6])==TRUE)

norploidysi.solanaceae<-simplified.solanaceae[-c(nor.ploidysi),]
matched<-make.treedata(solanaceae.tree,norploidysi.solanaceae) #636 species!!!!
diploidsi<-which(matched$dat[,4]==2 & matched$dat[,5]=="SI")# 120
diploidsc<-which(matched$dat[,4]==2 & matched$dat[,5]=="SC")#159
diploidsna<-which(matched$dat[,4]==2 & is.na(matched$dat[,5])==TRUE) #192
polyploids<-which(matched$dat[,4]>2) #102
polyploidsc<-which(matched$dat[,4]>2 & matched$dat[,5]=="SC")#82
polyploidsna<-which(matched$dat[,4]>2 & is.na(matched$dat[,5])==TRUE)#19
si<-which(is.na(matched$dat[,4])==TRUE & matched$dat[,5]=="SI")#12
sc<-which(is.na(matched$dat[,4])==TRUE & matched$dat[,5]=="SC")#50

threestate.matrix<-matrix(rep(0,3*(dim(matched$dat)[1])),ncol=3)
threestate.matrix[diploidsc,1]<- 0
threestate.matrix[diploidsc,2]<- 1
threestate.matrix[diploidsc,3]<- 0

threestate.matrix[polyploids,1]<-0
threestate.matrix[polyploids,2]<- 0
threestate.matrix[polyploids,3]<- 1

threestate.matrix[diploidsi,1]<- 1
threestate.matrix[diploidsi,2]<- 0
threestate.matrix[diploidsi,3]<- 0

threestate.matrix[si,1]<- 1
threestate.matrix[si,2]<- 0
threestate.matrix[si,3]<- 0

threestate.matrix[sc,1]<- 0
threestate.matrix[sc,2]<- 1
threestate.matrix[sc,3]<- 1

threestate.matrix[diploidsna,1]<- 1
threestate.matrix[diploidsna,2]<- 1
threestate.matrix[diploidsna,3]<- 0


threestate.matrix<-as.data.frame(threestate.matrix)
threestate.matrix<-cbind(matched$phy$tip.label,threestate.matrix)
### Corrections start here Columns are ID CD CP
# 1. Some challenging species
# Solanum_chacoense # ID
# Solanum_candolleanum # ID
# Solanum_elaeagnifolium # CD
# Solanum_chenopodinum # CD
# Lycium_gariepense  Solution CD and CP
aux2<-which(threestate.matrix[,1]=="Solanum_chacoense")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Solanum_candolleanum")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Solanum_elaeagnifolium")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Solanum_chenopodinum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 0

#4. Species that are informative for the three state but are not for neither of the two-state models (ID and CP simultaneously)= Solution leave them as is
#Lycium_californicum 
#Solanum_bulbocastanum
aux2<-which(threestate.matrix[,1]=="Lycium_californicum")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_bulbocastanum")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 1
aux2<-which(threestate.matrix[,1]=="Lycium_gariepense")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1


#5. Species that Rosana needs to correct by hand (automatic processing fails)
#Lycium_exsertum # This I need to correct by hand it is D/P but always SC
#Datura_metel    # This I need to correct by hand it is D/P but always SC
#Withania_somnifera  # This I need to correct by hand it is D/P but always SC
#Physalis_angulata   # This I need to correct by hand it is D/P but always SC
#Solanum_nigrum       # This I need to correct by hand it is D/P but always SC
#Solanum_etuberosum   # This I need to correct by hand it is D/P but always SC
#Solanum_andreanum    # This I need to correct by hand it is D/P but always SC
#Solanum_colombianum  # This I need to correct by hand it is D/P but always SC
#Solanum_torvum       # This I need to correct by hand it is D/P but always SC
#Solanum_capsicoides  # This I need to correct by hand it is D/P but always SC
#Solanum_pectinatum   # This I need to correct by hand it is D/P but always SC
#Solanum_dioicum      # This I need to correct by hand it is D/P but always SC
#Solanum_cinereum       # This I need to correct by hand it is D/P but always SC

aux2<-which(threestate.matrix[,1]=="Lycium_exsertum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Datura_metel")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Withania_somnifera")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Physalis_angulata")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_nigrum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1


aux2<-which(threestate.matrix[,1]=="Solanum_etuberosum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1


aux2<-which(threestate.matrix[,1]=="Solanum_andreanum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_colombianum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1


aux2<-which(threestate.matrix[,1]=="Solanum_torvum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_colombianum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_capsicoides")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_pectinatum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_dioicum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_cinereum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 1

## 6. Weird processing error
#Lycianthes_rantonnetii  # Weird processing error, need to correct to ID
#Capsicum_cardenasii     # Weird processing error, need to correct to ID
#Witheringia_meiantha    # Weird processing error, need to correct to ID
#Solanum_stipuloideum    # Weird processing error need to correct to ID
#Solanum_cyaneopurpureum # Weird processing error should be CD
#Solanum_leopoldensis    # Weird processing error should be CD

aux2<-which(threestate.matrix[,1]=="Lycianthes_rantonnetii")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Capsicum_cardenasii")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Witheringia_meiantha")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Solanum_stipuloideum")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Solanum_cyaneopurpureum")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Solanum_leopoldensis")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 0

aux2<-which(threestate.matrix[,1]=="Hyoscyamus_turcomanicus")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 0
threestate.matrix[aux2,4]= 1

aux2<-which(threestate.matrix[,1]=="Solanum_fructu-tecto")
threestate.matrix[aux2,2]= 0
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 0


#Solanum_proteanthum - 24csomes but ID CD
aux2<-which(threestate.matrix[,1]=="Solanum_proteanthum")
threestate.matrix[aux2,2]= 1
threestate.matrix[aux2,3]= 1
threestate.matrix[aux2,4]= 0


#### This is the checking statge with statesICDP
dat1 <- read.csv("~/Dropbox/solploidy/basicdata/statesICDP.csv", as.is=T)

rzfaddition<-0
notmatching<-0
for(i in 1:length(matched$phy$tip.label)){
	speciestomatch<-which(dat1$Species==matched$phy$tip.label[i])
	if(length(speciestomatch)==0){
		rzfaddition<-c(rzfaddition,i)
	}else{
		aux1<-as.numeric(dat1[speciestomatch,2:4])
		aux2<-as.numeric(threestate.matrix[i,2:4])
		matches<-rep(0,3)
		matches[1]<-((aux1[1]>0 & aux2[1]>0)|(aux1[1]==0 & aux2[1]==0))
		matches[2]<-((aux1[2]>0 & aux2[2]>0)|(aux1[2]==0 & aux2[2]==0))
		matches[3]<-((aux1[3]>0 & aux2[3]>0)|(aux1[3]==0 & aux2[3]==0))
		differences<-sum(matches)
		if(differences!=3){
			notmatching<-c(notmatching,i)
		}
	}
}

threestate.matrix[rzfaddition[-1],]
threestate.matrix[notmatching[-1],]

write.csv(threestate.matrix,file="~/Dropbox/solploidy/basicdata/newstatesICDP.csv")


