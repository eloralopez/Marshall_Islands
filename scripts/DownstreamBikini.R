setwd("~/Documents/GitHub/Bikini/mappedtoahya")
library("reshape2")
library(ggplot2)
library(stringr)
library(sciplot)
library(sinaplot)
library(ggforce)
library(gridExtra)
library(dplyr)
library(patchwork)
library(ggpubr)
meltfunc<-function(file) {
  AHB125<-read.delim(file, check.names = FALSE)
  header<-colnames(AHB125)
  samplenames<-header[4:length(header)]
  AHB125mdf<-melt(AHB125, id.vars="chrom.pos", measure.vars=c(samplenames), value.name="genotype",variable.name="sample")
  base<-basename(file)
  colony<-strsplit(base, "\\_")[[1]][1]
  mf<-strsplit(base, "\\_")[[1]][2]
  type<-strsplit(base, "\\_")[[1]][3]
  #print(type)
  pasted<-paste(colony,"_",mf,"_",type,"_melted.txt", sep="")
  print(pasted)
  write.table(AHB125mdf, file=pasted,sep="\t",quote=FALSE, row.name=FALSE)
}

#files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini", pattern="*26.txt", full.names=T, recursive=FALSE)
#files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini", pattern="*POP.txt", full.names=T, recursive=FALSE)
#files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini", pattern="*_svPOP.txt", full.names=T, recursive=FALSE)
files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/redo", pattern="*.txt", full.names=T, recursive=FALSE)

files
for (i in 1:length(files)) {
  #meltfunc(i)
  file =files[i]
  #print(file)
  meltfunc(file)
}
#meltfunc("AHB70-90_mf_svPOP.txt")
meltfunc()
#denominators:
AHB125denom<-123410268-(895*14)
AHB145denom<-108352551-(895*14)
AHB151denom<-112353652-(895*14)
AHB176denom<-115474965-(895*14)
AHB195denom<-122257120-(895*14)
AHB70denom<-119454920-(895*14)
AHB90denom<-117673579-(895*14)

AHAS46denom<-110743263-(895*14)
AHAS51denom<-87211176-(895*14)
AHAS56denom<-119919534-(895*14)
AHAS61denom<-110418543-(895*14)

AHP01denom<-107082003-(895*14)
AHP06denom<-89088343-(895*14)
#AHP11FOURSAMPLEdenom<-95288287-(895*14)
AHP16denom<-103199261-(895*14)
AHP21denom<-112945213-(895*14)

AHE01denom<-155603504-(895*14)
#AHE05denom<-10367803-(895*14)
AHE07denom<-134287918-(895*14)
AHE09denom<-151875529-(895*14)
AHE11denom<-131785459-(895*14)
mean(c(AHB125denom, AHB145denom, AHB151denom,AHB176denom,AHB195denom, AHB70denom,AHB90denom,
       AHAS46denom, AHAS51denom, AHAS56denom, AHAS61denom, 
       AHP01denom, AHP06denom, AHP16denom, AHP21denom))
se(c(AHB125denom, AHB145denom, AHB151denom,AHB176denom,AHB195denom, AHB70denom,AHB90denom,
       AHAS46denom, AHAS51denom, AHAS56denom, AHAS61denom, 
       AHP01denom, AHP06denom, AHP16denom, AHP21denom))
#surface areas: 
data<-read.delim("~/Documents/BikiniPipeline/sizebymutationfreq20191201.txt")
bikini<-data[which(data$location=="Bikini"),]
ofu<-data[which(data$location=="Ofu"),]
ofusizes<-ofu$distcm
bikinisize<-bikini$distcm
palausizes<-c(55,139,67,420)
enewetaksizes<-c(69,64,81,108)  
soma_func<- function(files) {
  #files<-list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)
  metadata= NULL
  for (i in 1:length(files)) { 
    file =files[i]
    data<-read.delim(file) #read in each file in "files"
    data<-data.frame(data) # transform the data from each file into a dataframe
    base<-basename(file)
    colony<-strsplit(base, "\\_")[[1]][1]
    print(colony)
    len<-nrow(data) 
    colonyrep<-rep(colony, len)
    withcolony<-data.frame(data, colonyrep) #combines the colonyname column with each dataframe
    metadata <- rbind(metadata, withcolony) #adds each dataframe to the overall metatadata, so that the information from all of the files are now in "metadata"
  }
  
  genoanddepth<-(metadata$genotype) #names the column
  split<-str_split_fixed(genoanddepth, ",", 4) #split the genotype, depths, and GQ score into their own separate strings
  
  genotypes<-split[,1] #defines the genotype as the first string in "split"
  position<-metadata$chrom.pos #names the column
  positionsplit<-str_split_fixed(position, "[.]", 2) #split the chromosome number and the position on the chromosome into their own separate strings
  
  chr<-positionsplit[,1] #defines the chromosome as the first string in "positionsplit"
  pos<-positionsplit[,2] #defines the position as the second string in "positionsplit"
  #totaldepth<-as.numeric(split[,2])
  refdepth<-as.numeric(split[,2]) #number of reference reads
  altdepth<-as.numeric(split[,3]) #number of alternate reads
  what<-metadata$WhattoWhat
  allelesplit<-str_split_fixed(what, "to", 2) #splits the normal and mutant alleles into different strings
  normalallele<-allelesplit[,1]
  mutantallele<-allelesplit[,2]
  totaldepth<-refdepth+altdepth #sum of all reference and alternate reads
  GQscore<-as.numeric(split[,4])
  mutationtype<-metadata$MutationType #eg intergenic_region 3_prime_UTR_variant 5_prime_UTR_premature_start_codon_gain_variant 5_prime_UTR_variant
  mutationstrength<-metadata$MutationStrength #eg HIGH LOW MODERATE MODIFIER
  mutant_alleledepth = rep("A", nrow(metadata))
  for (i in 1:nrow(metadata)){
    if (mutantallele[i] == metadata$alt[i]) {
      mutant_alleledepth[i] = altdepth[i]
      #print(mutant_alleledepth[i])
    } else {
      mutant_alleledepth[i] = refdepth[i]
      #print(mutant_alleledepth[i])
    }  #print("ALT", mutant_alleledepth, refdepth)
  }
  #print(mutant_alleledepth[1:10])
  #print(refdepth[1:10])
  #print(altdepth[1:10])
  
  metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt,
                         "normal_allele"= normalallele, "mutant_allele" = mutantallele, "mutant_allele_depth" = as.numeric(mutant_alleledepth), 
                         "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	
                         "GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat, "MutantSample1"=metadata$MutantSample1, "MutantSample2"=metadata$MutantSample2,
                         "MutationType"=mutationtype, "MutationStrength"=mutationstrength, "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
  
  uniquesomaticmetadatadf<-metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
  
  i <- sapply(metadatadf, is.factor) #gives TRUE or FALSE for whether each column is a factor
  metadatadf[i] <- lapply(metadatadf[i], as.character) #sets columns to be characters
  metadatadf[which(metadatadf$sample==metadatadf$MutantParent1),]
  DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean) #calculate mean depth per locus
  
  DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min) #colculate min depth per locus
  
  GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean) #calculate average GQ perlocus
  
  GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min) #calculate min GQ per locus
  
  #return(uniquesomaticmetadatadf)
  goh<-subset(uniquesomaticmetadatadf, GoH_or_LoH=="DeNovo")
  loh<-subset(uniquesomaticmetadatadf, GoH_or_LoH=="LoH")
  #return(list("goh"=goh,"loh"=loh))
  metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
  metadatadf.1<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-subset(metadatadf.1, startsWith(metadatadf.1$chrom, 'chr')==TRUE)
  #metadatadf.0$MutantSample1 <- factor(metadatadf.0$MutantSample1, levels=levels(metadatadf.0$sample))
  #metadatadf.0$MutantSample2 <- factor(metadatadf.0$MutantSample2, levels=levels(metadatadf.0$sample))
  
  uniquesomaticmetadatadf<-metadatadf.0[match(unique(metadatadf.0$chrom.pos), 					metadatadf.0$chrom.pos),]
  #return(uniquesomaticmetadatadf) #do not use for comparing against enewetak
  #somatic denovos:
  DeNovos<-subset(metadatadf.0, GoH_or_LoH=="DeNovo")
  withoutmutsample<-subset(DeNovos, (sample != MutantSample1 & sample != MutantSample2))# & sample != MutantParent2)) #subset the samples that are not the mutantparents
  
  uniquetoonesample<-subset(withoutmutsample, refdepth == 0 | altdepth ==0) #subsets just the nonmutant parent samples that have either refdepth=0 or altdepth=0. this eliminates the nonmutants for which the putatively mutant allele is present in "nonmutants" at low frequency
  uniquetoonesample_freq<-as.data.frame(table(uniquetoonesample$chrom.pos)) #gives  how many nonmutant parent samples are present per locus once they have been subsetted
  all_homozyg<-subset(uniquetoonesample_freq, Freq==5)
  uniquetoonesample_fulldf<-DeNovos[match(all_homozyg$Var1,DeNovos$chrom.pos),]
  #return(uniquetoonesample_fulldf)
  
  #somatic loh
  LoHonly<-subset(metadatadf.0, GoH_or_LoH=="LoH")
  LOHsample1<-subset(LoHonly, sample == MutantSample1)# & sample != MutantSample2))
  LOHsample2<-subset(LoHonly, sample == MutantSample2)   
  trueLoH1<-subset(LOHsample1, refdepth =="0" | altdepth=="0") #subsets to just the MutantParent1s that have zero minor allele frequency
  trueLoH2<-subset(LOHsample2, refdepth =="0" | altdepth=="0") #subsets to just the MutantParent2s that have zero minor allele frequency
  trueLoHtotal<- rbind(trueLoH1, trueLoH2) #combines trueLoHp2 and trueLoHp2.1
  removeuniques_LoH<-subset(trueLoHtotal,duplicated(chrom.pos)==TRUE | duplicated(chrom.pos, fromLast=TRUE)==TRUE) #subsets just the loci that appear in both MutantParent1 and MutantParent2 from trueLoHp2 and trueLoHp2.1
  u_removeuniques_LoH<-removeuniques_LoH[match(unique(removeuniques_LoH$chrom.pos), 					removeuniques_LoH$chrom.pos),]
  return(list("GOH"=goh, "LOH"=loh)) #multisample list
  #return(list("singlesampleGOH"=uniquetoonesample_fulldf,"singlesampleLOH"=u_removeuniques_LoH))
  allsinglesample<-rbind(uniquetoonesample_fulldf,u_removeuniques_LoH) 
  #return(allsinglesample) #compare to enewetak
} 
AHE_soma_func<- function(files) {
  #files<-list.files(path="~/Documents/GitHub/Bikini", pattern="AHE01-02_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)
  metadata= NULL
  for (i in 1:length(files)) { 
    file =files[i]
    data<-read.delim(file) #read in each file in "files"
    data<-data.frame(data) # transform the data from each file into a dataframe
    base<-basename(file)
    colony<-strsplit(base, "\\_")[[1]][1]
    print(colony)
    len<-nrow(data) 
    colonyrep<-rep(colony, len)
    withcolony<-data.frame(data, colonyrep) #combines the colonyname column with each dataframe
    metadata <- rbind(metadata, withcolony) #adds each dataframe to the overall metatadata, so that the information from all of the files are now in "metadata"
  }
  
  genoanddepth<-(metadata$genotype) #names the column
  split<-str_split_fixed(genoanddepth, ",", 4) #split the genotype, depths, and GQ score into their own separate strings
  
  genotypes<-split[,1] #defines the genotype as the first string in "split"
  position<-metadata$chrom.pos #names the column
  positionsplit<-str_split_fixed(position, "[.]", 2) #split the chromosome number and the position on the chromosome into their own separate strings
  
  chr<-positionsplit[,1] #defines the chromosome as the first string in "positionsplit"
  pos<-positionsplit[,2] #defines the position as the second string in "positionsplit"
  #totaldepth<-as.numeric(split[,2])
  refdepth<-as.numeric(split[,2]) #number of reference reads
  altdepth<-as.numeric(split[,3]) #number of alternate reads
  #what<-metadata$WhattoWhat
  #allelesplit<-str_split_fixed(what, "to", 2) #splits the normal and mutant alleles into different strings
  #normalallele<-allelesplit[,1]
  #mutantallele<-allelesplit[,2]
  totaldepth<-refdepth+altdepth #sum of all reference and alternate reads
  GQscore<-as.numeric(split[,4])
  mutationtype<-metadata$MutationType #eg intergenic_region 3_prime_UTR_variant 5_prime_UTR_premature_start_codon_gain_variant 5_prime_UTR_variant
  mutationstrength<-metadata$MutationStrength #eg HIGH LOW MODERATE MODIFIER
  mutant_alleledepth = rep("A", nrow(metadata))
  #for (i in 1:nrow(metadata)){
  # if (mutantallele[i] == metadata$alt[i]) {
  #  mutant_alleledepth[i] = altdepth[i]
  #print(mutant_alleledepth[i])
  #} else {
  #  mutant_alleledepth[i] = refdepth[i]
  #print(mutant_alleledepth[i])
  #}  #print("ALT", mutant_alleledepth, refdepth)
  #}
  #print(mutant_alleledepth[1:10])
  #print(refdepth[1:10])
  #print(altdepth[1:10])
  
  metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, 
                         
                         "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	
                         
                         "MutationType"=mutationtype, "MutationStrength"=mutationstrength, "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
  i <- sapply(metadatadf, is.factor) #gives TRUE or FALSE for whether each column is a factor
  metadatadf[i] <- lapply(metadatadf[i], as.character) #sets columns to be characters
  metadatadf[which(metadatadf$sample==metadatadf$MutantParent1),]
  DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean) #calculate mean depth per locus
  
  DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min) #colculate min depth per locus
  
  GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean) #calculate average GQ perlocus
  
  GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min) #calculate min GQ per locus
  
  
  #return(list("goh"=goh,"loh"=loh))
  metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
  metadatadf.1<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-subset(metadatadf.1, startsWith(metadatadf.1$chrom, 'chr')==TRUE)
  
  uniquesomaticmetadatadf<-metadatadf.0[match(unique(metadatadf.0$chrom.pos), 					metadatadf$chrom.pos.0),]
  
  justhomozyg<-subset(metadatadf.0, genotype =="T/T" | genotype =="C/C" | genotype =="G/G" | genotype =="A/A")
  truejusthomozyg<-subset(justhomozyg, refdepth == 0 | altdepth ==0)
  truejusthomozyg_freq<-as.data.frame(table(truejusthomozyg$chrom.pos))
  all_homozyg<-subset(truejusthomozyg_freq, Freq==2)
  all_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  uniqueall_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  return(uniqueall_homozyg_fulldf)
  
}
#GOH vs LOH for multisamples
ahb125<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB125denom
ahb145<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB145-149_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB145denom
ahb151<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB151-155_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB151denom
ahb176<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB176-180_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB176denom
ahb195<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB195-199_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB195denom
ahb70<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB70-74_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB70denom
ahb90<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB90-94_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHB90denom
bikrates_GOH<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS46<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS46-50_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHAS46denom
AHAS51<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS51-55_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHAS51denom
AHAS56<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHAS56denom
AHAS61<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHAS61denom

ofurates_GOH<-c(AHAS46, AHAS51, AHAS56, AHAS61)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHP01denom
AHP06<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP06-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHP06denom
AHP16<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP16-20_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHP16denom
AHP21<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP21-25_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$GOH)/AHP21denom

palrates_GOH<-c(AHP01,AHP06,AHP16,AHP21)

ahb125<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB125denom
ahb145<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB145-149_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB145denom
ahb151<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB151-155_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB151denom
ahb176<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB176-180_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB176denom
ahb195<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB195-199_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB195denom
ahb70<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB70-74_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB70denom
ahb90<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB90-94_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHB90denom
bikrates_LOH<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS46<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS46-50_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHAS46denom
AHAS51<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS51-55_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHAS51denom
AHAS56<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHAS56denom
AHAS61<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHAS61denom

ofurates_LOH<-c(AHAS46, AHAS51, AHAS56, AHAS61)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHP01denom
AHP06<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP06-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHP06denom
AHP16<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP16-20_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHP16denom
AHP21<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP21-25_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))$LOH)/AHP21denom

palrates_LOH<-c(AHP01,AHP06,AHP16,AHP21)

gohrates<-c(bikrates_GOH, ofurates_GOH, palrates_GOH)
lohrates<-c(bikrates_LOH, ofurates_LOH, palrates_LOH)
percentloh<-lohrates/(gohrates+lohrates)
percentgoh<-gohrates/(gohrates+lohrates)

size<-c(bikinisize,ofusizes,palausizes)
location<-c(rep("Bikini",7),rep("Ofu",4),rep("Palau",4))
treatment<-c(rep("Irradiated",7),rep("Not Irradiated",8))
colonynames<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")
ratioavg<-c(mean(bikrates_GOH/bikrates_LOH))
AvgDepth<-c(AHB70avg, AHB90avg, AHB125avg, AHB145avg, AHB151avg, AHB176avg, AHB195avg,
            AHAS46avg, AHAS51avg, AHAS56avg, AHAS61avg, AHP01avg, AHP06avg, AHP16avg, AHP21avg)
typeframe<-data.frame(ratio,size,location, colonynames, lohrates, gohrates, treatment, AvgDepth, percentloh, percentgoh)
ratiobox<-ggplot(typeframe, aes(x=location, y=ratio, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
percentlohbox<-ggplot(typeframe, aes(x=location, y=percentloh, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
gohbox<-ggplot(typeframe, aes(x=location, y=gohrates, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
lohbox<-ggplot(typeframe, aes(x=location, y=lohrates, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
ratiobox | gohbox | lohbox
ratiofit<-lm(ratio~treatment)
anova(ratiofit)
ratiopoint<-ggplot(typeframe, (aes(x=size, y=ratio, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  theme_bw()
lohpoint<-ggplot(typeframe, (aes(x=size, y=lohrates, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  theme_bw()
#vdepth:
gohvdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=gohrates, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(gohrates)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()
lohvdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=lohrates, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(lohrates)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()
ratiovdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=ratio, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(AvgDepth)*1.1,label.y=max(ratio)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()
percentlohvdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=(percentloh)*100, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(percentloh)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()
#this is fig:
gohvdepth | lohvdepth | percentlohvdepth

gohvsloh<-ggplot(typeframe, (aes(x=gohrates, y=lohrates, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(gohrates)*1.1,label.y=max(lohrates)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()
gohvdepth | lohvdepth | ratiovdepth
fit_goh<-lm(gohrates~ AvgDepth, data=typeframe)
fit_loh<-lm(lohrates~ AvgDepth, data=typeframe)
fit_percentloh<-lm(percentloh~ AvgDepth, data=typeframe)

typeframe$goh.resid<-fit_goh$residuals
typeframe$loh.resid<-fit_loh$residuals
typeframe$lohpercent.resid<-fit_percentloh$residuals
#residuals v location:
percentlohbox.resid<-ggplot(typeframe, aes(x=location, y=loh.resid, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
#residuals v size
lohpoint.resid<-ggplot(typeframe, (aes(x=size, y=loh.resid, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(size)*0.5,label.y=max(fit_loh$residuals)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1))+
  theme_bw()
gohpoint.resid<-ggplot(typeframe, (aes(x=size, y=goh.resid, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(size)*0.5,label.y=max(fit_goh$residuals)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1))+
  theme_bw()
lohpercent.resid<-ggplot(typeframe, (aes(x=size, y=lohpercent.resid, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(size)*0.5,label.y=max(fit_percentloh$residuals)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1))+
  theme_bw()
#this is fig:
gohpoint.resid | lohpoint.resid | lohpercent.resid
#resids boxplots
percentlohbox.resid<-ggplot(typeframe, aes(x=location, y=lohpercent.resid, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
lohbox.resid<-ggplot(typeframe, aes(x=location, y=loh.resid, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
gohbox.resid<-ggplot(typeframe, aes(x=location, y=goh.resid, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
gohbox.resid | lohbox.resid | percentlohbox.resid
##COMPARE AGAINST ENEWETAK: return(uniqueall_homozyg_fulldf) return(allsinglesample)
AHE01<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHE01-02_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE01denom
AHE05<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHE05-06_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE05denom
AHE07<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHE07-08_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE07denom
AHE09<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHE09-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE09denom
AHE11<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHE11-12_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE11denom

enerates<-c(AHE01,AHE07,AHE09,AHE11)

ahb125<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB125denom
ahb145<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB145-149_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB145denom
ahb151<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB151-155_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB151denom
ahb176<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB176-180_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB176denom
ahb195<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB195-199_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB195denom
ahb70<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB70-74_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB70denom
ahb90<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB90-94_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB90denom
bikrates<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS46<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS46-50_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS46denom
AHAS51<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS51-55_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS51denom
AHAS56<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS56denom
AHAS61<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS61denom

ofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP01denom
AHP06<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP06-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP06denom
AHP16<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP16-20_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP16denom
AHP21<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP21-25_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP21denom

palrates<-c(AHP01,AHP06,AHP16,AHP21)

irradiated<-c(enerates,bikrates)
notirradiated<-c(ofurates,palrates)

averages<-c(mean(bikrates),mean(enerates),mean(ofurates),mean(palrates))
rates<-c(bikrates, enerates, ofurates,palrates)
ses<-c(se(bikrates),se(enerates),se(ofurates),se(palrates))
location<-c(rep("Bikini",7),rep("Enewetak",4),rep("Ofu",4),rep("Palau",4))
sizes<-c(bikinisize,enewetaksizes,ofusizes,palausizes)
colony<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "ahe01","ahe07","ahe09","ahe11","AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")

ratedf<-data.frame("Rate"=rates, "Location"=location, "Distance"=sizes,"Colony"=colony)
#rate vs location:
boxun<-ggplot(ratedf, aes(x=Location, y=Rate))+
  geom_boxplot()+
  theme_bw()
#ratevs size for muts UNIQUE TO ONE SAMPLE:
un<-ggplot(ratedf, (aes(x=Distance, y=Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()
#rate vs size for muts NOT uniqueto one sample: return(uniquesomatic) #do not include enewetak
ahb125<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB125denom
ahb145<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB145-149_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB145denom
ahb151<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB151-155_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB151denom
ahb176<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB176-180_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB176denom
ahb195<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB195-199_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB195denom
ahb70<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB70-74_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB70denom
ahb90<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB90-94_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB90denom
bikrates_allsom<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS46<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS46-50_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS46denom
AHAS51<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS51-55_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS51denom
AHAS56<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS56denom
AHAS61<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS61denom

ofurates_allsom<-c(AHAS46, AHAS51, AHAS56, AHAS61)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP01denom
AHP06<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP06-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP06denom
AHP16<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP16-20_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP16denom
AHP21<-nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP21-25_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP21denom

palrates_allsom<-c(AHP01,AHP06,AHP16,AHP21)
rates_all<-c(bikrates_allsom, ofurates_allsom, palrates_allsom)
size<-c(bikinisize,ofusizes,palausizes)
location<-c(rep("Bikini",7),rep("Ofu",4),rep("Palau",4))
colonynames<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")
rate_alldf<-data.frame(rates_all,size,location,colonynames)
all<-ggplot(rate_alldf, (aes(x=size, y=rates_all, group=location, colour=location)))+
  #geom_point(size=3)+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  theme_bw()
#1a.	What is the average mutation frequency (differences between the two most distant samples from a colony) for each population? Is there a difference between Bikini/Enewetak and Palau/Ofu?

#1b how many from #1a are unique (only seen in one sample)?
#return(list("singlesampleGOH"=uniquetoonesample_fulldf,"singlesampleLOH"=u_removeuniques_LoH))

#spectra
spectrum<- function(dataframe, dataframecount) {
  #Amil_fasta<-read.fasta("~/Documents/GitHub/CoralGermline/dndscv/Amil.v2.01.chrs.fasta",forceDNAtolower = FALSE)
  WhattoWhat<-rep("Placeholder",nrow(dataframe))
  TiorTv<-rep("Placeholder",nrow(dataframe))
  
  ref<-as.character(dataframe$ref)
  alt<-as.character(dataframe$alt)
  for (i in 1:nrow(dataframe)) {
    if (ref[i] == "C" & alt[i] == "A") {
      WhattoWhat[i] <- 'CtoA'
      TiorTv[i]<-"Tv"
      #gh[i] <- 'CtoA'
      #print("ya")
      #print(i)
      #print(AtG[i])
    } else if (ref[i] == "A" & alt[i] == "C")  {
      WhattoWhat[i] <- "AtoC"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "A" & alt[i] == "G") {
      WhattoWhat[i] <- "AtoG"
      TiorTv[i]<-"Ti"
      
    } else if (ref[i] == "G" & alt[i] == "A") {
      WhattoWhat[i] <- 'GtoA'
      TiorTv[i]<-"Ti"
    } else if (ref[i] == "G" & alt[i] == "T") {
      WhattoWhat[i] <- "GtoT"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "T" & alt[i] == "G") {
      WhattoWhat[i] <- "TtoG"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "T" & alt[i] == "C") {
      WhattoWhat[i] <- "TtoC"
      TiorTv[i]<-"Ti"
    } else if (ref[i] == "C" & alt[i] == "T") {
      WhattoWhat[i] <-"CtoT"
      TiorTv[i]<-"Ti"
    } else if (ref[i] == "C" & alt[i] == "G") {
      WhattoWhat[i] <- "CtoG"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "G" & alt[i] == "C") {
      WhattoWhat[i] <-"GtoC"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "T" & alt[i] == "A") {
      WhattoWhat[i] <- "TtoA"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "A" & alt[i] == "T") {
      WhattoWhat[i] <-"AtoT"
      TiorTv[i]<-"Tv"
    }
  }
  
  together<-cbind(dataframe, WhattoWhat,TiorTv)				 
  togethercount<-nrow(together)
  
  verATGClist<-together[ which(together$WhattoWhat=="AtoG" | together$WhattoWhat=="TtoC"),]
  verATGCcount<-nrow(verATGClist)
  verATGCprop<-verATGCcount/togethercount
  SEverATGCprop<-se(verATGCprop)
  
  verGCATlist<-together[ which(together$WhattoWhat=="GtoA" | together$WhattoWhat=="CtoT"),]
  verGCATcount<-nrow(verGCATlist)
  verGCATprop<-verGCATcount/togethercount
  
  CpGdf<-data.frame(Date=as.Date(character()))
  CpHdf<-data.frame(Date=as.Date(character()))
  all<-NULL
  for (i in 1:nrow(together)) {
    if (together$WhattoWhat[i] == "CtoT") {
      character<-as.character(together$pos[i])
      integer<-as.integer(character)
      next_position<-integer + 1
      chr_num<-as.character(together$chrom[i])
      
      chroms<-Amil_fasta[[chr_num]]
      ref_nuc<- chroms[next_position]
      ref_nuc<-toupper(ref_nuc)
      allnucs<-data.frame(ref_nuc)
      all<-rbind(all,allnucs)
      if (ref_nuc == "G") {
        #print(ref_nuc)
        Gnucs<-data.frame(ref_nuc)
        CpGdf<-rbind(CpGdf, Gnucs)
      } else if (ref_nuc == "A"| ref_nuc == "C"| ref_nuc == "T") {
        othernucs<-data.frame(ref_nuc)
        CpHdf<-rbind(CpHdf,othernucs)
      }
    } else if (together$WhattoWhat[i] == "GtoA") {   
      character<-as.character(together$pos[i])
      integer<-as.integer(character)
      next_position<-integer + 1
      chr_num<-as.character(together$chrom[i])
      
      chroms<-Amil_fasta[[chr_num]]
      ref_nuc<- chroms[next_position]
      ref_nuc<-toupper(ref_nuc)
      allnucs<-data.frame(ref_nuc)
      all<-rbind(all,allnucs)
      
      if (ref_nuc == "C") {
        #print(ref_nuc)
        Gnucs<-data.frame(ref_nuc)
        CpGdf<-rbind(CpGdf, Gnucs)
      } else if (ref_nuc == "A"| ref_nuc == "G"| ref_nuc == "T") {
        othernucs<-data.frame(ref_nuc)
        CpHdf<-rbind(CpHdf,othernucs)
      }
    }
  }    
  CpGcount<-nrow(CpGdf)
  CpHcount<-nrow(CpHdf)
  
  CpGprop<- CpGcount/togethercount
  CpHprop<- CpHcount/togethercount
  
  verATTAlist<-together[ which(together$WhattoWhat=="AtoT" | together$WhattoWhat=="TtoA"),]
  verATTAcount<-nrow(verATTAlist)
  verATTAprop<-verATTAcount/togethercount
  
  verACTGlist<-together[ which(together$WhattoWhat=="AtoC" | together$WhattoWhat=="TtoG"),]
  verACTGcount<-nrow(verACTGlist)
  verACTGprop<-verACTGcount/togethercount
  
  verGCCGlist<-together[ which(together$WhattoWhat=="CtoG" | together$WhattoWhat=="GtoC"),]
  verGCCGcount<-nrow(verGCCGlist)
  verGCCGprop<-verGCCGcount/togethercount
  
  verGCTAlist<-together[ which(together$WhattoWhat=="GtoT" | together$WhattoWhat=="CtoA"),]
  verGCTAcount<-nrow(verGCTAlist)
  verGCTAprop<-verGCTAcount/togethercount
  
  binomATGC<-dbinom(verATGCcount, togethercount, verATGCprop)
  sdATGC<-togethercount*verATGCprop*(1-verATGCprop)
  seATGC<- sdATGC/sqrt(togethercount)
  
  binomGCAT<-dbinom(verGCATcount, togethercount, verGCATprop)
  sdGCAT<-togethercount*verGCATprop*(1-verGCATprop)
  seGCAT<- sdGCAT/sqrt(togethercount)
  
  binomCpG<-dbinom(CpGcount, togethercount, CpGprop)
  sdCpG<-togethercount*CpGprop*(1-CpGprop)
  seCpG<-sdCpG/sqrt(togethercount)
  
  binomCpH<-dbinom(CpHcount, togethercount, CpHprop)
  sdCpH<-togethercount*CpHprop*(1-CpHprop)
  seCpH<-sdCpH/sqrt(togethercount)
  
  binomATTA<-dbinom(verATTAcount, togethercount, verATTAprop)
  sdATTA<-togethercount*verATTAprop*(1-verATTAprop)
  seATTA<- sdATTA/sqrt(togethercount)
  
  binomACTG<-dbinom(verACTGcount, togethercount, verACTGprop)
  sdACTG<-togethercount*verACTGprop*(1-verACTGprop)
  seACTG<- sdACTG/sqrt(togethercount)
  
  binomGCCG<-dbinom(verGCCGcount, togethercount, verGCCGprop)
  sdGCCG<-togethercount*verGCCGprop*(1-verGCCGprop)
  seGCCG<- sdGCCG/sqrt(togethercount)
  
  binomGCTA<-dbinom(verGCTAcount, togethercount, verGCTAprop)
  sdGCTA<-togethercount*verGCTAprop*(1-verGCTAprop)
  seGCTA<- sdGCTA/sqrt(togethercount)
  
  standarderrors<-c(seATGC, seCpG, seCpH, seATTA, seACTG, seGCCG, seGCTA)
  se_forprops<-standarderrors/togethercount
  vertypesDF<-data.frame(Types=c("A>G/T>C","CpG", "CpH", "A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), coralProportion=c(verATGCprop, CpGprop, CpHprop, verATTAprop,verACTGprop, verGCCGprop, verGCTAprop),se_forprops)
  #pmdf$MutType<- factor(pmdf$MutType,levels=c("synonymous","missense","nonsense"))
  vertypesDF$Types<-factor(vertypesDF$Types,levels=c("A>G/T>C","CpG", "CpH", "A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"))
  #MUTATION SPECTRUM PLOT:
  #vertypesDFplot<-barplot(vertypesDF$coralProportion,names.arg=vertypesDF$Types, ylim=c(0,0.7), main="allmuts", las=2) #Figure 4a
  z<- ggplot(vertypesDF, aes(x=Types, y=coralProportion)) +
    #scale_fill_manual(value=c("white"))+
    geom_bar(stat="identity") +
    ylim(0, 0.35)+
    theme_bw()+
    ylab("Proportion of SNVs") +
    theme(axis.text.x = element_text(angle = 90, size=18), axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+#, axis.text.x=element_text(size=25))+
    geom_errorbar(aes(ymin=coralProportion-(se_forprops), ymax=coralProportion+(se_forprops)))#, width=.2,
  
  Transtions<- verATGCprop + verGCATprop
  Transversions<- verATTAprop+verACTGprop+verGCCGprop+verGCTAprop
  TiTv<- Transtions/Transversions
  counts<-c(verATGCcount,CpGcount, CpHcount, verATTAcount, verACTGcount, verGCCGcount, verGCTAcount)
  #return(z)
  #return(list(vertypesDF, TiTv))
  #return(list(z, vertypesDF))
  #return(z)
  return(counts)
}

AHAS61spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))

ahb125spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB125-129_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahb145spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB145-149_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB145-149_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahb151spec<-spectrum( soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB151-155_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB151-155_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahb176spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB176-180_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB176-180_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahb195spec<-spectrum( soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB195-199_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB195-199_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahb70spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB70-74_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB70-74_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahb90spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB90-94_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHB90-94_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
c1<-mean(c(ahb125spec[1], ahb145spec[1], ahb151spec[1], ahb176spec[1], ahb195spec[1],ahb70spec[1],ahb90spec[1]))
c2<-mean(c(ahb125spec[2], ahb145spec[2], ahb151spec[2], ahb176spec[2], ahb195spec[2],ahb70spec[2],ahb90spec[2]))
c3<-mean(c(ahb125spec[3], ahb145spec[3], ahb151spec[3], ahb176spec[3], ahb195spec[3],ahb70spec[3],ahb90spec[3]))
c4<-mean(c(ahb125spec[4], ahb145spec[4], ahb151spec[4], ahb176spec[4], ahb195spec[4],ahb70spec[4],ahb90spec[4]))
c5<-mean(c(ahb125spec[5], ahb145spec[5], ahb151spec[5], ahb176spec[5], ahb195spec[5],ahb70spec[5],ahb90spec[5]))
c6<-mean(c(ahb125spec[6], ahb145spec[6], ahb151spec[6], ahb176spec[6], ahb195spec[6],ahb70spec[6],ahb90spec[6]))
c7<-mean(c(ahb125spec[7], ahb145spec[7], ahb151spec[7], ahb176spec[7], ahb195spec[7],ahb70spec[7],ahb90spec[7]))
ahbmeancounts<-c(c1,c2,c3,c4,c5,c6,c7)

ahas46spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS46-50_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS46-50_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahas51spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS51-55_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS51-55_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahas56spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahas61spec<-spectrum( soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS61-65_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHAS56-60_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
c1<- mean(c(ahas46spec[1], ahas51spec[1], ahas56spec[1],ahas61spec[1]))
c2<- mean(c(ahas46spec[2], ahas51spec[2], ahas56spec[2],ahas61spec[2]))
c3<- mean(c(ahas46spec[3], ahas51spec[3], ahas56spec[3],ahas61spec[3]))
c4<- mean(c(ahas46spec[4], ahas51spec[4], ahas56spec[4],ahas61spec[4]))
c5<- mean(c(ahas46spec[5], ahas51spec[5], ahas56spec[5],ahas61spec[5]))
c6<- mean(c(ahas46spec[6], ahas51spec[6], ahas56spec[6],ahas61spec[6]))
c7<- mean(c(ahas46spec[7], ahas51spec[7], ahas56spec[7],ahas61spec[7]))
ahasmeancounts<-c(c1,c2,c3,c4,c5,c6,c7)

ahp01spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahp06spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP06-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahp16spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP16-20_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
ahp21spec<-spectrum(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP21-25_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)), nrow(soma_func(list.files(path="~/Documents/GitHub/Bikini", pattern="AHP01-05_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE))))
c1<- mean(c(ahp01spec[1], ahp06spec[1], ahp16spec[1],ahp21spec[1]))
c2<- mean(c(ahp01spec[2], ahp06spec[2], ahp16spec[2],ahp21spec[2]))
c3<- mean(c(ahp01spec[3], ahp06spec[3], ahp16spec[3],ahp21spec[3]))
c4<- mean(c(ahp01spec[4], ahp06spec[4], ahp16spec[4],ahp21spec[4]))
c5<- mean(c(ahp01spec[5], ahp06spec[5], ahp16spec[5],ahp21spec[5]))
c6<- mean(c(ahp01spec[6], ahp06spec[6], ahp16spec[6],ahp21spec[6]))
c7<- mean(c(ahp01spec[7], ahp06spec[7], ahp16spec[7],ahp21spec[7]))
ahpmeancounts<-c(c1,c2,c3,c4,c5,c6,c7)
chisq.test(ahbmeancounts,ahasmeancounts)
chisq.test(ahbmeancounts,ahpmeancounts)
chisq.test(ahasmeancounts,ahpmeancounts)
