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
files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/", pattern="*svPOP.txt", full.names=T, recursive=FALSE)
files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/", pattern="AHAS56-61_mf_svPOP.txt", full.names=T, recursive=FALSE)

for (i in 1:length(files)) {
  #meltfunc(i)
  file =files[i]
  #print(file)
  meltfunc(file)
}
# SEVEN SAMPLE denominators:
AHB125denom<-168738945-(949*14)
AHB145denom<-156842817-(949*14)
AHB151denom<-159349164-(949*14)
AHB176denom<-158135967-(949*14)
AHB195denom<-167199712-(949*14)
AHB70denom<-162467953-(949*14)
AHB90denom<-165400558-(949*14)

AHAS46denom<-153597774-(949*14)
AHAS51denom<-103151602-(949*14)
AHAS56denom<-147623534-(949*14)
AHAS61denom<-155221674-(949*14)

AHP01denom<-179162251-(949*14)
AHP06denom<-121528791-(949*14)
#AHP11FOURSAMPLEdenom<-95288287-(949*14)
AHP16denom<-180257116-(949*14)
#AHP21denom<-145573426-(949*14)

AHE01denom<-201273381-(949*14)
#AHE05denom<-10367803-(949*14)
AHE07denom<-179021266-(949*14)
AHE09denom<-198095283-(949*14)
AHE11denom<-167479025-(949*14)
sevensampledenomaverage<-mean(c(AHB125denom, AHB145denom, AHB151denom, AHB176denom, AHB195denom, AHB70denom, AHB90denom,
                                AHAS46denom, AHAS51denom, AHAS56denom, AHAS61denom,
                                AHP01denom, AHP06denom, AHP16denom))
se(c(AHB125denom, AHB145denom, AHB151denom, AHB176denom, AHB195denom, AHB70denom, AHB90denom,
          AHAS46denom, AHAS51denom, AHAS56denom, AHAS61denom,
          AHP01denom, AHP06denom, AHP16denom))
  #FOURSAMPLEdenoms
AHAS41FOURSAMPLEdenom<-189029139-(949*14)
AHAS46FOURSAMPLEdenom<-191386179-(949*14)
AHAS51FOURSAMPLEdenom<-179872303-(949*14)
AHAS56FOURSAMPLEdenom<-198843433-(949*14)
AHAS61FOURSAMPLEdenom<-186053791-(949*14)

AHB125FOURSAMPLEdenom<-187709580-(949*14)
AHB145FOURSAMPLEdenom<-171350522-(949*14)
AHB151FOURSAMPLEdenom<-171020741-(949*14)
AHB176FOURSAMPLEdenom<-179562285-(949*14)
AHB195FOURSAMPLEdenom<-184848275-(949*14)
AHB70FOURSAMPLEdenom<-182671963-(949*14)
AHB90FOURSAMPLEdenom<-170777305-(949*14)

AHP01FOURSAMPLEdenom<-191803202-(949*14)
AHP06FOURSAMPLEdenom<-175377087-(949*14)
AHP16FOURSAMPLEdenom<-192511356-(949*14)
AHP21FOURSAMPLEdenom<-158795991-(949*14)

AHE01FOURSAMPLEdenom<-AHE01denom
AHE07FOURSAMPLEdenom<-AHE07denom
AHE09FOURSAMPLEdenom<-AHE09denom
AHE11FOURSAMPLEdenom<-AHE11denom

foursampledenomaverage<-mean(c(AHB125FOURSAMPLEdenom, AHB145FOURSAMPLEdenom, AHB151FOURSAMPLEdenom, AHB176FOURSAMPLEdenom, AHB195FOURSAMPLEdenom, AHB70FOURSAMPLEdenom, AHB90FOURSAMPLEdenom,
                                AHAS41FOURSAMPLEdenom, AHAS46FOURSAMPLEdenom, AHAS51FOURSAMPLEdenom, AHAS56FOURSAMPLEdenom, AHAS61FOURSAMPLEdenom,
                                AHP01FOURSAMPLEdenom, AHP06FOURSAMPLEdenom, AHP16FOURSAMPLEdenom,
                                AHE01FOURSAMPLEdenom, AHE07FOURSAMPLEdenom, AHE09FOURSAMPLEdenom, AHE11FOURSAMPLEdenom))
se(c(AHB125FOURSAMPLEdenom, AHB145FOURSAMPLEdenom, AHB151FOURSAMPLEdenom, AHB176FOURSAMPLEdenom, AHB195FOURSAMPLEdenom, AHB70FOURSAMPLEdenom, AHB90FOURSAMPLEdenom,
     AHAS41FOURSAMPLEdenom, AHAS46FOURSAMPLEdenom, AHAS51FOURSAMPLEdenom, AHAS56FOURSAMPLEdenom, AHAS61FOURSAMPLEdenom,
     AHP01FOURSAMPLEdenom, AHP06FOURSAMPLEdenom, AHP16FOURSAMPLEdenom,
     AHE01FOURSAMPLEdenom, AHE07FOURSAMPLEdenom, AHE09FOURSAMPLEdenom, AHE11FOURSAMPLEdenom))
#pop denoms:
AHB125_145denom<-137219243-(949*14)
AHB125_151denom<-139303183-(949*14)
AHB125_176denom<-156341294-(949*14)
AHB125_195denom<-155887883-(949*14)
AHB125_70denom<-157333067-(949*14)
AHB125_90denom<-153343052-(949*14)
AHB145_151denom<-125861592-(949*14)
AHB145_176denom<-141957728-(949*14)
AHB145_195denom<-141083191-(949*14)
AHB145_70denom<-142181478-(949*14)
AHB145_90denom<-139198956-(949*14)
AHB151_176denom<-127012695-(949*14)
AHB151_195denom<-126717716-(949*14)
AHB151_70denom<-127377724-(949*14)
AHB151_90denom<-124132191-(949*14)
AHB176_195denom<-151893837-(949*14)
AHB176_70denom<-152893560-(949*14)
AHB176_90denom<-148600540-(949*14)
AHB195_70denom<-156759000-(949*14)
AHB195_90denom<-152681837-(949*14)
AHB70_90denom<-153545402-(949*14)

AHAS41_46denom<-160974368-(949*14)
AHAS41_51denom<-151557577-(949*14)
AHAS41_56denom<-156712982-(949*14)
AHAS41_61denom<-148479253-(949*14)
AHAS46_51denom<-143969188-(949*14)
AHAS46_56denom<-157701115-(949*14)
AHAS46_61denom<-150732543-(949*14)
AHAS51_56denom<-144770085-(949*14)
AHAS51_61denom<-136488309-(949*14)
AHAS56_61denom<-183096388-(949*14)

AHE01_07denom<-169071677-(949*14)
AHE01_09denom<-182067405-(949*14)
AHE01_11denom<-158749609-(949*14)
AHE07_09denom<-166304695-(949*14)
AHE07_11denom<-149653591-(949*14)
AHE09_11denom<-154905420-(949*14)

AHP01_06denom<-147234306-(949*14)
AHP01_16denom<-160271792-(949*14)
AHP01_21denom<-135455632-(949*14)
AHP06_16denom<-165552705-(949*14)
AHP06_21denom<-147619237-(949*14)
AHP16_21denom<-147087908-(949*14)
#surface areas: 
data<-read.delim("~/Documents/BikiniPipeline/sizebymutationfreq20191201.txt")
bikini<-data[which(data$location=="Bikini"),]
ofu<-data[which(data$location=="Ofu"),]
ofusizes<-ofu$distcm
bikinisize<-bikini$distcm
palausizes<-c(55,139,67)
enewetaksizes<-c(69,64,81,108)  
#ahb125<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya", pattern="AHB70-74_mT_snps.txt", full.names=T, recursive=FALSE))$GOH#/AHBdenom
#files<-"AHB70-74_mT_snps.txt"
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
  #mutationtype<-metadata$MutationType #eg intergenic_region 3_prime_UTR_variant 5_prime_UTR_premature_start_codon_gain_variant 5_prime_UTR_variant
  #mutationstrength<-metadata$MutationStrength #eg HIGH LOW MODERATE MODIFIER
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
                        "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
  
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
  return(list("GOH"=nrow(goh), "LOH"=nrow(loh))) #multisample list
  #return(list("singlesampleGOH"=uniquetoonesample_fulldf,"singlesampleLOH"=u_removeuniques_LoH))
  allsinglesample<-rbind(uniquetoonesample_fulldf,u_removeuniques_LoH) 
  #return(allsinglesample) #compare to enewetak
} 

AHE_soma_func<- function(files) {
  #files<-list.files(path="~/Documents/GitHub/Bikini/mappedtoahya", pattern="AHB125-129_mf_FOURsnps.txt_melted.txt", full.names=T, recursive=FALSE)
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
  #mutationtype<-metadata$MutationType #eg intergenic_region 3_prime_UTR_variant 5_prime_UTR_premature_start_codon_gain_variant 5_prime_UTR_variant
  #mutationstrength<-metadata$MutationStrength #eg HIGH LOW MODERATE MODIFIER
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
                         
                         "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
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
  metadatadf.0<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  #metadatadf.0<-subset(metadatadf.1, startsWith(metadatadf.1$chrom, 'chr')==TRUE)
  
  uniquesomaticmetadatadf<-metadatadf.0[match(unique(metadatadf.0$chrom.pos), 					metadatadf$chrom.pos.0),]
  
  justhomozyg<-subset(metadatadf.0, genotype =="T/T" | genotype =="C/C" | genotype =="G/G" | genotype =="A/A")
  truejusthomozyg<-subset(justhomozyg, refdepth == 0 | altdepth ==0)
  truejusthomozyg_freq<-as.data.frame(table(truejusthomozyg$chrom.pos))
  all_homozyg<-subset(truejusthomozyg_freq, Freq==2)
  all_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  uniqueall_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  ###return this:::
  return(nrow(uniqueall_homozyg_fulldf))
  
  somdf=NULL
  for (i in 1:nrow(all_homozyg_fulldf)){
    submdf<-subset(metadatadf.0, chrom.pos== all_homozyg_fulldf$chrom.pos[i])
    somdf<- rbind(somdf, submdf)
  }
  #return(somdf)
  #VAFplot:
  somaticmetadatadf<-somdf
  
  justhet<-subset(somaticmetadatadf, genotype !="A/A" & genotype !="T/T" & genotype !="C/C" & genotype !="G/G") 
  orderedjusthet<- justhet[order(justhet$sample),]
  justhet$sample<-as.factor(justhet$sample)
  #sample_levels<-c(levels(justhet$sample))
  data<-as.data.frame(table(orderedjusthet$sample))
  orderedjusthet$generic<-c(rep("sample1",data$Freq[1]),rep("sample2",data$Freq[2]), rep("sample3",data$Freq[3]), rep("sample4",data$Freq[4]))
  tx <- split(orderedjusthet, orderedjusthet$generic)
  #one<-noquote(paste("`",sample_levels[1],"`", sep=""))
  sample1<-tx$sample1#`AHE01-1_S1` #sample_levels[1]
  sample2<-tx$sample2
  sample3<-tx$sample3
  sample4<-tx$sample4
  props1<- c(sample1$altdepth/sample1$totaldepth.x,  sample3$altdepth/sample3$totaldepth.x)
  props2<- c(sample2$altdepth/sample2$totaldepth.x, sample4$altdepth/sample4$totaldepth.x)
  VAFavg<- (props1+props2)/2
  #VAFavg<-as.data.frame(VAFavg)
  td1<-c(sample1$totaldepth.x,  sample3$totaldepth.x)
  td2<-c(sample2$totaldepth.x,  sample4$totaldepth.x)
  average_depth_at_locus<-(td1+td2)/2
  VAFvsdepth<-data.frame(VAFavg, average_depth_at_locus)
  
  vvd<-ggplot(data=VAFvsdepth, aes(y=VAFavg, x=average_depth_at_locus))+
    geom_point()+
    theme_bw()
  VAFavg<-as.data.frame(VAFavg)
  p<-ggplot(data=VAFavg, aes(x=VAFavg, fill="red"))+
    geom_density(alpha=0.5)+
    theme_bw()
  #return(vvd | p)
}
SV_soma_func<- function(file) {
  #file<-"AHE01-02_mf_sv.txt_melted.txt"
  #file="AHP01-06_mT_svPOP.ann.txt"
  metadata= NULL
  
  data<-read.delim(file) #read in each file in "files"
  data<-data.frame(data) # transform the data from each file into a dataframe
  base<-basename(file)
  colony<-strsplit(base, "\\_")[[1]][1]
  print(colony)
  len<-nrow(data) 
  colonyrep<-rep(colony, len)
  withcolony<-data.frame(data, colonyrep) #combines the colonyname column with each dataframe
  metadata <- rbind(metadata, withcolony) #adds each dataframe to the overall metatadata, so that the information from all of the files are now in "metadata"
  
  
  genotyp<-(metadata$genotype) #names the column
  split<-str_split_fixed(genotyp, ",", 2) #split the genotype, depths, and GQ score into their own separate strings
  
  genotypes<-split[,1] #defines the genotype as the first string in "split"
  mantatype<-split[,2]
  position<-metadata$chrom.pos #names the column
  positionsplit<-str_split_fixed(position, "[.]", 2) #split the chromosome number and the position on the chromosome into their own separate strings
  
  chr<-positionsplit[,1] #defines the chromosome as the first string in "positionsplit"
  pos<-positionsplit[,2] #defines the position as the second string in "positionsplit"
  #totaldepth<-as.numeric(split[,2])
  #refdepth<-as.numeric(split[,2]) #number of reference reads
  #altdepth<-as.numeric(split[,3]) #number of alternate reads
  #what<-metadata$WhattoWhat
  allelesplit<-str_split_fixed(genotyp, "/", 2) #splits the normal and mutant alleles into different strings
  allele1<-allelesplit[,1]
  allele2<-allelesplit[,2]
  #totaldepth<-refdepth+altdepth #sum of all reference and alternate reads
  #GQscore<-as.numeric(split[,4])
  #mutationtype<-metadata$MutationType #eg intergenic_region 3_prime_UTR_variant 5_prime_UTR_premature_start_codon_gain_variant 5_prime_UTR_variant
  #mutationstrength<-metadata$MutationStrength #eg HIGH LOW MODERATE MODIFIER
  #mutant_alleledepth = rep("A", nrow(metadata))
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
                         
                         "genotype"= genotypes,"allele1"=allele1, "allele2"=allele2, "MantaType"=mantatype,	
                         
                        "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
  i <- sapply(metadatadf, is.factor) #gives TRUE or FALSE for whether each column is a factor
  metadatadf[i] <- lapply(metadatadf[i], as.character) #sets columns to be characters
  metadatadf[which(metadatadf$sample==metadatadf$MutantParent1),]
  #DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean) #calculate mean depth per locus
  
  #DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min) #colculate min depth per locus
  
  #GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean) #calculate average GQ perlocus
  
  #GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min) #calculate min GQ per locus
  
  
  #return(list("goh"=goh,"loh"=loh))
  #metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
  #metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
  #metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
  #metadatadf.1<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  #metadatadf.0<-subset(metadatadf.1, startsWith(metadatadf.1$chrom, 'chr')==TRUE)
  
  uniquesomaticmetadatadf<-metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
  
  #justhomozyg<-subset(metadatadf, allele1==allele2)
  #truejusthomozyg<-subset(justhomozyg, refdepth == 0 | altdepth ==0)
  #truejusthomozyg_freq<-as.data.frame(table(truejusthomozyg$chrom.pos))
  #all_homozyg<-subset(truejusthomozyg_freq, Freq==2)
  #all_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  #uniqueall_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  #return(uniquesomaticmetadatadf)
  deletions_insertions<-subset(uniquesomaticmetadatadf, MantaType=="MantaDEL" | MantaType=="MantaINS")
  duplications<-subset(uniquesomaticmetadatadf, MantaType=="MantaDUP")
  translocations<-subset(uniquesomaticmetadatadf, MantaType=="MantaBND")
  allmantatypes<-list("Indels"=nrow(deletions_insertions),"Duplications"=nrow(duplications), "Translocations"=nrow(translocations))
  return(allmantatypes)
  typestable<-as.data.frame(table(uniquesomaticmetadatadf$MutationType))
  genefusioncount<-subset(typestable, Var1=="gene_fusion")$Freq
  #return(genefusioncount)
  duplicationcount<-subset(typestable, Var1=="duplication")$Freq
  #return(duplicationcount)
  gf_fvcount<-subset(typestable, Var1=="gene_fusion&frameshift_variant")$Freq
  #return(gf_fvcount)
  strengthstable<-as.data.frame(table(uniquesomaticmetadatadf$MutationStrength))
  highprop<-subset(strengthstable, Var1=="HIGH")$Freq/sum(strengthstable$Freq)
  #return(highprop)
}
scatterplot_func<-function(somaticmetadatadf) {
  #########comparison for all somatic mutations:##########
  somaticmetadatadf<-somdf

  justhet<-subset(somaticmetadatadf, genotype !="A/A" & genotype !="T/T" & genotype !="C/C" & genotype !="G/G") 
  orderedjusthet<- justhet[order(justhet$sample),]
  justhet$sample<-as.factor(justhet$sample)
  #sample_levels<-c(levels(justhet$sample))
  data<-as.data.frame(table(orderedjusthet$sample))
  orderedjusthet$generic<-c(rep("sample1",data$Freq[1]),rep("sample2",data$Freq[2]), rep("sample3",data$Freq[3]), rep("sample4",data$Freq[4]))
  tx <- split(orderedjusthet, orderedjusthet$generic)
  #one<-noquote(paste("`",sample_levels[1],"`", sep=""))
  sample1<-tx$sample1#`AHE01-1_S1` #sample_levels[1]
  sample2<-tx$sample2
  sample3<-tx$sample3
  sample4<-tx$sample4
  props1<- c(sample1$altdepth/sample1$totaldepth.x,  sample3$altdepth/sample3$totaldepth.x)
  props2<- c(sample2$altdepth/sample2$totaldepth.x, sample4$altdepth/sample4$totaldepth.x)
  VAFavg<- (props1+props2)/2
  VAFavg<-as.data.frame(VAFavg)
  ggplot(data=VAFavg, aes(x=VAFavg, fill="red"))+
    geom_density(alpha=0.5)+
    theme_bw()
}





#GOH vs LOH for multisamples
ahb125<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB125-129_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB125denom
ahb145<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB145-149_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB145denom
ahb151<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB151-155_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB151denom
ahb176<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB176-180_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB176denom
ahb195<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB195-199_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB195denom
ahb70<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB70-74_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB70denom
ahb90<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB90-94_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHB90denom
bikrates_GOH<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS46<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS46-50_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHAS46denom
AHAS51<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS51-55_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHAS51denom
AHAS56<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS56-60_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHAS56denom
AHAS61<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS61-65_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHAS61denom

ofurates_GOH<-c(AHAS46, AHAS51, AHAS56, AHAS61)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP01-05_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHP01denom
AHP06<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP06-10_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHP06denom
AHP16<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP16-20_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHP16denom
#AHP21<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP21-25_mT_snps.txt", full.names=T, recursive=FALSE))$GOH/AHP21denom

palrates_GOH<-c(AHP01,AHP06,AHP16)#,AHP21)

ahb125<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB125-129_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB125denom
ahb145<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB145-149_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB145denom
ahb151<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB151-155_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB151denom
ahb176<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB176-180_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB176denom
ahb195<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB195-199_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB195denom
ahb70<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB70-74_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB70denom
ahb90<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHB90-94_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHB90denom
bikrates_LOH<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS46<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS46-50_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHAS46denom
AHAS51<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS51-55_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHAS51denom
AHAS56<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS56-60_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHAS56denom
AHAS61<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHAS61-65_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHAS61denom

ofurates_LOH<-c(AHAS46, AHAS51, AHAS56, AHAS61)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP01-05_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHP01denom
AHP06<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP06-10_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHP06denom
AHP16<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP16-20_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHP16denom
#AHP21<-soma_func(list.files(path="~/Documents/GitHub/Bikini/mappedtoahya/forR", pattern="AHP21-25_mT_snps.txt", full.names=T, recursive=FALSE))$LOH/AHP21denom

palrates_LOH<-c(AHP01,AHP06,AHP16)#,AHP21)

gohrates<-c(bikrates_GOH, ofurates_GOH, palrates_GOH)
lohrates<-c(bikrates_LOH, ofurates_LOH, palrates_LOH)
percentloh<-lohrates/(gohrates+lohrates)
percentgoh<-gohrates/(gohrates+lohrates)

size<-c(bikinisize,ofusizes,palausizes[1:3])
location<-c(rep("Bikini",7),rep("Ofu",4),rep("Palau",3))
treatment<-c(rep("Irradiated",7),rep("Not Irradiated",7))
colonynames<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16") #AHP21
ratioavg<-c(mean(bikrates_GOH/bikrates_LOH))
#AvgDepth<-c(AHB70avg, AHB90avg, AHB125avg, AHB145avg, AHB151avg, AHB176avg, AHB195avg,
            AHAS46avg, AHAS51avg, AHAS56avg, AHAS61avg, AHP01avg, AHP06avg, AHP16avg)#, AHP21avg)
AvgDepth<-c(AHB70avgsevensamples, AHB90avgsevensamples, AHB125avgsevensamples, AHB145avgsevensamples, AHB151avgsevensamples, AHB176avgsevensamples, AHB195avgsevensamples,
            AHAS46avgsevensamples, AHAS51avgsevensamples, AHAS56avgsevensamples, AHAS61avgsevensamples, AHP01avgsevensamples, AHP06avgsevensamples, AHP16avgsevensamples)
typeframe<-data.frame(size,location, colonynames, lohrates, gohrates, treatment, AvgDepth, percentloh, percentgoh)
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
gohbox | lohbox
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
  geom_point() +#geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(gohrates)*1.3, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("# GOH SNVs/bp")+ xlab("Average Depth")+
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
lohvdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=lohrates, group=location, colour=location)))+
  geom_point() +#geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(lohrates)*1.3, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("#LOH SNVs/bp")+ xlab("Average Depth")+
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
ratiovdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=ratio, group=location, colour=location)))+
  geom_point() +geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(AvgDepth)*1.1,label.y=max(ratio)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
percentlohvdepth<-ggplot(typeframe, (aes(x=AvgDepth, y=(percentloh)*100, group=location, colour=location)))+
  geom_point() +#geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(percentloh*100)*0.5, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("% of SNVs that are LOH")+ xlab("Average Depth")+
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
#this is fig:
gohvdepth | lohvdepth | percentlohvdepth
Figure6<-gohvdepth + theme(legend.position = "none")| lohvdepth + theme(legend.position = "none") | percentlohvdepth +scale_color_discrete(name="Location",labels=c("Bikini","Ofu","Palau")) + 
     theme(legend.text=element_text(size=15))
Figure6 + plot_annotation(tag_levels = 'a')
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
#rates v location:
percentlohbox.resid<-ggplot(typeframe, aes(x=location, y=loh.resid, group=location, colour=treatment))+
  geom_boxplot()+
  theme_bw()
Anova(aov(I(loh.resid*1e6)~location, typeframe), type="III")
Anova(aov(I(lohrates*1e6)~location+AvgDepth, typeframe), type="III")

#residuals v size
lohpoint.resid<-ggplot(typeframe, (aes(x=size, y=loh.resid, group=location, colour=location)))+
  geom_point() +#geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  ylab("LOH SNVs ~ Depth Residuals") + xlab("Colony Diameter (cm)")+
  stat_cor(label.x=max(size)*0.5,label.y=max(fit_loh$residuals)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1))+
  theme_bw()
gohpoint.resid<-ggplot(typeframe, (aes(x=size, y=goh.resid, group=location, colour=location)))+
  geom_point() +#geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  stat_cor(label.x=max(size)*0.5,label.y=max(fit_goh$residuals)*1.3, aes(group=1)) +
  ylab("GOH SNVs ~ Depth Residuals") + xlab("Colony Diameter (cm)")+
  geom_smooth(method=lm, aes(group=1))+
  theme_bw()
lohpercent.resid<-ggplot(typeframe, (aes(x=size, y=lohpercent.resid, group=location, colour=location)))+
  geom_point() +#geom_text(aes(label=colonynames),hjust=0, vjust=0)+
  ylab(" % LOH ~ Depth Residuals") + xlab("Colony Diameter (cm)")+
  stat_cor(label.x=max(size)*0.5,label.y=max(fit_percentloh$residuals)*1.3, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1))+
  theme_bw()
#this is fig:
gohpoint.resid | lohpoint.resid | lohpercent.resid
#resids boxplots
percentlohbox.resid<-ggplot(typeframe, aes(x=location, y=lohpercent.resid, group=location, colour=treatment))+
  geom_boxplot()+
  ylab(" % LOH ~ Depth residuals")+
  theme_bw()
lohbox.resid<-ggplot(typeframe, aes(x=location, y=loh.resid, group=location, colour=treatment))+
  geom_boxplot()+
  ylab("LOH SNVs ~ Depth residuals")+
  theme_bw()
gohbox.resid<-ggplot(typeframe, aes(x=location, y=goh.resid, group=location, colour=treatment))+
  geom_boxplot()+
  ylab("GOH SNVs ~ Depth residuals")+
  theme_bw()
gohbox.resid | lohbox.resid | percentlohbox.resid
#rates boxplots
percentlohbox<-ggplot(typeframe, aes(x=location, y=percentloh*100, group=location, colour=treatment))+
  geom_boxplot()+
  ylab(" % of SNVs that are LOH")+
  geom_point(aes(color=treatment, group=location), position = position_dodge(width = 0.75))+
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
percentlohbox_treatment<-ggplot(typeframe, aes(x=treatment, y=percentloh*100, group=treatment, colour=treatment))+
  geom_boxplot()+
  ylab(" % of SNVs that are LOH")+
  geom_point(aes(color=treatment, group=treatment), position = position_dodge(width = 0.75))+
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
##this is fig:
Figure7<-percentlohbox | percentlohbox_treatment
Figure7 + plot_annotation(tag_levels = 'a')
leveneTest(percentloh~location, typeframe)
leveneTest(percentloh~treatment, typeframe)
summary.lm(aov(percentloh~location+AvgDepth, typeframe))
summary.lm(aov(percentloh~treatment+AvgDepth, typeframe))

##ANCOVA for goh/loh
library(car)
leveneTest(gohrates~location, typeframe)
fit_goh <- aov(gohrates~location, typeframe)
fit_goh2 <- aov(gohrates~location+AvgDepth, typeframe)
summary.lm(fit_goh2)
anova(fit_goh, fit_goh2)
fit_goh3<-aov(I(gohrates*1e6)~location+AvgDepth,typeframe)
Anova(fit_goh3, type="III")

leveneTest(lohrates~location, typeframe)
fit_loh <- aov(lohrates~location, typeframe)
fit_loh2 <- aov(lohrates~location+AvgDepth, typeframe)
summary.lm(fit_loh2)
anova(fit_loh, fit_loh2)
fit_loh3<-aov(I(lohrates*1e6)~location+AvgDepth,typeframe)
Anova(fit_loh3, type="III")

leveneTest(percentloh~location, typeframe)
fit_percentloh <- aov(percentloh~location, typeframe)
fit_percentloh2 <- aov(percentloh~location+AvgDepth, typeframe)
fit_percentloh3<-aov(I(percentloh)~location+AvgDepth,typeframe)
Anova(fit_percentloh3, type="III")
summary.lm(fit_percentloh2)
anova(fit_percentloh, fit_percentloh2)

#treatment ancova
leveneTest(gohrates~treatment, typeframe)
fit_goh <- aov(gohrates~treatment, typeframe)
fit_goh2 <- aov(gohrates~treatment+AvgDepth, typeframe)
summary.lm(fit_goh2)
anova(fit_goh, fit_goh2)
fit_goh3<-aov(I(gohrates*1e6)~treatment+AvgDepth,typeframe)
Anova(fit_goh3, type="III")

leveneTest(lohrates~treatment, typeframe)
fit_loh <- aov(lohrates~treatment, typeframe)
fit_loh2 <- aov(lohrates~treatment+AvgDepth, typeframe)
summary.lm(fit_loh2)
anova(fit_loh, fit_loh2)
fit_loh3<-aov(I(lohrates*1e6)~treatment+AvgDepth,typeframe)
Anova(fit_loh3, type="III")

leveneTest(percentloh~treatment, typeframe)
fit_percentloh <- aov(percentloh~treatment, typeframe)
fit_percentloh2 <- aov(percentloh~treatment+AvgDepth, typeframe)
fit_percentloh3<-aov(I(percentloh)~treatment+AvgDepth,typeframe)
Anova(fit_percentloh3, type="III")
summary.lm(fit_percentloh2)
anova(fit_percentloh, fit_percentloh2)

##ANCOVA with size covariate
leveneTest(gohrates~location, typeframe)
fit_goh <- aov(gohrates~location, typeframe)
fit_goh2 <- aov(gohrates~location+size, typeframe)
summary.lm(fit_goh2)
anova(fit_goh, fit_goh2)
fit_goh3<-aov(I(gohrates*1e6)~location+size,typeframe)
Anova(fit_goh3, type="III")

leveneTest(lohrates~location, typeframe)
fit_loh <- aov(lohrates~location, typeframe)
fit_loh2 <- aov(lohrates~location+size, typeframe)
summary.lm(fit_loh2)
anova(fit_loh, fit_loh2)
fit_loh3<-aov(I(lohrates*1e6)~location+size,typeframe)
Anova(fit_loh3, type="III")

leveneTest(percentloh~location, typeframe)
fit_percentloh <- aov(percentloh~location, typeframe)
fit_percentloh2 <- aov(percentloh~location+size, typeframe)
fit_percentloh3<-aov(I(percentloh)~location+size,typeframe)
Anova(fit_percentloh3, type="III")
summary.lm(fit_percentloh2)
anova(fit_percentloh, fit_percentloh2)

################
##FOURSAMPLE indivs:
#SNPS
snps_AHE01<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE01-02_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE01FOURSAMPLEdenom
#AHE05<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE05-06_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE05FOURSAMPLEdenom
snps_AHE07<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE07-08_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE07FOURSAMPLEdenom
snps_AHE09<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE09-10_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE09FOURSAMPLEdenom
snps_AHE11<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE11-12_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE11FOURSAMPLEdenom

snps_enerates<-c(snps_AHE01,snps_AHE07,snps_AHE09,snps_AHE11)

snps_ahb125<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB125-129_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB125FOURSAMPLEdenom
snps_ahb145<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB145-149_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB145FOURSAMPLEdenom
snps_ahb151<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB151-155_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB151FOURSAMPLEdenom
snps_ahb176<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB176-180_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB176FOURSAMPLEdenom
snps_ahb195<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB195-199_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB195FOURSAMPLEdenom
snps_ahb70<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB70-74_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB70FOURSAMPLEdenom
snps_ahb90<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB90-94_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB90FOURSAMPLEdenom
snps_bikrates<-c(snps_ahb70, snps_ahb90, snps_ahb125, snps_ahb145, snps_ahb151, snps_ahb176, snps_ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
snps_AHAS41<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS41-45_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS41FOURSAMPLEdenom
snps_AHAS46<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS46-50_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS46FOURSAMPLEdenom
snps_AHAS51<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS51-55_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS51FOURSAMPLEdenom
snps_AHAS56<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS56-60_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS56FOURSAMPLEdenom
snps_AHAS61<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS61-65_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS61FOURSAMPLEdenom

snps_ofurates<-c(snps_AHAS41,snps_AHAS46, snps_AHAS51, snps_AHAS56, snps_AHAS61)
ofusizes_with41<-c(206,ofusizes)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

snps_AHP01<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP01-05_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP01FOURSAMPLEdenom
snps_AHP06<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP06-10_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP06FOURSAMPLEdenom
snps_AHP16<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP16-20_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP16FOURSAMPLEdenom
snps_AHP21<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP21-25_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP21FOURSAMPLEdenom

snps_palrates<-c(snps_AHP01,snps_AHP06,snps_AHP16)

##################
##pop_pairwise snps
#files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya", pattern="*mf_snpsPOP.txt_melted.txt", full.names=T, recursive=FALSE)
files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya", pattern="*mf_FOURsnps.txt_melted.txt", full.names=T, recursive=FALSE)
#svfiles<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya", pattern="*svPOP.txt_melted.txt", full.names=T, recursive=FALSE)
svfiles<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya", pattern="*sv.txt_melted.txt", full.names=T, recursive=FALSE)
svfiles<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya", pattern="*svPOP.txt_melted.txt", full.names=T, recursive=FALSE)

for (i in 1:length(files)) {
  file =files[i]
  base<-basename(file)
  colonies<-strsplit(base, "\\_")[[1]][1]
  colony1<-strsplit(colonies, "\\-")[[1]][1]
  colony2<-strsplit(colonies, "\\-")[[1]][2]
  
  mf<-strsplit(base, "\\_")[[1]][2]
  type<-strsplit(base, "\\_")[[1]][3]
  typ1<-strsplit(type, "\\.")[[1]][1]
  
  #print(type)
  pasted<-paste(colony1,"_",colony2, "_", typ1, sep="")
  assign(pasted,AHE_soma_func(file))
  #AHE_soma_func(file)
}
for (i in 1:length(svfiles)) {
  file =svfiles[i]
  base<-basename(file)
  colonies<-strsplit(base, "\\_")[[1]][1]
  colony1<-strsplit(colonies, "\\-")[[1]][1]
  colony2<-strsplit(colonies, "\\-")[[1]][2]
  
  mf<-strsplit(base, "\\_")[[1]][2]
  type<-strsplit(base, "\\_")[[1]][3]
  #typ1<-strsplit(type, "\\.")[[1]][1]
  
  #print(type)
  pasted_dup<-paste(colony1,"_",colony2, "_", "dup", sep="")
  assign(pasted_dup,SV_soma_func(file)$Duplications)
  
  pasted_transloc<-paste(colony1,"_",colony2, "_", "transloc", sep="")
  assign(pasted_transloc,SV_soma_func(file)$Translocations)
  
  pasted_mantaindel<-paste(colony1,"_",colony2, "_", "mantaindel", sep="")
  assign(pasted_mantaindel,SV_soma_func(file)$Indels)
  #AHE_soma_func(file)
}

AHAS51_61_snpsPOP<-1148477
AHAS41_46_snpsPOP<-1123973
AHAS41_51_snpsPOP<-1220052
AHAS41_56_snpsPOP<-1174711
AHAS41_61_snpsPOP<-1208701
AHAS46_51_snpsPOP<-759513
AHAS46_56_snpsPOP<-752767
AHAS46_61_snpsPOP<-807759
AHAS51_56_snpsPOP<-1121200
AHAS56_61_snpsPOP<-11031
  
AHB125_145_snpsPOP<-1085211
  AHB125_151_snpsPOP<-1638430
  AHB125_176_snpsPOP<-1362629
  AHB125_195_snpsPOP<-1376582
  AHB125_70_snpsPOP<-1346644
  AHB125_90_snpsPOP<-1306949
  AHB145_151_snpsPOP<-1419792
  AHB145_176_snpsPOP<-1181821
  AHB145_195_snpsPOP<-1195813
AHB145_70_snpsPOP<-1178082
AHB145_90_snpsPOP<-1158624
AHB151_176_snpsPOP<-1465459
AHB151_195_snpsPOP<-1492570
AHB151_70_snpsPOP<-1457266
AHB151_90_snpsPOP<-1434647
AHB176_195_snpsPOP<-1333387
AHB176_70_snpsPOP<-1316516
AHB176_90_snpsPOP<-1284876
AHB195_70_snpsPOP<-535673
AHB195_90_snpsPOP<-500924
AHB70_90_snpsPOP<-1310171

AHP01_06_snpsPOP<-1207671
AHP01_16_snpsPOP<-1584864
AHP01_21_snpsPOP<-1483824
AHP06_16_snpsPOP<-1388612
AHP06_21_snpsPOP<-1008007
AHP16_21_snpsPOP<-1652136

AHE01_07_snpsPOP<-28307
AHE01_09_snpsPOP<-12221
AHE01_11_snpsPOP<-252081
AHE07_09_snpsPOP<-13825
AHE07_11_snpsPOP<-192268
AHE09_11_snpsPOP<-154586

AHB_snps_pairwise<-c(AHB125_145_snpsPOP, AHB125_151_snpsPOP, AHB125_176_snpsPOP, AHB125_195_snpsPOP, AHB125_70_snpsPOP, AHB125_90_snpsPOP,
                     AHB145_151_snpsPOP, AHB145_176_snpsPOP, AHB145_195_snpsPOP, AHB145_70_snpsPOP, AHB145_90_snpsPOP,
                     AHB151_176_snpsPOP, AHB151_195_snpsPOP, AHB151_70_snpsPOP, AHB151_90_snpsPOP,
                     AHB176_195_snpsPOP, AHB176_70_snpsPOP, AHB176_90_snpsPOP, 
                     AHB195_70_snpsPOP, AHB195_90_snpsPOP,
                     AHB70_90_snpsPOP)
AHAS_snps_pairwise<-c(AHAS41_46_snpsPOP, AHAS41_51_snpsPOP, AHAS41_56_snpsPOP, AHAS41_61_snpsPOP,
                      AHAS41_51_snpsPOP, AHAS46_56_snpsPOP, AHAS46_61_snpsPOP, 
                      AHAS51_56_snpsPOP, AHAS51_61_snpsPOP,
                      AHAS56_61_snpsPOP)
AHE_snps_pairwise<-c(AHE01_07_snpsPOP, AHE01_09_snpsPOP, AHE01_11_snpsPOP,
                     AHE07_09_snpsPOP, AHE07_11_snpsPOP,
                     AHE09_11_snpsPOP)
AHP_snps_pairwise<-c(AHP01_06_snpsPOP, AHP01_16_snpsPOP,
                     AHP06_16_snpsPOP
                     )
AHB_transloc_pairwise<-c(AHB125_145_transloc, AHB125_151_transloc, AHB125_176_transloc, AHB125_195_transloc, AHB125_70_transloc, AHB125_90_transloc,
                     AHB145_151_transloc, AHB145_176_transloc, AHB145_195_transloc, AHB145_70_transloc, AHB145_90_transloc,
                     AHB151_176_transloc, AHB151_195_transloc, AHB151_70_transloc, AHB151_90_transloc,
                     AHB176_195_transloc, AHB176_70_transloc, AHB176_90_transloc, 
                     AHB195_70_transloc, AHB195_90_transloc)#,
                     #AHB70_90_transloc)
AHAS_transloc_pairwise<-c(AHAS41_46_transloc, AHAS41_51_transloc, AHAS41_56_transloc, AHAS41_61_transloc,
                      AHAS41_51_transloc, AHAS46_56_transloc, AHAS46_61_transloc, 
                      AHAS51_56_transloc, AHAS51_61_transloc, AHAS56_61_transloc)
AHE_transloc_pairwise<-c(AHE01_07_transloc, AHE01_09_transloc, AHE01_11_transloc,
                     AHE07_09_transloc, AHE07_11_transloc,
                     AHE09_11_transloc)
AHP_transloc_pairwise<-c(AHP01_06_transloc, AHP01_16_transloc,
                     AHP06_16_transloc
                     )
#
AHB_dup_pairwise<-c(AHB125_145_dup, AHB125_151_dup, AHB125_176_dup, AHB125_195_dup, AHB125_70_dup, AHB125_90_dup,
                         AHB145_151_dup, AHB145_176_dup, AHB145_195_dup, AHB145_70_dup, AHB145_90_dup,
                         AHB151_176_dup, AHB151_195_dup, AHB151_70_dup, AHB151_90_dup,
                         AHB176_195_dup, AHB176_70_dup, AHB176_90_dup, 
                         AHB195_70_dup, AHB195_90_dup)#,
#AHB70_90_dup)
AHAS_dup_pairwise<-c(AHAS41_46_dup, AHAS41_51_dup, AHAS41_56_dup, AHAS41_61_dup,
                          AHAS41_51_dup, AHAS46_56_dup, AHAS46_61_dup, 
                          AHAS51_56_dup, AHAS51_61_dup, AHAS56_61_dup)
AHE_dup_pairwise<-c(AHE01_07_dup, AHE01_09_dup, AHE01_11_dup,
                         AHE07_09_dup, AHE07_11_dup,
                         AHE09_11_dup)
AHP_dup_pairwise<-c(AHP01_06_dup, AHP01_16_dup,
                         AHP06_16_dup
                         )
#
AHB_mantaindel_pairwise<-c(AHB125_145_mantaindel, AHB125_151_mantaindel, AHB125_176_mantaindel, AHB125_195_mantaindel, AHB125_70_mantaindel, AHB125_90_mantaindel,
                         AHB145_151_mantaindel, AHB145_176_mantaindel, AHB145_195_mantaindel, AHB145_70_mantaindel, AHB145_90_mantaindel,
                         AHB151_176_mantaindel, AHB151_195_mantaindel, AHB151_70_mantaindel, AHB151_90_mantaindel,
                         AHB176_195_mantaindel, AHB176_70_mantaindel, AHB176_90_mantaindel, 
                         AHB195_70_mantaindel, AHB195_90_mantaindel)#,
#AHB70_90_mantaindel)
AHAS_mantaindel_pairwise<-c(AHAS41_46_mantaindel, AHAS41_51_mantaindel, AHAS41_56_mantaindel, AHAS41_61_mantaindel,
                          AHAS41_51_mantaindel, AHAS46_56_mantaindel, AHAS46_61_mantaindel, 
                          AHAS51_56_mantaindel, AHAS51_61_mantaindel, AHAS56_61_mantaindel)
AHE_mantaindel_pairwise<-c(AHE01_07_mantaindel, AHE01_09_mantaindel, AHE01_11_mantaindel,
                         AHE07_09_mantaindel, AHE07_11_mantaindel,
                         AHE09_11_mantaindel)
AHP_mantaindel_pairwise<-c(AHP01_06_mantaindel, AHP01_16_mantaindel, 
                         AHP06_16_mantaindel
                         )
#
AHB_POP_denoms<-c(AHB125_145denom, AHB125_151denom, AHB125_176denom, AHB125_195denom, AHB125_70denom, AHB125_90denom,
                  AHB145_151denom, AHB145_176denom, AHB145_195denom, AHB145_70denom, AHB145_90denom,
                  AHB151_176denom, AHB151_195denom, AHB151_70denom, AHB151_90denom,
                  AHB176_195denom, AHB176_70denom, AHB176_90denom, 
                  AHB195_70denom, AHB195_90denom,
                  AHB70_90denom)
AHAS_POP_denoms<-c(AHAS41_46denom, AHAS41_51denom, AHAS41_56denom, AHAS41_61denom,
                   AHAS41_51denom, AHAS46_56denom, AHAS46_61denom, 
                   AHAS51_56denom, AHAS51_61denom,
                   AHAS56_61denom)
AHE_POP_denoms<-c(AHE01_07denom, AHE01_09denom, AHE01_11denom,
                  AHE07_09denom, AHE07_09denom,
                  AHE09_11denom)
AHP_POP_denoms<-c(AHP01_06denom, AHP01_16denom,
                  AHP06_16denom
                  )
pop_depths<-c(AHB125_145_depth, AHB125_151_depth, AHB125_176_depth, AHB125_195_depth, AHB125_70_depth, AHB125_90_depth,
                  AHB145_151_depth, AHB145_176_depth, AHB145_195_depth, AHB145_70_depth, AHB145_90_depth,
                  AHB151_176_depth, AHB151_195_depth, AHB151_70_depth, AHB151_90_depth,
                  AHB176_195_depth, AHB176_70_depth, AHB176_90_depth, 
                  AHB195_70_depth, AHB195_90_depth,
                  AHAS41_46_depth, AHAS41_51_depth, AHAS41_56_depth, AHAS41_61_depth,
                   AHAS41_51_depth, AHAS46_56_depth, AHAS46_61_depth, 
                   AHAS51_56_depth, AHAS51_61_depth,
                   AHAS56_61_depth,AHE01_07_depth, AHE01_09_depth, AHE01_11_depth,
                  AHE07_09_depth, AHE07_09_depth,
                  AHE09_11_depth,AHP01_06_depth, AHP01_16_depth,
                  AHP06_16_depth)

AHB_snps_pairwise_rate<-AHB_snps_pairwise/AHB_POP_denoms
#AHB_indels_pairwise_rate<-AHB_indels_pairwise/AHB_POP_denoms
AHB_transloc_pairwise_rate<-AHB_transloc_pairwise/AHB_POP_denoms[1:20] #need AHB70-90
AHAS_snps_pairwise_rate<-AHAS_snps_pairwise/AHAS_POP_denoms
#AHAS_indels_pairwise_rate<-AHAS_indels_pairwise/AHAS_POP_denoms
AHAS_transloc_pairwise_rate<-AHAS_transloc_pairwise/AHAS_POP_denoms
AHE_snps_pairwise_rate<-AHE_snps_pairwise/AHE_POP_denoms
#AHE_indels_pairwise_rate<-AHE_indels_pairwise/AHE_POP_denoms
AHE_transloc_pairwise_rate<-AHE_transloc_pairwise/AHE_POP_denoms
AHP_snps_pairwise_rate<-AHP_snps_pairwise/AHP_POP_denoms
#AHP_indels_pairwise_rate<-AHP_indels_pairwise/AHP_POP_denoms
AHP_transloc_pairwise_rate<-AHP_transloc_pairwise/AHP_POP_denoms

AHB_dup_pairwise_rate<-AHB_dup_pairwise/AHB_POP_denoms[1:20]
AHAS_dup_pairwise_rate<-AHAS_dup_pairwise/AHAS_POP_denoms
AHE_dup_pairwise_rate<-AHE_dup_pairwise/AHE_POP_denoms
AHP_dup_pairwise_rate<-AHP_dup_pairwise/AHP_POP_denoms
AHB_mantaindel_pairwise_rate<-AHB_mantaindel_pairwise/AHB_POP_denoms[1:20]
AHAS_mantaindel_pairwise_rate<-AHAS_mantaindel_pairwise/AHAS_POP_denoms
AHE_mantaindel_pairwise_rate<-AHE_mantaindel_pairwise/AHE_POP_denoms
AHP_mantaindel_pairwise_rate<-AHP_mantaindel_pairwise/AHP_POP_denoms
##for indivs
bikdenoms<-c(AHB125FOURSAMPLEdenom, AHB145FOURSAMPLEdenom, AHB151FOURSAMPLEdenom, AHB176FOURSAMPLEdenom, AHB195FOURSAMPLEdenom, AHB70FOURSAMPLEdenom, AHB90FOURSAMPLEdenom)
mantaindel_bikrates<-c(AHB125_129_mantaindel, AHB145_149_mantaindel, AHB151_155_mantaindel, AHB176_180_mantaindel, AHB195_199_mantaindel, AHB70_74_mantaindel, AHB90_94_mantaindel)/bikdenoms
dup_bikrates<-c(AHB125_129_dup, AHB145_149_dup, AHB151_155_dup, AHB176_180_dup, AHB195_199_dup, AHB70_74_dup, AHB90_94_dup)/bikdenoms
transloc_bikrates<-c(AHB125_129_transloc, AHB145_149_transloc, AHB151_155_transloc, AHB176_180_transloc, AHB195_199_transloc, AHB70_74_transloc, AHB90_94_transloc)/bikdenoms
snps_bikrates<-c(AHB125_129_FOURsnps, AHB145_149_FOURsnps, AHB151_155_FOURsnps, AHB176_180_FOURsnps, AHB195_199_FOURsnps, AHB70_74_FOURsnps, AHB90_94_FOURsnps)/bikdenoms
#
paldenoms<-c(AHP01FOURSAMPLEdenom, AHP06FOURSAMPLEdenom, AHP16FOURSAMPLEdenom)#, AHB195FOURSAMPLEdenom, AHB70FOURSAMPLEdenom, AHB90FOURSAMPLEdenom)
mantaindel_palrates<-c(AHP01_05_mantaindel, AHP06_10_mantaindel, AHP16_20_mantaindel)/paldenoms
transloc_palrates<-c(AHP01_05_transloc, AHP06_10_transloc, AHP16_20_transloc)/paldenoms
dup_palrates<-c(AHP01_05_dup, AHP06_10_dup, AHP16_20_dup)/paldenoms
snps_palrates<-c(AHP01_05_FOURsnps, AHP06_10_FOURsnps, AHP16_20_FOURsnps)/paldenoms

#
ofudenoms<-c(AHAS41FOURSAMPLEdenom, AHAS46FOURSAMPLEdenom, AHAS51FOURSAMPLEdenom, AHAS56FOURSAMPLEdenom, AHAS61FOURSAMPLEdenom)
mantaindel_ofurates<-c(AHAS41_45_mantaindel, AHAS46_50_mantaindel, AHAS51_55_mantaindel, AHAS56_60_mantaindel, AHAS61_65_mantaindel)/ofudenoms  
dup_ofurates<-c(AHAS41_45_dup, AHAS46_50_dup, AHAS51_55_dup, AHAS56_60_dup, AHAS61_65_dup)/ofudenoms    
transloc_ofurates<-c(AHAS41_45_transloc, AHAS46_50_transloc, AHAS51_55_transloc, AHAS56_60_transloc, AHAS61_65_transloc)/ofudenoms  
snps_ofurates<-c(AHAS41_45_FOURsnps, AHAS46_50_FOURsnps, AHAS51_55_FOURsnps, AHAS56_60_FOURsnps, AHAS61_65_FOURsnps)/ofudenoms    

#
enedenoms<-c(AHE01FOURSAMPLEdenom, AHE07FOURSAMPLEdenom, AHE09FOURSAMPLEdenom, AHE11FOURSAMPLEdenom)
mantaindel_enerates<-c(AHE01_02_mantaindel, AHE07_08_mantaindel, AHE09_10_mantaindel, AHE11_12_mantaindel)/enedenoms
transloc_enerates<-c(AHE01_02_transloc, AHE07_08_transloc, AHE09_10_transloc, AHE11_12_transloc)/enedenoms
dup_enerates<-c(AHE01_02_dup, AHE07_08_dup, AHE09_10_dup, AHE11_12_dup)/enedenoms
snps_enerates<-c(AHE01_02_FOURsnps, AHE07_08_FOURsnps, AHE09_10_FOURsnps, AHE11_12_FOURsnps)/enedenoms

dup_rates<-c(dup_bikrates, dup_enerates, dup_ofurates,dup_palrates, 
             AHB_dup_pairwise_rate, AHE_dup_pairwise_rate, AHAS_dup_pairwise_rate, AHP_dup_pairwise_rate)
transloc_rates<-c(transloc_bikrates, transloc_enerates, transloc_ofurates,transloc_palrates, 
                  AHB_transloc_pairwise_rate, AHE_transloc_pairwise_rate, AHAS_transloc_pairwise_rate, AHP_transloc_pairwise_rate)
mantaindel_rates<-c(mantaindel_bikrates, mantaindel_enerates, mantaindel_ofurates,mantaindel_palrates, 
                    AHB_mantaindel_pairwise_rate, AHE_mantaindel_pairwise_rate, AHAS_mantaindel_pairwise_rate, AHP_mantaindel_pairwise_rate)
#indels_rates<-c(indels_bikrates, indels_enerates, indels_ofurates,indels_palrates, 
    #            AHB_indels_pairwise_rate, AHE_indels_pairwise_rate, AHAS_indels_pairwise_rate, AHP_indels_pairwise_rate)
snps_rates<-c(snps_bikrates, snps_enerates, snps_ofurates,snps_palrates, 
              AHB_snps_pairwise_rate[1:20], AHE_snps_pairwise_rate, AHAS_snps_pairwise_rate, AHP_snps_pairwise_rate)
#depthsfull<-rep(depths,2)
locationsfull<-c(rep("Bikini",7), rep("Enewetak",4), rep("Ofu",5), rep("Palau",3), rep("Bikini",20), rep("Enewetak",6), rep("Ofu",10), rep("Palau",3))
indorpop<-c(rep("ind", 19), rep("pop", 39))
treatment<-c(rep("irradiated",11), rep("not irradiated", 8),rep("irradiated",26), rep("not irradiated", 13))
avgdepth<-c()
names<-c("AHB70",  "AHB90",  "AHB125", "AHB145", "AHB151" ,"AHB176" ,"AHB195",
         "AHE01","AHE07","AHE09","AHE11",
         "AHAS41", "AHAS46", "AHAS51" ,"AHAS56", "AHAS61" ,
         "AHP01"  ,"AHP06", "AHP16", 
         "AHB125_145", "AHB125_151", "AHB125_176", "AHB125_195", "AHB125_70", "AHB125_90",
         "AHB145_151", "AHB145_176", "AHB145_195", "AHB145_70", "AHB145_90",
         "AHB151_176", "AHB151_195", "AHB151_70", "AHB151_90",
         "AHB176_195", "AHB176_70", "AHB176_90", 
         "AHB195_70", "AHB195_90",
         "AHAS41_46", "AHAS41_51", "AHAS41_56", "AHAS41_61",
         "AHAS41_51", "AHAS46_56", "AHAS46_61", 
         "AHAS51_56", "AHAS51_61",
         "AHAS56_61","AHE01_07", "AHE01_09", "AHE01_11",
         "AHE07_09", "AHE07_09",
         "AHE09_11","AHP01_06", "AHP01_16",
         "AHP06_16")
df<-data.frame(names, snps_rates, dup_rates, mantaindel_rates, transloc_rates, locationsfull, indorpop, treatment)# indels_rates, 
df$AvgDepth<-c(AHB70avg, AHB90avg, AHB125avg, AHB145avg, AHB151avg, AHB176avg, AHB195avg,
                    AHE01avg, AHE07avg, AHE09avg, AHE11avg, 
                    AHAS41avg, AHAS46avg, AHAS51avg, AHAS56avg, AHAS61avg,
                    AHP01avg, AHP06avg, AHP16avg, pop_depths
               )
df_boxun_loc_dup<-ggplot(df, aes(x=locationsfull, y=dup_rates,label=names))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("# duplicates/bp")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151"))+# | 
                          #justind$names == "AHAS61" | justind$names == "AHE11"))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_loc_transloc<-ggplot(df, aes(x=locationsfull, y=transloc_rates,label=names))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("# translocations/bp")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151"))+# | 
                          #justind$names == "AHAS61" | justind$names == "AHE11"))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_loc_mantaindel<-ggplot(df, aes(x=locationsfull, y=mantaindel_rates,label=names))+
  geom_boxplot(aes(color=indorpop))+
  ylab("# indels/bp")+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151"))+# | 
                          #justind$names == "AHAS61" | justind$names == "AHE11"))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_loc_snps<-ggplot(df, aes(x=locationsfull, y=snps_rates,label=names))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("# SNVs/bp")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151" ))+#| 
                          #justind$names == "AHAS61" | justind$names == "AHE11"))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.title.x=element_blank())
#this is figure:
Figure4<-(df_boxun_loc_snps + theme(legend.position = "none")| df_boxun_loc_dup + 
    scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + 
    theme(legend.text=element_text(size=15)))/( df_boxun_loc_mantaindel + theme(legend.position = "none")|
                                                  df_boxun_loc_transloc + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) +
                                                  theme(legend.text=element_text(size=15)))
Figure4 + plot_annotation(tag_levels = 'a')
## irr vs not irr
df_boxun_dup<-ggplot(df, aes(x=treatment, y=dup_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  ylab("# of duplications per bp")+
  theme(text = element_text(size=25),
                                         axis.text.x = element_text(size=20,angle=90, hjust=1),
                                         axis.title.x=element_blank())
df_boxun_transloc<-ggplot(df, aes(x=treatment, y=transloc_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  ylab("# of translocations per bp")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_mantaindel<-ggplot(df, aes(x=treatment, y=mantaindel_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  ylab("# of indels per bp")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_snps<-ggplot(df, aes(x=treatment, y=snps_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  ylab("# of SNVs per bp")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
#this is figure:
df_boxun_snps + theme(legend.position = "none")| df_boxun_dup + theme(legend.position = "none")| df_boxun_mantaindel + theme(legend.position = "none")| df_boxun_transloc + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + theme(legend.text=element_text(size=15))
#
#snps vs depth
justind<-df[c(1:19),]
justpop<-df[c(20:58),]
justind_minusAHE07<-justind[c(1:8,10:19),]
justind_minusAHE07_minusAHB195<-justind[c(1:4,6:8,10:19),]

#justind$AvgDepth<-c(AHB70avg, AHB90avg, AHB125avg, AHB145avg, AHB151avg, AHB176avg, AHB195avg,
#                    AHE01avg, AHE07avg, AHE09avg, AHE11avg, 
 #                   AHAS41avg, AHAS46avg, AHAS51avg, AHAS56avg, AHAS61avg,
#                    AHP01avg, AHP06avg, AHP16avg)
#snps vs depth
snpsvdepth<-ggplot(justind, (aes(x=AvgDepth, y=snps_rates, group=locationsfull, colour=locationsfull,label=names)))+
  geom_point() +
  #geom_text(data=subset(justind, justind$snps_rates>mean(justind$snps_rates)+ sd(justind$snps_rates)))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151" | 
                          justind$names == "AHAS61" | justind$names == "AHE11"))+
  stat_cor(label.x=5,label.y=max(justind$snps_rates)*1.2, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  ylab("# SNVs/bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
#SVs vs depth
dupvdepth<-ggplot(justind, (aes(x=AvgDepth, y=dup_rates, group=locationsfull, colour=locationsfull, label=names)))+
  geom_point() +
  #geom_text(data=subset(justind, justind$dup_rates>mean(justind$dup_rates)+ sd(justind$dup_rates)))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151" | 
                          justind$names == "AHAS61" | justind$names == "AHE11"))+
  stat_cor(label.x=5,label.y=max(justind$dup_rates)*0.94, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# duplications/bp")+  
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
mantaindelvdepth<-ggplot(justind, (aes(x=AvgDepth, y=mantaindel_rates, group=locationsfull, colour=locationsfull,label=names)))+
  geom_point() +
  #geom_text(data=subset(justind, justind$mantaindel_rates>mean(justind$mantaindel_rates)+ sd(justind$mantaindel_rates)))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151" | 
                          justind$names == "AHAS61" | justind$names == "AHE11"))+
  stat_cor(label.x=5,label.y=max(justind$mantaindel_rates)*0.94, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# indels/bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
translocvdepth<-ggplot(justind, (aes(x=AvgDepth, y=transloc_rates, group=locationsfull, colour=locationsfull,label=names)))+
  geom_point() +
  #geom_text(data=subset(justind, justind$transloc_rates>mean(justind$transloc_rates)+ sd(justind$transloc_rates)))+
  geom_text(data=subset(justind, justind$names == "AHE07" | justind$names == "AHB151" | 
                          justind$names == "AHAS61" | justind$names == "AHE11"))+
  stat_cor(label.x=5,label.y=max(justind$transloc_rates)*0.94, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# translocations/bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
################
#justindvdepth<-snpsvdepth + theme(legend.position = "none")| dupvdepth+ theme(legend.position = "none") | mantaindelvdepth+ theme(legend.position = "none") | translocvdepth + scale_color_discrete(name="Location") + theme(legend.text=element_text(size=15))
Figure2<-(snpsvdepth + theme(legend.position = "none")| dupvdepth + 
              scale_color_discrete(name="Location",labels=c("Bikini","Enewetak","Ofu","Palau")) + 
              theme(legend.text=element_text(size=15)))/( mantaindelvdepth + theme(legend.position = "none")|
                                                            mantaindelvdepth + scale_color_discrete(name="Location",labels=c("Bikini","Enewetak","Ofu","Palau")) +
                                                            theme(legend.text=element_text(size=15)))
Figure2 + plot_annotation(tag_levels = 'a')
################
#snps vs depth
snpsvdepth<-ggplot(justind_minusAHE07_minusAHB195, (aes(x=AvgDepth, y=snps_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justind_minusAHE07_minusAHB195$snps_rates)*1.2, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  ylab("# of SNVs per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
#SVs vs depth
dupvdepth<-ggplot(justind_minusAHE07_minusAHB195, (aes(x=AvgDepth, y=dup_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justind_minusAHE07_minusAHB195$dup_rates)*0.94, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of duplications per bp")+  
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
mantaindelvdepth<-ggplot(justind_minusAHE07_minusAHB195, (aes(x=AvgDepth, y=mantaindel_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justind_minusAHE07_minusAHB195$mantaindel_rates)*0.94, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of indels per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
translocvdepth<-ggplot(justind_minusAHE07_minusAHB195, (aes(x=AvgDepth, y=transloc_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justind_minusAHE07_minusAHB195$transloc_rates)*0.94, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of translocations per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
################

#snps vs depth
snpsvdepth<-ggplot(justpop, (aes(x=AvgDepth, y=snps_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justpop$snps_rates)*1.2, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  ylab("# of SNVs per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
#SVs vs depth
dupvdepth<-ggplot(justpop, (aes(x=AvgDepth, y=dup_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justpop$dup_rates)*0.94, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of duplications per bp")+  
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
mantaindelvdepth<-ggplot(justpop, (aes(x=AvgDepth, y=mantaindel_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justpop$mantaindel_rates)*0.94, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of indels per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
translocvdepth<-ggplot(justpop, (aes(x=AvgDepth, y=transloc_rates, group=locationsfull, colour=locationsfull)))+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justpop$transloc_rates)*0.94, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of translocations per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
#########
#dupvdepth | mantaindelvdepth | translocvdepth
####this is fig : 
snpsvdepth + theme(legend.position = "none")| dupvdepth+ theme(legend.position = "none") | mantaindelvdepth+ theme(legend.position = "none") | translocvdepth + scale_color_discrete(name="Location") + theme(legend.text=element_text(size=15))
#df_boxun_snps + theme(legend.position = "none")| df_boxun_dup + theme(legend.position = "none")| df_boxun_mantaindel + theme(legend.position = "none")| df_boxun_transloc + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + theme(legend.text=element_text(size=15))

##just pop:
snpsvdepth<-ggplot(df, (aes(x=AvgDepth, y=snps_rates, group=locationsfull, colour=locationsfull)))+
  geom_point(aes(shape=indorpop)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()+
  ylab("# of SNVs per bp")+
  xlab("Average Depth")+
  scale_y_continuous(trans='log10')+
  stat_cor(label.x=10,label.y=(log10(0.01)), aes(group=1)) +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
#SVs vs depth
dupvdepth<-ggplot(df, (aes(x=AvgDepth, y=dup_rates, group=locationsfull, colour=locationsfull)))+
  geom_point(aes(shape=indorpop)) +
  stat_cor(label.x=5,label.y=log10(max(df$dup_rates)*0.94), aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of duplications per bp")+ 
  scale_y_continuous(trans='log10')+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
mantaindelvdepth<-ggplot(df, (aes(x=AvgDepth, y=mantaindel_rates, group=locationsfull, colour=locationsfull)))+
  geom_point(aes(shape=indorpop)) +
  scale_y_continuous(trans='log10')+
  stat_cor(label.x=5,label.y=log10(max(df$mantaindel_rates)*0.94), aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of indels per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
translocvdepth<-ggplot(df, (aes(x=AvgDepth, y=transloc_rates, group=locationsfull, colour=locationsfull)))+
  geom_point(aes(shape=indorpop)) +
  scale_y_continuous(trans='log10')+
  stat_cor(label.x=5,label.y=log10(max(df$transloc_rates)*0.94), aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  ylab("# of translocations per bp")+
  xlab("Average Depth")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20))
#dupvdepth | mantaindelvdepth | translocvdepth
####this is fig : 
snpsvdepth + theme(legend.position = "none")| dupvdepth+ theme(legend.position = "none") | mantaindelvdepth+ theme(legend.position = "none") | translocvdepth + scale_color_discrete(name="Location") + theme(legend.text=element_text(size=15))
#df_boxun_snps + theme(legend.position = "none")| df_boxun_dup + theme(legend.position = "none")| df_boxun_mantaindel + theme(legend.position = "none")| df_boxun_transloc + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + theme(legend.text=element_text(size=15))
fit_snps_pop<-lm(snps_rates~ AvgDepth, data=justpop)
fit_dup_pop<-lm(dup_rates~ AvgDepth, data=justpop)
fit_mantaindel_pop<-lm(mantaindel_rates~ AvgDepth, data=justpop)
fit_transloc_pop<-lm(transloc_rates~ AvgDepth, data=justpop)
fit_snps_all<-lm(snps_rates~ AvgDepth, data=df)
fit_dup_all<-lm(dup_rates~ AvgDepth, data=df)
fit_mantaindel_all<-lm(mantaindel_rates~ AvgDepth, data=df)
fit_transloc_all<-lm(transloc_rates~ AvgDepth, data=df)
df$fit_snps.resid <-fit_snps_all$residuals
df$fit_dup.resid <-fit_dup_all$residuals
df$fit_mantaindel.resid <-fit_mantaindel_all$residuals
df$fit_transloc.resid <-fit_transloc_all$residuals
##box with resids:
df_boxun_dup<-ggplot(df, aes(x=treatment, y=fit_dup.resid))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("Duplications~Average Depth Residuals")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  #scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_transloc<-ggplot(df, aes(x=treatment, y=fit_transloc.resid))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("Translocations~Average Depth Residuals")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  #scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_mantaindel<-ggplot(df, aes(x=treatment, y=fit_mantaindel.resid))+
  geom_boxplot(aes(color=indorpop))+
  ylab("Indels~Average Depth Residuals")+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  #scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_snps<-ggplot(df, aes(x=treatment, y=fit_snps.resid))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("SNV~Average Depth Residuals")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  #scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
#this is figure:
df_boxun_snps + theme(legend.position = "none")| df_boxun_dup + theme(legend.position = "none")| df_boxun_mantaindel + theme(legend.position = "none")| df_boxun_transloc + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + theme(legend.text=element_text(size=15))

##box with treatment (rates):
df_boxun_dup<-ggplot(df, aes(x=treatment, y=dup_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("# duplications/bp")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_transloc<-ggplot(df, aes(x=treatment, y=transloc_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("# translocations/bp")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_mantaindel<-ggplot(df, aes(x=treatment, y=mantaindel_rates))+
  geom_boxplot(aes(color=indorpop))+
  ylab("# indels/bp")+
  theme_bw()+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
df_boxun_snps<-ggplot(df, aes(x=treatment, y=snps_rates))+
  geom_boxplot(aes(color=indorpop))+
  theme_bw()+
  ylab("# SNVs/bp")+
  geom_point(aes(color=indorpop, group=indorpop), position = position_dodge(width = 0.75))+
  scale_y_continuous(trans='log10')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=20,angle=90, hjust=1),
        axis.title.x=element_blank())
#this is figure:
Figure5<-(df_boxun_snps + theme(legend.position = "none")| df_boxun_dup + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + theme(legend.text=element_text(size=15)))/ (df_boxun_mantaindel + theme(legend.position = "none")| df_boxun_transloc + scale_color_discrete(name="Type",labels=c("Individual","Pairwise")) + theme(legend.text=element_text(size=15)))
Figure5+ plot_annotation(tag_levels = 'a')
##aovs by treatment
leveneTest(snps_rates ~ treatment, justind)
leveneTest(dup_rates ~ treatment, justind)
leveneTest(mantaindel_rates ~ treatment, justind)
leveneTest(transloc_rates ~ treatment, justind)

fit = aov(snps_rates ~ treatment, justind)
summary(fit)
fit2_snps=aov(snps_rates~treatment+AvgDepth,justind)
summary.lm(fit2_snps)
fit3_snps=aov(I(snps_rates*1e6)~treatment+AvgDepth,justind)
Anova(fit3_snps, type="III")
fit2_dup=aov(dup_rates~treatment+AvgDepth,justind)
summary.lm(fit2_dup)
fit3_dup=aov(I(dup_rates*1e6)~treatment+AvgDepth,justind)
Anova(fit3_dup, type="III")

fit2_mantaindel=aov(mantaindel_rates~treatment+AvgDepth,justind)
summary.lm(fit2_mantaindel)
fit3_mantaindel=aov(I(mantaindel_rates*1e6)~treatment+AvgDepth,justind)
Anova(fit3_mantaindel, type="III")

fit2_transloc=aov(transloc_rates~treatment+AvgDepth,justind)
summary.lm(fit2_transloc)
fit3_transloc=aov(I(transloc_rates*1e6)~treatment+AvgDepth,justind)
Anova(fit3_transloc, type="III")
########
fit_snps<-lm(snps_rates~ AvgDepth, data=justind)
justind$fit_snps.resid <-fit_snps$residuals
justind$size<-c(bikinisize, enewetaksizes,206,ofusizes,palausizes)
#justind$colony<-colony
#df$fit_indels.resid <-fit_indels$residuals

snpsvdist<-ggplot(justind, aes(size, snps_rates, colour=locationsfull)) + geom_point()+
  geom_point() +#geom_text(aes(label=colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=max(justind$snps_rates)*1.1, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("# SNVs/bp")+ xlab("Colony Diameter (cm)")+
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
  
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()#+
  #xlim(0,220)
dupvdist<-ggplot(justind, aes(size, dup_rates, colour=locationsfull)) + geom_point()+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justind$dup_rates)*1.1, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("# duplications/bp")+ xlab("Colony Diameter (cm)")+
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

translocvdist<-ggplot(justind, aes(size, transloc_rates, colour=locationsfull)) + geom_point()+
  geom_point() +
  stat_cor(label.x=5,label.y=max(justind$transloc_rates)*1.1, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("# translocations/bp")+ xlab("Colony Diameter (cm)")+
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

mantaindelvdist<-ggplot(justind, aes(size, mantaindel_rates, colour=locationsfull)) + geom_point()+
  geom_point() +
  stat_cor(aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x=5,label.y=max(justind$mantaindel_rates)*1.1) +
  geom_smooth(method=lm, aes(group=1)) +
  ylab("# indels/bp")+xlab("Colony Diameter (cm)")+
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#this is fig:
snpsvdist | dupvdist | mantaindelvdist | mantaindelvdist 
Figure3<-(snpsvdist + theme(legend.position = "none")| dupvdist + 
            scale_color_discrete(name="Location",labels=c("Bikini","Enewetak","Ofu","Palau")) + 
            theme(legend.text=element_text(size=15)))/( mantaindelvdist + theme(legend.position = "none")|
                                                          mantaindelvdist + scale_color_discrete(name="Location",labels=c("Bikini","Enewetak","Ofu","Palau")) +
                                                          theme(legend.text=element_text(size=15)))
Figure3 + plot_annotation(tag_levels = 'a')
###VAF
AHE11vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHE11-12_mf_FOURsnps.txt_melted.txt")
AHE09vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHE09-10_mf_FOURsnps.txt_melted.txt")
AHE07vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHE07-08_mf_FOURsnps.txt_melted.txt")
AHE01vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHE01-02_mf_FOURsnps.txt_melted.txt")
AHE01vaf / AHE07vaf / AHE09vaf / AHE11vaf
AHAS41vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHAS41-45_mf_FOURsnps.txt_melted.txt")
AHAS46vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHAS46-50_mf_FOURsnps.txt_melted.txt")
AHAS51vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHAS51-55_mf_FOURsnps.txt_melted.txt")
AHAS56vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHAS56-60_mf_FOURsnps.txt_melted.txt")
AHAS61vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHAS61-65_mf_FOURsnps.txt_melted.txt")
#(AHE11vaf | AHE09vaf | AHE07vaf | AHE01vaf)/(AHAS41vaf | AHAS46vaf | AHAS51vaf | AHAS56vaf | AHAS61vaf)
AHB125vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB125-129_mf_FOURsnps.txt_melted.txt")
AHB145vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB145-149_mf_FOURsnps.txt_melted.txt")
AHB151vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB151-155_mf_FOURsnps.txt_melted.txt")
AHB176vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB176-180_mf_FOURsnps.txt_melted.txt")
AHB195vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB195-199_mf_FOURsnps.txt_melted.txt")
AHB70vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB70-74_mf_FOURsnps.txt_melted.txt")
AHB90vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHB90-94_mf_FOURsnps.txt_melted.txt")

AHP01vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHP01-05_mf_FOURsnps.txt_melted.txt")
AHP06vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHP06-10_mf_FOURsnps.txt_melted.txt")
AHP16vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHP16-20_mf_FOURsnps.txt_melted.txt")
AHP21vaf<-AHE_soma_func("/Users/eloralopez/Documents/GitHub/Bikini/mappedtoahya/AHP21-25_mf_FOURsnps.txt_melted.txt")
(AHE11vaf | AHE09vaf | AHE07vaf | AHE01vaf)/(AHAS41vaf | AHAS46vaf | AHAS51vaf | AHAS56vaf | AHAS61vaf) /(AHB125vaf | AHB145vaf |AHB151vaf |AHB176vaf |AHB195vaf |AHB70vaf |AHB90vaf)/(AHP01vaf | AHP06vaf | AHP16vaf | AHP21vaf)
(AHB125vaf | AHB145vaf |AHB151vaf |AHB176vaf |AHB195vaf |AHB70vaf |AHB90vaf)/(AHP01vaf | AHP06vaf | AHP16vaf | AHP21vaf)
(AHB125vaf / AHB145vaf /AHB151vaf /AHB176vaf /AHB195vaf /AHB70vaf /AHB90vaf)
(AHP01vaf / AHP06vaf / AHP16vaf / AHP21vaf)
(AHAS41vaf / AHAS46vaf / AHAS51vaf / AHAS56vaf / AHAS61vaf)
#justind stats
bikini_ji<-subset(justind, locationsfull=="Bikini")
ene_ji<-subset(justind, locationsfull=="Enewetak")
palau_ji<-subset(justind, locationsfull=="Palau")
ofu_ji<-subset(justind, locationsfull=="Ofu")
snpsAnova <- lm(snps_rates ~ locationsfull, data = justind)
anova(snpsAnova)
dupAnova <- lm(dup_rates ~ locationsfull, data = justind)
anova(dupAnova)
translocAnova <- lm(transloc_rates ~ locationsfull, data = justind)
anova(translocAnova)
mantaindelAnova <- lm(mantaindel_rates ~ locationsfull, data = justind)
anova(mantaindelAnova)
#t-test by treatment
irradiatedset<-subset(justind, treatment=="irradiated")
notirradiatedset<-subset(justind, treatment=="not irradiated")
t.test(irradiatedset$snps_rates, notirradiatedset$snps_rates)
t.test(irradiatedset$dup_rates, notirradiatedset$dup_rates)
t.test(irradiatedset$transloc_rates, notirradiatedset$transloc_rates)
t.test(irradiatedset$mantaindel_rates, notirradiatedset$mantaindel_rates)
#ancova
library(car)
leveneTest(snps_rates~locationsfull, justind)
leveneTest(dup_rates~locationsfull, justind)
leveneTest(mantaindel_rates~locationsfull, justind)
leveneTest(transloc_rates~locationsfull, justind)

fit = aov(snps_rates ~ locationsfull, justind)
summary(fit)
fit2_snps=aov(snps_rates~locationsfull+AvgDepth,justind)
summary.lm(fit2_snps)
fit3_snps=aov(I(snps_rates*1e6)~locationsfull+AvgDepth,justind)
Anova(fit3_snps, type="III")
fit2_dup=aov(dup_rates~locationsfull+AvgDepth,justind)
summary.lm(fit2_dup)
fit3_dup=aov(I(dup_rates*1e6)~locationsfull+AvgDepth,justind)
Anova(fit3_dup, type="III")

fit2_mantaindel=aov(mantaindel_rates~locationsfull+AvgDepth,justind)
summary.lm(fit2_mantaindel)
fit3_mantaindel=aov(I(mantaindel_rates*1e6)~locationsfull+AvgDepth,justind)
Anova(fit3_mantaindel, type="III")

fit2_transloc=aov(transloc_rates~locationsfull+AvgDepth,justind)
summary.lm(fit2_transloc)
fit3_transloc=aov(I(transloc_rates*1e6)~locationsfull+AvgDepth,justind)
Anova(fit3_transloc, type="III")

leveneTest(snps_rates~locationsfull, justpop)
leveneTest(dup_rates~locationsfull, justpop)
leveneTest(mantaindel_rates~locationsfull, justpop)
leveneTest(transloc_rates~locationsfull, justpop)
fit2_snps=aov(snps_rates~locationsfull+AvgDepth,justpop)
summary.lm(fit2_snps)
fit3_snps=aov(I(snps_rates*1e6)~locationsfull+AvgDepth,justpop)
Anova(fit3_snps, type="III")
fit2_dup=aov(dup_rates~locationsfull+AvgDepth,justpop)
summary.lm(fit2_dup)
fit3_dup=aov(I(dup_rates*1e6)~locationsfull+AvgDepth,justpop)
Anova(fit3_dup, type="III")

fit2_mantaindel=aov(mantaindel_rates~locationsfull+AvgDepth,justpop)
summary.lm(fit2_mantaindel)
fit3_mantaindel=aov(I(mantaindel_rates*1e6)~locationsfull+AvgDepth,justpop)
Anova(fit3_mantaindel, type="III")

fit2_transloc=aov(transloc_rates~locationsfull+AvgDepth,justpop)
summary.lm(fit2_transloc)
fit3_transloc=aov(I(transloc_rates*1e6)~locationsfull+AvgDepth,justpop)
Anova(fit3_transloc, type="III")

##ancova rate vs size justind
leveneTest(snps_rates~locationsfull, justind)
leveneTest(dup_rates~locationsfull, justind)
leveneTest(mantaindel_rates~locationsfull, justind)
leveneTest(transloc_rates~locationsfull, justind)

fit = aov(snps_rates ~ locationsfull, justind)
summary(fit)
fit2_snps=aov(snps_rates~locationsfull+size,justind)
summary.lm(fit2_snps)
fit3_snps=aov(I(snps_rates*1e6)~locationsfull+size,justind)
Anova(fit3_snps, type="III")
fit2_dup=aov(dup_rates~locationsfull+size,justind)
summary.lm(fit2_dup)
fit3_dup=aov(I(dup_rates*1e6)~locationsfull+size,justind)
Anova(fit3_dup, type="III")

fit2_mantaindel=aov(mantaindel_rates~locationsfull+size,justind)
summary.lm(fit2_mantaindel)
fit3_mantaindel=aov(I(mantaindel_rates*1e6)~locationsfull+AvgDepth,justind)
Anova(fit3_mantaindel, type="III")

fit2_transloc=aov(transloc_rates~locationsfull+size,justind)
summary.lm(fit2_transloc)
fit3_transloc=aov(I(transloc_rates*1e6)~locationsfull+AvgDepth,justind)
Anova(fit3_transloc, type="III")
##with vs between colony differences
Bikinisnps<-aov(subset(df, locationsfull=="Bikini")$snps_rates~subset(df, locationsfull=="Bikini")$indorpop)
Bikinisnps2<-aov(subset(df, locationsfull=="Bikini")$snps_rates~subset(df, locationsfull=="Bikini")$indorpop+subset(df,locationsfull=="Bikini")$AvgDepth)

summary(Biksnps)
Anova(Bikinisnps2, type="III")
anova(Bikinisnps, Biksnps2)
leveneTest(subset(df, locationsfull=="Bikini")$snps_rates~subset(df, locationsfull=="Bikini")$indorpop)
leveneTest(subset(df, locationsfull=="Bikini")$dup_rates~subset(df, locationsfull=="Bikini")$indorpop)
leveneTest(subset(df, locationsfull=="Bikini")$mantaindel_rates~subset(df, locationsfull=="Bikini")$indorpop)
leveneTest(subset(df, locationsfull=="Bikini")$transloc_rates~subset(df, locationsfull=="Bikini")$indorpop)

Bikinidup<-aov(subset(df, locationsfull=="Bikini")$dup_rates~subset(df, locationsfull=="Bikini")$indorpop)
Bikinidup2<-aov(I(subset(df, locationsfull=="Bikini")$dup_rates*1e6)~subset(df, locationsfull=="Bikini")$indorpop+subset(df,locationsfull=="Bikini")$AvgDepth)
summary(Bikinidup)
Anova(Bikinidup2, type="III")
anova(Bikinidup, Bikinidup2)

Bikinisnps2<-aov(subset(df, locationsfull=="Bikini")$snps_rates~subset(df, locationsfull=="Bikini")$indorpop+subset(df,locationsfull=="Bikini")$AvgDepth)

summary(Biksnps)
Anova(Bikinisnps2, type="III")
anova(Bikinisnps, Biksnps2)

Bikinisnps<-aov(subset(df, locationsfull=="Bikini")$snps_rates~subset(df, locationsfull=="Bikini")$indorpop)
Bikinisnps2<-aov(subset(df, locationsfull=="Bikini")$snps_rates~subset(df, locationsfull=="Bikini")$indorpop+subset(df,locationsfull=="Bikini")$AvgDepth)

summary(Biksnps)
Anova(Bikinisnps2, type="III")
anova(Bikinisnps, Biksnps2)
#t.test(bikini_ji$snps_rates, subset(justpop, locationsfull=="Bikini")$snps_rates)
t.test(palau_ji$snps_rates, subset(justpop, locationsfull=="Palau")$snps_rates)
t.test(ofu_ji$snps_rates, subset(justpop, locationsfull=="Ofu")$snps_rates)
t.test(ene_ji$snps_rates, subset(justpop, locationsfull=="Enewetak")$snps_rates)

t.test(bikini_ji$dup_rates, subset(justpop, locationsfull=="Bikini")$dup_rates)
t.test(palau_ji$dup_rates, subset(justpop, locationsfull=="Palau")$dup_rates)
t.test(ofu_ji$dup_rates, subset(justpop, locationsfull=="Ofu")$dup_rates)
t.test(ene_ji$dup_rates, subset(justpop, locationsfull=="Enewetak")$dup_rates)

t.test(bikini_ji$mantaindel_rates, subset(justpop, locationsfull=="Bikini")$mantaindel_rates)
t.test(palau_ji$mantaindel_rates, subset(justpop, locationsfull=="Palau")$mantaindel_rates)
t.test(ofu_ji$mantaindel_rates, subset(justpop, locationsfull=="Ofu")$mantaindel_rates)
t.test(ene_ji$mantaindel_rates, subset(justpop, locationsfull=="Enewetak")$mantaindel_rates)

t.test(bikini_ji$transloc_rates, subset(justpop, locationsfull=="Bikini")$transloc_rates)
t.test(palau_ji$transloc_rates, subset(justpop, locationsfull=="Palau")$transloc_rates)
t.test(ofu_ji$transloc_rates, subset(justpop, locationsfull=="Ofu")$transloc_rates)
t.test(ene_ji$transloc_rates, subset(justpop, locationsfull=="Enewetak")$transloc_rates)

##two-way ancova
#linearity assumption (already shown in all the other plots)
#homogeneity of regression slopes:
homg<-aov(snps_rates~AvgDepth+indorpop+locationsfull+indorpop*locationsfull+AvgDepth*indorpop+
      AvgDepth*locationsfull+AvgDepth*indorpop*locationsfull, df)
summary(homg)
df$group<-factor(paste0(df$indorpop, df$locationsfull))
homg2<-aov(snps_rates~group*AvgDepth,df)
summary(homg2)
model<-lm(snps_rates~AvgDepth + indorpop * locationsfull, data=df)
summary(model)

mean(subset(justpop, locationsfull=="Bikini")$AvgDepth)
mean(subset(justind, locationsfull=="Bikini")$AvgDepth)
t.test(subset(justpop, locationsfull=="Bikini")$AvgDepth, 
       subset(justind, locationsfull=="Bikini")$AvgDepth, alternative = c("two.sided"))
t.test(subset(justpop, locationsfull=="Palau")$AvgDepth, 
       subset(justind, locationsfull=="Palau")$AvgDepth, alternative = c("two.sided"))
t.test(subset(justpop, locationsfull=="Ofu")$AvgDepth, 
       subset(justind, locationsfull=="Ofu")$AvgDepth, alternative = c("two.sided"))
t.test(subset(justpop, locationsfull=="Enewetak")$AvgDepth, 
       subset(justind, locationsfull=="Enewetak")$AvgDepth, alternative = c("two.sided"))
