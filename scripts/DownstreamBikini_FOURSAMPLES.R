setwd("~/Documents/GitHub/Bikini/forDownstream")
#FOURSAMPLEdenoms
AHAS41FOURSAMPLEdenom<-139034873-(895*14)
AHAS46FOURSAMPLEdenom<-150062856-(895*14)
AHAS51FOURSAMPLEdenom<-133379132-(895*14)
AHAS56FOURSAMPLEdenom<-153220786-(895*14)
AHAS61FOURSAMPLEdenom<-135529784-(895*14)

AHB125FOURSAMPLEdenom<-142568826-(895*14)
AHB145FOURSAMPLEdenom<-120296114-(895*14)
AHB151FOURSAMPLEdenom<-123029507-(895*14)
AHB176FOURSAMPLEdenom<-134003833-(895*14)
AHB195FOURSAMPLEdenom<-141499547-(895*14)
AHB70FOURSAMPLEdenom<-135674761-(895*14)
AHB90FOURSAMPLEdenom<-131210271-(895*14)

AHP01FOURSAMPLEdenom<-120996343-(895*14)
AHP06FOURSAMPLEdenom<-133809818-(895*14)
AHP16FOURSAMPLEdenom<-130687954-(895*14)
AHP21FOURSAMPLEdenom<-126197489-(895*14)

AHE01FOURSAMPLEdenom<-AHE01denom
AHE07FOURSAMPLEdenom<-AHE07denom
AHE09FOURSAMPLEdenom<-AHE09denom
AHE11FOURSAMPLEdenom<-AHE11denom
mean(c(AHB125FOURSAMPLEdenom, AHB145FOURSAMPLEdenom, AHB151FOURSAMPLEdenom,AHB176FOURSAMPLEdenom,AHB195FOURSAMPLEdenom, AHB70FOURSAMPLEdenom,AHB90FOURSAMPLEdenom,
       AHAS41FOURSAMPLEdenom, AHAS46FOURSAMPLEdenom, AHAS51FOURSAMPLEdenom, AHAS56FOURSAMPLEdenom, AHAS61FOURSAMPLEdenom, 
       AHP01FOURSAMPLEdenom, AHP06FOURSAMPLEdenom, AHP16FOURSAMPLEdenom, AHP21FOURSAMPLEdenom, 
       AHE01FOURSAMPLEdenom, AHE07FOURSAMPLEdenom, AHE09FOURSAMPLEdenom, AHE11FOURSAMPLEdenom))
se(c(AHB125FOURSAMPLEdenom, AHB145FOURSAMPLEdenom, AHB151FOURSAMPLEdenom,AHB176FOURSAMPLEdenom,AHB195FOURSAMPLEdenom, AHB70FOURSAMPLEdenom,AHB90FOURSAMPLEdenom,
       AHAS41FOURSAMPLEdenom, AHAS46FOURSAMPLEdenom, AHAS51FOURSAMPLEdenom, AHAS56FOURSAMPLEdenom, AHAS61FOURSAMPLEdenom, 
       AHP01FOURSAMPLEdenom, AHP06FOURSAMPLEdenom, AHP16FOURSAMPLEdenom, AHP21FOURSAMPLEdenom, 
       AHE01FOURSAMPLEdenom, AHE07FOURSAMPLEdenom, AHE09FOURSAMPLEdenom, AHE11FOURSAMPLEdenom))
#average depths
AHAS41avg<-
  
indels_files<-list.files(path="/Users/eloralopez/Documents/GitHub/Bikini/forDownstream", pattern="_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)
for (i in 1:length(indels_files)) {
  #meltfunc(i)
  file =files[i]
  print(file)
  indels_df<-AHE_soma_func(file)
  #return(indels_df)
}
AHE_soma_func<- function(file) {
  #files<-list.files(path="~/Documents/GitHub/Bikini", pattern="AHE01-02_mT_20201215.txt.ann.txt", full.names=T, recursive=FALSE)
  #file="AHE11-12_mT_indels.txt.ann.txt"
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
  allelesplit<-str_split_fixed(genotypes, "/", 2) #splits the normal and mutant alleles into different strings
  allele1<-allelesplit[,1]
  allele2<-allelesplit[,2]
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
                         
                         "genotype"= genotypes,"allele1"=allele1, "allele2"=allele2, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	
                         
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
  
  uniquesomaticmetadatadf<-metadatadf.0[match(unique(metadatadf.0$chrom.pos), 					metadatadf.0$chrom.pos),]
  
  justhomozyg<-subset(metadatadf.0, allele1==allele2)
  truejusthomozyg<-subset(justhomozyg, refdepth == 0 | altdepth ==0)
  truejusthomozyg_freq<-as.data.frame(table(truejusthomozyg$chrom.pos))
  all_homozyg<-subset(truejusthomozyg_freq, Freq==2)
  all_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  uniqueall_homozyg_fulldf<-metadatadf.0[match(all_homozyg$Var1,metadatadf.0$chrom.pos),]
  return(uniqueall_homozyg_fulldf)
  
}

##COMPARE AGAINST ENEWETAK: return(uniqueall_homozyg_fulldf) return(allsinglesample)
indels_AHE01<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE01-02_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE01FOURSAMPLEdenom
#AHE05<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE05-06_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE05FOURSAMPLEdenom
indels_AHE07<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE07-08_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE07FOURSAMPLEdenom
indels_AHE09<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE09-10_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE09FOURSAMPLEdenom
indels_AHE11<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE11-12_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHE11FOURSAMPLEdenom

indels_enerates<-c(indels_AHE01,indels_AHE07,indels_AHE09,indels_AHE11)

indels_ahb125<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB125-129_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB125FOURSAMPLEdenom
indels_ahb145<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB145-149_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB145FOURSAMPLEdenom
indels_ahb151<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB151-155_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB151FOURSAMPLEdenom
indels_ahb176<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB176-180_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB176FOURSAMPLEdenom
indels_ahb195<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB195-199_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB195FOURSAMPLEdenom
indels_ahb70<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB70-74_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB70FOURSAMPLEdenom
indels_ahb90<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB90-94_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB90FOURSAMPLEdenom
indels_bikrates<-c(indels_ahb70, indels_ahb90, indels_ahb125, indels_ahb145, indels_ahb151, indels_ahb176, indels_ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
indels_AHAS41<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS41-45_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS41FOURSAMPLEdenom
indels_AHAS46<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS46-50_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS46FOURSAMPLEdenom
indels_AHAS51<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS51-55_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS51FOURSAMPLEdenom
indels_AHAS56<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS56-60_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS56FOURSAMPLEdenom
indels_AHAS61<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS61-65_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS61FOURSAMPLEdenom

indels_ofurates<-c(indels_AHAS41,indels_AHAS46, indels_AHAS51, indels_AHAS56, indels_AHAS61)
ofusizes_with41<-c(206,ofusizes)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

indels_AHP01<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP01-05_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP01FOURSAMPLEdenom
indels_AHP06<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP06-10_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP06FOURSAMPLEdenom
indels_AHP16<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP16-20_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP16FOURSAMPLEdenom
indels_AHP21<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP21-25_mT_indels.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP21FOURSAMPLEdenom

indels_palrates<-c(indels_AHP01,indels_AHP06,indels_AHP16,indels_AHP21)

indels_irradiated<-c(indels_enerates,indels_bikrates)
indels_notirradiated<-c(indels_ofurates,indels_palrates)

indels_averages<-c(mean(indels_bikrates),mean(indels_enerates),mean(indels_ofurates),mean(indels_palrates))
indels_rates<-c(indels_bikrates, indels_enerates, indels_ofurates,indels_palrates)
indels_ses<-c(se(indels_bikrates),se(indels_enerates),se(indels_ofurates),se(indels_palrates))
location<-c(rep("Bikini",7),rep("Enewetak",4),rep("Ofu",5),rep("Palau",4))
sizes<-c(bikinisize,enewetaksizes,ofusizes,palausizes)
colony<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "ahe01","ahe07","ahe09","ahe11","AHAS41","AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")

indels_ratedf<-data.frame("Rate"=indels_rates, "Location"=location, "Distance"=sizes,"Colony"=colony)
#rate vs location:
indels_boxun<-ggplot(indels_ratedf, aes(x=Location, y=Rate))+
  geom_boxplot()+
  theme_bw()
#ratevs size for muts UNIQUE TO ONE SAMPLE:
indels_un<-ggplot(indels_ratedf, (aes(x=Distance, y=Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()
######SVs##########
SVs_AHE01<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE01-02_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHE01FOURSAMPLEdenom
#AHE05<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE05-06_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHE05FOURSAMPLEdenom
SVs_AHE07<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE07-08_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHE07FOURSAMPLEdenom
SVs_AHE09<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE09-10_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHE09FOURSAMPLEdenom
SVs_AHE11<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHE11-12_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHE11FOURSAMPLEdenom

SVs_enerates<-c(SVs_AHE01,SVs_AHE07,SVs_AHE09,SVs_AHE11)
SVs_ahb125<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB125-129_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB125FOURSAMPLEdenom
SVs_ahb145<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB145-149_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB145FOURSAMPLEdenom
SVs_ahb151<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB151-155_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB151FOURSAMPLEdenom
SVs_ahb176<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB176-180_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB176FOURSAMPLEdenom
SVs_ahb195<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB195-199_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB195FOURSAMPLEdenom
SVs_ahb70<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB70-74_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB70FOURSAMPLEdenom
SVs_ahb90<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB90-94_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHB90FOURSAMPLEdenom
SVs_bikrates<-c(SVs_ahb70, SVs_ahb90, SVs_ahb125, SVs_ahb145, SVs_ahb151, SVs_ahb176, SVs_ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
SVs_AHAS41<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS41-45_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHAS41FOURSAMPLEdenom
SVs_AHAS46<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS46-50_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHAS46FOURSAMPLEdenom
SVs_AHAS51<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS51-55_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHAS51FOURSAMPLEdenom
SVs_AHAS56<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS56-60_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHAS56FOURSAMPLEdenom
SVs_AHAS61<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS61-65_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHAS61FOURSAMPLEdenom

SVs_ofurates<-c(SVs_AHAS41,SVs_AHAS46, SVs_AHAS51, SVs_AHAS56, SVs_AHAS61)
ofusizes_with41<-c(206,ofusizes)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

SVs_AHP01<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP01-05_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHP01FOURSAMPLEdenom
SVs_AHP06<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP06-10_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHP06FOURSAMPLEdenom
SVs_AHP16<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP16-20_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHP16FOURSAMPLEdenom
SVs_AHP21<-SV_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP21-25_mT_sv20201226.txt.ann.txt", full.names=T, recursive=FALSE))/AHP21FOURSAMPLEdenom

SVs_palrates<-c(SVs_AHP01,SVs_AHP06,SVs_AHP16,SVs_AHP21)

SVs_irradiated<-c(SVs_enerates,SVs_bikrates)
SVs_notirradiated<-c(SVs_ofurates,SVs_palrates)
SVs_averages<-c(mean(SVs_bikrates),mean(SVs_enerates),mean(SVs_ofurates),mean(SVs_palrates))
SVs_rates<-c(SVs_bikrates, SVs_enerates, SVs_ofurates,SVs_palrates)
SVs_ses<-c(se(SVs_bikrates),se(SVs_enerates),se(SVs_ofurates),se(SVs_palrates))
location<-c(rep("Bikini",7),rep("Enewetak",4),rep("Ofu",5),rep("Palau",4))
sizes<-c(bikinisize,enewetaksizes,ofusizes,palausizes)
colony<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "ahe01","ahe07","ahe09","ahe11","AHAS41","AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")

SVs_ratedf<-data.frame("Rate"=SVs_rates, "Location"=location, "Distance"=sizes,"Colony"=colony)
#rate vs location:
SVs_boxun<-ggplot(SVs_ratedf, aes(x=Location, y=Rate))+
  geom_boxplot()+
  theme_bw()
#ratevs size for muts UNIQUE TO ONE SAMPLE:
SVs_un<-ggplot(SVs_ratedf, (aes(x=Distance, y=Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()
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

snps_palrates<-c(snps_AHP01,snps_AHP06,snps_AHP16,snps_AHP21)

snps_irradiated<-c(snps_enerates,snps_bikrates)
snps_notirradiated<-c(snps_ofurates,snps_palrates)

snps_averages<-c(mean(snps_bikrates),mean(snps_enerates),mean(snps_ofurates),mean(snps_palrates))
snps_rates<-c(snps_bikrates, snps_enerates, snps_ofurates,snps_palrates)
snps_ses<-c(se(snps_bikrates),se(snps_enerates),se(snps_ofurates),se(snps_palrates))
depths<-c(AHB70avg, AHB90avg, AHB125avg, AHB145avg, AHB151avg, AHB176avg, AHB195avg, AHE01avg, AHE07avg, AHE09avg, AHE11avg,
          AHAS41avg, AHAS46avg, AHAS51avg, AHAS56avg, AHAS61avg, AHP01avg, AHP06avg, AHP16avg, AHP21avg)
location<-c(rep("Bikini",7),rep("Enewetak",4),rep("Ofu",5),rep("Palau",4))
sizes<-c(bikinisize,enewetaksizes,ofusizes,palausizes)
colony<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "ahe01","ahe07","ahe09","ahe11","AHAS41","AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")
normdepths<-(depths - min(depths)) / (max(depths)-min(depths))
snps_ratedf<-data.frame("Rate"=snps_rates, "Location"=location, "Distance"=sizes,"Colony"=colony)
treatment<-c(rep("irradiated",11),rep("not irradiated",9))
siratio<-snps_rates/indels_rates
fit_snps<-lm(snps_rates~ depths)
#ratedf$fit_snps.resid <-fit_snps$residuals
fit_indels<-lm(indels_rates~ depths)
fit_SVs<-lm(SVs_rates~depths)
#ratedf$fit_indels.resid <-fit_indels$residuals
#resid.ratio<-fit_snps$residuals/fit_indels$residuals
ratedf<-data.frame("snps_Rate"=snps_rates,"indels_Rate"=indels_rates, "SVs_Rate" = SVs_rates, "Location"=location, 
                   "Distance"=sizes,"Colony"=colony, "AvgDepth"=depths, "Treatment"=treatment, 
                    "Indels.Resid"=fit_indels$residuals,"Snps.Resid"=fit_snps$residuals, "SVs.Resid"=fit_SVs$residuals)

#rate vs location:
snps_boxun<-ggplot(snps_ratedf, aes(x=Location, y=Rate))+
  geom_boxplot()+
  theme_bw()
#ratevs size for muts UNIQUE TO ONE SAMPLE:
snps_un<-ggplot(ratedf, (aes(x=Distance, y=SVs_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()+
  #
  stat_cor(label.x=100,label.y=0.00001, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  stat_regline_equation(label.x = 100, label.y = 0.000005)
snps_un1<-ggplot(ratedf, (aes(x=Distance, y=snps_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()

#ratios
siratiopoint<-ggplot(ratedf, (aes(x=Distance, y=SIratio, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()
siratiobox<-ggplot(ratedf, aes(x=Location, y=SIratio, group=Location, colour=Location))+
  geom_boxplot()+
  theme_bw()
residratiobox<-ggplot(ratedf, aes(x=Location, y=ResidRatio, group=Location, colour=Location))+
  geom_boxplot()+
  theme_bw()
summary(anova(lm(ResidRatio~Location,data=ratedf)))
indelbox<-ggplot(ratedf, aes(x=Location, y=Indels.Resid, group=Location, colour=Location))+
  geom_boxplot()+
  theme_bw()
snpsbox<-ggplot(ratedf, aes(x=Location, y=Snps.Resid, group=Location, colour=Location))+
  geom_boxplot()+
  theme_bw()
#snps vs indels rates
ggplot(ratedf, (aes(x=snps_Rate, y=indels_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()
#indels vs depth
indelsvdepth<-ggplot(ratedf, (aes(x=AvgDepth, y=indels_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=0.00004, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
#  stat_regline_equation(label.x = 5, label.y = 0.00003)+  
  theme_bw()#+
  xlim(0,210)
#snps vs depth
snpsvdepth<-ggplot(ratedf, (aes(x=AvgDepth, y=snps_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=0.000082, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
#  stat_regline_equation(label.x = 5, label.y = 0.000077)+    
  theme_bw()
#SVs vs depth
SVsvdepth<-ggplot(ratedf, (aes(x=AvgDepth, y=SVs_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=0.00003, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
#  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()
####this is fig : 
snpsvdepth | indelsvdepth | SVsvdepth
fit_snps<-lm(snps_Rate~ AvgDepth, data=ratedf)
ratedf$fit_snps.resid <-fit_snps$residuals
fit_indels<-lm(indels_Rate~ AvgDepth, data=ratedf)
ratedf$fit_indels.resid <-fit_indels$residuals

gg1<-ggplot(ratedf, aes(Distance, snps_Rate)) + geom_point()
gg2<- ggplot(ratedf, aes(AvgDepth, snps_Rate)) + geom_point()+
  stat_cor(label.x=20,label.y=0.00004, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  stat_regline_equation(label.x = 20, label.y = 0.00003)
snpsresidvdist<-ggplot(ratedf, aes(Distance, Snps.Resid, colour=Location)) + geom_point()+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=0.00003, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  xlim(0,220)
indelsresidvdist<-ggplot(ratedf, aes(Distance, Indels.Resid, colour=Location)) + geom_point()+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=0.00003, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  xlim(0,220)
SVsvdist<-ggplot(ratedf, aes(Distance, SVs_Rate, color=Location)) + geom_point()+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  stat_cor(label.x=5,label.y=0.00003, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  #  stat_regline_equation(label.x = 5, label.y = 0.000025)+      
  theme_bw()+
  xlim(0,220)
snpsresidvdist | indelsresidvdist | SVsvdist

gg1_snps<-ggplot(ratedf, aes(Location, snps_Rate, group=Location)) + geom_boxplot()
gg2<- ggplot(ratedf, aes(AvgDepth, snps_Rate)) + geom_point()
gg3_snps<-ggplot(ratedf, aes(Location, fit_snps.resid,  group=Location)) + geom_boxplot()
gg1 | gg2 | gg3
gg1_indels<-ggplot(ratedf, aes(Location, indels_Rate, group=Location)) + geom_boxplot()
gg2_indels<- ggplot(ratedf, aes(AvgDepth, indels_Rate)) + geom_point()+
  stat_cor(label.x=20,label.y=0.00004, aes(group=1)) +
  geom_smooth(method=lm, aes(group=1)) +
  stat_regline_equation(label.x = 20, label.y = 0.00003)
gg3_indels<-ggplot(ratedf, aes(Location, fit_indels.resid,  group=Location)) + geom_boxplot()
#snps vs dist when depth>20
ratedf20<-subset(ratedf, depths>20)
ggplot(ratedf20, (aes(x=Distance, y=snps_Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=Colony),hjust=0, vjust=0)+
  theme_bw()+
  xlim(0,210)

snps_boxun<-ggplot(ratedf20, aes(x=Location, y=snps_Rate))+
  geom_boxplot()+
  theme_bw()

indels_boxun<-ggplot(ratedf, aes(x=Location, y=snps_Rate/NormDepth))+
  geom_boxplot()+
  theme_bw()
#regressions
distance.snps<-lm(ratedf$snps_Rate~ratedf$Distance)
summary(distance.snps)
distance.depth.snps<-lm(formula = ratedf$snps_Rate~ ratedf$Distance + ratedf$AvgDepth)
summary(distance.depth.snps)
##anova
snpsAnova <- lm(snps_Rate/depths ~ Location, data = ratedf)
anova(snpsAnova)
snpsAnovaSummary <- summary(snpsAnova)
snpsAnovaSummary$r.squared
snpsAnova2<-lm(snps_Rate ~ depths*Location, data = ratedf)
anova(snpsAnova2)

residualsAnova_snps<-lm(fit_snps.resid~Location, data=ratedf)
anova(residualsAnova_snps)
residualsAnova_indels<-lm(fit_indels.resid~Location, data=ratedf)
anova(residualsAnova_indels)

residualsAnova_snps<-lm(fit_snps.resid~Treatment, data=ratedf)
summary(anova(residualsAnova_snps))
residualsAnova_indels<-lm(fit_indels.resid~Treatment, data=ratedf)
anova(residualsAnova_indels)

result1 <- aov(snps_Rate~Location*depths,data=ratedf)
summary(result1)
result2 <- aov(snps_Rate~Location+depths,data=ratedf)
summary(result2)
anova(result1,result2)

result1 <- aov(snps_Rate~Treatment*depths,data=ratedf)
summary(result1)
result2 <- aov(snps_Rate~Treatment+depths,data=ratedf)
summary(result2)
anova(result1,result2)
###########################
############################
###########################
ahb125<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB125-129_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB125FOURSAMPLEdenom
ahb145<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB145-149_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB145FOURSAMPLEdenom
ahb151<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB151-155_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB151FOURSAMPLEdenom
ahb176<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB176-180_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB176FOURSAMPLEdenom
ahb195<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB195-199_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB195FOURSAMPLEdenom
ahb70<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB70-74_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB70FOURSAMPLEdenom
ahb90<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHB90-94_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHB90FOURSAMPLEdenom
bikrates<-c(ahb70, ahb90, ahb125, ahb145, ahb151, ahb176, ahb195)
biksizes<-bikini$distcm
#ubikrates<-c(ahb125, ahb145, ahb151, ahb176, ahb195, ahb70, ahb90)
AHAS41<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS41-45_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS41FOURSAMPLEdenom
AHAS46<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS46-50_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS46FOURSAMPLEdenom
AHAS51<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS51-55_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS51FOURSAMPLEdenom
AHAS56<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS56-60_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS56FOURSAMPLEdenom
AHAS61<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHAS61-65_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHAS61FOURSAMPLEdenom

ofurates<-c(AHAS41,AHAS46, AHAS51, AHAS56, AHAS61)
ofusizes_with41<-c(206,ofusizes)
#uofurates<-c(AHAS46, AHAS51, AHAS56, AHAS61)

AHP01<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP01-05_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP01FOURSAMPLEdenom
AHP06<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP06-10_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP06FOURSAMPLEdenom
AHP16<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP16-20_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP16FOURSAMPLEdenom
AHP21<-nrow(AHE_soma_func(list.files(path="~/Documents/GitHub/Bikini/forDownstream", pattern="AHP21-25_mT_snps.txt.ann.txt", full.names=T, recursive=FALSE)))/AHP21FOURSAMPLEdenom

palrates<-c(AHP01,AHP06,AHP16,AHP21)

irradiated<-c(enerates,bikrates)
notirradiated<-c(ofurates,palrates)

averages<-c(mean(bikrates),mean(enerates),mean(ofurates),mean(palrates))
rates<-c(bikrates, enerates, ofurates,palrates)
ses<-c(se(bikrates),se(enerates),se(ofurates),se(palrates))
location<-c(rep("Bikini",7),rep("Enewetak",4),rep("Ofu",5),rep("Palau",4))
sizes<-c(bikinisize,enewetaksizes,ofusizes_with41,palausizes)
colony<-c("ahb70", "ahb90", "ahb125", "ahb145", "ahb151", "ahb176", "ahb195", "ahe01","ahe07","ahe09","ahe11","AHAS41","AHAS46", "AHAS51", "AHAS56", "AHAS61", "AHP01","AHP06","AHP16","AHP21")

ratedf<-data.frame("Rate"=rates, "Location"=location, "Distance"=sizes,"Colony"=colony)
#rate vs location:
boxun<-ggplot(ratedf, aes(x=Location, y=Rate))+
  geom_boxplot()+
  theme_bw()
#ratevs size for muts UNIQUE TO ONE SAMPLE:
un<-ggplot(ratedf, (aes(x=Distance, y=Rate, group=Location, colour=Location)))+
  geom_point() +geom_text(aes(label=colony),hjust=0, vjust=0)+
  theme_bw()