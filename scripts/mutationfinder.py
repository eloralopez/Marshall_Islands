#DESCRIPTION:this script calls two input files. Each input file contains genotype information for two replicate genomic DNA samples. The "fileinfo" function defined at the top of the script generates a list of all of the sites where the two replicate samples in a single file have the SAME GENOTYPE. Then, the two lists of matches are compared to each other and the final writeout gives just the sites where the genotypes differ between the two lists.

#USAGE: python3 mutationfinder20190820.py testAHAS57.txt testAHAS59.txt

import sys
import re


input1 = sys.argv[1]
numberofsamples = int(sys.argv[2])

def fileinfo(inputfile, rep1, rep2): 
    listofoutlist = [] #starts an empty list to which the matches between replicate libraries will be added later
    
    with open(inputfile, 'r') as f: #opens whichever file is specified
        
        for line in f:

            if line.startswith('##'): #ignores all the header lines
                continue

            if line.startswith('#CHROM'): #ignores all the header lines

                line = line.strip()

                items = line.split('\t')

                samplenames = items[9:]
                samplenames = '\t'.join(samplenames)
                #print(samplenames)
                header = ["chrom.pos", "ref", "alt", samplenames]
                #print(header)
                header = '\t'.join(header)
                listofoutlist.append(header)
                #print(header)
                #print(samplenames_out)
                

            else:

                line = line.strip()

                items = line.split('\t')

                contig = items[0] #the #CHROM column in the VCF

                position = items[1] #the POS column in the VCF

                concatenated = contig + "." + position #creates a string that contains both the chromosome and the position information

                ref_allele = items[3] #the REF column in the VCF


                alt_allele = items[4] #the ALT column in the VCF

                genotypes = items[9:] #each of the samplename columns in the VCF (e.g. AHAS57-1)

                geno_list=[] #empty list to which genotypes will be added later

                genoplusdepthlist=[] #empty list to which genotype depths will be added later
                #print(len(genotypes))
                for genotype in genotypes:

                    genos = genotype.split(':') #splits the genotype column up by ":"

                    alleles = re.split('; |, |\||\/',genos[0])#genos[0].split('/') #splits the first part of the genotype column into its two respective alleles
                    
                    alleledepths = genos[1].split(',') #depth per allele at locus
                    #FROM GATK: AD is the unfiltered allele depth, i.e. the number of reads that support each of the reported alleles. All reads at the position (including reads that did not pass the variant caller’s filters) are included in this number, except reads that were considered uninformative. Reads are considered uninformative when they do not provide enough statistical evidence to support one allele over another.

                    refdepth = alleledepths[0]
                    
                    altdepth = alleledepths[1]
                    
                    GQscore = genos[3]

#                    totaldepth_atlocus = genos[2] #totaldepth at locus

                    #FROM GATK: DP is the filtered depth, at the sample level. This gives you the number of filtered reads that support each of the reported alleles. You can check the variant caller’s documentation to see which filters are applied by default. Only reads that passed the variant caller’s filters are included in this number. However, unlike the AD calculation, uninformative reads are included in DP.
                    
                    if alleles[0] == ".": #ignore the sites where the depth is './.'
                        continue
                    if int(alleles[0]) == 0:

                        a1 = ref_allele #gives the corresponding character (A,C, G, or T) for the allele

        #
                    else:

                        a1 = alt_allele #gives the corresponding character (A,C, G, or T) for the allele

                    if int(alleles[1]) == 0:

                        a2 = ref_allele #gives the corresponding character (A,C, G, or T) for the allele

                    else:

                        a2 = alt_allele #gives the corresponding character (A,C, G, or T) for the allele
                    genostring = a1 + '/' + a2 #gives the genotype (e.g. A/T)
                    genoPLUSdepth = genostring + "," + refdepth + ","  + altdepth + "," + GQscore#gives a string of the genotype plus its corresponding depthsand GQscore

                    genoplusdepthlist.append(genoPLUSdepth)
                    genoPLUSdepth = '\t'.join(genoplusdepthlist)

                    geno_list.append(genostring)
                    genostring = '\t'.join(geno_list)
                    
                if geno_list[rep1] == geno_list[rep2]: #outputs just the sites where the two replicate libraries are the same genotype
                #if all(x==geno_list[0] for x in geno_list):
                    outlist = [concatenated, ref_allele, alt_allele, genostring, genoPLUSdepth] # collect each of these for each line where the genotypes are the same
                    outstring = '\t'.join(outlist) 

                    listofoutlist.append(outstring) #append each line to the list started up at the top
                    #print(concatenated, outlist)
                #else:
                    #continue #skips all sites that have different genotype calls betweeen the two technical replicates
                    #print(len(geno_list))
    return listofoutlist
    return dictionary #returns the full list of genotype matches in the file

    f.close()

if numberofsamples==7:
    genos1=fileinfo(input1, 1, 2)#, "geno_list[1]", "geno_list[2]") #calls the function for the first input file    
#print(genos1)
    genos2=fileinfo(input1, 4, 5) #calls the function for the second input file
else:   
    genos1=fileinfo(input1, 0, 1) 
    genos2=fileinfo(input1, 2, 3) 
    
#print(genos2)
# #
#
conc1list = []
geno1list = []
for line1 in genos1: #goes line by line in the first input file's list of matches
        #print(line1)
        line1 = line1.strip()

        items1 = line1.split('\t')
        #print(items1)
        conc1 = items1[0]

        genotypes1 = items1[3:]
        geno1 = genotypes1[1]
        #print(geno1)

        conc1list.append(conc1)
        geno1list.append(geno1)


conc1_and_geno1 = zip(conc1list, geno1list)

dictOfWords = dict(conc1_and_geno1)
#print(dictOfWords)
#
conc2list = []
geno2list = []


for line2 in genos2: #now goes line by line in the second input file's list of matches
    if line2.startswith('chrom.pos'):
                #print(line2)
                line2 = line2.strip()

                items2 = line2.split('\t')

                header = items2
                header = '\t'.join(header)
                #finalheader = [header, "TrueorFalse","TypeofMutation", "MutantParent1", "MutantParent2","MutantSperm1","MutantSperm2"]
                #finalheader = '\t'.join(finalheader)

                print(header)
    else:
            line2 = line2.strip()

            items2 = line2.split('\t')
            #print(items2)
            conc2 = items2[0]
            ref = items2[1]
            #print(ref)
            alt = items2[2]
    
            genotypes2= items2[3:]
            #print(genotypes2)
            if numberofsamples ==7:
                geno2= genotypes2[4]
            if numberofsamples ==4:
                geno2=genotypes2[2]
            #print(geno2)
            if conc2 in dictOfWords:
                geno1 = dictOfWords[conc2]

                if geno1 != geno2:
                    genolist = genotypes2[7:]
                    if numberofsamples==7:
                        writeout = [conc2, ref, alt, genotypes2[7],genotypes2[8], genotypes2[9], genotypes2[10], genotypes2[11], genotypes2[12], genotypes2[13]]
                    if numberofsamples==4:
                        writeout = [conc2, ref, alt, genotypes2[4],genotypes2[5], genotypes2[6], genotypes2[7]]
                    #writeout = [conc2, ref, alt, genolist]
                    writeout_string = '\t'.join(writeout)
                    #print(writeout_string)
                    header = "Chrom.pos","ref", "alt", "sample1","sample2","sample3","sample4","sample5","sample6","sample7"
                    header_string = '\t'.join(header)
                    #print(header)
                    print(writeout_string)    
                    

