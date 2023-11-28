import csv
import sys, getopt
import collections
import pandas as pd
from sys import exit

# Functions:

def readFeatures(filename): # internal function that saves input files as nested lists
  output = []
  with open(filename, 'r') as f:
      b = csv.reader(f, delimiter="\t")
      for line in b:
          #print(line)
          output.append(line)
      return(output)

def makeFeatureDict(path, filename): # makes a chromosome-keyed dict for each input feature file
    with open(path+"/chrom.txt") as f:
        TestDict = {}
        for line in f:
            line=line.strip("\n")
            #print(line)
            chrname = ''.join(line)
            TestDict[chrname] = readFeatures(path+"/"+chrname+"."+filename)
            #output = readFeatures(path+"/"+chrname+"."+filename)
            #print(path+"/"+chrname+"."+filename)
        return(TestDict)

def RollThru(input, query_pos, featurelist, mode): # loops through one entry in a feature file dictionary
    if mode == "gene_UTR": #genes and UTRs
        i = 5
    elif mode == "internal": #masked intron, exon
        i = 4
    elif mode == "weirdos": #splice site, TSS
        i = 6
    else:
        print("internal error!")
        i = 1000 #breaks the script
    
    
    # NOTE: find version for ordered genomes on gpfs, remake ALL genomes that come disordered. This script takes WAY too long otherwise...
    
    
    added=0
    for g in input:
        start = int(g[1])
        end = int(g[2])
        
        if query_pos >= start and query_pos <= end:
            featurelist.append(g[i])
            added=1
            break #if we find a gene, skip to the next position
    
    if added == 0:
        featurelist.append('NA')
    
    return(featurelist)

# input =  gene_Library[line_chrom]
# query_pos = line_pos
# featurelist =  gene


# Command Line Inputs:
def main(argv):
   inputfile = ''
   outputfile = ''
   genome = 'genome'
   chromLoc = 1
   posLoc = 2
   
   try:
      opts, args = getopt.getopt(argv, "hi:o:g:c:p", ["input=","output=","genome=","ChromCol=","PosCol="])
   except getopt.GetoptError:
      print('\nFind_Locus_Features.py\n-i <inputfile>\n-o <outputfile>\n-g <genome>\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Find_Locus_Features Required Inputs ---\n\n\t-i | --input\t<inputfile as a tab-delimited text file>\n\t-o | --output\t<outputfile>\n\n--- Optional Inputs ---\n\n\t-c | --ChromCol\t\t<CHROM column index in input data, first column=0>, default=1\n\t-p | --PosCol\t\t<POS column index in input data, first column=0>, default=2. *NOTE: This script seems to break if -p is used instead of --PosCol.\n\t-g | --genome\t\t<genome version to use>\tdefault=hg38.')
         sys.exit()
         sys.exit()
      elif opt in ("-i", "--input"):
         inputfile = arg
      elif opt in ("-o", "--output"):
         outputfile = arg
      elif opt in ("-g", "--genome"):
          genomePath = arg
      elif opt in ("-c", "--ChromCol"):
          chromLoc = int(arg)
      elif opt in ("-p", "--PosCol"):
          posLoc = int(arg)
          
   print('\n\n--------------------------------\nFinding locus features in', genomePath, '\nInput=', inputfile, "\nOutput=", outputfile, "\n--------------------------------\n\n")
   
   # Get genome source data:
   genomeDict = {genome: genomePath}
   # genome data paths must contain:
   # guide.txt = name of base file for each feature of: 'gene', 'exon', 'intron_masked', 'TSS', 'UTR', 'splice_site'
   # chrom.txt = name of each chromosome in the genome
   # for each base file in "guide.txt," parsed chromosome files in the format: chr1.basefile.bed, chr2.basefile.bed, chr3.basefile.bed...
   
   genomData = genomeDict[genome]
   with open(genomData+"/guide.txt") as tsv:
       FeatureFiles = dict(i.rstrip().split(None, 1) for i in tsv)
   # Base file names for individual features are now saved in a dictionary with the following keys:
   # 'gene', 'exon', 'intron_true', 'intron_masked', 'intergenic', 'TSS', 'UTR', 'splice_site'
   
   
   # So now, for each line of inputfile, we follow the logic tree below:
   #
   #                                             | true --> mark as TSS
   #               --> is in TSS of current gene?
   #              |                              | false
   #              |                                                    | true --> mark as masked intron|                      | false                                                              | false
   #          (parallel)            | true --> contains masked intron? |                               |--> is in splice_site?                                              |--> is in splice_site?
   #              |                 |                                  | false                         |                      | true --> mark as splice site, donor/acceptor                       | true --> UTR, pull out 3'/5'
   #              |                 |
   #          | true --> is in exon?                      | true --> UTR, pull out 3'/5'
   #          |                     | 
   #          |                     | false --> is in UTR?                              | false --> mark as intron (don't need to open file)
   #          |                                           | false --> is in splice_site?
   #  in gene?                                                                          | true --> mark as splice site, pull out donor/acceptor
   #          |
   #          |                     | true --> mark as non-genic TSS
   #          | false --> is in TSS? 
   #                                | false --> mark as intergenic (don't need to open file)
   #
   
   print("Making feature libraries...")
   # Make dictionaries for each input file, files can be called by chromosome!
   # 'gene', 'exon', 'intron_masked', 'TSS', 'UTR', 'splice_site'
   gene_Library = makeFeatureDict(genomData, FeatureFiles["gene"])
   exon_Library = makeFeatureDict(genomData, FeatureFiles["exon"])
   intron_masked_Library = makeFeatureDict(genomData, FeatureFiles["intron_masked"])
   TSS_Library = makeFeatureDict(genomData, FeatureFiles["TSS"])
   UTR_Library = makeFeatureDict(genomData, FeatureFiles["UTR"])
   splice_site_Library = makeFeatureDict(genomData, FeatureFiles["splice_site"])
   
   print("Scanning input for features...")
   with open(inputfile, 'r') as sites:
       
       # Columns for output file:
       varChrom = []
       varPos = []
       tss = []
       gene = []
       exon = []
       intron = []
       masked_intron = []
       intergenic = []
       utr = []
       splice_site = []
       
       s = csv.reader(sites, delimiter="\t")
       next(s)
       for line in s:
           # extract chrom and pos for single locus in input file:
           line_chrom = ''.join(line[chromLoc])
           
           if line_chrom in gene_Library.keys():#exclude weirdos
               line_pos = int(''.join(line[posLoc]))
               varChrom.append(line_chrom)
               varPos.append(line_pos)
               # Logic gate 1: is it a gene?
               # open the gene library for chromosome of interest and being scanning
               gene = RollThru(gene_Library[line_chrom], line_pos, gene, mode="gene_UTR")
               if gene[(len(gene)-1)] == 'NA': # if not gene
                   #print("FAIL")
                   # Append NA to genic feature lists:
                   intron.append('NA')
                   masked_intron.append('NA')
                   utr.append('NA')
                   splice_site.append('NA')
                   exon.append('NA')
               
                   #check if in TSS:
                   tss = RollThru(TSS_Library[line_chrom], line_pos, tss, mode="weirdos")
                   if tss[(len(tss)-1)] == 'NA': # if not TSS, sit is intergenic
                       intergenic.append('intergenic')
                   else:
                       intergenic.append('NA') #if is tss, not intergenic
               else:# if in gene, 
                   intergenic.append('NA') #append NA to intergenic list, then continue logic tree:
                   tss = RollThru(TSS_Library[line_chrom], line_pos, tss, mode="weirdos") #check if in TSS
               
                   # check if exon:
                   exon = RollThru(exon_Library[line_chrom], line_pos, exon, mode="internal")
                   if exon[(len(exon)-1)] == 'NA': # if not exon
                        masked_intron.append('NA') #can't be a masked intron
                        utr = RollThru(UTR_Library[line_chrom], line_pos, utr, mode="gene_UTR") # is it a utr?
                        if utr[(len(utr)-1)] == 'NA': # if not utr, check if it's a splice site:
                            splice_site = RollThru(UTR_Library[line_chrom], line_pos, splice_site, mode="weirdos")
                        
                            # if it's not a splice site or a utr, mark as intron:
                            if splice_site[(len(splice_site)-1)] == 'NA':
                                intron.append("intron")
                        else: #if it is in a UTR, mark splice site and intron as NA
                            splice_site.append('NA')
                            intron.append('NA')
                    
                   else: #if exon
                       #utr.append('NA') # edit mar 30 2023
                       intron.append('NA') #mark as not an intron
                       #check for masked introns:
                       masked_intron = RollThru(intron_masked_Library[line_chrom], line_pos, masked_intron, mode="internal") 
                       #check for splice sites:
                       splice_site = RollThru(UTR_Library[line_chrom], line_pos, splice_site, mode="weirdos")
                       utr = RollThru(UTR_Library[line_chrom], line_pos, utr, mode="gene_UTR") # is it a utr?
           
               #print('CHROM', 'POS', 'TSS', 'GENE', 'EXON', 'UNMASKED_INTRON', 'MASKED_INTRON', "INTERGENIC", "UTR", "SPLICE_SITE")
               #print(len(varChrom), len(varPos), len(tss), len(gene), len(exon), len(intron), len(masked_intron), len(intergenic), len(utr), len(splice_site))
           
               if not len(varChrom) == len(intron):
                   print("Problem at", varChrom[(len(varChrom)-1)], ":", varPos[(len(varPos)-1)])
                   intron.append('NA')
           #print(len(varChrom))
       #print('CHROM', 'POS', 'TSS', 'GENE', 'EXON', 'UNMASKED_INTRON', 'MASKED_INTRON', "INTERGENIC", "UTR", "SPLICE_SITE")
       #print(len(varChrom), len(varPos), len(tss), len(gene), len(exon), len(intron), len(masked_intron), len(intergenic), len(utr), len(splice_site))
       print("Writing output file with", len(varChrom),"lines...")
       df = pd.DataFrame(list(zip(varChrom, varPos, tss, gene, exon, intron, masked_intron, intergenic, utr, splice_site)),
                      columns =['CHROM', 'POS', 'TSS', 'GENE', 'EXON', 'UNMASKED_INTRON', 'MASKED_INTRON', "INTERGENIC", "UTR", "SPLICE_SITE"])
       df.to_csv(outputfile, sep="\t", chunksize=10000)
       print("Done :)")
   
if __name__ == "__main__":
   main(sys.argv[1:])














