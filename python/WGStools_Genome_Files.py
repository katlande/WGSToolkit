#!/usr/bin/env python
import csv
import sys, getopt
import collections
import pandas as pd

# Functions:

# Makes chromosome-specific files for an input feature
def OneFeature(inname, outname, filename, chri):
    subset = inname[inname[0] == chri]
    outname=("chr"+chri+"."+filename)
    subset.to_csv(outname,sep="\t", header=False, index=False)

# Command Line Inputs:
def main(argv):
   Name=''
   
   try:
      opts, args = getopt.getopt(argv, "hn", ["name="])
   except getopt.GetoptError:
      print('')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Script Required Inputs ---\n\n\t--name\t<genome name shorthand>\n\n')
         sys.exit()
      elif opt in ("-n", "--name"):
         Name = arg
   
   # Main code body
   exons=pd.read_csv((Name+".exons.bed"), sep="\t", header=None, low_memory=False)
   genes=pd.read_csv((Name+".genes.bed"), sep="\t", header=None, low_memory=False)
   intergenic=pd.read_csv((Name+".intergenic.bed"), sep="\t", header=None, low_memory=False)
   introns=pd.read_csv((Name+".introns.bed"), sep="\t", header=None, low_memory=False)
   masked_introns=pd.read_csv((Name+".masked_introns.bed"), sep="\t", header=None, low_memory=False)
   splice_site=pd.read_csv((Name+".splice_site.bed"), sep="\t", header=None, low_memory=False)
   tss_500=pd.read_csv((Name+".TSS_500.bed"), sep="\t", header=None, low_memory=False)
   utr=pd.read_csv((Name+".UTR.bed"), sep="\t", header=None, low_memory=False)
   
   chromList=(pd.unique(exons[0]))
   chromnamelist = [] # for chrom.txt
   Region_Sizelist = [] #remainder for region_sizes.txt
   Exon_Totallist = []
   Intron_Totallist = []
   TSS_Totallist = []
   Genes_Totallist= []
   Masked_Intron_Totallist= []
   Splice_Site_Totallist= []
   UTR_Totallist = [] 
   Intergenic_Totallist = []
   
   for chrom in chromList:
       chromnamelist.append("chr"+chrom) #chrom.txt
       print("Splitting features for chr"+chrom)
       
       OneFeature(exons, ("chr"+chrom+"."+Name+".exons.bed"), (Name+".exons.bed"), chrom)
       subset = exons[exons[0] == chrom] # for region sizes
       Exon_Totallist.append(sum(subset[2].astype(int) - subset[1].astype(int)))
       #print("Exons parsed")
       
       OneFeature(genes, ("chr"+chrom+"."+Name+".genes.bed"), (Name+".genes.bed"), chrom)
       subset = genes[genes[0] == chrom] # for region sizes
       genesize=sum(subset[2].astype(int) - subset[1].astype(int))
       Genes_Totallist.append(genesize)
       #print("Genes parsed")
       
       OneFeature(intergenic, ("chr"+chrom+"."+Name+".intergenic.bed"), (Name+".intergenic.bed"), chrom)
       subset = intergenic[intergenic[0] == chrom] # for region sizes
       intsize=sum(subset[2].astype(int) - subset[1].astype(int))
       Intergenic_Totallist.append(intsize)
       Region_Sizelist.append((int(intsize)+int(genesize)))#Region_Size = genic + intergenic
       #print("Intergenic regions parsed")
       
       OneFeature(introns, ("chr"+chrom+"."+Name+".introns.bed"), (Name+".introns.bed"), chrom)
       subset = introns[introns[0] == chrom] # for region sizes
       Intron_Totallist.append(sum(subset[2].astype(int) - subset[1].astype(int)))
       #print("Introns parsed")
       
       OneFeature(masked_introns, ("chr"+chrom+"."+Name+".masked_introns.bed"), (Name+".masked_introns.bed"), chrom)
       subset = masked_introns[masked_introns[0] == chrom] # for region sizes
       Masked_Intron_Totallist.append(sum(subset[2].astype(int) - subset[1].astype(int)))
       #print("Masked Introns parsed")
       
       OneFeature(splice_site, ("chr"+chrom+"."+Name+".splice_site.bed"), (Name+".splice_site.bed"), chrom)
       subset = splice_site[splice_site[0] == chrom] # for region sizes
       Splice_Site_Totallist.append(sum(subset[2].astype(int) - subset[1].astype(int)))
       #print("SSs parsed")
       
       OneFeature(tss_500, ("chr"+chrom+"."+Name+".TSS_500.bed"), (Name+".TSS_500.bed"), chrom)
       subset = tss_500[tss_500[0] == chrom] # for region sizes
       TSS_Totallist.append(sum(subset[2].astype(int) - subset[1].astype(int)))
       #print("TSSs parsed")
       
       OneFeature(utr, ("chr"+chrom+"."+Name+".UTR.bed"), (Name+".UTR.bed"), chrom)
       subset = utr[utr[0] == chrom] # for region sizes
       UTR_Totallist.append(sum(subset[2].astype(int) - subset[1].astype(int)))
       #print("UTRs parsed")
       
   chromdf = pd.DataFrame(chromnamelist)
   chromdf.to_csv("chrom.txt", sep="\t", header=False, index=False)
   
   # Region sizes file:
   Regiondf = pd.DataFrame(list(zip(chromnamelist, Region_Sizelist, Exon_Totallist, Intron_Totallist, TSS_Totallist, Genes_Totallist, Masked_Intron_Totallist, Splice_Site_Totallist, UTR_Totallist, Intergenic_Totallist)),
                  columns =["Region", "Region_Size", "Exon_Total", "Intron_Total", "TSS_Total", "Genes_Total", 
                  "Masked_Intron_Total", "Splice_Site_Total", "UTR_Total", "Intergenic_Total"])
   Regiondf.to_csv("region_sizes.txt", sep="\t", index=False)
   
if __name__ == "__main__":
   main(sys.argv[1:])

