# Python script that scans genomic bins for differential SNP enrichment

import csv
import sys, getopt
import collections
import pandas as pd
from scipy.stats import fisher_exact


# Command Line Inputs:
def main(argv):
   diffSNPs = ''
   allSNPs = ''
   outputfile = ''
   genome = 'hg38'
   binsize = 10000
   norm="chrom"
   
   try:
      opts, args = getopt.getopt(argv, "hi:a:o:g:b:n", ["diffSNPs=","allSNPs","output=","genome=", "binsize=", "NormalizeBy="])
   except getopt.GetoptError:
      print('\nFind_Locus_Features.py\n-i <diff SNPs>\n-a <all SNPs>\n-o <outputfile>\n-g <genome>\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Binned_Enrichment Required Inputs ---\n\n\t-i | --diffSNPs\t<differential SNP inputfile as a tab-delimited text file.>\n\t-a | --allSNPs\t<all SNP inputfile as a tab-delimited text file.>\n\t-o | --output\t<outputfile>\n\n--- Optional Inputs ---\n\n\t-g | --genome\t<reference genome name>, default=hg38; available genomes: hg38\n\t-n | --NormalizeBy\t<Value to normalize the bins against; either by genome or by chromosome.>, default=chrom; available options: genome, chrom\n\t-b | --binsize\t<query bin size (bp)>, default=10000\n\n')
         sys.exit()
      elif opt in ("-i", "--diffSNPs"):
         diffSNPs = arg
      elif opt in ("-a", "--allSNPs"):
         allSNPs = arg
      elif opt in ("-o", "--output"):
         outputfile = arg
      elif opt in ("-g", "--genome"):
          genome = arg
      elif opt in ("-b", "--binsize"):
          binsize = int(arg)
      elif opt in ("-n", "--NormalizeBy"):
          norm = arg
          
   print('\n\n--------------------------------\nLooking for SNP enrichment in', binsize, "bins against", genome, "\nInputs=", diffSNPs, "and", allSNPs, "\nOutput=", outputfile, "\n--------------------------------\n\n")
   
   regions = pd.read_csv(genome+'/region_sizes.txt', sep="\t",  low_memory=False)
   diff_sites = pd.read_csv(diffSNPs, sep="\t",  low_memory=False)
   all_sites = pd.read_csv(allSNPs, sep="\t",  low_memory=False)
   
   # universal values:
   genome_len = regions["Region_Size"].sum()
   #print(totalSNP_count, diffSNP_count, bkgdSNP_count, genome_len)
   
   # output columns:
   chromV = []
   stV = []
   eV = []
   actual_diff = []
   actual_bkgd = []
   expect_diff = []
   expect_bkdg = []
   ER_diff = []
   ER_bkgd = []
   sig2N =[]


# For chrom in unique chrom list
   for chrom in regions["Region"].unique():
       
       print("Scanning", chrom, "...")
       chrom_len = int(regions.loc[regions['Region'] == chrom, 'Region_Size'].iloc[0])
       
       #chrom specific DFs for all and diff SNPs:
       chrom_DIFF = diff_sites[diff_sites["CHROM"] == chrom] #not sure if syntax correct
       chrom_ALL = all_sites[all_sites["CHROM"] == chrom] #not sure if syntax correct
       
       if chrom_DIFF.empty:
           print("No variants detected in "+chrom)
           continue
       
       if norm == "chrom":
           compVal = chrom_len
           totalSNP_count = len(chrom_ALL["CHROM"]) # total # of SNPs in the chromosome
           diffSNP_count = len(chrom_DIFF["CHROM"]) # total # of diff SNPs in the chromosome
       elif norm == "genome":
           compVal = genome_len
           # genomic expected total SNP occurance (SNPs/bp) = total SNPs / genome length'
           totalSNP_count = len(all_sites["CHROM"])
           # genomic expected differential SNP occurance (differential SNPs/bp) = differential SNPs / genome length
           diffSNP_count = len(diff_sites["CHROM"])
       
       bkgdSNP_count = (totalSNP_count-diffSNP_count) #bkgd SNPs = all NS SNPs
       print(chrom, "has", totalSNP_count, "SNPs, of which", diffSNP_count, "are in the input (", ((diffSNP_count/totalSNP_count)*100), "%)\n")
       
       start=1
       while start <= chrom_len:
           end = start+(binsize-1)
           
           if end > chrom_len:
               # if the end is larger than the chromosome length:
               tbin =chrom_len-start
               end=chrom_len
           else:
               tbin=binsize
           
           #output df cols:
           chromV.append(chrom)
           stV.append(start)
           eV.append((end))
           
           #df for one iteration
           i_diff =  chrom_DIFF[(chrom_DIFF["POS"] >= start) & (chrom_DIFF["POS"] <= end)] # not sure if correct syntax
           i_all =  chrom_ALL[(chrom_ALL["POS"] >= start) & (chrom_ALL["POS"] <= end)] # not sure if correct syntax
           
           #print(i_diff.head())
           #print(i_all.head())
           
           #### Actual SNPs ###
           aDiff=len(i_diff["CHROM"]) # num differential SNPs in bin
           aTotal=len(i_all["CHROM"]) # total SNPs in bin
           aBkgd=(aTotal-aDiff) # background SNPs = (total num SNPs in bin - num differential SNPs in bin)
           #add to output:
           actual_diff.append(aDiff)
           actual_bkgd.append(aBkgd)
           
           ### Expected SNPs ###
           # expected diff SNPs in bin (diff SNPs/genome) * bin length
           eDiff = (diffSNP_count/compVal)*tbin
           # expected background SNPs in bin (background SNPs/genome) * bin length
           eBkgd = (bkgdSNP_count/compVal)*tbin
           #add to output:
           expect_diff.append(eDiff)
           expect_bkdg.append(eBkgd)
           
           ### Enrichment Ratios ###
           # background SNP enrichment ratio = (num background SNPs in bin/expected background SNPs in bin)
           # diff SNP enrichment ratio = (num SNP enrichment in bin/expected diff SNPs in bin)
           ER_diff.append((aDiff/eDiff))
           ER_bkgd.append((aBkgd/eBkgd))
           # Signal2Noise = diff SNP enrichment ratio/background SNP enrichment ratio
           
           if aBkgd == 0:
               sig2N.append("NA")
           else:
               sig2N.append(((aDiff/eDiff)/(aBkgd/eBkgd)))
           
           start = start + binsize
       
       #write output:
   outputDF = pd.DataFrame(list(zip(chromV, stV, eV, actual_diff, actual_bkgd, expect_diff, expect_bkdg, ER_diff, ER_bkgd, sig2N)),
                      columns =["Region", "start", "end", "dSNPs", "bkgdSNPs", "expected_dSNPs", "expected_bkgdSNPs", "dEnrichRatio", "bkgdEnrichRatio", "signal2noise"])
   outputDF.to_csv(outputfile, sep="\t", chunksize=10000)
   print("Done :)")
       
# write output df...
if __name__ == "__main__":
   main(sys.argv[1:])
