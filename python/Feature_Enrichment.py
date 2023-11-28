import csv
import sys, getopt
import collections
import pandas as pd
from scipy.stats import fisher_exact

# Functions

def findEnrichment(df, feature, regionsize, subsetsize):
    #print("total feature length:", len(df[feature]))
    #print("feature (no hits):", len(df[df[feature].isnull()]))
    #print("feature (hits):", len(df[df[feature].notna()]))
    
    act = len(df[df[feature].notna()])
    expect = (len(df[feature]) / regionsize)*subsetsize
    
    if expect == 0:
        # If a feature is not present in a chromosome/region:
        output = [str(df["CHROM"].unique()[0]), feature, "0", "0", "0", "1"]
    else:
        other_SNPs = len(df[df[feature].isnull()])
        # p-value with fisher's exact test:
        # asking: are frequencies pf SNPs in the feature significantly different from frequencies of SNPs not in the feature?
        #                   In Feature      Not In Feature
        # variable (SNP)      #bp [a]           #bp [b]
        # not variable        #bp [c]           #bp [d]
        
        # a = act                                   = SNPs annotated to feature
        # b = other_SNPs                            = SNPs not annotated to feature
        # c = (subsetsize-act)                      = all non-variable sites in feature
        # d = ((regionsize-subsetsize)-other_SNPs)  = all non-variable sites outside of feature
        oddsr, p = fisher_exact(table=pd.DataFrame({'InFeature':[act, (subsetsize-act)], 'NotInFeature':[other_SNPs, ((regionsize-subsetsize)-other_SNPs)]}, index=pd.Index(['Var', 'nonVar'])).to_numpy(), alternative='two-sided')
        output = [str(df["CHROM"].unique()[0]), feature, act, expect, (act/expect), p]
        
    return(output)

#df = tempdf (df of one region)
#feature = name of feature in input file; one of: ("TSS", "GENE", "EXON", "UNMASKED_INTRON", "MASKED_INTRON", "INTERGENIC", "UTR", "SPLICE_SITE")
#regionsize = chrom length (bp)
#subsetsize = feature length (bp)


# Command Line Inputs:
def main(argv):
   inputfile = ''
   outputfile = ''
   genome = 'hg38'
   
   try:
      opts, args = getopt.getopt(argv, "hi:o:g", ["input=","output=","genome="])
   except getopt.GetoptError:
      print('\nFind_Locus_Features.py\n-i <inputfile>\n-o <outputfile>\n-g <genome>\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Feature_Enrichment Required Inputs ---\n\n\t-i | --input\t<inputfile as a tab-delimited text file. Script is designed to take the output from Find_Locus_Features.py>\n\t-o | --output\t<outputfile>\n\n--- Optional Inputs ---\n\n\t-g | --genome\t<reference genome name>\n\n')
         sys.exit()
      elif opt in ("-i", "--input"):
         inputfile = arg
      elif opt in ("-o", "--output"):
         outputfile = arg
      elif opt in ("-g", "--genome"):
          genome = arg
          
   print('\n\n--------------------------------\nLooking for feature enrichment in', inputfile, "against", genome, "\nOutput=", outputfile, "\n--------------------------------\n\n")
   
   # genome data paths must contain:
   # guide.txt = name of base file for each feature of: 'gene', 'exon', 'intron_masked', 'TSS', 'UTR', 'splice_site'
   # chrom.txt = name of each chromosome in the genome
   # for each base file in "guide.txt," parsed chromosome files in the format: chr1.basefile.bed, chr2.basefile.bed, chr3.basefile.bed...
   
   regions = pd.read_csv(genome+'/region_sizes.txt', sep="\t")
   sampleData = pd.read_csv(inputfile, sep="\t")
   
   # Output columns:    Region  Feature     Actual  Expected    Enrichment Ratio    p-val
   
   # Region = each chr, total
   regionCol = []
   # Feature = Feature name (e.g., exon, intron, etc)
   featureCol = []
   # Actual = # annotated SNPs in region
   actualCol = []
   # Expected = total SNPs / total size (bp) * feature size (bp)
   expectedCol = []
   # Enrichment ratio = Actual/Expected
   erCol = []
   # p val from Chi or Fisher's Exact test
   pCol = []
   
   # for each input region:
   for chrom in regions["Region"].unique():
       
       tempdf = sampleData[sampleData["CHROM"] == chrom]
       sizedf = regions[regions["Region"] == chrom]
       region_size = int(sizedf["Region_Size"].unique())
       
       if tempdf.empty:
           print("No variants detected in "+chrom)
       else:
           print("Querying enrichment in", chrom, "...")
           # for each chromosome, get enrichment for each feature: 
           temp_genes = findEnrichment(tempdf, "GENE", region_size, int(sizedf["Genes_Total"].unique()))
           temp_tss = findEnrichment(tempdf, "TSS", region_size, int(sizedf["TSS_Total"].unique()))
           temp_exon = findEnrichment(tempdf, "EXON", region_size, int(sizedf["Exon_Total"].unique()))
           temp_uint = findEnrichment(tempdf, "UNMASKED_INTRON", region_size, int(sizedf["Intron_Total"].unique()))
           temp_mint = findEnrichment(tempdf, "MASKED_INTRON", region_size, int(sizedf["Masked_Intron_Total"].unique()))
           temp_intergenic = findEnrichment(tempdf, "INTERGENIC", region_size, int(sizedf["Intergenic_Total"].unique()))
           temp_ss = findEnrichment(tempdf, "SPLICE_SITE", region_size, int(sizedf["Splice_Site_Total"].unique()))
           temp_utr = findEnrichment(tempdf, "UTR", region_size, int(sizedf["UTR_Total"].unique()))
       
           regionCol.extend([temp_genes[0], temp_tss[0], temp_exon[0], temp_uint[0], temp_mint[0], temp_intergenic[0], temp_ss[0], temp_utr[0]])
           # Feature = Feature name (e.g., exon, intron, etc)
           featureCol.extend([temp_genes[1], temp_tss[1], temp_exon[1], temp_uint[1], temp_mint[1], temp_intergenic[1], temp_ss[1], temp_utr[1]])
           # Actual = # annotated SNPs in region
           actualCol.extend([temp_genes[2], temp_tss[2], temp_exon[2], temp_uint[2], temp_mint[2], temp_intergenic[2], temp_ss[2], temp_utr[2]])
           # Expected = ((feature size/region size) * total SNPs in region #)
           expectedCol.extend([temp_genes[3], temp_tss[3], temp_exon[3], temp_uint[3], temp_mint[3], temp_intergenic[3], temp_ss[3], temp_utr[3]])
           # Enrichment ratio = Actual/Expected
           erCol.extend([temp_genes[4], temp_tss[4], temp_exon[4], temp_uint[4], temp_mint[4], temp_intergenic[4], temp_ss[4], temp_utr[4]])
           # p val from Chi or Fisher's Exact test
           pCol.extend([temp_genes[5], temp_tss[5], temp_exon[5], temp_uint[5], temp_mint[5], temp_intergenic[5], temp_ss[5], temp_utr[5]])
   
   # region size = entire genome (bp)
   # subset size = feature size across all chrom
   print("Querying genome-wide enrichment...")
   full_genome_gene = findEnrichment(sampleData, "GENE", regions['Region_Size'].sum(), regions['Genes_Total'].sum())
   full_genome_tss = findEnrichment(sampleData, "TSS", regions['Region_Size'].sum(), regions['TSS_Total'].sum())
   full_genome_exon = findEnrichment(sampleData, "EXON", regions['Region_Size'].sum(), regions['Exon_Total'].sum())
   full_genome_uint = findEnrichment(sampleData, "UNMASKED_INTRON", regions['Region_Size'].sum(), regions['Intron_Total'].sum())
   full_genome_mint = findEnrichment(sampleData, "MASKED_INTRON", regions['Region_Size'].sum(), regions['Masked_Intron_Total'].sum())
   full_genome_intergenic = findEnrichment(sampleData, "INTERGENIC", regions['Region_Size'].sum(), regions['Intergenic_Total'].sum())
   full_genome_ss = findEnrichment(sampleData, "SPLICE_SITE", regions['Region_Size'].sum(), regions['Splice_Site_Total'].sum())
   full_genome_utr = findEnrichment(sampleData, "UTR", regions['Region_Size'].sum(), regions['UTR_Total'].sum())
   
   print("Combining enrichment data...")
   regionCol.extend(["Genomic", "Genomic", "Genomic", "Genomic", "Genomic", "Genomic", "Genomic", "Genomic"])
   featureCol.extend([full_genome_gene[1], full_genome_tss[1], full_genome_exon[1], full_genome_uint[1], full_genome_mint[1], full_genome_intergenic[1], full_genome_ss[1], full_genome_utr[1]])
   actualCol.extend([full_genome_gene[2], full_genome_tss[2], full_genome_exon[2], full_genome_uint[2], full_genome_mint[2], full_genome_intergenic[2], full_genome_ss[2], full_genome_utr[2]])
   expectedCol.extend([full_genome_gene[3], full_genome_tss[3], full_genome_exon[3], full_genome_uint[3], full_genome_mint[3], full_genome_intergenic[3], full_genome_ss[3], full_genome_utr[3]])
   erCol.extend([full_genome_gene[4], full_genome_tss[4], full_genome_exon[4], full_genome_uint[4], full_genome_mint[4], full_genome_intergenic[4], full_genome_ss[4], full_genome_utr[4]])
   pCol.extend([full_genome_gene[5], full_genome_tss[5], full_genome_exon[5], full_genome_uint[5], full_genome_mint[5], full_genome_intergenic[5], full_genome_ss[5], full_genome_utr[5]])
   
   print("Writing output...")
   outputDF = pd.DataFrame(list(zip(regionCol, featureCol, actualCol, expectedCol, erCol, pCol)),
                  columns =["Region", "ALL_REGIONS", "Actual_SNPs", "Expected_SNPs", "EnrichmentRatio", "p"])
   
   #print(outputDF)
   outputDF.to_csv(outputfile, sep="\t", chunksize=10000)
   print("Done :)")
   
if __name__ == "__main__":
   main(sys.argv[1:])
   
   
   
   
   
   
   
   
   
   
   