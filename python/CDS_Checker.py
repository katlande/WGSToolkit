import csv
import sys, getopt
import collections
import pandas as pd
from os import walk

def check_cds(gene, genome, pos):
    file_list = next(walk(genome+"/CDS/"+gene), (None, None, []))[2]# save each transcript filename to a list
    hit=False
    i=0
    while i < len(file_list) and hit == False:
        with open(genome+"/CDS/"+gene+"/"+file_list[i]) as f:
            cds = f.readline().strip('\n')
            cds = cds.split(', ')
            cds=list(map(int, cds))
            if int(pos) in cds:
                hit=True
        i=i+1
    if hit == True:
        return(True)
    else:
        return(False)

# Command Line Inputs:
def main(argv):
   inputfile = ''
   outputfile = ''
   genome = 'hg38'
   cCol = 0
   pCol = 1
   gCol = 2
   
   try:
      opts, args = getopt.getopt(argv, "hi:o:c:p:n:g", ["input=","output=","chromcol=", "poscol=", "genecol=", "genome=","base="])
   except getopt.GetoptError:
      print('\nCall_Codons.py\n-i <input>\n-o <output>\n-c <chromcol>\n-p <poscol>\n-n <genecol>\n-g <genome>\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\nCDS_Checker scans to see if a list of genomic loci are in the CDS regions of any protein coding genes.\n\n--- Required Inputs ---\n\n')
         print('\t-i | --input\t\t<Tab-separated text file containing chromosomes and positions of loci>\n\t-o | --output\t\t<outputfile>\n\t-c | --chromcol\t\t<index of chromosome column in input file, first index=0>, default=0\n\t-p | --poscol\t\t<index of position column in input file, first index=0>, default=1\n\t-n | --genecol\t\t<index of gene column in input file, first index=0>, default=2\n\t-g | --genome\t\t<genome version to use>\tdefault=hg38; available genomes:\n\t\t\t\t- hg38\n\t\t\t\t- hg37\n\t\t\t\t- BDGP6\n\t\t\t\t- mm10\n\t\t\t\t- TAIR10\n\t\t\t\t- WBcel235\n\n')
         sys.exit()
      elif opt in ("-i", "--input"):
         inputfile = arg
      elif opt in ("-o", "--output"):
         outputfile = arg
      elif opt in ("-g", "--genome"):
          genome = arg
      elif opt in ("-c", "--chromcol"):
          cCol = int(arg)
      elif opt in ("-p", "--poscol"):
          pCol = int(arg)
      elif opt in ("-n", "--genecol"):
          gCol = int(arg)
      elif opt in ("-b", "--base"):
          basename = arg
          
   print('\n\n--------------------------------\nFinding CDS Codon Alterations in', basename, ' using ', inputfile, "\nOutput=", outputfile, 
   "\n--------------------------------\n\n")
   
   print("Collecting strand data...")
   strand_dict = {}
   with open(genome+"/"+basename+"_ProtCoding_Gene_Strands.txt") as file:
    for line in file:
       strand_dict[line.split()[0]] = [i.replace(',','') for i in line.split()[1:]]
   
   print("Reading in additional gene information...")
   broken_genes = []
   with open(genome+"/"+basename+"_Broken_Gene_Models.txt") as file:
    for line in file:
       broken_genes.append(line.split()[0])
   
   chrom_list = []
   pos_list = []
   gene_list = []
   CDS_list = []
   
   with open(inputfile) as readfile:
        rf=csv.reader(readfile, delimiter="\t")
        next(rf)
        line_count=0
        CDS_count=0
        for line in rf:
            
            line_count=line_count+1
            gene=line[gCol]
            gene=gene.replace('/', "-") #I hate you arabidopsis
            pos=line[pCol]
            chrom=line[cCol]
            
            chrom_list.append(chrom)
            pos_list.append(pos)
            gene_list.append(gene)
            
            if gene in strand_dict.keys() and not gene in broken_genes:
                check=check_cds(gene, genome, pos)
                if check==True:
                    CDS_count=CDS_count+1
                    CDS_list.append("TRUE")
                else:
                    CDS_list.append("FALSE")
            else:
                CDS_list.append("FALSE")
   
   print(line_count, "lines scanned,", CDS_count, "CDS loci found.")
   print("Writing output file...")
   df = pd.DataFrame(list(zip(chrom_list, pos_list, gene_list, CDS_list)),
                  columns =['CHROM', 'POS', 'GENE', 'CDS'])
   df.to_csv(outputfile, sep="\t")
   
if __name__ == "__main__":
   main(sys.argv[1:])