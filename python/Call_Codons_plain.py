import csv
import sys, getopt
import collections
import pandas as pd
from os import walk
from Bio import SeqIO
from pyfaidx import Fasta

# RUNS IN CONDA ENV: WGS_Extended, *NOT* base!!

# Functions:

def translate(seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

def compliment(char):
    if char == "A":
        return("T")
    elif char == "T":
        return("A")
    elif char == "G":
        return("C")
    elif char == "C":
        return("G")
    else:
        print("Unrecognized base in codon!")

def test_readingframe(temp_order, strand, line_pos):
    if strand=="+":
        cds_pos = int((temp_order.index(line_pos)))
        if((cds_pos+1) % 3 == 0):
            return("POS 3")
        elif((cds_pos) % 3 == 0):
            return("POS 1")
        elif((cds_pos+2) % 3 == 0):
            return("POS 2")
    else:
        cds_pos_neg = (len(temp_order)-1)-int((temp_order.index(line_pos)))
        if((cds_pos_neg) % 3 == 0):
            return("POS 1")
        elif((cds_pos_neg+2) % 3 == 0):
            return("POS 2")
        elif((cds_pos_neg+1) % 3 == 0):
            return("POS 3")

def most_frequent(List):
    counter = 0
    num = List[0]
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
    return num

# Command Line Inputs:
def main(argv):
   inputfile = ''
   output = ''
   geneList = ''
   genome = 'hg38'
   parent = ''
   child = ''
   poscol = ''
   
   try:
      opts, args = getopt.getopt(argv, "ht:o:n:y:g:p:c:l:r:a:x", ["txt=","output=","genecol=","genome=","chromcol=","parent=","child=","poscol=","--refcol=","base="])
   except getopt.GetoptError:
      print('\nArgument Error\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Call_Codons Required Inputs ---\n\n\t-t | --txt\t\t<tab-delimited txt file containing SNPs of interest>\n\t-o | --output\t\t<outputfile>\n\t-n | --genecol\t\t<column index of gene names>\n\t-p | --parent\t\t<index of control genotype column in txt>\n\t-c | --child\t\t<index of treatment genotype column in txt>\n\t-l | --poscol\t\t<index of the position in the input txt file>\n\t-g | --genome\t\t<genome version to use>\tdefault=hg38;  available genomes:\n\t\t\t\t- hg38\n\t\t\t\t- hg37\n\t\t\t\t- BDGP6\n\t\t\t\t- mm10\n\t\t\t\t- TAIR10\n\t\t\t\t- WBcel235\n\n')
         sys.exit()
      elif opt in ("-t", "--txt"):
         inputfile = arg
      elif opt in ("-o", "--output"):
         outputfile = arg
      elif opt in ("-n", "--genecol"):
         genes = int(arg)
      elif opt in ("-y", "--chromcol"):
         chrom = int(arg)
      elif opt in ("-g", "--genome"):
          genome = arg
      elif opt in ("-p", "--parent"):
          pCol = int(arg)
      elif opt in ("-c", "--child"):
          cCol = int(arg)
      elif opt in ("-l", "--poscol"):
          loc = int(arg)
      elif opt in ("-r", "--refcol"):
          refcol = int(arg)
      elif opt in ("-x", "--base"):
          basename = arg
          
   print('\n\n--------------------------------\nFinding CDS Codon Alterations in', genome, ' using ', inputfile, "\nOutput=", outputfile, 
   "child = ", cCol, "parent = ", pCol, "\n--------------------------------\n\n")
   
   grantham_dict = {('A', 'A'): '0', ('A', 'C'): '195', ('A', 'D'): '126', ('A', 'E'): '107', ('A', 'F'): '113', ('A', 'G'): '60', ('A', 'H'): '86', 
   ('A', 'I'): '94', ('A', 'K'): '106', ('A', 'L'): '96', ('A', 'M'): '84', ('A', 'N'): '111', ('A', 'P'): '27', ('A', 'Q'): '91', ('A', 'R'): '112', 
   ('A', 'S'): '99', ('A', 'T'): '58', ('A', 'V'): '64', ('A', 'W'): '148', ('A', 'Y'): '112', ('C', 'A'): '195', ('C', 'C'): '0', ('C', 'D'): '154', 
   ('C', 'E'): '170', ('C', 'F'): '205', ('C', 'G'): '159', ('C', 'H'): '174', ('C', 'I'): '198', ('C', 'K'): '202', ('C', 'L'): '198', ('C', 'M'): '196', 
   ('C', 'N'): '139', ('C', 'P'): '169', ('C', 'Q'): '154', ('C', 'R'): '180', ('C', 'S'): '112', ('C', 'T'): '149', ('C', 'V'): '192', ('C', 'W'): '215', 
   ('C', 'Y'): '194', ('D', 'A'): '126', ('D', 'C'): '154', ('D', 'D'): '0', ('D', 'E'): '45', ('D', 'F'): '177', ('D', 'G'): '94', ('D', 'H'): '81', 
   ('D', 'I'): '168', ('D', 'K'): '101', ('D', 'L'): '172', ('D', 'M'): '160', ('D', 'N'): '23', ('D', 'P'): '108', ('D', 'Q'): '61', ('D', 'R'): '96', 
   ('D', 'S'): '65', ('D', 'T'): '85', ('D', 'V'): '152', ('D', 'W'): '181', ('D', 'Y'): '160', ('E', 'A'): '107', ('E', 'C'): '170', ('E', 'D'): '45', 
   ('E', 'E'): '0', ('E', 'F'): '140', ('E', 'G'): '98', ('E', 'H'): '40', ('E', 'I'): '134', ('E', 'K'): '56', ('E', 'L'): '138', ('E', 'M'): '126', 
   ('E', 'N'): '42', ('E', 'P'): '93', ('E', 'Q'): '29', ('E', 'R'): '54', ('E', 'S'): '80', ('E', 'T'): '65', ('E', 'V'): '121', ('E', 'W'): '152', 
   ('E', 'Y'): '122', ('F', 'A'): '113', ('F', 'C'): '205', ('F', 'D'): '177', ('F', 'E'): '140', ('F', 'F'): '0', ('F', 'G'): '153', ('F', 'H'): '100', 
   ('F', 'I'): '21', ('F', 'K'): '102', ('F', 'L'): '22', ('F', 'M'): '28', ('F', 'N'): '158', ('F', 'P'): '114', ('F', 'Q'): '116', ('F', 'R'): '97', 
   ('F', 'S'): '155', ('F', 'T'): '103', ('F', 'V'): '50', ('F', 'W'): '40', ('F', 'Y'): '22', ('G', 'A'): '60', ('G', 'C'): '159', ('G', 'D'): '94', 
   ('G', 'E'): '98', ('G', 'F'): '153', ('G', 'G'): '0', ('G', 'H'): '98', ('G', 'I'): '135', ('G', 'K'): '127', ('G', 'L'): '138', ('G', 'M'): '127', 
   ('G', 'N'): '80', ('G', 'P'): '42', ('G', 'Q'): '87', ('G', 'R'): '125', ('G', 'S'): '56', ('G', 'T'): '59', ('G', 'V'): '109', ('G', 'W'): '184', 
   ('G', 'Y'): '147', ('H', 'A'): '86', ('H', 'C'): '174', ('H', 'D'): '81', ('H', 'E'): '40', ('H', 'F'): '100', ('H', 'G'): '98', ('H', 'H'): '0', 
   ('H', 'I'): '94', ('H', 'K'): '32', ('H', 'L'): '99', ('H', 'M'): '87', ('H', 'N'): '68', ('H', 'P'): '77', ('H', 'Q'): '24', ('H', 'R'): '29', 
   ('H', 'S'): '89', ('H', 'T'): '47', ('H', 'V'): '84', ('H', 'W'): '115', ('H', 'Y'): '83', ('I', 'A'): '94', ('I', 'C'): '198', ('I', 'D'): '168', 
   ('I', 'E'): '134', ('I', 'F'): '21', ('I', 'G'): '135', ('I', 'H'): '94', ('I', 'I'): '0', ('I', 'K'): '102', ('I', 'L'): '5', ('I', 'M'): '10', 
   ('I', 'N'): '149', ('I', 'P'): '95', ('I', 'Q'): '109', ('I', 'R'): '97', ('I', 'S'): '142', ('I', 'T'): '89', ('I', 'V'): '29', ('I', 'W'): '61', 
   ('I', 'Y'): '33', ('K', 'A'): '106', ('K', 'C'): '202', ('K', 'D'): '101', ('K', 'E'): '56', ('K', 'F'): '102', ('K', 'G'): '127', ('K', 'H'): '32', 
   ('K', 'I'): '102', ('K', 'K'): '0', ('K', 'L'): '107', ('K', 'M'): '95', ('K', 'N'): '94', ('K', 'P'): '103', ('K', 'Q'): '53', ('K', 'R'): '26', 
   ('K', 'S'): '121', ('K', 'T'): '78', ('K', 'V'): '97', ('K', 'W'): '110', ('K', 'Y'): '85', ('L', 'A'): '96', ('L', 'C'): '198', ('L', 'D'): '172', 
   ('L', 'E'): '138', ('L', 'F'): '22', ('L', 'G'): '138', ('L', 'H'): '99', ('L', 'I'): '5', ('L', 'K'): '107', ('L', 'L'): '0', ('L', 'M'): '15', 
   ('L', 'N'): '153', ('L', 'P'): '98', ('L', 'Q'): '113', ('L', 'R'): '102', ('L', 'S'): '145', ('L', 'T'): '92', ('L', 'V'): '32', ('L', 'W'): '61', 
   ('L', 'Y'): '36', ('M', 'A'): '84', ('M', 'C'): '196', ('M', 'D'): '160', ('M', 'E'): '126', ('M', 'F'): '28', ('M', 'G'): '127', ('M', 'H'): '87', 
   ('M', 'I'): '10', ('M', 'K'): '95', ('M', 'L'): '15', ('M', 'M'): '0', ('M', 'N'): '142', ('M', 'P'): '87', ('M', 'Q'): '101', ('M', 'R'): '91', 
   ('M', 'S'): '135', ('M', 'T'): '81', ('M', 'V'): '21', ('M', 'W'): '67', ('M', 'Y'): '36', ('N', 'A'): '111', ('N', 'C'): '139', ('N', 'D'): '23', 
   ('N', 'E'): '42', ('N', 'F'): '158', ('N', 'G'): '80', ('N', 'H'): '68', ('N', 'I'): '149', ('N', 'K'): '94', ('N', 'L'): '153', ('N', 'M'): '142', 
   ('N', 'N'): '0', ('N', 'P'): '91', ('N', 'Q'): '46', ('N', 'R'): '86', ('N', 'S'): '46', ('N', 'T'): '65', ('N', 'V'): '133', ('N', 'W'): '174', 
   ('N', 'Y'): '143', ('P', 'A'): '27', ('P', 'C'): '169', ('P', 'D'): '108', ('P', 'E'): '93', ('P', 'F'): '114', ('P', 'G'): '42', ('P', 'H'): '77', 
   ('P', 'I'): '95', ('P', 'K'): '103', ('P', 'L'): '98', ('P', 'M'): '87', ('P', 'N'): '91', ('P', 'P'): '0', ('P', 'Q'): '76', ('P', 'R'): '103', 
   ('P', 'S'): '74', ('P', 'T'): '38', ('P', 'V'): '68', ('P', 'W'): '147', ('P', 'Y'): '110', ('Q', 'A'): '91', ('Q', 'C'): '154', ('Q', 'D'): '61', 
   ('Q', 'E'): '29', ('Q', 'F'): '116', ('Q', 'G'): '87', ('Q', 'H'): '24', ('Q', 'I'): '109', ('Q', 'K'): '53', ('Q', 'L'): '113', ('Q', 'M'): '101', 
   ('Q', 'N'): '46', ('Q', 'P'): '76', ('Q', 'Q'): '0', ('Q', 'R'): '43', ('Q', 'S'): '68', ('Q', 'T'): '42', ('Q', 'V'): '96', ('Q', 'W'): '130', 
   ('Q', 'Y'): '99', ('R', 'A'): '112', ('R', 'C'): '180', ('R', 'D'): '96', ('R', 'E'): '54', ('R', 'F'): '97', ('R', 'G'): '125', ('R', 'H'): '29', 
   ('R', 'I'): '97', ('R', 'K'): '26', ('R', 'L'): '102', ('R', 'M'): '91', ('R', 'N'): '86', ('R', 'P'): '103', ('R', 'Q'): '43', ('R', 'R'): '0', 
   ('R', 'S'): '110', ('R', 'T'): '71', ('R', 'V'): '96', ('R', 'W'): '101', ('R', 'Y'): '77', ('S', 'A'): '99', ('S', 'C'): '112', ('S', 'D'): '65', 
   ('S', 'E'): '80', ('S', 'F'): '155', ('S', 'G'): '56', ('S', 'H'): '89', ('S', 'I'): '142', ('S', 'K'): '121', ('S', 'L'): '145', ('S', 'M'): '135', 
   ('S', 'N'): '46', ('S', 'P'): '74', ('S', 'Q'): '68', ('S', 'R'): '110', ('S', 'S'): '0', ('S', 'T'): '58', ('S', 'V'): '124', ('S', 'W'): '177', 
   ('S', 'Y'): '144', ('T', 'A'): '58', ('T', 'C'): '149', ('T', 'D'): '85', ('T', 'E'): '65', ('T', 'F'): '103', ('T', 'G'): '59', ('T', 'H'): '47', 
   ('T', 'I'): '89', ('T', 'K'): '78', ('T', 'L'): '92', ('T', 'M'): '81', ('T', 'N'): '65', ('T', 'P'): '38', ('T', 'Q'): '42', ('T', 'R'): '71', 
   ('T', 'S'): '58', ('T', 'T'): '0', ('T', 'V'): '69', ('T', 'W'): '128', ('T', 'Y'): '92', ('V', 'A'): '64', ('V', 'C'): '192', ('V', 'D'): '152', 
   ('V', 'E'): '121', ('V', 'F'): '50', ('V', 'G'): '109', ('V', 'H'): '84', ('V', 'I'): '29', ('V', 'K'): '97', ('V', 'L'): '32', ('V', 'M'): '21', 
   ('V', 'N'): '133', ('V', 'P'): '68', ('V', 'Q'): '96', ('V', 'R'): '96', ('V', 'S'): '124', ('V', 'T'): '69', ('V', 'V'): '0', ('V', 'W'): '88', 
   ('V', 'Y'): '55', ('W', 'A'): '148', ('W', 'C'): '215', ('W', 'D'): '181', ('W', 'E'): '152', ('W', 'F'): '40', ('W', 'G'): '184', ('W', 'H'): '115', 
   ('W', 'I'): '61', ('W', 'K'): '110', ('W', 'L'): '61', ('W', 'M'): '67', ('W', 'N'): '174', ('W', 'P'): '147', ('W', 'Q'): '130', ('W', 'R'): '101', 
   ('W', 'S'): '177', ('W', 'T'): '128', ('W', 'V'): '88', ('W', 'W'): '0', ('W', 'Y'): '37', ('Y', 'A'): '112', ('Y', 'C'): '194', ('Y', 'D'): '160', 
   ('Y', 'E'): '122', ('Y', 'F'): '22', ('Y', 'G'): '147', ('Y', 'H'): '83', ('Y', 'I'): '33', ('Y', 'K'): '85', ('Y', 'L'): '36', ('Y', 'M'): '36', 
   ('Y', 'N'): '143', ('Y', 'P'): '110', ('Y', 'Q'): '99', ('Y', 'R'): '77', ('Y', 'S'): '144', ('Y', 'T'): '92', ('Y', 'V'): '55', ('Y', 'W'): '37', 
   ('Y', 'Y'): '0'}
   
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
       
   print("Reading reference genome...")
   gFile = genome+"/"+basename+".fa"
   #gChr = [i for i in SeqIO.parse(gFile, 'fasta')]
   gChr = Fasta(gFile)
   
   pos_gene_list = []
   chrom_list = []
   pos_list = []
   child_geno_list = []
   parent_geno_list = []
   ref_codon_list = []
   alt_codon_list = []
   ref_peptide_list = []
   alt_peptide_list = []
   grantham_dist_list = []
   transcript_list = []
   
   print("Checking input loci...")
   genes_checked=[]
   with open(inputfile, 'r') as sites:
       s = csv.reader(sites, delimiter="\t")
       next(s)
       for line in s:
           line_pos = int(''.join(line[loc]))
           line_gene = ''.join(line[genes])
           line_gene=line_gene.replace('/', "-")
           geneChrom = ''.join(line[chrom])
           ctl_geno = ''.join(line[pCol])
           trt_geno = ''.join(line[cCol])
           
           pos_gene_list.append(line_gene)
           chrom_list.append(geneChrom)
           pos_list.append(line_pos)
           parent_geno_list.append(ctl_geno)
           child_geno_list.append(trt_geno)
           
           if not geneChrom in gChr.keys():
               geneChrom = geneChrom.replace('chr', '')
               geneChrom = geneChrom.replace('Chr', '')
               geneChrom = geneChrom.replace('CHR', '')
               if not geneChrom in gChr.keys():
                   if "chr"+geneChrom in gChr.keys():
                       geneChrom="chr"+geneChrom
                   elif "Chr"+geneChrom in gChr.keys():
                       geneChrom="Chr"+geneChrom
                   elif "CHR"+geneChrom in gChr.keys():
                       geneChrom="CHR"+geneChrom
                   else:
                       print("Error! chromosome name not recognized. Try reformating input file.")
           #print(line_gene)
           #print(strand_dict.keys()[1])
           #print(line_gene in strand_dict.keys())
           if line_gene in strand_dict.keys() and not line_gene in broken_genes:
               strand = strand_dict[line_gene][0]
               file_list = next(walk(genome+"/CDS/"+line_gene), (None, None, []))[2]# save each transcript filename to a list
               hit=False
               i=0
               reading_frame_list=[]
               transcript_name_list=[]
               temp_order_list=[]
               while i < len(file_list):
                   with open(genome+"/CDS/"+line_gene+"/"+file_list[i]) as f:
                       cds = f.readline().strip('\n')
                       cds = cds.split(', ')
                       cds=list(map(int, cds))
                       if line_pos in cds:# open each file, check if CDS position is in the file
                           hit=True
                           temp_order=list(sorted(cds))
                           temp_order_list.append(temp_order)
                           reading_frame_list.append(test_readingframe(temp_order, strand, line_pos))
                           transcript_name=file_list[i]
                           transcript_name=transcript_name.replace('.txt','')
                           transcript_name_list.append(transcript_name)
                   i=i+1
               
               if hit == True:
                   if not (len(temp_order) % 3 == 0):
                       if not line_gene in genes_checked:
                           print("Error: "+line_gene+" has a coding region length not divisible by 3! ("+str(len(temp_order))+")") # This error should not be possible if the genome was set up properly.
                           genes_checked.append(line_gene)
                   else:
                       if not line_gene in genes_checked:
                           genes_checked.append(line_gene)
                           
                       # Extract the most common reading frame:
                       common=most_frequent(reading_frame_list)
                       v = reading_frame_list.index(common)
                       transcript_name=transcript_name_list[v]
                       temp_order=temp_order_list[v]
                       cds_pos = int((temp_order.index(line_pos)))
                       
                       if strand == "+":
                           if((cds_pos+1) % 3 == 0):
                               ref_codon = str(gChr[geneChrom][temp_order[cds_pos-2]-1])+str(gChr[geneChrom][temp_order[cds_pos-1]-1])+ctl_geno
                               alt_codon = str(gChr[geneChrom][temp_order[cds_pos-2]-1])+str(gChr[geneChrom][temp_order[cds_pos-1]-1])+trt_geno
                           elif((cds_pos) % 3 == 0):
                               ref_codon = ctl_geno+str(gChr[geneChrom][temp_order[cds_pos+1]-1])+str(gChr[geneChrom][temp_order[cds_pos+2]-1])
                               alt_codon = trt_geno+str(gChr[geneChrom][temp_order[cds_pos+1]-1])+str(gChr[geneChrom][temp_order[cds_pos+2]-1])
                           elif((cds_pos+2) % 3 == 0):
                               ref_codon = str(gChr[geneChrom][temp_order[cds_pos-1]-1])+ctl_geno+str(gChr[geneChrom][temp_order[cds_pos+1]-1])
                               alt_codon = str(gChr[geneChrom][temp_order[cds_pos-1]-1])+trt_geno+str(gChr[geneChrom][temp_order[cds_pos+1]-1])
                           else:
                               print("oh no.")
                       else:
                           cds_pos_neg = (len(temp_order)-1)-int((temp_order.index(line_pos))) #previous: len(temp_order))-1
                           if((cds_pos_neg) % 3 == 0):
                               ref_codon = ctl_geno+str(gChr[geneChrom][temp_order[cds_pos-1]-1])+str(gChr[geneChrom][temp_order[cds_pos-2]-1])
                               alt_codon = trt_geno+str(gChr[geneChrom][temp_order[cds_pos-1]-1])+str(gChr[geneChrom][temp_order[cds_pos-2]-1])
                           elif((cds_pos_neg+2) % 3 == 0):
                               ref_codon = str(gChr[geneChrom][temp_order[cds_pos+1]-1])+ctl_geno+str(gChr[geneChrom][temp_order[cds_pos-1]-1])
                               alt_codon = str(gChr[geneChrom][temp_order[cds_pos+1]-1])+trt_geno+str(gChr[geneChrom][temp_order[cds_pos-1]-1])
                           elif((cds_pos_neg+1) % 3 == 0):
                               ref_codon = str(gChr[geneChrom][temp_order[cds_pos+2]-1])+str(gChr[geneChrom][temp_order[cds_pos+1]-1])+ctl_geno
                               alt_codon = str(gChr[geneChrom][temp_order[cds_pos+2]-1])+str(gChr[geneChrom][temp_order[cds_pos+1]-1])+trt_geno
                           else:
                               print("oh no.")
                           # take direct compliment of reversed codon:
                           ref_codon = list(ref_codon)
                           alt_codon = list(alt_codon)
                           i=0
                           while i < (len(ref_codon)):
                               ref_codon[i] = compliment(ref_codon[i])
                               i=i+1
                           i=0
                           while i < (len(alt_codon)):
                               alt_codon[i] = compliment(alt_codon[i])
                               i=i+1
                           ref_codon = "".join(ref_codon)
                           alt_codon = "".join(alt_codon)
                           
                       # translate codons into peptides:
                       if(not len(ref_codon)==3 or not len(alt_codon)==3):# if indel:
                           
                           ref_true = ''.join(line[refcol])
                           rf_ref=(len(ref_true)%3)/3 #checks if we are in the same reading frame as the reference allele
                           rf_ctl=(len(ctl_geno)%3)/3
                           rf_trt=(len(trt_geno)%3)/3
                           
                           if(len(ctl_geno) == len(ref_true) or rf_ctl == rf_ref):
                               ref_pep = "In Frame"
                           else:
                               ref_pep = "Frameshift"
                           if(len(trt_geno) == len(ref_true) or rf_trt == rf_ref):
                               alt_pep = "In Frame"
                           else:
                               alt_pep = "Frameshift"
                           
                           granth="NA"
                           
                       else:# if SNP:
                           ref_pep = translate(ref_codon)
                           alt_pep = translate(alt_codon)
                           if(ref_pep == alt_pep): # ID mutation types by comparing peptides:
                               granth = 0
                           else:
                               if(ref_pep == "_"):
                                   granth = "nonstop"
                               elif(alt_pep == "_"):
                                   granth = "nonsense"
                               else:
                                   granth = grantham_dict[ref_pep, alt_pep]
                       ref_codon_list.append(ref_codon) # save positional info for output
                       alt_codon_list.append(alt_codon)
                       ref_peptide_list.append(ref_pep)
                       alt_peptide_list.append(alt_pep)
                       grantham_dist_list.append(granth)
                       transcript_list.append(transcript_name)
               else:
                   # NOT IN A CDS
                   ref_codon_list.append("Not in CDS")
                   alt_codon_list.append("Not in CDS")
                   ref_peptide_list.append("NA")
                   alt_peptide_list.append("NA")
                   grantham_dist_list.append("NA")
                   transcript_list.append('NA')
           else:
               if line_gene in broken_genes:
                   ref_codon_list.append("Gene CDS length error")
                   alt_codon_list.append("Gene CDS length error")
                   ref_peptide_list.append("NA")
                   alt_peptide_list.append("NA")
                   grantham_dist_list.append("NA")
                   transcript_list.append('NA')
               else:
                   ref_codon_list.append("Gene not recognized")
                   alt_codon_list.append("Gene not recognized")
                   ref_peptide_list.append("NA")
                   alt_peptide_list.append("NA")
                   grantham_dist_list.append("NA")
                   transcript_list.append('NA')
   #print(chrom_list, pos_list, pos_gene_list, parent_geno_list, child_geno_list, ref_codon_list, alt_codon_list, ref_peptide_list, alt_peptide_list, grantham_dist_list)
   print("Writing output file...")
   df = pd.DataFrame(list(zip(chrom_list, pos_list, pos_gene_list, parent_geno_list, child_geno_list, ref_codon_list, alt_codon_list, ref_peptide_list, alt_peptide_list, grantham_dist_list, transcript_list)),
                  columns =['CHROM', 'POS', 'Gene', 'Control_Genotype', "Treat_Genotype", "Control_Codon", "Treat_Codon", 'Control_Peptide', 'Treat_Peptide', 'Grantham', 'Transcript'])
   df.to_csv(outputfile, sep="\t", index=False) #print(df)
   
if __name__ == "__main__":
   main(sys.argv[1:])
