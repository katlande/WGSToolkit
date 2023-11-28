from Bio import pairwise2
import random
from Bio.pairwise2 import align
from pyfaidx import Fasta
import sys, getopt
import csv
import pandas as pd

def One_Comp(index, c, s, e, p, l, Verbose=False):
    e=int(e)
    s=int(s)
    l=int(l)
    if l > (e-s)/1.5:
        while l > (e-s)/1.5:
            l = int(l/2)
            print("Warning: Read length is too long for this region size. Shrinking read length to: "+str(l)+"bp.")
    cushion = int(0.2*l) # amount of padding to place on either side of the region when pulling the starting index
    start_val = int(s)-(l-cushion)
    end_val = int(e)+(l-cushion)
    fasta = str(index[str(c)][start_val:end_val])
    
    score_list = []
    i = 0
    if Verbose==True:
        print("Running "+str(p)+" permutations...")
    
    while i < p:
        # Random number as a to get read position:
        i_start = random.randrange(0, (len(fasta)-l-1), 1)
        i_end = i_start+l
        i_start2 = random.randrange(0, (len(fasta)-l-1), 1)
        i_end2 = i_start2+l
        # pull these regions from the fasta:
        seq1 = fasta[i_start:i_end]
        seq2 = fasta[i_start2:i_end2]
        # align two sequences:
        score = pairwise2.align.globalxx(seq1, seq2, score_only=True, penalize_end_gaps=False)
        score=(score/l) # roughly equal to percent of the read in agreement
        score_list.append(score)
        i=i+1
    
    out_score=sum(score_list)/len(score_list)
    variance = sum([((x - out_score) ** 2) for x in score_list]) / len(score_list)
    sd = variance ** 0.5

    sd=sd*100
    out_score=out_score*100
    return([out_score, sd])

def main(argv):
   
   try:
      opts, args = getopt.getopt(argv, "hp:l:f:c:s:e:b:j:i:o", ["permutes=","readlen=","fasta=","chrom=","start=","end=","batch=","flank_jitter=","input=","output="])
   except getopt.GetoptError:
      print('error')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('error')
         sys.exit()
      elif opt in ("-p", "--permutes"):
         perm = int(arg)
      elif opt in ("-l", "--readlen"):
         readlen = int(arg)
      elif opt in ("-f", "--fasta"):
         fasta_path = str(arg)
      elif opt in ("-c", "--chrom"):
         chrom = str(arg)
      elif opt in ("-s", "--start"):
         start = int(arg)
      elif opt in ("-e", "--end"):
          end = int(arg)
      elif opt in ("-b", "--batch"):
          batch = str(arg)
      elif opt in ("-j", "--flank_jitter"):
          flank = int(arg)
      elif opt in ("-i", "--input"):
          infile = arg
      elif opt in ("-o", "--output"):
          outfile = arg
   
   if batch=="b":
       print("Batch mode...")
   elif batch=="s":
        next
   else:
       print("Error: Mode option "+batch+" not recognized.")
       exit()
   
   # Read fasta in:
   # whatever command
   # from int(start)-(readlen-cushion) to end+(readlen-cushion)
   print("Reading reference genome...")
   gChr = Fasta(fasta_path)
   print("Done!")
   
   if batch=="s":
       out=One_Comp(gChr, chrom, start, end, perm, readlen, Verbose=True)
       print("")
       print("")
       print("In "+str(perm)+" permutations of "+str(readlen)+"bp, the mean complexity score is: "+str(out[0])+" with st dev "+str(out[1]))
       print("A complexity score of 100 means no variability whatsoever. Higher deviations in complexity score indicate a mix of low and high complexity sequences.\n")
       print("--------- Legend ---------")
       print("(From 100,000 permutations)")
       print("------------------------------------------------------------------------")
       print("readlen:          10bp           50bp           100bp          150bp")
       print("------------------------------------------------------------------------")
       print("3bp repeat:    93.3+/-4.70    98.7+/-0.94    99.3+/-0.47    99.6+/-0.31")
       print("4bp repeat:    90.0+/-7.08    98.0+/-1.42    99.0+/-0.71    99.3+/-0.47")
       print("5bp repeat:    88.0+/-7.49    97.6+/-1.50    98.8+/-0.75    99.2+/-0.50")
       print("10bp repeat:   75.1+/-15.0    95.0+/-3.00    97.5+/-1.50    98.3+/-1.00")
       print("------------------------------------------------------------------------")
       print("Random:        52.8+/-10.2    60.5+/-4.22    62.3+/-3.81    63.2+/-4.18")
       print("Semi-Random:   53.8+/-11.5    60.5+/-7.24    62.0+/-6.26    62.8+/-6.26")
       print("hg38 (chr2):   51.4+/-11.2    58.9+/-5.22    60.8+/-4.55    61.7+/-4.61")
       print("mm10 (chr2):   51.5+/-11.7    59.0+/-5.67    60.9+/-4.65    61.9+/-4.64")
       print("------------------------------------------------------------------------")
       print("* ~4kb regions. This is NOT a reference tool. Compare regions of interest to flanking regions of comparable size.")
   elif batch=="b":
       chrom_list=[]
       start_list=[]
       end_list=[]
       cRegion_list=[]
       sdRegion_list=[]
       cFlank_list=[]
       sdFlank_list=[]
       
       # read in infile, for each row, chr, start, end
       with open(infile) as gList:    # for each gene
          readfile=csv.reader(gList, delimiter="\t")
          for line in readfile:
              chrom=str(''.join(line[0]))
              start=int(''.join(line[1]))
              end=int(''.join(line[2]))
              chrom_list.append(chrom)
              start_list.append(start)
              end_list.append(end)
              
              out_region=One_Comp(gChr, chrom, start, end, perm, readlen, Verbose=False)
              out_ctl=One_Comp(gChr, chrom, (start+flank), (end+flank), perm, readlen, Verbose=False)
              
              cRegion_list.append(out_region[0])
              sdRegion_list.append(out_region[1])
              cFlank_list.append(out_ctl[0])
              sdFlank_list.append(out_ctl[1])
          
          df = pd.DataFrame(list(zip(chrom_list, start_list, end_list, cRegion_list, cFlank_list, sdRegion_list, sdFlank_list)),
                         columns =['CHROM', 'START', 'END', 'Complexity', 'Complexity_Flank', 'sd', "sd_Flank"])
          df.to_csv(outfile, sep="\t", index=False)
   else:
       print("Script error: contact klande@salk.edu")
   
if __name__ == "__main__":
   main(sys.argv[1:])

