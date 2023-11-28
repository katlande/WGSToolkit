# This script takes the SNP_Genotypes.txt file output from R and finds distinct sites between two user-defined groups

import csv
import sys, getopt
import collections
import pandas as pd

# Command Line Inputs:
def main(argv):
   inputfile = ''
   outputfile = ''
   Group1 = ''
   Group2 = ''
   Group1_name = "Group1"
   Group2_name = "Group2"
   max_missing = 0.5
   consensus_freq = 0.75
   anal_mode="SNP"
   try:
      opts, args = getopt.getopt(argv, "hi:o:a:b:c:d:m:f:v", ["input=","output=","group1=","group2=","name1=","name2=","max_missing=","consensus_freq=","mode="])
   except getopt.GetoptError:
      print('\nDifferential_Variants.py\n-i <inputfile>\n-o <outputfile>\n-g1 <group 1 (as column vector)>\n-g2 <group 2 (as column vector)>\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Differential Variants Required Inputs ---\n\n\t-i | --input\t<inputfile as a tab-delimited text file>\n\t-o | --output\t<outputfile>\n\t-a | --group1\t<group 1 (as column vector)>\n\t-b | --group2\t<group 2 (as column vector)>\n\n--- Differential Variants Optional Inputs ---\n\n\t-c | --name1\t\t<group 1 name>, default="Group1"\n\t-d | --name2\t\t<group 2 name>, default="Group2"\n\t-m | --max_missing\t<max allowed missing sample frequency per group, a value between 0 and 1>, default=0.5\n\t-f | --consensus_freq\t<min genotype frequency per group required for consensus, a value > 0.5 and ≤ 1>, default=0.75\n\t-v | --mode\t<SNP or INDEL>, default=SNP\n\n')
         sys.exit()
      elif opt in ("-i", "--input"):
         inputfile = arg
      elif opt in ("-o", "--output"):
         outputfile = arg
      elif opt in ("-a", "--group1"):
          Group1 = arg
      elif opt in ("-b", "--group2"):
          Group2 = arg
      elif opt in ("-c", "--name1"):
          Group1_name = arg
      elif opt in ("-d", "--name2"):
          Group2_name = arg
      elif opt in ("-m", "--max_missing"):
          max_missing = float(arg)
      elif opt in ("-f", "--consensus_freq"):
          consensus_freq = float(arg)
      elif opt in ("-v", "--mode"):
          anal_mode = arg
          
   print('\n\n--------------------------------\nFinding differential variants with parameters:\nInput=', inputfile, "\nOutput=", outputfile, "\nGroup 1 ( named", Group1_name,") spans input columns ", Group1, "\nGroup 2 ( named", Group2_name,") spans input columns ", Group2, "\nMax missing variants allowed in each group = ", max_missing, "\nMinimum genotype frequency to form a group consensus = ", consensus_freq,"\n--------------------------------\n\n")
   
   if anal_mode == "SNP":
       print("Running SNP mode. Zygosity predictions will not be accurate for variants longer than 1bp")
   else:
       print("Running INDEL mode. ")
   
   if consensus_freq <= 0.5:
       print("ERROR: consensus_freq must be above 0.5")
       sys.exit(2)
   
   # set up empty objects for output file:
   chrom = []
   pos = []
   g1_freqMISSING = []
   g1_freqCONSENSUS = []
   g1_GENOTYPE = []
   g1_ZYGOSITY = []
   g2_freqMISSING = []
   g2_freqCONSENSUS = []
   g2_GENOTYPE = []
   g2_ZYGOSITY = []   
   
   Group1 = Group1.strip('[]').split(',')
   Group1 = [int(i)-1 for i in Group1]
   Group2 = Group2.strip('[]').split(',')
   Group2 = [int(i)-1 for i in Group2]
   
   with open(inputfile) as tsv:
       df=csv.reader(tsv, delimiter="\t")
       next(df)
       print("Reading input file...")
       for line in df:
           #print(line)
           lineg1 = [line[i] for i in Group1] # For one line, gives the genotypes of samples in group 1 as a list
           lineg2 = [line[i] for i in Group2]
           g1freq = dict(collections.Counter(lineg1)) # For one line, gives the frequencies of group 1 genotypes
           g2freq = dict(collections.Counter(lineg2))
           
           # check for missing samples in the line:
           if "./." in g1freq:
               g1_missing = int(g1freq["./."])/len(lineg1)
               lineg1.remove('./.') #only deletes the first entry for some reason? this is trash. fixed lines 98-99
           else:
               g1_missing = 0
               
           if "./." in g2freq:
               g2_missing = int(g2freq["./."])/len(lineg2)
               lineg2.remove('./.')
           else:
               g2_missing = 0
           
           # continue analyzing the line only if the missing sample number is below the max missing threshold
           if g2_missing <= max_missing and g1_missing <= max_missing:
               # Next check to see if the consensus genotype of each group is above the threshold
               g1freq_filt = dict(collections.Counter(lineg1)) 
               g2freq_filt = dict(collections.Counter(lineg2))
               
               if './.' in g1freq_filt: del g1freq_filt['./.']# remake genotype dicts without ./. samples
               if './.' in g2freq_filt: del g2freq_filt['./.']# remake genotype dicts without ./. samples
               
               g1CF = int(list(g1freq_filt.values())[list(g1freq_filt.values()).index(max(list(g1freq_filt.values())))])/sum(g1freq_filt.values())
               g2CF = int(list(g2freq_filt.values())[list(g2freq_filt.values()).index(max(list(g2freq_filt.values())))])/sum(g2freq_filt.values())
               
               # If the genotype frequencies are above the threshold in both groups, continue.
               if g1CF >= consensus_freq and g2CF >= consensus_freq:
                   #Check if the groups have the same consensus genotype:
                   g1_allele = list(g1freq_filt.keys())[list(g1freq_filt.values()).index(max(list(g1freq_filt.values())))]
                   g2_allele = list(g2freq_filt.keys())[list(g2freq_filt.values()).index(max(list(g2freq_filt.values())))]
                   #print("Group 1:", g1_allele, " | Group 2:", g2_allele)
                   
                   # If the consensus sites are different alleles, make outputs:
                   if not g1_allele == g2_allele:
                       #print("Group1", g1freq_filt, g1CF)
                       #print("Group2", g2freq_filt, g2CF)
                       #print("Differential site!")
                       
                       # extract genotypes as nucleotides:
                       geno_g1 = str(line[int(g1_allele[0])+2])+"/"+str(line[int(g1_allele[2])+2])
                       geno_g2 = str(line[int(g2_allele[0])+2])+"/"+str(line[int(g2_allele[2])+2]) 
                       
                       chrom.append(line[0])
                       pos.append(line[1])
                       #print(chrom, pos)
                       g1_freqMISSING.append(g1_missing)
                       g1_freqCONSENSUS.append(g1CF)
                       g1_GENOTYPE.append(geno_g1)
                       if geno_g1[0] == geno_g1[2]:
                           g1_ZYGOSITY.append("HOMO")
                       else:
                           g1_ZYGOSITY.append("HETERO")
                       
                       #print(g1_GENOTYPE, g1_ZYGOSITY)
                       g2_freqMISSING.append(g2_missing)
                       g2_freqCONSENSUS.append(g2CF)
                       g2_GENOTYPE.append(geno_g2)
                       if geno_g2[0] == geno_g2[2]:
                           g2_ZYGOSITY.append("HOMO")
                       else:
                           g2_ZYGOSITY.append("HETERO")
               
               else: # if samples are not above the threshold, check with a fine toothed comb;
                   # if one sample is strictly reference/homozygous, and the other sample is a mix of heterozygous + homozyous/reference, it should still be considered differential
                   if g1CF == 1 or g2CF == 1: # if group 1 OR group 2 is 100% 0/0 or 1/1:
                       g1_allele = list(g1freq_filt.keys())[list(g1freq_filt.values()).index(max(list(g1freq_filt.values())))]
                       g2_allele = list(g2freq_filt.keys())[list(g2freq_filt.values()).index(max(list(g2freq_filt.values())))]
                       
                       # if the 100% sample is homozygous, and the other sample is at least 50% one genotype and contains NO samples with the other group's genotype, we will also consider this a "polymorphic hit"
                       #print(g1CF, g1_allele, g2CF, g2_allele)
                       #quit()
                       if g1CF == 1 and not g2CF == 1:
                           if g1_allele[0] == g1_allele[2] and g2CF >= 0.5: #if homozygous....
                              if not g1_allele in g2freq_filt.keys():
                                  chrom.append(line[0])
                                  pos.append(line[1])
                                  g1_freqMISSING.append(g1_missing)
                                  g1_freqCONSENSUS.append(g1CF)
                                  g1_GENOTYPE.append(g1_allele + " | " + str(line[int(g1_allele[0])+2])+"/"+str(line[int(g1_allele[2])+2]))
                                  g1_ZYGOSITY.append("HOMO")
                                  g2_freqMISSING.append(g2_missing)
                                  g2_freqCONSENSUS.append(str(g2freq_filt))
                                  g2_GENOTYPE.append("NA")
                                  g2_ZYGOSITY.append("Polymorphic")
                       elif g2CF == 1 and not g1CF == 1:
                           if g2_allele[0] == g2_allele[2] and g1CF >= 0.5: #if homozygous....
                               if not g2_allele in g1freq_filt.keys():
                                   chrom.append(line[0])
                                   pos.append(line[1])
                                   g1_freqMISSING.append(g1_missing)
                                   g1_freqCONSENSUS.append(str(g1freq_filt))
                                   g1_GENOTYPE.append("NA")
                                   g1_ZYGOSITY.append("Polymorphic")
                                   g2_freqMISSING.append(g2_missing)
                                   g2_freqCONSENSUS.append(g2CF)
                                   g2_GENOTYPE.append(g2_allele + " | " + str(line[int(g2_allele[0])+2])+"/"+str(line[int(g2_allele[2])+2]))
                                   g2_ZYGOSITY.append("HOMO")
                       else:
                           print("ERROR: CONSENSUS FREQUENCY SETTING IS ABOVE 100%\nTERMINATING...")
                           quit()
                    
                   
                   
       # Write output file:
       print("Writing output file...")
       
       if anal_mode == "SNP":
           df = pd.DataFrame(list(zip(chrom, pos, g1_freqMISSING, g1_freqCONSENSUS, g1_GENOTYPE, g1_ZYGOSITY, g2_freqMISSING, g2_freqCONSENSUS, g2_GENOTYPE, g2_ZYGOSITY)),
                          columns =['CHROM', 'POS', Group1_name+"_percMissing", Group1_name+"_ConsensusFreq", Group1_name+"_Genotype", Group1_name+"_Zygosity",  Group2_name+"_percMissing", Group2_name+"_ConsensusFreq", Group2_name+"_Genotype", Group2_name+"_Zygosity"])
       else:
           df = pd.DataFrame(list(zip(chrom, pos, g1_freqMISSING, g1_freqCONSENSUS, g1_GENOTYPE, g2_freqMISSING, g2_freqCONSENSUS, g2_GENOTYPE)),
                          columns =['CHROM', 'POS', Group1_name+"_percMissing", Group1_name+"_ConsensusFreq", Group1_name+"_Genotype", Group2_name+"_percMissing", Group2_name+"_ConsensusFreq", Group2_name+"_Genotype"])
       
       
       #print(df)
       df.to_csv(outputfile, sep="\t")
       
       
       

if __name__ == "__main__":
   main(sys.argv[1:])




