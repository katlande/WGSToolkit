#!/usr/bin/env python
import csv
import sys, getopt
import collections
import pandas as pd
from sys import exit
import os
pd.options.mode.chained_assignment = None

# Functions:

def check_transcript(CDSonly, t):
    test = CDSonly[CDSonly[7].isin([t])]
    test[4] = test[4] + 1
    if((sum(test[4] - test[3])) % 3 == 0):
        return(True)
    else:
        return(False)

# Command Line Inputs:
def main(argv):
   Name=''
   Gtf=''
   
   try:
      opts, args = getopt.getopt(argv, "hn:g", ["name=", "gtf="])
   except getopt.GetoptError:
      print('')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Script Required Inputs ---\n\n\t--name\t<genome name shorthand>\n\n\t--gtf\t<gtf file>\n\n')
         sys.exit()
      elif opt in ("-n", "--name"):
         Name = arg
      elif opt in ("-g", "--gtf"):
         Gtf = arg
   
   # Main code body
   print("Reading CDS sites...")
   gtf_file=pd.read_csv((Gtf), sep="\t", header=None, low_memory=False, comment='#')
   CDSonly = gtf_file[gtf_file[2] == 'CDS']
   
   total_rows=len(CDSonly[8])
   name_rows = len([x for x in CDSonly[8] if "transcript_name" in x])
   id_rows = len([x for x in CDSonly[8] if "transcript_id" in x])
   print("Separating transcripts...")
   transcript_aware=True
   if name_rows == total_rows:
       print("Parsing by transcript name...")
       CDSonly[7] = CDSonly[8].str.replace('.+transcript_name "', "")
       CDSonly[7] = CDSonly[7].str.replace('";.+', "")
   elif id_rows == total_rows:
       print("Parsing by transcript ID...")
       CDSonly[7] = CDSonly[8].str.replace('.+transcript_id "', "")
       CDSonly[7] = CDSonly[7].str.replace('";.+', "")
   else:
       if name_rows > 0 or id_rows > 0:
           print ("Incomplete transcript information.")
           print(str(name_rows)+"/"+str(total_rows)+" CDS regions with transcript names.")
           print(str(id_rows)+"/"+str(total_rows)+" CDS regions with transcript IDs.")
           print("Proceeding as if all gene models represent singular transcripts. This may causes reading frame errors in tools that translate mutant gene products such as Call_Codons.")
       else:
           print("Error: No transcript information found in GTF. Proceeding as if all gene models represent singular transcripts. This may causes reading frame errors in tools that translate mutant gene products such as Call_Codons.")
       transcript_aware=False
    
   #print(CDSonly.head(10))
   print("Reformatting CDS data...")
   CDSonly[8] = CDSonly[8].str.replace('.+gene_name "', "")
   CDSonly[8] = CDSonly[8].str.replace('";.+', "")
   CDSonly[8] = CDSonly[8].str.replace('gene_id "', "") # format weirdo gene names too, arabidopsis is breaking my heart
   CDSonly[8] = CDSonly[8].str.replace('/', "-") #I cant believe I have to include this line. Arabidopsis scientists are satanists.
   genes=list(set(list(CDSonly[8])))
   if(len(genes) > 5000):
       print(len(genes), " genes to query... (this may take a while)")
   else:
       print(len(genes), " genes to query...")
   
   if transcript_aware == True:
       print("Parsing CDS regions by transcript and writing indices...")
       false_cds = []
       for g in genes: # For each gene:
           gdf = CDSonly[CDSonly[8].isin([g])]
           os.mkdir("./CDS/"+g) # Make folder in CDS for gene name
           transcripts=list(set(list(gdf[7])))
           # for each transcript:
           t_count=0
           for t in transcripts:
               reading_frame=check_transcript(CDSonly, t) # check if transcript length is divisible by three, if True:
               if reading_frame == True:
                   
                   tdf = CDSonly[CDSonly[7] == t]
                   list_one_gene=[]
                   for index, row in tdf.iterrows(): # iterate over rows
                       rangeOBJ=[*range(int(row[3]), int(row[4])+1,1)] # fill values in between start and end
                       list_one_gene.append(', '.join(str(x) for x in (rangeOBJ)))
                   
                   rangeOBJ=', '.join(str(x) for x in (list_one_gene))
                   #print(rangeOBJ.split(",")) # once this is fixed the other 2 conditions need to be reformated
                   rangeOBJ = rangeOBJ.split(",")
                   rangeOBJ = list(sorted(list(map(int, set(list(rangeOBJ))))))
                   rangeOBJ=', '.join(str(x) for x in (list_one_gene))
                   val=[]
                   val.append(rangeOBJ)
                   df = pd.DataFrame(val, columns =[" "])
                   df.to_csv("./CDS/"+g+"/"+t+".txt", sep="\t", index=False, header=False)
                   
                   # take ranges, sort, unique, order, and write to file named transcript.txt
                   t_count=t_count+1
           
           if t_count == 0:
               # if all individual transcripts fail, check the entire gene
               list_one_gene=[]
               for index, row in gdf.iterrows(): # iterate over rows
                   rangeOBJ=[*range(int(row[3]), int(row[4])+1,1)] # fill values in between start and end
                   list_one_gene.append(', '.join(str(x) for x in (rangeOBJ)))
               
               rangeOBJ=', '.join(str(x) for x in (list_one_gene))
               rangeOBJ = rangeOBJ.split(",")
               #print(rangeOBJ)
               rangeOBJ = list(sorted(list(map(int, set(rangeOBJ)))))
               if not len(rangeOBJ) % 3 == 0:
                   false_cds.append(g)
                   #print("ERROR: "+g+" has NO transcripts with a total CDS length divisible by 3! Skipping this gene.")
               else:
                   rangeOBJ=', '.join(str(x) for x in (list_one_gene))
                   val=[]
                   val.append(rangeOBJ)
                   df = pd.DataFrame(val, columns =[" "])
                   df.to_csv("./CDS/"+g+"/"+g+".txt", sep="\t", index=False, header=False)
   else:
       
       for g in genes:
           gdf = CDSonly[CDSonly[8].isin([g])]
           os.mkdir("./CDS/"+g)
           # check if gene length is divisible by three:
           list_one_gene=[] 
           for index, row in gdf.iterrows(): # iterate over rows
               rangeOBJ=[*range(int(row[3]), int(row[4])+1,1)] # fill values in between start and end
               list_one_gene.append(', '.join(str(x) for x in (rangeOBJ)))
           
           rangeOBJ=', '.join(str(x) for x in (list_one_gene))
           rangeOBJ = rangeOBJ.split(",")
           rangeOBJ = list(sorted(list(map(int, set(rangeOBJ)))))
           if not len(rangeOBJ) % 3 == 0:
               false_cds.append(g)
               #print("ERROR: "+g+" has a total CDS length non-divisible by 3! Skipping this gene.")
           else:
               rangeOBJ=', '.join(str(x) for x in (list_one_gene))
               val=[]
               val.append(rangeOBJ)
               df = pd.DataFrame(val, columns =[" "])
               df.to_csv("./CDS/"+g+"/"+g+".txt", sep="\t", index=False, header=False)
           
   #PCgenes = pd.unique(CDSonly[8])
   
   print("Finding gene strands...")
   strands = CDSonly.drop(CDSonly.columns[[0, 1, 2, 3, 4, 5, 7]], axis=1) #should remove all except gene name and strand
   strands = strands.drop_duplicates()
   #print(strands.head())
   
   print("Writing strands to file...")
   strands = strands[strands.columns[[1,0]]]
   strands.to_csv(Name+"_ProtCoding_Gene_Strands.txt", sep="\t", header=False, index=False)
   
   print("Writing broken protein models to file...")
   df = pd.DataFrame(false_cds, columns =["Broken_Transcripts"])
   df.to_csv(Name+"_Broken_Gene_Models.txt", sep="\t", index=False, header=False) #print(df)
   
   print(str(len(false_cds))+"/"+str(len(genes))+" skipped genes had CDS region lengths that were not divisible by 3.")
if __name__ == "__main__":
   main(sys.argv[1:])

