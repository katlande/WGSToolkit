#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# This is a nested script
cat("\n----- Extracting Genotype Information to Plain Text -----\n\n")
cat("Loading libraries...\n")
suppressMessages(library("tidyr", warn.conflicts = F, quietly = T))
suppressMessages(library("dplyr", warn.conflicts = F, quietly = T))
suppressMessages(library("reshape2", warn.conflicts = F, quietly = T))
suppressMessages(library("stringr", warn.conflicts = F, quietly = T))
setwd(args[1])

# main apply function:
apply_SNP <- function(x){
  temp <- unlist(str_split(unlist(x), "/"))# note - cannot simply use x[1], x[3] because this needs to work for any possible ploidy; 
  REF_AF <- length(temp[temp=="0"])/length(temp) # % of alleles = ref
  ALT_AF <- length(temp[temp=="1"])/length(temp) # % of alleles = alt
  NB_AF <- length(temp[! temp%in%c("0","1",".")])/length(temp) # % of alleles != ref, alt, or missing (i.e., % multiallelic alleles)
  
  NB_AF <- ifelse(NB_AF == 0, "biallelic", paste0("multiAlleleFreq:",NB_AF)) # label biallelic or multiallelic with multiallelic AF
  
  return(ifelse((length(temp[temp=="."])/length(temp)) > 0.5, 
                paste0("MISSING;MISSING;MISSING"), #ignore allele is >half of alleles are missing
                paste0(REF_AF, ";", ALT_AF, ";", NB_AF)))# otherwise returns ref, alt, allele type label for one sample~locus
}

# Clean up sample names:
sample_names <- read.delim("tmp.list", header=F)$V1
sample_names <- str_remove_all(sample_names, "./gvcf/")
sample_names <- str_remove_all(sample_names, ".g.vcf.gz")
sample_names <- str_remove_all(sample_names, "_R1")
sample_names <- str_remove_all(sample_names, "_R2")
sample_names <- str_remove_all(sample_names, "_001")
#gsub("_CSF.*", "", sample_names)

cat("Reading filtered VCFs...\n")
VCF_SNP <- read.delim("Filtered_SNPs.recode.vcf", comment.char = "#", header=F, stringsAsFactors = F)
geno_SNP <- VCF_SNP[c(1,2,4,5,10:ncol(VCF_SNP))]
colnames(geno_SNP) <- c("CHROM",	"POS",	"REF",	"ALT", sample_names[!sample_names=="All.Samples"])
VCF_indel <- read.delim("Filtered_Indels.recode.vcf", comment.char = "#", header=F, stringsAsFactors = F)
geno_indel <- VCF_indel[c(1,2,4,5,10:ncol(VCF_indel))]
colnames(geno_indel) <- c("CHROM",	"POS",	"REF",	"ALT", sample_names[!sample_names=="All.Samples"])

# make sure chromosome names start with "chr" - will cause issues later if they don't
geno_SNP$CHROM <- paste0("chr",geno_SNP$CHROM)
geno_SNP$CHROM <- str_replace_all(geno_SNP$CHROM, "chrchr", "chr")
geno_indel$CHROM <- paste0("chr",geno_indel$CHROM)
geno_indel$CHROM <- str_replace_all(geno_indel$CHROM, "chrchr", "chr")


cat(paste0(args[2],"/region_sizes.txt"))
t<-paste0(args[2],"/region_sizes.txt")
if(file.exists(t)){
  usefulChrom <- read.delim(t, stringsAsFactors = F)
  usefulChrom <- usefulChrom$Region
  cat("Keeping the following chromosomes:\n")
  cat(usefulChrom)
  geno_SNP <- geno_SNP[geno_SNP$CHROM %in% usefulChrom,]
  geno_indel <- geno_indel[geno_indel$CHROM %in% usefulChrom,] 
}

cat("\nRemoving non-essential information from VCF...\n")
geno_SNP[c(5:ncol(geno_SNP))] <- apply(geno_SNP[c(5:ncol(geno_SNP))], 2, function(y) (gsub("\\:.*", "", y)))
geno_indel[c(5:ncol(geno_indel))] <- apply(geno_indel[c(5:ncol(geno_indel))], 2, function(y) (gsub("\\:.*", "", y)))

cat("\nCalculating frequency of SNP file...\n")
cat(paste0(length(5:ncol(geno_SNP)), " Samples to reformat... "))
geno_SNP[c(5:ncol(geno_SNP))] <- apply(geno_SNP[c(5:ncol(geno_SNP))], c(1,2), FUN=apply_SNP)
cat("Done!\n Reformatting SNP file for output...\n")
runV <- 5:ncol(geno_SNP)
runV <- ((seq_along(runV)-1)*2)+runV
for(c in runV){
  separate(geno_SNP, col=c, into=c(paste0("REF_Freq_",colnames(geno_SNP)[c]),
                                   paste0("ALT_Freq",colnames(geno_SNP)[c]),
                                   paste0("Other_Freq",colnames(geno_SNP)[c])), sep="\\;") -> geno_SNP
}
rm(c, runV)
geno_SNP<-geno_SNP[rowSums(geno_SNP=='MISSING')==0,]
geno_SNP[c(5:ncol(geno_SNP))] <- apply(geno_SNP[c(5:ncol(geno_SNP))], 2, function(x){sub("^$", 0, x)})

cat("\nCalculating frequency of indel file...\n")
cat(paste0(length(5:ncol(geno_indel)), " Samples to reformat..."))
geno_indel[c(5:ncol(geno_indel))] <- apply(geno_indel[c(5:ncol(geno_indel))], c(1,2), FUN=apply_SNP)
cat("Done!\n Reformatting indel file for output...\n")
runV <- 5:ncol(geno_indel)
runV <- ((seq_along(runV)-1)*2)+runV
for(c in runV){
  separate(geno_indel, col=c, into=c(paste0("REF_Freq_",colnames(geno_indel)[c]),
                                     paste0("ALT_Freq",colnames(geno_indel)[c]),
                                     paste0("Other_Freq",colnames(geno_indel)[c])), sep="\\;") -> geno_indel
}
rm(c, runV)
geno_indel<-geno_indel[rowSums(geno_indel=='MISSING')==0,]
geno_indel[c(5:ncol(geno_indel))] <- apply(geno_indel[c(5:ncol(geno_indel))], 2, function(x){sub("^$", 0, x)})

cat("\nFinal filtering and writing output files...")
write.table(geno_SNP, "SNP_Genotypes_Freq.txt", sep="\t", quote=F, row.names = F, col.names = T)
write.table(geno_indel, "indel_Genotypes_Freq.txt", sep="\t", quote=F, row.names = F, col.names = T)

geno_SNP %>% filter(! if_any(everything(), ~  grepl("multiAlleleFreq", .))) -> geno_SNP_bi
geno_SNP %>% filter(if_any(everything(), ~  grepl("multiAlleleFreq", .))) -> geno_SNP_multi
write.table(geno_SNP_bi, "SNP_Genotypes_Freq_Biallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)
write.table(geno_SNP_multi, "SNP_Genotypes_Freq_Multiallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)

geno_indel %>% filter(! if_any(everything(), ~  grepl("multiAlleleFreq", .))) -> geno_indel_bi
geno_indel %>% filter(if_any(everything(), ~  grepl("multiAlleleFreq", .))) -> geno_indel_multi
write.table(geno_indel_bi, "indel_Genotypes_Freq_Biallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)
write.table(geno_indel_multi, "indel_Genotypes_Freq_Multiallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)
cat(" Done!\n\n")














