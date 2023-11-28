args = commandArgs(trailingOnly=TRUE)
# This is a nested script
"Extracting genotype information to plain text..."
"Loading libraries..."
suppressMessages(library("tidyr", warn.conflicts = F, quietly = T))
suppressMessages(library("dplyr", warn.conflicts = F, quietly = T))
suppressMessages(library("reshape2", warn.conflicts = F, quietly = T))
suppressMessages(library("stringr", warn.conflicts = F, quietly = T))
setwd(args[1])

# Clean up sample names:
sample_names <- read.delim("tmp.list", header=F)$V1
sample_names <- str_remove_all(sample_names, "./gvcf/")
sample_names <- str_remove_all(sample_names, ".g.vcf.gz")
sample_names <- str_remove_all(sample_names, "_R1")
sample_names <- str_remove_all(sample_names, "_R2")
sample_names <- str_remove_all(sample_names, "_001")

# Rscript --vanilla /gpfs/genomes/KatWGS/GenoFreq.r $PWD

"Reading filtered VCFs..."
VCF_SNP <- read.delim("Filtered_SNPs.recode.vcf", comment.char = "#", header=F, stringsAsFactors = F)
geno_SNP <- VCF_SNP[c(1,2,4,5,10:ncol(VCF_SNP))]
#head(geno_SNP)
#sample_names
colnames(geno_SNP) <- c("CHROM",	"POS",	"REF",	"ALT", sample_names[!sample_names=="All.Samples"])
VCF_indel <- read.delim("Filtered_Indels.recode.vcf", comment.char = "#", header=F, stringsAsFactors = F)
geno_indel <- VCF_indel[c(1,2,4,5,10:ncol(VCF_indel))]
colnames(geno_indel) <- c("CHROM",	"POS",	"REF",	"ALT", sample_names[!sample_names=="All.Samples"])

# make sure chromosome names start with "chr" - will cause issues later if they don't
geno_SNP$CHROM <- paste0("chr",geno_SNP$CHROM)
geno_SNP$CHROM <- str_replace_all(geno_SNP$CHROM, "chrchr", "chr")
geno_indel$CHROM <- paste0("chr",geno_indel$CHROM)
geno_indel$CHROM <- str_replace_all(geno_indel$CHROM, "chrchr", "chr")

t<-paste0(args[2],"/region_sizes.txt")
if(file.exists(t)){
  "Filtering chromosomes..."}
if(file.exists(t)){
  usefulChrom <- read.delim(t, stringsAsFactors = F)
  usefulChrom <- usefulChrom$Region
  geno_SNP <- geno_SNP[geno_SNP$CHROM %in% usefulChrom,]
  geno_indel <- geno_indel[geno_indel$CHROM %in% usefulChrom,] 
  usefulChrom
}


# Remove non-genotype information:
"Extracting Genotypes..."
geno_SNP[c(5:ncol(geno_SNP))] <- apply(geno_SNP[c(5:ncol(geno_SNP))], 2, function(y) (gsub("\\:.*", "", y)))
geno_indel[c(5:ncol(geno_indel))] <- apply(geno_indel[c(5:ncol(geno_indel))], 2, function(y) (gsub("\\:.*", "", y)))

"Splitting SNP Alleles..."
geno_SNP %>% separate(ALT, c("ALT_1", "ALT_2", "ALT_3"), sep=",", extra = "drop") -> geno_SNP
"Splitting indel Alleles..."
geno_indel %>% separate(ALT, c("ALT_1", "ALT_2", "ALT_3"), sep=",", extra = "drop") -> geno_indel

"*Note: 'Expected 3 pieces' warning message is expected."

# write output files:
"Writing output files..."
write.table(geno_SNP, "SNP_Genotypes.txt", sep="\t", quote=F, row.names = F, col.names = T)
write.table(geno_indel, "indel_Genotypes.txt", sep="\t", quote=F, row.names = F, col.names = T)

geno_SNP[is.na(geno_SNP$ALT_2) & is.na(geno_SNP$ALT_3),] -> geno_SNP_bi
geno_SNP[! is.na(geno_SNP$ALT_2) | ! is.na(geno_SNP$ALT_3),] -> geno_SNP_multi
write.table(geno_SNP_bi, "SNP_Genotypes_Biallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)
write.table(geno_SNP_multi, "SNP_Genotypes_Multiallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)

geno_indel[is.na(geno_indel$ALT_2) & is.na(geno_indel$ALT_3),] -> geno_indel_bi
geno_indel[! is.na(geno_indel$ALT_2) | ! is.na(geno_indel$ALT_3),] -> geno_indel_multi
write.table(geno_indel_bi, "indel_Genotypes_Biallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)
write.table(geno_indel_multi, "indel_Genotypes_Multiallelic.txt", sep="\t", quote=F, row.names = F, col.names = T)

