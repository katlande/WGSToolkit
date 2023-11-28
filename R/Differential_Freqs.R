#!/usr/bin/env Rscript

# Converts GenoFreq.r into differential SNP files
args = commandArgs(trailingOnly=TRUE)
# This is a nested script
"Looking for differential variants by frequency..."
"Loading libraries..."
suppressMessages(library("tidyr", warn.conflicts = F, quietly = T))
suppressMessages(library("dplyr", warn.conflicts = F, quietly = T))
suppressMessages(library("reshape2", warn.conflicts = F, quietly = T))
suppressMessages(library("stringr", warn.conflicts = F, quietly = T))
setwd(args[1])

# Args:
# (1) <- $PWD
# (2) <- Allele diff cutoff
# (3) <- multi or bi mode
# (4) <- input file
# (5) <- .txt with treatment sample names
# (6) <- .txt with control sample names
# (7) <- output name
# (8) <- treatment name
# (9) <- ctl name

trtNAME <- args[8]
ctlNAME <- args[9]

trt.names <- read.delim(args[5], header=F)
trt.names <- as.character(trt.names$V1)
ctl.names <- read.delim(args[6], header=F)
ctl.names <- as.character(ctl.names$V1)
FreqCut <- args[2]

"Querying for differential sites..."
if(args[3] == "BI"){
  read.delim(args[4], stringsAsFactors = F) -> geno
  
  trt_sample_cols <- which(colnames(geno)%in%unique(grep(paste(trt.names,collapse="|"), colnames(geno), value=TRUE)))
  ctl_sample_cols <- which(colnames(geno)%in%unique(grep(paste(ctl.names,collapse="|"), colnames(geno), value=TRUE)))
  trt.cols <- which(colnames(geno)%in%colnames(geno)[trt_sample_cols][grepl("REF_Freq", colnames(geno)[trt_sample_cols])])
  ctl.cols <- which(colnames(geno)%in%colnames(geno)[ctl_sample_cols][grepl("REF_Freq", colnames(geno)[ctl_sample_cols])])
  
  geno$TRT_MEAN <- rowMeans(geno[trt.cols]) # for biallelic, take the mean of the samples' reference alleles
  geno$CTL_MEAN <- rowMeans(geno[ctl.cols])
  
  geno$dif <- geno$TRT_MEAN-geno$CTL_MEAN
  geno$abs_dif <- abs(geno$dif)
  geno <- geno[geno$abs_dif > FreqCut,]
  
  genoCall <- function(x, col){
    ifelse(x[col] > 0.5, x['REF'], x['ALT'])
  }  # Apply to add genotype columns
  geno$TRT_genotype <- apply(geno, 1, genoCall,'TRT_MEAN')
  geno$CTL_genotype <- apply(geno, 1, genoCall,'CTL_MEAN')
  
  geno <- geno[c(1:4,(ncol(geno)-5),(ncol(geno)-4),(ncol(geno)-2):ncol(geno))]
  colnames(geno) <- str_replace_all(colnames(geno), "TRT", trtNAME)
  colnames(geno) <- str_replace_all(colnames(geno), "CTL", ctlNAME)
  write.table(geno, args[7], quote=F, sep="\t", row.names = F, col.names = T) # This output is optimized for Call_Codons
  
} else if(args[3] == "MULTI"){
  read.delim(args[4], stringsAsFactors = F) -> multi
  multi.array <- function(multi){
    
    for(i in grep("Other_Freq",colnames(multi))){
      multi[i] <- str_remove_all(unlist(multi[i]), "multiAlleleFreq:")
      multi[i] <- str_replace_all(unlist(multi[i]), "biallelic", "0")
    }
    multi[c(5:ncol(multi))] <- apply(multi[c(5:ncol(multi))], 2, as.numeric)
    
    trt_sample_cols <- which(colnames(multi)%in%unique(grep(paste(trt.names,collapse="|"), colnames(multi), value=TRUE)))
    ctl_sample_cols <- which(colnames(multi)%in%unique(grep(paste(ctl.names,collapse="|"), colnames(multi), value=TRUE)))
    
    # TRT SAMPLES:
    trt.ref.cols <- which(colnames(multi)%in%colnames(multi)[trt_sample_cols][grepl("REF_Freq", colnames(multi)[trt_sample_cols])])
    trt.alt.cols <- which(colnames(multi)%in%colnames(multi)[trt_sample_cols][grepl("ALT_Freq", colnames(multi)[trt_sample_cols])])
    trt.multi.cols <- which(colnames(multi)%in%colnames(multi)[trt_sample_cols][grepl("Other_Freq", colnames(multi)[trt_sample_cols])])
    multi$TRT_MEAN_REF <- rowMeans(multi[trt.ref.cols]) 
    multi$TRT_MEAN_ALT <- rowMeans(multi[trt.alt.cols]) 
    multi$TRT_MEAN_MULTI <- rowMeans(multi[trt.multi.cols]) 
    
    # CTL SAMPLES:
    ctl.ref.cols <- which(colnames(multi)%in%colnames(multi)[ctl_sample_cols][grepl("REF_Freq", colnames(multi)[ctl_sample_cols])])
    ctl.alt.cols <- which(colnames(multi)%in%colnames(multi)[ctl_sample_cols][grepl("ALT_Freq", colnames(multi)[ctl_sample_cols])])
    ctl.multi.cols <- which(colnames(multi)%in%colnames(multi)[ctl_sample_cols][grepl("Other_Freq", colnames(multi)[ctl_sample_cols])])
    multi$CTL_MEAN_REF <- rowMeans(multi[ctl.ref.cols]) 
    multi$CTL_MEAN_ALT <- rowMeans(multi[ctl.alt.cols]) 
    multi$CTL_MEAN_MULTI <- rowMeans(multi[ctl.multi.cols]) 
    
    multi$abs_dif <- apply(multi, 1, function(x){
      ((abs(as.numeric(x['CTL_MEAN_REF'])-as.numeric(x['TRT_MEAN_REF']))+abs(as.numeric(x['CTL_MEAN_ALT'])-as.numeric(x['TRT_MEAN_ALT']))+abs(as.numeric(x['CTL_MEAN_MULTI'])-as.numeric(x['TRT_MEAN_MULTI'])))/2)
    })
    return(multi)
  } # Allele difference array function
  multi <- multi.array(multi)
  multi <- multi[multi$abs_dif > FreqCut,]
  multi <- multi[c(1:4,(ncol(multi)-6):ncol(multi))]
  colnames(multi) <- str_replace_all(colnames(multi), "TRT", trtNAME)
  colnames(multi) <- str_replace_all(colnames(multi), "CTL", ctlNAME)
  write.table(multi, args[7], quote=F, sep="\t", row.names = F, col.names = T) # This output is not optimized for Call_Codons; Call_Codons cannot take a multiallelic input.
} else {
  "Error: Mode option not recognized. Options: 'BI' or 'MULTI'"
  quit()
}

"Done!"













