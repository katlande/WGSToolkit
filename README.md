### Overview

WGSToolkit is a set of command line wrapper scripts and original tools for short read WGS data. Wrapper scripts are provided for alignment with BWA and small variant (SNP and INDEL) calling with GATK. WGSToolkit works for any genome available on ensembl.

Most functionalities of this toolkit work for non-diploid samples.

### Dependencies
R-4.0.0≥, with the following packages:
- tidyr
- dplyr
- reshape2
- stringr

Python 3, with the following packages:
- pandas
- scipy.stats
- csv
- Fasta (pyfaidx)
- vcf
- SeqIO (Bio)

BWA
Picard
samtools
gatk
vcftools

Note: while the paths to R and picard can be specified, the scripts expect BWA, python, samtools, vcftools, and gatk to be in your $PATH.



### Set-up scripts
| Script | Description | Output |
|-------|--------|----------|
| Align.sh | Align is a wrapper script for BWA. It takes a set of FASTQs and aligns them to a reference genome with BWA. Resulting bams are sorted, indexed, and marked for PCR duplicates. | Sorted and duplicate marked bam files.|
| HardCall.sh | HardCall takes the Align output and runs GATK's haplotypecaller. Variants are hard filtered for quality. | SNP and indel VCFs for all variants and biallelic variants only. If all samples are diploid, Genotype and Genotype_Freq text files will be produced as well. If a mix of ploidies or non-diploid samples are input, Genotype_Freq (but not Genotype) text files will be produced. These files are used for downstream analyses.|
| WGStools_Genome_Setup.sh | Adds a new genome to the WGSToolkit genome presets using a .gtf and a .fa file from ensembl. In addition to making a BWA index, this script will also pull out genomic features and CDS regions, and make additional supporting files for downstream analysis. This step needs only to be run once per genome. | Adds a genome to the preset genomes in WGSToolkit.|


### Analysis scripts 
| Script | Description | Output |
|-------|--------|----------|
| Differential_Variants.sh | For identifying variants that differ between two samples or sets of samples; diploid only. Accepts SNPs and indels. Differential_Variants.sh offers better functionality for identifying heterozygous variants than Differential_Freqs.sh. | tab-delimited txt file |
| Differential_Freqs.sh | For identifying variants that differ between two samples or sets of samples; accepts non-diploid and mixed ploidy samples, SNPs, and indels. Differential_Variants.sh is recommended for diploid-only samples as it is better at detecting heterozygous variation. | tab-delimited txt file |
| Find_Locus_Features.sh | Annotates a set of variants with genomic loci. | tab-delimited txt file|
| binned_enrichment.sh | Checks genomic regions in bins to see if a subset of variants are enriched compared to all variants. | tab-delimited txt file |
| Feature_Enrichment.sh | Checks if a Find_Locus_Features.sh output file is enriched or depleted for genomic features. | tab-delimited txt file |
| CDS_Checker.sh | Checks if variants are inside the CDS region of a protein-coding gene. Works best with a Find_Locus_Features.sh output file. | tab-delimited txt file |
| ExtractPCVariants.sh | Preps differential variants file for Call_Codons.sh automatically by combining the Differential freq and Find_Locus_Feature outputs. | tab-delimited txt file|
| Call_Codons.sh | Identifies protein coding changes introduced by variants. For SNPs, identifies type (missense, stop, silent, nonstop); for indels, identifies whether or not a frameshift is introduced. Optimized to take a Differential_Freqs.sh output file, but can take any text file as long as the requisite columns are present. | tab-delimited txt file |



