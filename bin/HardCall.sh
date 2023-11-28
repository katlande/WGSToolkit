#!/bin/bash
# From raw fastq, this script makes:
# - Bam files
# - Coverage files
# - SNP VCFs
# - indel VCFs

Help()
{
   # Display Help
   echo ""
   echo "This script calls variants from aligned WGS data. It is optimized for downstream variant calling."
   echo "It works in tandem with 'WGS_Align' and must be run from thae same working directory from which WGS_Align was run."
   echo "Separately aligned bams can be used if they are in a directory called 'MarkDup' in the current working directory"
   echo "and all have the extension *MarkDupRG.bam. These bam files must be indexed, sorted, duplicate marked, and"
   echo "must have read groups."
   echo
   echo "Syntax: [-h|p|r|n|t|m|g]"
   echo ""
   echo "h     Help"
   echo "p     Path to the WGSToolkit folder."
   echo "r     Path to R."
   echo "n     Sample ploidy as integer. E.g.., for diploid '-p 2'."
   echo "t     Filtering threshold; 'S' for strict filtering, 'L' for loose filtering."
   echo "m     Mixed ploidies. If 'True', script handles samples with different ploidies. If True,"
   echo "      Raw gvcfs must be supplied by placing them in a folder named 'gvcf' within the current working directory."
   echo "g     Reference genome. May input a name of a pre-set genome or input a path to a local genome directory."
   echo "      See 'WGStools_Genome_Setup' or 'WGStools_Genome_Setup_Local' for details."
   echo ""
}

# Read in arguments:
while getopts g:p:t:m:h flag
do
    case "${flag}" in
        p) toolkitpath=${OPTARG};;
        p) rpath=${OPTARG};;
        g) genome=${OPTARG};;
        n) ploid=${OPTARG};;
        t) thresh=${OPTARG};;
        m) mix=${OPTARG};;
        h) Help; exit 2;;
        ?) Help; exit 2;;
    esac
done

if [ -d $toolkitpath"/Genomes/BWA/"$genome ]; then
	genome_store=$toolkitpath"/Genomes/"$genome
	genome=$toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa"
else
	genome_store=$genome
	basename=${genome##*/}
	genome=$genome"/BWA/"$basename".fa" #LOCAL INPUT
fi

# Some sanity checks before we dive in:
if [ ! -f $genome ]; then
	echo "Error: Genome directory or fasta does not exist!"
	exit
fi

if [ $thresh == "L" ]; then
	echo "Loose filtering selected..."
	MM=0.3
	MM_b=0.24
	MQ=30
	MQ_b=28
elif [ $thresh == "S" ]; then
	echo "Strict filtering selected..."
	MM=0.5
	MM_b=0.4
	MQ=32
	MQ_b=30
else 
	echo "Filtering stringency option not recognized!"
	exit
fi

# if ploidies are not mixed:
if [ ! $mix == 'True' ]; then
	if [ -d "/gpfs/genomes/KatWGS/Indel_Goldsets/KnownIndels_"$genome".vcf.gz" ]; then
		known_sites="/gpfs/genomes/KatWGS/Indel_Goldsets/KnownIndels_"$genome".vcf.gz"
		echo "Running BQSR..."
		mkdir BQSR
		for bam in ./MarkDup/*MarkDupRG.bam
		do
		string=${bam//"_MarkDupRG.bam"}
		string=${string//"./MarkDup/"}
		echo $string"..."
		report="./BQSR/"$string".BQSRreport"
		out="./BQSR/"$string".BSQR.bam"
		# run BQSR:
		gatk BaseRecalibrator \
		-I $bam \
		--reference $genome/*.fa \
		--known-sites $known_sites \
		-O $report
		# Apply BQSR recalibration: 
		gatk ApplyBQSR \
		-I $bam \
		--reference $genome/*.fa \
		--bqsr-recal-file $report \
		-O $out
		done
	fi

	echo "Calling haplotypes..."
	mkdir gvcf
	if [ ! -d ./BQSR ]; then
		for bam in ./MarkDup/*MarkDupRG.bam
		do
			string=${bam//"_MarkDupRG.bam"}
			string=${string//"./MarkDup/"}
			echo $string"..."
			out="./gvcf/"$string".g.vcf.gz"
			gatk --java-options "-Xmx4g" HaplotypeCaller \
			-R $genome \
			--sample-ploidy $ploid \
			-I $bam \
			-O $out \
			-ERC GVCF
		done
	else 
		for bam in ./BQSR/*bam
		do
			string=${bam//".BSQR.bam"}
			string=${string//"./BQSR/"}
			echo $string"..."
			out="./gvcf/"$string".g.vcf.gz"
			gatk --java-options "-Xmx4g" HaplotypeCaller \
			-R $genome \
			--sample-ploidy $ploid \
			-I $bam \
			-O $out \
			-ERC GVCF
		done
	fi
fi

# Indel realignment is depreciated if we are using haplotype caller...
# If no goldset is available we can skip BQSR and hard filter more aggressively...

# combine gvcfs:
echo "Combining GVCFs..."
find ./gvcf -type f -name "*.vcf.gz" > tmp.list
gatk --java-options "-Xmx7g" CombineGVCFs \
  -R $genome \
  -V tmp.list \
  -O ./gvcf/All.Samples.g.vcf.gz

# call variants:
# Need to add the --sample-ploidy / -ploidy option in GenotypeGVCFs if included in HaplotypeCaller
echo "Genotyping GVCFs"
gatk --java-options "-Xmx7g" GenotypeGVCFs \
  -R $genome \
  -V ./gvcf/All.Samples.g.vcf.gz \
  -O All.Samples.vcf.gz


if [ ! $mix == 'True' ] && [ $ploid -eq 2 ]; then 
	# Filtering
	echo "Extracting biallelic variants..."
    gatk SelectVariants \
        -R $genome \
        -V All.Samples.vcf.gz \
        -restrict-alleles-to BIALLELIC \
        -O All.Samples.Biallelic.vcf.gz
	
	echo "Simple variants detected..."
	echo "Hard filtering variants with vcftools..."
	if [ -d ./BQSR ]; then
		echo "Extracting SNPs..."
		vcftools --gzvcf All.Samples.vcf.gz --max-missing $MM_b --mac 2 --min-meanDP 15 --minQ $MQ_b --remove-indels --recode --recode-INFO-all --out Filtered_SNPs
		echo "Extracting indels..."
		vcftools --gzvcf All.Samples.vcf.gz --max-missing $MM_b --mac 2 --min-meanDP 15 --minQ $MQ_b --keep-only-indels --recode --recode-INFO-all --out Filtered_Indels
		echo "Filtering biallelic variants..."
		vcftools --gzvcf All.Samples.Biallelic.vcf.gz --max-missing $MM_b --mac 2 --min-meanDP 15 --minQ $MQ_b --remove-indels --recode --recode-INFO-all --out Filtered_SNPs_Biallelic
		vcftools --gzvcf All.Samples.Biallelic.vcf.gz --max-missing $MM_b --mac 2 --min-meanDP 15 --minQ $MQ_b --keep-only-indels --recode --recode-INFO-all --out Filtered_Indels_Biallelic
	else
		echo "Extracting SNPs..."
		vcftools --gzvcf All.Samples.vcf.gz --max-missing $MM --mac 2 --min-meanDP 15 --minQ $MQ --remove-indels --recode --recode-INFO-all --out Filtered_SNPs
		echo "Extracting indels..."
		vcftools --gzvcf All.Samples.vcf.gz --max-missing $MM --mac 2 --min-meanDP 15 --minQ $MQ --keep-only-indels --recode --recode-INFO-all --out Filtered_Indels
		echo "Filtering biallelic variants..."
		vcftools --gzvcf All.Samples.Biallelic.vcf.gz --max-missing $MM --mac 2 --min-meanDP 15 --minQ $MQ --remove-indels --recode --recode-INFO-all --out Filtered_SNPs_Biallelic
		vcftools --gzvcf All.Samples.Biallelic.vcf.gz --max-missing $MM --mac 2 --min-meanDP 15 --minQ $MQ --keep-only-indels --recode --recode-INFO-all --out Filtered_Indels_Biallelic
	fi
	
	echo "Making simplified genotype file..."
	$rpath"/Rscript" --vanilla $toolkitpath"/R/GenoDiploid.r" $PWD $genome_store
	echo "Making frequency-based genotype file..."
	$rpath"/Rscript" --vanilla $toolkitpath"/R/GenoFreq.r" $PWD $genome_store
	
else
	echo "Complex variants detected..."
	echo "Hard filtering variants with bcftools..."
	bcftools view -i 'TYPE="snp" && QUAL>'$MQ' && DP>15 && AF > 0.02 & AF < 0.98' -O z -I All.Samples.vcf.gz -o Filtered_SNPs.vcf.gz
	bcftools view -i 'TYPE="indel" && QUAL>'$MQ' && DP>15 && AF > 0.02 & AF < 0.98' -O z -I All.Samples.vcf.gz -o Filtered_Indels.vcf.gz
	bcftools view -i 'TYPE="snp" && QUAL>'$MQ' && DP>15 && AF > 0.02 & AF < 0.98' -O z -I All.Samples.Biallelic.vcf.gz -o Filtered_SNPs_Biallelic.vcf.gz
	bcftools view -i 'TYPE="indel" && QUAL>'$MQ' && DP>15 && AF > 0.02 & AF < 0.98' -O z -I All.Samples.Biallelic.vcf.gz -o Filtered_Indels_Biallelic.vcf.gz	
	echo "Making frequency-based genotype file..."
	$rpath"/Rscript" --vanilla $toolkitpath"/R/GenoFreq.r" $PWD $genome_store
fi

rm tmp.list













