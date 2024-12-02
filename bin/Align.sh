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
   echo "This script aligns raw WGS fastqs. It is optimized for downstream variant calling."
   echo
   echo "Syntax: [-h|f|t|g|p|j]"
   echo ""
   echo "h     Help"
   echo "p     Path to the WGSToolkit folder."
   echo "j     Path to picard.jar."
   echo "f     Directory containing all input fastq files. All fastq in this directory will be aligned."
   echo "t     Sequencing type; 's' for single-end, 'p' for paired-end."
   echo "g     Reference genome. May input a name of a pre-set genome or input a path to a pre-made genome directory."
   echo "      See 'WGStools_Genome_Setup' or 'WGStools_Genome_Setup_Local' for details."
   echo ""
}

while getopts f:g:t:h flag
do
    case "${flag}" in
        p) toolkitpath=${OPTARG};;
        j) picard=${OPTARG};;
        f) fastq_dir=${OPTARG};;
        g) genome=${OPTARG};;
        t) type=${OPTARG};;
        h) Help; exit 2;;
        ?) Help; exit 2;;
    esac
done


# Read in the genome:
if [ -d $toolkitpath"/Genomes/BWA/"$genome ]; then
	genome=$toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa"
else
	basename=${genome##*/}
	genome=$genome"/BWA/"$basename".fa"
fi

# Some sanity checks before we dive in:
if [ ! -f $genome ]; then
	echo "Error: Genome directory does not exist!"
	exit
fi

if [ ! -d $fastq_dir ]; then
	echo "Error: fastq directory does not exist!"
	exit
fi


echo "Setting up alignment for variant calling pipeline..."
if [ ! -d tmp_fastq ]; then #if the tmp_fastq is reserved from a fail run no need to copy again
	mkdir tmp_fastq
	scp $fastq_dir/*fastq* ./tmp_fastq # Copy fastq to tmpdir
fi

echo "Aligning samples..."
# Align fastq
mkdir bam
# Calls BWA index in WGStools Genomes OR manual input:
if [ $type == "s" ] ; then
	for fq in ./tmp_fastq/*
	do
		string=${fq//".fastq.gz"}
		string=${string//".fastq"}
		string=${string//"./tmp_fastq/"}
		outBAM="./bam/"$string"_SORTED.bam"
		echo "aligning "$string"..."
		bwa mem $genome $fq -t 32 -v 1 | samtools sort -@32 -o $outBAM -
		echo "Sorting BAM file..."
		samtools index $outBAM -@ 32
	done
elif [ $type == "p" ] ; then
	for fq in ./tmp_fastq/*R1*
	do
		rev=${fq/R1/"R2"}
		string=${fq//".fastq.gz"}
		string=${string//".fastq"}
		string=${string//"./tmp_fastq/"}
		outBAM="./bam/"$string"_SORTED.bam"
		echo "aligning "${string//"R1_"}"..."
		bwa mem $genome $fq $rev -t 32 -v 1 | samtools sort -@32 -o $outBAM -
		echo "Sorting BAM file..."
		samtools index $outBAM -@ 32
	done
else
	echo "type option not specified; must be 's' (single-end) or 'p' (paired-end)"
	exit
fi
# if paired... if single...
# replace R1 with R2, str replace regex:
	#firstString="I love Suzi and Marry"
	#secondString="Sara"
	#echo "${firstString/Suzi/"$secondString"}"    
	# prints 'I love Sara and Marry'

rm -r ./tmp_fastq #remove tmp dir

echo "Extracting sample coverage..."
# Make coverage file:
multiBamSummary bins -bs 100000 --bamfiles ./bam/*bam -p 16 --outRawCounts coverageALL_100kb.txt 

echo "Marking duplicate reads..."
mkdir MarkDup
mkdir ./MarkDup/tmp
for bam in ./bam/*
do
	string=${bam//"SORTED.bam"}
	string=${string//"./bam/"}
	echo $string"..."
	string="./MarkDup/"$string
	output=$string"MarkDup.bam"
	met=$string"DupMetrics.txt"
	java -jar $picard"/picard.jar" MarkDuplicates \
	      I=$bam \
	      O=$output \
	      M=$met \
	      MAX_RECORDS_IN_RAM=2000000 \
	      TMP_DIR=./MarkDup/tmp
done

echo "Adding read groups..."
for bam2 in ./MarkDup/*MarkDup.bam
do
	string=${bam2//"_MarkDup.bam"}
	echo ${string//"./MarkDup/"}"..."
	string2=${string//"./MarkDup/"}
	output=$string"_MarkDupRG.bam"
	java -jar $picard"/picard.jar" AddOrReplaceReadGroups \
	      I= $bam2 \
	      O= $output \
	      RGID=1 \
	      RGLB=lib1 \
	      RGPL=tmp \
	      RGPU= ${string2: -2} \
	      RGSM= $string2
done

rm -r ./MarkDup/tmp

for file in ./MarkDup/*_MarkDupRG.bam
do
	echo "Indexing "$file
	samtools index $file -@ 32
done
