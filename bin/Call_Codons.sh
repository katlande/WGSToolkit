#!/bin/bash
#conda activate /gpfs/tools/anaconda/envs/WGStools/

Help()
{
   # Display Help
   echo ""
   echo "Call_Codons finds codon and peptide differences between a set of samples. It has two functionalities:"
   echo "VCF mode and plaintext mode. VCF mode compares all variable sites in a supplied list of genes and takes"
   echo "a VCF as input. Plaintext mode accepts a txt file as input and queries all loci included in the file."
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|i|o|c|p|r|a|b|t|g]"
   echo ""
   echo "h     Help"
   echo ""
   echo "Required Inputs:"
   echo "i     Input file: a txt file with chromosome, locus, gene, genotype, and reference allele columns."
   echo "o     Output file name"
   echo "l     Gene input: a column index (first index=0)."
   echo "c     Index of chromosome column in input file, first index=0."
   echo "p     Index of position column in input file, first index=0."
   echo "r     Index of reference allele column in input file, first index=0."
   echo "a     Control sample column index(first index=0)."
   echo "b     Treatment sample column index (first index=0)."
   echo "t     Path to the WGSToolkit folder."
   echo "g     Genome Directory. Choose from a list of preset genomes or input the"
   echo "      path to a local genome directory made with 'WGStools_Genome_Setup_Local'"
   echo ""
   echo "Currently available genomes:"
   find /gpfs/genomes/KatWGS/ -maxdepth 1 -mindepth 1 -type d -not -name BWA -not -name WGS_PYTHON -not -name Indel_Goldsets | sed 's|^/gpfs/genomes/KatWGS/||'
   echo ""
}


mode='VCF'
genome="hg38"
# Read in arguments:
while getopts i:o:g:l:c:p:r:a:b:h flag
do
    case "${flag}" in
        i) inf=${OPTARG};;
        o) outf=${OPTARG};;
        g) genome=${OPTARG};;
        l) gene=${OPTARG};;
        c) chrom=${OPTARG};;
        p) pos=${OPTARG};;
        r) ref=${OPTARG};;
        a) parent=${OPTARG};;
        b) child=${OPTARG};;
        t) toolkitpath=${OPTARG};;
        h) Help; exit 2;;
        ?) Help; exit 2;;
    esac
done

if [ -d $toolkitpath"/Genomes/"$genome ]; then
	genome=$toolkitpath"/Genomes/"$genome
	if [ $genome == $toolkitpath"/Genomes/" ]; then
		echo "Error: No genome specified"
		exit
	fi
else
	genome=$genome
fi

if [ ! -d $genome ]; then
	echo "Error: Genome does not exist"
	exit
fi

basename=${genome##*/}

python $toolkitpath"/python/Call_Codons_plain.py" -t $inf -o $outf -n $gene -y $chrom -p $parent -c $child -l $pos -r $ref -g $genome --base $basename


