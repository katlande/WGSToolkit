#!/bin/bash

Help()
{
   # Display Help
   echo ""
   echo "CDS_Checker returns whether a list of loci are in CDS regions."
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|i|o|c|p|n|t|g]"
   echo ""
   echo "h     Help"
   echo "i     Input file as a tab-delimited text file. Must have chromosome, locus, and gene columns."
   echo "o     Output file name"
   echo "c     Index of chromosome column in input file, first index=0. Auto=0"
   echo "p     Index of position column in input file, first index=0. Auto=1"
   echo "n     index of gene column in input file, first index=0, Auto=2"
   echo "t     Path to the WGSToolkit folder."
   echo "g     Genome Directory. Choose from a list of preset genomes or input the"
   echo "      path to a local genome directory made with 'WGStools_Genome_Setup_Local'"
   echo ""
   echo "Currently available genomes:"
   find /gpfs/genomes/KatWGS/ -maxdepth 1 -mindepth 1 -type d -not -name BWA -not -name WGS_PYTHON -not -name Indel_Goldsets | sed 's|^/gpfs/genomes/KatWGS/||'
   echo ""
}

chrom=0
pos=1
gene=2
genome="hg38"
# Read in arguments:
while getopts i:o:c:p:n:g:t:h flag
do
    case "${flag}" in
        g) genome=${OPTARG};;
        i) inf=${OPTARG};;
        o) outf=${OPTARG};;
        c) chrom=${OPTARG};;
        p) pos=${OPTARG};;
        n) gene=${OPTARG};;
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
echo $basename

python $toolkitpath"/python/CDS_Checker.py" -i $inf -o $outf -c $chrom -p $pos -n $gene --genome $genome --base $basename

