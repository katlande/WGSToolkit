#!/bin/bash

Help()
{
   # Display Help
   echo ""
   echo "Feature_Enrichment identifies if any features are enriched in a set of input variable sites"
   echo "compared to the genomic background."
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|i|o|p|g]"
   echo ""
   echo "h     Help"
   echo "i     Input file as a tab-delimited text file. Script is designed to take the output from Find_Locus_Features"
   echo "o     Output file name."
   echo "p     Path to the WGSToolkit folder."
   echo "g     Genome Directory. Choose from a list of preset genomes or input the"
   echo "      path to a local genome directory made with 'WGStools_Genome_Setup_Local'; default='hg38'"
   echo ""
   echo "Currently available genomes:"
   find /gpfs/genomes/KatWGS/ -maxdepth 1 -mindepth 1 -type d -not -name BWA -not -name WGS_PYTHON -not -name Indel_Goldsets | sed 's|^/gpfs/genomes/KatWGS/||'
   echo ""
}

genome="hg38"
# Read in arguments:
while getopts i:o:c:p:g:h flag
do
    case "${flag}" in
        g) genome=${OPTARG};;
        i) inf=${OPTARG};;
        o) outf=${OPTARG};;
        p) toolkitpath=${OPTARG};;
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

python $toolkitpath"/python/Feature_Enrichment.py" -i $inf -o $outf --genome $genome

