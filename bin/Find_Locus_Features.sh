#!/bin/bash

Help()
{
   # Display Help
   echo ""
   echo "Find_Locus_Features identifies the genomic features at an input set of variable sites."
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|i|o|t|p|g]"
   echo ""
   echo "h     Help"
   echo "i     Input file as a tab-delimited text file."
   echo "o     Output file name."
   echo "t     Path to the WGSToolkit folder."
   echo "c     Chromosome column index in input data, first column=0, default=1"
   echo "p     Position column index in input data, first column=0, default=2."
   echo "g     Genome Directory. Choose from a list of preset genomes or input the"
   echo "      path to a local genome directory made with 'WGStools_Genome_Setup_Local'"
   echo ""
   echo "Currently available genomes:"
   find /gpfs/genomes/KatWGS/ -maxdepth 1 -mindepth 1 -type d -not -name BWA -not -name WGS_PYTHON -not -name Indel_Goldsets | sed 's|^/gpfs/genomes/KatWGS/||'
   echo ""
}

poscol=2
chromcol=1
genome="hg38"
# Read in arguments:
while getopts i:o:c:p:t:g:h flag
do
    case "${flag}" in
        g) genome=${OPTARG};;
        p) poscol=${OPTARG};;
        c) chromcol=${OPTARG};;
        t) toolkitpath=${OPTARG};;
        i) inf=${OPTARG};;
        o) outf=${OPTARG};;
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

python $toolkitpath"/python/Find_Locus_Features.py" -i $inf -o $outf -g $genome -c $chromcol --PosCol $poscol 









