#!/bin/bash

Help()
{
   # Display Help
   echo ""
   echo "Binned_Enrichment bins a set of variants and looks to see if bins are enriched compared to a background set of variants."
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|i|o|a|n|b|g|p]"
   echo ""
   echo "h     Help"
   echo "i     Input file as a tab-delimited text file. Must have chromosome and locus columns named 'CHROM' and 'POS'"
   echo "a     All variant file as a tab-delimited text file. Must have chromosome and locus columns named 'CHROM' and 'POS'"
   echo "o     Output file name"
   echo "n     NormalizeBy, options = 'genome' or 'chrom', whether to normalize the bins by genome or by chromosome. Auto='chrom'"
   echo "b     Bin size in bp; auto=10000"
   echo "p     Path to the WGSToolkit folder."
   echo "g     Genome Directory. Choose from a list of preset genomes or input the"
   echo "      path to a local genome directory made with 'WGStools_Genome_Setup_Local'"
   echo ""
   echo "Currently available genomes:"
   find /gpfs/genomes/KatWGS/ -maxdepth 1 -mindepth 1 -type d -not -name BWA -not -name WGS_PYTHON -not -name Indel_Goldsets | sed 's|^/gpfs/genomes/KatWGS/||'
   echo ""
}

normby="chrom"
bin=10000
genome="hg38"
# Read in arguments:
while getopts i:o:a:n:g:b:h flag
do
    case "${flag}" in
        g) genome=${OPTARG};;
        i) inf=${OPTARG};;
        a) allf=${OPTARG};;
        n) normby=${OPTARG};;
        o) outf=${OPTARG};;
        b) bin=${OPTARG};;
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

python $toolkitpath"/python/binned_enrichment.py" -i $inf -o $outf -a $allf --NormalizeBy $normby --binsize $bin --genome $genome









