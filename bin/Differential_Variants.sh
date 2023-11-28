#!/bin/bash

Help()
{
   # Display Help
   echo ""
   echo "Differential_Variants identifies variable loci that are substantially different between two"
   echo "samples or sets of samples. This script can only take fully diploid sample sets."
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|i|o|a|b|c|d|p|m|f|v]"
   echo ""
   echo "h     Help"
   echo "i     Input file as a tab-delimited text file. Optimized for the 'Genotypes' files from WGS_HardCall"
   echo "o     Output file name."
   echo "a     Group 1; a vector of input file column indices representing one group of samples to compare."
   echo "b     Group 2; second vector of column indices."
   echo "c     Optional name for group 1, otherwise will be named 'Group 1'"
   echo "d     Optional name for group 2, otherwise will be named 'Group 2'"
   echo "m     Maximum allowed missing samples per group between 0 & 1; auto=0.5"
   echo "f     Minimum genotype frequency for a consensus between 0.5 & 1; auto=0.75"
   echo "v     Mode: 'SNP' or 'INDEL'. Auto='SNP'"
   echo "p     Path to the WGSToolkit folder."
   echo ""
}

g1n="Group1"
g2n="Group2"
maxmiss=0.5
minfreq=0.75
mode="SNP"
# Read in arguments:
while getopts i:o:a:b:c:d:m:f:v:p:h flag
do
    case "${flag}" in
        i) inf=${OPTARG};;
        o) outf=${OPTARG};;
        a) g1=${OPTARG};;
        b) g2=${OPTARG};;
        c) g1n=${OPTARG};;
        d) g2n=${OPTARG};;
        m) maxmiss=${OPTARG};;
        f) minfreq=${OPTARG};;
        v) mode=${OPTARG};;
        p) toolkitpath=${OPTARG};;
        h) Help; exit 2;;
        ?) Help; exit 2;;
    esac
done

python $toolkitpath"/python/Differential_Variants.py" --input $inf --output $outf --group1 $g1 --group2 $g2 --name1 $g1n --name2 $g2n --max_missing $maxmiss --consensus_freq $minfreq -v $mode

