#!/bin/bash
#conda activate /gpfs/tools/anaconda/envs/WGStools/

Help()
{
   # Display Help
   echo ""
   echo "Differential_Freqs looks for differential genotypes between sample sets from frequency-based genotype files."
   echo "This script is optimized for 'Genotypes_Freq' files output from WGS_HardCall"
   echo "Contact klande@salk.edu or igc@salk.edu for help."
   echo ""
   echo "Syntax: [-h|m|p|r|i|o|f|t|c|T|C]"
   echo ""
   echo "h     Help"
   echo "m     Mode: 'BI' or 'MULTI' for biallelic or multiallelic inputs respectively, auto='BI'"
   echo "p     Path to the WGSToolkit folder."
   echo "r     Path to R."
   echo "i     Input file: A Genotypes_Freq file from WGS_HardCall"
   echo "o     Output file name"
   echo "f     Allelic Frequency differential cutoff from 0.5 to 1. I.e., 1-% genotype frequency overlap between."
   echo "      input groups. 0 = 100% genotype overlap, 1 = 0% genotype overlap between groups. Auto=0.8"
   echo "t     Treatment samples: .txt file of sample names representing one group of samples to compare."
   echo "c     Control samples: .txt file of sample names representing another group of samples to compare."
   echo "T     Treatment sample name (optional)"
   echo "C     Control sample name (optional)"
   echo ""
}

mode='BI'
freq=0.8
tname="TRT"
cname="CTL"
# Read in arguments:
while getopts m:p:i:o:f:t:c:T:C:h flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
        r) rpath=${OPTARG};;
        p) toolkitpath=${OPTARG};;
        i) inf=${OPTARG};;
        o) outf=${OPTARG};;
        f) freq=${OPTARG};;
        t) trt=${OPTARG};;
        c) ctl=${OPTARG};;
        T) tname=${OPTARG};;
        C) cname=${OPTARG};;
        h) Help; exit 2;;
        ?) Help; exit 2;;
    esac
done

$rpath"/Rscript" --vanilla $toolkitpath"/R/Differential_Freqs.R" $PWD $freq $mode $inf $trt $ctl $outf $tname $cname

