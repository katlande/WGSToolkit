#!/bin/bash
# This script takes as its input a GTF and an output directory and formats the genome for use with the WGStoolkit

Help()
{
   # Display Help
   echo ""
   echo "This script sets up a genome for WGStools in the main WGStools path using a genome fasta and GTF."
   echo ""
   echo "Syntax: [-h|f|g|c|p|j|m]"
   echo ""
   echo "h     Help"
   echo "f     Genome fasta. Genomes from http://ftp.ensembl.org/pub/ are the best for this pipeline."
   echo "g     Genome name. E.g., 'mm10' or 'hg38'. This is the name all scripts will use to access this genome."
   echo "b     Genome GTF. Genomes from http://ftp.ensembl.org/pub/ are the best for this pipeline."
   echo "c     Chromosomes to include in feature parsing. Generally suggested to include all autosomes and sex chromosomes. Format example for mouse: '1-20,X,Y'"
   echo "p     Path to the WGSToolkit folder."
   echo "j     Path to picard.jar."
   echo "m     Mode. This argument can generally be ignored. If set to 'l', WGStools will partially overwrite extant genomes with the same name."
   echo "      This option is intended ONLY for use when a genome setup run fails halfway through as it skip some memory-intensive steps if they are already complete."
}

# Read in arguments:
while getopts b:f:g:c:m: flag
do
    case "${flag}" in
        b) gtf=${OPTARG};;
        f) fasta=${OPTARG};;
        g) genome=${OPTARG};;
        c) chroms=${OPTARG};;
        p) toolkitpath=${OPTARG};;
        j) picard=${OPTARG};;
        m) mode=${OPTARG};;
        h) Help; exit 2;;
        ?) Help; exit 2;;
    esac
done

# make a directory in which to store the genome:
# check if dir exists:
outdir= $toolkitpath"/Genomes/"$genome
if [ -d $outdir ] ; then
	
	if [ ! $mode == "l" ] ; then
	    echo -e "\n" $genome "already exists"
		echo -e "use '"$genome"' or give this genome a different version name\n"
		exit
	fi
else
	mkdir $outdir 
fi

scp $fasta $outdir"/"$genome".fa"
scp $gtf $outdir"/"$genome".gtf"

# make a BWA index for the genome if one is not present:
if [ ! -d $toolkitpath"/Genomes/BWA/"$genome ]; then
	echo "Making BWA index"
	mkdir $toolkitpath"/Genomes/BWA/"$genome
	scp $fasta $toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa"
	bwa index $toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa"
else
	echo "BWA index detected. Skipping index generation..."
fi

if [ ! -f $toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa.fai" ]; then
	echo "Indexing genome fasta..."
	samtools faidx $toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa"
else
	echo "Fasta index detected. Skipping genome indexing..."
fi

if [ ! -f $toolkitpath"/Genomes/BWA/"$genome"/"$genome".dict" ]; then
	echo "Making genome fasta dict..."
	java -jar $picard"/picard.jar" CreateSequenceDictionary \
	      R=$toolkitpath"/Genomes/BWA/"$genome"/"$genome".fa" \
	      O=$toolkitpath"/Genomes/BWA/"$genome"/"$genome".dict"
else
	echo "Fasta dictionary detected. Skipping dictionary generation..."
fi

# Run gtftools
echo -e "\n#############################################################################################\n This script formats genomes for use with WGStools scripts.\n Genomes from http://ftp.ensembl.org/pub/ are the best for this pipeline.\n For genomes not supported by ensembl, errors may occur. Contact klande@salk.edu for support.\n#############################################################################################\n"
echo -e "Running modified gtftools.py, developed by Hongdong Li at Central South University, Changsha, P.R. China...\n"
echo "Extracting Exons..."
python $toolkitpath"/python/gtftools.py" $gtf  -e $outdir"/"$genome".exons.bed" -c $chroms
echo "Extracting Genes..."
python $toolkitpath"/python/gtftools.py" -g $outdir"/"$genome".genes.bed" $gtf -c $chroms
echo "Extracting Intergenic Regions..."
python $toolkitpath"/python/gtftools.py" -b $outdir"/"$genome".intergenic.bed" $gtf -c $chroms
echo "Extracting UTRs..."
python $toolkitpath"/python/gtftools.py" -u $outdir"/"$genome".UTR.bed" $gtf -c $chroms
echo "Extracting regions 500bp upstream of TSSs..."
python $toolkitpath"/python/gtftools.py" -t $outdir"/"$genome".TSS_500.bed" $gtf -c $chroms -w 500-0
echo "Extracting Masked Introns..."
python $toolkitpath"/python/gtftools.py" -k $outdir"/"$genome".masked_introns.bed" $gtf -c $chroms
echo "Extracting Unmasked Introns (this may take a while)..."
python $toolkitpath"/python/gtftools.py" -d $outdir"/"$genome".introns.bed" $gtf -c $chroms
echo "Extracting Splice Sites..."
python $toolkitpath"/python/gtftools.py" -q $outdir"/"$genome".splice_site.bed" $gtf -c $chroms
echo "Genomic Feature Creation Complete"

# make the guide file:
echo -e "gene\t"$genome".genes.bed" > $outdir"/guide.txt"
echo -e "exon\t"$genome".exons.bed" >> $outdir"/guide.txt"
echo -e "intron_true\t"$genome".introns.bed" >> $outdir"/guide.txt"
echo -e "intron_masked\t"$genome".masked_introns.bed" >> $outdir"/guide.txt"
echo -e "intergenic\t"$genome".intergenic.bed" >> $outdir"/guide.txt"
echo -e "TSS\t"$genome".TSS_500.bed" >> $outdir"/guide.txt"
echo -e "UTR\t"$genome".UTR.bed" >> $outdir"/guide.txt"
echo -e "splice_site\t"$genome".splice_site.bed" >> $outdir"/guide.txt"

retdir=$PWD; cd $outdir
# parse gtftools files by chromosome:
echo "Parsing Feature Files by chromosome..."
python $toolkitpath"/python/WGStools_Genome_Files.py" --name $genome
# extract CDS indicies for
mkdir $outdir"/"$genome"/CDS"
echo "Formatting CDS Sites for WGSTools..."
python $toolkitpath"/python/Format_GTF.py" --name $genome --gtf $genome".gtf"
find $outdir"/"$genome"/CDS/" -mindepth 1 -type d -empty -delete
cd $retdir

chmod -R 777 $outdir
chmod -R 777 $toolkitpath"/Genomes/BWA/"$genome

