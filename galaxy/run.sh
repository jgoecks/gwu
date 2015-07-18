#!/bin/bash

# John M. Gaspar
# Dec. 2014

# Analyze a sample for variants.
# External software requirements:
#   - bowtie2 (2.2.3)
#   - samtools (0.1.19)
#   - VarScan (2.3.7)
#
# bowtie2 and samtools are assumed to be in PATH. Set VarScan location here:
VARSCAN="./VarScan.v2.3.7.jar"

# Arguments check.
if [ $# -ne "6" ]
then
  echo "Usage: `basename $0` <forward_reads> <reverse_reads> <primers_BED> <ref_genome_fasta> <ref_genome_bowtie2_prefix> <output_directory>"
  exit -1
fi


# input files
file1=$1          # Input FASTQ file #1
file2=$2          # Input FASTQ file #2
bed=$3            # BED file listing locations of primers
gen=$4            # reference genome (FASTA)
idx=$5            # bowtie2 index prefix (indexes will be generated if necessary)
dir=$6            # output directory

# Set home directory for script + analysis scripts.
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# check input files
if [[ ! -f $file1 || ! -f $file2 ]]; then
  echo "Input FASTQ files not found"
  exit -1
elif [ ! -f $bed ]; then
  echo "Input primer BED file not found"
  exit -1
elif [ ! -f $gen ]; then
  echo "Input reference genome not found"
  exit -1
fi

# retrieve primer-target sequences
prim=primers.txt
if [ ! -f $prim ]; then
  perl ${HOME_DIR}/getPrimers.pl $bed $gen $prim
fi

# stitch together reads
echo "Stitching reads"
tr1=join.fastq
tr2=un1.fastq
tr3=un2.fastq
stParam="-m 20 -p 0.1 -d"  # min overlap 20, 10% allowed mismatches, dovetailing
${HOME_DIR}/stitch -1 $file1 -2 $file2 -o $tr1 -u1 $tr2 -u2 $tr3 $stParam

# remove primers, with -rq
echo "Removing primers"
log1=joinlog.txt
tr0=join-pr.fastq
tr4=join-nopr.fastq
rpParam="-fp -1,1 -rp -1,1 -ef 2 -er 2"  # allowing 2 subs, can start at +/- 1
${HOME_DIR}/removePrimer -i $tr1 -p $prim -o $tr0 $rpParam -rq -l $log1 -w $tr4  # require both primers

# retrieve reads whose primers weren't found
echo "Getting failure reads"
tr5=nopr1.fastq
tr6=nopr2.fastq
perl ${HOME_DIR}/getReads.pl $tr4 $file1 $tr5
perl ${HOME_DIR}/getReads.pl $tr4 $file2 $tr6
# cat with unjoined reads
cat $tr2 >> $tr5
cat $tr3 >> $tr6

# removePrimer individually
echo "Removing primers individually"
tr7=nopr1-pr.fastq
tr8=nopr2-pr.fastq
log2=nopr1log.txt
log3=nopr2log.txt
rpParam2="-rl 16 -el 1 -b $bed -bp -1,1"  # more options to find second primer
${HOME_DIR}/removePrimer -i $tr5 -p $prim -o $tr7 $rpParam $rpParam2 -l $log2
${HOME_DIR}/removePrimer -i $tr6 -p $prim -o $tr8 $rpParam $rpParam2 -l $log3

# filter singletons
tr9=noprcomb.fastq
fsParam="-b -q -c"  # prefer both primers removed, higher quality read, no chimeras
perl ${HOME_DIR}/filterSingle.pl $tr7 $tr8 $tr9 $fsParam

# quality trim
echo "Quality filtering"
out1=joined.fastq
qtParam="-t 30 -n 20"  # min avg qual 30; min len 20; no window filtering
${HOME_DIR}/qualTrim -i $tr0 -o $out1 $qtParam
tr10=noprcomb-qt.fastq
${HOME_DIR}/qualTrim -i $tr9 -o $tr10 $qtParam

# cat joined and singletons
out2=combined.fastq
cat $out1 $tr10 > $out2

# check for bowtie2 indexes
if [[ ! -f $idx.1.bt2 || ! -f $idx.2.bt2 ||
    ! -f $idx.3.bt2 || ! -f $idx.4.bt2 ||
    ! -f $idx.rev.1.bt2 || ! -f $idx.rev.2.bt2 ]]; then
  echo "Building bowtie2 indexes"
  bowtie2-build $gen $idx
fi

# map with bowtie2
echo "Mapping with bowtie2"
out3=combined.sam
bwtParam="-D 200 -N 1 -L 18 -i S,1,0.50 -k 20"
proc=4   # number of processors
bowtie2 -x $idx -U $out2 -S $out3 $bwtParam -p $proc

# find length variants
echo "Finding length variants"
tr11=len1.txt
perl ${HOME_DIR}/findLengthVars.pl $out1 $bed $tr11
tr12=len2.fastq
perl ${HOME_DIR}/getLengthVars.pl $out1 $tr11 $tr12
len3=realign.txt
tr13=len3v.txt
perl ${HOME_DIR}/alignLengthVars.pl $tr12 $prim $bed $len3 $tr13 $gen

# check alternative mapping sites
out4=altMapping.txt
if [ ! -f $out4 ]; then
  echo "Checking alternative mapping sites"
  perl ${HOME_DIR}/checkAltMapping.pl $out2 $out3 $prim $bed $gen $out4
fi

# filter SAM -- multi-mapping and realignment of length variants
echo "Filtering SAM"
out5=combinedFiltered.sam
log4=realign.log
perl ${HOME_DIR}/filterSAM.pl $out2 $bed $out4 $out3 $out5 $len3 $log4

# convert SAM to sorted BAM
echo "Converting SAM to sorted BAM"
out6=combinedFiltered.bam
samtools view -b -S $out5 | samtools sort -f - $out6

# call variants
echo "Calling variants"
out7=combinedFiltered.pileup
samtools mpileup -B -Q 0 -d 100000 -f $gen $out6 > $out7
qual=30  # min quality score to count a base
tr15=combFil.snp
java -jar ${VARSCAN} pileup2snp --min-avg-qual $qual --min-coverage 0 --min-var-freq 0.01 --variants < $out7 > $tr15
tr16=combFil.indel
java -jar ${VARSCAN} pileup2indel --min-avg-qual $qual --min-coverage 0 --min-var-freq 0.01 --variants < $out7 > $tr16

# make VCF, filter
echo "Producing VCF"
tr17=temp.vcf
perl ${HOME_DIR}/makeVCF.pl $tr15 $tr16 $out7 $tr17 $qual
tr18=temp2.vcf
perl ${HOME_DIR}/filterVCF.pl $tr17 $bed $tr18
out8=combinedFiltered.vcf
perl ${HOME_DIR}/addHPtoVCF.pl $tr18 $gen $out8

# remove extra files
rm $tr0 $tr1 $tr2 $tr3 $tr4 $tr5 $tr6 $tr7 $tr8 $tr9 $tr10 \
  $tr11 $tr12 $tr13 $tr15 $tr16 $tr17 $tr18

# move files to output directory
if [ ! -d $dir ]; then
  mkdir $dir
fi
mv $out1 $out2 $out3 $out4 $out5 $out6 $out7 $out8 \
  $log1 $log2 $log3 $log4 $len3 $dir
