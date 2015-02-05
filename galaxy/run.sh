#!/bin/bash

# John M. Gaspar
# Dec. 2014

# Analyze a sample for variants.

# inputs
file1=          # Input FASTQ file #1
file2=          # Input FASTQ file #2
bed=            # BED file listing locations of primers
gen=            # reference genome (FASTA)
idx=            # bowtie2 indexes (will be generated if necessary)
dir=.           # output directory

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
  perl getPrimers.pl $bed $gen $prim
fi

# stitch together reads
echo "Joining reads"
overlap=20    # minimum overlap for joining
tr1=join.fastq
tr2=un1.fastq
tr3=un2.fastq
fastq-join -m $overlap $file1 $file2 -o %.fastq

# remove primers, with -rq
echo "Removing primers"
log1=joinlog.txt
out1=joined.fastq
tr4=join-nopr.fastq
rpParam="-fp -1,1 -rp -1,1 -ef 2 -er 3 -rl 16 -el 1 -b $bed -bp -1,1";
removePrimer -i $tr1 -p $prim -o $out1 $rpParam -rq -l $log1 -w $tr4

# retrieve reads whose primers weren't found
echo "Getting failure reads"
tr5=nopr1.fastq
tr6=nopr2.fastq
perl getReads.pl $tr4 $file1 $tr5
perl getReads.pl $tr4 $file2 $tr6
# cat with unjoined reads
cat $tr2 >> $tr5
cat $tr3 >> $tr6

# removePrimer individually
echo "Removing primers individually"
tr7=nopr1-pr.fastq
tr8=nopr2-pr.fastq
log2=nopr1log.txt
log3=nopr2log.txt
removePrimer -i $tr5 -p $prim -o $tr7 $rpParam -l $log2
removePrimer -i $tr6 -p $prim -o $tr8 $rpParam -l $log3

# get common reads and singletons
echo "Getting common reads"
out2=single.fastq
tr9=gc1.fastq
tr10=gc2.fastq
tr11=chim.fastq
perl getCommon.pl $tr7 $tr8 $tr9 $tr10 $out2 $tr11

# join matches
echo "Joining reads again"
tr12=gcjoin.fastq
tr13=gcun1.fastq
tr14=gcun2.fastq
fastq-join -m $overlap $tr9 $tr10 -o gc%.fastq

# cat joined together, and singletons together
cat $tr12 >> $out1
cat $tr13 $tr14 >> $out2

# cat joined and singletons
tr15=temp.fastq
cat $out1 $out2 > $tr15

# quality trim
echo "Quality filtering"
out3=comb.fastq
qtParam="-t 30 -n 20"  # min. avg. qual. of 30; min. len. 20; no window filtering
qualTrim -i $tr15 -o $out3 $qtParam

# check for bowtie2 indexes
if [[ ! -f $idx.1.bt2 || ! -f $idx.2.bt2 ||
    ! -f $idx.3.bt2 || ! -f $idx.4.bt2 ||
    ! -f $idx.rev.1.bt2 || ! -f $idx.rev.2.bt2 ]]; then
  echo "Building bowtie2 indexes"
  bowtie2-build $gen $idx
fi

# map with bowtie2
echo "Mapping with bowtie2"
out4=comb.sam
bwtParam="-D 200 -N 1 -L 18 -i S,1,0.50 -k 20"
proc=1   # number of processors
bowtie2 -x $idx -U $out3 -S $out4 $bwtParam -p $proc

# find length variants
echo "Finding length variants"
len1=len1.txt
perl findLengthVars.pl $out1 $bed $len1
len2=len2.fastq
perl getLengthVars.pl $out1 $len1 $len2
len3=len3.txt
len3v=len3v.txt
perl alignLengthVars.pl $len2 $prim $bed $len3 $len3v $gen

# check alternative mapping sites
out5=altMapping.txt
log4=altMapping.log
if [ ! -f $out5 ]; then
  echo "Checking alternative mapping sites"
  perl checkAltMapping.pl $out3 $out4 $prim $bed $gen $out5 $log4
fi

# filter SAM -- multi-mapping and realignment of length variants
echo "Filtering SAM"
out6=combFil.sam
log5=realign.log
perl filterSAM.pl $out3 $bed $out5 $out4 $out6 $len3 $log5

# convert SAM to sorted BAM
echo "Converting SAM to sorted BAM"
tr16=temp.bam
samtools view -b -S $out6 > $tr16
out7=combsort.bam
samtools sort $tr16 combsort

# call variants
echo "Calling variants"
tr17=temp.vcf
fbParam="-K -F 0.01 -m 0 -q 30"
freebayes $fbParam -f $gen $out7 -v $tr17";
tr18=temp2.vcf
vcfbreakmulti $tr17 > $tr18

# filter variants by region
tr19=temp3.vcf
perl filterVCF.pl $tr18 $bed $tr19

# filter variants by proximity to homopolymers
out8=varFil.vcf
log6=varFil.log
perl filterVars.pl $tr19 $gen $out8 $log6

# remove extra files
rm $tr1 $tr2 $tr3 $tr4 $tr5 $tr6 $tr7 $tr8 $tr9 $tr10 \
  $tr11 $tr12 $tr13 $tr14 $tr15 $tr16 $tr17 $tr18 $tr19

# move files to output directory
if [ ! -d $dir ]; then
  mkdir $dir
fi
mv $out1 $out2 $out3 $out4 $out5 $out6 $out7 $out8 \
  $log1 $log2 $log3 $log4 $log5 $log6 $len1 $len2 $len3 $len3v $dir
