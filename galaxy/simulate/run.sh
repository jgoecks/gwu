#!/bin/bash

# John M. Gaspar
# June 2015

# Simulate amplicon reads.
# External software requirements:
#   - ipcress (from exonerate 2.4.0)

# input files
bed=$1            # BED file listing locations of primers
gen=$2            # reference genome (FASTA)
var=$3            # file listing variants to make in reads

# check input files
if [ ! -f $bed ]; then
  echo "Input primer BED file not found"
  exit -1
elif [ ! -f $gen ]; then
  echo "Input reference genome not found"
  exit -1
fi

# retrieve primer-target sequences
prim=primers.txt
if [ ! -f $prim ]; then
  echo "Retrieving primers"
  perl getPrimers.pl $bed $gen $prim
fi

# convert sequences to ipcress-friendly format
ipc=primers.ipc
min=70       # minimum amplicon length
max=180      # maximum amplicon length
perl formatPrimers.pl $prim $ipc $min $max

# run ipcress
echo "Running ipcress"
mis=7  # maximum number of allowed mismatches
tr1=temp.ipcout
ipcress -i $ipc -s $gen -m $mis -p 0 > $tr1

# filter ipcress output
echo "Filtering ipcress output"
ipcout=primers.ipcout
score=0.85   # minimum primer matching score
perl filterIpc.pl $tr1 $ipc $gen $ipcout $score
rm $tr1

# simulate reads
echo "Simulating reads"
out1=reads_1.fastq  # output FASTQ file #1
out2=reads_2.fastq  # output FASTQ file #2
num=1000     # number of PE reads for each perfect primer match score
len=100      # length of PE reads to create
if [ ! -f $var ]; then
  perl readSim.pl $ipcout $ipc $gen $out1 $out2 $num \
    $len $score $min,$max $bed
else
  pct=10       # percent of reads to make variants
  perl readSim.pl $ipcout $ipc $gen $out1 $out2 $num \
    $len $score $min,$max $bed $var $pct
fi
