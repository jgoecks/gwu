#!/bin/bash

#
# Script reads SAM from standard in and adds a read group header and read group tag to each read
# entry.
#
# Example usage:
#	samtools view -h input.bam | ./add_read_group.sh myID myLibrary myPlatform mySample | samtools view -hS - > output.bam
#
# This command produces output.bam, which includes the read group for all reads in input.bam
#

# Parameter testing and evaluation.
if [ $# -ne "4" ]
then
  echo "Usage: `basename $0` <id> <library> <platform> <sample>"
  exit -1
fi

# Set read group attributes.
ID=$1
LIBRARY=$2
PLATFORM=$3
SAMPLE=$4

# Output read group header.
echo -e "@RG\tID:$ID\tSM:$SAMPLE\tLB:$LIBRARY\tPL:$PLATFORM"

# Add read group to reads from standard in.
awk -v id=$ID '{ if (substr($1,1,1)=="@") print; else printf "%s\tRG:Z:%s\n",$0,id; }' /dev/stdin
