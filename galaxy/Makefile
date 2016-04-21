SHELL=bash

# Sequence data to analyze:
FASTQ1=../test/JPE10900-4296-0036-VHL-Control_S36_L001_R1_001.fastq
FASTQ2=../test/JPE10900-4296-0036-VHL-Control_S36_L001_R2_001.fastq

# PCR primer locations:
PRIMERS=../test/rcc_primers.bed

# Output directory:
OUTPUT_DIR=test_output

# Reference genome fasta and bowtie2 indices.
FASTA=../data/hg19/hg19.fa
BOWTIE2_INDICES=../data/hg19/bowtie2/hg19

# Run complete amplicon analysis for a single sample.
run: all target_regions.bed
	source activate amplicon-analysis-env && \
	bash -x run.sh $(FASTQ1) $(FASTQ2) $(PRIMERS) $(FASTA) $(BOWTIE2_INDICES) $(OUTPUT_DIR)

# Compile programs and set up analysis environment. Creates empty target for record keeping.
all: programs env
	touch all

# Create sequencing targets from primers, removing "bad" regions that have too high a repeat density.
target_regions.bed:
	python primers_to_target_regions.py < $(PRIMERS) | python remove_bad_regions.py > sort -k1,1 -k2,2 > target_regions.bed

# Compile C programs. Creates empty target for record keeping.
programs: removePrimer qualTrim stitch
	touch programs

# Set up analysis environment. Creates empty target for record keeping.
env:
	conda create -n amplicon-analysis-env -y python
	# Special installation necessary for samtools to link correctly on linux.
	source activate amplicon-analysis-env && \
		conda install -c bioconda -y bowtie2 varscan vcflib && \
		conda install -c r -y ncurses && conda install -c bioconda -y samtools
	touch env

# Compile removePrimer program.
removePrimer: removePrimer.c removePrimer.h
	gcc -g -Wall -O3 -std=c99 -o removePrimer removePrimer.c

# Compile qualTrim program.
qualTrim: qualTrim.c qualTrim.h
	gcc -g -Wall -O3 -std=c99 -o qualTrim qualTrim.c

# Compile stich program.
stitch: stitch.c stitch.h
	gcc -g -Wall -O3 -std=c99 -o stitch stitch.c

# Remove compiled programs and analysis environment.
clean:
	rm -f removePrimer qualTrim stitch all programs env
	conda remove -n amplicon-analysis-env -y --all
