SHELL=bash

#
# General usage:
# % make env
# % make run
#


# For setting up a conda environment.
ENV_NAME=rcc-env
ACTIVATE_ENV=source activate $(ENV_NAME)

#
# For running the analyses.
#
DATA_DIR=~/projects/gwu
PRIMERS=$(DATA_DIR)/test/rcc_primers.bed
# Reference genome fasta and bowtie2 indices.
FASTA=$(DATA_DIR)/data/hg19/hg19.fa
BOWTIE2_INDICES=$(DATA_DIR)/data/hg19/bowtie2/hg19

#
# Commands.
#

# Process all matched tumor (P = plasma)-normal (BC = buffy coat) pairs in the current
# directory to create GEMINI databases and do somatic variant detection. This command
# should be run after run_dir has completed successfully; it uses the BAMs created via run_dir
# as inputs.
run_paired_somatic: env
	$(ACTIVATE_ENV); parallel -j 1 --xapply "bash -x run_paired_somatic.sh {1} {2} {=1 s/-P_.*$///g =}" ::: *-P_* ::: *BC_*

# Process all fastq files in a directory.
run_dir: all env target_primers.bed
	$(ACTIVATE_ENV) && bash -x ./run_dir.sh ~/projects/gwu/galaxy/target_primers.bed $(FASTA) $(BOWTIE2_INDICES) ~/projects/rcc/Plate7

run: all env
	#$(ACTIVATE_ENV) && bash -x ./run.sh $(PRIMERS) $(FASTA) $(BOWTIE2_INDICES)
	$(ACTIVATE_ENV) && bash -x /home/jgoecks/projects/gwu/galaxy/run.sh \
						/home/jgoecks/projects/rcc/Plate6/JPE10900-12028-VA006-P_S24_L001_R1_001.fastq \
						/home/jgoecks/projects/rcc/Plate6/JPE10900-12028-VA006-P_S24_L001_R2_001.fastq \
						/home/jgoecks/projects/gwu/test/rcc_primers.bed \
						/home/jgoecks/projects/gwu/data/hg19/hg19.fa \
						/home/jgoecks/projects/gwu/data/hg19/bowtie2/hg19 \
						test

# Compile programs and set up analysis environment. Creates empty target for record keeping.
all: programs env
	touch all

# Create list of good primers.
target_primers.bed:
	python remove_bad_regions.py < $(PRIMERS) | sort -k1,1 -k2,2 > target_primers.bed

# Create sequencing targets from primers, removing "bad" regions that have too high a repeat density.
target_regions.bed:
	python primers_to_target_regions.py < $(PRIMERS) | python remove_bad_regions.py | sort -k1,1 -k2,2 > target_regions.bed

# Set up conda virtual environment with needed packages and tools.
env:
	# Special installation necessary for samtools to link correctly on linux.
	conda create -n $(ENV_NAME) -y python
	$(ACTIVATE_ENV); pip install bioblend; \
		conda install -c bioconda -y bowtie2 varscan vcflib datamash freebayes parallel vt snpeff gemini; \
		conda install -c r -y ncurses && conda install -c bioconda -y samtools; \
	touch env

## Remove environment.
clean-env:
	conda remove --name $(ENV_NAME) --all -y
	rm env

# Compile C programs. Creates empty target for record keeping.
programs: removePrimer qualTrim stitch
	touch programs

# Compile removePrimer program.
removePrimer: removePrimer.c removePrimer.h
	gcc -g -Wall -O3 -std=c99 -o removePrimer removePrimer.c

# Compile qualTrim program.
qualTrim: qualTrim.c qualTrim.h
	gcc -g -Wall -O3 -std=c99 -o qualTrim qualTrim.c

# Compile stich program.
stitch: stitch.c stitch.h
	gcc -g -Wall -O3 -std=c99 -o stitch stitch.c
