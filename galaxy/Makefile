#
# General usage:
# % make env
# % make run
#

# For setting up a conda environment.
ENV_NAME=rcc-env
ACTIVATE_ENV=bash -c "source activate $(ENV_NAME)"

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

# Create sequencing targets from primers, removing "bad" regions that have too high a repeat density.
target_regions.bed:
	python primers_to_target_regions.py < $(PRIMERS) | python remove_bad_regions.py | sort -k1,1 -k2,2 > target_regions.bed

# Set up conda virtual environment with needed packages and tools.
env:
	# Special installation necessary for samtools to link correctly on linux.
	conda create -n $(ENV_NAME) -y python && \
		pip install bioblend && \
		conda install -c bioconda -y bowtie2 varscan vcflib && \
		conda install -c r -y ncurses && conda install -c bioconda -y samtools
	touch env

clean:
	conda remove --name $(ENV_NAME) --all -y

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
