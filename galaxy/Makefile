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
TARGET_DIR=~/projects/rcc/Plate3
DATA_DIR=~/projects/gwu
PRIMERS=$(DATA_DIR)/test/rcc_primers.bed
# Reference genome fasta and bowtie2 indices.
FASTA=$(DATA_DIR)/data/hg19/hg19.fa
BOWTIE2_INDICES=$(DATA_DIR)/data/hg19/bowtie2/hg19

#
# Commands.
#

# Run FFPE analysis to identify variants shared b/t plasma and FFPE samples. This was used to provide results for Plate 7 datasets.
ffpe_reports:
	$(ACTIVATE_ENV) && bash -x ./run_ffpe_reports.sh

# Process FFPE samples to create GEMINI tumor-normal databases using correct tumor-normal pairs.
ffpe_samples:
	$(ACTIVATE_ENV) && bash -x ./run_plate7_ffpe.sh

# Query all databases for somatic mutations.
# NOTE: this was used to provide results for Plate 6 and 8 datasets.
# NOTE: cannot use tumor_lod because Freebayes often assigns highest genotype likelihood to reference when tumor allele frequency (AF) is low,
# which in turn causes low AF tumor variants to be filtered out if tumor_lod is used as a criterion.
query_all:
	$(ACTIVATE_ENV) && find . -maxdepth 2 -name tumor_normal.db | sort | parallel -j 1 \
	"echo {}; gemini query -q 'select chrom, start, ref, alt, gene, transcript, \
							   impact, tumor_lod, normal_lod, tumor_aaf, normal_aaf, \
							   gt_alt_depths, gt_ref_depths, gt_depths from variants \
							   where normal_lod > 3.5 and tumor_aaf >= 0.05 and tumor_mabq >= 30 and \
							   aaf_exac_all <= 0.01 and aaf_1kg_all <= 0.01 and aaf_esp_all <= 0.01 and \
							   (normal_aaf <= 0.001 OR tumor_aaf/normal_aaf > 5)' \
							   --gt-filter 'gt_alt_depths.tumor >= 15 AND gt_alt_depths.normal <= 10' {} \
							   | sort -k1,1 -k2,2n"

# Process all matched tumor (P = plasma)-normal (BC = buffy coat) pairs in the current
# directory to create GEMINI databases and do somatic variant detection. This command
# should be run after run_dir has completed successfully; it uses the BAMs created via run_dir
# as inputs.
run_paired_somatic: all env target_regions.bed
	$(ACTIVATE_ENV); parallel -j 1 --xapply "bash -x run_paired_somatic.sh {1} {2} {=1 s/-Plasma_.*$///g =} $(FASTA)" ::: *-Plasma_* ::: *BC_*

# Process all fastq files in a directory.
run_dir: all env target_primers.bed
	$(ACTIVATE_ENV) && bash -x ./run_dir.sh ~/projects/gwu/galaxy/target_primers.bed $(FASTA) $(BOWTIE2_INDICES) $(TARGET_DIR)

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
	# Fix GEMINI version b/c it is shared with other projects.
	conda create -n $(ENV_NAME) -y python pandas
	$(ACTIVATE_ENV); pip install bioblend; \
		conda install -c bioconda -y bowtie2 varscan vcflib datamash freebayes parallel vt snpeff gemini=0.18.3; \
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
