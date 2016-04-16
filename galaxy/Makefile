# Sequence data to analyze:
FASTQ1=test/JPE10900-4296-0036-VHL-Control_S36_L001_R1_001.fastq
FASTQ2=test/JPE10900-4296-0036-VHL-Control_S36_L001_R2_001.fastq

# PCR primer locations:
PRIMERS=test/primers.bed

# Output directory:
OUTPUT_DIR=test_output

# Reference genome fasta and bowtie2 indices.
FASTA=~/projects/galaxy/data/hg19/hg19.fa
BOWTIE2_INDICES=~/projects/galaxy/data/hg19/bowtie2/hg19



# Run complete amplicon analysis for a single sample.
run: all
	source activate amplicon-analysis-env && \
	bash run.sh $(FASTQ1) $(FASTQ2) $(PRIMERS) $(FASTA) $(BOWTIE2_INDICES) $(OUTPUT_DIR)

# Compile programs and set up analysis environment. Creates empty target for record keeping.
all: programs env
	touch all

# Compile C programs. Creates empty target for record keeping.
programs: removePrimer qualTrim stitch
	touch programs

# Set up analysis environment. Creates empty target for record keeping.
env:
	conda create -n amplicon-analysis-env -y python
	source activate amplicon-analysis-env && \
		conda install -c bioconda -y bowtie2 samtools varscan
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
