# check command-line arguments
if [ $# -ne 3 ]; then
  echo "Usage: `basename $0`  <tumor_dir>  <normal_dir> <output_dir>" 1>&2
  exit -1
fi

TUMOR_DIR=$1
NORMAL_DIR=$2
OUTPUT_DIR=$3

FASTA=/home/jgoecks/projects/gwu/data/hg19/hg19.fa

mkdir ${OUTPUT_DIR}

# set home directory
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ADD_READ_GROUPS_SCRIPT=${HOME_DIR}/add_read_group.sh

# Add read groups for tumor.
samtools view -h ${TUMOR_DIR}/${TUMOR_DIR}.final.bam | \
    ${ADD_READ_GROUPS_SCRIPT} tumor tumor Illumina tumor | \
    samtools view -hSb - > ${OUTPUT_DIR}/tumor.bam
samtools index ${OUTPUT_DIR}/tumor.bam

# Add read groups for normal.
samtools view -h ${NORMAL_DIR}/${NORMAL_DIR}.final.bam | \
    ${ADD_READ_GROUPS_SCRIPT} normal normal Illumina normal | \
    samtools view -hSb - > ${OUTPUT_DIR}/normal.bam
samtools index ${OUTPUT_DIR}/normal.bam

# Call variants using Freebayes.
PATH=${PATH}:.
${HOME_DIR}/fasta_generate_regions.py ${FASTA}.fai 10000000 > regions.txt; \
freebayes-parallel regions.txt 6 --bam ${OUTPUT_DIR}/tumor.bam --bam ${OUTPUT_DIR}/normal.bam \
    --fasta-reference ${FASTA} --theta "0.001" \
    --ploidy "2" -J -K -n "0" --haplotype-length "3" --min-repeat-size "5" \
    --min-repeat-entropy "1" -m "1" -q "0" -R "0" -Y "0" -e "1000" -F "0.05" -C "2" \
    --min-alternate-qsum "0" -G "1" --min-coverage "0" -a --base-quality-cap "0" \
    --prob-contamination "1e-08" --report-genotype-likelihood-max \
    -B "1000" --genotyping-max-banddepth "6" -W "1,3" -D "0.9" --genotype-qualities \
    > ${OUTPUT_DIR}/called.vcf

# Keep only variants in target regions.
vcfintersect -b target_regions.bed ${OUTPUT_DIR}/called.vcf > ${OUTPUT_DIR}/called_in_rois.vcf

# NOTE: the following steps are taken from cherry-analysis and should be updated as that repo is updated.

# Annotate with SOMATIC to create final set of variants.
# TODO: may need to set --tumor-thresh to -50 to pick up low-frequency variants.
python annotate_somatic.py ${OUTPUT_DIR}/called_in_rois.vcf > ${OUTPUT_DIR}/final.vcf

# Create GEMINI database.
${HOME_DIR}/create_gemini_db.sh ${OUTPUT_DIR}/final.vcf ${OUTPUT_DIR}/tumor_normal.db ${FASTA}
egrep '(^#)|SOMATIC' ${OUTPUT_DIR}/final.vcf > ${OUTPUT_DIR}/somatic_only.vcf
vt decompose -s ${OUTPUT_DIR}/somatic_only.vcf | vt normalize -r ${FASTA} - > ${OUTPUT_DIR}/somatic_only.norm.vcf
bgzip ${OUTPUT_DIR}/somatic_only.norm.vcf && tabix -p vcf ${OUTPUT_DIR}/somatic_only.norm.vcf.gz
gemini annotate -f ${OUTPUT_DIR}/somatic_only.norm.vcf.gz -a boolean -c SOMATIC ${OUTPUT_DIR}/tumor_normal.db

# Working query, but want to replace this with better query using annotated LODs and allele frequencies.
find . -name tumor_normal.db | sort | parallel -j 1 \
"echo {}; gemini query -q 'select chrom, start, ref, alt, gene, impact, SOMATIC, gt_alt_depths, gt_ref_depths, gt_depths from variants' --gt-filter 'gt_alt_depths.tumor >= 15 AND gt_alt_depths.normal <= 10' {}"
