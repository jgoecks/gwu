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

# Annotate with additional to create final set of variants.
python annotate_tumor_normal.py ${OUTPUT_DIR}/called_in_rois.vcf > ${OUTPUT_DIR}/final.vcf

# Create GEMINI database and annotate with somatic information.
# FIXME: create_gemini_db.sh actually creates final.norm.vcf.gz in the parent directory, but its difficult to use.
${HOME_DIR}/create_gemini_db.sh ${OUTPUT_DIR}/final.vcf ${OUTPUT_DIR}/tumor_normal.db ${FASTA}
vt decompose -s ${OUTPUT_DIR}/final.vcf | vt normalize -r ${FASTA} - > ${OUTPUT_DIR}/final.norm.vcf
bgzip ${OUTPUT_DIR}/final.norm.vcf && tabix -p vcf ${OUTPUT_DIR}/final.norm.vcf.gz
gemini annotate -f ${OUTPUT_DIR}/final.norm.vcf.gz -a extract -o last -c tumor_aaf -t float -e TUMOR_AAF ${OUTPUT_DIR}/tumor_normal.db
gemini annotate -f ${OUTPUT_DIR}/final.norm.vcf.gz -a extract -o last -c normal_aaf -t float -e NORMAL_AAF ${OUTPUT_DIR}/tumor_normal.db
gemini annotate -f ${OUTPUT_DIR}/final.norm.vcf.gz -a extract -o last -c tumor_lod -t float -e TUMOR_LOD ${OUTPUT_DIR}/tumor_normal.db
gemini annotate -f ${OUTPUT_DIR}/final.norm.vcf.gz -a extract -o last -c normal_lod -t float -e NORMAL_LOD ${OUTPUT_DIR}/tumor_normal.db
gemini annotate -f ${OUTPUT_DIR}/final.norm.vcf.gz -a extract -o last -c tumor_mabq -t float -e TUMOR_MABQ ${OUTPUT_DIR}/tumor_normal.db
