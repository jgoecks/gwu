#
# Process a directory of sequence data.
#
# set home directory
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PRIMERS=$1
FASTA=$2
BOWTIE2_INDICES=$3
TARGET_DIR=$4

parallel -j 1 --xapply "bash -x ${HOME_DIR}/run.sh {1} {2} ${PRIMERS} ${FASTA} ${BOWTIE2_INDICES} test" ::: ${TARGET_DIR}/*R1*.fastq ::: ${TARGET_DIR}/*R2*.fastq
