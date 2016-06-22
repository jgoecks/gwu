#
# Run plate7 FFPE analysis.
#

SAMPLE_IDS=( 14 18 47 48 50 53 )

# check command-line arguments
# if [ $# -ne 3 ]; then
#   echo "Usage: `basename $0`  <tumor_dir>  <normal_dir> <output_dir>" 1>&2
#   exit -1
# fi

# Create tumor-normal databases for plasma and for each FFPE sample.
for id in "PD0${SAMPLE_IDS[@]}"
do
    NORMAL=$(ls -d *${id}WB*)

    # Create tumor-normal databases for each FFPE sample.
    for tumor in $(ls -d *${id}* -I ${NORMAL})
    do
        bash -x ./run_paired_somatic.sh ${tumor} ${NORMAL} ${tumor}
    done
done
