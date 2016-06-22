#
# Run plate7 reports.
#

SAMPLE_IDS=( PD014 PD018 PD047 PD048 PD050 PD053 )

# Create tumor-normal databases for plasma and for each FFPE sample.
for id in "${SAMPLE_IDS[@]}"
do
    rm -f *${id}*results.txt ${id}_final.txt
    ls -d *${id}* | parallel -j 1 \
    "gemini query --header -q 'select chrom, start, ref, alt, gene, transcript, \
                               impact, tumor_lod, normal_lod, tumor_aaf, normal_aaf, \
                               gt_alt_depths, gt_ref_depths, gt_depths from variants \
                               where normal_lod > 3.5 and tumor_aaf >= 0.01 and \
                               (normal_aaf <= 0.001 OR tumor_aaf/normal_aaf > 2.5)' \
                                {}/tumor_normal.db > {}_results.txt"
    python ffpe_report.py ${id} > ${id}_final.txt
done
