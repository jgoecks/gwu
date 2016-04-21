#
# Remove MUC4_42-124 because they target highly repetitive regions that will be
# excluded from the analysis.
#
# Usage:
#   python remove_bad_regions.py < all_target_regions.bed > final_target_regions.BED

import sys

for line in sys.stdin:
    chrom, start, end, name = line.split("\t")
    try:
        gene, primerno = name.split("_")
        primerno = int(primerno)
        if gene == "MUC4" and primerno >= 42 and primerno <= 124:
            continue
        print line,
    except:
        print line,
