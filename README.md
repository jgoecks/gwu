# Analysis steps for RCC paired tumor (plasma) - normal (buffy coat) amplicon datasets:

1. Download datasets from Galaxy histories on main:
   ```
   $ python ~/gwu/galaxy/download_history_datasets.py
   ```

2. Run analyses for a plate of data. This create many files; we care most about the BAM for each sample.
   ```
   $ make run_dir TARGET_DIR=~/projects/rcc/Plate8
   // Move results to Plate8 directory
   $ make run_dir TARGET_DIR=~/projects/rcc/Plate9
   // Move results to Plate9 directory
   ```

3. Run matched tumor-normal analysis for plasma-buffy coat for each sample
   ```
   // Move Plate8 results back to main directory
   $ make run_paired_somatic

   ```

    NOTES: needed to fix freebayes-parallel script to point to correct location for vcffirstheader and vcfstreamsort in conda. vcfstreamsort is in vcflib directory, but vcffirstheader needs to be added: https://github.com/vcflib/vcflib/blob/8ac9fd517579134ef3b9797714d20c9c99c18ec6/bin/vcffirstheader

4. Do final query
   ```
   $ make query_all
   ```

5. NOTE: for FFPE samples matched to a plasma/buffy coat sample (e.g. plate7), steps 3 and 4 must be replaced with:
   ```
   $ make ffpe_samples
   $ make ffpe_reports
   ```
