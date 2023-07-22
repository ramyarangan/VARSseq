# Overview
Analysis of RNA-seq data from VARS-seq to compare retained intron fractions and normalized mRNA levels between variant, rescue, and wildtype sequences for stem and junction sets of interest.

* `analyze_splice/*py`: Scripts to process and plot RNA-seq data
* `barcode_umi_tables_tophat/`: Per-barcode per-UMI read counts for spliced, unspliced, and other reads
* `var_barcode_spliced/`: Spliced and unspliced read counts organized by barcode and variant sequence
* `var_group_files/`: Variant sequences and ID's grouped by stem / junction disruption sets.
* `figures/`: Sample violin plots and heatmaps
* `stats_tables/`: Tables of per-barcode spliced, unspliced and other read count statistics grouped by variant sequence and construct. 
* `ref/`: Reference sequences used for tophat alignment

# Plotting
* Generate violin plots for stem / junction sets of interest, either comparing the retained intron fraction (RI) or normalized mRNA levels (mRNA). These violin plots are used in Fig. 6 and Fig. S18
	* `python plot_signif_variants.py var_barcode_spliced/S5_rd1_var_barcode_spliced.txt var_group_files/S5_var_groups.txt RPS9A RI`
	* `python plot_signif_variants.py var_barcode_spliced/S5_rd1_var_barcode_spliced.txt var_group_files/S5_var_groups.txt RPS9A mRNA`
* Generate heatmaps comparing variant, rescue, and wildtype sequences based on retained intron fractions or normalized mRNA levels for all constructs: `python plot_heatmap_signif_variants.py`. These heatmaps are used in Fig. S19.

