# Overview
Scripts used for analyzing gDNA sequencing data from VARS-seq to identify consensus library sequences for each barcode. These scripts should be used to process consensus sequences downstream of fgBio and bbmerge (as described in GEO accession ID: GSE209857). This directory includes scripts for collecting barcodes that correspond to desired variant sequences and for generating simple plots summarizing their coverage.

# Plotting
* Generating plots for the size distribution of consensus sequences: plot_consensus_seq_sizes.py consensus_size/ihist1.txt QCR9
* Generating barplots for the coverage of designed variants by barcodes: python plot_variant_coverage.py

Example plots are in `figures`.
