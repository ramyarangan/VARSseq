# VARS-seq
Designing sequence libraries for VARS-seq and analyzing gDNA and RNA sequencing data from VARS-seq.

# Overview
This repository contains scripts for: 
* Designing structure variant sequences for VARS-seq in `library_design/` 
	* Sampling random variants and rescue sequences for candidate stem disruption and rescue sequences
	* Predicting secondary structures for sampled sequences
	* Scoring secondary structures for desired properties for disrupting and rescuing stem sets
	* Utilities for processing RNA secondary structures
* Analyzing gDNA sequencing data for VARS-seq in `gdna_analysis/`
	* Pairing barcodes with consensus sequences
	* Obtaining transformation frequencies per barcode
	* Obtaining barcodes for designed variant sequences
	* Plotting designed variant sequence coverage and length distributions for consensus sequences
	* Filtering out barcodes with multiple consensus sequences
* Analyzing RNA sequencing data for VARS-seq in `rnaseq_analysis/`
	* Obtaining per-barcode, per-variant, and per-stem set counts for spliced and unspliced reads after deduplicating UMI's
	* Displaying heatmaps showing the significance of differences between wildtype, variant, and rescue sequence sets
	* Assembling stats tables for the number of spliced and unspliced reads per barcode, grouped by variant sequence

# Designing sequences for VARS-seq

Sample commands for designing variant library sequences: 
* First sample candidate variant sequences and predict their secondary structures: 
	* `python library_design/generate_variant_mfes.py library_design/intron_info/rps9b/rps9b.fasta library_design/intron_info/rps9b/rps9b_secstruct.txt library_design/intron_info/rps9b/rps9b_bpp.csv library_design/output_mfes/rps9b_varmfe_ 238 10000 library_design/intron_info/rps9b/rps9b_var_list.txt`
* Then score candidate variant sequences:
	* `python library_design/score_variants.py library_design/intron_info/rps9b/rps9b.fasta library_design/intron_info/rps9b/rps9b_secstruct.txt library_design/intron_info/rps9b/rps9b_bpp.csv library_design/intron_info/rps9b/rps9b_var_list.txt library_design/output_mfes/rps9b_varmfe_ library_design/barcodes/barcodes_n8_k3.20220511_104258.csv library_design/intron_info/rps9b/rps9b_fullseq.txt library_design/variants/rps9b_variants.csv 4`

Detailed descriptions for the input and output variables and data types for these commands are included in the `library_design/` directory's README.

We recommend carrying out sequence design with the help of a high performance cluster. To accelerate sequence design, we carried out secondary structure prediction and scored secondary structures for each stem set in each intron on separate cores in parallel.

# Producing figures 

Each subdirectory includes a figures directory with example analysis figures and a README detailing how to produce these figures.

# Requirements
* scipy, numpy, matplotlib, seaborn, pandas
* arnie, Contrafold
