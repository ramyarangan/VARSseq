# Overview
Scripts used for analyzing gDNA sequencing data from VARS-seq to identify consensus library sequences for each barcode. These scripts should be used to process consensus sequences downstream of fgBio and bbmerge (as described in GEO accession ID: GSE209857). This directory includes scripts for collecting barcodes that correspond to desired variant sequences and for generating simple plots summarizing their coverage.
* `barcode_consensus/*py` Scripts to process data from fgBio and bbmerge to obtain barcode-consensus sequence correspondences.
* `barcode_consensus/barcode_sets/`: Final barcodes for each designed variant sequence
* `consensus_fa/`: Consensus sequences for each barcode, with and without read coverage values
* `consensus_size/`: Histograms for the length of consensus sequences across all barcodes for each sample
* `consensus_tags/`: `S*_consensus_tags.txt` is assembled from fgBio's CallMolecularConsensusReads. `S*_tag_blacklist.txt` is a list of barcodes that have multiple high-coverage consensus sequences which we will ignore. `S*_tags_collated.txt` is a list of all barcode sequences with information on the consensus read IDs and coverage per barcode.
* `ref/`: Reference files for each sequence, including the intron sequence (e.g. `qcr9.fa`), list of full variant sequences excluding barcodes (e.g. `qcr9_complete.fa`), and variant sequences in the intron region (e.g. `qcr9_vars.fa`. 

# Plotting
* Generating plots for the size distribution of consensus sequences (Fig. S17A): `plot_consensus_seq_sizes.py consensus_size/ihist1.txt QCR9`
* Generating barplots for the coverage of designed variants by barcodes (Fig. S17B): `python plot_variant_coverage.py`.

Example plots are in `figures`.
