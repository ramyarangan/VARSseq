# Overview
Variant and compensatory variant sequence design for VARS-seq to generate structure variant libraries. Library sequences will be designed to specifically alter specified stems and junctions in wildtype sequence secondary structures while leaving the remaining sequence in-tact. 

Contents: 
* `generate_variant_mfes.py` and `score_variants.py`: scripts to generate library sequences 
* `intron_info/`: intron annotation files and stem / junction set specifications required to generate variant sequences
* `variants/`: code for final checks for variant sequences along with the final intron variant sequences designed for this paper
* `barcodes/`: barcode sets that were designed to have at least 2 or 3 pairwise edit distance
* `utrs/`: sequence annotations that were used to generate sequences with broader genomic context for library design

# Sample commands

Sample commands for designing variant library sequences: 
* First sample candidate variant sequences and predict their secondary structures: 
`python generate_variant_mfes.py intron_info/rps9b/rps9b.fasta intron_info/rps9b/rps9b_secstruct.txt intron_info/rps9b/rps9b_bpp.csv output_mfes/rps9b_varmfe_ 238 10000 intron_info/rps9b/rps9b_var_list.txt`
* Then score candidate variant sequences:
`python score_variants.py intron_info/rps9b/rps9b.fasta intron_info/rps9b/rps9b_secstruct.txt intron_info/rps9b/rps9b_bpp.csv intron_info/rps9b/rps9b_var_list.txt output_mfes/rps9b_varmfe_ barcodes/barcodes_n8_k3.20220511_104258.csv intron_info/rps9b/rps9b_fullseq.txt variants/rps9b_variants.csv 4`

The file-types and inputs required are: 
* Wildtype sequence interval in a fasta file (in our case, the intron sequence)
* Dot-bracket secondary structure for wildtype sequence - this secondary structure can be obtained through any method: secondary structure prediction, structure probing experiments, or homology to prior structures
* Predicted base-pair probability matrix for the wildtype sequence
* Variant list file: Each line in this file will specify a stem / junction set for which variants will be designed. For each line, after generating sample sequences and scoring them, a small number of variants will be designed meant to vary the specified stem / junction set. On each line of the file, list one nucleotide per stem or junction that will be mutated in variant sequences followed by the fraction of nucleotides in these stems / junctions that should be mutated. For instance, the line `31 41 54 0.5` implies that the stems / junctions including nucleotides 31, 41, and 54 will be mutated by variant sequences, and 50% of the nucleotides in these intervals will be mutated. 
* Barcode file: We will add one barcode from this file to each sequence variant, with a unique barcode per variant. 
* Full sequence file: A full sequence file including extra sequence context will be used to generate the final library sequences. The first line in this file will have three numbers: 1. The position for the barcode to be inserted in this full sequence. 2. The position for the start of the wildtype sequence (in our case the intron sequence). 3. The position for the start of the library in the full sequence. The second line in this file includes the wildtype sequence with extended sequence context. 

For more detailed input formatting, see examples in the `intron_info` folder. 