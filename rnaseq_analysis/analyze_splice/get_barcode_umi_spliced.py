"""
Example file to generate per-barcode per-UMI stats on spliced, unspliced, and other read counts
This is run on the compute cluster directly after Tophat alignment and generating fasta files
"""
import sys
from Bio import SeqIO

fa_filename = sys.argv[1]
barcode_umi_counts_filename = sys.argv[2]

SPLICED="CTAAAATGTCTTTT"
UNSPLICED="CTAAAATGGTATGT"

barcode_umi_spliced = {}

with open(fa_filename) as f:
	for record in SeqIO.parse(f, "fasta"):
		seq_id = record.id
		# Get barcode and UMI sequence based on read location
		barcode = seq_id.split('_')[-1][0:12]
		umi = seq_id.split('_')[-1][12:].replace('\n', '')
		seq = record.seq

		idx = 2
		if SPLICED in seq: 
			idx = 0
		if UNSPLICED in seq:
			idx = 1

		new_vals = [0, 0, 0]
		new_vals[idx] += 1

		# Accumulate spliced, unspliced, and other counts per-barcode and per-UMI
		if barcode in barcode_umi_spliced.keys():
			if umi in barcode_umi_spliced[barcode].keys():
				barcode_umi_spliced[barcode][umi][idx] += 1
			else:
				barcode_umi_spliced[barcode][umi] = new_vals
		else:
			barcode_umi_spliced[barcode] = {umi: new_vals}

# Write collated spliced, unspliced, other read count statistics along with barcode and UMI sequences
f = open(barcode_umi_counts_filename, 'w')
total_spliced = 0
total_unspliced = 0 
total_other = 0
for barcode in barcode_umi_spliced.keys():
	for umi in barcode_umi_spliced[barcode].keys():
		[spliced, unspliced, other] = barcode_umi_spliced[barcode][umi]
		f.write("%s %s %d %d %d\n" % (barcode, umi, spliced, unspliced, other))
		total_spliced += spliced
		total_unspliced += unspliced
		total_other += other

f.close()

print("Spliced: %d, Unspliced: %d, Other: %d" % \
	(total_spliced, total_unspliced, total_other))
