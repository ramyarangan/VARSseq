"""
Using the barcode-variant correspondences from gDNA sequencing and the per-barcode per-UMI 
splicing data, collect and write out splicing statistics per variant.
Also prints out stats on the number of variants covered by at least 10 reads and at least 3 barcodes.

Example:
python collect_designed_variant_barcode_splicing.py ../../gdna_analysis/barcode_consensus/S8_designed_var_barcodes.txt barcode_umi_tables_tophat/S8_rd1_barcode_umi_spliced.txt var_barcode_spliced/S8_rd1_var_barcode_spliced.txt
"""
import sys
import numpy as np 

designed_var_barcode_file = sys.argv[1] # Input: Barcode-variant correspondences from gDNA sequencing
barcode_umi_spliced_file = sys.argv[2] # Input: Per-barcode per-UMI splicing statistics from RNA sequencing
per_var_barcode_spliced_file = sys.argv[3] # Output: Per-variant per-barcode splicing statistics 

f = open(designed_var_barcode_file)
designed_var_lines = f.readlines()
f.close()

designed_var_dict = {}
barcode_stats_dict = {}

# Collect all barcodes corresponding to each variant sequence and prepare to
# read in splicing statistics for each one
ii = 0
cur_var = 1
while (ii < len(designed_var_lines)):
	num_barcodes = designed_var_lines[ii].split(' ')[-1]
	num_barcodes = int(num_barcodes.replace('\n', ''))
	barcode_list = []
	for jj in range(num_barcodes):
		ii += 1
		# Also get the coverage from gDNA sequencing - transformation frequency
		barcode_cov = designed_var_lines[ii].replace('\n', '')
		barcode = barcode_cov.split(' ')[0]
		cov = barcode_cov.split(' ')[1]
		barcode_list += [barcode]
		barcode_stats_dict[barcode] = np.array([0, 0, 0, 0, float(cov)])
	designed_var_dict[cur_var] = barcode_list
	ii += 1
	cur_var += 1

# Compute the total spliced, unspliced and other UMI's per barcode
with open(barcode_umi_spliced_file) as f:
	for line in f:
		barcode = line.split(' ')[0]
		if barcode in barcode_stats_dict:
			spliced = int(line.split(' ')[2])
			unspliced = int(line.split(' ')[3])
			other = int(line.split(' ')[4].replace('\n', ''))
			
			total_reads = float(spliced + unspliced + other)
			spliced = spliced/total_reads
			unspliced = unspliced/total_reads
			other = other/total_reads

			barcode_stats = np.array([1, spliced, unspliced, other, 0])
			barcode_stats_dict[barcode] += barcode_stats

# Collect stats on coverage and splicing rates
total_spliced = 0
total_unspliced = 0 
num_covered = 0 # Num variants with at least 10 reads
num_covered_barcodes = 0 # Num variants with at least 3 barcodes

f = open(per_var_barcode_spliced_file, 'w')

# For each designed variant, write out the per-barcode 
# UMI coverage (copies), spliced vs unspliced vs other read counts, and transformation frequency
for designed_var, barcodes in designed_var_dict.items():
	num_barcodes = 0
	num_barcode_copies = 0
	num_spliced = 0
	num_unspliced = 0
	num_other = 0

	f.write("Var: %d\n" % designed_var)
	for barcode in barcodes:
		(cur_copies, cur_spliced, cur_unspliced, cur_other, cov) = barcode_stats_dict[barcode]
		if cur_copies > 0: 
			f.write("%d %d %d %d %d %s\n" % (cur_copies, cur_spliced, cur_unspliced, cur_other, cov, barcode))
		if cur_copies == 0:
			continue
		num_barcodes += 1
		
		num_barcode_copies += cur_copies
		num_spliced += cur_spliced
		num_unspliced += cur_unspliced
		num_other += cur_other

	num_spliced = int(num_spliced)
	num_unspliced = int(num_unspliced)
	num_other = int(num_other)
	if num_barcodes > 2:
		num_covered_barcodes += 1
	if num_barcode_copies > 10:
		num_covered += 1

	total_spliced += num_spliced
	total_unspliced += num_unspliced
	print("Var: %d: %d %d %d %d %d \n" % \
		(designed_var, num_barcodes, num_barcode_copies, num_spliced, num_unspliced, num_other))

f.close()

print("Total spliced: %d" % total_spliced)
print("Total unspliced: %d" % total_unspliced)
print("Num with at least 10 unique reads: %d" % num_covered)
print("Num with at least 3 unique barcodes: %d" % num_covered_barcodes)
