"""
Make per-variant tables with stats on the spliced, unspliced, other read counts per barcode.

python make_barcode_stats_table.py var_barcode_spliced/S1_rd1_var_barcode_spliced.txt var_group_files/S1_var_groups.txt QCR9 stats_tables/qcr9.csv  
"""
import sys

barcode_stats_file = sys.argv[1] # like var_barcode_spliced/S1_rd1_var_barcode_spliced.txt
var_groups_file = sys.argv[2] # like var_group_files/S1_var_groups.txt
gene_name = sys.argv[3]
stats_table_file = sys.argv[4]

f = open(barcode_stats_file)
barcode_stats_lines = f.readlines()
f.close()

f = open(var_groups_file)
var_groups_lines = f.readlines()
f.close()

# Read in per-barcode spliced, unspliced, and other read counts
var_to_stats = {}
cur_var = None
cur_var_list = []
for line in barcode_stats_lines:
	if line[0] == 'V': 
		if cur_var is not None:
			var_to_stats[cur_var] = cur_var_list
			cur_var_list = []
		cur_var = int(line.split(' ')[1])
		continue
	line_items = line.split(' ')
	spliced = int(line_items[1])
	unspliced = int(line_items[2])
	other = int(line_items[3])
	barcode= line_items[-1].replace('\n', '')
	cur_var_list += [(barcode, spliced, unspliced, other)]

# Write table as CSV
f = open(stats_table_file, 'w')

ii = 0
# Go through stem / junction sets
while ii < len(var_groups_lines):
	var_info = var_groups_lines[ii].split(' ')
	num_vars = int(var_info[2])
	ii += 1
	# Go through each variant in the stem / junction set
	for jj in range(num_vars):
		var_seq_num = var_groups_lines[ii].split(' ')
		ii += 1
		var_num = int(var_seq_num[0])
		if var_num not in var_to_stats:
			continue
		var_seq = var_seq_num[1].replace('\n', '')
		num_barcodes = len(var_to_stats[var_num])

		# Write out the variant sequence and per-barcode stats
		f.write("%s,%s,%s,%s,%d\n" % \
			(gene_name, var_info[0], var_info[1], var_seq, num_barcodes))
		for (barcode, spliced, unspliced, other) in var_to_stats[var_num]:
			f.write("%s,%d,%d,%d\n" % (barcode, spliced, unspliced, other))

f.close()


