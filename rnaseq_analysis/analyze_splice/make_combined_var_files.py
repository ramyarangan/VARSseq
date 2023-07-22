"""
Used to make files in var_group_files/
Groups variants by stem disruption set

python make_combined_var_files.py ../../library_design/variants/final_variants/rps9a_variants.csv ../../gdna_analysis/ref/rps9a_complete.fa var_group_files/rps9a_var_groups.txt AGTACGT
"""
import sys

var_detail_file = sys.argv[1] # The designed library variants with info on different stem sets
var_file = sys.argv[2] # Full library sequences including regions external to intron
var_group_file = sys.argv[3] # Output file for variant sequences grouped by stem / junction disruption set
fivess_seq = sys.argv[4] # 5' splice site sequence to find start of intron 

f = open(var_detail_file)
var_detail_lines = f.readlines()
f.close()

f = open(var_file)
var_file_lines = f.readlines()
f.close()

var_id_seqs = []
for ii in range(int(len(var_file_lines)/2)):
	var_id = var_file_lines[ii * 2].split("_")[1].replace('\n', '')
	seq = var_file_lines[ii * 2 + 1].replace('\n', '')
	var_id_seqs += [(var_id, seq)]

group_id_to_var_list = {}

for var_detail_line in var_detail_lines:
	var = var_detail_line.split(',')[0]
	fivess_idx_var = var.index(fivess_seq)
	compens = var_detail_line.split(',')[2]
	fivess_idx_compens = -1
	if len(compens) > 0:
		fivess_idx_compens = compens.index(fivess_seq)

	group_id = '_'.join(var_detail_line.split(',')[4].split('_')[:-1])

	var_match = None
	for (var_id, seq) in var_id_seqs:
		if var[fivess_idx_var:] in seq:
			var_match = (var_id, seq)
			break

	compens_match = None
	if len(compens) > 0: 
		for (var_id, seq) in var_id_seqs:
			if compens[fivess_idx_compens:] in seq:
				compens_match = (var_id, seq)
				break

	if var_match is None:
		raise RuntimeError("No match for: %s" % var[fivess_idx_var:])
	if compens_match is None and len(compens) > 0:
		raise RuntimeError("No match for: %s" % compens[fivess_idx_compens:])

	var_key = (group_id, "variant")
	if var_key in group_id_to_var_list.keys():
		group_id_to_var_list[var_key] += [var_match]
	else: 
		group_id_to_var_list[var_key] = [var_match]

	compens_key = (group_id, "compens")
	if len(compens) > 0:
		if compens_key in group_id_to_var_list.keys():
			group_id_to_var_list[compens_key] += [compens_match]
		else: 
			group_id_to_var_list[compens_key] = [compens_match]


# Group idx ID, variant vs compens, num variants
# variant ID's sequence
f = open(var_group_file, 'w')
for group_key in group_id_to_var_list.keys():
	var_list = group_id_to_var_list[group_key]
	f.write("%s %s %d\n" % (group_key[0], group_key[1], len(var_list)))
	for var_item in var_list:
		f.write("%s %s\n" % (var_item[0], var_item[1]))
f.close()
