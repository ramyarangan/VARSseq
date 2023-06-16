import sys
from arnie.bpps import bpps 
from arnie.mfe import mfe
from secstruct_util import *
import random 

intron_seq_fasta = sys.argv[1]
intron_dotbracket_file = sys.argv[2]
intron_bpp_matrix = sys.argv[3]
node_lists_file = sys.argv[4]
mfe_dir = sys.argv[5]
barcode_file = sys.argv[6]
full_seq_file = sys.argv[7] # Should include 5'UTR and 50nt past intron
variant_outfile = sys.argv[8]
num_top_variants = int(sys.argv[9])

LIB_SIZE = 300

def run_mfe_checks(mfe_variant, mfe_wt, bp_matrix_wt, unpaired=[]):
	bp_matrix_variant = get_bp_matrix(mfe_variant)

	mfe_wt = np.array(list(mfe_wt))
	mfe_variant = np.array(list(mfe_variant))
	unpaired = np.array(unpaired)

	# How many positions in unpaired set are paired to something?
	num_fp_paired_var_region = 0	
	if len(unpaired) > 0:
		num_fp_paired_var_region = sum(mfe_variant[unpaired] != '.')

	# How many base-pairs are present in the unpaired set?
	num_fp_bp_var_region = 0
	if len(unpaired) > 0: 
		num_fp_bp_var_region = np.sum(bp_matrix_variant[unpaired,:][:,unpaired])/2

	var_region_stats = (num_fp_paired_var_region, num_fp_bp_var_region)

	# Assemble base-pairs outside the unpaired set
	full_pos_set = set(range(len(mfe_variant)))
	unpaired_set = set(unpaired)
	remaining_set = full_pos_set.difference(unpaired_set)
	remaining_pos = np.array(list(remaining_set))

	bp_subset_variant = bp_matrix_variant[remaining_pos,:][:,remaining_pos]
	bp_subset_variant = bp_subset_variant.astype(int)
	bp_subset_wt = bp_matrix_wt[remaining_pos,:][:,remaining_pos]
	bp_subset_wt = bp_subset_wt.astype(int)
	total_bp = np.sum(bp_subset_wt)/2

	# How many base-pairs outside the unpaired set are maintained between WT vs variant?
	num_tp_bp = np.sum(bp_subset_variant & bp_subset_wt)/2

	# How many base-pairs outside the unpaired set are lost from WT to variant?
	num_fn_bp = np.sum(~bp_subset_variant & bp_subset_wt)/2

	# How many base-pairs outside the unpaired set are gained from WT to variant?
	num_fp_bp = np.sum(bp_subset_variant & ~bp_subset_wt)/2

	# print(num_fp_bp)
	# print(num_fn_bp)

	cons_region_stats = (num_tp_bp, num_fn_bp, num_fp_bp, total_bp)

	return var_region_stats, cons_region_stats

def get_score(var_region_stats, cons_region_stats):
	var_region_errors = var_region_stats[0] + var_region_stats[1]
	(_, num_fn_bp, num_fp_bp, _) = cons_region_stats
	cons_region_errors = num_fn_bp + 0.5 * num_fp_bp
	return var_region_errors + cons_region_errors

def get_score_var_compens(score_variant, score_compens):
	return score_variant + score_compens

def get_variant_scores(all_variants, all_compens_variants, mfe_wt, \
	bp_matrix_wt, unpaired):
	has_compens = len(all_compens_variants) > 0

	mfe_check_scores = []
	total_scores = []
	for ii, variant in enumerate(all_variants):
		if ii %  100 == 0:
			print("Scoring variant %d of %d" % (ii, len(all_variants)))
		mfe_variant = all_variants[ii]
		stats = run_mfe_checks(mfe_variant, mfe_wt, bp_matrix_wt, unpaired=unpaired)
		var_score = get_score(stats[0], stats[1])
		cur_stats = list(stats[0]) + list(stats[1])
		total_score = var_score

		if has_compens:
			mfe_compens = all_compens_variants[ii]
			stats = run_mfe_checks(mfe_compens, mfe_wt, bp_matrix_wt)
			compens_score = get_score(stats[0], stats[1])
			cur_stats += list(stats[0]) + list(stats[1])
			total_score = get_score_var_compens(var_score, compens_score)
		else:
			cur_stats += [-1] * len(cur_stats)

		# print(total_score)
		mfe_check_scores += [cur_stats]
		total_scores += [total_score]

	return mfe_check_scores, total_scores

def check_bpp(bpp_matrix_wt, variant, bpp_perc_change, unpaired=[]):
	return True

def has_restriction_sites(seq):
	# BglII: AGATCT
	for ii in range(len(seq) - 6):
		if seq[ii:(ii + 6)] == 'AGATCT':
			print(seq[ii:(ii+ 6)])
			print(seq)
			return True

	# XCMI: CCANNNNNNNNNTGG
	for ii in range(len(seq) - 15):
		if seq[ii:(ii + 3)] == 'CCA':
			if seq[(ii+12):(ii+15)] == 'TGG':
				print(seq[ii:(ii+ 15)])
				print(seq)
				return True

	return False

def get_top_variants(node_list, perc_change, stems, junctions, 
	mfe_wt, bp_matrix_wt, bpp_matrix_wt, bpp_perc_change=0.2):
	print("Getting top variants for:")
	print(node_list)
	nodelist_filename = get_nodelist_filename(node_list, perc_change, mfe_dir)

	f = open(nodelist_filename)
	nodelist_lines = f.readlines()
	f.close()

	all_variants = []
	all_secstructs = []
	all_compens_variants = []
	all_compens_secstructs = []
	
	for ii in range(int(len(nodelist_lines)/2)):
		seq_line = nodelist_lines[ii * 2]
		struct_line = nodelist_lines[ii * 2 + 1]

		seq = seq_line.split(" ")[1].replace('\n', '')
		struct = struct_line.replace('\n', '')
		
		if seq_line[0] == "V":
			all_variants += [seq]
			all_secstructs += [struct]
		else:
			all_compens_variants += [seq]
			all_compens_secstructs += [struct]

	has_compens = len(all_compens_variants) > 0
	if has_compens and len(all_compens_variants) != len(all_variants):
		print("Warning: Not all variants had compensatory variant")

	unpaired = []
	for stem in stems:
		unpaired += stem.strand1_nts + stem.strand2_nts
	for junction in junctions:
		unpaired += junction.nts

	mfe_check_scores, total_scores = get_variant_scores(all_secstructs, \
		all_compens_secstructs, mfe_wt, bp_matrix_wt, unpaired)


	# Sort variants by score
	total_scores = np.array(total_scores)
	sort_idxs = np.argsort(total_scores)
	total_scores = total_scores[sort_idxs]

	mfe_check_scores = np.array(mfe_check_scores)[sort_idxs]
	all_variants = np.array(all_variants)[sort_idxs]
	all_secstructs = np.array(all_secstructs)[sort_idxs]
	if has_compens:
		all_compens_variants = np.array(all_compens_variants)[sort_idxs]
		all_compens_secstructs = np.array(all_compens_secstructs)[sort_idxs]

	# Do final checks until enough variants have been collected
	final_variants = []
	final_structs = []
	final_variant_set = set()
	final_compens_variants = []
	final_compens_structs = []
	final_compens_variant_set = set()
	final_scores = []
	# Get 2N variants if compensatory variants are not present, otherwise
	# get N of each.
	num_variants = num_top_variants
	if not has_compens:
		num_variants = 2 * num_top_variants

	for ii, variant in enumerate(all_variants):
		if len(final_variants) == num_variants:
			break

		# Don't allow for repeat variants that were randomly generated
		if variant in final_variant_set:
			continue
		final_variant_set.add(variant)
		if has_compens:
			if all_compens_variants[ii] in final_compens_variant_set:
				continue
			final_compens_variant_set.add(all_compens_variants[ii])

		# Check base-pair probability matrix to minimize changes (try keeping <10%)
		# For now skipping this
		if not check_bpp(bpp_matrix_wt, variant, bpp_perc_change, unpaired=unpaired):
			print("Variant failed BPP check 1")
			continue
		if has_compens:
			if not check_bpp(bpp_matrix_wt, all_compens_variants[ii], bpp_perc_change):
				print("Compensatory variant failed BPP check 2")
				continue

		# Check for any new restriction sites
		if has_restriction_sites(variant):
			print("Variant failed restriction site check 1")
			continue
		if has_compens:
			if has_restriction_sites(all_compens_variants[ii]):
				print("Compensatory variant failed restriction site check 1")
				continue
		
		# Collect final passing variants and compensatory variants
		final_variants += [variant]
		final_structs += [all_secstructs[ii]]
		if has_compens:
			final_compens_variants += [all_compens_variants[ii]]
			final_compens_structs += [all_compens_secstructs[ii]]
		final_scores += [tuple(list(mfe_check_scores[ii]) + [total_scores[ii]])]
		print(mfe_wt)
		print(all_secstructs[ii])
		if has_compens:
			print(all_compens_secstructs[ii])
		print(mfe_check_scores[ii])
		print(total_scores[ii])

	return final_variants, final_structs, final_compens_variants, final_compens_structs, final_scores


def score_mfe_barcode(mfe_variant, var_intron_pos, mfe_barcode, bc_intron_pos, intron_len):
	bp_matrix_variant = get_bp_matrix(mfe_variant)
	bp_matrix_barcode = get_bp_matrix(mfe_barcode)

	score = 0

	# Number of base-pairs within intron changed
	var_intron_end = var_intron_pos+intron_len
	var_in_intron = bp_matrix_variant[var_intron_pos:var_intron_end,:][:,var_intron_pos:var_intron_end]
	bc_intron_end = bc_intron_pos + intron_len
	bc_in_intron = bp_matrix_barcode[bc_intron_pos:bc_intron_end,:][:,bc_intron_pos:bc_intron_end]

	score += np.sum(np.not_equal(var_in_intron, bc_in_intron)/2)

	# Number of base-pairs external to intron changed 
	# (bp partner identity outside intron doesnt matter)
	var_out_intron_idxs = list(np.arange(var_intron_pos)) + list(np.arange(var_intron_end, len(mfe_variant)))
	var_out_intron = bp_matrix_variant[var_intron_pos:var_intron_end,:][:,var_out_intron_idxs]
	var_external_bps = np.sum(var_out_intron, axis=1)
	bc_out_intron_idxs = list(np.arange(bc_intron_pos)) + list(np.arange(bc_intron_end, len(mfe_barcode)))
	bc_out_intron = bp_matrix_barcode[bc_intron_pos:bc_intron_end,:][:,bc_out_intron_idxs]
	bc_external_bps = np.sum(bc_out_intron, axis=1)

	score += np.sum(np.not_equal(var_external_bps, bc_external_bps))

	return score


def add_unique_barcode(variant, barcodes, full_seq, intron_pos, barcode_pos, variant_pos):
	variant = variant[variant_pos:]
	full_seq[intron_pos:(intron_pos + len(variant))] = list(variant)
	full_mfe_variant = mfe(''.join(full_seq[:barcode_pos] + full_seq[(barcode_pos+8):]), package='contrafold')

	not_found = True
	ii = 0
	idx_order = list(range(len(barcodes)))
	random.shuffle(idx_order)

	best_score = 1000
	new_barcodes = barcodes
	new_full_seq = full_seq

	while ii < 25 and not_found:
		idx = idx_order[ii]
		full_seq[barcode_pos:(barcode_pos + 8)] = list(barcodes[idx].replace('T', 'U'))
		full_mfe_barcode = mfe(''.join(full_seq), package='contrafold')

		score = score_mfe_barcode(full_mfe_variant, intron_pos-8, \
			full_mfe_barcode, intron_pos, len(variant))
		print("Score: %d" % score)

		if score < best_score:
			new_barcodes = barcodes[:idx] + barcodes[(idx+1):]
			new_full_seq = full_seq
			best_score = score

		# Verify that structure with new barcode is relatively unchanged
		if best_score < 5:
			print("Found unique barcode iteration: %d" % ii)
			not_found = False
			print(barcodes[idx])

		ii += 1
	
	return ''.join(new_full_seq), new_barcodes

# Add randomized barcode, verify that secondary structure is constant
def get_full_variants(all_variants):
	f = open(barcode_file)
	barcode_lines = f.readlines()[1:]
	f.close()

	barcodes = []
	for barcode_line in barcode_lines:
		barcodes += [barcode_line.split(',')[0]]

	f = open(full_seq_file)
	full_seq_lines = f.readlines()
	f.close()
	full_seq = list(full_seq_lines[1].replace('\n', '').replace('T', 'U'))
	barcode_pos = int(full_seq_lines[0].split(" ")[0])
	intron_pos = int(full_seq_lines[0].split(" ")[1])
	library_pos = int(full_seq_lines[0].split(" ")[2])
	variant_pos = 0
	if len(full_seq_lines[0].split(" ")) > 3:
		variant_pos = int(full_seq_lines[0].split(" ")[3])
	all_full_variants = []

	print("Adding barcodes")
	for ii, (variants, structs, compens_variants, compens_structs, scores, \
		node_list, perc_change) in enumerate(all_variants):
		print("Adding barcodes for: %s" % str(node_list))

		full_variants = []
		for variant in variants: 
			full_variant, barcodes = add_unique_barcode(variant, barcodes, full_seq, \
				intron_pos, barcode_pos, variant_pos)
			full_variants += [full_variant[library_pos:(library_pos+LIB_SIZE)]]

		full_compens_variants = []
		for variant in compens_variants: 
			full_compens_variant, barcodes = add_unique_barcode(variant, barcodes, full_seq, \
				intron_pos, barcode_pos, variant_pos)
			full_compens_variants += [full_compens_variant[library_pos:(library_pos+LIB_SIZE)]]

		all_full_variants += [(full_variants, structs, full_compens_variants, compens_structs, scores, \
			node_list, perc_change)]

	return all_full_variants

def write_variants(all_variants):
	f = open(variant_outfile, 'w')

	for (full_variants, structs, full_compens_variants, compens_structs, scores, \
		node_list, perc_change) in all_variants:
		tag = get_nodelist_tag(node_list, perc_change)

		for ii, variant in enumerate(full_variants):
			variant = variant.replace('U', 'T')
			struct = structs[ii]
			(num_fp_paired_var_region, num_fp_bp_var_region, num_tp_bp, num_fn_bp, \
				num_fp_bp, _, _, _, num_tp_bp_compens, num_fn_bp_compens, \
				num_fp_bp_compens, _, total_score) = scores[ii]
			compens_variant = ""
			compens_struct = ""
			if len(full_compens_variants) > 0:
				compens_variant = full_compens_variants[ii]
				compens_variant = compens_variant.replace('U', 'T')
				compens_struct = compens_structs[ii]
			f.write("%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n" % (variant, struct, \
				compens_variant, compens_struct, tag, num_fp_paired_var_region, \
				num_fp_bp_var_region, num_tp_bp, num_fn_bp, num_fp_bp, num_tp_bp_compens, \
				num_fn_bp_compens, num_fp_bp_compens, total_score))

	f.close()

# Get basic intron info
intron_seq = get_seq(intron_seq_fasta)
mfe_wt = mfe(intron_seq, package='contrafold')
bp_matrix_wt = get_bp_matrix(mfe_wt)
bpp_matrix_wt = bpps(intron_seq, package='contrafold')

dotbracket = get_dotbracket(intron_dotbracket_file)
bpp_matrix = get_bpp_matrix(intron_bpp_matrix)
G, nt_to_node_dict, all_stems = get_stem_graph(dotbracket, bpp_matrix)

# Get region list for evaluating variants in
node_lists = get_node_lists(node_lists_file)

# Should be a list of ((variant, compens_variant), (node_list, perc_change))
all_variants = []
for (node_list, perc_change) in node_lists:
	stems, junctions = get_stems_junctions_from_nodelist(node_list, nt_to_node_dict, bpp_matrix)
	(variants, structs, compens_variants, compens_structs, scores) = get_top_variants(node_list, perc_change, \
		stems, junctions, mfe_wt, bp_matrix_wt, bpp_matrix_wt)
	all_variants += [(variants, structs, compens_variants, compens_structs, scores, node_list, perc_change)]

# Add 8 copies of wildtype
all_variants += [([intron_seq] * 8, [mfe_wt] * 8, [], [], [[-1] * 13] * 8, [0], 0)]

full_variants = get_full_variants(all_variants)
write_variants(full_variants)
