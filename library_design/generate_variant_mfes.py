"""
Generate sample intron variant sequences, predict, and write out their secondary structures
"""

import sys
from arnie.bpps import bpps 
from arnie.mfe import mfe
import numpy as np
import random
from secstruct_util import *

intron_seq_fasta = sys.argv[1] # E.g. intron_info/rps9b/rps9b.fasta
intron_dotbracket_file = sys.argv[2] # E.g. intron_info/rps9b/rps9b_secstruct.txt
intron_bpp_matrix = sys.argv[3] # E.g. intron_info/rps9b/rps9b_bpp.csv
mfe_output_pref = sys.argv[4] # Output prefix for secondary structure MFE files
library_end = int(sys.argv[5]) # Library's 3' end in the fasta sequence
num_structs = int(sys.argv[6]) # Number of secondary structures to generate per variant
node_list_file = sys.argv[7] # E.g. intron_info/rps9b/rps9b_var_list.txt

# Randomly generate a new sequence replacing n nucleotides at pos_list positions in seq
def vary_pos_random(pos_list, seq, n):
	matches_curseq = True
	pos_sample = []
	new_nts = []
	while matches_curseq:
		pos_sample = random.sample(pos_list, n)
		pos_sample = np.array(pos_sample)
		cur_nts = seq[pos_sample]
		nt_list = ["A", "C", "G", "U"]
		new_nts = np.array(random.choices(nt_list, k=n))

		if sum(cur_nts != new_nts) > 0:
			matches_curseq = False

	return pos_sample, new_nts

# Randomly scramble n nucleotides in pos_list positions in seq
def vary_pos_scramble(pos_list, seq, n):
	matches_curseq = True
	pos_sample = []
	shuffle_nts = []

	while matches_curseq:
		pos_sample = random.sample(pos_list, n)
		pos_sample = np.array(pos_sample)
		cur_nts = seq[pos_sample]
		cur_nt_list = list(cur_nts)
		shuffle_nts = random.sample(cur_nt_list, len(cur_nt_list))
		shuffle_nts = np.array(shuffle_nts)
		
		if sum(cur_nts != shuffle_nts) > 0:
			matches_curseq = False

	return pos_sample, shuffle_nts

# Get base-pair partners for the variant positions and compensatory mutations
def get_compens_muts(stem_list, var_pos, var_nts):
	comp_nts = np.copy(var_nts)
	comp_nts[var_nts == 'A'] = ['U']
	comp_nts[var_nts == 'U'] = ['A']
	comp_nts[var_nts == 'C'] = ['G']
	comp_nts[var_nts == 'G'] = ['C']

	comp_pos = []
	for pos in var_pos:
		for stem in stem_list:
			seq_pos_dict = stem.get_seq_pos_to_strand1_pos()
			if pos in stem.strand1_nts:
				comp_pos += [stem.strand2_nts[seq_pos_dict[pos]]]
	comp_pos = np.array(comp_pos)

	return comp_pos, comp_nts

# given a stem list, generate randomized variants that: 
# destabilize all the stems
# have compensatory mutants if available
# generate random variants through scrambles, or by mutating a max of N residues per stem
# predict and save MFE structures
def write_stem_var_mfes(stem_list, seq, library_end, out_file, \
	perc_mutate=1, max_iter=1000, do_scramble=False):
	f = open(out_file, 'a')

	seq = seq.replace('T', 'U')
	seq = np.array(list(seq))

	all_variants = []

	ii = 0
		
	# Build variant and compensatory variant over all stems that will be modified
	variant = np.copy(seq)
	compens_variant = np.copy(seq)
	do_compens = True

	# To do compensatory mutations, library must include at least some 
	# nucleotides in all stems in this set
	for stem in stem_list:
		if library_end < min(stem.strand2_nts):
			do_compens = False

	stem_res = []
	for stem in stem_list:
		# Get range of values where mutations should be attempted
		max_3p = max(stem.strand1_nts) + 1
		mut_range = max(stem.strand2_nts) - min(stem.strand2_nts) + 1
		if do_compens:
			mut_range = min(library_end - min(stem.strand2_nts), mut_range)
		if mut_range < 0: 
			continue
		min_3p = max_3p - mut_range
		if library_end > max(stem.strand1_nts):
			stem_res += list(range(min_3p, max_3p))
		else:
			stem_res += list((range(min_3p, library_end)))

	print(stem_res)

	num_var_res = int(perc_mutate * len(stem_res))

	for jj in range(max_iter):
		if jj % 20 == 0:
			print("Doing structures %d of %d \n" % (jj, max_iter))

		# Get sample variant and compensatory variant positions / nts
		if do_scramble:
			(var_pos, var_vals) = vary_pos_scramble(stem_res, seq, num_var_res)
		else:
			(var_pos, var_vals) = vary_pos_random(stem_res, seq, num_var_res)

		(comp_pos, comp_vals) = get_compens_muts(stem_list, var_pos, var_vals)

		test_variant = np.copy(variant)
		test_variant[var_pos] = var_vals
		test_compens_variant = np.copy(compens_variant)
		
		if do_compens:
			test_compens_variant[var_pos] = var_vals
			test_compens_variant[comp_pos] = comp_vals

		# Write out variant candidate sequences and MFEs
		mfe_variant = mfe(test_variant, package='contrafold')
		f.write("Variant: %s\n" % ''.join(test_variant))
		f.write("%s\n" % mfe_variant)

		# Write out compensatory variant candidate sequences and MFEs 
		if do_compens:
			mfe_compens_variant = mfe(test_compens_variant, package='contrafold')
			f.write("Compensatory: %s\n" % ''.join(test_compens_variant))
			f.write("%s\n" % mfe_compens_variant)

	f.close()


# given a loop, generate random scrambles
def write_junction_var_mfes(junctions, seq, library_end, out_file, \
		perc_mutate=1, max_iter=1000, do_scramble=False):
	f = open(out_file, 'a')

	seq = seq.replace('T', 'U')
	seq = np.array(list(seq))

	all_variants = []

	ii = 0
		
	# Build variant and compensatory variant over all stems that will be modified
	variant = np.copy(seq)

	junc_res = []
	for junction in junctions:
		# Get range of values where mutations should be attempted
		for nt in junction.nts:
			if nt < library_end:
				junc_res += [nt]	
		
	num_var_res = int(perc_mutate * len(junc_res))

	for jj in range(max_iter):
		if jj % 20 == 0:
			print("Doing structures %d of %d \n" % (jj, max_iter))

		# Get sample variant and compensatory variant positions / nts
		if do_scramble:
			(var_pos, var_vals) = vary_pos_scramble(junc_res, seq, num_var_res)
		else:
			(var_pos, var_vals) = vary_pos_random(junc_res, seq, num_var_res)

		test_variant = np.copy(variant)
		test_variant[var_pos] = var_vals

		mfe_variant = mfe(test_variant, package='contrafold')
		f.write("Variant: %s\n" % ''.join(test_variant))
		f.write("%s\n" % mfe_variant)

	f.close()

# Iterate through all node lists and generate / write variant sequences and MFEs
def write_vars_from_lists(node_lists, nt_to_node_dict, bpp_matrix, intron_seq):
	for (node_list, perc_change) in node_lists:
		stems, junctions = get_stems_junctions_from_nodelist(node_list, nt_to_node_dict, bpp_matrix)

		out_file = get_nodelist_filename(node_list, perc_change, mfe_output_pref)

		if len(stems) > 0:
			write_stem_var_mfes(stems, intron_seq, library_end, out_file, \
				perc_mutate=perc_change, max_iter=num_structs)
			write_stem_var_mfes(stems, intron_seq, library_end, out_file, \
				perc_mutate=perc_change, max_iter=num_structs, do_scramble=True)

		if len(junctions) > 0:
			write_junction_var_mfes(junctions, intron_seq, library_end, out_file, \
				perc_mutate=perc_change, max_iter=num_structs)
			write_junction_var_mfes(junctions, intron_seq, library_end, out_file, \
				perc_mutate=perc_change, max_iter=num_structs, do_scramble=True)

# Get basic intron info
intron_seq = get_seq(intron_seq_fasta)

# Get stem and junction nodes
dotbracket = get_dotbracket(intron_dotbracket_file)
bpp_matrix = get_bpp_matrix(intron_bpp_matrix)
G, nt_to_node_dict, all_stems = get_stem_graph(dotbracket, bpp_matrix)

for stem in all_stems:
	print(stem)

# generate variants for specified sets of stems and loops
node_lists = get_node_lists(node_list_file)

write_vars_from_lists(node_lists, nt_to_node_dict, bpp_matrix, intron_seq)

