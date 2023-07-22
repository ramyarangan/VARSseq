"""
Make heatmaps showing p-values for all comparisons between wildtype, variant, and rescue groups 
Can compare using RI fractions or normalized mRNA spliced levels

python plot_heatmap_signif_variants.py
"""

import numpy as np 
import seaborn as sns
from matplotlib import pyplot as plt 
import scipy.stats as stats
import pandas as pd

MIN_UMI_COV = 10


# Get the per variant and wildtype lists of per-barcode RI or mRNA stats
def get_variant_stats(var_barcode_spliced_file, do_RI=True):

	f = open(var_barcode_spliced_file)
	var_barcode_spliced_lines = f.readlines()
	f.close()

	var_stats = {}

	cur_barcode_stats = []
	cur_id = -1
	for var_line in var_barcode_spliced_lines:
		if var_line[0] == "V":
			if cur_id > -1:
				var_stats[cur_id] = np.array(cur_barcode_stats)
			cur_id = int(var_line.split(" ")[-1].replace('\n', ''))
			cur_barcode_stats = []
			continue
		total_umis = int(var_line.split(' ')[0])
		if total_umis < MIN_UMI_COV:
			continue
		spliced = int(var_line.split(' ')[1])
		unspliced = int(var_line.split(' ')[2])
		if spliced + unspliced < MIN_UMI_COV:
			continue
		barcode_cov = int(var_line.split(' ')[4])
		new_stats = spliced/barcode_cov
		if do_RI:
			new_stats = unspliced/(spliced + unspliced)
		cur_barcode_stats += [new_stats]

	wt_stats = np.array(cur_barcode_stats)

	var_stats["wt"] = wt_stats

	return var_stats


# Get lists of per-barcode stats by stem set / junction set
def get_group_stats(var_stats, var_group_file):
	var_groups = {}
	group_stats = {}
	group_labels = {}

	f = open(var_group_file)
	var_group_lines = f.readlines()
	f.close()

	ii = 0
	while ii < len(var_group_lines):
		# Collect all the variant sequence ID's for a given stem / junction set
		group_id_line = var_group_lines[ii]
		group_id = (group_id_line.split(' ')[0], group_id_line.split(' ')[1])
		num_vars = int(group_id_line.split(' ')[2].replace('\n', ''))
		label = group_id_line.split(' ')[3]
		if 'comb' not in label:
			label += '_' + group_id_line.split(' ')[4]
		var_ids = []
		for jj in range(num_vars):
			ii += 1
			var_id = int(var_group_lines[ii].split(' ')[0])
			var_ids += [var_id]
		unique_var_ids = list(set(var_ids))

		# Get all the RI or mRNA stats per stem / junction set
		cur_group_stats = []
		for var_id in unique_var_ids:
			if var_id in var_stats.keys():
				cur_group_stats += list(var_stats[var_id])
		var_groups[group_id] = unique_var_ids
		group_stats[group_id] = np.array(cur_group_stats) 
		group_labels[group_id] = label
		var_ids = []
		ii += 1

	return group_stats, group_labels


def mean_statistic(x, y, axis):
	return np.mean(x, axis=axis) - np.mean(y, axis=axis)

# Compare two lists using the permutation test with the mean statistic
def compare_stats_signif(stats1, stats2):
	signif_test = stats.permutation_test((stats1, stats2), mean_statistic, n_resamples=100000, \
			vectorized=True)
	pval = signif_test.pvalue

	is_greater = 1
	if np.mean(stats1) > np.mean(stats2): 
		is_greater = -1

	min_size = min(len(stats1), len(stats2))

	return is_greater, pval, min_size

# Compare stats for variant vs wildtype, variant vs rescue, and rescue vs wildtype
def compare_groups(group_stats, wt_stats, group_labels):
	all_group_comps = []
	all_group_labels = []

	num_groups = 0
	for group_id in group_stats.keys():
		if group_id[1] != 'variant':
			continue

		# For each variant, get stats for the rescue variant if available
		var_stats = group_stats[group_id]
		compens_id = (group_id[0], 'compens')
		compens_stats = None
		if compens_id in group_stats.keys():
			compens_stats = group_stats[compens_id]

		# Compare variant vs wildtype
		wt_var_compare = (np.nan, np.nan, 0)
		if len(wt_stats) > 1 and len(var_stats) > 1:
			is_greater, p_val, min_size = compare_stats_signif(wt_stats, var_stats)
			wt_var_compare = (is_greater, p_val, min_size)
			if np.isnan(wt_var_compare[1]):
				wt_var_compare = (is_greater, 1, min_size)

		wt_compens_compare = (np.nan, np.nan, 0)
		var_compens_compare = (np.nan, np.nan, 0)
		if len(wt_stats) > 1 and (compens_stats is not None and len(compens_stats) > 1):
			# Compare wildtype vs rescue
			is_greater, p_val, min_size = compare_stats_signif(wt_stats, compens_stats)
			wt_compens_compare = (is_greater, p_val, min_size)
			if np.isnan(wt_compens_compare[1]):
				wt_compens_compare = (is_greater, 1, min_size)

			# Compare variant vs rescue
			is_greater, p_val, min_size = compare_stats_signif(var_stats, compens_stats)
			var_compens_compare = (is_greater, p_val, min_size)
			if np.isnan(var_compens_compare[1]):
				var_compens_compare = (is_greater, 1, min_size)

		# Collect all comparisons' p-values
		all_group_comps += [(wt_var_compare, wt_compens_compare, var_compens_compare)]
		all_group_labels += [group_labels[group_id]]

	return all_group_comps, all_group_labels


# Get a count of all significant comparisons (comparing to pval_cutoff) 
# between variant and wildtype sequences
def get_signif_stats(var_barcode_spliced_file, var_group_file, do_RI=True, do_both=True, pval_cutoff = 0.05):
	# Get all stats grouped by stem / junction set 
	var_stats = get_variant_stats(var_barcode_spliced_file, do_RI=do_RI)
	group_stats, group_labels = get_group_stats(var_stats, var_group_file)

	# Get p-values for comparisons between variant, rescue, and wildtype sequences
	all_group_comps, all_group_labels = \
		compare_groups(group_stats, var_stats['wt'], group_labels)

	# Compare sequences using the other metric 
	# (count a stem / junction set as significant if there is a signif change in either RI or mRNA)
	all_group_comps_2 = None 
	all_group_labels_2 = None
	if do_both:
		var_stats = get_variant_stats(var_barcode_spliced_file, do_RI=(not do_RI))
		group_stats, group_labels = get_group_stats(var_stats, var_group_file)

		all_group_comps_2, all_group_labels_2 = \
			compare_groups(group_stats, var_stats['wt'], group_labels)

	# Count the number of significant shifts along with the totals for stems and junctions
	stem_signif = 0
	loop_signif = 0
	stem_cnt = 0
	loop_cnt = 0
	for ii, group_comp in enumerate(all_group_comps):
		label = all_group_labels[ii].split('_')[-1].replace('\n', '')
		group_comp_2 = None
		if all_group_labels_2 is not None:
			for jj, group_label in enumerate(all_group_labels_2):
				if group_label == all_group_labels[ii]:
					group_comp_2 = all_group_comps_2[jj]

		node_is_signif = False
		if group_comp[0][1] < pval_cutoff:
			node_is_signif = True 
		if group_comp_2 is not None and group_comp_2[0][1] < pval_cutoff:
			node_is_signif = True

		if label == 'stem':
			stem_cnt += 1
			if node_is_signif:
				stem_signif += 1

		if label == 'loop':
			loop_cnt += 1
			if node_is_signif:
				loop_signif += 1

	return stem_signif, loop_signif, stem_cnt, loop_cnt


# Create a heatmap showing the comparisons between variant, rescue, and wildtype sequences
def make_heatmap(var_barcode_spliced_file, var_group_file, plt_title, do_RI=True):
	var_stats = get_variant_stats(var_barcode_spliced_file, do_RI=do_RI)
	group_stats, group_labels = get_group_stats(var_stats, var_group_file)

	all_group_comps, all_group_labels = \
		compare_groups(group_stats, var_stats['wt'], group_labels)

	# Keep track of direction of comparison, significance
	heatmap_arr = []
	# Keep track of minimum number of barcodes included in the comparison
	size_arr = []
	# Keep track of stem / junction set label
	labels = []

	for ii, group_comp in enumerate(all_group_comps):
		if "wt" in all_group_labels[ii]:
			continue
		labels += [all_group_labels[ii]]
		sign = 1
		if do_RI:
			sign = -1
		# Color heatmap by log p-value, with sign indicating direction of effect
		signed_log_p = [sign * x[0] * np.log10(x[1]) for x in group_comp] 

		heatmap_arr += [signed_log_p]
		size_arr += [[str(x[2]) for x in group_comp]]
	heatmap_arr = np.array(heatmap_arr)
	size_arr = np.array(size_arr)

	# Generate heatmap
	comparison_labels=['WT vs Var', 'WT vs Compens', 'Var vs Compens']

	cmap = sns.color_palette("vlag", as_cmap=True)
	if not do_RI:
		cmap = sns.color_palette('BrBG_r', as_cmap=True)

	fig, ax = plt.subplots(figsize=(8.5,8.5))
	sns.heatmap(heatmap_arr, annot=size_arr, mask=np.isnan(heatmap_arr), \
		yticklabels=labels, xticklabels=comparison_labels, \
		vmin=-4, vmax=4, fmt='', cmap=cmap, ax=ax)
	ax.tick_params(axis="y", labelsize=9)
	ax.set_facecolor("gray")
	plt.title(plt_title)

	fig_name = 'figures/heatmap/' + plt_title + '_heatmap_perm_mean'
	if do_RI: 
		fig_name += '_RI'
	else:
		fig_name += '_normmRNA'
	# plt.savefig(fig_name + '.png', dpi=300)
	plt.show()


# Plot all heatmaps for all constructs
def plot_all_heatmaps(do_RI=True):
	gene_names = {
		1: "QCR9", 
		2: "RPL28", 
		3: "RPS9B", 
		4: "RPL7A", 
		5: "RPS9A", 
		6: "ACT1", 
		7: "RPS14B",
		8: "RPL36B"
	}

	gene_names = {
		1: "QCR9", 
		2: "RPL28"
	}

	for ii in gene_names.keys():
		idx = str(ii)
		make_heatmap("var_barcode_spliced/S" + idx + "_rd1_var_barcode_spliced.txt", \
			"var_group_files/S" + idx + "_var_groups.txt", gene_names[ii], do_RI=do_RI)


# Count all significantly different variant vs wildtype sequences across all constructs
def get_all_signif_stats(do_RI=True):
	gene_names = {
		1: "QCR9", 
		2: "RPL28", 
		3: "RPS9B", 
		4: "RPL7A", 
		5: "RPS9A", 
		6: "ACT1", 
		7: "RPS14B",
		8: "RPL36B"
	}

	gene_names = {
		1: "QCR9", 
		2: "RPL28"
	}

	all_stem_signif = 0
	all_loop_signif = 0
	all_stem_cnt = 0
	all_loop_cnt = 0

	for ii in gene_names.keys():
		idx = str(ii)
		print(gene_names[ii])
		stem_signif, loop_signif, stem_cnt, loop_cnt = \
			get_signif_stats("var_barcode_spliced/S" + idx + "_rd1_var_barcode_spliced.txt", \
				"var_group_files/S" + idx + "_var_groups.txt", do_RI=do_RI)
		all_stem_signif += stem_signif
		all_loop_signif += loop_signif
		all_stem_cnt += stem_cnt
		all_loop_cnt += loop_cnt

	print("%d significant effects in %d stems" % (all_stem_signif, all_stem_cnt))
	print("%d significant effects in %d loops" % (all_loop_signif, all_loop_cnt))


# plot_all_heatmaps(do_RI=True)
plot_all_heatmaps()

get_all_signif_stats()
