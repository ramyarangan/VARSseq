"""
Compare wildtype, variant, and compensatory variant retained intron fraction and spliced mRNA
levels. Print out significance values by permutation test for significant comparisons and make violin plots.

To compare retained intron fractions: 
python plot_signif_variants.py var_barcode_spliced/S5_rd1_var_barcode_spliced.txt var_group_files/S5_var_groups.txt RPS9A RI

To compare normalized spliced mRNA levels: 
python plot_signif_variants.py var_barcode_spliced/S5_rd1_var_barcode_spliced.txt var_group_files/S5_var_groups.txt RPS9A mRNA

"""

import sys
import scipy.stats as stats
import pandas as pd 
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt

# Only include barcodes that have at least 10 UMI's
MIN_UMI_COV = 10

var_barcode_spliced_file = sys.argv[1] # Per-variant per-barcode spliced / unspliced / other stats
var_group_file = sys.argv[2] # Sequences for variants that belong to the same stem set disruption / rescue set
construct_tag = sys.argv[3] # Current construct analyzed (used to select violin plots to show)
do_RI_tag = sys.argv[4] # Whether to compare RI or mRNA (must be either "RI" or "mRNA")

do_RI = None
if do_RI_tag == "RI":
	do_RI = True
elif do_RI_tag == "mRNA":
	do_RI = False
else:
	raise RuntimeError("Must compare with mRNA or RI")

# Get the per variant and wildtype lists of per-barcode RI or mRNA stats
def get_var_wt_stats(var_barcode_spliced_file, do_RI=True):
	f = open(var_barcode_spliced_file)
	var_barcode_spliced_lines = f.readlines()
	f.close()

	all_var_stats = {}

	cur_barcode_stats = []
	cur_id = -1
	for var_line in var_barcode_spliced_lines:
		if var_line[0] == "V":
			if cur_id > -1:
				all_var_stats[cur_id] = np.array(cur_barcode_stats)
			cur_id = int(var_line.split(" ")[-1].replace('\n', ''))
			cur_barcode_stats = []
			continue

		# Skip barcodes with low coverage
		total_umis = int(var_line.split(' ')[0])
		if total_umis < MIN_UMI_COV:
			continue
		
		# Get retained intron fraction and normalized spliced mRNA level
		spliced = int(var_line.split(' ')[1])
		unspliced = int(var_line.split(' ')[2])
		if spliced + unspliced < MIN_UMI_COV:
			continue
		RI_frac = unspliced/(spliced + unspliced)
		spliced_norm_gDNA = spliced/int(var_line.split(' ')[4])
		if do_RI:
			cur_barcode_stats += [RI_frac]
		else:
			cur_barcode_stats += [spliced_norm_gDNA]

	wt_stats = np.array(cur_barcode_stats)

	return all_var_stats, wt_stats


# Get lists of per-barcode stats by stem set / junction set
def get_group_stats(all_var_stats, wt_stats, var_group_file):
	group_stats = {}

	f = open(var_group_file)
	var_group_lines = f.readlines()
	f.close()

	ii = 0
	while ii < len(var_group_lines):
		# Collect all the variant sequence ID's for a given stem / junction set
		group_id_line = var_group_lines[ii]
		group_id = (group_id_line.split(' ')[0], group_id_line.split(' ')[1])
		num_vars = int(group_id_line.split(' ')[2].replace('\n', ''))
		var_ids = []
		for jj in range(num_vars):
			ii += 1
			var_id = int(var_group_lines[ii].split(' ')[0])
			var_ids += [var_id]

		unique_var_ids = list(set(var_ids))

		# Get all the RI or mRNA stats per stem / junction set
		cur_group_stats = []
		for var_id in unique_var_ids:
			if var_id in all_var_stats.keys():
				cur_group_stats += list(all_var_stats[var_id])
		group_stats[group_id] = np.array(cur_group_stats) 
		var_ids = []
		ii += 1

	return group_stats


def mean_statistic(x, y, axis):
	return np.mean(x, axis=axis) - np.mean(y, axis=axis)

# Permutation test between two arrays with the mean statistic
def get_perm_stat(arr1, arr2):
	perm = stats.permutation_test((arr1, arr2), mean_statistic, n_resamples=100000, \
		vectorized=True)
	return perm

def compare_sets(arr1, arr2):
	return get_perm_stat(arr1, arr2)


# Compare variants with wildtype sequences and compensatory rescue
# sequences and print out significance of comparisons
def compare_group_stats(group_stats, wt_stats):
	group_id_pvals = {}
	collected_pvals = {}
	for group_id in group_stats.keys():
		barcode_stats = group_stats[group_id]
		if len(barcode_stats) < 2:
			continue
		perm = compare_sets(wt_stats, barcode_stats)

		# print(group_id)
		# print("Num variant barcodes: %d" % len(barcode_stats))

		if group_id[1] != 'variant':
			continue

		collected_pvals[('wt', group_id)] = perm.pvalue

		group_id_pvals[group_id] = perm.pvalue

		compens_group_id = (group_id[0], 'compens')
		if compens_group_id not in group_stats.keys():
			continue

		compens_stats = group_stats[compens_group_id]
		if len(compens_stats) < 2:
			continue		
		perm = compare_sets(barcode_stats, compens_stats)
		collected_pvals[(group_id, compens_group_id)] = perm.pvalue
		
		# print("Has compens")
		# print("Num compens barcodes: %d" % len(compens_stats))
		# print(compens_group_id)
		print_str = "(%.4g)" % perm.pvalue
		group_id_pvals[compens_group_id] = perm.pvalue

	for group_id in group_stats.keys():
		if group_id[1] == "compens":
			continue
		do_print = False
		if group_id in group_id_pvals.keys():
			print_str = group_id[0] + ": "
			if group_id_pvals[group_id] < 0.05: 
				do_print = True
			print_str += "%.4g" % group_id_pvals[group_id]
			compens_group_id = (group_id[0], 'compens')
			if compens_group_id in group_id_pvals.keys():
				print_str += "; " 
				print_str += "%.4g" % group_id_pvals[compens_group_id]
				if group_id_pvals[compens_group_id] < 0.05: 
					do_print = True
			if do_print: 
				print(print_str)


# Make violin plots comparing wildtype, variant, and compensatory rescue stats
def make_violin_plot(wt_stats, var_stats, violin_title, compens_stats=[], do_RI=True):
	variant_type = ['wt'] * len(wt_stats) + ['var'] * len(var_stats) + ['compens'] * len(compens_stats)
	RI_frac = list(wt_stats) + list(var_stats) + list(compens_stats)

	label = 'RI_frac'
	if not do_RI:
		label = 'norm_mRNA'
	data_dict = {'variant_type': variant_type, label: RI_frac}

	df = pd.DataFrame(data_dict)
	my_palette = {"wt": sns.color_palette()[8], 
		"var": sns.color_palette()[9], 
		"compens": sns.color_palette()[6]
	}
	sns.violinplot(x='variant_type', y=label, data=df, scale='width', cut=0, gridsize=1000, palette=my_palette)#, color='goldenrod')
	sns.stripplot(x='variant_type', y=label, data=df, color='black', size=3, \
		jitter=0.2, edgecolors='black')
	plt.savefig('figures/violin/' + violin_title + '.png', dpi=300)
	plt.show()


# Make all violin plots of interest across tested constructs
def make_violin_plots(wt_stats, group_stats, construct_tag, do_RI=True):
	tag_to_plot_lists = {
		"QCR9": [('8_199', True), ('113', True), ('33_37', True)], 
		"RPL7A": [('63_67', False), ('47_53_57_63_67', False), ('75', False), \
			('100', True), ('135', True), ('135_140_new', True), ('135_140', True), \
			('47_53_57_63_67_75', False)], 
		"RPS9B": [('69', False), ('80', False), ('69_80', False),
			('89_95', False), ('69_80_89_95', False), 
			('149', False), ('201_180', True), ('167_201_180', True), ('186', False)], 
		"RPL36B": [('54_58', True), ('31_41_54_58_82', True), 
			('31_41_54_58', True), ('31_41_54_58_82_comb', True)],
		"RPS9A": [('55', False), ('77_80_86_88', False), ('77', False), ('77_80', False), \
			('103_112', False), ('193', True), ('185_193', True), \
			('182_185_193', True), ('174_182_185_193', True), ('200', False)],
		"RPL28": [('29_34', False), ('29_34_44_55_63', False), ('29_34_44_55', False), 
			('44_55', False), ('44', False), ('55', False), ('40', False), ('76_81_93', False),
			('63', False), ('105_109', False), ('138_141', False),
			('76_81_93_100', False), ('76_81_93_100_105_109_116_120', False)],
		"RPS14B": [('40', False), ('57', False), ('71_77', False), \
			('57_71_77', False), ('88_99', True)]
	}

	for tag_to_plot in tag_to_plot_lists[construct_tag]:
		violin_title = construct_tag + '_' + tag_to_plot[0]
		if tag_to_plot[1]:
			make_violin_plot(wt_stats, group_stats[(tag_to_plot[0], 'variant')], \
				violin_title, compens_stats=group_stats[(tag_to_plot[0], 'compens')], do_RI=do_RI)
		else:
			make_violin_plot(wt_stats, group_stats[(tag_to_plot[0], 'variant')], violin_title, do_RI=do_RI)	

all_var_stats, wt_stats = get_var_wt_stats(var_barcode_spliced_file, do_RI=do_RI)
group_stats = get_group_stats(all_var_stats, wt_stats, var_group_file)
compare_group_stats(group_stats, wt_stats)
make_violin_plots(wt_stats, group_stats, construct_tag, do_RI=do_RI)

