"""
Pairwise compare RI fraction between barcode sets in barcode_compare_set_file in violin plots, and 
get significance of each pairwise comparison.
"""
import sys
import numpy as np 
import pandas as pd 
import seaborn as sns
from matplotlib import pyplot as plt 
import scipy.stats as stats

MIN_UMI_COV = 10

barcode_umi_spliced_file = sys.argv[1]
barcode_compare_set_file = sys.argv[2]

def get_barcode_stats(all_barcodes, barcode_umi_spliced_file, verbose=False):
	barcode_stats_dict = {}
	for barcode in all_barcodes:
		barcode_stats_dict[barcode] = np.zeros(4)

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

				barcode_stats = np.array([1, spliced, unspliced, other])
				barcode_stats_dict[barcode] += barcode_stats

	cur_barcode_stats = []

	for barcode in barcode_stats_dict.keys():
		(total_umis, spliced, unspliced, other) = barcode_stats_dict[barcode]
		if total_umis < MIN_UMI_COV:
			continue

		if spliced + unspliced < MIN_UMI_COV:
			continue

		if verbose: 
			print("%s %d %d %d %d" % (barcode, total_umis, spliced, unspliced, other))

		RI_frac = unspliced/(spliced + unspliced)

		cur_barcode_stats += [RI_frac]

	return cur_barcode_stats

def get_all_barcode_stats(barcode_umi_spliced_file, barcode_compare_set_file):
	f = open(barcode_compare_set_file)
	all_barcode_lines = f.readlines()
	f.close()

	var_name_to_barcodes = {}
	var_name = ""
	barcodes = []
	for barcode_line in all_barcode_lines:
		if barcode_line[0] == "V":
			if len(barcodes) > 0:
				var_name_to_barcodes[var_name] = barcodes
			var_name = barcode_line.split(" ")[-1].replace('\n', '')
			barcodes = []
			continue
		barcodes += [barcode_line.replace('\n', '')]
	if len(barcodes) > 0:
		var_name_to_barcodes[var_name] = barcodes

	var_names = []
	all_barcode_stats = []
	for var_name in var_name_to_barcodes:
		verbose = False
		if var_name != 'wt':
			verbose = True
			print(var_name)
		barcode_stats = get_barcode_stats(var_name_to_barcodes[var_name], \
			barcode_umi_spliced_file, verbose=verbose)
		var_names += [var_name]
		all_barcode_stats += [barcode_stats]

	return (var_names, all_barcode_stats)

def make_violin_plot(var_names, all_barcode_stats):
	variant_type = []
	RI_frac = []
	for ii, var_name in enumerate(var_names): 
		variant_type += [var_name] * len(all_barcode_stats[ii])
		RI_frac += all_barcode_stats[ii]

	data_dict = {'variant_type': variant_type, 'RI_frac': RI_frac}

	df = pd.DataFrame(data_dict)

	my_palette = {"wt": sns.color_palette()[8], 
		"bp_mutant": sns.color_palette()[9], 
		"compens": sns.color_palette()[6]
	}
	sns.violinplot(x='variant_type', y='RI_frac', data=df, scale='width', cut=0, gridsize=1000, palette=my_palette)#, color='goldenrod')
	sns.stripplot(x='variant_type', y='RI_frac', data=df, color='black', size=3, \
		jitter=0.2, edgecolors='black')
	plt.savefig('figures/violin/qcr9_branchpoint.png', dpi=300)
	plt.show()

def mean_statistic(x, y, axis):
	return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def get_perm_stat(arr1, arr2):
	perm = stats.permutation_test((arr1, arr2), mean_statistic, n_resamples=100000, \
		vectorized=True)
	return perm

def compare_sets(stats1, stats2):
	stats1 = np.array(stats1)
	stats2 = np.array(stats2)
	mw_bucket = stats.mannwhitneyu(x=(stats1 > 0.2), y=(stats2 > 0.2), alternative='two-sided')
	mw = stats.mannwhitneyu(x=stats1, y=stats2, alternative='two-sided')
	ttest = stats.ttest_ind(stats1, stats2, alternative='two-sided')

	print(mw_bucket)
	print(mw)
	print(ttest)
	print(get_perm_stat(stats1, stats2))

def compare_all_sets(var_names, all_barcode_stats):
	for ii in range(len(var_names)):
		for jj in range(ii + 1, len(var_names)):
			stats1 = all_barcode_stats[ii]
			stats2 = all_barcode_stats[jj]
			print("Var1: %s, Var2: %s" % (var_names[ii], var_names[jj]))
			compare_sets(stats1, stats2)
			print()

(var_names, all_barcode_stats) = get_all_barcode_stats(barcode_umi_spliced_file, barcode_compare_set_file)
print(all_barcode_stats)
compare_all_sets(var_names, all_barcode_stats)
make_violin_plot(var_names, all_barcode_stats)

