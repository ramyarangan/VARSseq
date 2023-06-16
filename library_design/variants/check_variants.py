import sys
from arnie.mfe import mfe
from arnie.bpps import bpps 
from matplotlib import pyplot as plt 
import numpy as np

variant_file = sys.argv[1]
intron_seq_fasta = sys.argv[2]
intron_full_seq_fasta = sys.argv[3]


def get_seq(fasta_file):
	f = open(fasta_file)
	fasta_lines = f.readlines()
	f.close()

	return fasta_lines[1].replace('\n', '')

def get_intron(variant, library_pos, full_seq, intron_pos, intron_len):
	new_full_seq = list(full_seq)
	new_full_seq[library_pos:(library_pos+len(variant))] = list(variant)
	new_full_seq = ''.join(new_full_seq)
	return new_full_seq[intron_pos:(intron_pos + intron_len)]

intron_seq = get_seq(intron_seq_fasta)
mfe_wt = mfe(intron_seq, package='contrafold')

f = open(intron_full_seq_fasta)
full_seq_lines = f.readlines()
f.close()
library_pos = int(full_seq_lines[0].split(" ")[2])
intron_pos = int(full_seq_lines[0].split(" ")[1])
full_seq = full_seq_lines[1].replace("\n", "")

f = open(variant_file)
variant_lines = f.readlines()
f.close()

# QCR9
check_tags = ["21_28_33_1.0", "21_28_33_37_1.0", "8_199_21_28_33_37_55_58_60_63_69_1.0", "21_28_33_37_55_58_60_63_69_1.0", "55_58_60_63_69_1.0"]
# RPL36B
check_tags = ["31_41_54_58_82_1.0", "31_41_54_58_1.0", "31_41_54_1.0"]
# RPL7A
check_tags = ["47_53_57_63_67_1.0", "47_53_57_63_67_75_1.0", "180_191_194_199_300_1.0", "180_191_194_199_300_209_220_226_235_1.0"]
# RPS14B
check_tags = ["57_71_77_1.0", "99_113_121_126_130_134_0.5", "113_121_126_130_134_1.0", "40_57_71_77_1.0"]

num_variants = 0
for variant_line in variant_lines:
	print(variant_line.split(",")[4])

	if variant_line.split(",")[4] in check_tags:
		bp_matrix = bpps(intron_seq, package='contrafold')
		intron_var = get_intron(variant_line.split(",")[0], library_pos, full_seq, intron_pos, len(intron_seq))
		bp_matrix_var = bpps(intron_var, package='contrafold')
		if len(variant_line.split(",")[2]) > 0:
			intron_var = get_intron(variant_line.split(",")[2], library_pos, full_seq, intron_pos, len(intron_seq))
			bp_matrix_var_compens = bpps(intron_var, package='contrafold')

		min_bpp = min(np.min(bp_matrix), np.min(bp_matrix_var))
		if len(variant_line.split(",")[2]) > 0:
			min_bpp = min(min_bpp, np.min(bp_matrix_var_compens))
		max_bpp = max(np.max(bp_matrix), np.max(bp_matrix_var))
		if len(variant_line.split(",")[2]) > 0:
			max_bpp = max(max_bpp, np.max(bp_matrix_var_compens))

		if len(variant_line.split(",")[2]) > 0:
			fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
			ax1.imshow(bp_matrix, origin='lower', cmap='gist_heat_r', vmin=min_bpp, vmax=max_bpp)
			ax2.imshow(bp_matrix_var, origin='lower', cmap='gist_heat_r', vmin=min_bpp, vmax=max_bpp)
			ax3.imshow(bp_matrix_var_compens, origin='lower', cmap='gist_heat_r', vmin=min_bpp, vmax=max_bpp)
		else:
			fig, (ax1, ax2) = plt.subplots(1, 2)
			ax1.imshow(bp_matrix, origin='lower', cmap='gist_heat_r', vmin=min_bpp, vmax=max_bpp)
			ax2.imshow(bp_matrix_var, origin='lower', cmap='gist_heat_r', vmin=min_bpp, vmax=max_bpp)			
		plt.show()

	print(variant_line.split(",")[0])
	num_variants += 1
	if len(variant_line.split(",")[2]) > 0:
		print(variant_line.split(",")[2])
		num_variants += 1
	print(mfe_wt)
	print(variant_line.split(",")[1])
	if len(variant_line.split(",")[2]) > 0:
		print(variant_line.split(",")[3])
	print()

print(num_variants)

