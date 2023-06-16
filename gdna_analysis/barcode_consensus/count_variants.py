import sys
from matplotlib import pyplot as plt
import numpy as np 

SEQ_MIN = 285
SEQ_MAX = 305

fa_file = sys.argv[1]
designed_var_file = ""
if len(sys.argv) > 2:
	designed_var_file = sys.argv[2]
	designed_var_barcode_file = sys.argv[3]

f = open(fa_file)
fa_lines = f.readlines()
f.close()

designed_var_counts = {}
designed_var_barcodes = {}
if designed_var_file != "": 
	f = open(designed_var_file)
	designed_var_lines = f.readlines()
	f.close()

	for ii in range(int(len(designed_var_lines)/2)):
		cur_designed_var = designed_var_lines[ii * 2 + 1].replace('\n', '')
		designed_var_counts[cur_designed_var] = 0
		designed_var_barcodes[cur_designed_var] = []

seq_dict = {}
for ii in range(int(len(fa_lines)/2)):
	barcode_cov_line = fa_lines[ii * 2].replace('\n', '').replace('>', '')
	barcode = barcode_cov_line.split(' ')[0]
	cov = barcode_cov_line.split(' ')[1]
	cur_seq = fa_lines[ii * 2 + 1].replace('\n', '')
	if len(cur_seq) < SEQ_MIN or len(cur_seq) > SEQ_MAX:
		continue
	if cur_seq in seq_dict:
		seq_dict[cur_seq] += 1
	else:
		seq_dict[cur_seq] = 1
	if cur_seq in designed_var_counts.keys():
		designed_var_counts[cur_seq] += 1
		designed_var_barcodes[cur_seq] += [(barcode, cov)]

freqs = np.array(list(seq_dict.values()))

ge_1 = sum(freqs > 0)
ge_2 = sum(freqs > 1)
ge_3 = sum(freqs > 2)
ge_5 = sum(freqs > 4)
ge_8 = sum(freqs > 7)
ge_10 = sum(freqs > 9)
print("At least (1, 2, 3, 5, 8, 10) barcodes: (%d, %d, %d, %d, %d, %d)" % \
	(ge_1, ge_2, ge_3, ge_5, ge_8, ge_10))

f = open(designed_var_barcode_file, 'w')

ii = 1
for designed_var, barcode_covs in designed_var_barcodes.items():
	f.write("Variant %d, total barcodes: %d\n" % (ii, len(barcode_covs)))
	for barcode_cov in barcode_covs:
		f.write("%s %s\n" % (barcode_cov[0], barcode_cov[1]))
	ii += 1

f.close()

designed_var_freqs = np.array(list(designed_var_counts.values()))

ge_0 = sum(designed_var_freqs > -1)
ge_1 = sum(designed_var_freqs > 0)
ge_3 = sum(designed_var_freqs > 2)
ge_5 = sum(designed_var_freqs > 4)
ge_10 = sum(designed_var_freqs > 9)
print("At least (0, 1, 3, 5, 10) barcodes (designed vars): (%d, %d, %d, %d, %d)" % \
	(ge_0, ge_1, ge_3, ge_5, ge_10))
