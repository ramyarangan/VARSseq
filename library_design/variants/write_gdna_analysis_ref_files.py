import sys
from arnie.mfe import mfe
from arnie.bpps import bpps 
from matplotlib import pyplot as plt 
import numpy as np

variant_tag = sys.argv[1]

variant_file = "final_variants/variant_files_sherlock/" + variant_tag + "_variants.csv"
full_seq_file = "../intron_info/" + variant_tag + "/" + variant_tag + "_fullseq.txt"

variant_fasta = "../../gdna_analysis/ref/" + variant_tag + "_vars.fa"

f = open(full_seq_file)
full_seq_lines = f.readlines()
f.close()

intron_pos = int(full_seq_lines[0].split(" ")[1])
library_pos = int(full_seq_lines[0].split(" ")[2])

f = open(variant_file)
variant_lines = f.readlines()
f.close()


def write_line(f, idx, cur_var):
	f.write('>' + variant_tag + '_' + str(idx))
	f.write('\n')
	f.write('%s\n' % cur_var[(intron_pos - library_pos):])	

f = open(variant_fasta, 'w')

var_count = 1
for variant_line in variant_lines:
	write_line(f, var_count, variant_line.split(",")[0])
	var_count += 1
	if len(variant_line.split(",")[2]) > 0:
		write_line(f, var_count, variant_line.split(",")[2])
		var_count += 1

f.close()