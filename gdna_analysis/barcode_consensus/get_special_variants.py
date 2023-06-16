import sys

SEQ_MIN = 285
SEQ_MAX = 305

fa_file = sys.argv[1]
output_file = sys.argv[2]

QCR9_SEQ = 'CTAAAATGGTATGT'

def get_seq_to_barcodes_dict(fa_file):
	f = open(fa_file)
	fa_lines = f.readlines()
	f.close()

	seq_to_barcodes = {}
	for ii in range(int(len(fa_lines)/2)):
		barcode = fa_lines[ii * 2][1:].replace('\n', '')
		seq = fa_lines[ii * 2 + 1].replace('\n', '')
		if seq in seq_to_barcodes.keys():
			seq_to_barcodes[seq] += [barcode]
		else:
			seq_to_barcodes[seq] = [barcode]
	return seq_to_barcodes

def check_seq(seq, check_condition):
	if check_condition == 'bp':
		if (QCR9_SEQ in seq) and ('TACTAAC' not in seq):
			return True
		return False
	if check_condition == 'fivess':
		if ('GTATGT' not in seq) and ('GTACGT' not in seq):
			return True
		return False
	return False

def get_mutants(seq_to_barcodes, barcode_file, check_condition):
	variants = []
	num_barcodes = 0

	f = open(barcode_file, 'w')
	for seq in seq_to_barcodes.keys():
		if check_seq(seq, check_condition):
			variants += [seq]
			num_barcodes += len(seq_to_barcodes[seq])
			for barcode in seq_to_barcodes[seq]: 
				f.write("%s\n" % barcode)
	f.close()

	print("Number of barcodes: %d" % num_barcodes)
	print("Number of unique variants: %d" % len(variants))

seq_to_barcodes = get_seq_to_barcodes_dict(fa_file)
get_mutants(seq_to_barcodes, output_file, 'bp')