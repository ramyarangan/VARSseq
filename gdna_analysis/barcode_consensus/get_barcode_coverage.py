import sys

collated_tags_file = sys.argv[1]
consensus_fa_file = sys.argv[2]
consensus_cov_file = sys.argv[3]

f = open(collated_tags_file)
collated_tag_lines = f.readlines()
f.close()

f = open(consensus_fa_file)
consensus_fa_lines = f.readlines()
f.close()

barcode_cov_dict = {}
ii = 0
cur_total = 0
while ii < len(collated_tag_lines):
	cur_line = collated_tag_lines[ii]
	cur_barcode = cur_line.split(' ')[0]
	num_seqs = int(cur_line.split(' ')[1].replace('\n', ''))
	ii += 1

	cur_total = 0
	for jj in range(num_seqs):
		cur_total += int(collated_tag_lines[ii].split(' ')[-1].replace('\n', ''))
		ii += 1

	barcode_cov_dict[cur_barcode] = cur_total

f = open(consensus_cov_file, 'w')

for ii in range(int(len(consensus_fa_lines)/2)):
	barcode_line = consensus_fa_lines[ii * 2].split(">")[1].replace('\n', '')
	coverage = barcode_cov_dict[barcode_line]
	f.write('>%s %d\n' % (barcode_line, coverage))
	f.write('%s' % consensus_fa_lines[ii * 2 + 1])

f.close()