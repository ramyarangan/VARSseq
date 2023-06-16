import sys

var_file = sys.argv[1]
complete_var_file = sys.argv[2]
add_idx = sys.argv[3]

add_5prime = {
	1: 'AACAATAGCAAGCGGACGATACGGACTAAAATG', 
	2: 'AACTAATCAAAGCGGACGGAATGCCTTCCAGATTCACTAAGACTAGAAAGCACAGAGGTCACGTCTCAG',
	3: 'AACTAAGACTAACTGCCCAAGCAACAATGCCAA',
	4: 'AAGCAAATTAAGAGAAGAGCATCACAATGGCCGCTGA',
	5: 'AATACAAAAGAGCGGACGTATACAACATGCCAA', 
	6: 'ATCGAAAATTAGACGGCGTACTGAATTAACAATGGATTCTG', 
	7: 'AACAACTCCCAGAAGACGAATAACAATTAAGAATGGCTAACG',
	72: 'AACAACTCCCAGAAGACG',
	8: 'GGATAATACAAGCGGACGAAATGGCTGTCAAGACTG'
}

add_3prime = {
	1: '',
	2: '', 
	3: '',
	4: '',
	5: '', 
	6: 'T',
	7: 'G', 
	72: 'G',
	8: ''
}

f = open(var_file)
orig_var_lines = f.readlines()
f.close()

add_idx = int(add_idx)

f = open(complete_var_file, 'w')

for ii in range(int(len(orig_var_lines)/2)):
	f.write("%s" % orig_var_lines[ii * 2])
	new_str = orig_var_lines[ii * 2 + 1].replace('\n', '')
	new_str = add_5prime[add_idx] + new_str + add_3prime[add_idx]
	f.write("%s\n" % new_str)

f.close()

