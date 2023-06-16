import sys

fq_file = sys.argv[1]
tag_blacklist_file = sys.argv[2]
id_to_tag_file = sys.argv[3]
fa_file= sys.argv[4]

f = open(fq_file)
fq_lines = f.readlines()
f.close()

f = open(tag_blacklist_file)
blacklist_tags = f.readlines()
f.close()
blacklist_tags = [x.replace('\n', '') for x in blacklist_tags]

f = open(id_to_tag_file)
id_to_tag_lines = f.readlines()
f.close()

tag_to_id_max = {}
id_to_tag = {}
for id_to_tag_line in id_to_tag_lines:
	cur_id = id_to_tag_line.split(' ')[0]
	cur_tag = id_to_tag_line.split(' ')[1].split(":")[-1].replace('\n', '')
	cur_cov = int(id_to_tag_line.split(' ')[2].split(":")[-1].replace('\n', ''))
	id_to_tag[cur_id] = cur_tag
	if cur_tag in tag_to_id_max.keys():
		cur_max = tag_to_id_max[cur_tag][1]
		if cur_cov > cur_max:
			tag_to_id_max[cur_tag] = (cur_id, cur_cov)
	else:
		tag_to_id_max[cur_tag] = (cur_id, cur_cov)

f = open(fa_file, 'w')

for ii in range(int(len(fq_lines)/4)):
	tag_line = fq_lines[4 * ii]
	seq_line = fq_lines[4 * ii + 1].replace('\n', '')

	id_val = tag_line[1:].replace('\n', '')
	cur_tag = id_to_tag[id_val]

	if cur_tag in blacklist_tags:
		continue
	if tag_to_id_max[cur_tag][0] != id_val:
		continue

	f.write('>%s\n' % cur_tag)
	f.write('%s\n' % seq_line)

f.close()

