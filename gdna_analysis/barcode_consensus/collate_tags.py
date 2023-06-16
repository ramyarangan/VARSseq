import sys

FILTER_AMOUNT=0.05

consensus_tag_file = sys.argv[1]
collated_tag_file = sys.argv[2]
tag_blacklist_file = sys.argv[3]

f = open(consensus_tag_file)
consensus_tag_lines = f.readlines()
f.close()

consensus_tag_lines = consensus_tag_lines[3:]

tag_to_info = {}

for tag_line in consensus_tag_lines:
	tag_items = tag_line.split(' ')
	
	if tag_items[1] == '141':
		continue

	tag = tag_items[-1].split(':')[-1].replace('\n', '')
	max_val = tag_items[2].split(':')[-1]
	min_val = tag_items[3].split(':')[-1]
	id_val = tag_items[0][1:]

	if tag in tag_to_info.keys():
		tag_to_info[tag] += [(id_val, min_val, max_val)]
	else:
		tag_to_info[tag] = [(id_val, min_val, max_val)]

num_tags = len(tag_to_info.keys())

f = open(collated_tag_file, 'w')
blacklist_tags = []

for tag, tag_list in tag_to_info.items():
	f.write("%s %d\n" % (tag, len(tag_list)))
	max_vals = []
	for tag_item in tag_list:
		max_vals += [int(tag_item[-1])]
		f.write("%s %s %s\n" % tag_item)
	max_max_vals = max(max_vals)
	max_vals.remove(max_max_vals)
	num_other_reads = sum(max_vals)
	if num_other_reads/max_max_vals > FILTER_AMOUNT:
		blacklist_tags += [tag]
f.close()

f = open(tag_blacklist_file, 'w')
for tag in blacklist_tags:
	f.write("%s\n" % tag)
f.close()

print("Fraction blacklisted: %d of %d = %f" % (len(blacklist_tags), num_tags, \
	len(blacklist_tags)/num_tags))

