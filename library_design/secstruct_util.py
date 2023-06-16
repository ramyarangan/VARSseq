import os
import networkx as nx 
import numpy as np

def get_seq(fasta_file):
	f = open(fasta_file)
	fasta_lines = f.readlines()
	f.close()

	return fasta_lines[1].replace('\n', '')

def get_dotbracket(dotbracket_filename):
	if not os.path.exists(dotbracket_filename):
		return None

	f = open(dotbracket_filename)
	lines = f.readlines()
	f.close()
	return lines[0].replace('\n', '')

def get_bpp_matrix(bpp_matrix_filename):
	if not os.path.exists(bpp_matrix_filename):
		return None 

	f = open(bpp_matrix_filename)
	lines = f.readlines()
	f.close()
	
	bpp_matrix = []
	for line in lines: 
		bpp_matrix_line = line.split(',')
		bpp_matrix_line = [float(x) for x in bpp_matrix_line]
		bpp_matrix += [bpp_matrix_line]
	
	return bpp_matrix

def get_bp_matrix(mfe):
	bp_matrix = np.zeros((len(mfe), len(mfe)))

	bp_open_stack = []
	for ii, cur_char in enumerate(mfe):
		if cur_char == ".":
			continue
		if cur_char == "(":
			bp_open_stack += [ii]
		if cur_char == ")":
			bp_matrix[bp_open_stack[-1]][ii] = 1
			bp_matrix[ii][bp_open_stack[-1]] = 1
			bp_open_stack = bp_open_stack[:-1]

	return bp_matrix

def get_node_lists(node_list_file):
	f = open(node_list_file)
	node_list_items = f.readlines()
	f.close()

	node_lists = []
	for node_list_item in node_list_items:
		node_list = node_list_item.replace('\n', '').split(' ')
		perc_change = float(node_list[-1])

		node_list = [int(x) for x in node_list[:-1]]
		node_lists += [(node_list, perc_change)]

	return node_lists

def get_nodelist_tag(node_list, perc_change):
	node_list = [str(x) for x in node_list]
	perc_change_str = '%s' % float('%.1g' % perc_change)
	nodelist_tag = 	'_'.join(node_list) + '_' + perc_change_str
	return nodelist_tag

def get_nodelist_filename(node_list, perc_change, mfe_output_pref):
	nodelist_tag = get_nodelist_tag(node_list, perc_change)
	nodelist_filename = mfe_output_pref + nodelist_tag + '.txt'
	return nodelist_filename

class BigNode():
	def get_type(self):
		return 'base'

class Stem(BigNode):
	def __init__(self, strand1_nts, strand2_nts, bpp_matrix):
		if len(strand1_nts) != len(strand2_nts):
			raise RuntimeError("Stem needs to have equal # nts in strands")

		self.strand1_nts = strand1_nts
		self.strand2_nts = strand2_nts
		self.bpp = self.get_bpp(bpp_matrix)

	def get_bpp(self, bpp_matrix):
		bpp = 0
		for ii, nt1 in enumerate(self.strand1_nts):
			nt2 = self.strand2_nts[ii]
			bpp = max(bpp, bpp_matrix[nt1][nt2])
		return bpp

	def get_type(self):
		return 'stem'

	def len(self):
		return len(self.strand1_nts)

	def get_seq_pos_to_strand1_pos(self):
		seq_pos_to_strand_pos = {}
		for ii, pos in enumerate(self.strand1_nts):
			seq_pos_to_strand_pos[pos] = ii
		return seq_pos_to_strand_pos

	def get_pairs(self):
		pairs = []
		for ii, nt in enumerate(self.strand1_nts):
			pairs += [(nt, self.strand2_nts[ii])]
		return pairs

	def __str__(self):
		print_str = "Stem of length %d containing base-pairs:\n" % len(self.strand2_nts)
		for ii, n1 in enumerate(self.strand1_nts):
			print_str += "%d %d\n" % (n1, self.strand2_nts[ii])
		print_str += "Base pair probability: %f\n" % self.bpp
		return print_str

class Junction(BigNode):
	def __init__(self, nts):
		self.nts = nts

	def get_type(self):
		return 'junction'

class External(BigNode):
	def __init__(self, nts):
		self.nts = nts

	def get_type(self):
		return 'external'


# Update graph with nodes for (cur_start, cur_end - 1) (inclusive)
# nt_to_node_dict: for each nucleotide, the stem, junction, or external node it belongs to
# dotbracket: dot bracket notation for the full structure
# G: current graph
# is_internal: False for a portion of the structure that is not in any stems 
#              or internal to any stems
def get_stem_graph_rec(cur_start, cur_end, nt_to_node_dict, \
	dotbracket, G, is_internal, verbose=False):
	# Add node for current portion of stem graph
	G.add_node(cur_start)

	if dotbracket[cur_start] == ')':
		raise RuntimeError("Unexpected recursion architecture: start index is a closing base-pair")

	# Case 1: stems. Here cur_start and cur_end are the ends of a stem; recursively update internal loops
	if dotbracket[cur_start] == '(':
		if verbose: 
			print("Doing stem at position: %d" % cur_start)
		if cur_start not in nt_to_node_dict.keys() or \
			nt_to_node_dict[cur_start].get_type() != 'stem':
			raise RuntimeError("Expected base-paired residue to be in a stem.")

		cur_stem = nt_to_node_dict[cur_start]
		
		stem_start = min(cur_stem.strand1_nts)
		junc_start = max(cur_stem.strand1_nts)
		junc_end = min(cur_stem.strand2_nts)
		stem_end = max(cur_stem.strand2_nts)
		
		if stem_start != cur_start:
			raise RuntimeError("Unexpected recursion architecture")

		# Update graph for nts between the 5' and 3' strands of stem
		# E.g. when Called on (((....)))... from (((((....)))...)), get node for ....
		node_id = get_stem_graph_rec(junc_start + 1, junc_end, nt_to_node_dict, \
			dotbracket, G, True, verbose=verbose)
		# Connect node for .... to ((()))
		G.add_edge(stem_start, node_id)
		
		node_id_2 = cur_start
		if cur_end > stem_end + 1:
			# This happens for internal loops with only 3' nucleotides or
			# external ssRNA outside of stems 
			# E.g. when Called on (((....)))... from (((((....)))...)), get node for ...
			node_id_2 = get_stem_graph_rec(stem_end + 1, cur_end, nt_to_node_dict, \
				dotbracket, G, is_internal, verbose=verbose)
			# Connect node for ... to ((()))
			G.add_edge(stem_start, node_id_2)

		if verbose:
			print("Finished stem at position: %d" % cur_start)

		if is_internal:
			# Called on (((....)))... from (((((....)))...))
			# Should return node for ... 
			# Called on (((....))) from ((...(((....)))...)) 
			# Should return node for ((()))
			return node_id_2
		# Called on (((....)))..... from ...(((....))).....
		# Should return node for ((()))
		return cur_start

	# all ssRNA in the loop / junctions / external nucleotides
	junc_nts = []

	# Case 2: external ssRNA
	if not is_internal:
		if verbose:
			print("Doing external ssRNA at position: %d" % cur_start)
		ssrna_end = cur_end
		cur_dotbracket = dotbracket[cur_start:cur_end]
		if '(' in cur_dotbracket: 
			ssrna_end = cur_dotbracket.index('(') + cur_start
			if ')' not in cur_dotbracket: 
				raise RuntimeError("Unexpected stem architecture")
			if dotbracket.index(')') < dotbracket.index('('):
				raise RuntimeError("Unexpected stem architectures")

		junc_nts = list(range(cur_start, ssrna_end))

		if cur_end > ssrna_end:
			# E.g. when called on ....(((...))) from ((....))....(((...)))
			node_id = get_stem_graph_rec(ssrna_end, cur_end, nt_to_node_dict, \
				dotbracket, G, False, verbose=verbose)
			# Connect .... to ((()))
			G.add_edge(cur_start, node_id)

		if verbose:
			print("Finished external ssRNA at position: %d" % cur_start)

		# Add new junction to datastructure
		new_external = External(junc_nts)
		for ii in junc_nts:
			nt_to_node_dict[ii] = new_external	

	# Case 3: internal ssRNA - loops and junctions
	else:
		if verbose:
			print("Doing internal loop / junction at position: %d" % cur_start)

		junc_nts = [] 
		neighbor_stems = [] # all stems internal to the junction
		cur_idx = cur_start
		while (cur_idx < cur_end):
			if dotbracket[cur_idx] == '.':
				junc_nts += [cur_idx]
				cur_idx += 1
			else:
				if cur_idx not in nt_to_node_dict.keys() or \
					nt_to_node_dict[cur_idx].get_type() != 'stem':
					raise RuntimeError("Expected base-paired residue to be in a stem.")
				neighor_stem = nt_to_node_dict[cur_idx]
				neighbor_stems += [neighor_stem]
				
				# Skip over all nucleotides in the new stem
				cur_idx = max(neighor_stem.strand2_nts) + 1

		# Recursion over all stems in the junction
		for neighbor_stem in neighbor_stems:
			start_idx = min(neighbor_stem.strand1_nts)
			end_idx = max(neighbor_stem.strand2_nts)
			# E.g. when called on ....((...))...(....).. from ((....((...))...(....)..))
			node_id = get_stem_graph_rec(start_idx, end_idx + 1, nt_to_node_dict, \
				dotbracket, G, True, verbose=verbose)
			# Connect ......... to (()) and ......... to ()
			G.add_edge(cur_start, node_id)

		if verbose:
			print("Finished internal loop / junction at position: %d" % cur_start)

		# Add new junction to datastructure
		new_junction = Junction(junc_nts)
		for ii in junc_nts:
			nt_to_node_dict[ii] = new_junction

	# When called on ....((...))...(....).. from ((....((...))...(....)..)) 
	# return start of .........
	# When called on ....(((...))) from ((....))....(((...)))
	# return start of ....
	return cur_start

# Node types: 
#      Stem (residues, length, bootstrapping probability) 
#      Loop (all loop nts condensed into 1)
# Used for: NWJ accounting, longest stem
# Return type: networkx graph, dictionary connecting nucleotides to the stems/ junctions
#              they belong in 
def get_stem_graph(dotbracket, bpp_matrix, stem_verbose=False, graph_verbose=False):
	stems = []

	# First build list of stems
	base1_list = []
	cur_stem = []
	nt_to_stem_dict = {}
	for ii, curchar in enumerate(dotbracket):
		in_stem = False
		if len(cur_stem) > 0 and curchar == ')':
			# We could be continuing a stem
			if base1_list[-1] == cur_stem[-1][0] - 1:
				in_stem = True

		if not in_stem and len(cur_stem) > 0:
			strand1_nts = [x[0] for x in cur_stem]
			strand2_nts = [x[1] for x in cur_stem]
			new_stem = Stem(strand1_nts, strand2_nts, bpp_matrix)
			if stem_verbose:
				print(new_stem)
			stems += [new_stem]
			for nt1 in strand1_nts:
				nt_to_stem_dict[nt1] = new_stem
			for nt2 in strand2_nts:
				nt_to_stem_dict[nt2] = new_stem
			cur_stem = []

		if curchar == '(':
			base1_list += [ii]

		if curchar == ')':
			cur_stem += [(base1_list[-1], ii)]
			base1_list = base1_list[:-1]

	if len(cur_stem) > 0:
		strand1_nts = [x[0] for x in cur_stem]
		strand2_nts = [x[1] for x in cur_stem]
		new_stem = Stem(strand1_nts, strand2_nts, bpp_matrix)
		stems += [new_stem]		
		for nt1 in strand1_nts:
			nt_to_stem_dict[nt1] = new_stem
		for nt2 in strand2_nts:
			nt_to_stem_dict[nt2] = new_stem

	# Initiate stem graph building recursion
	G = nx.Graph()
	nt_to_node_dict = nt_to_stem_dict
	get_stem_graph_rec(0, len(dotbracket), nt_to_node_dict, dotbracket, G, \
		False, verbose=graph_verbose)
	return G, nt_to_node_dict, stems

def get_stems_junctions_from_nodelist(node_list, nt_to_node_dict, bpp_matrix):
	stems = []
	junctions = []
	seen_nodes = []
	
	for node_nt in node_list:
		cand_node = nt_to_node_dict[node_nt]
		if cand_node.get_type() == "junction":
			if cand_node.nts[0] in seen_nodes:
				continue
			seen_nodes += [cand_node.nts[0]]
			junctions += [cand_node]

		if cand_node.get_type() == "stem":
			if cand_node.strand1_nts[0] in seen_nodes:
				continue
			if cand_node.get_bpp(bpp_matrix) < 0.7:
				print("Warning: stem is low confidence.")
			seen_nodes += [cand_node.strand1_nts[0]]
			stems += [cand_node]

	return stems, junctions

