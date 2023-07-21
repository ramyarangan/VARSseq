"""
Barcode sampling used to generate Fig. S16D
"""

import numpy as np 
import random

# Bounded Levenshtein distance
def levenshtein_dist(str1, str2, max_diff=3):
	arr = np.full((len(str1) + 1, len(str2) + 1), max_diff)

	arr[0][0] = 0

	for ii in range(1, len(str1)):
		arr[ii][0] = ii

	for jj in range(1, len(str2)):
		arr[0][jj] = jj

	for ii in range(1, len(str1) + 1):
		min_range = max(1, ii - max_diff)
		max_range = min(ii + max_diff + 1, len(str2) + 1)
		for jj in range(min_range, max_range):
			cost = 1
			if (str1[ii - 1] == str2[jj - 1]):
				cost = 0
			arr[ii][jj] =  min(
				min(arr[ii - 1][jj] + 1, arr[ii][jj - 1] + 1), 
				arr[ii - 1][jj - 1] + cost)
			if arr[ii][jj] > max_diff: 
				arr[ii][jj] = max_diff

	return arr[len(str1)][len(str2)]

def levenshtein_tests():
	# 3
	print(levenshtein_dist("kitten", "sitting"))
	# 3 (max dist)
	print(levenshtein_dist("testkitten", "sitting"))
	# 2
	print(levenshtein_dist("GUAGU", "GAAU"))
	# 3 (max dist)
	print(levenshtein_dist("GUAGUGUAGU", "GAAUGAAU"))

# 30:26:24:20
# T>G>C>A
def get_kmers(k, n):
	rand_vals = np.random.uniform(size=(n, k))
	kmers = np.full((n, k), '')
	kmers[rand_vals <= 0.3] = 'T'
	kmers[(rand_vals > 0.3) & (rand_vals <= 0.56)] = 'G'
	kmers[(rand_vals > 0.56) & (rand_vals <= 0.8)] = 'C'
	kmers[rand_vals > 0.8] = 'A'

	kmers = [''.join(kmer) for kmer in kmers]
	return kmers

def get_all_pairs_dist(k, n):
	kmers = get_kmers(k, n)
	min_dist = k * 2
	for ii in range(len(kmers)):
		print(ii)
		for jj in range(ii + 1, len(kmers)):
			best_dist = levenshtein_dist(kmers[ii], kmers[jj])
			if best_dist < 3:
				print(kmers[ii])
				print(kmers[jj])
			min_dist = min(min_dist, best_dist)
	return min_dist

def stochastic_dist_test(k, n, N):
	kmers = get_kmers(k, n)

	num_dist = [0, 0, 0, 0]
	for ii in range(N):
		if (ii % 10000) == 0:
			print(ii)
		rand1 = int(random.random() * n)
		rand2 = int(random.random() * n)
		while rand2 == rand1: 
			rand2 = int(random.random() * n)
		str1 = kmers[rand1]
		str2 = kmers[rand2]
		cur_dist = levenshtein_dist(str1, str2, max_diff=4)
		if cur_dist < 4:
			num_dist[cur_dist] += 1
	return num_dist

def stochastic_barcode_check(k, n, N):
	kmers = get_kmers(k, n)

	num_passing = {2: 0, 3: 0, 4: 0}

	for ii in range(N):
		print(ii)
		rand_idx = int(random.random() * n)
		str1 = kmers[rand_idx]

		# All passes?
		min_dist = 2 * k
		for jj in range(n):
			if jj == rand_idx: 
				continue
			str2 = kmers[jj]
			min_dist = min(min_dist, levenshtein_dist(str1, str2, max_diff=4))

		if min_dist > 1:
			num_passing[2] += 1
		if min_dist > 2:
			num_passing[3] += 1
		if min_dist > 3: 
			num_passing[4] += 1

	return num_passing

def run_stochastic_dist_test():
	for ii in range(8, 17):
		print("Barcode length: %d" % ii)
		print(stochastic_dist_test(ii, 100000, 5000000))


def run_stochastic_barcode_check():
	for ii in range(8, 17):
		print("Barcode length: %d" % ii)
		num_barcodes = 10000
		print("Num barcodes: %d" % num_barcodes)
		print(stochastic_barcode_check(ii, num_barcodes, 20))


run_stochastic_barcode_check()
