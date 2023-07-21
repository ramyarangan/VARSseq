"""
Make barplots showing the number of variant sequences with different barcode coverage levels
python plot_variant_coverage.py
"""
from matplotlib import pyplot as plt
import numpy as np

all_variants = {
	'QCR9': [57800, 354, 225, 187], 
	'RPL28': [2600, 216, 197, 183], 
	'RPS9B': [4549, 286, 220, 186], 
	'RPL7A': [5745, 305, 219, 203], 
	'RPS9A': [2175, 213, 189, 150], 
	'ACT1': [184, 9, 5, 2], 
	'RPS14B': [1799, 209, 193, 174], 
	'RPL36B': [2953, 207, 192, 187]
}

designed_variants = {
	'QCR9': [181, 180, 180, 179], 
	'RPL28': [190, 190, 190, 180], 
	'RPS9B': [192, 189, 188, 172], 
	'RPL7A': [191, 191, 191, 190], 
	'RPS9A': [193, 190, 186, 149], 
	'ACT1': [10, 3, 3, 2], 
	'RPS14B': [185, 185, 185, 184], 
	'RPL36B': [192, 191, 190, 187]
}

num_variants = {
	'QCR9': 187, 
	'RPL28': 193, 
	'RPS9B': 193, 
	'RPL7A': 193, 
	'RPS9A': 193, 
	'ACT1': 179, 
	'RPS14B': 185, 
	'RPL36B': 193
}

name_list = all_variants.keys()

fig, axs = plt.subplots(2, 4, figsize=(20,7))
for ii, name in enumerate(name_list):
	axs[int(ii/4), ii%4].bar(["1", "3", "5", "10"], all_variants[name])
	axs[int(ii/4), ii%4].set_title(name)
	#axs[int(ii/4), ii%4].set_xlim((0, 11))	
	#axs[int(ii/4), ii%4].set_xticks([1, 3, 5, 10])
plt.show()

fig, axs = plt.subplots(2, 4, figsize=(20,7))
for ii, name in enumerate(name_list):
	axs[int(ii/4), ii%4].bar(["1", "3", "5", "10"], designed_variants[name])
	axs[int(ii/4), ii%4].axhline(num_variants[name], color='red', linestyle='--')
	axs[int(ii/4), ii%4].set_ylim((0, 200))
	axs[int(ii/4), ii%4].set_title(name)
plt.savefig('figures/designed_variant_coverage.png', dpi=300)
plt.show()