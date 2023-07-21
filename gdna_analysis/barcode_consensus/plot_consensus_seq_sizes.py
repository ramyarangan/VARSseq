"""
Plots histogram of consensus sequence lengths using the ihist files output by bbmerge
plot_consensus_seq_sizes.py consensus_size/ihist1.txt QCR9
"""

import sys
from matplotlib import pyplot as plt

hist_file = sys.argv[1]
plt_title = sys.argv[2]

f = open(hist_file)
hist_lines = f.readlines()
f.close()

hist_lines = hist_lines[6:]

x_vals = [int(x.split('\t')[0]) for x in hist_lines]
y_vals = [int(x.split('\t')[1]) for x in hist_lines]

plt.plot(x_vals, y_vals, '-o', color='black', markersize=3)
plt.title(plt_title)
plt.savefig('figures/'+ plt_title + '_hist.png', dpi=300)
plt.show()