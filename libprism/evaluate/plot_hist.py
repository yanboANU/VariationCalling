#########################################################################
# File Name: plot_gap_distribution.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Fri 22 Mar 2019 15:13:10 AEDT
#########################################################################
#!/bin/bash
import matplotlib
#import numpy as np
import matplotlib.pyplot as plt
import sys
import math

#use one vector, real histogram

def read_data(filename):
    x = []
    with open(filename) as f:
        for line in f:
            words = line.strip().split()
            x.append(int(words[0]))
    return x



filename=sys.argv[1] # TP_coverage
x = read_data(filename)
filename=sys.argv[2] # FP_coverage
y = read_data(filename)



plt.hist(x, 50, histtype='step', label='TP')
plt.hist(y, 50, histtype='step', label='FP')
plt.legend()

plt.xlabel('21-mer Coverage', fontsize=18)
plt.ylabel('Count', fontsize=18)

#plt.xlim(0,200)
#plt.ylim(10,1000000)
# add a 'best fit' line
#y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
     #np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
#ax.plot(bins, y, '--')
#ax1.set_title('(a) Gap length distribution for PacBio data (18)', fontsize=24)
#ax1.tick_params(labelsize=22)
#ax1.text()

#ax1.text(0, -0.05, '(a)', transform=ax1.transAxes,
#     fontsize=24, fontweight='bold', va='top', ha='right')


'''
pos =[]
for p in patches:
    height = p.get_height()
    #ax.text(p.get_x() + p.get_width() / 2, height, str(int(height*len(data))), ha="center", va="bottom", fontsize=8)
    pos.append(p.get_x() + p.get_width() /2) 
    #print (height)
'''

# Tweak spacing to prevent clipping of ylabel
#fig.tight_layout()
#plt.xticks([],('0','1',))
#plt.xticks(pos, label_list)
plt.subplots_adjust(bottom=0.1, left=0.1)
plt.show()

