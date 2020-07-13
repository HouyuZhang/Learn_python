# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:01:08 2018

@author: Administrator
"""
import numpy as np
import getopt,sys

if len(sys.argv) < 2:
	print("python bin_peak_overlap_1~8.py -i WT_12_CO_motif.bed -g LSC_group1~8.bed")
	print("better use bedtools intersect instead!!!")
	sys.exit() 
opts, args = getopt.getopt(sys.argv[1:], 'i:g:')
for op, value in opts:
	if op == '-i':
		file = value
	if op == '-g':
		peak = value


seg_dict = {}

def overlap(a0,a1,a2,b0,b1,b2):
	if a0 == b0:
		if not (int(b1) + 100 > int(a2) or int(b2) - 100 < int(a1)):
			return True

seg_dict = {}
with open (peak) as peak_file:
	next(peak_file)
	for pf_line in peak_file:
		tmp_pf = pf_line.strip().split()
		seg_dict[tmp_pf[0] + '\t' + tmp_pf[1] + '\t' +  tmp_pf[2]] = tmp_pf[3]

peak_dict = {}
with open(file) as seg:
	for seg_line in seg:  
		tmp_seg = seg_line.strip().split()
		tmp_key = tmp_seg[0] + '\t' + tmp_seg[1] + '\t' +  tmp_seg[2]
		peak_dict[tmp_key] = tmp_seg[3]



#out = open(peak + file,"w+")
#res = np.zeros((30,1))
group_dict = {}
#key = group_dict.keys()
for k,v in seg_dict.items():
	for peak_key,peak_num in peak_dict.items():
		a0,a1,a2 = peak_key.split()[0:3]
		b0,b1,b2 = k.split()[0:3]
		if overlap(a0,a1,a2,b0,b1,b2):
			if v in group_dict.keys():
				group_dict[v] +=1
			else:
				group_dict[v] = 1 
			#out.write(v + '\t' + a0 + '\t' + a1 + '\t' + a2 + '\t' + peak_num + '\n')
			break
for k,v in sorted(group_dict.items()):
	print(k,v)
print(file + " done!\n")

