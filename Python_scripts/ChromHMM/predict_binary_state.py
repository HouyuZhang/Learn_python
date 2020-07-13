# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 15:51:10 2018

@author: Administrator
"""
##############get pattern
scan_region = {}
seg_dict = {}
with open("GMP_CO_motif_9_all_segments.bed") as seg:
    for seg_line in seg:  
        tmp_seg = seg_line.strip().split()
        tmp_key = tmp_seg[0] + '\t' + tmp_seg[1] + '\t' +  tmp_seg[2]
        seg_dict[tmp_key] = tmp_seg[3]
        if tmp_seg[3] == 'E9':
            scan_region[tmp_key] = tmp_seg[3]

binary_dict = {}
with open("GMP_CO_binary_ref.txt") as binary:
    next(binary)
    for seg_line in binary:  
        tmp_seg = seg_line.strip().split()
        tmp_key = tmp_seg[1].rstrip('"').lstrip('"') + '\t' + tmp_seg[2] + '\t' +  tmp_seg[3]
        binary_dict[tmp_key] = str(tmp_seg[4:])

pattern_dict = {}
for k,v in binary_dict.items():
    if k not in pattern_dict.keys():
        pattern_dict[v] = seg_dict.get(k)

print('model done!')
###scan group9
out = open("_GMP_all_signal_segments.bed","w+")
     
group9_bin = 0
with open("GMP_all_binary_ref.txt") as all_binary:
    next(all_binary)
    for line in all_binary:
        tmp = line.strip().split()
        chr = tmp[1].rstrip('"').lstrip('"')
        tmp_key = chr + '\t' + tmp[2] + '\t' +  tmp[3]
        #print(tmp_key)
        if tmp_key in scan_region.keys():
            #print('h')
            group9_bin += 1
            out.write(chr + '\t' + tmp[2] + '\t' +  tmp[3] + '\t' + pattern_dict[str(tmp[4:])] + '\n')
print("Work Done !!!  there are " + str(group9_bin) + " bins in group9")










