# -*- coding: utf-8 -*-
import getopt,sys

#
if len(sys.argv) < 2:
    print("This script is for seperate ChromHMM result to 200bp bins, because ChromHMM will merge adjacent bins in same state.")
    print("USAGE: python bin_seperate.py -i LSC_CO_motif_9_segments.bed -o ***")
    sys.exit() 
    
opts, args = getopt.getopt(sys.argv[1:], 'i:o:')
for op, value in opts:
    if op == '-i':
        file = value
    if op == '-o':
        out = value

out = open(out,"w+")
seg_dict = {}
with open(file) as seg:
    for seg_line in seg:  
        tmp_seg = seg_line.strip().split()
        for i in range(int((int(tmp_seg[2])-int(tmp_seg[1]))/200)):
            out.write(tmp_seg[0] + '\t' +str(int(tmp_seg[1]) + i*200) + '\t' +  str(int(tmp_seg[1]) + i*200 + 200) + '\t' + tmp_seg[3] + '\n')
            
print("DONE!")