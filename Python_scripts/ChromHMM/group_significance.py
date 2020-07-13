# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 11:45:44 2018

@author: Administrator
"""
import numpy as np
import scipy.stats as st
import re,random,sys,getopt
import math

if len(sys.argv) < 2:
    print("python sig.py -r MA9_12group_geno_motifID.bed -m MA9_12group_motifID -p MA9_12h_CO.bed -s 10 -g 12 -o ./result.txt")
    sys.exit() 
opts, args = getopt.getopt(sys.argv[1:], 'r:o:m:p:s:g:')
for op, value in opts:
    if op == '-r':
        ref_bins = value
    if op == '-o':
        out = value
    if op == '-m':
        Motifs_ID = value
    if op == '-p':
        Peak = value
    if op == '-s':
        Shuffle = int(value)
    if op == '-g':
        Group = int(value)

ref = np.zeros((Group,9),dtype = int)

motifID_dict = {}
with open(Motifs_ID) as motifID:
    for ID_line in motifID:
        motifID_dict[ID_line.strip()] = 1
                     

ref_dict = {}
with open(ref_bins) as genome_ref:
    for ref_line in genome_ref:
        tmp_ref = ref_line.strip().split()
        key = motifID_dict.get(tmp_ref[3],0)
        tmp_ref.append(key)
        group = int(re.match('E[0-9]{1,2}',tmp_ref[3]).group().lstrip('E'))-1
        ref[group,0] += 1 
        ref[group,1] += key
        ref_dict[tmp_ref[0] + '\t' + tmp_ref[1] + '\t' + tmp_ref[2]] = tmp_ref[3:5]
        
with open(Peak) as peak_motif:
    for peak_line in peak_motif:
        tmp_peak = peak_line.strip().split()
        header = tmp_peak[0] + '\t' + tmp_peak[1] + '\t' + tmp_peak[2]
        group = int(tmp_peak[3].lstrip('E'))-1
        ref[group,2] += 1 
        ref[group,3] += ref_dict[header][1]
        
#shuffle
np_shuffle = np.zeros((Group,int(Shuffle)))
for i in [ 'E' + str(x + 1) + '-' for x in range(Group)]:
    tmp_dict = {}
    group = int(i.lstrip('E').rstrip('-')) - 1
    for k,v in ref_dict.items():
        if v[0].startswith(i):
            tmp_dict[k] = v

    ref[group,4] = ref[group,2]

    for shuff in range(int(Shuffle)) :
        print("Shuffling group" + str(group + 1) + " times " + str(shuff + 1) + '...')
        keys = random.sample(list(ref_dict), ref[group,2])
        values = [ref_dict[k] for k in keys]
        for value in values:
            np_shuffle[group,shuff] += value[1]
    con = st.t.interval(0.999, len(np_shuffle[group,:])-1, loc=np.mean(np_shuffle[group,:]), scale=st.sem(np_shuffle[group,:]))
    ci = list(con)
    if math.isnan(ci[0]):
        ci[0] = 0
    if math.isnan(ci[1]):
        ci[1] = 0
    ave = (ci[0] + ci[1])/2
    if ave == 0:
        ave = 0.00001
    ratio = ((ref[group,3] - ave)/ave)*100
    ref[group,5] = ci[0]
    ref[group,6] = ci[1]
    ref[group,7] = 1000
    ref[group,8] = ratio
np.savetxt(out,ref,delimiter='\t',fmt='%d',header = 'Total group bins\tbin motif\ttotl peak bins\tbin motif\tsample bins\tCI upper\tCIlower\tp-value(rev)\tDFC' )
np.savetxt(ref_bins + Motifs_ID + Peak + "Shuffle_matrix", np_shuffle, delimiter='\t')