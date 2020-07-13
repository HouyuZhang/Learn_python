# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 15:51:10 2018

@author: Administrator
"""
##############get pattern
import numpy as np

predict_matrix = np.zeros((9,6,20000))
times = np.zeros((9,1),dtype=int)
seg_dict = {}
with open("LSC_CO_motif_9_all_segments.bed") as seg:
    for seg_line in seg:  
        tmp_seg = seg_line.strip().split()
        tmp_key = tmp_seg[0] + '\t' + tmp_seg[1] + '\t' +  tmp_seg[2]
        seg_dict[tmp_key] = tmp_seg[3]

signal_dict = {}
with open("LSC_CO_motif_ref.txt") as signal:
    next(signal)
    for seg_line in signal:  
        tmp_seg = seg_line.strip().split()
        tmp_key = tmp_seg[1].rstrip('"').lstrip('"') + '\t' + tmp_seg[2] + '\t' +  tmp_seg[3]
        signal_dict[tmp_key] = tmp_seg[4:]


for k,v in signal_dict.items():
    if k in seg_dict.keys():
        group = int(seg_dict[k].strip('E')) -1 
        predict_matrix[group,:,times[group,0]] = np.array(v,dtype=int)
        times[group,0] += 1

##############build model
model = np.zeros((9,6,2))
for i in range(9):
    for j in range(6):
        model[i,j,0] = np.mean(predict_matrix[i,j,0:times[i,0]]) - 2*np.std(predict_matrix[i,j,0:times[i,0]])
        model[i,j,1] = np.mean(predict_matrix[i,j,0:times[i,0]]) + 2*np.std(predict_matrix[i,j,0:times[i,0]])
        
def traning(v):
    groups = []
    
    for i in range(9):
        if v[0] >= model[i,0,0] and model[i,0,1] >= v[0]:
            if v[1] >= model[i,1,0] and model[i,1,1] >= v[1]:
                if v[2] >= model[i,2,0] and model[i,2,1] >= v[2]:
                    if v[3] >= model[i,3,0] and model[i,3,1] >= v[3]:
                        if v[4] >= model[i,4,0] and model[i,4,1] >= v[4]:
                                groups.append(i+1)
    return(groups)


out = open("MA9_all_signal_predicted_segments.bed","w+") 
for k,v in signal_dict.items():
    v = np.array(v[0:5],dtype = int)
    groups = traning(v)
    tmp = k.strip().split()
    if k in seg_dict.keys():
        state = seg_dict[k]
    out.write(tmp[0] + '\t' + tmp[1] + '\t' +  tmp[2] +'\t'+str(groups)+ '\t' + state +'\n')

###predict genome
'''
scan_dict = {}
with open("LSC_CO_motif_ref.txt") as scan:
    next(scan)
    for seg_line in scan:  
        tmp_seg = seg_line.strip().split()
        tmp_key = tmp_seg[1].rstrip('"').lstrip('"') + '\t' + tmp_seg[2] + '\t' +  tmp_seg[3]
        scan_dict[tmp_key] = tmp_seg[4:]

''' 










