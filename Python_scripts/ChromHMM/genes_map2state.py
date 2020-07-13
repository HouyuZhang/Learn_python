# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 19:30:12 2018

@author: Administrator
"""

import getopt,sys
import numpy as np
state = np.zeros((3,20,30))

state_file = "all_chr.state_test"
gene_list = "genelist.TXT"

def cal_state(gene_list,state_file):

	state_dict = {}
	f= open(state_file)
	for bin in f:
		line = bin.strip().split()
		state_dict["\t".join(line[0:3])] = line[3:6]

	gene_dict = {}
	g = open(gene_list)
	for gene_info in g:
		gene = gene_info.strip().split()
		gene_dict[gene[3]] = gene[0:3]

	for num,key in enumerate(gene_dict.keys()):
		coor = gene_dict[key]
		start = int(coor[1]) - int(coor[1])%200
		end = int(coor[2]) - int(coor[2])%200
		bin_num = int((end - start)/200) -1
		for i in range(bin_num):
			new_term = coor[0] + "\t" + str(start + i*200) + "\t" + str(start + i*200 + 200)
			state[0,num,int(state_dict[new_term][0])] += 1
			state[1,num,int(state_dict[new_term][1])] += 1
			state[2,num,int(state_dict[new_term][2])] += 1
		state[0,num,:] = state[0,num,:]/sum(state[0,num,:])
		state[1,num,:] = state[1,num,:]/sum(state[1,num,:])
		state[2,num,:] = state[2,num,:]/sum(state[2,num,:])


	np.savetxt("B.txt",state[0], delimiter="\t", fmt='%.2f', header='\t'.join(['group' + str(x) for x in list(range(1,31))]))
	np.savetxt("CD4.txt",state[1], delimiter="\t", fmt='%.2f',header='\t'.join(['group' + str(x) for x in list(range(1,31))]))
	np.savetxt("CD8.txt",state[2], delimiter="\t", fmt='%.2f',header='\t'.join(['group' + str(x) for x in list(range(1,31))]))

def main():
	opts, args = getopt.getopt(sys.argv[1:], 'i:s:')
	for op, value in opts:
		if op == '-i':
			gene_list = value
		if op == '-s':
			state_file = value  

	cal_state(gene_list,state_file)


if __name__ == '__main__':
	main()