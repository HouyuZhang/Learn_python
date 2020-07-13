# -*- coding: utf-8 -*-
"""
Created on Tue May  1 12:05:43 2018

@author: Administrator
"""  

import getopt,sys

def overlap(s1,e1,s2,e2):
    if not (e1<s2 or e2<s1):
        return (True)
    
def main():
    parameter_list = 'hi:o:e:a:'
    opts, args = getopt.getopt(sys.argv[1:], parameter_list)
    input_file = ''
    output_file = './result.txt'
    for op, value in opts:
        if op == '-i':
            input_file = value
        elif op == '-o':
            output_file = value
        elif op == '-e':
            ext_dis = value
        elif op == '-a':
            anno = value
                                  
    f = open(anno)
    t = open(input_file)
    o = open(output_file,"w+")
    next(f)
    next(t)
    next(t)
    
    print("*" * 100)
    print("Gene list processed...")
    gene_dict = {}
    for line in f:
        line = line.strip()
        line = line.split('\t')
        if line[4]==line[3]:
            line[3] = int(line[3]) + int(ext_dis)
        elif line[4]==line[2]:
            line[2] = int(line[2]) - int(ext_dis)
        gene_dict[line[0]] = line[1]+'\t'+str(line[2])+'\t'+str(line[3])+'\t'+line[4]+'\t'+line[5]+'\t'+'0'
    
    print("Peak list processed...")
    peak_dict = {}
    for line in t:
        line = line.strip()
        line = line.split('\t')
        peak_dict[line[0]] = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\t'+line[6]
    
    print("Start computing overlap...")    
    for key in peak_dict:
        for g_key in gene_dict:
            if peak_dict[key].split('\t')[0] == gene_dict[g_key].split('\t')[0]:
                new_line = gene_dict[g_key].split('\t')
                if overlap(int(new_line[1]),int(new_line[2]),int(peak_dict[key].split('\t')[1]),int(peak_dict[key].split('\t')[2])):
                    new_line[5] = int(new_line[5])+int(peak_dict[key].split('\t')[5])
                    gene_dict[g_key] = new_line[0]+'\t'+new_line[1]+'\t'+new_line[2]+'\t'+new_line[3]+'\t'+new_line[4]+'\t'+str(new_line[5])
                else:
                    pass
            else:
                pass
    o.write("GeneID" + '\t' + input_file + '\n')           
    for k in sorted(gene_dict.keys()):
        o.write(str(k)+'\t'+str(gene_dict[k].split('\t')[5])+'\n')
    o.close()
    print(input_file + "Computed done!  :) ")

if __name__ == '__main__':
	main()






