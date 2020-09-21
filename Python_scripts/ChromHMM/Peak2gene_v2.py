#!/usr/bin/env python
# -*- coding: utf-8 -*-

Description = '''
This scripts can map Peak tags to given region on genome.
**NOTE:
1. Peak count can be calculated using FeatureCount.
Args: featureCounts -a file.saf -F SAF -p --minOverlap 15 -o file_peakcount file.bam  
2. Gene annotation must in bed format follows: Chr\tStart\tEnd\tTSS\tName\t...
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2019-04-11.
Issue report on Hughiez047@gmail.com
Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import argparse,sys

def overlap(s1,e1,s2,e2):
    if not (e1<s2 or e2<s1):
        return (True)
    
def run(options):
    a = open(options.Gene_annotation)
    p = open(options.Peak_count)
    o = open(options.gene_count_file,"w+")
    next(a)
    next(p)
    next(p)
    
    print(">" * 50)
    print("Gene list processing...")
    gene_dict = {}
    for line in a:
        line = line.strip().split('\t')
        if(options.Ext_feature == "TSS"):
            if line[3]==line[2]:
                ss = int(line[2]) - options.Ext_size
                ee = int(line[2]) + options.Ext_size
            elif line[3]==line[1]:
                ss = int(line[1]) - options.Ext_size
                ee = int(line[1]) + options.Ext_size
            gene_dict[line[4]] = line[0]+'\t'+str(ss)+'\t'+str(ee)+ '\t' + '0'
        elif(options.Ext_feature == "TSS_genebody"):
            if line[3]==line[2]:
                ss = int(line[1])
                ee = int(line[2]) + options.Ext_size
            elif line[3]==line[1]:
                ss = int(line[1]) - options.Ext_size
                ee = int(line[2])
            gene_dict[line[4]] = line[0]+'\t'+str(ss)+'\t'+str(ee) + '\t' +'0'
            
    print("Peak list processing...")
    peak_dict = {}
    for line in p:
        line = line.strip().split('\t')
        peak_dict[line[0]] = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[6]
    
    print("Start computing overlap...")    
    for key in peak_dict:
        pl = peak_dict[key].split('\t')
        for g_key in gene_dict:
            #compare chrosome 
            gl = gene_dict[g_key].split('\t')
            if pl[0] == gl[0]:
                if overlap(int(gl[1]),int(gl[2]),int(pl[1]),int(pl[2])):
                    gl[3] = int(gl[3]) + int(pl[3])
                    gene_dict[g_key] = ('\t').join(gl[0:3]) + '\t' + str(gl[3])
                else:
                    pass
            else:
                pass
    
    #o.write("GeneID" + '\t' + options.Peak_count + '\n')           
    for k in sorted(gene_dict.keys()):
        o.write(str(k)+'\t'+str(gene_dict[k].split('\t')[3])+'\n')
    o.close()
    print(options.Peak_count + " Computed done! :) ")
    

def main():   
    parser = argparse.ArgumentParser(prog = 'Peak2genes.py', usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-a', '--Gene_annotation',action = 'store', type=str, dest = 'Gene_annotation',required=True,
                        help = 'Gene annotation file in bed format')
    parser.add_argument('-p', '--Peak_count',action = 'store', type=str, dest = 'Peak_count', required=True,
                        help = 'Peak counts file [FeatureCount result is used here]')
    parser.add_argument('-f', '--Ext_feature',action = 'store', type=str, dest = 'Ext_feature', required=True,
                        help = 'Specify extend features [TSS|TSS_genebody]')
    parser.add_argument('-e', '--Ext_size',action = 'store', type=int, dest = 'Ext_size', required=True,
                        help = 'Specify extend size (bp) from feature')
    parser.add_argument('-o', action = 'store', type=str, dest = 'gene_count_file',
                        help = 'Destination of reult gene count file')    

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)

    run(options)

if __name__ == '__main__':
    main()
