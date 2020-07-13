# -*- coding: utf-8 -*-

import getopt,sys,os
import numpy as np

if len(sys.argv) < 2:
    print("This script will calculate the corresponding chromatin state (-r) for your offered regions (-p) under specified resolution (-b)")
    print("NOTE: multi-sample tmp files shouldn't specified same name (-t)!!!")
    print("USAGE: python group_composit.py -p 12_groups_3.txt -r MA9_12group_geno.bed -b 1 -t tmp.txt -w 2000 -s 12 -o test.txt")
    sys.exit() 

opts, args = getopt.getopt(sys.argv[1:], 'o:p:r:b:t:w:s:')
for op, value in opts:
    if op == '-o': out = value
    if op == '-p': Peak = value
    if op == '-r': Ref = value
    if op == '-b': bin_size = int(value)
    if op == '-t': tmp_sp = value
    if op == '-w': width = int(value)
    if op == '-s': state = int(value)

line = 0
with open(Peak) as peak:
    for peak_line in peak:
        line += 1

matrix = np.zeros((line,int(width/bin_size)))

with open(Peak) as peak:
    peak_rows = 0
    for peak_line in peak:
        peak_rows += 1
        if peak_rows % 10 == 0:
            print("Processing the " + str(peak_rows) + "th peak...")

        tmpo = open(tmp_sp,"w+")
        tmp_peak = peak_line.strip().split()
        for i in range(int((int(tmp_peak[2])-int(tmp_peak[1]))/bin_size)):
            tmpo.write(tmp_peak[0] + '\t' + str(int(tmp_peak[1]) + i*bin_size) + '\t' +  str(int(tmp_peak[1]) + i*bin_size + 1) + '\n')
        tmpo.close()

        os.system("bedtools intersect -a " + Ref + " -b " + tmp_sp + " -wa -wb > " + tmp_sp + "_")
        os.system("rm -f " + tmp_sp)

        with open(tmp_sp + "_") as tmp_file:
            row = 0
            for line in tmp_file:
                tmp_line = line.strip().split()
                group = int(tmp_line[3].lstrip('E'))
                matrix[peak_rows,row] += group
                row +=1
        tmp_file.close()

        os.system("rm -f " + tmp_sp + "_")
    #matrix = matrix/peak_rows
    np.savetxt(out, matrix, delimiter="\t", header='\t'.join(['E' + str(x) for x in list(range(1,state + 1))]))
    print("Work Done!")