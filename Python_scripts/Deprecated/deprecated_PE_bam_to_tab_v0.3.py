#!/usr/bin/env python
# -*- coding: utf-8 -*-

Description = '''
This scripts can help you calculate the tag coverage acorss a genome from a bam file (for PE sequencing currently).
**ATTENTION**: 
1. The input bam file should be only proper paired reads left and without filter low qualtiyreads.
   (samtools view -b -f 2 -F 2304 file.sam > file.bam)
2. Please don't sort the bam file using command like (samtools sort file.bam)
3. I will discard all scaffolds or contigs (i.e. chr1 ...chrY will be used).
'''.lstrip()

Copyright = '''
***
Inspired by Pinding's script -- 'bam_to_tab.py'
Created by Houyu Zhang on 2018-05-18.
Copyright (c) 2018 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import numpy as np
import argparse,sys,os,re
import datetime

def write_tab(genome_array,options):
    out = open(options.tab_file,'w+')
    #out.write('chr' + '\t' + 'position' + '\t' + 'number' + '\n')
    [rows, cols] = genome_array.shape  
    for i in range(rows):  
        for j in range(cols): 
            if genome_array[i,j] != 0:
                #the + = strand tags are assigned to 0 
                if i == 19:
                    out.write('chrX' + '\t' + str(j+1) + '\t' + '0' + '\t' + '0' + '\t' + str(genome_array[i, j]) + '\n')
                elif i == 20:
                    out.write('chrY' + '\t' + str(j+1) + '\t' + '0' + '\t' + '0' + '\t' + str(genome_array[i, j]) + '\n')
                elif i == 21:
                    out.write('chrM' + '\t' + str(j+1) + '\t' + '0' + '\t' + '0' + '\t' + str(genome_array[i, j]) + '\n')
                else:
                    out.write('chr'+str(i+1) + '\t' + str(j+1) + '\t' + '0' + '\t' + '0' + '\t' + str(genome_array[i, j]) + '\n')
    
    out.close()
    
def check_pair(tmp_line1, tmp_line2):
        return(tmp_line1[3][:-2] == tmp_line2[3][:-2])
 

def map_bam2genome(options):
    
    if not os.path.exists(options.bam_file):
	print(options.bam_file + " doesn't exist, please check your file!")
	sys.exit(1)
    bed_file = options.bam_file + '.bed'
    print("Converting bam to bed ...")
    start_time = datetime.datetime.now()
    os.system("bamToBed -i " + options.bam_file + " > " + bed_file )
    end_time = datetime.datetime.now()
    print("Done after " + str((end_time - start_time).seconds) + " seconds!")
    
    bed = open(bed_file)
    #bed = open("_f2_F2048.bed")

    genome_array = np.zeros((22,200000000), dtype=int)

    start_time = datetime.datetime.now()
    print("Starting mapping tags to genome...")
    
    while True:
        pair_range = []
        tmp_line1 = bed.readline().split('\t')
        tmp_line2 = bed.readline().split('\t')
       # print(tmp_line1 ,tmp_line2)

        if tmp_line1 == [''] or tmp_line2 == ['']:
            break
        if tmp_line1[4] != '0' and tmp_line2[4] != '0':
            if not check_pair(tmp_line1,tmp_line2):
                print("Please check your bed file and make sure the pair reads are adjcent! [maybe you have sorted your bam file]")
                sys.exit(1)
            if re.match('chr[0-9]{1,2}|chr[X|Y|M]',tmp_line1[0]):
                if tmp_line1[5] == '+':
                    pair_range.append(int(tmp_line1[1]) + 1)
                    pair_range.append(int(tmp_line2[2]))
                else:
                    pair_range.append(int(tmp_line2[1]) + 1)
                    pair_range.append(int(tmp_line1[2]))
            
                chr_index = re.search('[0-9]{1,2}|[X|Y|M]',tmp_line1[0]).group()
                if chr_index == 'X':
                    genome_array[19,pair_range[0]-1:pair_range[1]] += 1
                elif chr_index == 'Y':
                    genome_array[20,pair_range[0]-1:pair_range[1]] += 1
                elif chr_index == 'M':
                    genome_array[21,pair_range[0]-1:pair_range[1]] += 1
                else:
                    genome_array[int(chr_index)-1,pair_range[0]-1:pair_range[1]] += 1
    bed.close()        
    end_time = datetime.datetime.now()
    print("Done after " + str((end_time - start_time).seconds) + " seconds!")
    
    start_time = datetime.datetime.now()
    print("Writing result to " + options.tab_file + '...')
    write_tab(genome_array,options)
    end_time = datetime.datetime.now()
    print("Done after " + str((end_time - start_time).seconds) + " seconds!")
    print('I have converted your ' + options.bam_file + ' to ' + options.tab_file )

def main():   
    parser = argparse.ArgumentParser(prog = 'PE_bam_to_tab.py', usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-i', action = 'store', type=str, dest = 'bam_file',required=True,
                        help = 'A bam file which filtered by proper paired [i.e. only FLAG == 0x2 left].')
    parser.add_argument('-o', action = 'store', type=str, dest = 'tab_file', default = 'result.tab',
                        help = 'Destination of result tab file.')

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)
                   
    map_bam2genome(options)

if __name__ == '__main__':
    main()
