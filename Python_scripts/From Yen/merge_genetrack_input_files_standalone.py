#!/usr/bin/env python
# encoding: utf-8
usage = """
Usage: 
It merges all the genetrack files into one genetrack files. 

***
<required>:
-m: model organisms: yeast, human, mouse, fly
-g: all the genetrack files you wish to merge into one
-o: output file

Example: python merge_genetrack_input_files_standalone.py -m model -g infile1 infile2 infile3 infile4 .... -o outfile
***

Created by Zhenhai Zhang on 2010-02-02.
Modified by Kuangyu Yen on 2012-04-20
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
"""

import sys, getopt, os, csv, re
from genome_info import *
from itertools import ifilter
from socket import gethostname

#
# define functions
#
def which_organism(organism):
  if organism == "cerevisiae": chr_len_organism = chr_len_yeast
  elif organism == "human": chr_len_organism = chr_len_hg18
  elif organism == "mouse": chr_len_organism = chr_len_mm9
  elif organism == "fly": chr_len_organism = chr_len_dm3
  
  return chr_len_organism


def get_chrom_tags(chrom, chrom_len):
    global infiles
    f_list, r_list = [0] * (chrom_len + 1), [0] * (chrom_len + 1)
    
    def check(l):
        return l[0] == chrom
    
    print "populating tags on chromosome: %s" %chrom,
    for infile in infiles:
        print "%s... " %infile, 
        reader = csv.reader(open(infile, "rU"), delimiter = "\t")
        for row in ifilter(check, reader):
            index, forward, reverse = int(row[1]), int(row[2]), int(row[3])
            try:
                f_list[index] += forward
                r_list[index] += reverse
            except IndexError:
                print "index out of range; index: %d chrom_len: %d" %(index, chrom_len)
            
    print "done!"
    
    return f_list, r_list
    
def write_to_file(chrom, chrom_len, f_list, r_list):
    global writer
    print "writing result on chrom: %s to file..." %chrom,
    for index in range(chrom_len + 1):
        forward, reverse = f_list[index], r_list[index]
        if forward + reverse > 0:
            writer.writerow([chrom, index, forward, reverse])
            
print "done!"
    
    
if __name__ == "__main__":
  if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)
  
  organism = sys.argv[2]
  infiles, outfile = sys.argv[ 4 : -2], sys.argv[-1]
  print infiles, outfile
  
  
  chr_len_organism = which_organism(organism)
  
  writer = csv.writer(open(outfile, "w"), delimiter = "\t")
  writer.writerow(["chrom", "index", "forward", "reverse"])
  
  for chrom, chrom_len in chr_len_organism.items():
    f_list, r_list = get_chrom_tags(chrom, chrom_len)
    write_to_file(chrom, chrom_len, f_list, r_list)
        
        
