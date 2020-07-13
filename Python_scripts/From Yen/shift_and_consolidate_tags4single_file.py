#!/usr/bin/env python
# encoding: utf-8

usage = """
Usage:  
Shift the tags. 

***
<required>:
-f: genetrack file
-s: how much to shift
-m: model organisms: cerevisiae, pombe, human, mouse, fly

Example: shift_and_consolidate_tags4single_file.py -g genetrack.tab -s 0 -m yeast
***

Created by Zhenhai Zhang on 2009-08-11.
Modified by Kuangyu Yen on 2012-04-16
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
"""
import sys, getopt, os, csv, re
from genome_info import *
from socket import gethostname
from itertools import ifilter
from datetime import datetime

#
# define functions
#
def which_organism(organism):
  if organism == "cerevisiae": chr_len_organism = chr_len_yeast
  elif organism == "pombe": chr_len_organism = chr_len_sp11
  elif organism == "human": chr_len_organism = chr_len_hg18
  elif organism == "mouse": chr_len_organism = chr_len_mm9
  elif organism == "fly": chr_len_organism = chr_len_dm3
  
  return chr_len_organism
  

def get_tag_dict(infile, chrom, chrom_len, shift):
	
	print "retrieving tag informatoin on %s" %(chrom)
	result = dict()		# dict[coor] = [forward, reverse]
	
	def check_chrom(l):
		return l[0] == chrom
	
	reader = csv.reader(open(infile, "rU"), delimiter = "\t")
	for row in ifilter(check_chrom, reader):
		coor, forward, reverse = int(row[1]), int(row[2]), int(row[3])
		f_coor, r_coor = coor + shift, coor - shift
		
		if f_coor <= chrom_len:
			if f_coor not in result:
				result[f_coor] = [0, 0]

			result[f_coor][0] += forward

		else:
			pass; 		# index exceeds chromosome length

		if r_coor > 0:
			if r_coor not in result:
				result[r_coor] = [0, 0]

			result[r_coor][1] += reverse

		else:
			pass;		# index reaches left end of chromosome
			
	return result
		
if __name__ == '__main__':
  if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)
  
  optlist, alist = getopt.getopt(sys.argv[1:], "hg:s:m:")
  for opt in optlist:
    if opt[0] == "-h": sys.exit(usage)
    elif opt[0] == "-g": genetrack_file = opt[1]
    elif opt[0] == "-s": shift = int(opt[1])
    elif opt[0] == "-m": organism = opt[1]
  
  chr_len_organism = which_organism(organism)
  
  
  print genetrack_file, ":", shift
  
  outfile = genetrack_file[:-4] + "_shift" + str(shift) + ".tab"
  
  print "processing file %s, shift %d, result will be saved in %s...." %(genetrack_file, shift, outfile)
  writer = csv.writer(open(outfile, "w"), delimiter = "\t")
  writer.writerow(["chrom", "index", "forward", "reverse"])
  
  for chrom in chr_len_organism:
    chrom_len = chr_len_organism
    tag_index_dict = get_tag_dict(genetrack_file, chrom, chrom_len, shift)
    
    keys = tag_index_dict.keys()
    keys.sort()
    
    for key in keys:
      writer.writerow([chrom, key] + tag_index_dict[key])
      
  print "processing %s finished...." %(genetrack_file)		
	

