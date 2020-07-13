#!/usr/bin/env python
# encoding: utf-8
"""
ye_plusone_mapping.py -s suffix -u up_lim -d dn_lim

This program map tags into Tsukiyama et. al. defined plus one dyads.

Usage:
Put all experimental data files (genetrack shifted tag index file in format desceibed below) into a folder; 
Make sure they have same extentioin (.txt or .tab).

input file format:
#chrom	index	forward	reverse

***
suffix:		.txt or .tab
up_lim:		-300 if you wanna map to 300bp upstream; 100 if you wanna map from downstream 100bp
dn_lim: 	downstream limitation.

example:	ye_plusone_mapping.py -s .tab -u -300 -d 500
This command will map your tags to upstream 300 bp to downstream 500 bp from plus one dyad. (strand information of features has been taken care)

***

Created by Zhenhai Zhang on 2009-11-11.
Copyright (c) 2009 The Pennsylvania State Univ.. All rights reserved.
"""

import sys
from socket import gethostname

# append code path to python path
hostname = gethostname()
if hostname == "apollo":
	sys.path.append("/export/share/zhenhai")
elif hostname == "zyz5012.bx.psu.edu":
	pass
elif hostname == "zhenhai-macpro.bx.psu.edu":
	sys.path.append("/Users/zhenhaizhang/Dropbox/codes")
else:
	sys.path.append("/home/zzhang/codes")

from mytools import *
from datetime import datetime

def get_current_list(count_list, coor, strand, chrom_len):
	global up_lim, dn_lim
	
	start, end, left, right = 0, 0, [], []
	if strand in "WwFf+":
		start, end = coor + up_lim, coor + dn_lim
	else:
		start, end = coor - dn_lim + 1, coor - up_lim + 1
	
	if start < 0:
		start, left = 0, [0.0] * abs(start)
	
	if end > chrom_len:
		end, right = chrom_len + 1, [0.0] * (end - chrom_len)
		
	result = left + count_list[start : end] + right
	
	if strand in "CcRr-":
		result.reverse()
	
	return result
	
	

if __name__ == "__main__":
	
	# get arguments
	if len(sys.argv) <  7:
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, s="suffix", u="up_lim", d="dn_lim")
	suffix, up_lim, dn_lim = getParas(dict_args, "suffix", "up_lim", "dn_lim")
	
	infiles = getFiles("", suffix)
	for infile in infiles:
		outfile = infile[ : - len(suffix)] + "_plus_one" + swapExt(suffix)
		writer = csv.writer(open(outfile, "w"), delimiter = sep)
		for chrom, chrom_len in chr_len_yeast.items():
			print "mapping reads on chromosome: %s to features" %chrom,
			chrom_count_list = populateTagOnChrom(infile, chrom)
			
			for gene, strand, coor in generate_plus_one_coor(chrom):
				curr_list = get_current_list(chrom_count_list, coor, strand, chrom_len)
				writer.writerow([gene] + curr_list)
			
			print "done!"
		print "done processing file %s" %infile

			
	
