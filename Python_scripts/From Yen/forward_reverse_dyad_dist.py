#!/usr/bin/env python
# encoding: utf-8
usage = """
Usage:
It calculates the median distance between c-w peak pairs, and reports the half distance into shift.config
The peak filename needs to be looked like "genetrack_index.tab.s520.txt"

***
<required>:
-s: suffix
-m: model organisms: cerevisiae, pombe, human, mouse, fly 

Example: python forward_reverse_dyad_dist.py -s .txt -m cerevisiae
***

Created by Zhenhai Zhang on 2009-08-05.
Modified by Kuangyu on 2012-06-07, 2013-03-07
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
"""

import sys, getopt, os, csv, re
from genome_info import *
from itertools import ifilter
from socket import gethostname
from numpy import array, zeros
#from pylab import *
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

def getFiles(suffix):
	all_files = os.listdir(os.getcwd())
	return [x for x in all_files if x.endswith(suffix)]

def swapExt(s):
	if s == ".txt": return ".tab"
	if s == ".tab": return ".txt"

def populate_peaks(infile, chrom, chrom_len):
	print "retrieving dyad information on %s from %s" %(chrom, infile)
	print chrom_len
	def check_chrom(l):
		return l[0] == chrom and float(l[4]) > 1.0 # this removes all the singleton

	reader, f_list, r_vec, total_f, total_r = csv.reader(open(infile, "rU"), delimiter = "\t"), [], zeros(chrom_len + 1), 0, 0
	for peak in ifilter(check_chrom, reader):
		strand, dyad = peak[6], ( int(peak[3]) + int(peak[4]) ) / 2
		if strand in "Ww+":
			f_list.append(dyad)
			total_f += 1
		elif strand in "Cc-":
			try: 
				r_vec[dyad] += 1
				total_r += 1
			except IndexError:
				pass
	print "finished..... %d peaks in forward strand and %d peaks in reverse strand" %(total_f, total_r)
	return f_list, r_vec

def pop_dist(my_list, my_vec):
	print "calculating distance for adjacent dyads..."
	result, chrom_len = zeros(2 * region_len), len(my_vec) - 1

	for coor in my_list:
		if coor < region_len or coor + region_len > chrom_len: pass
		else: result += my_vec[coor - region_len : coor + region_len]
	return result

def validateList(l):
	# validate the list according to its type
	return not (None in set(l))
	
def binList(l, bin_size):
	result = []
	if len(l) % bin_size <> 0:
		raise Exception("Please select bin size appropriately")
	if not validateList(l):
		raise Exception("Invalid value in the list: None")
	for i in range(0, len(l) / bin_size):
		result.append(sum(l[i * bin_size : (i + 1) * bin_size]))
	return result

def smoothList(l, smooth_factor):
	#smooth the list value;
	#calculate the sum of surrounding smooth_factor values together with ith value
	result = []
	l_len = len(l)
	if not validateList(l):
		raise Exception("invalid value in the list: None")
	for i in range(0, l_len):
		l_lim = max(0, i - smooth_factor/2)
		r_lim = min(l_len, i + smooth_factor/2 + 1)
		result.append( float(sum(l[l_lim : r_lim])) / float(r_lim - l_lim))
	return result



def find_shift(my_vec):
	my_list = list(my_vec)
	most_present = my_list[region_len :].index(max(my_list[region_len :]))
	shift = most_present / 2
	if most_present % 2 > 0:
		shift += 1
	return shift

def parse_dash(l):
	result = []
	for ind, item in enumerate(l):
		if l[ind:ind+5] == ".tab.":
			result.append(ind)
	return result[0]

def process_name(org_name): # this function will trace back to original input file that led to the genetrack peak file
	start, end = 0, parse_dash(org_name)+4
	#print org_name, start, end
	return org_name[start : end]



if __name__ == '__main__':
	if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)
	
	# get argument
	optlist, alist = getopt.getopt(sys.argv[1:], 'hs:m:')
	for opt in optlist:
		if opt[0] == "-h": sys.exit(usage)
		elif opt[0] == "-s": suffix = opt[1]
		elif opt[0] == "-m": organism = opt[1]
	
	chr_len_organism = which_organism(organism)
	
	infiles = getFiles(suffix)
	region_len, bin, smooth, master_list = 200, 5, 2, [] 		# 200 bp each side
	
	config_file, plot_file = "../a genetrack_index/shift.config", "shift" + swapExt(".txt")
	config_handle = open(config_file, "w")
	
	titles = [process_name(x) for x in infiles]
	print titles
	

	for file_index, infile in enumerate(infiles):
		dist_vec = zeros(2 * region_len)
		
		for chrom in chr_len_organism:
			f_list, r_vec = populate_peaks(infile, chrom, chr_len_organism[chrom])
			dist_vec += pop_dist(f_list, r_vec)
			
		print "finished processing file", infile
		smooth_vec = smoothList(binList(dist_vec, bin), smooth)
		master_list.append(smooth_vec)
		curr_shift = find_shift(dist_vec)
		print "suggested shift for %s is %d" %(infile, curr_shift)
		config_handle.write(titles[file_index] + "\t" + str(curr_shift) + os.linesep)
		
	
		#plot(range(-region_len, region_len, bin), smooth_vec)
	#show()
	
	config_handle.close()

