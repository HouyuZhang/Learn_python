#!/usr/bin/env python
# encoding: utf-8
"""
Modified by Kuangyu Yen on 2012-10-26
Created by Zhenhai Zhang on 2009-07-29.
Copyright (c) 2009 Bioinformatics and Genomics@PSU. All rights reserved.
Usage: bin_smooth_single_file.py -i infile -b bin -s smooth
"""

#
# import modules start here
#

import sys, os, csv, re
import random as rand
from itertools import ifilter
from genome_info import *
from Bio import SeqIO
from Bio.Seq import Seq
from socket import gethostname
from numpy import mean, array, zeros, ones
from math import sqrt, pi, exp
from time import time

sep="\t"

#
# argument process functions from here
#
def processParas(para_list, **keywords):
	# process the parameter information, all parameter start with "-paraname"
	# return a dictionary (paraName, paraValue)
	# remove the first parameter which is the program name
	para_list = para_list[1 :]
	kwgs, values = para_list[ :: 2], para_list[1 :: 2]
	if len(kwgs) != len(values):
		print "number of keywords and values does not equal"
		sys.exit(0)
	
	kwgs = map(lambda x : keywords[x[1 :]], kwgs)
	values = map(evalValues, values)
	return dict(zip(kwgs,values))
	
def evalValues(v):
	# Evaluate strings and return a value corresponding to its real type (int, float, list, tuple)
	try:	return eval(v)
	except:	return v

def getParas(my_dict, *args):
		if len(args) == 1:	return my_dict[args[0]]
		else:	return (my_dict[arg] for arg in args)


def swapExt(s):
	if s == ".txt":
		return ".tab"
	if s == ".tab":
		return ".txt"


def get_title_content(f_csv):
	reader, master_list = csv.reader(open(f_csv, "rU"), delimiter = sep), []

	title = reader.next()
	for i in range(len(title)):
		master_list.append([])

	for row in reader:
		row = map(float, row)
		for ind, item in enumerate(row):
			master_list[ind].append(item)
	return title, master_list

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



if __name__ == '__main__':
	#get arguments
	if len(sys.argv) < 7:
		print __doc__
		sys.exit(0)
	dict_args = processParas(sys.argv, i="infile", b="bin", s="smooth")
	infile, bin, smooth = getParas(dict_args, "infile", "bin", "smooth")
	head, tail = os.path.splitext(infile)
	outfile = head + "_b" + str(bin) + "s" + str(smooth) + swapExt(tail)
	
	title, master_list = get_title_content(infile)
	
	min_coor, max_coor = int(min(master_list[0])), int(max(master_list[0]) + 1)
	coor_list = range(min_coor, max_coor, bin)
	
	master_list = [smoothList(binList(x, bin), smooth) for x in master_list]
	
	writer = csv.writer(open(outfile, "w"), delimiter = sep)
	writer.writerow(title)
	
	for i in range(len(coor_list)):
		aline = [ coor_list[i] ]
		for j in range(1, len(title)):
			aline.append(master_list[j][i])
		writer.writerow(aline)
		
	print "finished... result saved in ", outfile
	


