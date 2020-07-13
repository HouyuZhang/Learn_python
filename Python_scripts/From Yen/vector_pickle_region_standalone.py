#!/usr/bin/env python
# encoding: utf-8
"""
Created by Zhenhai Zhang on 2009-09-10.
Modified by Kuangyu Yen on 2011-11-05
Copyright (c) 2009 Bioinformatics and Genomics@PSU. All rights reserved.
Usage: python vector_pickle_region_standalone.py -suffix suffix -bi base_index -u uplim -d dn_lim -b bin -s smooth

###
This program will pick the defined region (uplim ~ dnlim) from vector files and generate cluster-enabled input file.
Title column starts from gene, then followed by "c" and indices of each column.
Note: base_index, normally, if you use my script, it will generate vectors from -5000 to +5000, [-5000, +5000) semiopen, regarding reference point.
In this case, base_index is 5000. 

bin: 1, do not bin; other than 1, you need to specify an aliquot of vector length
smooth: 0, do not smooth. Other than 0, you can choose any even number.
"""


import sys, os, csv, re
import random as rand
import operator
from itertools import ifilter
from socket import gethostname
from numpy import mean, array, zeros, ones
from math import sqrt, pi, exp
from time import time
from socket import gethostname
from datetime import datetime

sep="\t"


# define functions
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

def getFiles(my_folder = os.getcwd(), suffix = ""):
	"""
	getting all files under current directory;
	if suffix is not empty, only files with given suffix will be return
	"""
	tmp_result, result = [], []
	def checkSuffix(f):
		return f.endswith(suffix) and os.path.isfile(f)

	for dirpath, dirname, filename in os.walk(my_folder):		#os.walk generate 3-tuple (dirpath, dirnames, filenames)
		tmp_result = filename

	#print "tmp_files:", tmp_result

	if suffix != "":
		result = [x for x in ifilter(checkSuffix, tmp_result)]
		return result
	else:
		return tmp_result
	print "getFile1"
		
def getFiles(suffix):
	all_files = os.listdir(os.getcwd())
	return [x for x in all_files if x.endswith(suffix)]
	print "getFile2"

def getFiles(folder, suffix):
	all_files = os.listdir(os.getcwd() + "/" + folder)
	return [x for x in all_files if x.endswith(suffix)]
	print "getFile3"


def swapExt(s):
	if s == ".txt":
		return ".tab"
	if s == ".tab":
		return ".txt"


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

def validateList(l):
	# validate the list according to its type

	return not (None in set(l))





if __name__ == "__main__":
  #get arguments
  if len(sys.argv) < 13:
    print __doc__
    sys.exit(0)
    
  dict_args = processParas(sys.argv, suffix="suffix", bi="base_index", u="uplim", d="dnlim", b="bin", s="smooth")
  suffix, base_index, uplim, dnlim, bin, smooth = getParas(dict_args, "suffix", "base_index", "uplim", "dnlim", "bin", "smooth")
  
  up, dn = base_index + uplim, base_index + dnlim
  
  infiles = getFiles("", suffix)
  
  bin_list = range(uplim, dnlim, bin)
  title = ["gene"] + ["c" + str(x) for x in bin_list]
  
  for infile in infiles:
    print "processing file", infile
    
    fname, ext = os.path.splitext(infile)
    outfile = fname + "_u" + str(uplim) + "_d" + str(dnlim) + swapExt(suffix)
    reader, writer = csv.reader(open(infile, "rU"), delimiter = sep), csv.writer(open(outfile, "w"), delimiter = sep)
    writer.writerow(title)
    
    for no, row in enumerate(reader):
      if no > 0:
        gene, vector = row[0], map(float, row[1 :][up : dn])
                
        vector = smoothList(binList(vector, bin), smooth)
        writer.writerow([gene] + vector)
      
      else:
        pass
        
    print "finished processing ", infile
	


