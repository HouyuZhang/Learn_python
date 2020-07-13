#!/usr/bin/env python
# encoding: utf-8
"""
Created by Zhenhai Zhang on 2009-09-05.
Modified by Kuangyu Yen  on 2011-11-04, 2013-01-17
Remodified by Houyu Zhang on 2018-12-13 for deleting the 'Margin Center features'
Copyright (c) 2009 Bioinformatics and Genomics@PSU. All rights reserved.

Usage: compose_multi_vector_files.py -s suffix -n name -c normalization -bi base_index 

You can choose if you want to normalize your tag counts or not by changing the parameter -c (yes, then it will be normalized)
"""

import sys, os, csv, re
import random as rand
import operator
from itertools import ifilter
from Bio import SeqIO
from Bio.Seq import Seq
from socket import gethostname
from numpy import mean, array, zeros, ones, shape
from math import sqrt, pi, exp
from time import time
from socket import gethostname
from datetime import datetime

sep="\t"

#
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
    try:    return eval(v)
    except:    return v

def getParas(my_dict, *args):
        if len(args) == 1:    return my_dict[args[0]]
        else:    return (my_dict[arg] for arg in args)


def getFiles(folder, suffix):
    all_files = os.listdir(os.getcwd() + "/" + folder)
    return [x for x in all_files if x.endswith(suffix)]


def swapExt(s):
    if s == ".txt":
        return ".tab"
    if s == ".tab":
        return ".txt"


def group_infiles(l, name):
    global dict_infiles
    l.sort()
    
    dict_infiles[name] = [x for x in l if x.find(name) > 0]


def get_single_exp_vector(f):
    print "calculate count vecotrs in %s for all genes" %f
    global base_index
    reader, count_array, total = csv.reader(open(f, "rU"), delimiter = sep), zeros(2 * base_index), 0
    reader.next()
    for row in reader:
        gene, counts = row[0], array(map(float, row[1 :]))
        if shape(count_array) != shape(counts):
            pass    
        else:
            count_array += counts
            total += 1
        if total % 10 ** 3 == 0:
            print "%d finished..." %total
    print "finished... %d in total..." %total
    return count_array, sum(count_array)

def aggregate_count_vectors(grouped_infiles, normalization):
    global dict_infiles, base_index
    item_list, sum_list, result = [], [], dict()

    for ind, item in enumerate(grouped_infiles):
        result[item], curr_sum = get_single_exp_vector(item)
        sum_list.append(curr_sum)
    max_sum = max(sum_list)
    
    # normalize to maximum sum
    
    if normalization == "yes":
        for item in result:
            curr_sum = sum(result[item])
            factor = max_sum / curr_sum
            result[item] = result[item] * factor
    else: pass
    
    return result

def write_to_file(outfile, dict_list):
    global base_index
    
    exps = list(dict_list.keys())
    exps.sort()
    
    writer = csv.writer(open(outfile, "w"), delimiter = sep)
    writer.writerow(["coor"] + exps)
    coor_list = range(- base_index, base_index)
    
    for i in range(2 * base_index):
        aline = [coor_list[i]]
        for item in exps:
            aline.append(dict_list[item][i])
        writer.writerow(aline)
    

if __name__ == '__main__':
    #get arguments
    if len(sys.argv) < 7:
        print __doc__
        sys.exit(0)
    
    print sys.argv[1:]
    
    dict_args = processParas(sys.argv, s="suffix", n="name", c="normalization", bi="base_index")
    suffix, name, normalization, base_index = getParas(dict_args, "suffix", "name", "normalization", "base_index")
    print suffix, name, normalization, base_index
    
    infiles, dict_infiles = getFiles("", suffix), dict()
    
    group_infiles(infiles, name)
    
    for infile_group in dict_infiles:
        print "processing group : %s " %infile_group
        
        outfile = infile_group + "_" + str(base_index) + swapExt(suffix)
        dict_count_list = aggregate_count_vectors(dict_infiles[infile_group], normalization)
        
        print "writing result of %s to file %s" %(infile_group, outfile)
        write_to_file(outfile, dict_count_list)
        
        del dict_count_list
        
            
    


