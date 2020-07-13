#!/usr/bin/env python
# encoding: utf-8

usage = """
Usage:
It maps tags with consideration of orientation of your cooridnate

***
<required>:
-g: the genetrack_index file 
-r: the coorinate file for mapping
-m: model organisms: cerevisiae, pombe, human, mouse, fly 
-u: upstream limitation 
-d: downstream limitation 
-o: output file

[Options]:
for "-g" genetrack file:
  -d separator character used in genetrack file [Default = "\\t"]
  -c Index of column of chromosome [Default = 0]
  -i Index of column of index [Default = 1]
  -F Index of column of Forward read [Default = 2]
  -R Index of column of Reverse read  [Default = 3]


for "-r" coordinate file:
  -T separator character used in reference file [Default = "\\t"]
  -C Index of column of chromosome [Default = 3] (It has to be chr1 format)
  -I Index of column of index [Default = -1]
  -S Index of column of strand [Default = 5]
  -E Index of column of featureid [Default = 4]


Example: python tag_mapping_standalone_allorganism_v2.py -g genetrack.tab -r coor.txt -m model -u up_lim -d dn_lim -o outfile
***

Created by Zhenhai Zhang on 2009-11-11. (ye_plusone_mapping.py)
Modified by Kuangyu Yen on 2012-04-10
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
"""


import sys, getopt, os, csv, re
from genome_info import *
from itertools import ifilter


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


def populateTagOnChrom(f, chrom, chr_len_organism):
    print "populating reads from %s on %s" %(f, chrom)
    # f is a genetrack input file contain lines in format: chrom    index    forward    reverse
    # chrom is chromosome information in format "chr01" or "chr16"
    result = [0] * (chr_len_organism[chrom] + 1)
    def checkChrom(aList):
        return aList[ind_chrom_g] == chrom

    reader = csv.reader(open(f, "rU"), delimiter = sep_char_g)
    #reader.next()        #skip title line
    total_tags = 0
    for aLine in ifilter(checkChrom, reader):
        index, total = int(aLine[ind_index_g]), float(aLine[ind_forward_g]) + float(aLine[ind_reverse_g])
        try:
            result[index] += total
            total_tags += total
        except:
            print index, chr_len_organism[chrom], len(result)
            print aLine

    #print "finished %d reads in total" %total_tags
    return result


def generate_any_coors(chrom):
    reader = csv.reader(open(coor_file, "rU"), delimiter = sep_char_r)
    #reader.next()        # skip title line
    
    def checkChrom(aList):
      return aList[ind_chrom_r] == chrom

    for row in ifilter(checkChrom, reader):
        gene, strand, coor = row[ind_featureid_r], row[ind_strand_r], int(row[ind_index_r])
        yield gene, strand, coor


def get_current_list(count_list, coor, strand, chrom_len):
    global up_lim, dn_lim

    start, end, left, right = 0, 0, [], []
    if strand in "WwFf+":
        start, end = coor + up_lim, coor + dn_lim
    else:
        start, end = coor - dn_lim + 1, coor - up_lim + 1

    if start < 0:
        start, left = 0, [0] * abs(start)

    if end > chrom_len:
        end, right = chrom_len, [0] * (end - chrom_len)

    result = left + count_list[start : end] + right

    if strand in "CcRr-":
        result.reverse()

    return result



if __name__ == "__main__":
    if len(sys.argv) < 2 or not sys.argv[1].startswith('-'): sys.exit(usage)
    print sys.argv[1:]
    
    # get arguments
    sep_char_g, ind_chrom_g, ind_index_g, ind_forward_g, ind_reverse_g = "\t", 0, 1, 2, 3 
    sep_char_r, ind_chrom_r, ind_index_r, ind_strand_r, ind_featureid_r = "\t", 3, -1, 5, 4
    
    optlist, alist = getopt.getopt(sys.argv[1:], 'hg:t:c:i:F:R:r:T:C:I:S:E:m:u:d:o:')
    for opt in optlist:
      if opt[0] == "-h": sys.exit(usage)
      elif opt[0] == "-g": index_file = opt[1]
      elif opt[0] == "-t": sep_char_g = opt[1]
      elif opt[0] == '-c': ind_chrom_g = int(opt[1])
      elif opt[0] == '-i': ind_index_g = int(opt[1])
      elif opt[0] == '-F': ind_forward_g = int(opt[1])
      elif opt[0] == '-R': ind_reverse_g = int(opt[1])
      
      elif opt[0] == '-r': coor_file = opt[1]
      elif opt[0] == '-T': sep_char_r = opt[1]
      elif opt[0] == '-C': ind_chrom_r = int(opt[1])
      elif opt[0] == '-I': ind_index_r = int(opt[1])
      elif opt[0] == '-S': ind_strand_r = int(opt[1])
      elif opt[0] == '-E': ind_featureid_r = int(opt[1])
      
      elif opt[0] == '-m': organism = opt[1]
      elif opt[0] == '-u': up_lim = int(opt[1])
      elif opt[0] == '-d': dn_lim = int(opt[1])
      elif opt[0] == "-o": out_file = opt[1]
        
    chr_len_organism = which_organism(organism)

    title = range(up_lim, dn_lim)
    title.insert(0, "gene")
    

    writer = csv.writer(open(out_file, "w"), delimiter = "\t")
    writer.writerow(title)
    for chrom, chrom_len in chr_len_organism.items():
        print "mapping reads on chromosome: %s to features" %chrom,
        chrom_count_list = populateTagOnChrom(index_file, chrom, chr_len_organism)
        
        for gene, strand, coor in generate_any_coors(chrom):
            curr_list = get_current_list(chrom_count_list, coor, strand, chrom_len)
            writer.writerow([gene] + curr_list)
        
        print "done!"
    print "done processing file %s" %index_file

    
    
