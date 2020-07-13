#!/usr/bin/env python
# encoding: utf-8

usage = """
Usage:
It will arrange your dataset into the same gene      as the reference file

***
<required>:
-r Reference file. The required gene list 
-t Test file. The records are retrieved if the IDs are in reference file.
-o Output file

[Options]:
for "-r" reference file:
  -d separator character used in reference file [Default = "\\t"]
  -g Index of column of featureid [Default = 0]

for “-t” reference file:
  -G Index of column of featureid [Default = 0]


example: python     ed_by_genelist.py -r TATA-less_Genelist.txt -t tag_mapping.txt -o tag_mapping_TATA-less.txt
***

Created by Kuangyu Yen on 2012-04-11.
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
"""

import sys, getopt, csv, os, re 

#
# define functions
#
def get_genelist(infile, sep, ind):
  reader, genelist = csv.reader(open(infile, "rU"), delimiter=sep), []
  for no, row in enumerate(reader):
    if no > 0:
      genelist.append(row[ind])
    else:
      pass
  
  print "finished retrieving interested GeneList!"
  return genelist
      
  

if __name__ == '__main__':
  if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)
  print sys.argv[1:]
  
  # get arguments
  sep_char_r, ind_featureid_r, ind_featureid_t = "\t", 0, 0
  
  optlist, alist = getopt.getopt(sys.argv[1:], "hr:d:g:t:G:o:")
  print optlist
  for opt in optlist:
    if opt[0] == "-h": sys.exit(usage)
    elif opt[0] == "-r": reference_file = opt[1]
    elif opt[0] == "-d": sep_char_r = opt[1]
    elif opt[0] == "-g": ind_featureid_r = int(opt[1])
    elif opt[0] == "-t": test_file = opt[1]
    elif opt[0] == "-G": ind_featureid_t = int(opt[1])
    elif opt[0] == "-o": out_file = opt[1]


  # retrieve data from base file
  reader, data = csv.reader(open(test_file, "rU"), delimiter="\t"), {}
  title = reader.next()

  for row in reader:
    data[row[ind_featureid_t]] = row[1:]

  print "finished retrieving base data "

  # write out interested genelist from the base data to output file
  genelist = get_genelist(reference_file, sep_char_r, ind_featureid_r)
  
  
  handle = csv.writer(open(out_file, "w"), delimiter="\t")
  handle.writerow(title)
  
  for line in genelist:
    try: handle.writerow([line]+data[line])
    except KeyError: pass
  
  print "Done!"

