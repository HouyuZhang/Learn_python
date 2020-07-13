# -*- coding: utf-8 -*-
#!/usr/bin/env python

#***************************parse_FastQC.py*****************************
#Script description: This scripts can return you a PASS_FAIL_WARN signed file for each sample through parse FastQC's result
#Version:1.2
#Created by Houyu Zhang on 2018-01-10.
#Copyright (c) 2018 __YenLab@SCUT__. All rights reserved.

#Alterations compared to version 1.0:
#1.Improve the flexiblity when read zip files
#2.Increase the parse speed using only one fast iteration
#Further alteration:
#1.Output excel file 
#***********************************************************************

import os,sys
import getopt
import glob

def usage():
	help_info = '''
Usage: python parse_FastQC.py [options] -i FastQC_results_path

#This scripts can return you a PASS_FAIL_WARN signed file for each sample through parse FastQC's result.

Options:
	-h Get help information
	-i Path of input file
	-o Path of output file
***

Created by Houyu Zhang on 2018-01-10.
Copyright (c) 2018 __YenLab@SCUT__. All rights reserved.
'''
	print(help_info)

def parse_file(need_to_parse,parsed_result,unzipped_dir_name):
	f = open(need_to_parse).readlines()
	out_line = unzipped_dir_name[:-7] + '\t'
	for line in f:
		line = line.split('\t')
		out_line = out_line + str(line[0]) + '\t'
	parsed_result.write(str(out_line) + '\n')


def get_targetfile(input_path,output_path):
	parsed_result = open(output_path,'w+')
	headline = 'Sample	Basic Statistics	Per base sequence quality	Per tile sequence quality	Per sequence quality scores\
	Per base sequence content	Per sequence GC content	Per base N content	Sequence Length Distribution	Sequence Duplication Levels\
	Overrepresented sequences	Adapter Content	Kmer Content\n'
	parsed_result.write(headline)
	os.chdir(input_path)
	for zip_file in glob.glob('*.zip'):
		os.system('unzip -q ' + zip_file)
		unzipped_dir_name = zip_file[:-4]
		os.chdir(unzipped_dir_name)
		need_to_parse = 'summary.txt'
		parse_file(need_to_parse,parsed_result,unzipped_dir_name)
		os.chdir('../')
		os.system('rm -rf '+ unzipped_dir_name) #remove unzipped dirctories

def main():
	if len(sys.argv) < 2:
		usage() 
		sys.exit() 
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hi:o:')
	except getopt.GetoptError as err:
		print('*'*10 + 'ERROR!!!'+'*'*10)
		print(err)
		print('*'*28)
		usage()
		sys.exit()
	input_path = ''
	output_path = './FastQC_parse_result.txt'
	for op, value in opts:
		if op == '-i':
			input_path = value
		elif op == '-o':
			output_path = value
		elif op == '-h':
			usage()
			sys.exit()
	get_targetfile(input_path,output_path)
	print('*'*100 + '\nParse Done!!!')
	print('*'*100)
	
if __name__ == '__main__':
	main()
