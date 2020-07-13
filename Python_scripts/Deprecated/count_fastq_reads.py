# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 19:33:12 2017

@author: zhy
"""
import sys,getopt

def count_read(input_file,output_file):
	dic = {}
	total_reads = 0
	f = open(output_file,'w')
	for line in open(input_file,'r'):
		line = line.strip()
		if not (line.startswith('@') or line.startswith('+')):
			length = len(line)
			if length in dic:
				dic[length] += 1
			else:
				dic[length] = 0
				dic[length] += 1
		else:
			pass
	f.write('*'*100+'\n')
	f.write('Read_len\tRead_num\n')
	for key in sorted(dic):
		total_reads += int(dic[key]/2)
		f.write(str(key) +'\t\t'+str(int((dic[key]/2)))+'\n')
	
	f.write('*'*100+'\n')
	f.write(str(total_reads)+' reads have been input!\n')
	print('*'*100)
	print('| COUNT DONE!!!\n| The result has been written into '+output_file)
	print('*'*100)
	f.close()

def usage():
	help_info = '''
****************************************
#This script is for counting reads of fastq file
Usage:
-h Get help info
-i Path of input file
-o Path of output file
****************************************'''
	print(help_info)
	
def main():
	if len(sys.argv) < 2:
		usage() 
		sys.exit() 
	parameter_list = 'hi:o:'
	opts, args = getopt.getopt(sys.argv[1:], parameter_list)
	input_file = ''
	output_file = ''
	for op, value in opts:
		if op == '-i':
			input_file = value
		elif op == '-o':
			output_file = value
		elif op == '-h':
			usage()
			sys.exit()
	count_read(input_file,output_file)

if __name__ == '__main__':
	main()


