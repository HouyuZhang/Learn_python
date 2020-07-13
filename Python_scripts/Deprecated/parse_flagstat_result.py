import re
import os,sys
import getopt,glob

def usage():
	help_info = '''
*********************************************************************************************************
#This Scripts is for parsing mapping rate for your sam files
Usage:
-h Get help info
-i Path of input file
-o Path of output file
*********************************************************************************************************'''
	print(help_info)

def parse_flagstat_file(sample_name,parsed_result):
	f = open(sample_name).readlines()
	parsed_result.write(sample_name + '\t')
	for lines in f:
		if re.search('[0-9] mapped',lines):
			totally_mapped = re.search('[0-9]{1,3}\.[0-9]{1,3}%',lines).group()
			parsed_result.write(str(totally_mapped) +'\t')
		elif re.search('properly paired',lines):
			if re.search('[0-9]{1,3}\.[0-9]{1,3}%',lines):
				properly_paired = re.search('[0-9]{1,3}\.[0-9]{1,3}%',lines).group()
				parsed_result.write(str(properly_paired) +'\t\n')
			else:
				properly_paired = 'NA'
				parsed_result.write(str(properly_paired) +'\t\n')
			


def get_bam_file(input_path,output_path):
	parsed_result = open(output_path,'w+')
	headline = 'Samples	Totally mapped	properly paired\n'
	parsed_result.write(headline)
	os.chdir(input_path)
	for sam_file in glob.glob('*.sam'):
		sample_name = sam_file[:-4]
		os.system('samtools flagstat ' + sam_file + ' > ' + sample_name)
		parse_flagstat_file(sample_name,parsed_result)
		os.system('rm -rf ' + sample_name)



def main():
	if len(sys.argv) < 2:
		usage() 
		sys.exit() 
	parameter_list = 'hi:o:'
	opts, args = getopt.getopt(sys.argv[1:], parameter_list)
	input_path = ''
	output_path = './flagstat_parse_result.txt'
	for op, value in opts:
		if op == '-i':
			input_path = value
		elif op == '-o':
			output_path = value
		elif op == '-h':
			usage()
			sys.exit()
	get_bam_file(input_path,output_path)
	print('*'*100 + '\nParse sam_files Done!!!')
	print('*'*100)
if __name__ == '__main__':
	main()

































