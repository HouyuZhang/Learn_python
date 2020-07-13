# -*- coding: utf-8 -*-
Description = '''
This script is for splitting PE fastq files according to a supplied barcode file.
Often, the barcode is in the first 8bp of Read_2.fastq, and the sequential 8bp is UMI, 
After Read_2 is demultiplexed, the corresponding reads in Read_1 is matched by header information.
Each barcode name corresponds to a fastq file.
If the sample information is provided, the result fastq will be renamed.
*NOTE: A perl script named "fastx_barcode_splitter.pl" from FASTX-toolkit is used.*
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2018-12-15.
Copyright (c) 2018 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import gzip, subprocess, shutil
import argparse,sys,os

class BarcodeSplitter(object):
    '''
    Split a multiplexed pair of FASTQ files into corresponding demultiplexed FASTQ files.
    '''

    def __init__(self, options): 
        self.args = ["fastx_barcode_splitter.pl", "--suffix", "_R2.fastq"]
        self.log = { "not matched" : 0, "demultiplexed": 0 }
        self.file_1, self.file_2 = options.fastq1, options.fastq2     
        
        self.bcfile = options.barcode_file
        self.args.extend(["--bcfile", options.barcode_file])
        
        self.mismatch = options.mismatch
        self.args.extend(["--mismatch", options.mismatch])
        
        self.prefix = options.OUT_prefix
        
        self.args.extend(["--prefix", options.OUT_prefix])
        
        self.position = options.position
        self.args.append("--" + options.position)
            
    def run(self):
        
        lib = os.path.basename(self.prefix)
        path = os.path.dirname(self.prefix)
        
        # demultiplex fastq2
        openf2 = gzip.open if self.file_2.endswith(".gz") else open
        with openf2(self.file_2, "rb") as f:
            p = subprocess.Popen(self.args, stdin=subprocess.PIPE)
            shutil.copyfileobj(f, p.stdin)
            p.stdin.close()
        p.wait()
        
        # get output file names
        fnames = [ fname for fname in os.listdir(path) if fname.startswith(lib) and fname.endswith("_R2.fastq") ]
    
        # read all outputs and the fastq1
        openf1 = gzip.open if self.file_1.endswith(".gz") else open
        with openf1(self.file_1, "rb") as f_1:
            total = 1
            # to finally close the open file arrays
            try:
            
                # open the mate #2 demultiplexed files for reading
                demultiplexed_files_2 = [ open(os.path.join(path, fname)) for fname in fnames ]
                
                # open the mate #1 demultiplexed files for writing
                demultiplexed_files_1 = [ open(os.path.join(path, fname.replace("_R2", "_R1")), "w") for fname in fnames ]
                                          
                # read 4 lines (1 read info) of each of the mate #2 
                # demultiplexed files
                slots = [ [ f_2.readline().strip() for j in range(4) ] for f_2 in demultiplexed_files_2 ]
                
                # read 4 lines of the file to be demultiplexed
                read_1_lines = [ f_1.readline() for j in range(4) ]
                
                # continue while the lines are not null
                while all(read_1_lines):
                    
                    # get the read name before any whitspace
                    read_name = read_1_lines[0].split()[0]
                    #print(read_name)
                    # search the read in the demultiplex files, if found, write read_1_lines to the corresponding mate #1 
                    # demultiplexed file then read the next 4 lines of the corresponding mate #2 demultiplexed file.
                    found = False
                    for i in range(len(fnames)):
                        read_2_lines = slots[i]
                        
                        # pass when it has reach the end of the file
                        if not read_2_lines[0]: continue
                        
                        if read_2_lines[0].split()[0] == read_name: 
                            demultiplexed_files_1[i].writelines(read_1_lines)
                            slots[i] = [ demultiplexed_files_2[i].readline() for j in range(4) ]
                            self.log["demultiplexed"] += 1
                            found = True
                            break
                    
                    # report unmatched mates                    
                    if not found:
                        sys.stderr.write("mate missing in Read_1: {}\n".format(read_name))
                        self.log["not matched"] += 1
                    
                    if total % 10 ** 5 == 0:
                        print("%d reads in Read_1 finished serching..." %total)
                        
                    # read the next 4 lines of file_2
                    read_1_lines = [ f_1.readline() for j in range(4) ]
                    total +=1 
            finally:
                map(lambda x: x.close(), demultiplexed_files_2)
                map(lambda x: x.close(), demultiplexed_files_1)
        
def barcode_splitter(options):
    '''
    Run a BarcodeSplitter instance.
    '''
    
    # handle options
    go = BarcodeSplitter(options)    
    
    # run the process
    go.run()

def Change_name(options):
    #sample names dictionary
    sample_dict = {}
    with open(options.sample) as s:
        for line in s:
            line = line.strip().split("_")
            sample_dict["_".join(line[2:4])] = "_".join(line[0:2])
    
    #recode barcode in dictionary
    barcode_dict = {}
    with open(options.barcode_file) as f:
        for line in f:
            line = line.strip().split()
            barcode_dict[line[1]] = line[0]
    
    #rename files
    lib = os.path.basename(options.OUT_prefix)
    path = os.path.dirname(options.OUT_prefix)
    
    for value in barcode_dict.values():
        sample_key = lib + value
        if sample_key in sample_dict.keys():
            os.system("mv " + options.OUT_prefix + value + "_R1.fastq " + os.path.join(path, sample_dict[sample_key]) + "_" + sample_key + "_R1.fastq")
            os.system("mv " + options.OUT_prefix + value + "_R2.fastq " + os.path.join(path, sample_dict[sample_key]) + "_" + sample_key + "_R2.fastq")
            
    
def main():   
    parser = argparse.ArgumentParser(prog = 'Fastq_splitter.py', usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--f1', action = 'store', type=str, dest = 'fastq1',required=True,
                        help = 'A fastq or fastq.gz file of Read_1')
    parser.add_argument('--f2', action = 'store', type=str, dest = 'fastq2',required=True,
                        help = 'A fastq or fastq.gz file of Read_2 (Barcode in this file)')
    parser.add_argument('-B', action = 'store', type=str, dest = 'barcode_file',required=True,
                        help = 'A barcode file which the first col is brcode name and second is barcode')
    parser.add_argument('-S', action = 'store', type=str, dest = 'sample',
                        help = 'If provided, it serves as reference for revising splitted fastq file names ')  
    parser.add_argument('--mismatch', action = 'store', type=str, dest = 'mismatch', default = '0',
                        help = 'Number of tolerated mismatches for barcodes.')
    parser.add_argument('--position', action = 'store', type=str, dest = 'position', default = 'bol', choices=['bol', 'eol'],
                        help = "Position for mapping barcode, Choose 'bol' from the 5' of reads, or 'eol' from the 3'. ") 
    parser.add_argument('--zip', action = 'store', type=str, dest = 'zip', default = 'NO', choices=['YES', 'NO'],
                        help = "Whether or not zip your result fastq files, default using 10 threads for gzipping ")    
    parser.add_argument('--prefix', action = 'store', type=str, dest = 'OUT_prefix', default = '.',
                        help = 'Prefix for the output splitted fastq files, default output into the current dir')

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)

    
    barcode_splitter(options)
    
    if options.sample:
        Change_name(options)
        
    if options.zip == "YES":
        os.system("parallel -j 10 'gzip {}' ::: " + os.path.dirname(options.OUT_prefix) + "/*fastq")
    
    print(options.fastq1 + " and " + options.fastq2 + " Have beed successfully splitted!")
    
if __name__ == '__main__':
    sys.exit(main())