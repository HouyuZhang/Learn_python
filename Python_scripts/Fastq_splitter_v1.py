# -*- coding: utf-8 -*-
Description = '''
This script is for splitting PE fastq files according to a supplied barcode file,
each barcode name corresponds to a fastq file. The barcode is in the first 8bp of Read_2.fastq,
and the sequential 8bp is UMI, the corresponding reads in Read_1 is matched by header information.
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2018-12-15.
Copyright (c) 2018 __YenLab@SCUT__. All rights reserved.
'''.lstrip()


import gzip 
import argparse,sys,os
from itertools import islice

def Splitter(options):
    
    #recode barcode in dictionary
    barcode_dict = {}
    with open(options.barcode_file) as f:
        for line in f:
            line = line.strip().split()
            barcode_dict[line[1]] = line[0]
    
    #build every single file for barcodes        
    for value in barcode_dict.values():
        vars()[value + "_R1"] = gzip.open(options.OUT_prefix + value + "_R1.fastq.gz","wb")
        vars()[value + "_R2"] = gzip.open(options.OUT_prefix + value + "_R2.fastq.gz","wb")     
        
    
    #handle Read_2 and Read_1 file
    f1 = gzip.open(options.fastq1,"rb")
    
    with gzip.open(options.fastq2,"rb") as f2:
        while True:
            unit_R2 = list(islice(f2, 4))
            unit_R1 = list(islice(f1, 4))
            if len(unit_R1) > 0:
                sc = str(unit_R2[1][:8])[2:10]
                if sc in barcode_dict.keys():
                    sc_R1 = barcode_dict[sc] + "_R1"
                    sc_R2 = barcode_dict[sc] + "_R2"
                    vars()[sc_R1].writelines(unit_R1)
                    vars()[sc_R2].writelines(unit_R2)
            else:
                break
                
    f1.close()    
    for value in barcode_dict.values():
        vars()[value + "_R1"].close()
        vars()[value + "_R2"].close()     

    if options.sample:
            #sample names dictionary
            sample_dict = {}
            with open(options.sample) as s:
                for line in s:
                    line = line.strip().split("_")
                    sample_dict["_".join(line[2:4])] = "_".join(line[0:2])
            
            #rename files
            lib = os.path.basename(options.OUT_prefix)
            path = os.path.dirname(options.OUT_prefix) + "/"
            for value in barcode_dict.values():
                sample_key = lib + value
                if sample_key in sample_dict.keys():
                    os.system("mv " + options.OUT_prefix + value + "_R1.fastq.gz " + path + sample_dict[sample_key] + "_" + sample_key + "_R1.fastq.gz")
                    os.system("mv " + options.OUT_prefix + value + "_R2.fastq.gz " + path + sample_dict[sample_key] + "_" + sample_key + "_R2.fastq.gz")

def main():   
    parser = argparse.ArgumentParser(prog = 'Fastq_splitter.py', usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--f1', action = 'store', type=str, dest = 'fastq1',required=True,
                        help = 'A fastq or fastq.gz file of Read_1')
    parser.add_argument('--f2', action = 'store', type=str, dest = 'fastq2',required=True,
                        help = 'A fastq or fastq.gz file of Read_2')
    parser.add_argument('-B', action = 'store', type=str, dest = 'barcode_file',required=True,
                        help = 'A barcode file which the first col is brcode name and second is barcode')
    parser.add_argument('-S', action = 'store', type=str, dest = 'sample',
                        help = 'If provided, it serves as reference for revising splitted fastq file names ')  
    parser.add_argument('--prefix', action = 'store', type=str, dest = 'OUT_prefix', default = '.',
                        help = 'Prefix for the output splitted fastq files')

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)
                   
    Splitter(options)

if __name__ == '__main__':
    main()