# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 11:36:38 2018

@author: Administrator
"""
import argparse, sys, textwrap
import numpy as np


__DESCRIPTION__ = "This script is for comparing the segment results between IDEAS and ChromHMM"
__COPYRIGHT__ = textwrap.dedent('''***
Created by Houyu Zhang on 2018-11-1.
Copyright (c) 2018 __YenLab@SCUT__. All rights reserved.''')

       

def Compare(args):
    
    IDEAS_dict = {}
    with open(args.IDEAS_file) as f:
        #next(f)
        for bin in f:
            line = bin.strip().split()
            IDEAS_dict['\t'.join(line[1:4])] = int(line[4])
            
    ChromHMM_dict = {}
    with open(args.ChromHMM_file) as f:
        for bin in f:
            line = bin.strip().split()
            ChromHMM_dict['\t'.join(line[0:3])] = int(line[3].lstrip('E'))

    res = np.zeros((args.State,args.State))
    for key,value in IDEAS_dict.items():
        if key in ChromHMM_dict.keys():
            res[value,ChromHMM_dict[key] -1] += 1
    

    np.savetxt(args.Outfile, res, delimiter="\t", fmt='%d')

def run():
    parser = argparse.ArgumentParser(usage='python %(prog)s [options]', description=__DESCRIPTION__, epilog=__COPYRIGHT__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--IDEAS', metavar='Segment_file', dest='IDEAS_file', required=False,
                        help='Provide an IDEAS segment file')
    parser.add_argument('--ChromHMM', metavar='Segment_file', dest='ChromHMM_file', required=False,
                        help='Provide a ChromHMM segment file')
    parser.add_argument('-o', dest='Outfile', required=False,
                        help='The path for result file')    
    parser.add_argument('-s', metavar='INT', dest='State', required=False,type=int,
                        help='The state number')     
    args = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    Compare(args)  

    
if __name__ == '__main__':
    run()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    