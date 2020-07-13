# -*- coding: utf-8 -*-
"""
Created on Tue May 22 20:22:45 2018

@author: Administrator
"""
Description = '''
This scripts can help you calculate motif change given differential peaks.
NOTE:you should provide a file formatted as 'peak_name'	'chr'	'start'	'end' for reference peaks in the working dir.
'''
Copyright = ''


import glob,argparse,csv


def process_file(args):
            
    peak = open(args.strain_name+'_peakset.txt')
    peak_set_dict = {}
    next(peak)
    for peak_line in peak:
        peak_tmp = peak_line.rstrip().split()
        peak_set_dict[str(peak_tmp[1])+str(peak_tmp[2])+str(peak_tmp[3])] = peak_tmp[0]
    peak.close()
    
    out = open(args.strain_name + '_changed.txt','w+')
    
    #write header line
    out.write('motif' + '\t')
    for diffpeak_file in glob.glob(args.strain_name +'*_500.txt' ):
        dp_name =diffpeak_file.lstrip(args.strain_name).rstrip('_500.txt')
        out.write(dp_name + '_up' + '\t' + dp_name + '_down' + '\t')
    out.write('\n')
    
    #for each motif 
    for motif_file in glob.glob(args.motif_dir + '/' + args.strain_name + '_A_Con0_*.txt'):
        motif  = motif_file.replace(args.motif_dir + args.strain_name +'_A_Con0_','').rstrip('.txt')
        print('Start ' + motif + ' !')
        out.write(motif + '\t')
        
        mf = open(motif_file)
        next(mf)
        mf_dict = {}
        for mf_line in mf:
            mf_tmp = mf_line.rstrip().split()
            if mf_tmp[0] in mf_dict.keys():
                mf_dict[mf_tmp[0]]  += 1
            else:
                mf_dict[mf_tmp[0]]  = 0
        
        #calculate each diffpeak file           
        for diffpeak_file in glob.glob(args.strain_name +'*_500.txt' ):
            dp_name =diffpeak_file.lstrip(args.strain_name).rstrip('_500.txt')
            dp_dict = {}
            dp = open(diffpeak_file) 
            next(dp)
            
            for dp_line in dp:
                dp_tmp = dp_line.rstrip().split()
                peak_key = peak_set_dict[str(dp_tmp[1]).lstrip('"').rstrip('"')+str(dp_tmp[2])+str(dp_tmp[3])]
                if peak_key in mf_dict.keys():
                    if float(dp_tmp[10]) < 0.05:
                        if float(dp_tmp[9]) > 0:
                            dp_dict[peak_key] = 1
                        elif float(dp_tmp[9]) < 0:
                            dp_dict[peak_key] = -1
                        else:
                            dp_dict[peak_key] = 0
            up = sum(v > 0 for v in dp_dict.values())
            down = sum(v < 0 for v in dp_dict.values())
            out.write(str(up) + '\t' + str(down) + '\t')
            print('\t' + diffpeak_file)
        out.write('\n')
        print('Done ' + motif + ' !')
        mf.close()
    out.close()
    print('*' * 10 + args.strain_name + ' done! and the result has been output into ' + args.strain_name + '_changed.txt')
    
   
def main():   
    parser = argparse.ArgumentParser(prog = 'motif_multiply_occupancy.py', usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--motif_dir','-m', action = 'store', type=str, dest = 'motif_dir',required=True,default = './',
                        help = 'a directory contains homer calculated motif files')
    parser.add_argument('--diffpeak_dir','-p',action = 'store', type=str, dest = 'diffpeak_dir',required=True,default = './',
                        help = 'a directory contains DiffBind result diffpeak file')
    parser.add_argument('--strain','-s', action = 'store', type=str, dest = 'strain_name',required=True,
                        help = 'a directory contains featurecount result peakcount file')
    args = parser.parse_args()
    
                   
    process_file(args)

if __name__ == '__main__':
    main()
