# -*- coding: utf-8 -*-

import glob,getopt,os,sys
import datetime


def process_file(file):
    outname = "avg_" + file
    f = open(file)
    out = open(outname,"w+")
    next(f)
    gene_dict = {}
    times_dic = {}
    for gene in f:
        gene_pos = gene.rstrip().split("\t")  
        count_sum = sum(map(lambda x: float(x), gene_pos[1:len(gene_pos)]))
        if gene_pos[0] in gene_dict.keys():
            times_dic[gene_pos[0]] += 1
            gene_dict[gene_pos[0]] += count_sum
        else:
            times_dic[gene_pos[0]] = 1
            gene_dict[gene_pos[0]] = count_sum
    for k in sorted(gene_dict.keys()):
        out.write(str(k)+'\t'+str(int(gene_dict[k]/times_dic[k]))+'\n')
    out.close()
    print('Done averaging transcripts for ***' + file, end='')

def get_targetfile(input_path):
    os.chdir(input_path)
    for file in glob.glob('*TSS*bp.txt'):
        process_file(file)

def main():
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:')
    input_path = '.'
    for op, value in opts:
        if op == '-i':
            input_path = value
            
    start_time = datetime.datetime.now()
    get_targetfile(input_path)
    end_time = datetime.datetime.now()
    print('*** after ' + str((end_time - start_time).seconds) + ' seconds!')
          
if __name__ == '__main__':
    main()
    
    
    
    
    