# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import sys
import time
import tables

#change to expression matrix location
#by default the organization is as {path_to_run_folder}/{run_name}/{sample_name}/out/filtered_gene_bc_matrices_h5.h5
path_to_run_folder = '/media/sf_D/10x_data_analysis/'
matrix_file = 'filtered_gene_bc_matrices_h5.h5'

def main():
    
    
    
    """
    generate 00-external/CBC-list/{index}_CBC-list.fasta
    """
    with open('00-external/multiplex_index.txt', 'r') as rmi:
        milines = rmi.readlines()
    if not os.path.exists('00-external/CBC-list'):
        os.mkdir('00-external/CBC-list')
    for line in milines[1:]:
        line = line.split('\t')
        print line[3].strip(), '...' #sample
        out = '00-external/CBC-list/' + line[3].strip() + '_CBC-list.fasta'
        if os.path.isfile(out):
            print 'Already exists, skipped.'
            continue
        path = path_to_run_folder + '/' + line[4].strip() + '/' + line[3].strip() + '/out/' + matrix_file
        if not os.path.isfile(path):
            print "Cellranger H5 file not found, skipped."
            continue
        with tables.open_file(path, 'r') as f:
            try:
                group = f.get_node(f.root, 'mm10')
            except tables.NoSuchNodeError:
                try:
                    group = f.get_node(f.root, 'hg19')
                except tables.NoSuchNodeError:
                    print "That genome does not exist in this file."
                    sys.exit()
            barcodes = getattr(group, 'barcodes').read()
        with open(out, 'w') as wf:
            for j in xrange(len(barcodes)):
                wf.writelines('>CBC' + str(j+1) + '\r\n' + barcodes[j][:16] + '\r\n')

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))