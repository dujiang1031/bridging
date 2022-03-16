# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time

#change to clonal tracking matrix data directory
#clonal tracking data is stored as tab delimited txt file, where rows are clones, columns are time_points/cell_types/etc
clonal_tracking_dir = r'/media/sf_D/10x_data_analysis/run_12/blood'

def main():
    
    
    
    """
    generate 00-external/TBC-list/{sample}_TBC-list.fasta
    """
    loc = '00-external/TBC-list'
    if not os.path.exists(loc):
        os.mkdir(loc)
    with open('00-external/multiplex_index.txt', 'r') as rmi:
        milines = rmi.readlines()
    for line in milines[1:]:
        line = line.split('\t')
        sample = line[3].strip()
        print sample
        mice = line[5].strip().split(',')
        result = []
        for mouse in mice:
            bloodfile = ''
            for fil in os.listdir(clonal_tracking_dir):
                if mouse == fil.split('_')[0] and fil.endswith('.txt'):
                    bloodfile += os.path.join(clonal_tracking_dir, fil)
            if bloodfile == '':
                print "Blood file of %s not found, skipped." % mouse
                continue
            with open(bloodfile, 'r') as rf:
                bloodlines = rf.readlines()
            for i in range(1, len(bloodlines)):
                tbc = bloodlines[i].split('\t')[0].strip()
                result.append('>TBC' + mouse + '-' + str(i) + '\r\n' + tbc + '\r\n')
        if len(result) == 0:
            print 'No TBCs for %s, skipped.' % sample
            continue
        with open(os.path.join(loc, sample + '_TBC-list.fasta'), 'w') as wf:
            for ll in result:
                wf.writelines(ll)

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))