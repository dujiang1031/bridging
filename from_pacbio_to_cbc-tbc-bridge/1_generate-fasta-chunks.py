# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time
from Bio import SeqIO

#input PacBio fastq location
ccsfastq = '/media/sf_D/sequencing_results/txg/SMRT_library/Irvine1/data/CCS/ccs.fastq'

def main():
    
    
    
    """
    generate 00-data/chunk-{number}.fasta
    """
    loc = '00-data'
    if not os.path.exists(loc):
        os.mkdir(loc)
    n1,n2 = 0,0
    with open(ccsfastq, 'r') as rf:
        iterator = SeqIO.parse(rf, 'fastq')
        eof = False
        seqID = 0
        chunk = 1
        while True:
            wf = open(os.path.join(loc, 'chunk-'+'{:03}'.format(chunk)+'.fasta'), 'w')
            for ii in range(5000):
                try:
                    record = next(iterator)
                except StopIteration:
                    eof = True
                    break
                read = record.seq
                seqID += 1
                if len(read) <= 1000:
                    n1 += 1
                else:
                    n2 += 1
                    wf.writelines('>seqid' + str(seqID) + '\r\n' + read + '\r\n')
            wf.close()
            chunk += 1
            print 'Finshed one chunk...'
            if eof:
                break
    print "# reads shorter than 1000bp:", n1
    print "# reads longer than 1000bp:", n2

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))