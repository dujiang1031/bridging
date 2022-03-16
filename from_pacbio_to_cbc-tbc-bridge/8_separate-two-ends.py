# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time

def reverseComplement(seq):
    out = ''
    for i in range(len(seq)):
        item = seq[len(seq)-i-1]
        if item == 'A':
            out += 'T'
        elif item == 'T':
            out += 'A'
        elif item == 'C':
            out += 'G'
        elif item == 'G':
            out += 'C'
        else:
            out += item
    return out


def main():
    
    
    
    """
    generate 03-pairs/{sample}_CBC.fasta and 03-pairs/{sample}_TBC.fasta
    """
    loc = '03-pairs'
    if not os.path.exists(loc):
        os.mkdir(loc)
    locin = '02-demultiplex'
    for fil in os.listdir(locin):
        if not fil.endswith('.fasta'):
            continue
        sample = fil.split('.')[0]
        print sample, '...'
        wf1 = open(os.path.join(loc, sample+'_TBC.fasta'), 'w')
        wf2 = open(os.path.join(loc, sample+'_CBC.fasta'), 'w')
        with open(os.path.join(locin, fil), 'r') as rf:
            eof = False
            while True:
                head = rf.readline().strip()
                if not head:
                    eof = True
                if eof:
                    break
                line = rf.readline().strip()
                wf1.writelines(head + '\r\n' + line[:130] + '\r\n')
                wf2.writelines(head + '\r\n' + reverseComplement(line[-130:]) + '\r\n')
        wf1.close()
        wf2.close()

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))