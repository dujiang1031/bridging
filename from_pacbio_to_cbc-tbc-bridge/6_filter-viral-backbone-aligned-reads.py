# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time

def getScoreList(fil):
    """
    input: alignment result in water format
    return: a list of alignment scores
    """
    output = []
    with open(fil, 'r') as rf:
        eof = False
        while True:
            line = rf.readline()
            if not line:
                eof = True
            if eof:
                break
            if not 'Score' in line[:10]:
                continue
            output.append(float(line.split(':')[1]))
    return output

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
    genenerate 01-fix_orientation/assorted/{chunk}.fasta
        --save reads with expected 1kb sequence in consistent orientation
    """
    loc = '01-fix_orientation/assorted'
    if not os.path.exists(loc):
        os.mkdir(loc)
    n_tot,n_pos,n_neg = 0,0,0
    for fil in os.listdir('00-data'):
        chunk = fil.split('.')[0]
        scores1 = getScoreList(os.path.join('01-fix_orientation/inside', chunk+'.txt'))
        scores2 = getScoreList(os.path.join('01-fix_orientation/inside_rc', chunk+'.txt'))
        pos,neg = [],[]
        assert len(scores1)==len(scores2), "ERROR: unequal numbers of reads."
        for j in range(len(scores1)):
            if scores1[j] >= 2000 and scores1[j] > scores2[j]:
                pos.append(j)
            elif scores2[j] >= 2000 and scores2[j] > scores1[j]:
                neg.append(j)
        wf = open(os.path.join(loc, chunk+'.fasta'), 'w')
        with open(os.path.join('00-data', fil), 'r') as rf:
            eof = False
            click = 0
            while True:
                click += 1
                head = rf.readline().strip()
                if not head:
                    eof = True
                if eof:
                    break
                assert head[0] == '>', "ERROR: incorrect line counting."
                n_tot += 1
                read = rf.readline().strip()
                if click-1 in pos:
                    n_pos += 1
                elif click-1 in neg:
                    read = reverseComplement(read)
                    n_neg += 1
                else:
                    continue
                wf.writelines(head + '\r\n' + read + '\r\n')
        wf.close()
        print 'accumulated:', chunk, n_tot, n_pos, n_neg
    print 'Final:', n_tot, n_pos, n_neg

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))