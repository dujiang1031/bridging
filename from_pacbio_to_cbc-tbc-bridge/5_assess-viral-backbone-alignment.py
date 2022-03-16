# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time
import matplotlib.pyplot as plt

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

def main():
    
    
    
    """
    generate 01-fix_orientation/score_distribution.png
    """
    loc = '01-fix_orientation'
    cc = 'center'
    ff = 'figure fraction'
    fig = plt.figure(figsize=(12,12))
    for i,sub in enumerate(['inside', 'outside']):
        print "Working on", sub, ':'
        scores = []
        for fil1 in os.listdir(loc+'/'+sub):
            fil2 = os.path.join(loc+'/'+sub+'_rc', fil1)
            if not os.path.isfile(fil2):
                continue
            print fil1, '...'
            scores1 = getScoreList(os.path.join(loc+'/'+sub, fil1))
            scores2 = getScoreList(fil2)
            assert len(scores1)==len(scores2), "ERROR: unequal numbers of reads."
            for j in xrange(len(scores1)):
                scores.append(max([scores1[j], scores2[j]]))
        large = [item for item in scores if item >= 2000]
        percent = int(1000.0*len(large)/len(scores))/10.0
        s1 = str(100-percent) +'%'
        s2 = str(percent) + '%, n='+str(len(large))
        ax = fig.add_subplot(2,1,i+1)
        ax.hist(scores, bins=100, range=(0,5000), color='r')
        ax.axvline(x=2000, color='b', ls='--')
        ax.text(.1, .9, s1, transform=ax.transAxes, ha=cc, va=cc, fontsize=20, color='g')
        ax.text(.8, .9, s2, transform=ax.transAxes, ha=cc, va=cc, fontsize=20, color='g')
        ax.set_title('Aligned by '+sub+' fragment', fontsize=20)
        ax.tick_params(axis='both', labelsize=15)
    plt.xlabel('Alignment Score', fontsize=25)
    plt.annotate('Number of Reads', (.025,.5), xycoords=ff, ha=cc, va=cc, rotation=90, fontsize=25)
    plt.savefig(os.path.join(loc, 'score_distribution.png'), transparent=True)
    plt.close(fig)

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))