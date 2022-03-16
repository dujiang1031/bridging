# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from Levenshtein import distance

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

def identifyMultiplexIndex(header, indices):
    """
    input:
        header = 5'-12nt + rc(3'-12nt)
        indices = [primer IDs]
    """
    dist = {}
    for key in indices:
        perfect = key + key
        dist[key] = [distance(header[:12], perfect), distance(header[12:], perfect)]
    for key in dist.keys():
        if dist[key][0] > 4 and dist[key][1] > 4:
            del dist[key]
    if not dist:
        return -1
    else:
        temp = []
        for tup in dist.items():
            if bool(set(tup[1]) & set([0, 1, 2])):
                temp.append(tup)
        if not temp:
            return min(dist, key=lambda c:c[1]+c[2])
        else:
            return min(temp, key=lambda c:min(c[1]))[0]

def main():
    
    
    
    """
    generate 02-demultiplex/:
        {sample}.txt
        demultiplex_stats.xlsx
        length_distribution.png
    """
    cc = 'center'
    ff = 'figure fraction'
    loc = '02-demultiplex'
    if not os.path.exists(loc):
        os.mkdir(loc)
    with open('00-external/multiplex_index.txt', 'r') as ri:
        lines = ri.readlines()
    index = []
    wf = {}
    count = {}
    pooled = {}
    length = {}
    titles = {}
    for i in range(1, len(lines)):
        key = lines[i].split('\t')[1].strip()
        mouse = lines[i].split('\t')[3].strip()
        index.append(key)
        titles[key] = mouse
        wf[key] = open(os.path.join(loc, mouse+'.fasta'), 'w')
        count[key] = 0
        pooled[key] = float(lines[i].split('\t')[-1].strip())
        length[key] = []
    count['unknown'] = 0
    length['unknown'] = []
    for fil in os.listdir('01-fix_orientation/assorted'):
        chunk = fil.split('.')[0]
        print 'Demultiplexing', chunk, '...'
        with open(os.path.join('01-fix_orientation/assorted', fil), 'r') as rf:
            eof = False
            while True:
                head = rf.readline().strip()
                if not head:
                    eof = True
                if eof:
                    break
                line = rf.readline().strip()
                header = line[:12] + reverseComplement(line[-12:])
                destination = identifyMultiplexIndex(header, index)
                if destination == -1:
                    count['unknown'] += 1
                    length['unknown'].append(len(line))
                else:
                    wf[destination].writelines(head + '\r\n' + line + '\r\n')
                    count[destination] += 1
                    length[destination].append(len(line))
    for key in wf.keys():
        wf[key].close()
    print count, sum([count[key] for key in count.keys()])
    fig = plt.figure(figsize=(15,12))
    for i in range(len(index)):
        k = index[i]
        ax = fig.add_subplot(5, 4, i+1)
        ax.hist(length[k], bins=80, color='blue')
        ax.set_title(titles[k], fontsize=12)
        if not len(length[k]) == 0:
            ax.axvline(min(length[k]), color='green', ls='--')
            ax.axvline(max(length[k]), color='green', ls='--')
    fig.add_subplot(5, 4, 18).hist(length['unknown'], bins=80, color='orange')
    fig.subplots_adjust(left=.08, right=.96, bottom=.08, top=.92, hspace=.4)
    plt.annotate('# Reads', (.04,.5), xycoords=ff, ha=cc, va=cc, rotation=90, fontsize=20)
    plt.annotate('Size (bp)', (.5,.04), xycoords=ff, ha=cc, va=cc, fontsize=20)
    fig.savefig(os.path.join(loc, 'length_distrubtion.png'), transparent=True)
    plt.close(fig)
    total_ng = sum([pooled[key] for key in pooled.keys()])
    total_rd = sum([count[key] for key in count.keys()])
    total_rx = sum([count[key] for key in count.keys() if not key == 'unknown'])
    rows,numbers = [],[]
    for k in count.keys():
        if k == 'unknown':
            continue
        rows.append(k)
        pct1 = 100.0*count[k]/total_rd
        pct2 = 100.0*pooled[k]/total_ng
        pct3 = 100.0*count[k]/total_rx
        numbers.append([count[k], pct1, pct3, pct2])
    rows.append('unknown')
    numbers.append([count['unknown'], 100.0*count['unknown']/total_rd, None, 0])
    cols = ('reads#', 'reads%', "reads%-among-id'ed", 'pooling%')
    pd.DataFrame(numbers, index=rows, columns=cols).to_excel(os.path.join(loc, 'demultiplex_stats.xlsx'))
    

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))