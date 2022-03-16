# -*- coding: utf-8 -*-
"""
Created on Mon Jul 09 14:26:42 2018

@author: Du Jiang

"""

import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import WaterCommandline

def readChunkIterator(fil, ref):
    with open(fil, 'r') as rf:
        eof = False
        while True:
            result = []
            for iii in xrange(500):
                head = rf.readline().strip()
                if not head:
                    eof = True
                    break
                read = rf.readline().strip()
                result.append([head,read])
            yield (result,ref)
            if eof:
                break

def getCellBarcodeAlignment(read, fil):
    """
    use stdin and stdout to simplify water
        asequence: one SMRT read
        bsequence: {index}_CBC-list.fasta
        return: best matched CBC for this SMRT read and the corresponding score
    """
    water_cline = WaterCommandline(asequence='stdin', filter=True, bsequence=fil, gapopen=10.0, gapextend=.5, stdout=True)
    child = subprocess.Popen(str(water_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
    rec = SeqRecord(Seq(read), id="temp")
    SeqIO.write(rec, child.stdin, "fasta")
    child.stdin.close()
    seqs,scores = [],[]
    line = child.stdout.readline()
    eof = False
    while True:
        if not line:
            eof = True
        if eof:
            break
        if '2:' in line[:8]:
            seqs.append(line.strip().split(':')[1])
        elif 'Score' in line:
            scores.append(float(line.split(':')[1]))
        line = child.stdout.readline()
    assert len(seqs)==len(scores), "ERROR: incorrect alignment file line counting."
    return seqs[scores.index(max(scores))], max(scores)

def alignReadChunkToCBCList(PAIR):
    chunk, reference=PAIR
    result = []
    for read in chunk:
        head,seq = read
        cbcid,score = getCellBarcodeAlignment(seq, reference)
        result.append(head + ', ' + cbcid + ', ' + str(score) + '\r\n' + seq + '\r\n')
    return result

def main():
    
    loc = '05-align_TBC'
    if not os.path.exists(loc):
        os.mkdir(loc)
    with open('00-external/multiplex_index.txt', 'r') as rmi:
        milines = rmi.readlines()
    mice,fillist1,fillist2 = [],[],[]
    for miline in milines[1:]:
        k = miline.split('\t')[3].strip()
        fil1 = '00-external/TBC-list/'+k+'_TBC-list.fasta'
        fil2 = '03-pairs/'+k+'_TBC.fasta'
        mice.append(k)
        fillist1.append(fil1)
        fillist2.append(fil2)
    pairs = [(i,i) for i in range(len(mice))]
    for a,b in pairs:
        tail = '_aligned-TBC.txt'
        if not os.path.isfile(fillist1[a]) and os.path.isfile(fillist2[b]):
            continue
        if os.path.exists(os.path.join(loc, mice[b]+tail)):
            continue
        print 'Writing:', mice[b]+tail
        pool = mp.Pool(mp.cpu_count())
        results = pool.imap(alignReadChunkToCBCList, readChunkIterator(fillist2[b], fillist1[a]))
        pool.close()
        pool.join()
        with open(os.path.join(loc, mice[b]+tail), 'w') as wf:
            for result in results:
                for line in result:
                    wf.writelines(line)
    
    #####################################################################################################
    
    cut = 170
    cc = 'center'
    ff = 'figure fraction'
    tag = '_aligned'
    scores = {}
    for fil in os.listdir(loc):
        if not fil.endswith('TBC.txt'):
            continue
        if not tag in fil:
            continue
        sample = fil.split('_')[0]
        print sample, '...'
        scores[sample] = []
        with open(os.path.join(loc, fil), 'r') as rf:
            eof = False
            while True:
                head = rf.readline().strip()
                if not head:
                    eof = True
                if eof:
                    break
                if not head[0] == '>':
                    continue
                scores[sample].append(float(head.split(',')[-1]))
    assert len(scores)<=16, "ERROR: too many samples."
    fig = plt.figure(figsize=(15,12))
    i = 0
    for sample in scores.keys():
        i += 1
        temp = [scores[sample][j] for j in range(len(scores[sample])) if scores[sample][j] >= cut]
        st = str(int(1000.0*len(temp)/len(scores[sample]))/10.0) + '%'
        ax = fig.add_subplot(4, 4, i)
        ax.hist(scores[sample], bins=50, color='b')
        ax.set_title(sample, fontsize=18)
        ax.set_xlim(right=255)
        ax.tick_params(axis='both', labelsize=15)
        ax.axvline(x=cut-1, color='r', ls='--', lw=1)
        ax.text(.85, .9, st, ha=cc, va=cc, transform=ax.transAxes, color='r', fontsize=15)
    plt.annotate('# Reads', (.03,.5), xycoords=ff, ha=cc, va=cc, fontsize=25, rotation=90)
    plt.annotate('Alignment Score', (.5,.03), xycoords=ff, ha=cc, va=cc, fontsize=25)
    plt.subplots_adjust(left=.08, bottom=.08, right=.95, top=.95, hspace=.3, wspace=.3)
    fig.savefig(os.path.join(loc, 'score-distribution'+tag+'.png'), transparent=True)
    plt.close(fig)


if __name__ == '__main__':
    main()