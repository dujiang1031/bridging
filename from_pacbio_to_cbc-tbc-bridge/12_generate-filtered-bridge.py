# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time
import pandas as pd

def transformFastaToDic(fasta):
    """
    input: fasta-like file
    return: dict{id:seq}
    """
    dic = {}
    with open(fasta, 'r') as rf:
        eof = False
        while True:
            head = rf.readline().strip()
            if not head:
                eof = True
            if eof:
                break
            assert head[0]=='>', "Error: not begin with >."
            seq = rf.readline().strip()
            dic[head[1:]] = seq
    return dic

def collectBarcodeSequenceForSeqid(aligned, barcodes, cutoff):
    """
    input:
        aligned: file, {sample}_aligned-TBC/CBC.txt
        barcodes: dict, {barcode-id:barcode-sequence}
        cutoff: float, score threshold
    return:
        dict{seqid: barcode-sequence}
    """
    result = {}
    with open(aligned, 'r') as ra:
        eof = False
        while True:
            head = ra.readline().strip()
            if not head:
                eof = True
            if eof:
                break
            if not head[0] == '>':
                continue
            seqid = head[1:].split(',')[0].strip()
            bcid = head[1:].split(',')[1].strip()
            score = float(head[1:].split(',')[2].strip())
            if score < cutoff:
                continue
            result[seqid] = barcodes[bcid]
    return result

def main():
    
    
    
    """
    generate 07-cleaned_bridge/{sample}_bridge.txt and clean_bridge_stats.xlsx
        -- delete cells that are mapped to more than one clone in 06
    """
    loc = '07-cleaned_bridge'
    locin = '06-TBC_CBC_bridge'
    if not os.path.exists(loc):
        os.mkdir(loc)
    rows,numbers = [],[]
    for fil in os.listdir(locin):
        if not fil.endswith('bridge.txt'):
            continue
        sample = fil.split('_')[0]
        print sample, '...'
        n1,n2 = 0,0
        droplets = []
        with open(os.path.join(locin, fil), 'r') as rf:
            rawlines = rf.readlines()
        for line in rawlines:
            droplets.append(line.split('\t')[0].strip())
        wf = open(os.path.join(loc, fil), 'w')
        for line in rawlines:
            cell = line.split('\t')[0].strip()
            n1 += 1
            if droplets.count(cell) > 1:
                continue
            wf.writelines(line)
            n2 += 1
        wf.close()
        n3 = 100.0*len([c for c in set(droplets) if droplets.count(c)==1])/len(set(droplets))
        rows.append(sample)
        numbers.append([n1, n2, 100*float(n2)/float(n1), n3])
    df = pd.DataFrame(numbers, index=rows, columns=['#initial bridge', '#singlet bridge', '%singlet bridge', '%singlets'])
    df.to_excel(os.path.join(loc, 'clean_bridge_stats.xlsx'))

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))