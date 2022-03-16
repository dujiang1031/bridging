# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 22:08:45 2020

@author: Du Jiang
"""

import os
import time

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
    generate 06-TBC_CBC_bridge/{sample}_bridge.txt
    """
    loc = '06-TBC_CBC_bridge'
    cbc_cutoff = 79
    tbc_cutoff = 249
    if not os.path.exists(loc):
        os.mkdir(loc)
    with open('00-external/multiplex_index.txt', 'r') as rmi:
        milines = rmi.readlines()
    files = {}
    for miline in milines[1:]:
        mouse = miline.split('\t')[3].strip()
        fil1 = '04-align_CBC/' + mouse + '_aligned-CBC.txt'
        fil2 = '05-align_TBC/' + mouse + '_aligned-TBC.txt'
        if not (os.path.isfile(fil1) and os.path.isfile(fil2)):
            continue
        fil3 = '00-external/CBC-list/' + mouse + '_CBC-list.fasta'
        fil4 = '00-external/TBC-list/' + mouse + '_TBC-list.fasta'
        files[mouse] = [fil1, fil2, transformFastaToDic(fil3), transformFastaToDic(fil4)]
    for sample in files.keys():
        print sample, '...'
        molecule = {}
        seqid_cbc = collectBarcodeSequenceForSeqid(files[sample][0], files[sample][2], cbc_cutoff)
        seqid_tbc = collectBarcodeSequenceForSeqid(files[sample][1], files[sample][3], tbc_cutoff)
        for seqid in seqid_cbc.keys():
            if not seqid in seqid_tbc.keys():
                continue
            cbc = seqid_cbc[seqid]
            tbc = seqid_tbc[seqid]
            try:
                molecule[cbc, tbc] += 1
            except KeyError:
                molecule[cbc, tbc] = 1
        with open(os.path.join(loc, sample+'_bridge.txt'), 'w') as wb:
            for k1,k2 in molecule.keys():
                wb.writelines(k1 + '\t' + k2 + '\t' + str(molecule[k1, k2]) + '\r\n')

if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))