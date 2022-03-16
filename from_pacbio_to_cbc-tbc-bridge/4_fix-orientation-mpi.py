# -*- coding: utf-8 -*-
"""
Created on Mon Jul 09 14:42:14 2018

@author: Du Jiang

generate 01-fix_orientation/{fragment}/{chunk}.txt
    -- done on HPC interactive node:
        $ salloc --ntasks={number of chunks} --time=5:00:00 --constraint=avx
        $ source /usr/usc/python/2.7.8/setup.sh
        $ source /usr/usc/openmpi/default/setup.sh
        $ source /home/rcf-proj/ac2/dujiang/software/emboss/sourceme.sh
        $ mpirun python 01-fix_orientation_mpi.py

"""

import os
from mpi4py import MPI
from Bio.Emboss.Applications import WaterCommandline

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    loc = '01-fix_orientation'
    chunk = 'chunk-'+'{:03}'.format(rank+1)
    if not os.path.exists(loc):
        os.mkdir(loc)
    for fragment in ('inside', 'inside_rc', 'outside', 'outside_rc'):
        subloc = loc + '/' + fragment
        if not os.path.exists(subloc):
            os.mkdir(subloc)
        infile = os.path.join('00-external', fragment+'.txt')
        data = os.path.join('00-data', chunk +'.fasta')
        if not os.path.isfile(data):
            continue
        outfile = os.path.join(subloc, chunk+'.txt')
        water_cline = WaterCommandline(asequence=infile, bsequence=data, gapopen=10.0, gapextend=.5, outfile=outfile)
        stdout,stderr = water_cline()

if __name__ == '__main__':
    main()