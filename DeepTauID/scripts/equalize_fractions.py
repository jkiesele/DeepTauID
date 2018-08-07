#!/usr/bin/env python2



import multiprocessing as mp

from argparse import ArgumentParser
import os
import subprocess
import time


parser = ArgumentParser('equalize tau and non tau fraction. run in same dir as files are located in')
parser.add_argument('infile')
parser.add_argument('outdir')

args = parser.parse_args()

args.outdir = os.path.abspath(args.outdir)

rootfiles=[]
with open(args.infile) as f:
    rootfiles = f.readlines()



def worker(j):
    while not os.path.isfile(rootfiles[j][:-1]):
        counter+=1
        time.sleep(1)
        if counter > 120:
            print('file '+rootfiles[j][:-1]+' failed')
            return
    command='removeNonTauFraction '+rootfiles[j][:-1]+' '+args.outdir+'/'+rootfiles[j][:-1]+'.equal.root 0.2'
    print(command) #leave 4 times more background
    os.system(command)
    
    
pool = mp.Pool(processes=int(mp.cpu_count()*0.8), )
pool.map(worker, range(len(rootfiles)))