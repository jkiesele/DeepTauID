#!/usr/bin/env python2


from argparse import ArgumentParser
import os
import subprocess

print("This script is still experimental and not fully completed")

parser = ArgumentParser('merge samples')
parser.add_argument('nsamples')
parser.add_argument('outdir')
parser.add_argument('--batchdir',type=str,default="")
parser.add_argument('infiles', metavar='N', nargs='+',
                    help='sample list files')

args = parser.parse_args()

if not os.path.isdir(args.outdir):

    allins=''
    for l in args.infiles:
        allins+=' '+l
        
    os.system('createMergeList '+str(args.nsamples)+' '+args.outdir+' '+allins)
    
    
#read number of jobs
file=open(args.outdir+'/nentries','r')
nJobs=file.read()

listtoberun=[]
listsucc=[]

for j in range(int(nJobs)):
    
    if os.path.exists(args.outdir+'/'+str(j)+'.succ'):
        listsucc.append(j)
        continue
    
    listtoberun.append(j)

print('successful: ',listsucc)
print('remaining: ',int(nJobs)-len(listsucc))

cwd = os.getcwd()
cmssw_base = os.environ.get('CMSSW_BASE')

if len(args.batchdir):
    args.outdir = os.path.abspath(args.outdir)
    os.system('mkdir -p '+args.batchdir)
    args.batchdir = os.path.abspath(args.batchdir)
    for j in listtoberun:
        batchscript=''' sleep $(shuf -i1-120 -n1) ; \
cd {cmssw_dir}; \
eval \`scramv1 runtime -sh\`; \
cd {cwd} ; \
merge {outdir}/mergeconfig {jobno} ; '''.format(
    outdir=args.outdir,
    cmssw_dir=cmssw_base,
    cwd=cwd,
    jobno=j,
    batchdir=args.batchdir)
        command = "cd "+ args.batchdir +" ; echo \""+ batchscript +"\" | bsub -q 1nh -J merge"+str(j)
        print(command)
        os.system(command)
    exit()


import multiprocessing as mp


def worker(j):
    print('starting '+str(j))
    os.system('merge '+args.outdir+'/mergeconfig '+str(j))


pool = mp.Pool(processes=mp.cpu_count(),) 
pool.map(worker, listtoberun)


for j in range(int(nJobs)):
    if os.path.exists(args.outdir+'/'+str(j)+'.succ'):
        listsucc.append(j)
    
if len(listsucc) == int(nJobs):
    print('merge successful, creating file list')
    file=open(args.outdir+'/samples.txt','w')
    for filenumber in listsucc:
        file.write('ntuple_merged_'+str(filenumber)+'.root\n')
    file.close()






