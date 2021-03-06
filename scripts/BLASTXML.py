# -*- encoding: utf-8 -*-
#Generating BLAST: (v6): all parameters
#Using a .FASTA file, (generated by multiFasta.py) BLAST both files from both progenitors
#Each line contains...
#FORMATDV OPTIONS:
#-i: imput file
#-p: type of file: F=nucleotide
#-o: parsing options: F=not parsing SeqId not creating index

import sys
import os
import multiprocessing

############################# PARAMETERS #############################
num_args = 7

if len(sys.argv) != num_args:
    print("\nERROR calling the script. Number of arguments must be %d\nSyntax: python %s [progenitor1MultiFasta.FASTA] [progenitor2MultiFasta.FASTA] [mode] [p.value] [evalue] [blast mode]\n" %(num_args-1, sys.argv[0]))
    sys.exit(1)

progenitor1= sys.argv[1]
progenitor2 = sys.argv[2]
mode = sys.argv[3]
pvalue = sys.argv[4]
evalue = sys.argv[5]
blastMode = sys.argv[6]
#progenitor1= "PC9_multiFASTA.fasta" 
#progenitor2 = "PC15_multiFASTA.fasta"
#mode = '1'
#pvalue = '80'
#evalue = '0.00000000000000000001'
#blastMode = 1

############################# FUNCTIONS #############################

def worker(progenitor1,  progenitor2, modes, pvalue,  evalue):
    """Blasting muliFastas"""
    name_progenitor1 = progenitor1[0:progenitor1.find('_')]
    name_progenitor2 = progenitor2[0:progenitor2.find('_')]
    for mode in modes:
        cmd = ""
        if mode == 1: cmd = "blast2 -p blastn -m 7 -i %s -d %s -P %s -e %s -v 1 -b 1 > %s*%s.xml" %(progenitor1, progenitor2,  pvalue,  evalue,  name_progenitor1, name_progenitor2)
        elif mode == 2: cmd = "blast2 -p blastn -m 7 -i %s -d %s -P %s -e %s -v 1 -b 1 > %s*%s.xml" %(progenitor2, progenitor1,  pvalue,  evalue, name_progenitor2, name_progenitor1)
        elif mode == 3: cmd = "blast2 -p blastn -m 7 -i %s -d %s -P %s -e %s -v 1 -b 1 > %s*%s.xml" %(progenitor1, progenitor1,   pvalue,  evalue,  name_progenitor1, name_progenitor1)
        elif mode == 4: cmd = "blast2 -p blastn -m 7 -i %s -d %s -P %s -e %s -v 1 -b 1 > %s*%s.xml" %(progenitor2, progenitor2,  pvalue,  evalue,  name_progenitor2, name_progenitor2)
        os.system(cmd)
    return

############################# PROGRAM #############################

if __name__ == '__main__':
    if blastMode == '1':
        blastMode = 3
    elif blastMode == '2':
        blastMode = 5
    cores = multiprocessing.cpu_count()
    #Unless multicore, only one list (for one core)
    listToRun = 1
    if mode == '1': #Convervative. Uses N-1 cores, letting OS and other programs to run normaly
        listToRun = cores - 1
    elif mode == '2': #Boost. Uses ALL cores
        listToRun = cores
    cmd = "formatdb -i %s -p F -o F" % progenitor1
    os.system(cmd)
    cmd = "formatdb -i %s -p F -o F" % progenitor2
    os.system(cmd)
    i = 0
    #Cores-1 processes are launched at the same time
    listProcess = [None] * (listToRun)
    processes = []
    tasks = []
    #Preparing tasks
    for i in range(1, blastMode):
        tasks.append(i)
    #Planing tasks in listToRun lists/cores
    for i in range(listToRun):
        listProcess[i] = []
    for task in tasks:
        listProcess[i % listToRun].append(task)
        i = i + 1
    for list in listProcess:
        processes.append(multiprocessing.Process(target=worker,  args=(progenitor1,  progenitor2, list,  pvalue,  evalue)))
    for process in processes:
        process.start()
    for process in processes:
        process.join()
        
#    print 'Memory usage (Self): %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
#    print 'Memory usage (Children): %s (kb)' % resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

