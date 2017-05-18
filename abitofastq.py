  

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### scriptname.py
### Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 

import numpy as np
import itertools
import random
import os
import pandas as pd
from Bio import SeqIO
from StringIO import StringIO
from Bio import AlignIO

import glob, os


##Goal: to take snp data and design oligos for site directed mutagenesis on WT vectors. 

abilist = []

with open("example.fasta", "w") as output_handle:
    SeqIO.write("", output_handle, "fastq")

os.chdir("/media/sf_O_DRIVE/Usager/pelletie/scripts/sequencingalignment/mydir")
for file in glob.glob("*.ab1"):
    abilist.append(str(file))


for i in abilist:
    count = SeqIO.convert(i, "abi", (str(i)+".fastq"), "fastq")
    # count = SeqIO.parse(i, "abi", (str(i)+".fastq"), "fastq")
    #print("Converted %i records" % count)


fastqlist = [] 
for file in glob.glob("*.fastq"):
    fastqlist.append(str(file))

corid = []
for i in fastqlist:
    corid.append("corr/"+i)


for i in xrange(len(fastqlist)):
    with open(fastqlist[i]) as original, open(corid[i], 'w') as corrected:
        records = SeqIO.parse(fastqlist[i], 'fastq')
        for record in records:
                       
            record.id = fastqlist[i].replace(".fastq","") 
            print record.id 
        SeqIO.write(record, corrected, 'fastq')

with open('fastqSeq_date.txt', 'w') as outfile:
    for fname in corid:
        with open(fname) as infile:
            outfile.write(infile.read() + "\n")