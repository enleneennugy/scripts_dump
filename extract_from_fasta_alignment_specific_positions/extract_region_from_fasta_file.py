#!/usr/bin/env python2

import glob
import time
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#aligner
import subprocess


#make a new directory for the analyses, day/time

tm = (time.strftime("%H.%M.%S"))
dt = (time.strftime("%Y_%m_%d"))
now = (dt + "_" + tm)
os.mkdir(now)
os.chdir(now)
list_genes = open('list.gene.formated', 'r').read().splitlines()

for genes in list_genes:

	name = genes.split()[0]
	first = int(genes.split()[1])
	second = int(genes.split()[2])
	f = open((name+".fasta"), "w")
	for record in SeqIO.parse(open('plastomes.fasta', 'r'), 'fasta'):
		#print record.id
#		print record.seq[(first-1):(second-1)]
		f.write('>'+str(record.id) + "\n")
		f.write(str(record.seq[(first-1):(second-1)])+'\n')