#!/usr/bin/env python
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from collections import Counter
import collections
import itertools


z = {}

listfile = glob.glob("*fas")
list_sample = open("list", "r")
list_sample = list_sample.read().splitlines()
a=len(listfile)
for matrix in listfile:
	d = {}
	for record in SeqIO.parse(matrix, 'fasta'):
		length = len(record.seq)
		recordSEQ = record.seq
		recordid = record.id.split('.')[0]
		nCount = recordSEQ.lower().count('-') 
		ratio = ((nCount/(length*1.00)))
		
		d[recordid] = ratio	
	z = {x: d.get(x, 0) + z.get(x, 0) for x in set(d).union(z)}


print "percentage of missing data per sample"
for key, value in z.items():
    z[key] = float(format(value/a, '.4f'))
for kv in (sorted(z.items(), key=lambda x:x[1])):
	print kv[0],'\t',kv[1]
