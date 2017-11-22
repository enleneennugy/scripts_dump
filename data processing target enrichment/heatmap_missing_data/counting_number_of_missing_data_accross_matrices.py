#!/usr/bin/env python
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

import pandas as pd


z = {}
df = pd.DataFrame()

listfile = glob.glob("*fasta")

a = len(listfile)


for matrix in listfile:
	d = {}
	indexd = matrix.split('.')[0]
	for record in SeqIO.parse(matrix, 'fasta'):
		length = len(record.seq)
		recordSEQ = record.seq
		recordid = record.id.split('.')[0]
		nCount = recordSEQ.lower().count('-')
		ratio = ((nCount/(length*1.00)))
		
		d[recordid] = ratio	
	df1 = pd.DataFrame([d], index=[indexd])
	df = pd.concat([df,df1], axis=0)
	z = {x: d.get(x, 0) + z.get(x, 0) for x in set(d).union(z)}

df.to_csv("missing_data.csv", sep='\t')
print "percentage of missing data per sample"
for key, value in z.items():
	z[key] = float(format(value/a, '.4f'))
for kv in (sorted(z.items(), key=lambda x:x[1])):
	print kv[0],'\t',kv[1]
