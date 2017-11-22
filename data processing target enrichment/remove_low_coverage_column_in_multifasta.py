#!/usr/bin/env python
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from Bio import AlignIO
from pandas import Series,DataFrame
import pandas as pd
import numpy as np



df = pd.DataFrame()

listfile = glob.glob("*fasta")
n = 0
fasta = {}

for matrix in listfile:
	
	alignment = AlignIO.read(matrix, "fasta")
	n = n+1
	number_samples = len(alignment)
	f = open("cleaning/"+matrix, 'w')
	print matrix

	for record in SeqIO.parse(matrix, 'fasta'):
		sequence = str(record.seq)
		sequence_id = str(record.id)
		fasta[sequence_id] = (sequence)
	d = fasta
	for key, value in d.iteritems():
		d[key] = list(value)

	df = pd.DataFrame.from_dict(d, orient='index')
	print ("Before:" + str(df.shape))

	for col in df:
			val = (df[col].value_counts("-"))*number_samples
			val = val.to_frame()
			if "-" in val.index:
				b = str(val.loc['-'])
				b = float(b.split()[1])
				if (b) > (number_samples-4):
					df.drop(col, axis=1, inplace=True)

			else:
				continue
	print ("After: " + str(df.shape))

	df.replace(['None', 'nan', None], '-', inplace=True)
	print df
	df['new'] = df.values.sum(axis=1)
	dic = dict(zip(df.index, df.new))

	for key, value in dic.iteritems():
		seq_ID = str(key)
		seq = str(value)
		new_record = SeqRecord(Seq(seq), id= seq_ID,description="")
		SeqIO.write(new_record, f, "fasta")



print "number of matrices processed: " + str(n)