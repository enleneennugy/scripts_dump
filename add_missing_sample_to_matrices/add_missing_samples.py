#!/usr/bin/env python
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from collections import Counter

#in target enrichment when some marker are missing for some sample, in order to have the same sample in all the matrices, this script is generating empty sequences for these samples. 


#to make the list of all used samples
#cat *.fas | grep ">" | tr "." "\t" | cut -f1 | sort | uniq | sort |sed 's/>//'>list


list_sample = open("list", "r")
list_sample = list_sample.read().splitlines()
missing_samples_list = []
listfile = glob.glob("*fas")
for matrix in listfile:
	matrix_list_sample = []
	for record in SeqIO.parse(matrix, 'fasta'):
		length_matrix = len(record)
		record = record.id.split('.')[0]
		matrix_list_sample.append(record)
	missing_samples = list(set(list_sample) - set(matrix_list_sample))
	missing_samples_list.extend(missing_samples)

	for samplem in missing_samples:
		new_name = samplem + "." + matrix.split(".")[0]
		missing_sequence = length_matrix * "-"
		new_record = SeqRecord(Seq(missing_sequence), id=new_name, description='')
		
		with open(matrix, 'a') as fileout:
			SeqIO.write(new_record, fileout, "fasta")
		fileout.close()
d = Counter(missing_samples_list)
print "number of time a sample has been replaced:"
for kv in (sorted(d.items(), key=lambda x:x[1])):
	print kv[0],'\t',kv[1]