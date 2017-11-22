#!/usr/bin/env python2
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys

fasta_file = sys.argv[1]
#fasta_file = "270928382.alignment.fasta"
new_mane_list = open("new_name.list","r")
f = open("rename/"+fasta_file, 'w')

for seq_record in SeqIO.parse(fasta_file, "fasta"):
	sequences = {}
	sequence = str(seq_record.seq)
	print seq_record.id
	ID = seq_record.id.split("_")[0]
	print ID
	new_mane_list = open("new_name.list","r")
	for newname in new_mane_list:
		newname = newname.split("\n")[0]
		newname_split = newname.split("_")[-1]
		if ID == newname_split:
			new_record = SeqRecord(Seq(sequence), id= newname,description="")
			print new_record
			SeqIO.write(new_record, f, "fasta")
	print sequences
print "done!"






