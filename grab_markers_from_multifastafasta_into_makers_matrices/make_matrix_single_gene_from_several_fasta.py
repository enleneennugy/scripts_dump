#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import subprocess
import time

## the script is taking from several multifasta files the first sequence in each mutltifasta file and write it in a new fasta file.
##then perform an alignment with mafft.

#For several individuals, you map the targeted loci, you get a multifasta file per individual with all the targeted sequences. The script will make matrices per target sequence.
##From multifasta per individual to multifasta per target sequences. 


fasta_file = glob.glob('*.fasta')

print(fasta_file)

list_id = []
for seq in SeqIO.parse(fasta_file[0], 'fasta'):
	record = seq.id.split('_')[0:4]
	record = '_'.join(record)
	list_id.append(record)

for name in list_id:
	out = open(name + ".fa", "w")
out.close()

for files in fasta_file:
	fileid = (files[:2])
	for seq1 in SeqIO.parse(files, 'fasta'):
		dna_seq = str(seq1.seq)
		newname1 = seq1.id.split('_')[0:4]
		newname2 = '_'.join(newname1)
		newname3 = fileid + '_' + str(newname2)
		new_record = SeqRecord(Seq(dna_seq), id=newname3, description='')
		with open(newname2 + '.fa', 'a') as fileout:
			SeqIO.write(new_record, fileout, "fasta")
		fileout.close()

print '\n\nMatrices done!!\n\n'




print 'Start alignment with mafft\n\n'
print ((time.strftime("%Y_%m_%d")) + '__' + (time.strftime("%H.%M.%S")))


matrix_fasta = glob.glob('*.fa')

for file in matrix_fasta:
	print('\n\n Star alignment for ' + file + '\n\n')
	#cmd = 'mafft --thread 8 --maxiterate 100 --genafpair --quiet ' + file + '>' + file + '_mafft.fasta'
	#subprocess.check_output(cmd, shell=True)

print 'Done! :)'
