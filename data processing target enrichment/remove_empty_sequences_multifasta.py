import sys
from Bio import SeqIO

#Remove empty sequences in a fasta alignment

#usage: python clean_seq.py matrix.fasta

fasta_file = sys.argv[1]



f = open("formated/" + fasta_file, 'w')
sequences = {}
n = 0
for seq_record in SeqIO.parse(fasta_file, "fasta"):

    sequence = str(seq_record.seq)
    a = (len(sequence))
    if ((sequence.count("-")) != a) == True:
        n=n+1
        sequences[sequence] = seq_record.id
        SeqIO.write(seq_record, f, "fasta") ###change "gb" for "fasta"

    else:
        print seq_record.id
        print (len(sequence))
        print ((sequence.count("-")))

# Write the clean sequences

# Create a file in the same directory where you ran this script

print("\n\nCLEAN!!!\n\nPlease check clean_" + fasta_file)
print "Number of remaining sequences: " + str(n)