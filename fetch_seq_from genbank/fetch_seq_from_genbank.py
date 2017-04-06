#!/usr/bin/env python2

##Dependencies: Biopython
# to install biopython copy paste in your terminal: sudo pip install biopython

## Output:
# individual file per genera or family or species from the input list
# if your asking thousands of sequences to genbank, you might concider to divide your list into smaller ones.

##Usage
# python fetch_seq_from_genbank.py
# note: need to have the list file in the same folder.

##list file example: one query per line
#Aulosepalum
#Beloglottis
#Brachystele
#Coccineorchis
#Cotylolabium
#Cyclopogon
#Deiregyne
#Dichromanthus
#Dithyridanthus

##Fasta output
# if you want some fasta as output file:
# CHANGE line 47: f = open(genus+".gb", 'w') ###change ".gb" for ".fas"
# CHANGE line 52: SeqIO.write(record, f, "gb") ###change "gb" for "fasta"

from Bio import Entrez, SeqIO

Entrez.email = "vincent.manzanilla@nhm.uio.no" ###CHANGE email address

#handle = Entrez.einfo()
#record = Entrez.read(handle)
#record["DbList"]
list_genus = open("genus_list.txt") ###CHANGE list file if necessary
print "genus\t number_of_seq \t sequence_ID \t length"
for genus in list_genus:
    genus = genus.split("\n")[0]
    handle = Entrez.esearch(db="nucleotide",term=genus+"[Orgn] AND rbcL[Gene]") ### change gene name if necessary
    record = Entrez.read(handle)
    ID_list = record["IdList"]

    if len(ID_list) >=1:
        print(genus + ": \t" + str(len(ID_list)))
        f = open(genus+".gb", 'w') ###change ".gb" for ".fas"  
        for accessions in ID_list:
            handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            print ("\t\t\t" + record.name)+"\t" + str(len(record))
            SeqIO.write(record, f, "gb") ###change "gb" for "fasta"
            handle.close()
    else:
        print(genus + ": \t" + str(len(ID_list)))
        continue