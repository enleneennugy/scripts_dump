#!/usr/bin/env python2

# Download sequences from genbank according a list of genera or families or species.


# - Dependencies: Biopython
#     To install biopython copy paste in your terminal:
#     ```
#      sudo pip install biopython
#     ```

# - List file example: one query per line
#     ```
#     Aulosepalum
#     Beloglottis
#     Brachystele
#     Coccineorchis
#     Cotylolabium
#     Cyclopogon
#     Deiregyne
#     Dichromanthus
#     Dithyridanthus
#     ```
# - Usage:
#     note: need to have the list file in the same folder.
#     ```
#     python fetch_seq_from_genbank.py
#     ```

# - Output:
#     - individual file per genera or family or species from the input list
#     - if your asking thousands of sequences to genbank, you might concider to divide your list into smaller ones.
#     Example of terminal output:
#     ```
#     genus    number_of_seq   sequence_ID     length
#     Aulosepalum:    1
#                 AJ542433    1316
#     Beloglottis:    1
#                 AJ542432    1314
#     Brachystele:    3
#                 FJ571317    1331
#                 FN870776    1341
#                 JQ045462    699
#     Coccineorchis:  1
#                 AJ542422    1316
#     Cotylolabium:   0
#     Cyclopogon:     3
#                 GQ916978    834
#                 FJ571323    1331
#                 AJ542425    1316
#     Deiregyne:  1
#                 AJ542440    1316
#     Dichromanthus:  2
#                 AJ542439    1316
#                 AJ542438    1314
#     ```

# - Fasta output
#     If you want a fasta formated file as output file:
#     ```
#     CHANGE line 96: f = open(genus+".gb", 'w') ###change ".gb" for ".fas"
#     CHANGE line 106: SeqIO.write(record, f, "gb") ###change "gb" for "fasta"
#     ```



from Bio import Entrez, SeqIO
from urllib2 import HTTPError
import time

Entrez.email = "vincent.manzanilla@gmail.com" ###CHANGE email address



#handle = Entrez.einfo()
#record = Entrez.read(handle)
#record["DbList"]
list_genus = open("genus_list.txt") ###CHANGE list file if necessary
print "genus\t number_of_seq \t sequence_ID \t length"
for genus in list_genus:
    genus = genus.split("\n")[0]
    try:
        handle = Entrez.esearch(db="nucleotide", term="("+genus+"[Orgn]) AND (rbcL[Gene])", retmax=500000) ### change gene name if necessary
    except HTTPError:
        time.sleep(20)
        handle = Entrez.esearch(db="nucleotide", term="("+genus+"[Orgn]) AND (rbcL[Gene])", retmax=500000) ### change gene name if necessary
    record = Entrez.read(handle)
    handle.close()
    ID_list = record["IdList"]

    if len(ID_list) >= 1:
        print(genus + ": \t" + str(len(ID_list)))
        f = open(genus+".gb", 'w') ###change ".gb" for ".fas"
        for accessions in ID_list:
            try:
                handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text")
            except HTTPError:
                time.sleep(20)
                handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            dir(record)
            print ("\t\t\t" + record.name)+"\t" + str(len(record))
            SeqIO.write(record, f, "gb") ###change "gb" for "fasta"
            handle.close()
            time.sleep(1)
    else:
        print(genus + ": \t" + str(len(ID_list)))
        continue
list_genus.close()
