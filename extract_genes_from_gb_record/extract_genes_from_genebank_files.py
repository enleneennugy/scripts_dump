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

#    '''
#    Reads in the requested CDS name, searches for the CDS in the genbank file,
#    returns the DNA sequence, with organism name, gene name in the
#    sequence description. Create individual fasta file for each gene.
#    '''

#make a new directory for the analyses, day/time
tm = (time.strftime("%H.%M.%S"))
dt = (time.strftime("%Y_%m_%d"))
now = (dt + "_" + tm)
print 'Starting to extract the genes from the multiple genbank files.\nFasta files writen in the folder: '+(now) + '/01_fasta_file'
os.mkdir(now)

#directory where the new fasta file will be store
os.mkdir(now+'/01_fasta_file')
os.mkdir(now+'/02_fasta_file_warning')


#Get the file list and the list of genes to extract from the genbank file. 
filelist = glob.glob('/Users/vincem/Research/PhD/Research topics/panax/gb_files/*.gb')
list_genes = open('gene.list', 'r').read().splitlines()

#change dir and create the log file and the warning file
os.chdir(now)
w = open(('02_fasta_file_warning/genbank_warning.log'), 'w')
print 'Warning writen in the folder 02_fasta_file_warning/genbank_warning.log'
l = open('01_fasta_file/genbank.log', 'w')
print 'Log file writen in the folder 01_fasta_file/genbank.log'
found = str(len(filelist)) + ' files found.\n'
l.write(found)

print '****warnings****'

file_list = []
for genes in list_genes:
    gene_append = []
    #print genes
    for file_gb in filelist:
        for gb_record in SeqIO.parse(open(file_gb, 'r'), 'genbank'):
            n_hits = 0
            id_name = gb_record.id
            organism_name = gb_record.annotations['organism']
            #print id_name + ' ' + organism_name
            for index, feature in enumerate(gb_record.features):
                if feature.type == 'CDS':
                    gb_feature = gb_record.features[index]
                    if 'gene' in gb_feature.qualifiers:
                        if genes == gb_feature.qualifiers['gene'][0]:
                            n_hits = n_hits + 1
                            
                            if n_hits <= 1:
                                product = gb_feature.qualifiers['gene'][0]
                                DNAseq = gb_feature.extract(gb_record.seq)
                                record_new = SeqRecord(Seq(str(DNAseq), generic_dna), id=id_name, name=product, description=product + ' ' + organism_name)
                                record_new.annotations['organism'] = organism_name 
#                                print (record_new.format("genbank"))
                                gene_append.append(record_new)
                                SeqIO.write(gene_append, '01_fasta_file/' + genes + '.fas', "fasta")
#                                print (record_new.format("fasta"))                               
                            if n_hits > 1:
                                DNAseq = gb_feature.extract(gb_record.seq)
                                if record_new.seq == DNAseq:
                                    continue
                                elif record_new.seq != DNAseq:
                                    #warning printed and log file: genbank_warning.log
                                    warning = 'Warning, ' + str(n_hits) + ' found in ' + id_name + ' for ' + genes 
                                    print (warning)
                                    w.write(warning)
                                    product = gb_feature.qualifiers['gene']
                                    DNAseq = gb_feature.extract(gb_record.seq)
                                    record_new = SeqRecord(Seq(str(DNAseq), generic_dna), id=id_name, name=genes, description=product[0])
                                    gene_append.append(record_new)
                                    SeqIO.write(gene_append, '02_fasta_file_warning/' + genes + '.fas', "fasta")
    # print in a log file the genes recovered + number of files
    log = str(len(gene_append)) + ' sequences recovered for the gene ' + genes + '\n'
    #print log
    l.write(log)
w.close()
l.close()


print '\n****alignment started****'
#make alignments
#   dependencies:
#        - pasta aligner with an alias: python /usr/local/bin/pasta/pasta/pasta/run_pasta.py

os.mkdir('03_alignments')

filea_genes = glob.glob('01_fasta_file/*.fas')
for filea in filea_genes:
    filea_split = filea.split('01_fasta_file/')[1]
    print filea_split
    cmd = 'python /usr/local/bin/pasta/pasta/run_pasta.py  --auto -j ' + filea_split + ' --input=./' + filea
    subprocess.check_output(cmd, shell=True)



