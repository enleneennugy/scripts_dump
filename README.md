# scripts_dump
Where all the small scripts are stored


- change_filename_with_list
```
	rename the files in the fastq.gz format given a list of new names. 

	change_filename_list.py
	test/

	USAGE: run the script in the folder with your fasta file: python  change_filename_list.py

```	
- grab_markers_from_multifastafasta_into_makers_matrices.	
```
	Build loci specific matrices from several multi-fasta files that belong to different individuals. 
 	
	For several individuals, you map the targeted loci, you get a multifasta file per individual with the targeted loci. The script will make matrices per targeted loci.
	From a multifasta per individual to a multifasta per targeted loci. 

	The script is taking from several multifasta files (individual) the first sequence then second... in each mutltifasta file and write it in a new fasta file, name the new fasta file with the name of the the locus.

	Then perform an alignment with mafft.

	make_matrix_single_gene_from_several_fasta.py

	USAGE: run the script in the folder with your fasta file: python change_filename_list.py
```

- compare list.

```
	Compare list, give the common elements of all the lists. 
```

- Extract genes from GenBank file format given a list of genes, make new fasta files per genes and align. 

```
    Reads in the requested CDS name, searches for the CDS in the genbank file,
    returns the DNA sequence, with organism name, gene name in the
    sequence description. Create individual fasta file for each gene.


	USAGE: ./extract_genes_from_genebank_files.py -l gene.list -f ~/Downloads/

	Dependencies: pasta aligner with in the PATH, run_pasta.py.



	example folder: genes.list and genebank records
```




