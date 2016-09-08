# scripts_dump
Where all the little script are stored


- change_filename_with_list
```
	rename the file from a folder from a new list of file name. 

	change_filename_list.py
	test/
	
```	
- Build loci specific matrices from several multi-fasta files that belong to different individuals. 
 ```
	For several individuals, you map the targeted loci, you get a multifasta file per individual with the targeted loci. The script will make matrices per targeted loci.
	From a multifasta per individual to a multifasta per targeted loci. 

	The script is taking from several multifasta files (individual) the first sequence then second... in each mutltifasta file and write it in a new fasta file, name the new fasta file with the name of the the locus.

	Then perform an alignment with mafft.


```