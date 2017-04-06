Download sequences from genbank according a list of genea or families or species. 


- Dependencies: Biopython
	To install biopython copy paste in your terminal: 
	```
	 sudo pip install biopython
	```

- List file example: one query per line
	```
	Aulosepalum
	Beloglottis
	Brachystele
	Coccineorchis
	Cotylolabium
	Cyclopogon
	Deiregyne
	Dichromanthus
	Dithyridanthus
	```
- Usage:
	note: need to have the list file in the same folder.
	```
	python fetch_seq_from_genbank.py
	```

- Output:
	- individual file per genera or family or species from the input list
	- if your asking thousands of sequences to genbank, you might concider to divide your list into smaller ones.
	Example of terminal output: 
	```
	genus	 number_of_seq 	 sequence_ID 	 length
	Aulosepalum: 	1
				AJ542433	1316
	Beloglottis: 	1
				AJ542432	1314
	Brachystele: 	3
				FJ571317	1331
				FN870776	1341
				JQ045462	699
	Coccineorchis: 	1
				AJ542422	1316
	Cotylolabium: 	0
	Cyclopogon: 	3
				GQ916978	834
				FJ571323	1331
				AJ542425	1316
	Deiregyne: 	1
				AJ542440	1316
	Dichromanthus: 	2
				AJ542439	1316
				AJ542438	1314
	```

- Fasta output
	If you want a fasta formated file as output file: 
	```
	CHANGE line 47: f = open(genus+".gb", 'w') ###change ".gb" for ".fas"
	CHANGE line 52: SeqIO.write(record, f, "gb") ###change "gb" for "fasta"
	```

