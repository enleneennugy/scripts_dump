#!/usr/bin/env python2

from glob import glob
import os

#change_filename_with_list
#	rename the file from a folder from a new list of file name. 


oldfilelist = glob('*.fastq.gz')
print oldfilelist
new_name_list = open('newnames.list', 'r')

for newname in new_name_list:
	for oldfile in oldfilelist:
		new = newname[-27:].strip()
		if oldfile.endswith(new):
			os.rename(oldfile, newname.strip())

print  (glob('*.f*q.gz'))
