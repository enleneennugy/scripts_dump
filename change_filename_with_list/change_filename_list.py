#!/usr/bin/env python2

from glob import glob
import os

#change_filename_with_list
#	rename files from a folder with a new list of file name.
#	Add a prefix on each filename. The identification is based on the end of the filename. 


oldfilelist = glob('*.fastq.gz')
print oldfilelist
new_name_list = open('newnames.list', 'r')

for newname in new_name_list:
	for oldfile in oldfilelist:
		old = newname[-27:-12].strip()+'_r.fastq.gz'
		new = newname[-27:-12].strip()+'_f.fastq.gz'
		print old + ' ' + new
		
		if oldfile.endswith(old):
			print oldfile
			os.rename(oldfile, newname.strip())

print (glob('*.f*q.gz'))
