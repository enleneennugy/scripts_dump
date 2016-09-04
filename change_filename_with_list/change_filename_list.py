#!/usr/bin/env python2

from glob import glob
import os


oldfilelist = glob('*.fastq.gz')
print oldfilelist
new_name_list = open('newnames.list', 'r')

for newname in new_name_list:
	for oldfile in oldfilelist:
		new = newname[-27:].strip()
		if oldfile.endswith(new):
			os.rename(oldfile, newname.strip())

print  (glob('*.f*q.gz'))
