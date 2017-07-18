#!/usr/bin/env python2

listname = open("/Users/vincem/name.list")

for name in listname:
	
	name = name.replace('\n', '')
	print ("\n\n"+name)
	csvfile = open("/Users/vincem/Downloads/Specie_Mount.csv")
	for rowfile in csvfile:
		if name in str(rowfile):
			print 1
		else:
			print 0
			
		
	
