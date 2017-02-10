#!/usr/bin/env python3

from glob import glob
import os
from pandas import DataFrame, read_csv
import pandas as pd
import numpy as np
from collections import Counter
import itertools

'Compare file list and return the common rows, based on csv files.'

filelist = glob('*.list')
l = len(filelist)

d = []
for file in filelist:
	print('\n' + file)
	df = pd.read_csv(file, index_col =False, sep='\t', header=(0), )
	a = list(df.loc[(df['Reads']>= 30) & (df['Mismatch %']>0) & (df['Mismatch %']< 5), 'Contig'])
	d.append(a)


d = list(itertools.chain.from_iterable(d))

dict_list = (Counter(d))
for key, value in dict_list.items():
    if value == l:
        print key

	