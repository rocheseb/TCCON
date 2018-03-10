#!/usr/bin/env python2.7
 # -*- coding: ascii -*-

from __future__ import print_function # allows the use of Python 3.x print function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

'''
This program can be used the merge saved sessions of lft_app.
The resulting file can be loaded from the app.

A list of filenames from the lft_app/saved_session/ should be written in lft_app/utils/to_merge.dat

e.g.

saved_session_1.npy
saved_session_2.npy
etc. 

If two sessions contain tests with the same name, the data from the last file in the list will overwrite that test's data.

e.g. above if both saved_session_1.npy and saved_session_2.npy contain a test "HCl_45_eu_150916 reg=1.8", then the data in saved_session_2.npy will be the one in the merged session.

run with:
python merge_sessions.py arg1

arg1: name of the merged file (with no extensions)
'''

import os
import sys
import numpy as np

util_path = os.getcwd() # ... /linefit/lft145/lft_app/utils
app_path = (os.sep).join(util_path.split(os.sep)[:-1]) # the app should be in ... /linefit/lft145/lft_app
save_path = os.path.join(app_path,'saved_sessions') # ... /linefit/lft145/lft_app/saved_sessions

infile = open('to_merge.dat','r')
content = infile.readlines()
infile.close()

argu = sys.argv
try:
	save_name = argu[1]+'.npy'
except IndexError:
	print('Need one argument: name of the merged file (with no extensions)')
	sys.exit()

merged_session = {}
for line in content:
	file_path = os.path.join(save_path,line.split('\n')[0])
	session = np.load(file_path).item()

	merged_session.update(session)

np.save(os.path.join(save_path,save_name),merged_session)