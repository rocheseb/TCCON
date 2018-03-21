#!/var/lib/py27_sroche/bin/python
 # -*- coding: ascii -*-

from __future__ import print_function # allows the use of Python 3.x print function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# code description #
####################

'''
Functions to read GFIT files

Since some dictionnaries are nested I also included somes useful functions to inspect dictionaries.

Note: I still need to also store the units for most of them

'''

####################
# import libraries #
####################

import os

import sys

import subprocess

import numpy as np

from numpy.linalg import inv

import netCDF4

import collections

from datetime import datetime,timedelta

#############
# Functions #
#############

def list_to_matrix(mylist,n):
	'''
	function to convert a 1d array to a 2d array with n columns
	'''
	return [mylist[i:i+n] for i in xrange(0,len(mylist),n)]

# output could be very lengthy !
def show_keys(dic,indent=0):
	'''
	function that displays the keys of a nested dictionary using indents for different levels

	key 1
		sub_key 1
		sub_key 2
			sub_sub_key 1
	key 2
		sub_key 1
		sub_key 2
	etc.
	'''
	for key in sorted(dic.keys()):
		print('\t'*indent,key)
		if type(dic[key]) in [dict,collections.OrderedDict]:
			show_keys(dic[key],indent+1) #recursive function ! calling itself.

def descend_strings(obj):
	'''
	generator function to get all the strings in an iterable object
	This will not work in python 3 because string objects have an __iter__ attribute in 3.x, but not in 2.x
	'''
	if hasattr(obj,'__iter__'):
		if type(obj) in [dict,collections.OrderedDict]:
			for key in obj:
				for result in descend_strings(obj[key]):
					yield result
		else:
			for item in obj:
				for result in descend_strings(item):
					yield result
	if type(obj)==str:
		yield obj

def descend_keys(dic):
	'''
	gets all the keys in a dictionnary and put them in a 1d list
	'''
	for key in sorted(dic.keys()):
		yield key
		if type(dic[key]) in [dict,collections.OrderedDict]:
			for x in descend_keys(dic[key]):
				yield x

def descend_values(dic):
	'''
	gets all the values in a dictionnary and put them in a 1d list
	'''
	for key in sorted(dic.keys()):
		if type(dic[key]) in [dict,collections.OrderedDict]:
			for x in descend_values(dic[key]):
				yield x
		else:
			yield dic[key]

def flatten(obj,keep=[]):
	'''
	gets all the items in a nested iterable oject and put them in a 1d list
	In python 3.x, you need to include 'str' in the keep list if you do not want to explode all strings in individual letters
	'''
	if (type(obj) in [dict,collections.OrderedDict]) and ((dict not in keep) and (collections.OrderedDict not in keep)):
		for x in descend_values(obj):
			yield x
	else:	
		for item in obj:
			if hasattr(item,'__iter__') and (type(item) not in keep):
				for x in flatten(item,keep=keep):
					yield x
			else:
				yield item

def read_opus_header(igram_file_path):
	'''
	use the perl progran OpusHdr to read the acquisition parameters and return them in a dictionary

	igram_file_path: full path to Opus file
	'''

	GGGPATH = os.environ['GGGPATH']

	cwd = os.path.join(GGGPATH,'i2s','scripts') # path

	proc = subprocess.Popen(['./OpusHdr',igram_file_path], stdout=subprocess.PIPE,shell=False,cwd=os.path.join(GGGPATH,'i2s','scripts'))
	tmp = proc.stdout.read() # shell output of OpusHdr
	prompt = tmp.split('\n')

	data = {}

	for line in prompt:
		if (':' not in line) and (line!=''):
			key = line.strip()
			data[key] = {}
		elif line!='':
			try:
				subkey,val = [x.strip() for x in line.split(': ')]
			except ValueError:
				pass
			else:
				data[key][subkey] = val

	return data

def read_col(path):
	'''
	.col files contain the output of GFIT scaling retrievals
	'''

	DATA = {}

	infile = open(path,'r')
	content = infile.readlines()
	infile.close()

	for line in content:
		if 'Spectrum' in line:
			header = line
			break

	start = content.index(header) + 1

	header = header.split()

	content_T = np.array([line.split() for line in content[start:]]).T

	for var in header:
		try:
			DATA[var] = np.array([float(elem) for elem in content_T[header.index(var)]])
		except ValueError:
			DATA[var] = content_T[header.index(var)]

	return DATA

def read_mav(path):
	'''
	.mav files contain a priori information
	'''

	DATA = {}

	infile = open(path,'r')
	content = infile.readlines()
	infile.close()

	content_T = np.array([[float(elem) for elem in line.split()] for line in content[9:]]).T

	header=content[6].split()

	for var in header:
		DATA[var] = content_T[header.index(var)]

	DATA['tropalt'] = float(content[3].split(':')[1])

	return DATA

def read_map(path):
	'''
	.map files are more condensed a priori files
	'''

	DATA = {}

	if os.path.exists(path) == False:
		print('You need to run write_aux to create .map files.')
		return

	infile = open(path,'r')
	content = infile.readlines()
	infile.close()

	for line in content:
		if 'Height' in line:
			header = [elem.split()[0] for elem in line.split(',')]
			start = content.index(line)+2
			units = [elem.split()[0] for elem in content[content.index(line)+1].split(',')]
			break

	content_T = np.array([[float(elem) for elem in line.split(',')] for line in content[start:]]).T

	DATA['units'] = {}

	for var in header:
		DATA[var] = content_T[header.index(var)]
		DATA['units'][var] = units[header.index(var)]

	for i in range(5,start-2):
		var = content[i].split()[0]
		DATA[var] = content[i].split()[2]
		DATA['units'][var] = content[i].split()[1]

	return DATA

def read_spt(path):
	'''
	spt files are the spectrum files output by GFIT/GFIT2
	'''

	DATA = {}

	infile=open(path,'r')
	content=infile.readlines()
	infile.close()

	head = content[2].split()

	DATA['header'] = head

	content_T = np.array([[elem for elem in line.split()] for line in content[3:]],dtype=np.float).T

	DATA['columns'] = {}
	for var in head:
		DATA['columns'][var] = content_T[head.index(var)]

	DATA['sza'] = float(content[1].split()[3])

	resid = 100.0*(DATA['columns']['Tc']-DATA['columns']['Tm']) #the % residuals, tm and tc are transmittances so we just need to multiply by 100 to get %
	rms_resid = np.sqrt(np.mean(np.square(resid)))  #rms of residuals

	DATA['resid'] = resid
	DATA['rms_resid'] = rms_resid

	DATA['params'] = content[1].split()

	return DATA

def read_mod(path):
	'''
	mod files are interpolated NCEP profiles produced by mod_maker
	'''

	DATA = {}

	infile=open(path,'r')
	content=infile.readlines()
	infile.close()

	info = np.array(content[1].split(),dtype=np.float)

	DATA['radius'] = info[0]
	DATA['ecc2'] = info[1]
	DATA['lat'] = info[2]
	DATA['gravity'] = info[3]
	DATA['hold'] = info[4]
	DATA['pfact'] = info[5]

	header = content[3].split()
	units = content[2].split()

	content_T = np.array([[elem for elem in line.split()] for line in content[4:]],dtype=np.float).T

	for var in header:
		DATA[var] = content_T[header.index(var)]

	DATA['units'] = [header[i]+': '+units[i] for i in range(len(header))]

	return DATA

def read_vmr(path):
	'''
	vmr files are used to construct the a priori
	'''

	DATA = {}

	infile=open(path,'r')
	content=infile.readlines()
	infile.close()

	DATA['ztrop'] = content[1].split()[1]
	DATA['date'] = content[2].split()[1]
	DATA['lat'] = content[3].split()[1]

	header = content[4].split()

	content_T = np.array([[elem for elem in line.split()] for line in content[5:]],dtype=np.float).T

	for var in header:
		DATA[var] = content_T[header.index(var)]

	return DATA


def read_tccon(path,mode='eof',variables=[],key_variables=[],flag='all'):
	'''
	 can read tccon .eof, .eof.csv, or .nc files with specified variables of two types:
	 "variables" is a list of exact variables names to be read
	 "key_variables" is a list of keywords, every variable that includes one of the keywords in their name will be read
	 if those lists are not specified all variables (>1200) will be read
	 flag is "all" by default and data with any flag will be read, you can give the appropriate integer if you want to only read data with a specific flag
	 the mode can be set to 'eof', 'csv', or 'netcdf'
	'''

	DATA = {}

	if mode.lower() not in ['eof','netcdf']:
		print('in function read_tccon() : The file format was not recognized')
		return DATA

	if mode.lower() == 'eof':
		infile = open(path,'r')
		content = infile.readlines()
		infile.close()

		if 'csv' in path:
			header = content[2].split(',')
		else:
			header = content[2].split()

		if variables == []:
			variables = header

		# add all the variables that includes keys in key_variables, and make sure there are no repeat variables
		add_var = []
		if key_variables != []:
			for key in key_variables:
				add_var += [var for var in header if key in var]
			for var in add_var:
				if var not in variables:
					variables.append(var)

		if flag == 'all':
			if 'csv' in path:
				content_T = np.array([line.split(',') for line in content[3:]]).T
			else:
				content_T = np.array([line.split() for line in content[3:]]).T
		else:
			if 'csv' in path:
				content_T = np.array([line.split(',') for line in content[3:] if line.split(',')[0]==flag]).T
			else:
				content_T = np.array([line.split() for line in content[3:] if line.split()[0]==flag]).T

		for var in variables:
			if var not in header:
				print('in function read_tccon(): wrong variable input:',var,'not in',path)
				return {}				
			try:
				DATA[var] = np.array( content_T[header.index(var)], dtype=np.float )
			except ValueError:
				print('Skipping string variable:',var)

		DATA['xtime'] = [datetime(int(DATA['year'][i]),1,1)+timedelta(days=int(DATA['day'][i])-1,hours=float(DATA['hour'][i])) for i in range(len(DATA[var]))]

		del DATA['year']
		del DATA['day']
		del DATA['hour']
				
	if mode.lower() == 'netcdf':
		f = netCDF4.Dataset(path,'r')

		header = [v for v in f.variables]

		if variables == []:
			variables = header

		# add all the variables that includes keys in key_variables, and make sure there are no repeat variables
		add_var = []
		if key_variables != []:
			for key in key_variables:
				add_var += [var for var in header if key in var]
			for var in add_var:
				if var not in variables:
					variables.append(var)

		for var in variables:
			try:
				if flag == 'all':
					DATA[var] = np.array( f.variables[var][:], dtype=np.float )
				else:
					DATA[var] = np.array( [f.variables[var][ID] for ID,elem in enumerate(list(f.variables['flag'][:])) if int(elem)==int(flag)] ,dtype=np.float ) # this is kind of slow ...
			except ValueError:
				'Skipping string variable:',var

		DATA['xtime'] = [datetime(int(DATA['year'][i]),1,1)+timedelta(days=int(DATA['day'][i])-1,hours=float(DATA['hour'][i])) for i in range(len(DATA['hour']))]

		f.close()

		del DATA['year']
		del DATA['day']
		del DATA['hour']

	return DATA


def read_runlog(path):
	"""
	runlogs are part of the input files of GGG, they list the spectra as well as coincident atmospheric parameters
	"""

	DATA = {}

	infile = open(path,'r')
	content = infile.readlines()
	infile.close()

	for line in content:
		if 'Spectrum' in line:
			header = line
			break

	start = content.index(header) + 1

	header = header.split()

	nice_content = [line for line in content[start:] if len(line.split())==len(header)]

	content_T = np.array([line.split() for line in nice_content]).T

	for var in header:
		try:
			DATA[var] = [float(elem) for elem in content_T[header.index(var)]]
		except ValueError:
			DATA[var] = content_T[header.index(var)]	

	return DATA

def read_isotopologs(path):
	"""
	isotopologs.dat contains information on the isotopologues of each gas.
	They have a tropospheric and stratospheric isotopic fractionation that is used in gsetup to create the mav file from the vmr file
	"""

	DATA = {}

	infile = open(path,'r')
	content = infile.readlines()
	infile.close()

	gas = ''
	for line in content:

		split_line = line.split()

		gas = split_line[2]
		iso = split_line[1]+gas.lower()

		if iso == '0cirrus15cirrus':
			iso = '0cirrus15'

		DATA[iso] = {}

		for i,elem in enumerate(split_line):
			if i>2:
				if elem == '???':
					break

				try:
					int(elem)
					break
				except ValueError:
					pass

		# some special cases with inconsistencies between the gas name in the isotopologs.dat file and in the vmr file
		if iso == '0cirrus15':
			DATA[iso]['gas'] = 'cirrus15'
		elif iso == '0f141b':
			DATA[iso]['gas'] = 'f141b'
		elif 'c6h6' in iso:
			DATA[iso]['gas'] = 'c6h6'
		elif iso == '1air':
			DATA[iso]['gas'] = 'air'
		else:
			DATA[iso]['gas'] = split_line[2]

		DATA[iso]['kgas'] = int(split_line[0])
		DATA[iso]['kiso'] = int(split_line[1])

		try:
			DATA[iso]['id'] = int(split_line[i])
		except ValueError:
			DATA[iso]['id'] = '???'

		DATA[iso]['fia'] = float(split_line[i+1])
		DATA[iso]['delta'] = float(split_line[i+2])
		DATA[iso]['epsilon'] = float(split_line[i+3])

	return DATA

def read_windows(path):
	'''
	read a file in GGGPATH/windows
	'''

	infile = open(path,'r')
	content = infile.readlines()[1:]
	infile.close()

	for i,line in enumerate(content):
		if line[0]==':':
			del content[i]

	return content

def read_levels(path):
	'''
	read a file in GGGPATH/levels
	'''

	infile = open(path,'r')
	content = [float(line.split()[0]) for line in infile.readlines()]
	infile.close()

	return content

def read_ray(path):
	"""
	slant path distances, bend, fov asza, for each spectrum
	"""

	DATA = {}

	infile = open(path,'r')
	content = infile.readlines()[3:]
	infile.close()

	header = content[0].split()

	nice_content = [line for line in content[1:] if len(line.split())==len(header)]

	content_T = np.array([line.split() for line in nice_content]).T

	for var in header:
		try:
			DATA[var] = [float(elem) for elem in content_T[header.index(var)]]
		except ValueError:
			DATA[var] = content_T[header.index(var)]	

	return DATA




