"""
Code to convert an existing runlog to a "synthetic" runlog that can be used to generate synthetic spectra with GGG


Dependencies:
	- fortranformat=0.2.5: https://pypi.org/project/fortranformat/0.2.5/
	- parse=0.12.0: https://pypi.org/project/parse/1.12.0/

For usage, run:

python write_syn_runlog.py -h
"""


import os
import sys
import parse
import argparse
import fortranformat as ff


def write_syn_runlog(path,v0=4750,v1=8250,DELTA_NU=0.0111111111):
	'''
	Input:
		- path : full path to an existing runlog
		- v0 : starting wavenumber (cm-1)
		- v1 : last wavenumber (cm-1)
		- delta_nu : wavenumber spacing (cm-1)

	This will write a new runlog that can be used to generate synthetic spectra.
	The new runlog will be saved in the same place but the filename will finish with "_syn.grl"
	'''

	with open(path,'r') as infile:
		content = infile.readlines()

	syn_dict = {
		'BPW':7,
		'POINTER':0,
		'APF':'N1',
		'DELTA_NU':DELTA_NU,
		'IFIRST':int(round(v0/DELTA_NU)),
		'ILAST':int(round(v1/DELTA_NU)),
		'SNR':1000,
	}

	header = [' ']+content[3].split()

	ID_dict = {key:header.index(key) for key in syn_dict}

	fmt = list(parse.parse('format={}\n',content[2]))[0]

	reader = ff.FortranRecordReader(fmt)
	writer = ff.FortranRecordWriter(fmt)

	for i,line in enumerate(content):
		if i>=4:
			data = reader.read(line)
			for key in syn_dict:
				data[ID_dict[key]] = syn_dict[key]
			new_line = writer.write(data)+'\n'
			content[i] = new_line

	new_path = path.replace('.grl','_syn.grl')
	with open(new_path,'w') as outfile:
		outfile.writelines(content)


if __name__=="__main__":

	parser = argparse.ArgumentParser(description='Code to convert an existing runlog (.grl) to a "synthetic" runlog that can be used to generate synthetic spectra with GGG')

	parser.add_argument('path',help='full path to an existing runlog')
	parser.add_argument('--v0',type=int,default=4750,help='first wavenumber, default=4750')
	parser.add_argument('--v1',type=int,default=8250,help='last wavenumber, default=8250')
	parser.add_argument('--dnu',type=float,default=0.0111111111,help='wavenumber spacing, default=0.0111111111')

	args = parser.parse_args()

	if not os.path.exists(args.path):
		sys.exit('Invalid path: {}'.format(args.path))

	write_syn_runlog(args.path,v0=args.v0,v1=args.v1,DELTA_NU=args.dnu)