#!~/anaconda2/bin/python
 # -*- coding: utf-8 -*-

import os
from datetime import datetime, timedelta
import sys
import urllib
import requests
import shutil

####################
# Code Description #
####################
"""
Functions to create list of URLs and/or download them
"""

#############
# Functions #
#############

def URLlist_FP(start,end,timestep=timedelta(hours=3),outpath='',surf=False):
	"""
	GEOS5-FP data has one global file every 3 hours (from 00:00 to 21:00 UTC each day)
	start: datetime object, start of the desired date range
	end: datetime object, end of the desired date range
	timestep: use the model time resolution to get all files, or a multiple of it to get less files
	outpath: full path to the file in which the list of URLs will be written
	"""
	if surf:
		fmt = "https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/das/Y{}/M{:0>2}/D{:0>2}/GEOS.fp.asm.inst3_2d_asm_Nx.{}_{:0>2}00.V01.nc4\n"
	else:
		fmt = "https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/das/Y{}/M{:0>2}/D{:0>2}/GEOS.fp.asm.inst3_3d_asm_Np.{}_{:0>2}00.V01.nc4\n"
	
	if outpath=='': # if no specified full path to make the file, just write a file in the current directory 
		outpath = 'getFP.dat'

	print 'Writting URL list in:',outpath

	curdate = start
	with open(outpath,'w') as f:
		while curdate<end:
			f.write(fmt.format(curdate.year,curdate.month,curdate.day,datetime.strftime(curdate,'%Y%m%d'),curdate.hour))
			curdate += timestep

	return outpath

def URLlist_FPIT(start,end,timestep=timedelta(hours=3),outpath='',surf=False):
	"""
	GEOS5-FP-IT data has one global file every 3 hours (from 00:00 to 21:00 UTC each day)
	start: datetime object, start of the desired date range
	end: datetime object, end of the desired date range
	timestep: use the model time resolution to get all files, or a multiple of it to get less files
	outpath: full path to the file in which the list of URLs will be written
	"""
	if surf:
		fmt = "http://goldsfs1.gesdisc.eosdis.nasa.gov/data/GEOS5/DFPITI3NXASM.5.12.4/{}/{:0>3}/.hidden/GEOS.fpit.asm.inst3_2d_asm_Nx.GEOS5124.{}_{:0>2}00.V01.nc4\n"
	else:
		fmt = "http://goldsfs1.gesdisc.eosdis.nasa.gov/data/GEOS5/DFPITI3NPASM.5.12.4/{}/{:0>3}/.hidden/GEOS.fpit.asm.inst3_3d_asm_Np.GEOS5124.{}_{:0>2}00.V01.nc4\n"
	
	if outpath=='': # if no specified full path to make the file, just write a file in the current directory 
		outpath = 'getFPIT.dat'

	print 'Writting URL list in:',outpath

	curdate = start
	with open(outpath,'w') as f:
		while curdate<end:
			f.write(fmt.format(curdate.year,curdate.timetuple().tm_yday,datetime.strftime(curdate,'%Y%m%d'),curdate.hour))
			curdate += timestep

	return outpath

def URLlist_MERRA2(start,end,timestep=timedelta(days=1),outpath='',surf=False):
	"""
	MERRA-2 data has one global file every 1 days (from 00:00 to 21:00 UTC each day)
	start: datetime object, start of the desired date range
	end: datetime object, end of the desired date range
	timestep: use the model time resolution to get all files, or a multiple of it to get less files
	outpath: full path to the file in which the list of URLs will be written
	"""

	if surf:
		fmt = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/{}/{:0>2}/MERRA2_400.inst1_2d_asm_Nx.{}.nc4\n"
	else:
		fmt = "https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4/{}/{:0>2}/MERRA2_400.inst3_3d_asm_Np.{}.nc4\n"

	if outpath=='': # if no specified full path to make the file, just write a file in the current directory 
		outpath = 'getMERRA2.dat'
	
	print 'Writting URL list in:',outpath

	curdate = start
	with open(outpath,'w') as f:
		while curdate<end:
			YYYYMMDD = datetime.strftime(curdate,'%Y%m%d')
			f.write(fmt.format(curdate.year,curdate.month,YYYYMMDD))
			curdate += timestep

	return outpath

def download_file(url,outpath=''):
	"""
	url: url pointing to a file to download
	outpath: path to the directory in which the file will be saved

	Note: make a .netrc file in your home directory for requests to automatically handle authentifications
	"""
	filename = url.split('/')[-1]
	if outpath != '':
		outpath = os.path.join(outpath,filename)
	else:
		outpath = filename
	r = requests.get(url, stream=True)
	with open(outpath, 'wb') as f:
		shutil.copyfileobj(r.raw, f)
	print url,'downloaded to',outpath

def get_data(URL_list,outpath=''):
	"""
	URL_list: full path to a file contain a list of URLs to download files

	This will not work with MERRA2 lists unless you have your earthdata credentials in a .netrc file in your home directory (~/.netrc)

	For FP-IT it will only work on computers which have been given access to the data
	"""
	with open(URL_list) as f:
		content = f.read().splitlines()

	for URL in content:
		fail = 0
		while fail < 3:
			try:
				download_file(URL,outpath)
			except:
				fail += 1
				print "Fail",fail,"with",URL

########
# Main #
########

if __name__=="__main__":
	argu = sys.argv

	func_dict = {'FP':URLlist_FP,'FPIT':URLlist_FPIT,'MERRA':URLlist_MERRA2}

	mode = argu[1]
	start = datetime.strptime(argu[2],'%Y%m%d') # YYYYMMDD
	end = datetime.strptime(argu[3],'%Y%m%d') # YYYYMMDD

	surf = False
	if 'surf' in argu:
		surf = True

	URL_list = func_dict[mode](start,end,surf=surf)

	download = raw_input('Download? (Y/N)') in ['Y','y']

	if download:
		get_data(URL_list)



