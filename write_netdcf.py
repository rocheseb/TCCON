from __future__ import print_function
"""
Use the official output files (.oof) to write the actually official output files (.nc)
"""

import os
import sys
import platform
import subprocess
import netCDF4
import pandas as pd
import numpy as np
import hashlib
import argparse
from collections import OrderedDict
import time
from datetime import datetime, timedelta

def progress(i,tot,bar_length=20,word=''):
	"""
	a fancy loadbar to be displayed in the prompt while executing a time consuming loop
	"""
	if tot==0:
		tot=1
	percent=float(i+1)/tot
	hashes='#' * int(round(percent*bar_length))
	spaces=' ' * (bar_length - len(hashes))
	sys.stdout.write("\rPercent:[{0}] {1}%".format(hashes + spaces, int(round(percent * 100)))+"    "+str(i+1)+"/"+str(tot)+'  Spectrum: '+word)
	sys.stdout.flush()

def md5(file_name):
	"""
	Reads file_name and get its md5 sum, returns the hexdigest hash string

	file_name: full path to the file
	"""
	hash_md5 = hashlib.md5()
	with open(file_name, "rb") as f:
		for chunk in iter(lambda: f.read(4096), b""):
			hash_md5.update(chunk)
	return hash_md5.hexdigest()

def checksum(file_name,hexdigest):
	"""
	Compare the hash string hexdigest with the md5 sum of the file file_name

	hexdigest: hash string
	file_name: full path to the file
	"""

	return md5(file_name) == hexdigest

def file_info(file_name):
	"""
	Read the first line of a file and get the number of header lines and number of data columns

	file_name: full path to the file
	"""
	with open(file_name,'r') as infile:
		nhead,ncol = [int(i) for i in infile.readline().strip().split()[:2]]
	nhead = nhead-1

	return nhead,ncol

if __name__=='__main__': # execute only when the code is run by itself, and not when it is imported

	try:
		GGGPATH = os.environ['GGGPATH']
	except:
		try:
			GGGPATH = os.environ['gggpath']
		except:
			print('You need to set a GGGPATH (or gggpath) environment variable')
			sys.exit()

	wnc_version = ' WRITE_NC             Version 1.0     2019-10-04     SR '
	print(wnc_version)

	description = wnc_version + "\nThis writes TCCON EOF files in NETCDF"
	
	parser = argparse.ArgumentParser(description=description,formatter_class=argparse.RawTextHelpFormatter)
	
	def file_choices(choices,file_name):
		"""
		Function handler to check file extensions with argparse

		choices: tuple of accepted file extensions
		file_name: path to the file
		"""
		ext = os.path.splitext(file_name)[1][1:]
		if ext not in choices:
			parser.error("file doesn't end with one of {}".format(choices))
		return file_name
	
	parser.add_argument('tav-file',type=lambda file_name:file_choices(('tav'),file_name),help='The .tav file')

	args = parser.parse_args()

	# input and output file names
	tav_file = args.tav_file
	vav_file = tav_file.replace('.tav','.vav')
	asw_file = tav_file.replace('.tav','.asw')
	ada_file = vav_file+'.ada'
	aia_file = ada_file+'.aia'
	esf_file = aia_file+'.daily_error.out'
	oof_file = aia_file+'.oof'
	
	siteID = tav_file.split(os.sep)[-1][:2] # two letter site abbreviation
	qc_file = os.path.join(GGGPATH,'tccon','{}_qc.dat'.format(siteID))
	header_file = os.path.join(GGGPATH,'tccon','{}_oof_header.dat'.format(siteID))
	correction_file =  os.path.join(GGGPATH,'tccon','corrections.dat')
	lse_file = os.path.join(GGGPATH,'lse','gnd',tav_file.split(os.sep)[-1].replace('.tav','.lse'))
	nc_file = tav_file.replace('.tav','.nc') # the final output file

	col_file_list = sorted([i for i in os.listdir(os.getcwd()) if '.col' in i])
	map_file_list = sorted([i for i in os.listdir(os.getcwd()) if '.map' in i])

	if not col_file_list: # [] evaluates to False
		print('No .col files !')
		sys.exit()
	if not map_file_list:
		print('No .map files !')
		sys.exit()

	## read data, I add the file_name to the data dictionaries for some of them

	# multiggg.sh
	with open('multiggg.sh','r') as infile:
		content = [line for line in infile.readlines() if line[0]!=':' or line.strip()!=''] # the the file without blank lines or commented out lines starting with ':'
	ncol = len(content)
	if ncol!=len(col_file_list):
		print('/!\\ multiggg.sh has {} command lines but there are {} .col files'.format(ncol,len(col_file_list)))
		sys.exit()

	# header file
	with open(header_file,'r') as infile:
		header_content = infile.read()

	# correction file
	nhead, ncol = file_info(correction_file)
	correction_data = pd.read_csv(correction_file,delim_whitespace=True,skiprows=nhead)

	# qc file
	nhead, ncol = file_info(qc_file)
	qc_data = pd.read_fwf(qc_file,widths=[15,3,8,7,10,9,10,45],skiprows=nhead+1,names='Variable Output Scale Format Unit Vmin Vmax Description'.split())
	for key in ['Variable','Format','Unit']:
		qc_data[key] = [i.replace('"','') for i in qc_data[key]]

	# error scale factors
	nhead, ncol = file_info(esf_file)
	esf_data = pd.read_csv(esf_file,delim_whitespace=True,skiprows=nhead)

	# oof file
	nhead, ncol = file_info(oof_file)
	oof_data = pd.read_csv(oof_file,delim_whitespace=True,skiprows=nhead)
	oof_data['file'] = oof_file
	site_info = pd.read_csv(oof_file,delim_whitespace=True,skiprows=lambda x: x in range(nhead-3) or x>=nhead-1) # has keys ['Latitude','Longitude','Altitude','siteID']

	# lse file
	nhead, ncol = file_info(lse_file)
	lse_data = pd.read_csv(lse_file,delim_whitespace=True,skiprows=nhead)
	lse_data['file'] = lse_file
	lse_data.rename(index=str,columns={'Specname':'spectrum'},inplace=True) # the other files use 'spectrum'

	# tav file
	with open(tav_file,'r') as infile:
		nhead,ncol,nrow,naux = np.array(infile.readline().split()).astype(int)
	nhead = nhead-1
	tav_data = pd.read_csv(tav_file,delim_whitespace=True,skiprows=nhead)
	tav_data['file'] = tav_file
	nwin = int((ncol-naux)/2)

	# vav file
	nhead, ncol = file_info(vav_file)
	vav_data = pd.read_csv(vav_file,delim_whitespace=True,skiprows=nhead)
	vav_data['file'] = vav_file

	# ada file
	nhead, ncol = file_info(ada_file)
	ada_data = pd.read_csv(ada_file,delim_whitespace=True,skiprows=nhead)
	ada_data['file'] = ada_file
	
	# aia file
	nhead, ncol = file_info(aia_file)
	aia_data = pd.read_csv(aia_file,delim_whitespace=True,skiprows=nhead)
	aia_data['file'] = aia_file

	## check all files have the same spectrum list as the vav file
	check_spec = np.array([data['spectrum']==vav_data['spectrum'] for data in [tav_data,ada_data,aia_data,oof_data]]).flatten()
	if False in check_spec:
		print('Files have inconsistent spectrum lists !')
		for data in [vav_data,ada_data,aia_data,oof_data]:
			print(len(data['spectrum']),'spectra in',data['file'][0])
		sys.exit()

	# make all the column names consistent between the different files
	for dataframe in [correction_data,qc_data,esf_data,oof_data,lse_data,vav_data,ada_data,aia_data]:
		dataframe.rename(str.lower,axis='columns',inplace=True) # all lower case
		if 'doy' in dataframe.columns: # all use 'day' and not 'doy'
			dataframe.rename(index=str,columns={'doy':'day'},inplace=True)
		if 'lon' in dataframe.columns:
			dataframe.rename(index=str,columns={'lon':'long'},inplace=True)

	qc_data['rsc'] = qc_data['scale'].copy()

	gas_list = [i for i in tav_data.columns[naux:] if ('_error' not in i) and ('file' not in i)] # gas names

	mchar = 0
	if 'spectrum' in tav_data.columns:
		mchar = 1

	#sys.exit()

	# Let's try to be CF compliant: http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.pdf
	standard_name_dict = {
	'year':'year',
	'run':'run_number',
	'lat':'latitude',
	'long':'longitude',
	'hour':'decimal_hour',
	'azim':'solar_azimuth_angle',
	'asza':'astronomical_solar_zenith_angle',
	'day':'day_of_year',
	'wspd':'wind_speed',
	'wdir':'wind_direction',
	'graw':'spectrum_spectral_point_spacing',
	'tins':'instrument_internal_temperature',
	'tout':'atmospheric_temperature',
	'pins':'instrument_internal_pressure',
	'pout':'atmospheric_pressure',
	'hout':'atmospheric_humidity',
	'sia':'solar_intensity_average',
	'fvsi':'fractional_variation_in_solar_intensity',
	'zobs':'observation_altitude',
	'zmin':'pressure_altitude',
	'osds':'observer_sun_doppler_stretch',
	'gfit_version':'gfit_version',
	'gsetup_version':'gsetup_version',
	'fovi':'internal_field_of_view',
	'opd':'maximum_optical_path_difference',
	'rmsocl':'fit_rms_over_continuum_level',
	'nit':'number_of_iterations',
	'cl':'continuum_level',
	'ct':'continuum_tilt',
	'cc':'continuum_curvature',
	'fs':'frequency_shift',
	'sg':'solar_gas_shift',
	'zo':'zero_level_offset',
	'zpres':'pressure_altitude',
	}

	checksum_var_list = ['config','apriori','runlog','levels','mav','ray','isotopologs','windows','telluric_linelists','solar']

	standard_name_dict.update({var+'_checksum':var+'_checksum' for var in checksum_var_list})

	long_name_dict = {key:val.replace('_',' ') for key,val in standard_name_dict.iteritems()} # standard names without underscores

	if os.path.exists(nc_file):
		os.remove(nc_file)

	with netCDF4.Dataset(nc_file,'w',format='NETCDF4_CLASSIC') as eof:
		
		## global attributes
		eof.More_information = "https://tccon-wiki.caltech.edu"
		eof.TCCON_reference = "Wunch, D., G. C. Toon, J.-F. L. Blavier, R. A. Washenfelder, J. Notholt, B. J. Connor, D. W. T. Griffith, V. Sherlock, and P. O. Wennberg (2011), The total carbon column observing network, Philosophical Transactions of the Royal Society - Series A: Mathematical, Physical and Engineering Sciences, 369(1943), 2087-2112, doi:10.1098/rsta.2010.0240. Available from: http://dx.doi.org/10.1098/rsta.2010.0240"
		eof.Data_Use_Policy = "https://tccon-wiki.caltech.edu/Network_Policy/Data_Use_Policy"
		eof.Auxiliary_Data_Description = "https://tccon-wiki.caltech.edu/Network_Policy/Data_Use_Policy/Auxiliary_Data"
		eof.description = '\n'+header_content
		eof.history = "Created {} (UTC)".format(time.asctime(time.gmtime(time.time())))
		eof.source = "Created with Python {} and the library netCDF4 {}".format(platform.python_version(),netCDF4.__version__)
		eof.Flag_info = 'The Vmin and Vmax attributes of some variables indicate the range of values out of which the data would be flagged bad'
		eof.Number_of_species = str(nwin)
		eof.Number_of_spectral_windows = str(len(col_file_list))
		
		proc = subprocess.Popen(['hg','summary'],cwd=GGGPATH,stdout=subprocess.PIPE)
		out, err = proc.communicate()
		eof.GGGtip = "The output of 'hg summary' from the GGG repository:\n"+out

		## groups

		# if we want to get rid of write_aux, we could write the priors in the netcdf file in this code
		#eof.createGroup('apriori')
		#eof['apriori'].description = "A priori profiles, TCCON uses one apriori per day at local noon"
		#eof['apriori'].createDimension('level')
		#eof['apriori'].createDimension('time',None) # unlimited dimension
		#times = eof['apriori'].createVariable('time',np.float64,('time',))
		#times.description = 'fractional days since epoch (UTC)'
		#times.units = 'days since 1970-01-01 00:00:00'
		#times.calendar = 'gregorian'

		## create dimensions
		tim = eof.createDimension('time',None)
		fixstring20 = eof.createDimension('fixstring20',20)
		fixstring32 = eof.createDimension('fixstring32',32)

		## create coordinate variables
		eof.createVariable('time',np.float64,('time',))
		eof['time'].standard_name = "time"
		eof['time'].long_name = "time"
		eof['time'].description = 'fractional days since 1970-01-01 00:00:00 (UTC)'
		eof['time'].units = 'days since 1970-01-01 00:00:00'
		eof['time'].calendar = 'gregorian'

		## create variables

		# checksums
		for var in checksum_var_list:
			checksum_var = eof.createVariable(var+'_checksum','S1',('time','fixstring32'))
			checksum_var.standard_name = standard_name_dict[var+'_checksum']
			checksum_var.long_name = long_name_dict[var+'_checksum']
			checksum_var.description = 'hexdigest hash string of the md5 sum of the {} file'.format(var)

		# code versions
		eof.createVariable('gfit_version',np.float64,('time',))
		eof['gfit_version'].description = "version number of the GFIT code that generated the data"
		eof['gfit_version'].standard_name = standard_name_dict['gfit_version']
		eof['gfit_version'].long_name_dict = long_name_dict['gfit_version']

		eof.createVariable('gsetup_version',np.float64,('time',))
		eof['gsetup_version'].description = "version number of the GSETUP code that generated the priors"
		eof['gsetup_version'].standard_name = standard_name_dict['gsetup_version']
		eof['gsetup_version'].long_name_dict = long_name_dict['gsetup_version']

		eof.createVariable('flag',np.float64,('time',))
		eof['flag'].description = 'data quality flag, 0 = good'
		eof['flag'].standard_name = 'quality_flag'
		eof['flag'].long_name = 'quality flag'

		eof.createVariable('spectrum','S1',('time','fixstring20',))
		eof['spectrum'].standard_name = 'spectrum_file_name'
		eof['spectrum'].description = 'spectrum file name'
		eof['spectrum'][:] = [list(elem) for elem in aia_data['spectrum'].values]

		# auxiliary variables
		aux_var_list = [tav_data.columns[i] for i in range(1,naux)]
		for var in aux_var_list: 
			qc_id = list(qc_data['variable']).index(var)
			digit = int(qc_data['format'][qc_id].split('.')[-1])
			eof.createVariable(var,np.float64,('time',),zlib=True,least_significant_digit=digit)
			if var in standard_name_dict.keys():
				eof[var].standard_name = standard_name_dict[var]
				eof[var].long_name = long_name_dict[var]
			# set attributes using the qc.dat file
			eof[var].description = qc_data['description'][qc_id]
			eof[var].units = qc_data['unit'][qc_id].replace('(','').replace(')','').strip()
			eof[var].vmin = qc_data['vmin'][qc_id]
			eof[var].vmax = qc_data['vmax'][qc_id]

		# averaged variables (from the different windows of each species)
		main_var_list = [tav_data.columns[i] for i in range(naux,len(tav_data.columns)-1)]  # minus 1 because I added the 'file' column
		for var in main_var_list:
			xvar = 'x'+var
			qc_id = list(qc_data['variable']).index(xvar)

			digit = int(qc_data['format'][qc_id].split('.')[-1])
			eof.createVariable(xvar,np.float64,('time',),zlib=True,least_significant_digit=digit)
			eof[xvar].standard_name = xvar
			eof[xvar].long_name = xvar.replace('_',' ')
			eof[xvar].description = qc_data['description'][qc_id]
			eof[xvar].units = qc_data['unit'][qc_id].replace('(','').replace(')','').strip()
			eof[xvar].vmin = qc_data['vmin'][qc_id]
			eof[xvar].vmax = qc_data['vmax'][qc_id]
			eof[xvar].precision = qc_data['format'][qc_id]

			eof.createVariable('vsf_'+var,np.float64,('time',))
			eof['vsf_'+var].description = var+" Volume Scale Factor"
			eof['vsf_'+var][:] = vav_data[var].values
			
			eof.createVariable('column_'+var,np.float64,('time',))
			eof['column_'+var].description = var+' molecules per square meter'
			eof['column_'+var][:] = tav_data[var].values

			eof.createVariable('ada_x'+var,np.float64,('time',))
			eof['ada_x'+var].description = var+' column-average dry-air mole fraction'
			eof['ada_x'+var].units = qc_data['unit'][qc_id].replace('(','').replace(')','').strip()
			eof['ada_x'+var][:] = ada_data['x'+var].values

		# lse data
		lse_description = {'lst':'Laser sampling T','lse':'Laser sampling error','lsu':'Laser sampling U'}
		common_spec = np.intersect1d(aia_data['spectrum'],lse_data['spectrum'],return_indices=True)[2]
		for var in lse_description.keys():
			eof.createVariable(var,np.float64,('time',))
			eof[var].description = lse_description[var]
			eof[var][:] = lse_data[var][common_spec].values

		# corrections
		for var in correction_data['gas']:
			for key in correction_data.columns[1:]:
				varname = var+'_'+key
				eof.createVariable(varname,np.float64,('time',))
				eof[varname][:] = correction_data[key][list(correction_data['gas']).index(var)] # write directly

		## write data
		# update data with new scale factors and determine flags
		esf_id = 0
		nflag = 0
		for esf_id in range(esf_data['year'].size):
					
			# indices to slice the data for the concerned spectra
			start = np.sum(esf_data['n'][:esf_id])
			end = start + esf_data['n'][esf_id]	

			"""
			If new day, read in the daily error scale factors and compute
			new scale factors (RSC) as weighted averages of the a priori
			ESF factors from the pa_qc.dat file, and the daily values.
			A priori ESF values are the ratio of the xx_error/xxa scale factors
			read in from the pa_qc.dat file, with 100% uncertainties assumed.
			"""		
			for gas in gas_list:
				xgas = 'x'+gas
				qc_id = list(qc_data['variable']).index(xgas)
				apesf = qc_data['scale'][qc_id+1]/qc_data['scale'][qc_id]
				qc_data.loc[qc_data.index==qc_id+1,'rsc'] = qc_data['scale'][qc_id]*(1.0/apesf+esf_data[xgas][esf_id]/esf_data[xgas+'_error'][esf_id]**2)/(1.0/apesf**2+1.0/esf_data[xgas+'_error'][esf_id]**2)

			"""
			Look within each data record to see if any of the data values are
			outside their VMIN to VMAX range. If so, set eflag to the index of
			the variable that was furthest out of range. Then write out the data.
			"""
			eflag = np.zeros(end-start)
			kmax = np.zeros(end-start)
			dmax = np.zeros(end-start)
			for var_id,var in enumerate(aia_data.columns[mchar:-1]):

				if len(aia_data[var][aia_data[var]>=9e29]) >= 1:
					print('Missing value found (>=9e29) for variable {}.\nEnding Program'.format(var))
					print('You may need to remove missing .col files from multiggg.sh and rerun post_processing.sh')
					sys.exit()

				qc_id = list(qc_data['variable']).index(var)
				
				eof[var][start:end] = aia_data[var][start:end].values*qc_data['rsc'][qc_id]

				dev = abs( (qc_data['rsc'][qc_id]*aia_data[var][start:end].values-qc_data['vmin'][qc_id])/(qc_data['vmax'][qc_id]-qc_data['vmin'][qc_id]) -0.5 )
				
				kmax[dev>dmax] = qc_id
				dmax[dev>dmax] = dev[dev>dmax]

			eflag[dmax>0.5] = kmax[dmax>0.5]
			
			eof['flag'][start:end] = eflag

		nflag = np.count_nonzero(eof['flag'][:])

		# time		
		eof['year'][:] = np.round(aia_data['year'][:].values-aia_data['day'][:].values/365.25)
		eof['day'][:] = np.round(aia_data['day'][:].values-aia_data['hour'][:].values/24.0)

		specdate = np.array([datetime(int(aia_data['year'][i]),1,1)+timedelta(days=aia_data['day'][i]-1) for i in range(nrow)])
		eof['time'][:] = np.array([elem.total_seconds() for elem in (specdate-datetime(1970,1,1))])/(24.0*3600.0)

		# write data from col files
		for col_file in col_file_list:

			print(col_file)

			gas_XXXX = col_file.split('.')[0] # gas_XXXX, suffix for eof variable names corresponding to each .col file (i.e. VSF_h2o from the 6220 co2 window becomes co2_6220_VSF_co2)

			nhead,ncol = file_info(col_file)

			# read col_file headers
			with open(col_file,'r') as infile:
				content = infile.readlines()[1:nhead]
			gfit_version, gsetup_version = content[:2]
			gfit_version = gfit_version.strip().split()[2]
			gsetup_version = gsetup_version.strip().split()[2]

			if col_file == col_file_list[0]:
				# check that the checksums are right for the files listed in the .col file header
				checksum_dict = OrderedDict((key+'_checksum',None) for key in checksum_var_list)
				for i,line in enumerate([line for line in content if len(line.split())==2]):
					csum,fpath = line.split()
					check = checksum(fpath,csum)
					if not check:
						print('/!\\ Checksum mismatch for',fpath)
					
					checksum_dict[checksum_var_list[i]+'_checksum'] = csum

				eof['gfit_version'][:] = gfit_version
				eof['gsetup_version'][:] = gsetup_version
				for var in checksum_var_list:
					checksum_var = var+'_checksum'
					eof[checksum_var][:] = checksum_dict[checksum_var]				

			# read col_file data
			col_data = pd.read_csv(col_file,delim_whitespace=True,skiprows=nhead)
			col_data.rename(str.lower,axis='columns',inplace=True)
			col_data.rename(index=str,columns={'rms/cl':'rmsocl'},inplace=True)
			if not all(col_data['spectrum'].values == vav_data['spectrum'].values):
				print('\nMismatch between .col file spectra and .vav spectra')
				print('col file:',col_file)
				continue # contine or exit here ? Might not need to exit if we can add in the results from the faulty col file afterwards

			# create window specific variables
			for var in col_data.columns[1:]: # skip the first one ("spectrum")
				varname = '_'.join([gas_XXXX,var])
				eof.createVariable(varname,np.float64,('time',))
				if var in standard_name_dict.keys():
					eof[varname].standard_name = standard_name_dict[var]
					eof[varname].long_name = long_name_dict[var]

				eof[varname][:] = col_data[var].value