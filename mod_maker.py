#!~/anaconda2/bin/python
 # -*- coding: utf-8 -*-

"""
mod_maker10.f translated into python

What's different?

This is run with:

python mod_maker.py arg1 arg2 ar3 arg4 arg5

arg1: two letter site abbreviation (e.g. 'oc' for Lamont, Oklahoma; see the "site_dict" dictionary)
arg2: date range (YYYYMMDD-YYYYMMDD, second one not inclusive, so you don't have to worry about end of months) or single date (YYYYMMDD)
arg3: mode ('ncep', 'merradap42', 'merradap72', 'merraglob', 'fpglob', 'fpitglob'), ncep and 'glob' modes require local files
the 'merradap' modes require a .netrc file in your home directory with credentials to connect to urs.earthdata.nasa.gov
arg4: (optional, default=12:00)  hour:minute (HH:MM) for the starting time, default is local time, add 'UT' to use UTC time (15:30 will be local, 15:30UTC will be UTC)
arg5: (optional, default=24) time step in hours (can be decimal)

In GGGPATH/models/gnd it will write one mod file per day.
The merra modes require an internet connection and EarthData credentials
The ncep mode requires the global NCEP netcdf files of the given year to be present in GGGPATH/ncdf

There is dictionary of sites with their respective lat/lon, so this works for all TCCON sites, lat/lon values were taken from the wiki page of each site.

note: still need to make it work for Darwin, or any site that changed location at some point and must use different lat/lon for different periods
"""

import os, sys
import numpy as np
import netCDF4 # netcdf I/O
import re # used to parse strings
import time
import netrc # used to connect to earthdata
from datetime import datetime, timedelta
from astropy.time import Time # this is essentialy like datetime, but with better methods for conversion of datetime to / from julian dates, can also be converted to datetime
from pydap.cas.urs import setup_session # used to connect to the merra opendap servers
from pydap.client import open_url
import xarray
from urllib2 import HTTPError

def svp_wv_over_ice(temp):
	"""	
	Uses the Goff-Gratch equation to calculate the saturation vapor
	pressure of water vapor over ice at a user-specified temperature.
		Input:  temp (K)
		Output: svp (mbar)
	"""
	t0 = 273.16	# triple point temperature
	tr = t0/temp 
	yy = -9.09718*(tr-1)-3.56654*np.log10(tr)+0.876793*(1-1/tr)
	svp = 6.1173*10**yy # saturation vapor pressure over ice (mbar)

	return svp

def write_mod(mod_path,version,site_lat,data=0,surf_data=0):
	"""
	Creates a GGG-format .mod file
	INPUTS:
		mod_path: full path to write the .mod file
		version: the mod_maker version
		site_lat: site latitude (-90 to 90)
		data: dictionary of the inputs
		surf_data: dictionary of the surface inputs (for merra/geos5)
	"""

	# Define US Standard Atmosphere (USSA) for use above 10 mbar
	p_ussa=[10.0,  5.0,   2.0,   1.0,   0.1,   0.01,  0.001, 0.0001]
	t_ussa=[227.7, 239.2, 257.9, 270.6, 231.6, 198.0, 189.8, 235.0]
	z_ussa=[31.1,  36.8,  42.4,  47.8,  64.9,  79.3,  92.0,  106.3]

	if type(surf_data)==int: # ncep mode

		# The head of the .mod file	
		fmt = '{:8.3f} {:11.4e} {:7.3f} {:5.3f} {:8.3f} {:8.3f} {:8.3f}\n'
		mod_content = []
		mod_content+=[	'5  6\n',
						fmt.format(6378.137,6.000E-05,site_lat,9.81,data['H'][0],1013.25,data['TROPP']),
						version+'\n',
						' mbar        Kelvin         km      g/mole      DMF       %\n',
						'Pressure  Temperature     Height     MMW        H2O      RH\n',	]

		fmt = '{:9.3e}    {:7.3f}    {:7.3f}    {:7.4f}    {:9.3e}{:>6.1f}\n' # format for writting the lines

		# Export the Pressure, Temp and SHum for lower levels (1000 to 300 mbar)
		for k,elem in enumerate(data['H2O_DMF']):
			svp = svp_wv_over_ice(data['T'][k])
			h2o_wmf = data['H2O_DMF'][k]/(1+data['H2O_DMF'][k]) # wet mole fraction of h2o
			frh = h2o_wmf*data['T'][k]/svp # Fractional relative humidity

			# Relace H2O mole fractions that are too small
			if (frh < 30./data['T'][k]):
				print 'Replacing too-small H2O ',mod_path, data['lev'][k],h2o_wmf,svp*30./data['lev'][k]/data['lev'][k],frh,30./data['lev'][k]				
				frh = 30./data['lev'][k]
				h2o_wmf = svp*frh/data['lev'][k]
				data['H2O_DMF'][k] = h2o_wmf/(1-h2o_wmf)

			# Relace H2O mole fractions that are too large (super-saturated)  GCT 2015-08-05
			if (frh > 1.0):
				print 'Replacing too-large H2O ',mod_path,data['lev'][k],h2o_wmf,svp/data['lev'][k],frh,1.0
				frh=1.0
				h2o_wmf = svp*frh/data['lev'][k]
				data['H2O_DMF'][k] = h2o_wmf/(1-h2o_wmf)

			mmw = 28.964*(1-h2o_wmf)+18.02*h2o_wmf

			mod_content += [fmt.format(data['lev'][k], data['T'][k],data['H'][k],mmw,data['H2O_DMF'][k],100*frh)]

		# Export Pressure and Temp for middle levels (250 to 10 mbar)
		# which have no SHum reanalysis
		ptop = data['lev'][k] # Top pressure level
		frh_top = frh  # remember the FRH at the top (300 mbar) level

		for k in range(len(data['H2O_DMF']),len(data['T'])): 
			zz = np.log10(data['lev'][k])  # log10[pressure]
			strat_wmf = 7.5E-06*np.exp(-0.16*zz**2)
			svp = svp_wv_over_ice(data['T'][k])
			trop_wmf = frh_top*svp/data['lev'][k]
			wt = (data['lev'][k]/ptop)**3
			avg_wmf = trop_wmf*wt + strat_wmf*(1-wt)
			avg_frh = avg_wmf*data['lev'][k]/svp
			if (avg_frh > 1.0):
				print 'Replacing super-saturated H2O ',mod_path, data['lev'][k],avg_wmf,svp*avg_frh/data['lev'][k],avg_frh,1.0
				avg_frh = 1.0
				avg_wmf = svp*avg_frh/data['lev'][k]

			mmw = 28.964*(1-avg_wmf)+18.02*avg_wmf
			
			mod_content += [fmt.format(data['lev'][k],data['T'][k],data['H'][k],mmw,avg_wmf/(1-avg_wmf),100*avg_frh)]

		# Get the difference between the USSA and given site temperature at 10 mbar,
		Delta_T=data['T'][16]-t_ussa[0]

		# Export the P-T profile above 10mbar
		for k in range(1,len(t_ussa)):
			Delta_T=Delta_T/2
			zz = np.log10(p_ussa[k])  # log10[pressure]
			strat_wmf = 7.5E-06*np.exp(-0.16*zz**2)
			svp = svp_wv_over_ice(data['T'][k])
			mmw = 28.964*(1-strat_wmf)+18.02*strat_wmf
			mod_content += [fmt.format(p_ussa[k],t_ussa[k]+Delta_T,z_ussa[k],mmw,strat_wmf,100*strat_wmf*p_ussa[k]/svp)]

	else: # merra/geos mode

		# The head of the .mod file	
		fmt1 = '{:8.3f} {:11.4e} {:7.3f} {:5.3f} {:8.3f} {:8.3f} {:8.3f}\n'
		fmt2 = '{:9.3e}    {:7.3f}    {:7.3f}    {:7.4f}    {:9.3e}{:>6.1f}    {:9.3e}    {:9.3e}    {:9.3e}    {:9.3e}    {:7.3f}\n'
		mod_content = []
		mod_content+=[	'7  8\n',
						fmt1.format(6378.137,6.000E-05,site_lat,9.81,data['H'][0],1013.25,surf_data['TROPPB']),
						'Pressure  Temperature     Height     MMW        H2O      RH         SLP        TROPPB        TROPPV      TROPPT       TROPT\n',
						fmt2.format(*[surf_data[key] for key in ['PS','T2M','H','MMW','H2O_DMF','RH','SLP','TROPPB','TROPPV','TROPPT','TROPT']]),
						version+'\n',
						' mbar        Kelvin         km      g/mole      DMF       %       k.m+2/kg/s   kg/kg\n',
						'Pressure  Temperature     Height     MMW        H2O      RH          EPV         O3\n',	]

		fmt = '{:9.3e}    {:7.3f}    {:7.3f}    {:7.4f}    {:9.3e}{:>6.1f}    {:9.3e}    {:9.3e}\n' # format for writting the lines

		# not sure if merra needs all the filters/corrections used for ncep data?

		# Export the Pressure, Temp and SHum
		for k,elem in enumerate(data['H2O_DMF']):
			svp = svp_wv_over_ice(data['T'][k])
			h2o_wmf = data['H2O_DMF'][k]/(1+data['H2O_DMF'][k]) # wet mole fraction of h2o

			# Relace H2O mole fractions that are too large (super-saturated)  GCT 2015-08-05
			if (data['RH'][k] > 1.0):
				print 'Replacing too-large H2O ',mod_path,data['lev'][k],h2o_wmf,svp/data['T'][k],data['RH'][k],1.0
				data['RH'][k] = 1.0
				h2o_wmf = svp*data['RH'][k]/data['T'][k]
				data['H2O_DMF'][k] = h2o_wmf/(1-h2o_wmf)

			mmw = 28.964*(1-h2o_wmf)+18.02*h2o_wmf

			mod_content += [fmt.format(data['lev'][k], data['T'][k],data['H'][k],mmw,data['H2O_DMF'][k],100*data['RH'][k],data['EPV'][k],data['O3'][k])]

	with open(mod_path,'w') as outfile:
		outfile.writelines(mod_content)

def trilinear_interp(DATA,varlist,site_lon_360,site_lat,site_tim):
	"""
	Evaluates  fout = fin(xx,yy,*,tt) 
	Result is a 1-vector
	"""
	INTERP_DATA = {}

	dx = DATA['lon'][1]-DATA['lon'][0]
	dy = DATA['lat'][1]-DATA['lat'][0]
	dt = DATA['time'][1]-DATA['time'][0]

	xx = (site_lon_360-DATA['lon'][0])/dx
	yy = (site_lat-DATA['lat'][0])/dy
	tt = (site_tim-DATA['time'][0])/dt

	nxx =  len(DATA['lon'])
	nyy =  len(DATA['lat'])
	ntt =  len(DATA['time'])

	index_xx = int(xx)
	ixpomnxx = (index_xx+1) % nxx
	fr_xx = xx-index_xx	

	index_yy = int(yy)
	if index_yy > nyy-2:
		index_yy = nyy-2  #  avoid array-bound violation at SP
	fr_yy = yy-index_yy

	index_tt = int(tt)
	
	if index_tt < 0:
		index_tt = 0          # Prevent Jan 1 problem
	if index_tt+1 > ntt-1:
		index_tt = ntt-2  # Prevent Dec 31 problem

	fr_tt=tt-index_tt  #  Should be between 0 and 1 when interpolating in time

	if (fr_tt < -1) or (fr_tt > 2):
	   print 'Excessive time extrapolation:',fr_tt,' time-steps   =',fr_tt*dt,' days'
	   print ' tt= ',tt,'  index_tt=',index_tt,'  fr_tt=',fr_tt
	   print 'input file does not cover the full range of dates'
	   print 'site_tim',site_tim
	   print 'tim_XX',DATA['time']
	   raw_input() # will hold the program until something is typed in commandline

	if (fr_xx < 0) or (fr_xx > 1):
	   print 'Excessive longitude extrapolation:',fr_xx,' steps   =',fr_xx*dx,' deg'
	   print ' xx= ',xx,'  index_xx=',index_xx,'  fr_xx=',fr_xx
	   print 'input file does not cover the full range of longitudes'
	   raw_input() # will hold the program until something is typed in commandline

	if (fr_yy < 0) or (fr_yy > 1):
	   print 'Excessive latitude extrapolation:',fr_yy-1,' steps   =',(fr_yy-1)*dy,' deg'
	   print ' yy= ',yy,'  index_yy=',index_yy,'  fr_yy=',fr_yy
	   print 'input file does not cover the full range of latitudes'
	   raw_input() # will hold the program until something is typed in commandline

	if (fr_tt < 0) or (fr_tt > 1):
		print ' Warning: time extrapolation of ',fr_tt,' time-steps'
	if (fr_xx < 0) or (fr_xx > 1):
		print ' Warning: longitude extrapolation of ',fr_xx,' steps'
	if (fr_yy < 0) or (fr_yy > 1):
		print ' Warning: latitude extrapolation of ',fr_yy,' steps'

	for varname in varlist:

		fin = DATA[varname]
		
		if fin.ndim==4:
			fout =	((fin[index_tt,:,index_yy,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt,:,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
			+ (fin[index_tt,:,index_yy+1,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt,:,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*(1.0-fr_tt) \
			+ ((fin[index_tt+1,:,index_yy,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt+1,:,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
			+ (fin[index_tt+1,:,index_yy+1,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt+1,:,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*fr_tt
		elif fin.ndim==3: # for data that do not have the vertical dimension
			fout =	((fin[index_tt,index_yy,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
			+ (fin[index_tt,index_yy+1,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*(1.0-fr_tt) \
			+ ((fin[index_tt+1,index_yy,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt+1,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
			+ (fin[index_tt+1,index_yy+1,index_xx]*(1.0-fr_xx) \
			+ fin[index_tt+1,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*fr_tt
		else:
			print 'Data has unexpected dimensions, ndim =',fin.ndim
			sys.exit()

		INTERP_DATA[varname] = fout*DATA['scale_factor_'+varname] + DATA['add_offset_'+varname]

	return INTERP_DATA

def read_data(dataset, varlist, lat_lon_box=0):
	"""
	for ncep files "dataset" is the full path to the netcdf file

	for merra files "dataset" is a pydap.model.DatasetType object
	"""
	DATA = {}

	opendap = type(dataset)!=netCDF4._netCDF4.Dataset

	if lat_lon_box == 0: # ncep mode

		varlist += ['level','lat','lon','time']

		for varname in varlist:
			DATA[varname] = dataset[varname][:]

			for attribute in ['add_offset','scale_factor']:
				try:
					DATA['add_offset_'+varname] = dataset[varname].getncattr('add_offset')
				except:
					DATA['add_offset_'+varname] = 0.0
					DATA['scale_factor_'+varname] = 1.0

	else: # merra/geos5 mode
		
		min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID = lat_lon_box

		DATA['lat'] = dataset['lat'][min_lat_ID:max_lat_ID] 	# Read in variable 'lat'
		DATA['lon'] = dataset['lon'][min_lon_ID:max_lon_ID] 	# Read in variable 'lon'
		DATA['time'] = dataset['time'][:]	# Read in variable 'time'

		try:
			dataset['lev']
		except:
			pass
		else:
			if dataset['lev'].shape[0] == 72:
				DATA['lev'] = dataset['PL'][:,:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID]	# the merra 72 mid level pressures are not fixed
			elif dataset['lev'].shape[0] == 42:
				DATA['lev'] = dataset['lev'][:]	# the 42 levels data is on a fixed pressure grid
			else:
				DATA['lev'] = dataset['PS'][:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID] # surface data doesn't have a 'lev' variable

		if opendap:
			for varname in ['time','lev','lat','lon']:
				try:
					DATA[varname] = DATA[varname].data	
				except KeyError,IndexError:
					pass

		# get longitudes as 0 -> 360 instead of -180 -> 180, needed for trilinear_interp
		for i,elem in enumerate(DATA['lon']):
			if elem < 0:
				DATA['lon'][i] = elem + 360.0

		for varname in varlist:

			if dataset[varname].ndim == 4:
				DATA[varname] = dataset[varname][:,:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID] 	# Read in variable varname
			else:
				DATA[varname] = dataset[varname][:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID] 	# Read in variable varname
			if opendap or ('Masked' in str(type(DATA[varname]))):
				DATA[varname] = DATA[varname].data

			for attribute in ['add_offset','scale_factor']:
				try:
					DATA[attribute+'_'+varname] = dataset[varname].getncattr(attribute)
				except:
					DATA['add_offset_'+varname] = 0.0
					DATA['scale_factor_'+varname] = 1.0

			print varname, DATA[varname].shape

	time_units = dataset['time'].units  # string containing definition of time units

	# two lines to parse the date (no longer need to worry about before/after 2014)
	date_list = re.findall(r"[\w]+",time_units.split(' since ')[1])
	common_date_format = '{:0>4}-{:0>2}-{:0>2} {:0>2}:{:0>2}:{:0>2}'.format(*date_list)

	start_date = datetime.strptime(common_date_format,'%Y-%m-%d %H:%M:%S')
	astropy_start_date = Time(start_date)

	DATA['julday0'] = astropy_start_date.jd # gives same results as IDL's JULDAY function

	return DATA

def querry_indices(dataset,site_lat,site_lon_180,box_lat_half_width,box_lon_half_width):
	"""	
	Set up a lat-lon box for the data querry

	Unlike with ncep, this will only use daily files for interpolation, so no time box is defined
	
	NOTE: merra lat -90 -> +90 ;  merra lon -180 -> +179.375

	To be certain to get two points on both side of the site lat and lon, use the grid resolution
	"""
	# define 2xbox_lat_half_width°x2xbox_lon_half_width° lat-lon box centered on the site lat-lon
	min_lat, max_lat = site_lat-box_lat_half_width, site_lat+box_lat_half_width
	min_lon, max_lon = site_lon_180-box_lon_half_width, site_lon_180+box_lon_half_width
	
	# handle edge cases
	if min_lat < -90:
		min_lat = -(90-abs(90-min_lat))
	if max_lat > 90:
		max_lat = 90-abs(90-max_lat)

	if max_lat < min_lat:
		swap = max_lat
		max_lat = min_lat
		min_lat = swap
	
	if min_lon < -180:
		min_lon = 180 - abs(180-min_lon)
	if max_lon > 180:
		max_lon = - (180 - abs(180-max_lon))

	if max_lon < min_lon:
		swap = max_lon
		max_lon = min_lon
		min_lon = swap

	# read the latitudes and longitudes from the merra file
	if type(dataset)==netCDF4._netCDF4.Dataset:
		merra_lon = dataset['lon'][:]
		merra_lat = dataset['lat'][:]
	else: # for opendap datasets
		merra_lon = dataset['lon'][:].data
		merra_lat = dataset['lat'][:].data

	# get the indices of merra longitudes and latitudes that fit in the lat-lon box
	merra_lon_in_box_IDs = np.where((merra_lon>=min_lon) & (merra_lon<=max_lon))[0]
	merra_lat_in_box_IDs = np.where((merra_lat>=min_lat) & (merra_lat<=max_lat))[0]

	min_lat_ID, max_lat_ID = merra_lat_in_box_IDs[0], merra_lat_in_box_IDs[-1]+1
	min_lon_ID, max_lon_ID = merra_lon_in_box_IDs[0], merra_lon_in_box_IDs[-1]+1
	# +1 because ARRAY[i:j] in python will return elements i to j-1

	return [min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID]

# ncep has geopotential height profiles, not merra(?, only surface), so I need to convert geometric heights to geopotential heights
# the idl code uses a fixed radius for the radius of earth (6378.137 km), below the gravity routine of gsetup is used
# also the surface geopotential height of merra is in units of m2 s-2, so it must be divided by surface gravity
def gravity(gdlat,altit):
	"""
	copy/pasted from fortran routine comments
	This is used to convert

	Input Parameters:
	    gdlat       GeoDetric Latitude (degrees)
	    altit       Geometric Altitude (km)
	
	Output Parameter:
	    gravity     Effective Gravitational Acceleration (m/s2)
	    radius 		Radius of earth at gdlat
	
	Computes the effective Earth gravity at a given latitude and altitude.
	This is the sum of the gravitational and centripital accelerations.
	These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
	The Earth is assumed to be an oblate ellipsoid, with a ratio of the
	major to minor axes = sqrt(1+con) where con=.006738
	This eccentricity makes the Earth's gravititational field smaller at
	the poles and larger at the equator than if the Earth were a sphere
	of the same mass. [At the equator, more of the mass is directly
	below, whereas at the poles more is off to the sides). This effect
	also makes the local mid-latitude gravity field not point towards
	the center of mass.
	
	The equation used in this subroutine agrees with the International
	Gravitational Formula of 1967 (Helmert's equation) within 0.005%.
	
	Interestingly, since the centripital effect of the Earth's rotation
	(-ve at equator, 0 at poles) has almost the opposite shape to the
	second order gravitational field (+ve at equator, -ve at poles),
	their sum is almost constant so that the surface gravity could be
	approximated (.07%) by the simple expression g=0.99746*GM/radius^2,
	the latitude variation coming entirely from the variation of surface
	r with latitude. This simple equation is not used in this subroutine.
	"""

	d2r=3.14159265/180.0	# Conversion from degrees to radians
	gm=3.9862216e+14  		# Gravitational constant times Earth's Mass (m3/s2)
	omega=7.292116E-05		# Earth's angular rotational velocity (radians/s)
	con=0.006738       		# (a/b)**2-1 where a & b are equatorial & polar radii
	shc=1.6235e-03  		# 2nd harmonic coefficient of Earth's gravity field 
	eqrad=6378178.0   		# Equatorial Radius (m)

	gclat=np.arctan(np.tan(d2r*gdlat)/(1.0+con))  # radians

	radius=1000.0*altit+eqrad/np.sqrt(1.0+con*np.sin(gclat)**2)
	ff=(radius/eqrad)**2
	hh=radius*omega**2
	ge=gm/eqrad**2                      # = gravity at Re

	gravity=(ge*(1-shc*(3.0*np.sin(gclat)**2-1)/ff)/ff-hh*np.cos(gclat)**2)*(1+0.5*(np.sin(gclat)*np.cos(gclat)*(hh/ge+2.0*shc/ff**2))**2)

	return gravity, radius

def read_merradap(username,password,mode,site_lon_180,gravity_at_lat,date,end_date,time_step,varlist,surf_varlist):
	"""
	Read MERRA2 data via opendap.

	This has to connect to the daily netcdf files, and then concatenate the subsetted datasets.

	This is EXTREMELY slow, should probably make that use separate to generate local files, and then use to files in mod_maker
	"""
	DATA = {}
	SURF_DATA = {}

	if '42' in mode:
		letter = 'P'
	elif '72' in mode:
		letter = 'V'
		varlist += ['PL']

	old_UTC_date = ''
	urllist = []
	surface_urllist = []
	print '\n\t-Making lists of URLs'
	while date < end_date:
		UTC_date = date + timedelta(hours = -site_lon_180/15.0) # merra times are in UTC, so the date may be different than the local date, make sure to use the UTC date to querry the file
		if (UTC_date.strftime('%Y%m%d') != old_UTC_date):
			print '\t\t',UTC_date.strftime('%Y-%m-%d')
			urllist += ['https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2I3N{}ASM.5.12.4/{:0>4}/{:0>2}/MERRA2_400.inst3_3d_asm_N{}.{:0>4}{:0>2}{:0>2}.nc4'.format(letter,UTC_date.year,UTC_date.month,letter.lower(),UTC_date.year,UTC_date.month,UTC_date.day)]
			surface_urllist += ['https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2I1NXASM.5.12.4/{:0>4}/{:0>2}/MERRA2_400.inst1_2d_asm_Nx.{:0>4}{:0>2}{:0>2}.nc4'.format(UTC_date.year,UTC_date.month,UTC_date.year,UTC_date.month,UTC_date.day)]
			if old_UTC_date == '':
				session = setup_session(username,password,check_url=urllist[0]) # just need to setup the authentication session once
		old_UTC_date = UTC_date.strftime('%Y%m%d')
		date = date + time_step

	# multi-level data
	print '\nNow doing multi-level data'
	print '\t-Connecting to datasets ...'
	store_list = [xarray.backends.PydapDataStore.open(url,session) for url in urllist]
	dataset_list = [xarray.open_dataset(store) for store in store_list]
	print '\t-Datasets opened'
	min_lat_ID,max_lat_ID,min_lon_ID,max_lon_ID = querry_indices(dataset_list[0],site_lat,site_lon_180,2.5,2.5) # just need to get the lat/lon box once
	subsest_dataset_list = [dataset[{'lat':range(min_lat_ID,max_lat_ID+1),'lon':range(min_lon_ID,max_lon_ID+1)}] for dataset in dataset_list]
	print '\t-Datasets subsetted'
	print '\t-Merging datasets (time consuming)'
	merged_dataset = xarray.concat(subsest_dataset_list,'time')
	merged_dataset = merged_dataset.fillna(1e15)

	# single-level data
	print '\nNow doing single-level data'
	print '\t-Connecting to datasets ...'
	surface_store_list = [xarray.backends.PydapDataStore.open(url,session) for url in surface_urllist]
	surface_dataset_list = [xarray.open_dataset(store) for store in surface_store_list]
	print '\t-Datasets opened'
	subsest_surface_dataset_list = [dataset[{'lat':range(min_lat_ID,max_lat_ID+1),'lon':range(min_lon_ID,max_lon_ID+1)}] for dataset in surface_dataset_list]
	print '\t-Datasets subsetted'
	print '\t-Merging datasets (time consuming)'
	merged_surface_dataset = xarray.concat(subsest_surface_dataset_list,'time')
	merged_surface_dataset = merged_surface_dataset.fillna(1e15)

	for varname in varlist:
		DATA[varname] = merged_dataset[varname].data	
		DATA['add_offset_'+varname] = 0.0
		DATA['scale_factor_'+varname] = 1.0
	for varname in surf_varlist:
		SURF_DATA[varname] = merged_surface_dataset[varname].data
		SURF_DATA['add_offset_'+varname] = 0.0
		SURF_DATA['scale_factor_'+varname] = 1.0

	for varname in ['time','lat','lon']:
		DATA[varname] = merged_dataset[varname].data
		SURF_DATA[varname] = merged_surface_dataset[varname].data
	DATA['lev'] = merged_dataset['lev'].data

	DATA['PHIS'] = DATA['PHIS'] / gravity_at_lat # convert from m2 s-2 to m

	delta_time = [(i-DATA['time'][0]).astype('timedelta64[h]') / np.timedelta64(1,'h') for i in DATA['time']] # hours since base time
	surf_delta_time = [(i-SURF_DATA['time'][0]).astype('timedelta64[h]') / np.timedelta64(1,'h') for i in SURF_DATA['time']] # hours since base time

	DATA['julday0'] = Time(str(DATA['time'][0]),format="isot").jd
	SURF_DATA['julday0'] = Time(str(SURF_DATA['time'][0]),format="isot").jd

	DATA['time'] = delta_time
	SURF_DATA['time'] = surf_delta_time

	# get longitudes as 0 -> 360 instead of -180 -> 180, needed for trilinear_interp
	for i,elem in enumerate(DATA['lon']):
		if elem < 0:
			DATA['lon'][i] = elem + 360.0
	for i,elem in enumerate(SURF_DATA['lon']):
		if elem < 0:
			SURF_DATA['lon'][i] = elem + 360.0

	return DATA, SURF_DATA

def read_ncep(ncdf_path,year):
	"""
	Read data from yearly NCEP netcdf files and return it in one dictionary
	"""

	# path to the netcdf files
	ncdf_AT_file = os.path.join(ncdf_path,'.'.join(['air','{:0>4}'.format(year),'nc']))
	ncdf_GH_file = os.path.join(ncdf_path,'.'.join(['hgt','{:0>4}'.format(year),'nc']))
	ncdf_SH_file = os.path.join(ncdf_path,'.'.join(['shum','{:0>4}'.format(year),'nc']))

	print 'Read global',year,'NCEP data ...'
	# Air Temperature
	DATA = read_data(netCDF4.Dataset(ncdf_AT_file,'r'), ['air'])
	if len(DATA['air']) < 17:
		print 'Need 17 levels of AT data: found only ',len(lev_AT)

	# Specific Humidity
	SHUM_DATA = read_data(netCDF4.Dataset(ncdf_SH_file,'r'), ['shum'])
	if len(SHUM_DATA['level']) <  8:
		print 'Need  8 levels of SH data: found only ',len(lev_SH)

	if list(SHUM_DATA['level'])!=list(DATA['level'][:len(SHUM_DATA['level'])]):
		print 'Warning: air and shum do not share the same lower pressure levels'
	
	DATA.update(SHUM_DATA)
	
	# Geopotential Height
	GH_DATA = read_data(netCDF4.Dataset(ncdf_GH_file,'r'), ['hgt'])
	if len(GH_DATA['level']) < 17:
		print 'Need 17 levels of GH data: found only ',len(lev_GH)
	
	DATA.update(GH_DATA)

	for key in DATA:
		if 'air' in key:
			DATA[key.replace('air','T')] = DATA[key]
			del DATA[key]
		if 'hgt' in key:
			DATA[key.replace('hgt','H')] = DATA[key]
			del DATA[key]
		if 'shum' in key:
			DATA[key.replace('shum','QV')] = DATA[key]
			del DATA[key]
	
	DATA['lev'] = DATA['level']
	del DATA['level']

	return DATA

def read_global(ncdf_path,mode,site_lat,site_lon_180,gravity_at_lat,varlist,surf_varlist):
	"""
	Read data from GEOS5 and MERRA2 datasets

	This assumes those are saved locally in GGGPATH/ncdf with two files per dataset (inst3_3d_asm_np and inst3_2d_asm_nx)
	"""

	key_dict = {'merraglob':'MERRA','fpglob':'_fp_','fpitglob':'_fpit_'}

	# assumes only one file with all the data exists in the GGGPATH/ncdf folder
	# path to the netcdf file
	ncdf_list = [i for i in os.listdir(ncdf_path) if key_dict[mode] in i]
	
	ncdf_file = ncdf_list[0]
	dataset = netCDF4.Dataset(os.path.join(ncdf_path,ncdf_file),'r')
	
	surf_file = ncdf_list[1]
	surface_dataset = netCDF4.Dataset(os.path.join(ncdf_path,surf_file),'r')

	# get the min/max lat-lon indices of merra lat-lon that lies within a given box.
	if 'fpglob' in mode: # geos5-fp has a smaller grid than merra2 amd geos5-fp-it
		box_lat_half_width = 0.250001
		box_lon_half_width = 0.312501
	else:
		box_lat_half_width = 0.500001
		box_lon_half_width = 0.625001
	lat_lon_box = querry_indices(dataset,site_lat,site_lon_180,box_lat_half_width,box_lon_half_width)

	# multi-level data
	print 'Read global',mode,'multi-level data ...'
	DATA = read_data(dataset,varlist,lat_lon_box)
	DATA['PHIS'] = DATA['PHIS'] / gravity_at_lat # convert from m2 s-2 to m

	# single level data
	print 'Read global',mode,'single-level data ...'
	SURF_DATA = read_data(surface_dataset,surf_varlist,lat_lon_box)
	
	# merra/geos time is minutes since base time, need to convert to hours
	DATA['time'] = DATA['time'] / 60.0 
	SURF_DATA['time'] = SURF_DATA['time'] / 60.0

	return DATA,SURF_DATA

if __name__ == "__main__": # this is only executed when the code is used directly (e.g. not executed when imported from another python code)

	# dictionary mapping TCCON site abbreviations to their lat-lon-alt data, and full names
	site_dict = {
				'pa':{'name': 'Park Falls','loc':'Wisconsin, USA','lat':45.945,'lon':269.727,'alt':442},
				'oc':{'name': 'Lamont','loc':'Oklahoma, USA','lat':36.604,'lon':262.514,'alt':320},
				'wg':{'name': 'Wollongong','loc':'Australia','lat':-34.406,'lon':150.879,'alt':30},
				'db':{'name': 'Darwin','loc':'Australia','lat':-12.45606,'lon':130.92658,'alt':37},
				#,tuple([datetime(2005,8,1),datetime(2015,7,1)]):{'lat':-12.422445,'lon':130.89154,'alt':30},
				#										tuple([datetime(2015,7,1)]):{'lat':-12.45606,'lon':130.92658,'alt':37}
				#		},
				'or':{'name': 'Orleans','loc':'France','lat':47.97,'lon':2.113,'alt':130},
				'bi':{'name': 'Bialystok','loc':'Poland','lat':53.23,'lon':23.025,'alt':180},
				'br':{'name': 'Bremen','loc':'Germany','lat':53.1037,'lon':8.849517,'alt':30},
				'jc':{'name': 'JPL 01','loc':'California, USA','lat':34.202,'lon':241.825,'alt':390},
				'jf':{'name': 'JPL 02','loc':'California, USA','lat':34.202,'lon':241.825,'alt':390},
				'ra':{'name': 'Reunion Island','loc':'France','lat':-20.901,'lon':55.485,'alt':87},
				'gm':{'name': 'Garmisch','loc':'Germany','lat':47.476,'lon':11.063,'alt':743},
				'lh':{'name': 'Lauder 01','loc':'New Zealand','lat':-45.038,'lon':169.684,'alt':370},
				'll':{'name': 'Lauder 02','loc':'New Zealand','lat':-45.038,'lon':169.684,'alt':370},
				'tk':{'name': 'Tsukuba 02','loc':'Japan','lat':63.0513,'lon':140.1215,'alt':31},
				'ka':{'name': 'Karlsruhe','loc':'Germany','lat':49.1002,'lon':8.4385,'alt':119},
				'ae':{'name': 'Ascenssion Island','loc':'United Kingdom','lat':-7.933333,'lon':345.583333,'alt':0},
				'eu':{'name': 'Eureka','loc':'Canada','lat':80.05,'lon':273.58,'alt':610},
				'so':{'name': 'Sodankyla','loc':'Finland','lat':67.3668,'lon':26.6310,'alt':188},
				'iz':{'name': 'Izana','loc':'Spain','lat':28.0,'lon':344.0,'alt':0},
				'if':{'name': 'Idianapolis','loc':'Indiana, USA','lat':39.861389,'lon':273.996389,'alt':270},
				'df':{'name': 'Dryden','loc':'California, USA','lat':34.959917,'lon':242.118931,'alt':700},
				'js':{'name': 'Saga','loc':'Japan','lat':33.240962,'lon':130.288239,'alt':7},
				'fc':{'name': 'Four Corners','loc':'USA','lat':36.79749,'lon':251.51991,'alt':1643},
				#'ci':{'name': 'Pasadena','loc':'California, USA','lat':34.13623,'lon':241.873103,'alt':230},
				'ci':{'name': 'Pasadena','loc':'California, USA','lat':34.136,'lon':241.873,'alt':230},
				'rj':{'name': 'Rikubetsu','loc':'Japan','lat':43.4567,'lon':143.7661,'alt':380},
				'pr':{'name': 'Paris','loc':'France','lat':48.846,'lon':2.356,'alt':60},
				'ma':{'name': 'Manaus','loc':'Brazil','lat':-3.2133,'lon':299.4017,'alt':50},
				'sp':{'name': 'Ny-Alesund','loc':'Norway','lat':78.92324,'lon':11.92298,'alt':0},
				'et':{'name': 'East Trout Lake','loc':'Canada','lat':54.353738,'lon':255.013333,'alt':501.8},
				'an':{'name': 'Anmyeondo','loc':'Korea','lat':36.5382,'lon':126.331,'alt':30},
				'bu':{'name': 'Burgos','loc':'Philippines','lat':18.5325,'lon':120.6496,'alt':35},
				'we':{'name': 'Jena','loc':'Austria','lat':50.91,'lon':11.57,'alt':211.6},
				}

	GGGPATH = os.environ['GGGPATH'] # reads the GGGPATH environment variable
	print 'GGGPATH =',GGGPATH

	argu = sys.argv # list of commandline arguments, argu[0] will be "mod_maker.py"

	site_abbrv = argu[1]	# two letter site abbreviation
	try:
		print 'Site:',site_dict[site_abbrv]['name'],site_dict[site_abbrv]['loc']
	except KeyError:
		print 'Wrong 2 letter site abbreviation (check the site_dict dictionary)'
		sys.exit()
	print 'lat,lon,masl:',site_dict[site_abbrv]['lat'],site_dict[site_abbrv]['lon'],site_dict[site_abbrv]['alt']

	# parse the selected range of dates for which .mod files will be generated
	date_range = argu[2].split('-')
	start_date = datetime.strptime(date_range[0],'%Y%m%d')
	try:
		end_date = datetime.strptime(date_range[1],'%Y%m%d')
	except IndexError: # if a single date is given set the end date as the next day
		end_date = start_date + timedelta(days=1)
	if start_date>=end_date:
		print 'Error: the second argument must be a date range YYYYMMDD-YYYYMMDD or a single date YYYYMMDD'
		sys.exit()
	print 'Date range: from',start_date.strftime('%Y-%m-%d'),'to',end_date.strftime('%Y-%m-%d')

	mode = argu[3].lower() # ncep or merra
	if False not in [elem not in mode for elem in ['merradap','merraglob','fpglob','fpitglob','ncep']]:
		print 'Wrong mode, must be one of [ncep, merradap42, merradap72, merraglob, fpglob, fpitglob]'
		sys.exit()
	if 'merradap' in mode: # get the earthdata credentials
		try:
			username,account,password = netrc.netrc().authenticators('urs.earthdata.nasa.gov')
		except:
			print 'When using MERRA mode, you need a ~/.netrc file to connect to urs.earthdata.nasa.gov'
			sys.exit()

	print 'Mode:',mode.upper()

	simple = {'merradap42':'merra','merradap72':'merra','merraglob':'merra','ncep':'ncep','fpglob':'fp','fpitglob':'fpit'}
	mod_path = os.path.join(GGGPATH,'models','gnd','comparison',simple[mode],site_abbrv)	# .mod files will be saved here
	if not os.path.exists(mod_path):
		os.makedirs(mod_path)
	print 'MOD files will be saved in:',mod_path

	# hour and minute for time interpolation, default is local noon
	UTC = False
	if len(argu)>4:
		time_input = re.search('([0-9][0-9]):([0-9][0-9])(UTC)?',argu[4]).groups()
		if 'UTC' in time_input:
			UTC = True
		HH = int(time_input[0])
		MM = int(time_input[1])
	else:
		HH = 12
		MM = 0

	# small checks
	if HH>=24 or HH<0:
		print 'Need 0<=H<24'
		sys.exit()
	if MM>=60 or MM<0:
		print 'Need 0<=MM<60'
		sys.exit()

	if len(argu)>5:
		time_step = int(argu[5])
	else:
		time_step = 24

	# will need to handle sites that changed location at some point
	site_moved = False
	try:
		site_lat = site_dict[site_abbrv]['lat']
	except KeyError:
		site_moved = True
	else:
		site_lon_360 = site_dict[site_abbrv]['lon']
		if site_lon_360 > 180:
			site_lon_180 = site_lon_360-360
		else:
			site_lon_180 = site_lon_360
		site_alt = site_dict[site_abbrv]['alt']

	# directions for .mod file name
	if site_lat > 0:
		ns = 'N'
	else:
		ns = 'S'

	if site_lon_180>0:
		ew = 'E'
	else:
		ew = 'W'

	rmm = 28.964/18.02	# Ratio of Molecular Masses (Dry_Air/H2O)
	gravity_at_lat, earth_radius_at_lat = gravity(site_lat,site_alt/1000.0) # used in merra/fp mode

	ncdf_path = os.path.join(GGGPATH,'ncdf')

	if 'ncep' in mode:
		DATA = read_ncep(ncdf_path,start_date.year)
		varlist = ['T','H','QV']
	elif 'glob' in mode:
		varlist = ['T','QV','RH','H','EPV','O3','PHIS']
		surf_varlist = ['T2M','QV2M','PS','SLP','TROPPB','TROPPV','TROPPT','TROPT']	
		DATA,SURF_DATA = read_global(ncdf_path,mode,site_lat,site_lon_180,gravity_at_lat,varlist,surf_varlist)

	if UTC:
		date = start_date + timedelta(hours=HH+site_lon_180/15.0,minutes=MM) # date with local time
	else:
		date = start_date + timedelta(hours=HH,minutes=MM) # date with local time
	astropy_date = Time(date)
	print 'Starting local time for interpolation:',date.strftime('%Y-%m-%d %H:%M')

	time_step = timedelta(hours=time_step) # time step between mod files; will need to change the mod file naming and gsetup to do sub-daily files
	print 'Time step:',time_step.total_seconds()/3600.0,'hours'

	if 'merradap' in mode: # read all the data first, this could take a while ...	
		print 'Reading MERRA2 data via opendap'
		varlist = ['T','QV','RH','H','EPV','O3','PHIS']
		surf_varlist = ['T2M','QV2M','PS','SLP','TROPPB','TROPPV','TROPPT','TROPT']	
		DATA,SURF_DATA = read_merradap(username,password,mode,site_lon_180,gravity_at_lat,date,end_date,time_step,varlist,surf_varlist)

	new_year = False
	while date<end_date:

		# use the local date for the name of the .mod file
		YYYYMMDD = date.strftime('%Y%m%d')
		HHMM = date.strftime('%H%M')
		if time_step < timedelta(days=1):
			mod_name = '{}_{}_{:0>2.0f}{:>1}_{:0>3.0f}{:>1}.mod'.format(YYYYMMDD,HHMM,round(abs(site_lat)),ns,round(abs(site_lon_180)),ew)
		else:
			mod_name = '{}_{:0>2.0f}{:>1}_{:0>3.0f}{:>1}.mod'.format(YYYYMMDD,round(abs(site_lat)),ns,round(abs(site_lon_180)),ew)
		mod_file_path = os.path.join(mod_path,mod_name)
		print '\n',mod_name
		
		"""
		Interpolation time:
			julday0 is the fractional julian day number of the base time of the dataset: dataset times are in UTC hours since base time
			astropy_date.jd is the fractional julian day number of the current local day
			(astropy_date.jd-julday0)*24.0 = local hours since julday0
		"""
		site_tim = (astropy_date.jd-DATA['julday0'])*24.0 - site_lon_180/15.0 # UTC hours since julday0
		# interpolate the data to the site's location and the desired time
		INTERP_DATA = trilinear_interp(DATA,varlist,site_lon_360,site_lat,site_tim)

		if 'ncep' in mode:
			INTERP_DATA['lev'] = np.copy(DATA['lev'])
			INTERP_DATA['TROPP'] = 0  # tropopause pressure not used with NCEP data
			INTERP_DATA['RH'] = 0 # won't be used, just to feed something to write_mod frh
		else: # merra/geos5
			if 'lev' not in varlist:
				INTERP_DATA['lev'] = np.copy(DATA['lev'])
			
			# get rid of fill values
			without_fill_IDs = np.where(INTERP_DATA['T']<1e10) # merra/geos fill value is 1e15
			for varname in list(set(varlist+['lev'])):
				try:
					INTERP_DATA[varname] = INTERP_DATA[varname][without_fill_IDs]
				except IndexError:
 					pass

		if ('merradap' in mode) or ('glob' in mode):
			site_tim = (astropy_date.jd-SURF_DATA['julday0'])*24.0 - site_lon_180/15.0 # UTC hours since julday0
			INTERP_SURF_DATA = trilinear_interp(SURF_DATA,surf_varlist,site_lon_360,site_lat,site_tim)
			for varname in ['PS','SLP','TROPPB','TROPPV','TROPPT']:
				INTERP_SURF_DATA[varname] = INTERP_SURF_DATA[varname] / 100.0 # convert Pa to hPa

			if 'merradap72' in mode: # merra42 and ncep go from high pressure to low pressure, but merra 72 does the reverse
				# reverse merra72 profiles
				INTERP_DATA['lev'] = INTERP_DATA['PL'] / 100.0
				for varname in list(set(varlist+['lev'])):
					try:
						INTERP_DATA[varname] = INTERP_DATA[varname][::-1]
					except IndexError:
						pass

		INTERP_DATA['H2O_DMF'] = rmm*INTERP_DATA['QV']/(1-INTERP_DATA['QV']) # Convert specific humidity, a wet mass mixing ratio, to dry mole fraction
		INTERP_DATA['H'] = INTERP_DATA['H']/1000.0	# Convert m to km

		if ('merradap' in mode) or ('glob' in mode):
			INTERP_SURF_DATA['H2O_DMF'] = rmm*INTERP_SURF_DATA['QV2M']/(1-INTERP_SURF_DATA['QV2M'])
			INTERP_DATA['PHIS'] = INTERP_DATA['PHIS']/1000.0
			# compute surface relative humidity
			svp = svp_wv_over_ice(INTERP_SURF_DATA['T2M'])
			INTERP_SURF_DATA['H2O_WMF'] = INTERP_SURF_DATA['H2O_DMF']/(1+INTERP_SURF_DATA['H2O_DMF']) # wet mole fraction of h2o
			INTERP_SURF_DATA['RH'] = 100*INTERP_SURF_DATA['H2O_WMF']*INTERP_SURF_DATA['PS']/svp # Fractional relative humidity
			INTERP_SURF_DATA['MMW'] = 28.964*(1-INTERP_SURF_DATA['H2O_WMF'])+18.02*INTERP_SURF_DATA['H2O_WMF']
			INTERP_SURF_DATA['H'] = INTERP_DATA['PHIS']

		# write the .mod file
		version = 'mod_maker_10.6   2017-04-11   GCT'
		if 'ncep' in mode:
			write_mod(mod_file_path,version,site_lat,data=INTERP_DATA)
		else:
			write_mod(mod_file_path,version,site_lat,data=INTERP_DATA,surf_data=INTERP_SURF_DATA)

		if ((date+time_step).year!=date.year):
			new_year = True

		date = date + time_step
		astropy_date = Time(date)
		
		if ('ncep' in mode) and new_year:
			DATA = read_ncep(ncdf_path,date.year)
			new_year = False
