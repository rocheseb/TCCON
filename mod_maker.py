#!~/anaconda2/bin/python
 # -*- coding: utf-8 -*-

"""
mod_maker10.f translated into python

What's different?

This is run with:

python mod_maker.py arg1 arg2 ar3 arg4 arg 5

arg1: two letter site abbreviation (e.g. 'oc' for Lamont, Oklahoma; see the "site_dict" dictionary)
arg2: date range (YYYYMMDD-YYYYMMDD, second one not inclusive, so you don't have to worry about end of months) or single date (YYYYMMDD)
arg3: mode ('ncep', 'merradap42', 'merradap72', 'merraglob', 'fpglob', 'fpitglob'), ncep and 'glob' modes require local files
the 'merradap' modes require a .netrc file in your home directory with credentials to connect to urs.earthdata.nasa.gov
arg4: (optional, default=12) local hour (HH) for time interpolation (12 is local noon)
arg5: (optional, default=0) minute (MM) for time interpolation

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

def write_mod(mod_path,version,site_lat,lev_AT,sat,sgh,site_TP,h2o_dmf,frh=0,epv=0,SLP=0,surf_P=0,surf_AT=0,surf_GH=0,surf_RH=0,surf_H2ODMF=0):
	"""
	Creates a GGG-format .mod file
	INPUTS:
	    site_lat    The latitude of the site
	    lev_AT      The pressure levels on which the data are tabulated
	    sat         site Atmospheric Temperature profile (vector)
	    sgh         site Geometric Height profile in km (vector)	# this is from the idl code; but geopotential height is taken from ncep, not geometric
	    site_TP        site Tropopause Pressure (scalar)
	    h2o_dmf       site H2O Dry Mole Fraction (VMR) profile (vector)
	"""

	# Define US Standard Atmosphere (USSA) for use above 10 mbar
	p_ussa=[10.0,  5.0,   2.0,   1.0,   0.1,   0.01,  0.001, 0.0001]
	t_ussa=[227.7, 239.2, 257.9, 270.6, 231.6, 198.0, 189.8, 235.0]
	z_ussa=[31.1,  36.8,  42.4,  47.8,  64.9,  79.3,  92.0,  106.3]

	if type(frh)==int: # ncep mode

		# The head of the .mod file	
		fmt = '{:8.3f} {:11.4e} {:7.3f} {:5.3f} {:8.3f} {:8.3f} {:8.3f}\n'
		mod_content = []
		mod_content+=[	'5  6\n',
						fmt.format(6378.137,6.000E-05,site_lat,9.81,sgh[0],1013.25,site_TP),
						version+'\n',
						' mbar        Kelvin         km      g/mole      DMF       %\n',
						'Pressure  Temperature     Height     MMW        H2O      RH\n',	]

		fmt = '{:9.3e}    {:7.3f}    {:7.3f}    {:7.4f}    {:9.3e}{:>6.1f}\n' # format for writting the lines

		# Export the Pressure, Temp and SHum for lower levels (1000 to 300 mbar)
		for k,elem in enumerate(h2o_dmf):
			svp = svp_wv_over_ice(sat[k])
			h2o_wmf = h2o_dmf[k]/(1+h2o_dmf[k]) # wet mole fraction of h2o
			frh = h2o_wmf*lev_AT[k]/svp # Fractional relative humidity

			# Relace H2O mole fractions that are too small
			if (frh < 30./lev_AT[k]):
				print 'Replacing too-small H2O ',mod_path, lev_AT[k],h2o_wmf,svp*30./lev_AT[k]/lev_AT[k],frh,30./lev_AT[k]				
				frh = 30./lev_AT[k]
				h2o_wmf = svp*frh/lev_AT[k]
				h2o_dmf[k] = h2o_wmf/(1-h2o_wmf)

			# Relace H2O mole fractions that are too large (super-saturated)  GCT 2015-08-05
			if (frh > 1.0):
				print 'Replacing too-large H2O ',mod_path,lev_AT[k],h2o_wmf,svp/lev_AT[k],frh,1.0
				frh=1.0
				h2o_wmf = svp*frh/lev_AT[k]
				h2o_dmf[k] = h2o_wmf/(1-h2o_wmf)

			mmw = 28.964*(1-h2o_wmf)+18.02*h2o_wmf

			mod_content += [fmt.format(lev_AT[k], sat[k],sgh[k],mmw,h2o_dmf[k],100*frh)]

		# Export Pressure and Temp for middle levels (250 to 10 mbar)
		# which have no SHum reanalysis
		ptop = lev_AT[k] # Top pressure level
		frh_top = frh  # remember the FRH at the top (300 mbar) level

		for k in range(len(h2o_dmf),len(lev_AT)): 
			zz = np.log10(lev_AT[k])  # log10[pressure]
			strat_wmf = 7.5E-06*np.exp(-0.16*zz**2)
			svp = svp_wv_over_ice(sat[k])
			trop_wmf = frh_top*svp/lev_AT[k]
			wt = (lev_AT[k]/ptop)**3
			avg_wmf = trop_wmf*wt + strat_wmf*(1-wt)
			avg_frh = avg_wmf*lev_AT[k]/svp
			if (avg_frh > 1.0):
				print 'Replacing super-saturated H2O ',mod_path, lev_AT[k],avg_wmf,svp*avg_frh/lev_AT[k],avg_frh,1.0
				avg_frh = 1.0
				avg_wmf = svp*avg_frh/lev_AT[k]

			mmw = 28.964*(1-avg_wmf)+18.02*avg_wmf
			
			mod_content += [fmt.format(lev_AT[k],sat[k],sgh[k],mmw,avg_wmf/(1-avg_wmf),100*avg_frh)]

		# Get the difference between the USSA and given site temperature at 10 mbar,
		Delta_T=sat[16]-t_ussa[0]

		# Export the P-T profile above 10mbar
		for k in range(1,len(t_ussa)):
			Delta_T=Delta_T/2
			zz = np.log10(p_ussa[k])  # log10[pressure]
			strat_wmf = 7.5E-06*np.exp(-0.16*zz**2)
			svp = svp_wv_over_ice(sat[k])
			mmw = 28.964*(1-strat_wmf)+18.02*strat_wmf
			mod_content += [fmt.format(p_ussa[k],t_ussa[k]+Delta_T,z_ussa[k],mmw,strat_wmf,100*strat_wmf*lev_AT[k]/svp)] # why not use p_ussa instead of lev_AT[k] ?

	else: # merra/geos mode

		# The head of the .mod file	
		fmt1 = '{:8.3f} {:11.4e} {:7.3f} {:5.3f} {:8.3f} {:8.3f} {:8.3f}\n'
		fmt2 = '{:9.3e}    {:9.3e}    {:7.3f}    {:7.3f}    {:9.3e}{:>6.1f}  {:9.3e}\n'
		mod_content = []
		mod_content+=[	'7  7\n',
						fmt1.format(6378.137,6.000E-05,site_lat,9.81,sgh[0],1013.25,site_TP),
						'   SLP          SP           ST         SGH        DMF       RH       TROPPB\n',
						fmt2.format(SLP,surf_P,surf_AT,surf_GH,surf_H2ODMF,surf_RH,site_TP),
						version+'\n',
						' mbar        Kelvin         km      g/mole      DMF       %       k.m+2/kg/s\n',
						'Pressure  Temperature     Height     MMW        H2O      RH          EPV\n',	]

		fmt = '{:9.3e}    {:7.3f}    {:7.3f}    {:7.4f}    {:9.3e}{:>6.1f}    {:9.3e}\n' # format for writting the lines

		# not sure if merra needs all the filters/corrections used for ncep data?

		# Export the Pressure, Temp and SHum
		for k,elem in enumerate(h2o_dmf):
			svp = svp_wv_over_ice(sat[k])
			h2o_wmf = h2o_dmf[k]/(1+h2o_dmf[k]) # wet mole fraction of h2o

			# Relace H2O mole fractions that are too large (super-saturated)  GCT 2015-08-05
			if (frh[k] > 1.0):
				print 'Replacing too-large H2O ',mod_path,lev_AT[k],h2o_wmf,svp/lev_AT[k],frh[k],1.0
				frh[k] = 1.0
				h2o_wmf = svp*frh[k]/lev_AT[k]
				h2o_dmf[k] = h2o_wmf/(1-h2o_wmf)

			mmw = 28.964*(1-h2o_wmf)+18.02*h2o_wmf

			mod_content += [fmt.format(lev_AT[k], sat[k],sgh[k],mmw,h2o_dmf[k],100*frh[k],epv[k])]

	with open(mod_path,'w') as outfile:
		outfile.writelines(mod_content)

def trilinear_interp(fin, fscale_factor, fadd_offset, site_lon_360, lon_XX, site_lat, lat_XX, site_tim, tim_XX):
	"""
	Evaluates  fout = fin(xx,yy,*,tt) 
	Result is a 1-vector
	"""

	dx = lon_XX[1]-lon_XX[0]
	dy = lat_XX[1]-lat_XX[0]
	dt = tim_XX[1]-tim_XX[0]

	xx = (site_lon_360-lon_XX[0])/dx
	yy = (site_lat-lat_XX[0])/dy
	tt = (site_tim-tim_XX[0])/dt

	nxx =  len(lon_XX)
	nyy =  len(lat_XX)
	ntt =  len(tim_XX)

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
	   print 'An NCEP file doesnt cover the full range of dates'
	   print 'site_tim',site_tim
	   print 'tim_XX',tim_XX
	   raw_input() # will hold the program until something is typed in commandline

	if (fr_xx < 0) or (fr_xx > 1):
	   print 'Excessive longitude extrapolation:',fr_xx,' steps   =',fr_xx*dx,' deg'
	   print ' xx= ',xx,'  index_xx=',index_xx,'  fr_xx=',fr_xx
	   print 'NCEP file doesnt cover the full range of longitudes'
	   raw_input() # will hold the program until something is typed in commandline

	if (fr_yy < 0) or (fr_yy > 1):
	   print 'Excessive latitude extrapolation:',fr_yy-1,' steps   =',(fr_yy-1)*dy,' deg'
	   print ' yy= ',yy,'  index_yy=',index_yy,'  fr_yy=',fr_yy
	   print 'NCEP file doesnt cover the full range of latitudes'
	   raw_input() # will hold the program until something is typed in commandline

	if (fr_tt < 0) or (fr_tt > 1):
		print ' Warning: time extrapolation of ',fr_tt,' time-steps'
	if (fr_xx < 0) or (fr_xx > 1):
		print ' Warning: longitude extrapolation of ',fr_xx,' steps'
	if (fr_yy < 0) or (fr_yy > 1):
		print ' Warning: latitude extrapolation of ',fr_yy,' steps'
	
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

	fout = fout*fscale_factor + fadd_offset

	return fout

def read_data(dataset, varname, min_lat_ID='', max_lat_ID='', min_lon_ID='',max_lon_ID=''):
	"""
	for ncep files "dataset" is the full path to the netcdf file

	for merra files "dataset" is a pydap.model.DatasetType object
	"""
	opendap = type(dataset)!=netCDF4._netCDF4.Dataset

	# Initialize add_offset and scale_factor
	data_add_offset_XX = 0
	data_scale_factor_XX = 1.0

	if min_lat_ID == '': # ncep mode
		try: 
			dataset = netCDF4.Dataset(dataset,'r')	# open the netcdf files
		except:
			print dataset,'not found'
			sys.exit()

		lev_XX = dataset.variables['level'][:]	# Read in variable 'level'
		lat_XX = dataset.variables['lat'][:]	# Read in variable 'lat'
		lon_XX = dataset.variables['lon'][:]	# Read in variable 'lon'
		tim_XX = dataset.variables['time'][:]	# Read in variable 'time'
		data_XX = dataset.variables[varname][:]	# Read in variable varname

		# Check variable for offset and scaling factor attributes
		for attribute in dataset.variables[varname].ncattrs():
			# Set data_add_offset_XX
			if attribute == 'add_offset':
				data_add_offset_XX = dataset.variables[varname].getncattr(attribute)
			# Set data_scale_factor_XX
			if attribute == 'scale_factor':
				data_scale_factor_XX = dataset.variables[varname].getncattr(attribute)

		time_units = dataset.variables['time'].units  # string containing definition of time units
			
		dataset.close()
	else: # merra/fp mode

		if dataset[varname].ndim == 4:
			if dataset['lev'].shape[0] == 72:
				# the merra 72 mid level pressures are not fixed
				lev_XX = dataset['PL'][:,:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID]
			elif dataset['lev'].shape[0] == 42:
				# the 42 levels data is on a fixed pressure grid
				lev_XX = dataset['lev'][:]
		else: # surface data doesn't have a 'lev' variable
			lev_XX = dataset['PS'][:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID]
			

		lat_XX = dataset['lat'][min_lat_ID:max_lat_ID] 	# Read in variable 'lat'
		lon_XX = dataset['lon'][min_lon_ID:max_lon_ID] 	# Read in variable 'lon'
		tim_XX = dataset['time'][:]	# Read in variable 'time'
		if opendap:
			tim_XX = tim_XX.data
			lev_XX = lev_XX.data
			lat_XX = lat_XX.data
			lon_XX = lon_XX.data
		# get longitudes as 0 -> 360 instead of -180 -> 180, needed for trilinear_interp
		for i,elem in enumerate(lon_XX):
			if elem < 0:
				lon_XX[i] = elem + 360.0
		#sort_IDs = np.argsort(lon_XX)
		#lon_XX = lon_XX[sort_IDs]

		if dataset[varname].ndim == 4:
			data_XX = dataset[varname][:,:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID] 	# Read in variable varname
		else:
			data_XX = dataset[varname][:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID] 	# Read in variable varname
		if opendap or ('Masked' in str(type(data_XX))):
			data_XX = data_XX.data

		time_units = dataset['time'].units  # string containing definition of time units

	print varname,': lev_shape,lat_shape,lon_shape,tim_shape =',lev_XX.shape,lat_XX.shape,lon_XX.shape,tim_XX.shape

	# two lines to parse the date (no longer need to worry about before/after 2014)
	date_list = re.findall(r"[\w]+",time_units.split(' since ')[1])
	common_date_format = '{:0>4}-{:0>2}-{:0>2} {:0>2}:{:0>2}:{:0>2}'.format(*date_list)

	start_date = datetime.strptime(common_date_format,'%Y-%m-%d %H:%M:%S')
	astropy_start_date = Time(start_date)

	julday0 = astropy_start_date.jd # gives same results as IDL's JULDAY function

	return lev_XX, lat_XX, lon_XX, tim_XX, data_XX, data_scale_factor_XX, data_add_offset_XX, julday0

def querry_indices(dataset,site_lat,site_lon_180,box_lat_half_width,box_lon_half_width):
	"""	
	Set up a lat-lon box for the data querry

	Unlike with ncep, this will only use daily files for interpolation, so no time box is defined
	
	NOTE: merra lat -90 -> +90 ;  merra lon -180 -> +179.375

	To be certain to get two points on both side of the site lat and lon, use the grid resolution
	"""
	# define 2xbox_lat_half_width°x2xbox_lon_half_width° lat-lon box centered on the site lat-lon and construct the url querry
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

	return min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID

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

	# hour and minute for time interpolation, default is noon
	if len(argu)>4:
		HH = int(argu[4])
		MM = int(argu[5])
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
	print 'Local time for interpolation: {:0>2}:{:0>2}'.format(HH,MM)

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

	if mode == 'ncep':

		# path to the netcdf files
		ncdf_AT_file = os.path.join(ncdf_path,'.'.join(['air','{:0>4}'.format(start_date.year),'nc']))
		ncdf_GH_file = os.path.join(ncdf_path,'.'.join(['hgt','{:0>4}'.format(start_date.year),'nc']))
		ncdf_SH_file = os.path.join(ncdf_path,'.'.join(['shum','{:0>4}'.format(start_date.year),'nc']))

		print 'Read global',start_date.year,'NCEP data ...'
		# Air Temperature
		lev_AT,lat_AT, lon_AT, tim_AT, data_AT, data_scale_factor_AT, data_add_offset_AT, julday0 = read_data(ncdf_AT_file, 'air')
		if len(lev_AT) < 17:
			print 'Need 17 levels of AT data: found only ',len(lev_AT)
		
		# Geopotential Height
		lev_GH,lat_GH, lon_GH, tim_GH, data_GH, data_scale_factor_GH, data_add_offset_GH, julday0 = read_data(ncdf_GH_file, 'hgt')
		if len(lev_GH) < 17:
			print 'Need 17 levels of GH data: found only ',len(lev_GH)
		
		# Specific Humidity
		lev_SH,lat_SH, lon_SH, tim_SH, data_SH, data_scale_factor_SH, data_add_offset_SH, julday0 = read_data(ncdf_SH_file, 'shum')
		if len(lev_SH) <  8:
			print 'Need  8 levels of SH data: found only ',len(lev_SH)
	elif 'glob' in mode:
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
		min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID = querry_indices(dataset,site_lat,site_lon_180,box_lat_half_width,box_lon_half_width)

		# multi-level data
		print 'Read global',start_date.year,mode,'data ...'
		# Air temperature
		lev_AT,lat_AT, lon_AT, tim_AT, data_AT, data_scale_factor_AT, data_add_offset_AT, julday0 = read_data(dataset,'T',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
		# Specific humidity
		lev_SH,lat_SH, lon_SH, tim_SH, data_SH, data_scale_factor_SH, data_add_offset_SH, julday0 = read_data(dataset,'QV',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
		# Relative humidity
		lev_RH,lat_RH, lon_RH, tim_RH, data_RH, data_scale_factor_RH, data_add_offset_RH, julday0 = read_data(dataset,'RH',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
		# Height
		lev_H,lat_H, lon_H, tim_H, data_H, data_scale_factor_H, data_add_offset_H, julday0 = read_data(dataset,'H',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
		# Potential vorticity
		lev_EPV,lat_EPV, lon_EPV, tim_EPV, data_EPV, data_scale_factor_EPV, data_add_offset_EPV, julday0 = read_data(dataset,'EPV',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)

		# single level data
		print 'Read global',start_date.year,mode,'surface data ...'
		# 2 meter Air Temperature
		lev_surf_AT,lat_surf_AT, lon_surf_AT, tim_surf_AT, data_surf_AT, data_scale_factor_surf_AT, data_add_offset_surf_AT, julday0 = read_data(surface_dataset,'T2M',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)		
		# 2 meter Specific humidity
		lev_surf_SH,lat_surf_SH, lon_surf_SH, tim_surf_SH, data_surf_SH, data_scale_factor_surf_SH, data_add_offset_surf_SH, julday0 = read_data(surface_dataset,'QV2M',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)		
		# Surface pressure
		lev_surf_P,lat_surf_P, lon_surf_P, tim_surf_P, data_surf_P, data_scale_factor_surf_P, data_add_offset_surf_P, julday0 = read_data(surface_dataset,'PS',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
		# Sea level pressure
		lev_SLP,lat_SLP, lon_SLP, tim_SLP, data_SLP, data_scale_factor_SLP, data_add_offset_SLP, julday0 = read_data(surface_dataset,'SLP',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)								
		# surface geopotential height
		lev_surf_GH,lat_surf_GH, lon_surf_GH, tim_surf_GH, data_surf_GH, data_scale_factor_surf_GH, data_add_offset_surf_GH, julday0 = read_data(dataset,'PHIS',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID) # surface geopotential height is in the levels datasets			
		data_surf_GH = data_surf_GH / gravity_at_lat # convert from m2 s-2 to m
		# Tropopause pressure (blended)
		lev_TP,lat_TP, lon_TP, tim_TP, data_TP, data_scale_factor_TP, data_add_offset_TP, julday0 = read_data(surface_dataset,'TROPPB',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)								
		
		# merra/geos time is minutes since base time, need to convert to hours
		tim_AT = tim_AT / 60.0
		tim_SH = tim_SH / 60.0
		tim_RH = tim_RH / 60.0
		tim_H = tim_H / 60.0
		tim_EPV = tim_EPV / 60.0
		tim_surf_AT = tim_surf_AT / 60.0
		tim_surf_SH = tim_surf_SH / 60.0
		tim_surf_GH = tim_surf_GH / 60.0
		tim_surf_P = tim_surf_P / 60.0
		tim_SLP = tim_SLP / 60.0
		tim_TP = tim_TP / 60.0

	date = start_date + timedelta(hours=HH,minutes=MM) # date with local time
	astropy_date = Time(date)

	time_step = timedelta(days=1) # time step between mod files; will need to change the mod file naming and gsetup to do sub-daily files
	print 'Time step:',time_step.total_seconds()/3600.0,'hours'
	new_year = False
	while date<end_date:

		# use the local date for the name of the .mod file
		YYYYMMDD = date.strftime('%Y%m%d')
		mod_name = '_'.join(['NCEP',YYYYMMDD,'{:0>2.0f}'.format(round(abs(site_lat)))+ns,'{:0>3.0f}'.format(round(abs(site_lon_180)))+ew+'.mod'])
		mod_file_path = os.path.join(mod_path,mod_name)
		print '\n',mod_name

		if 'merradap' in mode: # MERRA with openDAP

			# in the url below, M2I3NVASM is for every 3 hour instantaneous assimilated meteorological fields; VASM is for 72 levels; PASM is for 42 levels
			if 'merradap42' in mode:
				letter = 'P'
			elif 'merradap72' in mode:
				letter = 'V'

			UTC_date = date + timedelta(hours = -site_lon_180/15.0) # merra times are in UTC, so the date may be different than the local date, make sure to use the UTC date to querry the file

			# multi levels data
			url = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2I3N{}ASM.5.12.4/{:0>4}/{:0>2}/MERRA2_400.inst3_3d_asm_N{}.{:0>4}{:0>2}{:0>2}.nc4'.format(letter,UTC_date.year,UTC_date.month,letter.lower(),UTC_date.year,UTC_date.month,UTC_date.day)
			session = setup_session(username,password,check_url=url)
			try:
				dataset = open_url(url,session=session)
			except HTTPError:
				print 'url not valid\n',url
				sys.exit()

			# single level data, from M2I1NXASM ; it is on the goldsmr4 server and not goldsmr5 like above
			surface_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2I1NXASM.5.12.4/{:0>4}/{:0>2}/MERRA2_400.inst1_2d_asm_Nx.{:0>4}{:0>2}{:0>2}.nc4'.format(UTC_date.year,UTC_date.month,UTC_date.year,UTC_date.month,UTC_date.day)
			surface_session = setup_session(username,password,check_url=surface_url)
			try:
				surface_dataset = open_url(surface_url,session=surface_session)
			except HTTPError:
				print 'url not valid\n',url
				sys.exit()				

			# get the min/max lat-lon indices of merra lat-lon that lies within a lat-lon box centered on the site lat-lon.
			min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID = querry_indices(dataset,site_lat,site_lon_180,2.5,2.5)

			# Air temperature
			lev_AT,lat_AT, lon_AT, tim_AT, data_AT, data_scale_factor_AT, data_add_offset_AT, julday0 = read_data(dataset,'T',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Specific humidity
			lev_SH,lat_SH, lon_SH, tim_SH, data_SH, data_scale_factor_SH, data_add_offset_SH, julday0 = read_data(dataset,'QV',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Relative humidity
			lev_RH,lat_RH, lon_RH, tim_RH, data_RH, data_scale_factor_RH, data_add_offset_RH, julday0 = read_data(dataset,'RH',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Height
			lev_H,lat_H, lon_H, tim_H, data_H, data_scale_factor_H, data_add_offset_H, julday0 = read_data(dataset,'H',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Potential vorticity
			lev_EPV,lat_EPV, lon_EPV, tim_EPV, data_EPV, data_scale_factor_EPV, data_add_offset_EPV, julday0 = read_data(dataset,'EPV',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			

			# 2 meter Air Temperature
			lev_surf_AT,lat_surf_AT, lon_surf_AT, tim_surf_AT, data_surf_AT, data_scale_factor_surf_AT, data_add_offset_surf_AT, julday0 = read_data(surface_dataset,'T2M',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)		
			# 2 meter Specific humidity
			lev_surf_SH,lat_surf_SH, lon_surf_SH, tim_surf_SH, data_surf_SH, data_scale_factor_surf_SH, data_add_offset_surf_SH, julday0 = read_data(surface_dataset,'QV2M',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)		
			# Surface pressure
			lev_surf_P,lat_surf_P, lon_surf_P, tim_surf_P, data_surf_P, data_scale_factor_surf_P, data_add_offset_surf_P, julday0 = read_data(surface_dataset,'PS',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Sea level pressure
			lev_SLP,lat_SLP, lon_SLP, tim_SLP, data_SLP, data_scale_factor_SLP, data_add_offset_SLP, julday0 = read_data(surface_dataset,'SLP',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)								
			# surface geopotential height
			lev_surf_GH,lat_surf_GH, lon_surf_GH, tim_surf_GH, data_surf_GH, data_scale_factor_surf_GH, data_add_offset_surf_GH, julday0 = read_data(dataset,'PHIS',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID) # surface geopotential height is in the levels datasets			
			data_surf_GH = data_surf_GH / gravity_at_lat # convert from m2 s-2 to m
			# Tropopause pressure (blended)
			lev_TP,lat_TP, lon_TP, tim_TP, data_TP, data_scale_factor_TP, data_add_offset_TP, julday0 = read_data(surface_dataset,'TROPPB',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)								
			
			# merra/geos time is minutes since base time, need to convert to hours
			tim_AT = tim_AT / 60.0
			tim_SH = tim_SH / 60.0
			tim_RH = tim_RH / 60.0
			tim_H = tim_H / 60.0
			tim_EPV = tim_EPV / 60.0
			tim_surf_AT = tim_surf_AT / 60.0
			tim_surf_SH = tim_surf_SH / 60.0
			tim_surf_GH = tim_surf_GH / 60.0
			tim_surf_P = tim_surf_P / 60.0
			tim_SLP = tim_SLP / 60.0
			tim_TP = tim_TP / 60.0

		"""
		Interpolation time:
			julday0 is the fractional julian day number of the base time of the dataset: dataset times are in UTC hours since base time
			astropy_date.jd is the fractional julian day number of the current local day
			(astropy_date.jd-julday0)*24.0 = local hours since julday0
		"""
		site_tim = (astropy_date.jd-julday0)*24.0 - site_lon_180/15.0 # UTC hours since julday0
		
		# interpolate the data to the site's location and the desired time
		site_AT = trilinear_interp(data_AT, data_scale_factor_AT, data_add_offset_AT, site_lon_360, lon_AT, site_lat, lat_AT, site_tim, tim_AT) 
		site_SH = trilinear_interp(data_SH, data_scale_factor_SH, data_add_offset_SH, site_lon_360, lon_SH, site_lat, lat_SH, site_tim, tim_SH)

		if 'ncep' in mode:
			# ncep has geopotential height profiles, merra has heights
			site_GH = trilinear_interp(data_GH, data_scale_factor_GH, data_add_offset_GH, site_lon_360, lon_GH, site_lat, lat_GH, site_tim, tim_GH)
			site_TP = 0.0 # tropopause pressure not used with NCEP data
			site_RH = 0 # won't be used, just to feed something to write_mod frh
		else:
			# interpolate height levels
			site_H = trilinear_interp(data_H, data_scale_factor_H, data_add_offset_H, site_lon_360, lon_H, site_lat, lat_H, site_tim, tim_H)
			# convert to geopotential
			site_GH = (site_H)/(1.0+(site_H)/earth_radius_at_lat)		# Convert from geometric to geopotential
			# interpolate relative humidity
			site_RH = trilinear_interp(data_RH, data_scale_factor_RH, data_add_offset_RH, site_lon_360, lon_RH, site_lat, lat_RH, site_tim, tim_RH)
			# interpolate potential vorticity
			site_EPV = trilinear_interp(data_EPV, data_scale_factor_EPV, data_add_offset_EPV, site_lon_360, lon_EPV, site_lat, lat_EPV, site_tim, tim_EPV)

		if 'merradap72' in mode:
			# the merra 72 mid level pressures are not fixed, so need to interpolate to get just 1 array of levels
			# Mid-level pressure
			lev_PL,lat_PL, lon_PL, tim_PL, data_PL, data_scale_factor_PL, data_add_offset_PL, julday0 = read_data(dataset,'PL',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)			
			tim_PL = tim_PL/60.0 # convert minutes to hours
			# interpolate pressure levels
			lev_AT = trilinear_interp(data_PL, data_scale_factor_PL, data_add_offset_PL, site_lon_360, lon_PL, site_lat, lat_PL, site_tim, tim_PL) 
			lev_AT = lev_AT / 100.0 # convert from Pa to hPa

		if 'ncep' not in mode:
			# get rid of fill values
			without_fill_IDs = np.where(site_AT<1e10) # merra/geos fill value is 1e15

			site_AT = site_AT[without_fill_IDs]
			site_GH = site_GH[without_fill_IDs]
			site_SH = site_SH[without_fill_IDs]
			site_RH = site_RH[without_fill_IDs]
			site_EPV = site_EPV[without_fill_IDs]
			if 'glob' in mode:
				lev_AT = lev_SH[without_fill_IDs] # I use lev_SH because lev_AT is not redefined in the while loop contrary to the others, otherwise it would trigger index errors
			else:
				lev_AT = lev_AT[without_fill_IDs]

		if ('merradap' in mode) or ('glob' in mode):
			# interpolate the surface data
			site_surf_AT = trilinear_interp(data_surf_AT, data_scale_factor_surf_AT, data_add_offset_surf_AT, site_lon_360, lon_surf_AT, site_lat, lat_surf_AT, site_tim, tim_surf_AT) 
			site_surf_GH = trilinear_interp(data_surf_GH, data_scale_factor_surf_GH, data_add_offset_surf_GH, site_lon_360, lon_surf_GH, site_lat, lat_surf_GH, site_tim, tim_surf_GH)
			site_surf_SH = trilinear_interp(data_surf_SH, data_scale_factor_surf_SH, data_add_offset_surf_SH, site_lon_360, lon_surf_SH, site_lat, lat_surf_SH, site_tim, tim_surf_SH)
			site_surf_P = trilinear_interp(data_surf_P, data_scale_factor_surf_P, data_add_offset_surf_P, site_lon_360, lon_surf_P, site_lat, lat_surf_P, site_tim, tim_surf_P)
			site_surf_P = site_surf_P / 100.0 # convert Pa to hPa
			site_SLP = trilinear_interp(data_SLP, data_scale_factor_SLP, data_add_offset_SLP, site_lon_360, lon_SLP, site_lat, lat_SLP, site_tim, tim_SLP)
			site_SLP = site_SLP / 100.0 # convert Pa to hPa
			site_TP = trilinear_interp(data_TP, data_scale_factor_TP, data_add_offset_TP, site_lon_360, lon_TP, site_lat, lat_TP, site_tim, tim_TP)
			site_TP = site_TP / 100.0 # convert Pa to hPa

			# site_TP = site_TP/100.0 # convert from Pa to hPa
			# I think the TROPPB is wrongly indicated as having 'Pa' units: the values are of order 1-3 *10e4, if we assume they are hPa this is usually between 8-16 km
			# if we divide those by 100 we get values of order 1-3 *10e2 which is between 40-45 km !!

			if 'merradap72' in mode: # merra42 and ncep go from high pressure to low pressure, but merra 72 does the reverse
				# reverse merra72 profiles
				site_AT = site_AT[::-1]
				site_GH = site_GH[::-1]
				site_SH = site_SH[::-1]
				site_RH = site_RH[::-1]
				lev_AT = lev_AT[::-1]

		site_H2ODMF = rmm*site_SH/(1-site_SH) # Convert specific humidity, a wet mass mixing ratio, to dry mole fraction
		site_GH = site_GH/1000.0	# Convert m to km

		if ('merradap' in mode) or ('glob' in mode):
			site_surf_H2ODMF = rmm*site_surf_SH/(1-site_surf_SH)
			site_surf_GH = site_surf_GH/1000.0
			# compute surface relative humidity
			svp = svp_wv_over_ice(site_surf_AT)
			h2o_wmf = site_surf_H2ODMF/(1+site_surf_H2ODMF) # wet mole fraction of h2o
			site_surf_RH = h2o_wmf*site_surf_P/svp # Fractional relative humidity

		# write the .mod file
		version = 'mod_maker_10.6   2017-04-11   GCT'
		if 'ncep' in mode:
			write_mod(mod_file_path,version,site_lat,lev_AT,site_AT,site_GH,site_TP,site_H2ODMF,frh=site_RH)
		else:
			write_mod(mod_file_path,version,site_lat,lev_AT,site_AT,site_GH,site_TP,site_H2ODMF,frh=site_RH,epv=site_EPV,SLP=site_SLP,surf_P=site_surf_P,surf_AT=site_surf_AT,surf_GH=site_surf_GH,surf_RH=site_surf_RH,surf_H2ODMF=site_surf_H2ODMF)

		if ((date+time_step).year!=date.year):
			new_year = True

		date = date + time_step
		astropy_date = Time(date)
		
		if ('ncep' in mode) and new_year:
			# path to the netcdf files
			ncdf_AT_file = os.path.join(ncdf_path,'.'.join(['air','{:0>4}'.format(date.year),'nc']))
			ncdf_GH_file = os.path.join(ncdf_path,'.'.join(['hgt','{:0>4}'.format(date.year),'nc']))
			ncdf_SH_file = os.path.join(ncdf_path,'.'.join(['shum','{:0>4}'.format(date.year),'nc']))

			print '\nRead global',date.year,'ncep data ...'
			# Air Temperature
			lev_AT,lat_AT, lon_AT, tim_AT, data_AT, data_scale_factor_AT, data_add_offset_AT, julday0 = read_data(ncdf_AT_file, 'air')
			if len(lev_AT) < 17:
				print 'Need 17 levels of AT data: found only ',len(lev_AT)
			
			# Geopotential Height
			lev_GH,lat_GH, lon_GH, tim_GH, data_GH, data_scale_factor_GH, data_add_offset_GH, julday0 = read_data(ncdf_GH_file, 'hgt')
			if len(lev_GH) < 17:
				print 'Need 17 levels of GH data: found only ',len(lev_GH)
			
			# Specific Humidity
			lev_SH,lat_SH, lon_SH, tim_SH, data_SH, data_scale_factor_SH, data_add_offset_SH, julday0 = read_data(ncdf_SH_file, 'shum')
			if len(lev_SH) <  8:
				print 'Need  8 levels of SH data: found only ',len(lev_SH)

			new_year = False
