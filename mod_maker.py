#!~/anaconda2/bin/python
 # -*- coding: utf-8 -*-

"""
mod_maker10.f translated into python

What's different?

This is run with:

python mod_maker.py arg1 arg2 ar3 arg4

arg1: two letter site abbreviation (e.g. 'oc' for Lamont, Oklahoma; see the "site_dict" dictionary)
arg2: year (YYYY)
arg3: mode ('ncep', 'merra42,username,password', or 'merra72,username,password')
arg4: (optional, default=12) hour (HH) for time interpolation
arg2: (optional, default=0) minute (MM) for time interpolation

In GGGPATH/models/gnd it will write one mod file per day.
The merra modes require an internet connection and EarthData credentials
The ncep mode requires the global NCEP netcdf files of the given year to be present in GGGPATH/ncdf

There is dictionary of sites with their respective lat/lon, so this works for all TCCON sites, lat/lon values were taken from the wiki page of each site.

note: still need to make it work for Darwin, or any site that changed location at some point and must use different lat/lon for different periods
"""

import numpy as np
import os, sys
import netCDF4
from datetime import datetime, timedelta
import time
import re
import jdcal

from pydap.cas.urs import setup_session
from pydap.client import open_url

def svp_wv_over_ice(temp):
	"""	
	Uses the Goff-Gratch equation to calculate the saturation vapor
	pressure of water vapor over ice at a user-specified temperature.
		Input:  temp (K)
		Output: svp (mbar)
	"""
	t0 = 273.16
	tr = t0/temp # triple point temperature
	yy = -9.09718*(tr-1)-3.56654*np.log10(tr)+0.876793*(1-1/tr)
	svp = 6.1173*10**yy # saturation vapor pressure over ice

	return svp

def write_mod(mod_path,site_lat,lev_AT,sat,sgh,sntp,h2omf):
	"""
	Creates a GGG-format .mod file
	INPUTS:
	    site_lat    The latitude of the site
	    lev_AT      The pressure levels on which the data are tabulated
	    sat         site Noon Atmospheric Temperature profile (vector)
	    sgh         site Noon Geometric Height profile in km (vector)
	    sntp        site Noon Tropopause Pressure (scalar)
	    h2omf       site Noon H2O Mole Fraction (VMR) profile (vector)
	"""

	# Define US Standard Atmosphere (USSA) for use above 10 mbar
	p_ussa=[10.0,  5.0,   2.0,   1.0,   0.1,   0.01,  0.001, 0.0001]
	t_ussa=[227.7, 239.2, 257.9, 270.6, 231.6, 198.0, 189.8, 235.0]
	z_ussa=[31.1,  36.8,  42.4,  47.8,  64.9,  79.3,  92.0,  106.3]

	# The head of the .mod file	
	fmt = '{:8.3f} {:11.4e} {:7.3f} {:5.3f} {:8.3f} {:8.3f} {:8.3f}'
	mod_content = []
	mod_content+=[	'4  6',
					fmt.format(6378.137,6.000E-05,site_lat,9.81,sgh[0],1013.25,sntp/10**2),
					' mbar        Kelvin         km      g/mole      vmr       %',
					'Pressure  Temperature     Height     MMW        H2O      RH',	]

	fmt = '{:9.3e}    {:7.3f}    {:7.3f}    {:7.4f}    {:9.3e}{:>6.1f}'

	# Export the Pressure, Temp and SHum for lower levels (1000 to 300 mbar)
	for k,elem in enumerate(h2omf):
		svp = svp_wv_over_ice(sat[k])
		satvmr = svp/lev_AT[k]   # vmr of saturated WV at T/P
		
		frh = h2omf[k]/satvmr    # fractional RH

		# Relace H2O mole fractions that are too small
		if (frh < 30./lev_AT[k]):
			frh = 30./lev_AT[k]
			print 'Replacing too small H2O_MF ',lev_AT[k],h2omf[k],frh*satvmr,frh
			h2omf[k] = frh*satvmr

		# Relace H2O mole fractions that are too large (super-saturated)  GCT 2015-08-05
		if (frh > 1.0):
			frh=1.0
			print 'Replacing too large H2O_MF ',lev_AT[k],h2omf[k],frh*satvmr,frh
			h2omf[k] = frh*satvmr

		mod_content+=[fmt.format(lev_AT[k], sat[k],sgh[k],28.9640,h2omf[k],100*frh)]

	# Export Pressure and Temp for middle levels (250 to 10 mbar)
	# which have no SHum reanalysis.
	ptop = lev_AT[k] # Top pressure level covered by NCEP H2O
	frhtop = frh  # remember the FRH at the top (300 mbar) level

	for k in range(len(h2omf),len(lev_AT)):
		zz = np.log10(lev_AT[k])  # log10[pressure]
		strat_h2o = 7.5E-06*np.exp(-0.16*zz**2)
		svp = svp_wv_over_ice(sat[k])
		satvmr = svp/lev_AT[k]   # vmr of saturated WV at T/P
		trop_h2o = satvmr*frhtop
		wt = (lev_AT[k]/ptop)**3
		avg_h2o = trop_h2o*wt + strat_h2o*(1-wt)
		if (avg_h2o > satvmr):
		  print 'Replacing super-saturated H2O_MF ',lev_AT[k],avg_h2o,satvmr
		  avg_h2o = satvmr
		
		mod_content += [fmt.format(lev_AT[k],sat[k],sgh[k],28.9640,avg_h2o,100*avg_h2o/satvmr)]

	# Get the difference between the USSA and given site temperature at 10 mbar,
	Delta_T=sat[16]-t_ussa[0]

	# Export the P-T profile above 10mbar
	for k in range(1,len(t_ussa)):
		Delta_T=Delta_T/2
		zz = np.log10(p_ussa[k])  # log10[pressure]
		strat_h2o = 7.5E-06*np.exp(-0.16*zz**2)
		svp = svp_wv_over_ice(sat[k])
		satvmr = svp/lev_AT[k]   # vmr of saturated WV at T/P
		# print p_ussa[k],strat_h2o
		mod_content += [fmt.format(p_ussa[k],t_ussa[k]+Delta_T,z_ussa[k],28.9640,strat_h2o,100*strat_h2o/satvmr)]

	mod_content = [line+'\n' for line in mod_content]

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

	if(fr_tt < -1 or fr_tt > 2):
	   print 'Excessive time extrapolation:',fr_tt,' time-steps   =',fr_tt*dt,' days'
	   print ' tt= ',tt,'  index_tt=',index_tt,'  fr_tt=',fr_tt
	   print 'An NCEP file doesnt cover the full range of dates'
	   raw_input() # will hold the program until something is typed in commandline

	if(fr_xx < 0 or fr_xx > 1):
	   print 'Excessive longitude extrapolation:',fr_xx,' steps   =',fr_xx*dx,' deg'
	   print ' xx= ',xx,'  index_xx=',index_xx,'  fr_xx=',fr_xx
	   print 'NCEP file doesnt cover the full range of longitudes'
	   raw_input() # will hold the program until something is typed in commandline

	if(fr_yy < 0 or fr_yy > 1):
	   print 'Excessive latitude extrapolation:',fr_yy-1,' steps   =',(fr_yy-1)*dy,' deg'
	   print ' yy= ',yy,'  index_yy=',index_yy,'  fr_yy=',fr_yy
	   print 'NCEP file doesnt cover the full range of latitudes'
	   raw_input() # will hold the program until something is typed in commandline

	if(fr_tt < -0 or fr_tt > 1):
		print ' Warning: time extrapolation of ',fr_tt,' time-steps'
	if(fr_xx < -0 or fr_xx > 1):
		print ' Warning: longitude extrapolation of ',fr_xx,' steps'
	if(fr_yy < -0 or fr_yy > 1):
		print ' Warning: latitude extrapolation of ',fr_yy,' steps'
	
	if fin.ndim==4:
		fout =	((fin[index_tt,:,index_yy,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt,:,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
		+ (fin[index_tt,:,index_yy+1,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt,:,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*(1.0-fr_tt) \
		+ ((fin[index_tt+1,:,index_yy,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt+1,:,index_yy+1,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
		+ (fin[index_tt+1,:,index_yy+1,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt+1,:,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*fr_tt
	else: # surface geopotential from merra data will not have the levels dimension
		fout =	((fin[index_tt,index_yy,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
		+ (fin[index_tt,index_yy+1,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*(1.0-fr_tt) \
		+ ((fin[index_tt+1,index_yy,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt+1,index_yy+1,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
		+ (fin[index_tt+1,index_yy+1,index_xx]*(1.0-fr_xx) \
		+ fin[index_tt+1,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*fr_tt

	fout = fout*fscale_factor + fadd_offset

	return fout

def read_data(dataset, varname, min_lat_ID='', max_lat_ID='', min_lon_ID='',max_lon_ID=''):
	"""
	for ncep files "dataset" is the full path to the netcdf file

	for merra files "dataset" is a pydap.model.DatasetType object
	"""

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
	else: # merra mode

		if dataset[varname].ndim == 4:
			if dataset['lev'].shape[0] == 72:
				# the merra 72 mid level pressures are not fixed
				lev_XX = dataset['PL'][:,:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID].data
			elif dataset['lev'].shape[0] == 42:
				# the 42 levels data is on a fixed pressure grid
				lev_XX = dataset['lev'][:].data
		else: # surface data doesn't have a 'lev' variable
			lev_XX = dataset['PS'][:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID].data

		lat_XX = dataset['lat'][min_lat_ID:max_lat_ID].data 	# Read in variable 'lat'
		lon_XX = dataset['lon'][min_lon_ID:max_lon_ID].data 	# Read in variable 'lon'
		# get longitudes as 0 -> 360 instead of -180 -> 180, needed for trilinear_interp
		for i,elem in enumerate(lon_XX):
			if elem < 0:
				lon_XX[i] = elem + 360.0
		#sort_IDs = np.argsort(lon_XX)
		#lon_XX = lon_XX[sort_IDs]

		tim_XX = dataset['time'][:].data 	# Read in variable 'time'
		if dataset[varname].ndim == 4:
			data_XX = dataset[varname][:,:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID].data 	# Read in variable varname
		else:
			data_XX = dataset[varname][:,min_lat_ID:max_lat_ID,min_lon_ID:max_lon_ID].data 	# Read in variable varname

		time_units = dataset['time'].units  # string containing definition of time units

	print varname,': lev_shape,lat_shape,lon_shape,tim_shape =',lev_XX.shape,lat_XX.shape,lon_XX.shape,tim_XX.shape

	# two lines to parse the date (no longer need to worry about before/after 2014)
	date_list = re.findall(r"[\w]+",time_units.split(' since ')[1])
	common_date_format = '{:0>4}-{:0>2}-{:0>2} {:0>2}:{:0>2}:{:0>2}'.format(*date_list)

	start_date = time.strptime(common_date_format,'%Y-%m-%d %H:%M:%S')

	julday0 = int(round(sum(jdcal.gcal2jd(*start_date[:3])))) # gives same results as IDL's JULDAY function
	
	return lev_XX, lat_XX, lon_XX, tim_XX, data_XX, data_scale_factor_XX, data_add_offset_XX, julday0

def merra_querry_indices(dataset,site_lat,site_lon_180,box_lat_half_width,box_lon_half_width):
	"""	
	Set up a lat-lon box for the data querry

	Unlike with ncep, this will only use daily files for interpolation, so no time box is defined
	
	NOTE: merra lat -90 -> +90 ;  merra lon -180 -> +179.375
	"""
	# define 5°x5° lat-lon box centered on the site lat-lon and construct the url querry
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
	merra_lon = dataset['lon'][:].data
	merra_lat = dataset['lat'][:].data

	# get the indices of merra longitudes and latitudes that fit in the lat-lon box
	merra_lon_in_box_IDs = np.where((merra_lon>=min_lon) & (merra_lon<=max_lon))[0]
	merra_lat_in_box_IDs = np.where((merra_lat>=min_lat) & (merra_lat<=max_lat))[0]

	min_lat_ID, max_lat_ID = merra_lat_in_box_IDs[0], merra_lat_in_box_IDs[-1]+1
	min_lon_ID, max_lon_ID = merra_lon_in_box_IDs[0], merra_lon_in_box_IDs[-1]+1
	# +1 because ARRAY[i:j] in python will return elements i to j-1

	return min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID


def gravity(gdlat,altit):
	"""
	copy/pasted from fortran routine comments
	This is used to convert

	Input Parameters:
	    gdlat       GeoDetric Latitude (degrees)
	    altit       Geometric Altitude (km)
	
	Output Parameter:
	    gravity     Effective Gravitational Acceleration (m/s2)
	
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

if __name__ == "__main__":

	# dictionary mapping TCCON site abbreviations to their lat-lon-alt data, and full names
	site_dict = {
				'pa':{'name': 'Park Falls','loc':'Wisconsin, USA','lat':45.945,'lon':269.727,'alt':442},
				'oc':{'name': 'Lamont','loc':'Oklahoma, USA','lat':36.604,'lon':262.514,'alt':320},
				'wg':{'name': 'Wollongong','loc':'Australia','lat':-34.406,'lon':150.879,'alt':30},
				'db':{'name': 'Darwin','loc':'Australia',tuple([datetime(2005,8,1),datetime(2015,7,1)]):{'lat':-12.422445,'lon':130.89154,'alt':30},
														tuple([datetime(2015,7,1)]):{'lat':-12.45606,'lon':130.92658,'alt':37}
						},
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
				'iz':{'name': 'Izana','loc':'Spain','lat':28,'lon':344,'alt':0},
				'if':{'name': 'Idianapolis','loc':'Indiana, USA','lat':39.861389,'lon':273.996389,'alt':270},
				'df':{'name': 'Dryden','loc':'California, USA','lat':34959917,'lon':242.118931,'alt':700},
				'js':{'name': 'Saga','loc':'Japan','lat':33.240962,'lon':130.288239,'alt':7},
				'fc':{'name': 'Four Corners','loc':'USA','lat':36.79749,'lon':251.51991,'alt':1643},
				'ci':{'name': 'Pasadena','loc':'California, USA','lat':34.13623,'lon':241.873103,'alt':230},
				'rj':{'name': 'Rikubetsu','loc':'Japan','lat':43.4567,'lon':143.7661,'alt':380},
				'pr':{'name': 'Paris','loc':'France','lat':48.846,'lon':2.356,'alt':60},
				'ma':{'name': 'Manaus','loc':'Brazil','lat':-3.2133,'lon':299.4017,'alt':50},
				'sp':{'name': 'Ny-Alesund','loc':'Norway','lat':78.92324,'lon':11.92298,'alt':0},
				'et':{'name': 'East Trout Lake','loc':'Canada','lat':54.353738,'lon':255.013333,'alt':501.8},
				'an':{'name': 'Anmyeondo','loc':'Korea','lat':36.5382,'lon':126.331,'alt':30},
				'bu':{'name': 'Burgos','loc':'Philippines','lat':18.5325,'lon':120.6496,'alt':35},
				'we':{'name': 'Jena','loc':'Austria','lat':50.91,'lon':11.57,'alt':211.6},
				}

	GGGPATH = os.environ['GGGPATH']
	print 'GGGPATH =',GGGPATH

	mod_path = os.path.join(GGGPATH,'models','gnd')

	argu = sys.argv

	site_abbrv = argu[1]
	year = int(argu[2])

	mode = argu[3].lower() # ncep or merra
	print mode
	if ('merra' not in mode) and ('ncep' not in mode):
		print 'Wrong mode, you must specify either "ncep" or "merra"'
		sys.exit()
	if 'merra' in mode:
		try:
			username = mode.split(',')[1]
			password = mode.split(',')[2]
		except IndexError:
			print 'When using MERRA mode, the argument must be "merra,username,password"'
			sys.exit()

	try:
		print 'Site:',site_dict[site_abbrv]['name'],site_dict[site_abbrv]['loc']
	except KeyError:
		print 'Wrong 2 letter site abbreviation (check the site_dict dictionary)'
		sys.exit()
	print 'lat,lon,masl:',site_dict[site_abbrv]['lat'],site_dict[site_abbrv]['lon'],site_dict[site_abbrv]['alt']
	print 'Year:',year

	date = datetime(year,1,1)
	one_day = timedelta(days=1)

	if len(argu)>4:
		HH = int(argu[4])
		MM = int(argu[5])
	else:
		HH = 12
		MM = 0

	if HH>=24 or HH<0:
		print 'Need 0<=H<24'
		sys.exit()
	if MM>=60 or MM<0:
		print 'Need 0<=MM<60'
		sys.exit()

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

	if site_lat > 0:
		ns = 'N'
	else:
		ns = 'S'

	if site_lon_180>0:
		ew = 'E'
	else:
		ew = 'W'

	if mode == 'ncep':

		# path to the netcdf files
		ncdf_AT_file = os.path.join(GGGPATH,'ncdf','.'.join(['air','{:0>4}'.format(year),'nc']))
		ncdf_GH_file = os.path.join(GGGPATH,'ncdf','.'.join(['hgt','{:0>4}'.format(year),'nc']))
		ncdf_SH_file = os.path.join(GGGPATH,'ncdf','.'.join(['shum','{:0>4}'.format(year),'nc']))

		print 'Read global data ...'
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

	rmm = 28.964/18.02	# Ratio of Molecular Masses (Dry_Air/H2O)
	gravity_at_lat, earth_radius_at_lat = gravity(site_lat,site_alt/1000.0)

	while date.year==year:

		YYYYMMDD = time.strftime('%Y%m%d',date.timetuple())
		mod_name = '_'.join(['NCEP',YYYYMMDD,'{:0>2.0f}'.format(round(abs(site_lat)))+ns,'{:0>3.0f}'.format(round(abs(site_lon_180)))+ew+'.mod'])
		mod_file_path = os.path.join(mod_path,mod_name)
		print mod_name

		if 'merra' in mode:

			# in the url below, M2I3NVASM is for every 3 hour instantaneous assimilated meteorological fields; VASM is for 72 levels; PASM is for 42 levels
			if 'merra42' in mode:
				letter = 'P'
			elif 'merra72' in mode:
				letter = 'V'

			# levels data
			url = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2I3N{}ASM.5.12.4/{:0>4}/{:0>2}/MERRA2_400.inst3_3d_asm_N{}.{:0>4}{:0>2}{:0>2}.nc4'.format(letter,date.year,date.month,letter.lower(),date.year,date.month,date.day)
			session = setup_session(username,password,check_url=url)
			dataset = open_url(url,session=session)

			# surface data, from M2I1NXASM ; it is on the goldsmr4 server and not goldsmr5 like above
			surface_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2I1NXASM.5.12.4/{:0>4}/{:0>2}/MERRA2_400.inst1_2d_asm_Nx.{:0>4}{:0>2}{:0>2}.nc4'.format(date.year,date.month,date.year,date.month,date.day)
			surface_session = setup_session(username,password,check_url=surface_url)
			surface_dataset = open_url(surface_url,session=surface_session)

			# get the min/max lat-lon indices of merra lat-lon that lies within a 5x5 lat-lon box centered on the site lat-lon.
			min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID = merra_querry_indices(dataset,site_lat,site_lon_180,2.5,2.5)

			# surface data
			# 2 meter Air Temperature
			lev_surf_AT,lat_surf_AT, lon_surf_AT, tim_surf_AT, data_surf_AT, data_scale_factor_surf_AT, data_add_offset_surf_AT, julday0 = read_data(surface_dataset,'T2M',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)		
			# 2 meter Specific humidity
			lev_surf_SH,lat_surf_SH, lon_surf_SH, tim_surf_SH, data_surf_SH, data_scale_factor_surf_SH, data_add_offset_surf_SH, julday0 = read_data(surface_dataset,'QV2M',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)		
			# Surface pressure
			lev_surf_P,lat_surf_P, lon_surf_P, tim_surf_P, data_surf_P, data_scale_factor_surf_P, data_add_offset_surf_P, julday0 = read_data(surface_dataset,'PS',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)								
			# surface geopotential height
			lev_surf_GH,lat_surf_GH, lon_surf_GH, tim_surf_GH, data_surf_GH, data_scale_factor_surf_GH, data_add_offset_surf_GH, julday0 = read_data(dataset,'PHIS',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID) # surface geopotential height is in the levels datasets			
			data_surf_GH = data_surf_GH / gravity_at_lat # convert from m2 s-2 to m

			# levels data
			# Air temperature
			lev_AT,lat_AT, lon_AT, tim_AT, data_AT, data_scale_factor_AT, data_add_offset_AT, julday0 = read_data(dataset,'T',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Specific humidity
			lev_SH,lat_SH, lon_SH, tim_SH, data_SH, data_scale_factor_SH, data_add_offset_SH, julday0 = read_data(dataset,'QV',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)
			# Height
			lev_H,lat_H, lon_H, tim_H, data_H, data_scale_factor_H, data_add_offset_H, julday0 = read_data(dataset,'H',min_lat_ID, max_lat_ID, min_lon_ID, max_lon_ID)

		# The mod file time, Noon by default
		site_tim = (int(round(sum(jdcal.gcal2jd(date.year,date.month,date.day))))-julday0)*24.0 + HH + MM/60.0 - site_lon_180/15.0

		# interpolate the data to the site's location and the desired time
		site_AT = trilinear_interp(data_AT, data_scale_factor_AT, data_add_offset_AT, site_lon_360, lon_AT, site_lat, lat_AT, site_tim, tim_AT) 
		site_SH = trilinear_interp(data_SH, data_scale_factor_SH, data_add_offset_SH, site_lon_360, lon_SH, site_lat, lat_SH, site_tim, tim_SH)

		if 'ncep' in mode:
			site_GH = trilinear_interp(data_GH, data_scale_factor_GH, data_add_offset_GH, site_lon_360, lon_GH, site_lat, lat_GH, site_tim, tim_GH)

		if 'merra72' in mode:
			# the merra 72 mid level pressures are not fixed, so need to interpolate to get just 1 array of levels
			lev_AT = trilinear_interp(lev_AT, lev_scale_factor_AT, lev_add_offset_AT, site_lon_360, lon_AT, site_lat, lat_AT, site_tim, tim_AT) 

		if 'merra' in mode:

			# interpolate height levels
			site_H = trilinear_interp(data_H, data_scale_factor_H, data_add_offset_H, site_lon_360, lon_H, site_lat, lat_H, site_tim, tim_H)
			# convert to geopotential
			site_GH = (site_H)/(1.0+(site_H)/earth_radius_at_lat)		# Convert from geometric to geopotential

			# interpolate the surface data to the site's location and the desired time
			site_surf_AT = trilinear_interp(data_surf_AT, data_scale_factor_surf_AT, data_add_offset_surf_AT, site_lon_360, lon_surf_AT, site_lat, lat_surf_AT, site_tim, tim_surf_AT) 
			site_surf_GH = trilinear_interp(data_surf_GH, data_scale_factor_surf_GH, data_add_offset_surf_GH, site_lon_360, lon_surf_GH, site_lat, lat_surf_GH, site_tim, tim_surf_GH)
			site_surf_SH = trilinear_interp(data_surf_SH, data_scale_factor_surf_SH, data_add_offset_surf_SH, site_lon_360, lon_surf_SH, site_lat, lat_surf_SH, site_tim, tim_surf_SH)
			site_surf_P = trilinear_interp(data_surf_P, data_scale_factor_surf_P, data_add_offset_surf_P, site_lon_360, lon_surf_P, site_lat, lat_surf_P, site_tim, tim_surf_P)

			print site_AT.shape,site_SH.shape,site_GH.shape,lev_AT.shape

			without_fill_IDs = np.where(site_AT<1e10)

			# get rid of fill values
			site_AT = site_AT[without_fill_IDs]
			site_GH = site_GH[without_fill_IDs]
			site_SH = site_SH[without_fill_IDs]
			lev_AT = lev_AT[without_fill_IDs]

			print site_AT.shape,site_SH.shape,site_GH.shape,lev_AT.shape

			#insert surface values
			try:
				insert_ID = np.where(lev_AT>site_surf_P/100)[0][0]
			except IndexError:
				insert_ID = 0
			
			site_AT = np.insert(site_AT,insert_ID,site_surf_AT)
			site_SH = np.insert(site_SH,insert_ID,site_surf_SH)
			site_GH = np.insert(site_GH,insert_ID,site_surf_GH)
			lev_AT = np.insert(lev_AT,insert_ID,site_surf_P/100) # site_surf_P is in Pa, lev_AT is in hPa

		site_H2OMF = rmm*site_SH/(1-site_SH*(1-rmm)) 	# Convert Mass Mixing Ratio to Mole Fraction
		site_GH = site_GH/1000.0	# Convert m to km
		site_TP = 0.0

		write_mod(mod_file_path,site_lat,lev_AT,site_AT,site_GH,site_TP,site_H2OMF)

		date = date + one_day