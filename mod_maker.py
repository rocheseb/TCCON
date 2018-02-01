"""
mod_maker10.f translated into python

What's different?

This is run with:

python mod_maker.py arg1 arg2 ar3 arg4

arg1: two letter site abbreviation
arg2: year (YYYY)
arg3: (optional, default=12) hour (HH) for time interpolation
arg4: (optional, default=0) minute (MM) for time interpolation

In GGGPATH/models/gnd it will write one mod file per day if the given year.
This require the global NCEP netcdf files of the given year to be present in GGGPATH/ncdf

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
	Creates a GGG-format .mod file when called by mod_maker.pro.
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

def trilinear_interp1(fin, fscale_factor, fadd_offset, site_lon_360, lon_XX, site_lat, lat_XX, site_tim, tim_XX):
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
	
	fout =	((fin[index_tt,:,index_yy,index_xx]*(1.0-fr_xx) \
	+ fin[index_tt,:,index_yy,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
	+ (fin[index_tt,:,index_yy+1,index_xx]*(1.0-fr_xx) \
	+ fin[index_tt,:,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*(1.0-fr_tt) \
	+ ((fin[index_tt+1,:,index_yy,index_xx]*(1.0-fr_xx) \
	+ fin[index_tt+1,:,index_yy+1,ixpomnxx]*fr_xx)*(1.0-fr_yy) \
	+ (fin[index_tt+1,:,index_yy+1,index_xx]*(1.0-fr_xx) \
	+ fin[index_tt+1,:,index_yy+1,ixpomnxx]*fr_xx)*fr_yy)*fr_tt 
	
	fout = fout*fscale_factor + fadd_offset

	return fout

def read_global_data_file(ncdf_file, varname):

	dataset = netCDF4.Dataset(ncdf_file,'r')	# open the netcdf files

	lev_XX = dataset.variables['level'][:]	# Read in variable 'level'
	lat_XX = dataset.variables['lat'][:]	# Read in variable 'lat'
	lon_XX = dataset.variables['lon'][:]	# Read in variable 'lon'
	tim_XX = dataset.variables['time'][:]	# Read in variable 'lon'

	print varname,': Nlev,Nlat,Nlon,Ntim =',len(lev_XX),len(lat_XX),len(lon_XX),len(tim_XX)

	# Set julday0 => Julian date where time=0
	time_units = dataset.variables['time'].units  # string containing definition of time units

	# two lines to parse the date (no longer need to worry about before/after 2014)
	date_list = re.findall(r"[\w]+",time_units.split('since')[1])
	common_date_format = '{:0>4}-{:0>2}-{:0>2} {:0>2}:{:0>2}:{:0>2}'.format(*date_list)

	start_date = time.strptime(common_date_format,'%Y-%m-%d %H:%M:%S')

	julday0 = int(round(sum(jdcal.gcal2jd(*start_date[:3])))) # gives same results as IDL's JULDAY function

	# Read global data
	global_data_XX = dataset.variables[varname][:]	# Read in variable 'lon'

	# Initialize add_offset and scale_factor
	global_data_add_offset_XX = 0
	global_data_scale_factor_XX = 1.0

	# Check variable for offset and scaling factor attributes
	for attribute in dataset.variables[varname].ncattrs():
		# Set global_data_add_offset_XX
		if attribute == 'add_offset':
			global_data_add_offset_XX = dataset.variables[varname].getncattr(attribute)
		# Set global_data_scale_factor_XX
		if attribute == 'scale_factor':
			global_data_scale_factor_XX = dataset.variables[varname].getncattr(attribute)

	dataset.close()

	return lev_XX, lat_XX, lon_XX, tim_XX, global_data_XX, global_data_scale_factor_XX, global_data_add_offset_XX, julday0


if __name__ == "__main__":

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

	print 'Site:',site_dict[site_abbrv]['name'],site_dict[site_abbrv]['loc']
	print 'lat,lon,masl:',site_dict[site_abbrv]['lat'],site_dict[site_abbrv]['lon'],site_dict[site_abbrv]['alt']
	print 'Year:',year

	date = datetime(year,1,1)
	one_day = timedelta(days=1)

	if len(argu)>3:
		HH = int(argu[3])
		MM = int(argu[4])
	else:
		HH = 12
		MM = 0

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

	# path to the netcdf files
	ncdf_AT_file = os.path.join(GGGPATH,'ncdf','.'.join(['air','{:0>4}'.format(year),'nc']))
	ncdf_GH_file = os.path.join(GGGPATH,'ncdf','.'.join(['hgt','{:0>4}'.format(year),'nc']))
	ncdf_SH_file = os.path.join(GGGPATH,'ncdf','.'.join(['shum','{:0>4}'.format(year),'nc']))

	print 'Read global data ...'
	# Air Temperature
	lev_AT,lat_AT, lon_AT, tim_AT, global_data_AT, global_data_scale_factor_AT, global_data_add_offset_AT, julday0	= read_global_data_file(ncdf_AT_file, 'air')
	if len(lev_AT) < 17:
		print 'Need 17 levels of AT data: found only ',len(lev_AT)
	
	# Geopotential Height
	lev_GH,lat_GH, lon_GH, tim_GH, global_data_GH, global_data_scale_factor_GH, global_data_add_offset_GH, julday0	= read_global_data_file(ncdf_GH_file, 'hgt')
	if len(lev_GH) < 17:
		print 'Need 17 levels of GH data: found only ',len(lev_GH)
	
	# Specific Humidity
	lev_SH,lat_SH, lon_SH, tim_SH, global_data_SH, global_data_scale_factor_SH, global_data_add_offset_SH, julday0	= read_global_data_file(ncdf_SH_file, 'shum')
	if len(lev_SH) <  8:
		print 'Need  8 levels of SH data: found only ',len(lev_SH)

	rmm = 28.964/18.02	# Ratio of Molecular Masses (Dry_Air/H2O)
	while date.year==year:

		YYYYMMDD = time.strftime('%Y%m%d',date.timetuple())
		mod_name = '_'.join(['NCEP',YYYYMMDD,'{:0>2.0f}'.format(round(abs(site_lat)))+ns,'{:0>3.0f}'.format(round(abs(site_lon_180)))+ew+'.mod'])
		mod_file_path = os.path.join(mod_path,mod_name)
		print mod_name

		site_tim=(int(round(sum(jdcal.gcal2jd(date.year,date.month,date.day))))-julday0)*24.0 + HH + MM/60.0 - site_lon_180/15.0

		site_AT = trilinear_interp1(global_data_AT, global_data_scale_factor_AT, global_data_add_offset_AT, site_lon_360, lon_AT, site_lat, lat_AT, site_tim, tim_AT) 
		site_GH = trilinear_interp1(global_data_GH, global_data_scale_factor_GH, global_data_add_offset_GH, site_lon_360, lon_GH, site_lat, lat_GH, site_tim, tim_GH) 
		site_SH = trilinear_interp1(global_data_SH, global_data_scale_factor_SH, global_data_add_offset_SH, site_lon_360, lon_SH, site_lat, lat_SH, site_tim, tim_SH) 

		site_H2OMF = rmm*site_SH/(1-site_SH*(1-rmm))  # Convert Mass Mixing Ratio to Mole Fraction
		site_GH = site_GH/1000.0   # Convert m to km
		site_TP = 0.0                 # Convert m to km

		write_mod(mod_file_path,site_lat,lev_AT,site_AT,site_GH,site_TP,site_H2OMF)

		date = date + one_day