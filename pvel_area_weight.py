"""
Note:

Reading GEOS/MERRA data yields arrays with dimensions:
	- (lon,lat,lev,time) from IDL
	- (time,lev,lat,lon) from Python
"""
import numpy as np 


def pvel(pv,lon,lat):
	"""
	1.1 08/04/05

	Calculation of PV area equivalent latitude following the method described in Allen and Nakamura 2003, JAS 60, p 287-303

	Inputs: PV on 3D or 4D regular latlon grid, lon and lat vectors
	Returns: PVEL on same 3D or 4D grid

	Alan Geer 4/8/5 see 17-151
	"""

	nTimes,nLev,nLat,nLon = pv.shape

	if nLon != len(lon) or nLat != len(lat):
		print "lon and lat do not match PV field"

	# Calculate grid box area in m^3
	box_area = area_weights(lon, lat, normalise=False)
	earth_area = np.sum(box_area)

	# For each level and time, order PV and calculate PVEL
	R = 6371e3 # earth's radius, m, from back of Houghton

	pvel = np.full(pv.shape,np.nan)

	for iLev in range(nLev):
		for iTime in range(nTimes):
			iPVsort = pv[iTime,iLev,:,:].argsort(axis=None,kind='mergesort')
			pv_area = box_area.T.flatten()[iPVsort].cumsum()
			# Note it's important to use the calculated total earth area in the fraction to avoid floating point wobbles:
			rel_area = pv_area/(pv_area[nLon*nLat-1]/2.) - 1.
			iSortLons = iPVsort % nLon
			iSortLats = iPVsort/nLon
			pvel[iTime,iLev,iSortLats,iSortLons] = np.rad2deg(np.arcsin(rel_area)) 

	return pvel


def area_weights(longitudes,latitudes,normalise=True):
	"""
	1.3 03/31/05

	Calculates area weights based for a regular lat/lon grid
	based on the input vectors of longitude and latitude used on that
	grid.

	Area weights are normalised so that sum(area_weights) = 1
	Hence sum(area_weights*field) gives the area weighted mean of
	a lon/lat field.

	If "normalise" is true: gives the area in m^3 instead

	Longitudes and latitudes expected as arrays of for all grid points, in
	degrees

	Based on 11-136

	Alan Geer 11/11/2003
	"""
	earth_r = 6371e3 # km, from back of Houghton

	nLons = len(longitudes)
	nLats = len(latitudes)
	lon_width = longitudes[1:nLons]-longitudes[0:nLons-1]
	lat_width = latitudes[1:nLats]-latitudes[0:nLats-1]

	delta_lon = np.deg2rad(np.abs(lon_width[0]/2.0))
	delta_lat = np.deg2rad(np.abs(lat_width[0]/2.0))      

	# Calculate a latitude slice of weights (formula from 11-136) for
	# grid boxes not centred at the North or South pole  
	lat_weights = 4 * earth_r**2 * delta_lon * np.sin(delta_lat) * np.cos(np.deg2rad(latitudes))

	# Now check if any boxes are centred on the North or South pole and
	# use a separate formula (from 11-137) for a triangular-shaped area
	# with its tip at the pole 
	iPoles = np.where(np.abs(latitudes)==90)[0]
	count = len(iPoles)
	if count != 0:
		lat_weights[iPoles] = earth_r**2 * delta_lon * delta_lat**2

	# Convert to a lon/lat field and normalise (a waste of space but easy to use)
	area_weights = np.outer(np.ones(nLons),lat_weights) # NB total of this is earth's surface area

	if normalise:
		area_weights = area_weights/np.sum(area_weights)   

	return area_weights




