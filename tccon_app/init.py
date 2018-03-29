'''
This code includes objects that are used to set up main.py
Users can modify inputs here
'''

def setup():
	'''
	returns the variables needed in main.py
	This should be the only place where users should apply modifications
	'''

	layout_mode = 'comp' # can be set to 'simple' or 'comp'; in 'simple' mode there is just one plot with one site and one variable to select

	cache_max_size = 2E8 # maximum size of the cache file (in bytes), any new cached data after that will remove the oldest data following certain rules (see add_cache function)
	# set cache_max_size to 0 if you do not want to make a cache file (reading from netcdf files can be faster than from the cache file)

	# if you modify the plotting colors, you will need to remove your cache_dic.npy file; you will also need to edit the styles.css in tccon_app/templates/styles.css to match the new colors.
	main_color = 'yellowgreen' # this will be the color used for the flag=0 data; I use css 'YellowGreen' (#9ACD32) by default
	main_color2 = 'plum' # the main color for data from the second site in 'comp' mode; I use css 'Plum' (#DDA0DD) by default
	flag_color = 'grey' # this will be the color used for the flag!=0 data
	hover_color = 'red' # this will be the color used for data hovered by the mouse

	# this is the mode of the hovertool for the Figure 1 and Figure 2 in 'comp' mode
	# 'width' to have all y data select over a x range (you can only pan the box selection left-right)
	# 'height' can only pan the box selection up-down and it selects all the x data
	# 'both' you can select a specific rectangle of data
	boxselecttool_dimensions = 'both'

	# variables including the keywords in 'skip_list' won't be shown in the variable input dropdowns
	skip_list = [	'_Version','ak_','prio','checksum','graw','spectrum','year','ada','aicf','adcf','time',
					'_OVC_','_VSF_','_AM_','_ZO','_Zpres','opd_cm',] 

	# dictonnary with the full names of TCCON sites
	T_FULL = {
				'pa':'Park Falls',
				'oc':'Lamont',
				'wg':'Wollongong',
				'db':'Darwin',
				'or':'Orleans',
				'bi':'Bialystok',
				'br':'Bremen',
				'jc':'JPL 01',
				'jf':'JPL 02',
				'ra':'Reunion Island',
				'gm':'Garmisch',
				'lh':'Lauder 01',
				'll':'Lauder 02',
				'tk':'Tsukuba 02',
				'ka':'Karlsruhe',
				'ae':'Ascension Island',
				'eu':'Eureka',
				'so':'Sodankyla',
				'iz':'Izana',
				'if':'Indianapolis',
				'df':'Dryden',
				'js':'Saga',
				'fc':'Four Corners',
				'ci':'Pasadena',
				'rj':'Rikubetsu',
				'pr':'Paris',
				'ma':'Manaus',
				'sp':'Ny-Alesund',
				'et':'East Trout Lake',
				'an':'Anmyeondo',
				'bu':'Burgos',
				'we':'Jena',
				#to add a new site (like for em27 data), just make up a new 2 letter site abbreviation and add it to this dictionary like this:
				'zf':'site or instrument name',
			 }

	# dictonnary with the country/state of TCCON sites
	T_LOC =  {	
				'pa':'Wisconsin, USA',
				'oc':'Oklahoma, USA',
				'wg':'Australia',
				'db':'Australia',
				'or':'France',
				'bi':'Poland',
				'br':'Germany',
				'jc':'California, USA',
				'jf':'California, USA',
				'ra':'France',
				'gm':'Germany',
				'lh':'New Zealand',
				'll':'New Zealand',
				'tk':'Japan',
				'ka':'Germany',
				'ae':'United Kingdom',
				'eu':'Canada',
				'so':'Finland',
				'iz':'Spain',
				'if':'Indiana, USA',
				'df':'California, USA',
				'js':'Japan',
				'fc':'USA',
				'ci':'California, USA',
				'rj':'Japan',
				'pr':'France',
				'ma':'Brazil',
				'sp':'Norway',
				'et':'Canada',
				'an':'South Korea',
				'bu':'Philippines',
				'we':'Germany',
				#to add a new site (like for em27 data), just make up a new 2 letter site abbreviation and add it to this dictionary like this:
				'zf':'site or instrument location',
			 }

	return layout_mode, cache_max_size, main_color, main_color2, flag_color, hover_color, boxselecttool_dimensions, skip_list, T_FULL, T_LOC
