#!/var/lib/py27_sroche/bin/python
 # -*- coding: utf-8 -*-

from __future__ import print_function # allows the use of Python 3.x print(function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# Code description #
####################

'''
Read TCCON data from netCDF files.
Creates .html plots with bokeh

compare all tccon sites to a selected site for Xspecies
e.g. with all sites compared to Lamont  http://www.atmosp.physics.utoronto.ca/~sroche/FREQ_TCCON_1_days_Lamont.html

or just get a full time series for a specific variable
e.g.  http://www.atmosp.physics.utoronto.ca/~sroche/TCCON_datetime_xco2_ppm.html

LARGE (~60 MB for the full time series) resulting files that can take ~20 sec to load but then have fluid interactions

How to use:

python TCCON_comp.py

Just answer the questions.
All netcdf files must be in a user specified folder.

There are two "modes" for the program. If you answer "no" to "skip to plotting", you will be guided to write new plotting input files.
If you have already written plotting input files, you can answer "yes" and will be asked which variables you want to plot.

You will be asked if you want to use DATA or FREQ_DATA:
- DATA is for full time series
- FREQ_DATA contains data averaged with a certain frequency (e.g. weekly / daily) and will be used for comparison plots

TIPS:
the first time you write a plotting input file, it is going to take a while (especially if you do matching and averaging with a frequency smaller than 1 weeks).
make sure you write in all the variables that you may be using later so you don't need to run the writting mode multiple times.

'flag','asza_deg','time','lat_deg', and 'long_deg' are read by default, you don't need to specify them

'''

#############
# Libraries #
#############

# general
import os
import sys

# netcdf reading/writing
import netCDF4

# time handling
import calendar
import time
from datetime import datetime
from datetime import timedelta

# special arrays with special functions
import numpy as np

# round up
from math import ceil

# interactive plots
from bokeh.plotting import figure, output_file
from bokeh.models import Legend, Panel, Tabs, CustomJS, ColumnDataSource, CheckboxGroup, RadioGroup, Button, VBox, PreText, Range1d
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.layouts import gridplot, widgetbox
from bokeh.resources import CDN
from bokeh.embed import file_html

#############
# Functions #
#############

#a fancy loadbar to be displayed in the prompt while executing a time consuming loop.
def progress(i,tot,bar_length=20,char=''):
	percent=(float(i)+1.0)/tot
	hashes='#' * int(round(percent*bar_length))
	spaces=' ' * (bar_length - len(hashes))
	sys.stdout.write("\rPercent:[{0}] {1}%".format(hashes + spaces, int(round(percent * 100)))+"    "+str(i+1)+"/"+str(tot)+' Site: '+char+'                    ')
	sys.stdout.flush()


#########
# SETUP #
#########

######## /!\ important time handling to make sure times don't get shifted from UTC due to computer environment variables when using datetime objects ###########
os.environ['TZ'] = 'UTC'
time.tzset()
# This will not change your system timezone settings
################################################################################################################################################################

print('\n\nPut in a folder all the .ncdf files (containing the complete eof data) you want to process, give the path to that folder')

TCCON_path='1'
while os.path.isdir(TCCON_path)==False:
	TCCON_path=raw_input('Give the path to your folder /YOUR/PATH/TO/FILES  :\n')
	if os.path.isdir(TCCON_path)==False:
		print('/!\\ You gave a wrong path /!\\\n')

save_path=os.path.join(TCCON_path,'SAVE')
if not os.path.isdir(save_path):
	os.makedirs(save_path)

kelly_colors = {	
					'vivid_yellow':(255, 179, 0),
					'strong_purple':(128, 62, 117),
					'vivid_orange':(255, 104, 0),
					'very_light_blue':(166, 189, 215),
					'vivid_red':(193, 0, 32),
					'grayish_yellow':(206, 162, 98),
					'medium_gray':(129, 112, 102),

					# these aren't good for people with defective color vision:
					'vivid_green':(0, 125, 52),
					'strong_purplish_pink':(246, 118, 142),
					'strong_blue':(0, 83, 138),
					'strong_yellowish_pink':(255, 122, 92),
					'strong_violet':(83, 55, 122),
					'vivid_orange_yellow':(255, 142, 0),
					'strong_purplish_red':(179, 40, 81),
					'vivid_greenish_yellow':(244, 200, 0),
					'strong_reddish_brown':(127, 24, 13),
					'vivid_yellowish_green':(147, 170, 0),
					'deep_yellowish_brown':(89, 51, 21),
					'vivid_reddish_orange':(241, 58, 19),
					'dark_olive_green':(35, 44, 22),

					#filling in for the last sites
					'gold':'gold',
					'chartreuse':'chartreuse',
					'cyan':'cyan',
					'firebrick':'firebrick',
					'lightsalmon':'lightsalmon',
					'peru':'peru',
					'goldenrod':'goldenrod',
					'navy':'navy',
					'green':'green'
				}

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
		 }

# dictonnary with the country/state of TCCON sites
T_LOC =  {	
			'pa':' Wisconsin, USA',
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
		 }

site_files =[i for i in os.listdir(TCCON_path) if '.' in i]

# start a new analysis or use existing files
Q = ''
while Q not in ['Y','y','N','n']:
	Q=raw_input("\nSkip to plotting? (y/n):\n")
	if Q not in ['Y','y','N','n']:
		print('Your must type y for yes or n for no')

# plots tools
TOOLS = "pan,wheel_zoom,box_zoom,undo,redo,reset,save"

#############
# MAIN CODE #
#############

codestart = time.time() #To keep track of time while the code runs

if Q in ['N','n']:

	select_vars = (raw_input('\nGive the exact name of variables to read (.EOF headers), separate by commas (lat,long,flag,time,sza are read by default) (e.g. "xco2_ppm,xch4_ppm,xn2o_ppb,xco2_error_ppm"):\n')).split(',')
	
	sza_check = False
	while sza_check == False:
		SZA=raw_input('\nSelect data with SZA below (give a number or just press enter to take data with any SZA):\n')
		if SZA == '':
			sza_check = True
		try:
			float(SZA)
			sza_check = True
		except ValueError:
			pass

	QF = ''
	while QF not in ['Y','y','N','n']:
		QF=raw_input("\nSkip matching and averaging? (y/n):\n")
		if QF not in ['Y','y','N','n']:
			print('Your must type y for yes or n for no')

	if QF in ['N','n']:

		print('\nTCCON sites and two letter abbreviations:\n')
		for i in T_FULL.keys():
			print('\t-',T_FULL[i],':',i)

		SELECT = ''
		while SELECT not in T_FULL.keys():
			SELECT = raw_input('\nGive the two letter abreviation of the site you wish to compare with other TCCON sites (e.g. "eu" for Eureka):\n')
			if SELECT not in T_FULL.keys():
				print('Wrong entry')
		SELECT = T_FULL[SELECT]

		# Ask user if (s)he wants to set a custom time range; if not, use the time range of the selected site
		QS=''
		while QS not in ['Y','y','N','n']:
			QS=raw_input("Use a specific time range? (y/n) if no, "+SELECT+" full time range will be used:\n")
			if QS not in ['Y','y','N','n']:
				print('Your must type y for yes or n for no')

		if QS in ['Y','y']:		
			time_switch = True
		else:
			time_switch = False

		# if time_switch = True, get a user specified time range
		if time_switch == True:

			t0 = ''
			while type(t0) != datetime: 
				start = raw_input("Starting date? YYYY-MM-DD-HH (e.g. 2010-01-02-14 for 2 PM January 2nd 2010, -HH is optional):\n")
			
				if len(start) in [10,13]:	
					if len(start) == 13:
						t0 = datetime.strptime(start,'%Y-%m-%d-%H')
					else:
						t0 = datetime.strptime(start,'%Y-%m-%d')
				else:
					print("Invalid input")
					t0 = ''

			tf = ''	
			while type(tf) != datetime: 
				end = raw_input("End date? YYYY-MM-DD-HH (e.g. 2010-01-02-14 for 2 PM January 2nd 2010, -HH is optional):\n")
				if len(end) in [10,13]:
					if len(end) == 13:
						tf = datetime.strptime(end,'%Y-%m-%d-%H')
					else:
						tf = datetime.strptime(end,'%Y-%m-%d')
				else:
					print("Invalid input")
					tf = ''

				if type(tf) == datetime:
					if tf < t0:
						print("The starting date shall precede the end date")
						tf = ''

		# Frequency of data matching and averaging
		switch = False # for checking the validity of user input
		while switch == False:
			FREQ = raw_input("Frequency of data matching and averaging? (e.g. '1 hours'; '5.62 days','7.1 weeks', put an 's' even for 1):\n")
			if len(FREQ.split()) == 2: # there should be exactly two arguments
				if FREQ.split()[1] in ['hours','days','weeks']: # the second argument must be in that list
					try: # the first argument must be a number
						float(FREQ.split()[0])
						switch = True # if above conditions are met, get out of the while loop, otherwise output an error message and prompt user with the question again
					except ValueError:
						print("Invalid input (1)")
				else:
					print("Invalid input (2)")			
			else:
				print("Invalid input (3)")

		if FREQ.split()[1] == 'weeks':
			time_step = timedelta(weeks=float(FREQ.split()[0]))
		if FREQ.split()[1] == 'days':
			time_step = timedelta(days=float(FREQ.split()[0]))
		if FREQ.split()[1] == 'hours':
			time_step = timedelta(hours=float(FREQ.split()[0]))

		frequency = time_step.total_seconds()
	
	###########################
	# Read TCCON netCDF files #
	###########################

	milestone = time.time()
	print('\nRead TCCON netCDF files ...')

	ALL_DATA = {}
	prev_site = ''

	for file in site_files:

		site_path = os.path.join(TCCON_path,file)
		
		site = T_FULL[file[:2]] # filnames start with the site two letter abbreviation

		check = False
		if site == prev_site:
			check = True
		else:
			print('\n',site)
			data_site = {}

		print(file)

		f = netCDF4.Dataset(site_path,'r')

		try:
			for var in select_vars:
				if check == False:
					data_site[var] = f.variables[var][:]
					if ('time' not in select_vars) and (var == select_vars[0]):
						data_site['time'] = f.variables['time'][:]
						data_site['datetime'] = np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in f.variables['time'][:]])
					if ('asza_deg' not in select_vars) and (var == select_vars[0]):
						data_site['asza_deg'] = f.variables['asza_deg'][:]
					if ('flag' not in select_vars) and (var == select_vars[0]):
						data_site['flag'] = f.variables['flag'][:]
					if ('lat_deg' not in select_vars) and (var == select_vars[0]):
						data_site['lat_deg'] = f.variables['lat_deg'][:]	
					if ('long_deg' not in select_vars) and (var == select_vars[0]):
						data_site['long_deg'] = f.variables['long_deg'][:]									
				else:
					data_site[var] = np.append(data_site[var][:],f.variables[var][:])
					if ('time' not in select_vars) and (var == select_vars[0]):
						data_site['time'] = np.append(data_site['time'][:],f.variables['time'][:])
						data_site['datetime'] = np.append( data_site['datetime'][:] , np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in f.variables['time'][:]]) )		
					if ('asza_deg' not in select_vars) and (var == select_vars[0]):
						data_site['asza_deg'] = np.append( data_site['asza_deg'][:],f.variables['asza_deg'][:] )	
					if ('flag' not in select_vars) and (var == select_vars[0]):
						data_site['flag'] = np.append( data_site['flag'][:],f.variables['flag'][:] )
					if ('lat_deg' not in select_vars) and (var == select_vars[0]):
						data_site['lat_deg'] = np.append( data_site['lat_deg'][:],f.variables['lat_deg'][:] )
					if ('long_deg' not in select_vars) and (var == select_vars[0]):
						data_site['long_deg'] = np.append( data_site['long_deg'][:],f.variables['long_deg'][:] )

			ALL_DATA[site] = data_site

		except KeyError:
			print(site,'has no',var)
			pass

		prev_site = site

		f.close()

	del data_site # we will just need the ALL_DATA[site] from now on

	print('\nRead TCCON netCDF files DONE in',time.time()-milestone,'seconds\n')

	#####################
	#    Filter DATA    #
	#####################

	milestone = time.time()
	print('Get flag=0 data ...')
	DATA = {} # flag 0 data

	for site in ALL_DATA:
		DATA[site] = {}
		no_flag = [i for i in range(len(ALL_DATA[site]['flag'])) if ALL_DATA[site]['flag'][i]==0]
		tot=len([var for var in ALL_DATA[site]])
		count=0
		for var in ALL_DATA[site]:
			progress(count,tot,char=site)
			count+=1
			if len(ALL_DATA[site][var])==len(ALL_DATA[site]['flag']):
				DATA[site][var] = np.array([ALL_DATA[site][var][i] for i in no_flag])
	print('\nGet flag=0 data DONE in',time.time()-milestone,'seconds\n')

	del ALL_DATA # memory relief as this can represent several GB of data

	if SZA!='':
		milestone = time.time()
		print('Get sza<'+SZA+' data ...')
		for site in DATA:
			print('\n',site)
			small_sza = [i for i in range(len(DATA[site]['asza_deg'])) if DATA[site]['asza_deg'][i]<=float(SZA)]
			tot=len([var for var in DATA[site]])
			count=0
			for var in DATA[site]:

				print(var)
				#progress(count,tot,char=site)
				count+=1
				if len(DATA[site][var])==len(DATA[site]['asza_deg']):
					DATA[site][var] = np.array([DATA[site][var][i] for i in small_sza])
				print(len(DATA[site][var]))
		print('Get sza<'+SZA+' data DONE in',time.time()-milestone,'seconds\n')

	if time_switch == False:
		t0 = DATA[SELECT]['datetime'][0]
		tf = DATA[SELECT]['datetime'][-1]
	
	if QF in ['N','n']:
		print(SELECT,'time range:\nStart',t0.strftime('%d-%m-%Y %H:%M'),'\nEnd',tf.strftime('%d-%m-%Y %H:%M'))
		
	########################
	# FREQly averaged data #
	########################

	if QF in ['N','n']:
		FREQ_DATA = {}

		span = int(ceil((tf-t0).total_seconds()/frequency))

		milestone = time.time()
		print('Dividing',SELECT,'time range in',span,'intervals of',FREQ)
		tmin = t0
		temp = [[] for i in range(span)]
		for interval in range(span):
			progress(interval,span,char=SELECT)
			tmax = tmin+time_step
			times = DATA[SELECT]['datetime'][(DATA[SELECT]['datetime']>=tmin) & (DATA[SELECT]['datetime']<tmax)]
			times_ID=np.nonzero(np.in1d(DATA[SELECT]['datetime'],times))[0]
			temp[interval] = [j for j in times_ID]
			tmin = tmax

		print('\ntimes DONE in',time.time()-milestone,'seconds')
		print(SELECT,'has',len([i for i in temp if i!=[]]),'intervals of',FREQ,'with data within the time range\n')

		for site in DATA:

			freq_data = {}

			freq_site_data = {}
			freq_select_data = {}

			if site != SELECT:

				milestone = time.time()
				print(site+':\nMatching and averaging:')
				tmin = t0
				temp1 = [[] for i in range(span)]
				# this loop can be very time consuming if the frequency of matching/averaging is small
				for interval in range(span):
					progress(interval,span,char=site)
					tmax=tmin+time_step
					if len(temp[interval])>0:
						times = DATA[site]['datetime'][(DATA[site]['datetime'] >= tmin) & (DATA[site]['datetime'] < tmax)]
						if len(times)>0:
							times_ID = np.nonzero(np.in1d(DATA[site]['datetime'],times))[0]
							temp1[interval] = [j for j in times_ID]
						else:
							temp1[interval] = []
					else:
						temp1[interval] = []
					tmin=tmax
				try:
					if len([i for i in temp1 if i!=[]])==0:
						print('\nmatching DONE in',time.time()-milestone,'seconds')
						print('(1) Matching intervals of',FREQ,'within the time range: 0 /',len([i for i in temp if i!=[]]),'\n')
						continue
				except NameError:
					print('\n',site,' has no data within the time range\n')
					continue

				tp=[i for i in temp if i!=[] and temp1[temp.index(i)]!=[]] # Matching indices for data in DATA[SELECT]
				tp1=[i for i in temp1 if i!=[] and temp[temp1.index(i)]!=[]] # Matching indices for data in DATA[site]

				print('\nmatching and averaging DONE in',time.time()-milestone,'seconds') # the averaging is in fact done after that but it takes less than 0.01 seconds for each site

				if len(tp)!=len(tp1):
					print('WARNING: length of matches is different: tp1=',len(tp1),'; tp=',len(tp))

				comn_var = [var for var in DATA[site] if var in DATA[SELECT]]
				for var in comn_var:
					if len(DATA[site][var]) == len(DATA[site]['time']):		
						if False not in [elem not in var for elem in ['date','year','day','hour','lat','lon','km','Ver']]:
							freq_select_data[var] = np.array([sum([DATA[SELECT][var][j] for j in i])/len(i) for i in tp])
							freq_site_data[var] = np.array([sum([DATA[site][var][j] for j in i])/len(i) for i in tp1])

				freq_select_data['datetime'] = np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in freq_site_data['time'][:]])
				freq_site_data['datetime'] = np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in freq_site_data['time'][:]])
				
				print('(2) Matching','intervals of',FREQ,'within the time range: ',len(tp),'/',len([i for i in temp if i!=[]]),'\n')

				freq_data[site] = freq_site_data
				freq_data[SELECT] = freq_select_data

				FREQ_DATA[site] = freq_data
				

	for site in DATA:
		ID = DATA.keys().index(site)
		color = kelly_colors.keys()[ID]
		DATA[site]['color'] = kelly_colors[color]

	np.save(os.path.join(save_path,'DATA.npy'),DATA)
	
	if QF in ['N','n']:
		np.save(os.path.join(save_path,'FREQ_DATA_'+'_'.join(FREQ.split())+'_'+SELECT+'.npy'),FREQ_DATA)

else: 
	try:
		DATA = np.load(os.path.join(save_path,'DATA.npy')).item()
	except IOError:
		print("You can't skip to plots because the data files don't exist\nExiting now ...\n")
		sys.exit()

QP = ''
while QP not in ['D','d','F','f']:
	QP=raw_input("\nPlot with DATA or FREQ_DATA? (d/f):\n")
	if QP not in ['D','d','F','f']:
		print('Your must type d for DATA or f for FREQ_DATA')


# for DATA
all_latitudes = sorted([DATA[site]['lat_deg'][0] for site in DATA])[::-1]

ALL_site_lat = {} # latitude per site in an ordered dictionary

for lat in all_latitudes:
	ALL_site_lat[all_latitudes.index(lat)] = {}
	for site in DATA:
		if DATA[site]['lat_deg'][0] == lat:
			ALL_site_lat[all_latitudes.index(lat)][site] = lat

#########
# Plots #
#########

print('\nVariables:\n','  '.join([var for var in DATA[DATA.keys()[0]]]),'\n')

if QP in ['D','d']:

	which=(raw_input('What plot do you want? xaxis,yaxis ( e.g. datetime,xair ):\n')).split(',')

	sources = {}
	for site in DATA:
		sources[site] = ColumnDataSource(data = {'x':DATA[site][which[0]],'y':DATA[site][which[1]]})

	print('Plotting:')

	lat_ordered_sites = [[key for key in dic][0] for dic in [ALL_site_lat[ID] for ID in ALL_site_lat]] # sites ordered by decreasing latitude

	min_x = min([DATA[site][which[0]][0] for site in DATA])
	max_x = max([DATA[site][which[0]][-1] for site in DATA])

	min_y = min([min(DATA[site][which[1]]) for site in DATA])
	max_y = max([max(DATA[site][which[1]]) for site in DATA])

	# we will set the y axis range at +/- 10% of the data max amplitude
	ampli = max_y - min_y
	min_y = min_y - ampli*0.1/100
	max_y = max_y + ampli*0.1/100

	milestone = time.time()

	if type(min_x) == datetime:
		fig = figure(output_backend = "webgl", title = 'TCCON '+which[1]+' vs '+which[0], y_range=[min_y,max_y], plot_width = 900, plot_height = 650, tools = TOOLS, toolbar_location = 'above', x_axis_type='datetime', x_range = Range1d(min_x,max_x)) 
	else:
		fig = figure(output_backend = "webgl", title = 'TCCON '+which[1]+' vs '+which[0], y_range=[min_y,max_y], plot_width = 900, plot_height = 650, tools = TOOLS, toolbar_location = 'above', x_range = Range1d(int(min_x),ceil(max_x))) 

	plots=[]
	for site in lat_ordered_sites:
		plots.append( fig.scatter(x='x',y='y',color=DATA[site]['color'],alpha=0.5,source=sources[site]) )

	N_plots = range(len(plots))

	legend=Legend(items=[(lat_ordered_sites[i],[plots[i]]) for i in N_plots],location=(0,0))
	fig.add_layout(legend,'right')
	fig.yaxis.axis_label = which[1]
	fig.xaxis.axis_label = which[0]

	checkbox = CheckboxGroup(labels=[site+' '+str(round(all_latitudes[lat_ordered_sites.index(site)],2)) for site in lat_ordered_sites],active=[],width=200)

	iterable = [('p'+str(i),plots[i]) for i in N_plots]+[('checkbox',checkbox)]

	checkbox_code = """var indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; };"""
	checkbox_code += ''.join(['p'+str(i)+'.visible = indexOf.call(checkbox.active, '+str(i)+') >= 0;' for i in N_plots])
	checkbox.callback = CustomJS(args={key: value for key,value in iterable}, code=checkbox_code)

	clear_button = Button(label='Clear all',width=200)
	clear_button_code = """checkbox.set("active",[]);"""+checkbox_code
	clear_button.callback = CustomJS(args={key: value for key,value in iterable}, code=clear_button_code)

	check_button = Button(label='Check all',width=200)
	check_button_code = """checkbox.set("active","""+str(N_plots)+""");"""+checkbox_code
	check_button.callback = CustomJS(args={key: value for key,value in iterable}, code=check_button_code)

	group = widgetbox(checkbox,clear_button,check_button)

	grid = gridplot([[fig,group]])

	print(' -writting TCCON_'+which[0]+'_'+which[1]+'.html')
	outfile=open(os.path.join(save_path,'TCCON_'+which[0]+'_'+which[1]+'.html'),'w')
	outfile.write(file_html(grid,CDN,which[1]))
	outfile.close()
	print('TCCON_'+which[0]+'_'+which[1]+'.html DONE in', time.time()-milestone,'seconds\n')

else:
	# for FREQ_DATA

	freq_filenames = [i for i in os.listdir(save_path) if 'FREQ_DATA' in i]

	rangelist = range(len(freq_filenames))

	for i in rangelist:
		print('\t',i,'-',freq_filenames[i])

	Qid = ''
	while Qid not in rangelist:
		try:
			Qid = int(raw_input('\n\nWhich file to use?:\n'))
		except ValueError: 
			pass
		if Qid not in rangelist:
			print('You must enter the corresponding number')

	FREQ_DATA = np.load(os.path.join(save_path,freq_filenames[Qid])).item()
	SELECT = freq_filenames[Qid].split('_')[-1].split('.')[0]

	FREQ = ' '.join(freq_filenames[Qid].split('_')[2:4])

	arbit = FREQ_DATA.keys()[0]
	variabs = FREQ_DATA[arbit][arbit].keys()
	print('\n\nVariables:\n')

	for i in range(len(variabs)):
		print('\t',i,' - ',variabs[i])

	numtab = ''
	vartoplot = []
	while len(numtab)<1:
		numtab=raw_input('\n\nWhat variable do you want to use? (give a series of comma separated numbers corresponding to each variable):\n')

	vartoplot = [variabs[int(i)] for i in numtab.split(',')]

	latitudes = sorted([DATA[site]['lat_deg'][0] for site in FREQ_DATA])[::-1] # latitudes of sites in decreasing order.

	site_lat = {}

	for lat in latitudes:
		site_lat[latitudes.index(lat)] = {}
		for site in FREQ_DATA:
			if DATA[site]['lat_deg'][0] == lat:
				site_lat[latitudes.index(lat)][site] = lat

	lat_ordered_sites = [[key for key in dic][0] for dic in [site_lat[ID] for ID in site_lat]]

	tabs = []
	TOOLS = "pan,wheel_zoom,box_zoom,undo,redo,reset,box_select,save"

	temp = [0 for i in FREQ_DATA]

	columns = [ 
					TableColumn(field='site',title='Site'),
					TableColumn(field='latitude',title='Latitude'),
					TableColumn(field='N',title='N'),
					TableColumn(field='RMS',title='RMS'),
					TableColumn(field='Bias',title='Bias'),
					TableColumn(field='Scatter',title='Scatter'),
					TableColumn(field='R',title='R'),
					]

	milestone = time.time()

	freq_sources = {}
	freq_cor_sources = {}

	for var in vartoplot:

		table_source = ColumnDataSource( data = {'site':lat_ordered_sites,'latitude':[round(DATA[i]['lat_deg'][0],2) for i in lat_ordered_sites],'N':temp,'RMS':temp,'Bias':temp,'Scatter':temp,'R':temp} )
		
		data_table = DataTable(source=table_source, columns=columns, width= 700, height=400)

		if var == 'xair':
			prec = '4'
		else:
			prec = '2'

		txt = PreText(text='The table must be displayed with increasing # before you start a new selection',width = 1000)

		count = 0
		for site in lat_ordered_sites:
			freq_sources[site]={site:{},SELECT:{}}
			freq_sources[site][site] = ColumnDataSource(data={'x':FREQ_DATA[site][site]['datetime'],'y':FREQ_DATA[site][site][var],'nicetime':[dat.strftime('%d-%m-%Y %H:%M') for dat in FREQ_DATA[site][site]['datetime']]})
			freq_sources[site][SELECT] = ColumnDataSource(data={'x':FREQ_DATA[site][SELECT]['datetime'],'y':FREQ_DATA[site][SELECT][var]})

			freq_cor_sources[site] = ColumnDataSource(data={'x':[],'y':[]})

			freq_sources[site][site].callback = CustomJS(args = dict(s2=freq_sources[site][SELECT],dt=data_table,scor=freq_cor_sources[site]), code="""
			var inds = cb_obj.get('selected')['1d'].indices;
			var d1 = cb_obj.get('data');
			var d2 = s2.get('data');
			var tab = dt.get('source').get('data');
			var dcor = scor.get('data');

			var difm = 0;
			var difm2 = 0;
			var scat = 0;

			var ym1 = 0;
			var ym2 = 0;

			var T1 = 0;
			var T2 = 0;
			var T3 = 0;

			dcor['x'] = [];
			dcor['y'] = [];

			tab['N']["""+str(count)+"""] = inds.length;


			if (inds.length == 0) {
				tab['RMS']["""+str(count)+"""] = 0;
				tab['Bias']["""+str(count)+"""] = 0;
				tab['Scatter']["""+str(count)+"""] = 0;
				tab['R']["""+str(count)+"""] = 0;
				dt.change.emit();
				return;
			}

			for (i=0; i < inds.length; i++){
				difm += d1['y'][inds[i]] - d2['y'][inds[i]];
				difm2 += Math.pow(d1['y'][inds[i]] - d2['y'][inds[i]],2);
				ym1 += d1['y'][inds[i]];
				ym2 += d2['y'][inds[i]];

				dcor['x'].push(d2['y'][inds[i]]);
				dcor['y'].push(d1['y'][inds[i]]);
			}

			difm /= inds.length;
			difm2 /= inds.length;
			ym1 /= inds.length;
			ym2 /= inds.length;

			for (i=0; i < inds.length; i++){
				scat += Math.pow(d1['y'][inds[i]] - d2['y'][inds[i]] - difm,2);
			}

			for (i=0; i < inds.length; i++){
				T1 += (d1['y'][inds[i]] - ym1)*(d2['y'][inds[i]] - ym2);
				T2 += Math.pow(d1['y'][inds[i]] - ym1,2);
				T3 += Math.pow(d2['y'][inds[i]] - ym2,2);
			}

			tab['RMS']["""+str(count)+"""] = Math.sqrt(difm2).toFixed("""+prec+""");
			tab['Bias']["""+str(count)+"""] = difm.toFixed("""+prec+""");
			tab['Scatter']["""+str(count)+"""] = Math.sqrt(scat/(inds.length -1)).toFixed("""+prec+""");
			tab['R']["""+str(count)+"""] = (T1/Math.sqrt(T2*T3)).toFixed("""+prec+""");

			dt.change.emit();
			scor.change.emit();
			""")
			count += 1

		min_x = min([FREQ_DATA[site][site]['datetime'][0] for site in FREQ_DATA])
		max_x = max([FREQ_DATA[site][site]['datetime'][-1] for site in FREQ_DATA])

		min_y = min([min(abs(FREQ_DATA[site][site][var])) for site in FREQ_DATA])
		max_y = max([max(FREQ_DATA[site][site][var]) for site in FREQ_DATA])

		ampli = max_y - min_y

		# we will set the y axis range at +/- 10% of the data max amplitude
		min_y = min_y - ampli*0.1/100
		max_y = max_y + ampli*0.1/100

		fig = figure(output_backend = "webgl", title = 'TCCON '+var+' vs '+'datetime', y_range=[min_y,max_y], plot_width = 900, plot_height = 650, tools = TOOLS, x_axis_type='datetime', x_range = Range1d(min_x,max_x))

		fig.tools[-2].dimensions='width' # only allow the box select tool to select data along the X axis (will select all Y data in a given X range)

		fig.tools[-2].callback = CustomJS(args=dict(txt=txt),code="""
			var sel = cb_data["geometry"];
			
			var startsec = sel["x0"]/1000;
			var start = new Date(0);

			start.setUTCSeconds(startsec)

			var startstring = ("0" + start.getDate()).slice(-2) + "-" + ("0"+(start.getMonth()+1)).slice(-2) + "-" +start.getFullYear() + " " + ("0" + start.getHours()).slice(-2) + ":" + ("0" + start.getMinutes()).slice(-2);

			var finishsec = sel["x1"]/1000;
			var finish = new Date(0);

			finish.setUTCSeconds(finishsec)

			var finishstring = ("0" + finish.getDate()).slice(-2) + "-" + ("0"+(finish.getMonth()+1)).slice(-2) + "-" +finish.getFullYear() + " " + ("0" + finish.getHours()).slice(-2) + ":" + ("0" + finish.getMinutes()).slice(-2);

			txt.text = 'Selection range from '+startstring + ' to ' + finishstring;

			txt.change.emit(); 
			""")

		plots=[]
		for site in lat_ordered_sites:
			plots.append( fig.scatter(x='x',y='y',color=DATA[site]['color'],alpha=0.5,source=freq_sources[site][site]) )
			plots.append( fig.scatter(x='x',y='y',color='black',alpha=0.5,source=freq_sources[site][SELECT]) )

		N_plots = range(len(plots))

		N_plots2 = range(len(plots)/2)

		even_plots = [i for i in N_plots if i%2 == 0]

		legend=Legend(items=[(SELECT,[plots[1]])]+[(lat_ordered_sites[i],[plots[even_plots[i]]]) for i in range(len(lat_ordered_sites))],location=(0,0))
		fig.add_layout(legend,'right')
		fig.yaxis.axis_label = var
		fig.xaxis.axis_label = 'Time'

		# correlation plot
		corfig = figure(output_backend = "webgl", title = 'Correlations', plot_width = 400, plot_height = 400, x_range = [min_y,max_y], y_range = [min_y,max_y]) 
		corfig.toolbar.logo = None
		corfig.toolbar_location = None
		corfig.xaxis.axis_label = ' '.join([SELECT,var])
		corfig.yaxis.axis_label = ' '.join(['sites',var])

		corplots = []
		linerange = list(np.arange(0,int(max_y),ceil(max_y)/10.0))+[max_y]
		for site in lat_ordered_sites:
			corplots.append( corfig.scatter(x='x',y='y',color=DATA[site]['color'],alpha=0.5,source=freq_cor_sources[site]) )
		
		corfig.line(x=linerange,y=linerange,color='black')

		N_corplots = range(len(corplots))

		checkbox = CheckboxGroup(labels=lat_ordered_sites,active=[],width=200)

		iterable = [('p'+str(i),plots[i]) for i in N_plots]+[('pcor'+str(i),corplots[i]) for i in N_corplots]+[('checkbox',checkbox)]

		checkbox_code = """var indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; };"""
		checkbox_code += ''.join(['p'+str(i)+'.visible = indexOf.call(checkbox.active, '+str(i/2)+') >= 0; p'+str(i+1)+'.visible= indexOf.call(checkbox.active, '+str(i/2)+') >= 0; pcor'+str(i/2)+'.visible = indexOf.call(checkbox.active, '+str(i/2)+') >= 0;' for i in range(0,len(N_plots),2)])
		checkbox.callback = CustomJS(args={key: value for key,value in iterable}, code=checkbox_code)

		clear_button = Button(label='Clear all',width=200)
		clear_button_code = """checkbox.set("active",[]);"""+checkbox_code
		clear_button.callback = CustomJS(args={key: value for key,value in iterable}, code=clear_button_code)

		check_button = Button(label='Check all',width=200)
		check_button_code = """checkbox.set("active","""+str(N_plots2)+""");"""+checkbox_code
		check_button.callback = CustomJS(args={key: value for key,value in iterable}, code=check_button_code)

		download_button = Button(label='Save Table to CSV', width = 200)
		download_button.callback = CustomJS(args=dict(dt=data_table),code="""
		var tab = dt.get('source').get('data');
		var filetext = 'site,latitude,N,RMS,Bias,Scatter,R'+String.fromCharCode(10);
		for (i=0; i < tab['site'].length; i++) {
		    var currRow = [tab['site'][i].toString(),
		                   tab['latitude'][i].toString(),
		                   tab['N'][i].toString(),
		                   tab['RMS'][i].toString(),
		                   tab['Bias'][i].toString(),
		                   tab['Scatter'][i].toString(),
		                   tab['R'][i].toString()+String.fromCharCode(10)];

		    var joined = currRow.join();
		    filetext = filetext.concat(joined);
		}

		var filename = 'data_result.csv';
		var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

	    var link = document.createElement("a");
	    link = document.createElement('a')
	    link.href = URL.createObjectURL(blob);
	    link.download = filename
	    link.target = "_blank";
	    link.style.visibility = 'hidden';
	    link.dispatchEvent(new MouseEvent('click'))

		""")

		group = widgetbox(checkbox,clear_button,check_button)

		grid = gridplot([[fig,group],[download_button,txt],[data_table,corfig]],toolbar_location='left')

		tabs.append(Panel(child=grid,title=var))

	final=Tabs(tabs=tabs)

	print('\n -writting FREQ_TCCON_'+'_'.join(FREQ.split())+'_'+SELECT+'.html'+' ...')
	outfile=open(os.path.join(save_path,'FREQ_TCCON_'+'_'.join(FREQ.split())+'_'+SELECT+'.html'),'w')
	outfile.write(file_html(final,CDN,'FREQ_TCCON'))
	outfile.close()
	print('FREQ_TCCON.html DONE in', time.time()-milestone,'seconds\n')

###########
# THE END #
###########

print('TCCON_comp.py DONE in',time.time()-codestart,'seconds')
