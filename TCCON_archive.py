####################
# Code description #
####################
"""
This code will produce an interactive plots reading into TCCON netcdf files to display time series of the data.

In the same directory as this code, make a 'TCCON' folder with netcdf files from TCCON (http://tccon.ornl.gov/ ~566MB for the public files as of 2017-08).
You can also edit the 'tccon_path' variable if you already have tccon files in a different location.

You can also use .eof.csv files, but not together with netcdf files; only one type of files should be in the 'TCCON' folder.
It is much slower to read from the .eof.csv files than it is to read from the netcdf files !
And unlike the netcdf files, using the date input widget won't make loading of data subsets faster with the .eof.csv files.

All the file names must start with the format xxYYYYMMDD_YYYYMMDD , xx is the two letters site abbreviation

This code needs to be run by a bokeh server like this:

bokeh serve --show TCCON_archive.py --args A

Where A is an argument that should be equal to 'simple' or 'comp' (without quotes), this will produce a specific layout for each:

- 'simple': There will be just one plot to display time series of any variables
- 'comp' : there will be three plots, two figures for time series, and a third figure that will plot the y axis of the two other figures against each other

If the selected site's data is in different files it will take longer to load new variables.

The program will create a cache_dic.npy file in which it will save full time series of variables that correspond to previous inputs.
The size of this cache file will be kept under 200 MB. If you like you can increase that size by changing the 'cache_max_size' value (in bytes)

e.g.	Reading a new variable in 'comp' mode for Lamont for the entire time series can take several minutes when done for the first time.
		The same data read from the cache file will load in 20-50 seconds.
		For Lamont in 'comp' mode, a set of two variables (+time,flag,color, and spectrum arrays) is ~ 50 MB for the full time series
"""

import sys
import os
import netCDF4
from datetime import datetime, timedelta
import time
import calendar
import numpy as np
import pandas as pd
from functools import partial

import bokeh
from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, TextInput, Div, CustomJS, Button, TextInput, Select, HoverTool, BoxSelectTool, DataTable, TableColumn
from bokeh.layouts import gridplot, widgetbox

# to ignore warnings (like for the save_source which is just meant to "pass" python variables to a callback but will return warnings because of columns of different sizes)
import warnings
warnings.filterwarnings('ignore')

mode = sys.argv[1]
if mode not in ['simple','comp']:
	print 'bokeh serve --show TCCON_archive.py --args A'
	print 'A is either "simple" or "comp"'
	sys.exit()

## MODIFIABLE SECTION
tccon_path = os.path.join(os.getcwd(),'full_archive') ## this is the only line that may need editing; full path to the folder containing the tccon netcdf files

cache_max_size = 2E8 # maximum size of the cache file (in bytes), any new cached data after that will remove the oldest data following certain rules (see add_cache function)

# if you modify the plotting colors, you will need to remove your cache_dic.npy file
main_color = 'goldenrod' # this will be the color used for the flag=0 data
flagged_color = 'grey' # this will be the color used for the flag!=0 data
hover_color = 'red' # this will be the color used for data hovered by the mouse

# this is the mode of the hovertool for the Figure 1 and Figure 2 in 'comp' mode
# 'width' to have all y data select over a x range (you can only pan the box selection left-right)
# 'height' can only pan the box selection up-down and it selects all the x data
# 'both' you can select a specific rectangle of data
boxselecttool_dimensions = 'both'

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
			# the ones below are not TCCON sites; I made them up for Toronto's EM27s
			'ta':'ta',
			'tb':'tb',
			'nn':'nn',
			'dn':'dn',
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
			'sp':'Norway',
			'et':'Canada',
			'an':'South Korea',
			'bu':'Philippines',
			'we':'Germany',
			# the ones below are not TCCON sites; I made them up for Toronto's EM27s
			'ta':'',
			'tb':'',
			'nn':'',
			'dn':'',
		 }
## END OF MODIFIABLE SECTION
#########################################################################################################################################################################
## SETUP SECTION
T_site = {key:'https://tccon-wiki.caltech.edu/Sites/' for key in T_FULL} # dictionary mapping each site prefix to its webpage
for key in T_FULL:
	if any(char.isdigit() for char in T_FULL[key]):
		T_site[key] += '_'.join(T_FULL[key].split()[:-1])
	else:
		T_site[key] += '_'.join(T_FULL[key].split())

netcdf = True
tccon_file_list = [i for i in os.listdir(tccon_path) if '.nc' in i] # list of the files in the 'TCCON' folder
if len(tccon_file_list) == 0:
	netcdf = False
	tccon_file_list = [i for i in os.listdir(tccon_path) if '.eof.csv' in i] # list of the files in the 'TCCON' folder

# list of TCCON 2 letters abbreviations from the files in the 'TCCON' folder doing list(set(a)) prevents repeated elements in the final list
prefix_list = list(set([i[:2] for i in tccon_file_list])) 

# determine if the files are from the public or private archive
public = True
if netcdf:
	f = netCDF4.Dataset(os.path.join(tccon_path,tccon_file_list[0]),'r')
	if 'flag' in [var for var in f.variables]:
		public = False
	f.close()
else:
	public = False

TOOLS = "box_zoom,wheel_zoom,pan,box_select,redo,undo,hover,reset" # the tools that will be available in the figure's toolbar

plot_width = 800
if mode == 'simple':
	plot_height = 450
elif mode == 'comp':
	plot_height = 250

fig = figure(output_backend="webgl",plot_width=plot_width,plot_height=plot_height,x_axis_type='datetime',tools=TOOLS,active_inspect=[],active_drag="box_zoom") # Figure 1

if mode == 'comp':
	fig2 = figure(output_backend="webgl",plot_width=plot_width,plot_height=plot_height,x_axis_type='datetime',x_range=fig.x_range,tools=TOOLS,active_inspect=[],active_drag="box_zoom") # Figure 2
	fig2.xaxis[0].axis_label = 'Time'
	source = ColumnDataSource(data={'x':[],'y1':[],'y2':[],'colo':[]}) # this is the object that stores the data to be plotted; it will be filled during callbacks
	fig2.scatter(x='x',y='y2',color='colo',hover_color=hover_color,alpha=0.7,source=source) # this is the plot of the bottom figure, it will read data from the 'source' object
	fig2.select_one(BoxSelectTool).dimensions = boxselecttool_dimensions

	# Figure 3, will plot y axis of Figure 1 vs y axis of Figure 2
	fig3 = figure(output_backend="webgl",plot_width=400,plot_height=400,tools=TOOLS,active_inspect=[],active_drag="box_zoom")
	fig3.scatter(x='y2',y='y1',color='colo',hover_color=hover_color,alpha=0.7,source=source)

	var_input = Select(title='Figure 1 variable:',width=220) # dropdown widget to select the variable to plot in Figure 1
	var_input2 = Select(title='Figure 2 variable:',width=220) # dropdown widget to select the variable to plot in Figure 2
elif mode == 'simple':
	source = ColumnDataSource(data={'x':[],'y1':[],'colo':[]}) # this is the object that stores the data to be plotted; it will be filled during callbacks
	fig.xaxis[0].axis_label = 'Time'
	var_input = Select(title='Variable to plot:',width=220) # dropdown widget to select the variable to plot in Figure 1

fig.select_one(BoxSelectTool).dimensions = boxselecttool_dimensions
fig.scatter(x='x',y='y1',color='colo',hover_color=hover_color,alpha=0.7,source=source) # this is the plot of the top figure, it will read data from the 'source' object

# list of sites that will be available in the site_input selection widget; I add an empty string at the beginning so that the first site can be selected right away
ordered_site_list = ['']+sorted([T_FULL[i] for i in prefix_list]) 

site_input = Select(title='Site:',options = ordered_site_list,width=220) # dropdown widget to select the TCCON site
date_input = TextInput(title='start-end yyyymmdd:') # text input widget to specify dates between which data should be fetched

# Information text for a Div widget
notes = """
<font size=4><b>Notes:</b></font></br>
</br>
<font size=2>
Use the dropdown buttons to select a site and a variable</br>
</br>
Use the text input to specify dates between which data will be fetched</br>
The second date is optional (e.g. 20120101-20140101 or 20120101)</br>
</br>
You can explore the plotted data using the toolbar</br>
</br>
If axis labels do not update, use the "Reset" tool</br>
</br>
Links: <a href='https://tccon-wiki.caltech.edu'>TCCON</a></font>
"""

notes_div = Div(text=notes,width=600) # this will display the 'notes' above
load_div = Div(text='Select a site',width=400) # the text of this widget will be updated with information on the state of callbacks
select_text = Div(text='',width = 450) # text div that will be updated with the selected range of date within the BoxSelect tool
dumdiv = Div(text='',height=10) # dummy empty Div widget for spacing

skip_list = ['_Version','ak_','prio','checksum','graw','spectrum','year','ada'] # variables including those keywords won't be shown in the 'var_input' dropdown

# javascript code for a dummy (invisible) button, it starts and stops a timer that will be displayed in the 'load_div' widget
dum_button_code = """
if (cb_obj.button_type.includes('success')){
var start = new Date();	
var intervalID = setInterval(function(){var current = new Date(); var diff=((current-start)/1000.0).toFixed(1); load_div.text='Loading data: '+diff.toString()+' s';	}, 100)
cb_obj.button_type = 'danger';
} else {
var noIntervals = setInterval(function(){});
for (var i = 0; i<noIntervals; i++) { window.clearInterval(i);}
load_div.text='Data loaded';
cb_obj.button_type = 'success';
}
"""
dum_button = Button(label='dummy',button_type='success',width=200) # the dummy button itself
dum_button.callback = CustomJS(args={'load_div':load_div},code=dum_button_code) # the callback of the button

# this code will trigger when the source data is changed in order to update the data table that shows correlations of selected variables
correlation_code = """
var inds = source.selected['1d'].indices;
var data = source.data;
var tab = dt.source.data;

var ym1 = 0;
var ym2 = 0;

var T1 = 0;
var T2 = 0;
var T3 = 0;

tab['N'][0] = inds.length;

if (inds.length == 0) {
	tab['R'][0] = 0;
	dt.change.emit();
	return;
}

for (i=0; i < inds.length; i++){
	ym1 += data['y1'][inds[i]];
	ym2 += data['y2'][inds[i]];
}

ym1 /= inds.length;
ym2 /= inds.length;

for (i=0; i < inds.length; i++){
	T1 += (data['y1'][inds[i]] - ym1)*(data['y2'][inds[i]] - ym2);
	T2 += Math.pow(data['y1'][inds[i]] - ym1,2);
	T3 += Math.pow(data['y2'][inds[i]] - ym2,2);
}

tab['R'][0] = (T1/Math.sqrt(T2*T3)).toFixed(3);

dt.change.emit();
"""

# this code will trigger when the BoxSelectTool is used in 'comp' mode; it fills the data_table based on the selected data
# it also triggers the 'center_button' after the selection if the files are not public
box_select_code = """
var sel = cb_data["geometry"];

var startsec = sel["x0"]/1000;
var start = new Date(0);

start.setUTCSeconds(startsec)

var startstring = ("0" + start.getUTCDate()).slice(-2) + "-" + ("0"+(start.getUTCMonth()+1)).slice(-2) + "-" +start.getUTCFullYear() + " " + ("0" + start.getUTCHours()).slice(-2) + ":" + ("0" + start.getUTCMinutes()).slice(-2);

var finishsec = sel["x1"]/1000;
var finish = new Date(0);

finish.setUTCSeconds(finishsec)

var finishstring = ("0" + finish.getUTCDate()).slice(-2) + "-" + ("0"+(finish.getUTCMonth()+1)).slice(-2) + "-" +finish.getUTCFullYear() + " " + ("0" + finish.getUTCHours()).slice(-2) + ":" + ("0" + finish.getUTCMinutes()).slice(-2);

txt.text = 'Selection range from '+startstring + ' to ' + finishstring;

txt.change.emit(); 


setTimeout(function(){
	button_list = document.getElementsByTagName('button');
	
	for(i=0;i<button_list.length;i++){
		if(button_list[i].textContent.includes("Scale")){button_list[i].click()}
	}

	}, 1000); // I need the timeout because it sometimes misses without it
"""

# in 'comp' mode, each of the 3 plots has a different hover tool, this button can click them all at once for convenience.
hover_button_code="""
if(cb_obj.button_type.includes("success")){
cb_obj.button_type = 'danger';
cb_obj.label = 'Disable hover tools'
} else {
	
cb_obj.button_type = 'success';
cb_obj.label= 'Enable hover tools';
}

for(i=7;i<10;i++){document.getElementsByClassName('bk-toolbar-button')[i].click()}
"""
hover_button = Button(label='Enable hover',button_type='success',width=220)
hover_button.callback = CustomJS(code=hover_button_code)

dum_text = TextInput() # dummy text input widget; it will be used to trigger the dummy button callback when the dropdown widgets are used.
# callback of the dummy 'dum_text' TextInput widget to trigger a button click on the dummy button 'dum_button'
dum_text_code = """
button_list = document.getElementsByTagName('button');

for(i=0;i<button_list.length;i++){
	if(button_list[i].textContent.includes("dummy")){button_list[i].click()}
}
"""
dum_text.js_on_change('value',CustomJS(code=dum_text_code))

# if the 'timeout' callback is set too early, the dummy widgets won't be hidden.
# but if it is set too late, the user will have time to see the dummy widgets before they disappear.
dum_alert = Div(text='<b>OOPS !</b> the widgets below are not supposed to be hidden. Try refreshing the page',width=700)

# use a different cache file for each data type
if netcdf and public:
	cache_file = 'cache_dic_pub.npy'
elif netcdf:
	cache_file = 'cache_dic.npy'
else:
	cache_file = 'cache_dic_eof.npy'

# global variables to make the data reading faster, it will save entire time series corresponding to different inputs; size will be limited, see the add_cache function
try:
	cache_dic = np.load(cache_file).item()
except IOError:
	cache_dic = {}
	np.save(cache_file,cache_dic)
else:
	print 'Cache dictionnary found:\nIf you added or modified files since last time, remove cache_dic.npy and run the program again, load times will initially be longer'

# global variables to keep track of errors in the load_var function
broken_text = ''
broken = False 

## END OF SETUP SECTION
#########################################################################################################################################################################
## ADD_CACHE FUNCTION
def add_cache(date_val,site,datasource,max_size=cache_max_size):
	'''
	This function will add data to the cache_dic dictionnary in order to make load time shorter.
	It will keep the size of the cache file below 'max_size' (in bytes)
	Whenever that size is exceeded, it will remove the oldest cached data following certain priority rules.
	'''

	global cache_dic

	stamp = datetime.now()

	keep_loop = True
	it = 1
	while keep_loop:
		if os.path.getsize(cache_file)<max_size:
			try:
				type(cache_dic[date_val])
			except KeyError:
				print 'Caching new date input:',date_val
				cache_dic[date_val] = {}
				cache_dic[date_val]['stamp'] = stamp

			try:
				type(cache_dic[date_val][site])
			except KeyError:
				print 'Caching new site input:',date_val,site
				cache_dic[date_val][site] = {}
				cache_dic[date_val][site]['stamp'] = stamp

			try:
				type(cache_dic[date_val][site]['x'])
			except KeyError:
				print 'Caching new time data:',date_val,site,len(datasource['x']),'values'
				cache_dic[date_val][site]['x'] = datasource['x']

			for var in datasource:
				if var!='x':
					nice_var = var
					if var=='y1':
						nice_var = var_input.value
					elif var=='y2':
						nice_var = var_input2.value

					try:
						type(cache_dic[date_val][site][nice_var])
					except KeyError:
						print 'Caching new variable:',date_val,site,nice_var,len(datasource[var]),'values'
						cache_dic[date_val][site][nice_var] = {}
						cache_dic[date_val][site][nice_var]['stamp'] = stamp
						cache_dic[date_val][site][nice_var]['value'] = datasource[var]
		else:
			print 
			while os.path.getsize(cache_dic)>max_size:
				# if there are more than 5 date_val, remove the oldest date_val
				if len(cache_dic)>5:
					del_date = min([(cache_dic[date_val]['stamp'],date_val) for date_val in cache_dic])[1]
					del cache_dic[del_date]
					print 'Removing',del_date,'from cache'
				# if there are more than 5 sites, remove the oldest site
				site_list = list(set([[site for site in cache_dic[date_val] if site!='stamp'] for date_val in cache_dic][0])) # list(set()) removes duplicates
				if len(site_list)>5:
					del_site = min([[(cache_dic[date_val][site]['stamp'],(date_val,site)) for site in cache_dic[date_val] if site!='stamp'] for date_val in cache_dic][0])[1]
					del cache_dic[del_site[0]][del_site[1]]
					print 'Removing',del_site[0],del_site[1],'from cache'
				# if there are more than 20 variables, remove the oldest variable (except for times, flags, colors, and spectrum names; those will only be removed when a site or date_val is removed)
				var_list = list(set([[[var for var in cache_dic[date_val][site] if var not in ['x','flag','spectrum','colo','stamp']] for site in cache_dic[date_val] if site!='stamp'] for date_val in cache_dic][0][0]))
				if len(var_list)>20:
					del_var = min([[[(cache_dic[date_val][site][var]['stamp'],(date_val,site,cache_dic[date_val][site][var]['value'])) for var in cache_dic[date_val][site] if var not in ['x','flag','spectrum','colo','stamp']] for site in cache_dic[date_val] if site !='stamp'] for date_val in cache_dic][0][0])
					del cache_dic[del_var[0]][del_var[1]][del_var[2]]
					print 'Removing',del_var[0],del_var[1],del_var[2],'from cache'

		if os.path.getsize(cache_file)<max_size:
			keep_loop = False
		else:
			print 'Cache is over',max_size/1E6,'MB after caching iteration',it
			print 'If this does not stop iterating it means that the current selection alone is over the max_size'
			it+=1
	print 'Saving new cache file ...'
	np.save(cache_file,cache_dic)
	print cache_file,'size:',os.path.getsize(cache_file)/1E6,'MB'
## END OF ADD_CACHE FUNCTION
#########################################################################################################################################################################
## INITIALIZE FUNCTION
def initialize(all_var):
	'''
	Function that fills the var_input widgets with options and resets the data source
	Also sets up the hovertool tooltips according to the file type and mode
	'''
	if len(var_input.options)==0:
		var_list = sorted([var for var in all_var if True not in [elem in var for elem in skip_list]])
		var_input.options = var_list # fill the 'var_input' dropdown with options
		if mode == 'comp':
			var_input2.options = var_list # fill the 'var_input2' dropdown with options

	# configure the source and HoverTool based on the type of file
	if public and mode=='simple':
		source.data.update({'x':[],'y1':[],'colo':[]})
		fig.select_one(HoverTool).tooltips = [(fig.yaxis[0].axis_label,'@y1')]
	elif not public and mode=='simple':
		source.data.update({'x':[],'y1':[],'colo':[],'flag':[],'spectrum':[]})
		fig.select_one(HoverTool).tooltips = [('value','@y1'),('spectrum','@spectrum'),('flag','@flag')]
	elif public and mode=='comp':
		source.data.update({'x':[],'y1':[],'y2':[],'colo':[]})
		fig.select_one(HoverTool).tooltips = [(fig.yaxis[0].axis_label,'@y1'),(fig2.yaxis[0].axis_label,'@y2')]
		fig2.select_one(HoverTool).tooltips = [(fig2.yaxis[0].axis_label,'@y2'),(fig.yaxis[0].axis_label,'@y1')]	
		fig3.select_one(HoverTool).tooltips = [('y','@y1'),('x','@y2')]
	elif not public and mode=='comp':
		source.data.update({'x':[],'y1':[],'y2':[],'colo':[],'flag':[],'spectrum':[]})
		fig.select_one(HoverTool).tooltips = [(fig.yaxis[0].axis_label,'@y1'),('spectrum','@spectrum'),('flag','@flag'),(fig2.yaxis[0].axis_label,'@y2')]
		fig2.select_one(HoverTool).tooltips = [(fig2.yaxis[0].axis_label,'@y2'),('spectrum','@spectrum'),('flag','@flag'),(fig.yaxis[0].axis_label,'@y1')]	
		fig3.select_one(HoverTool).tooltips = [('y','@y1'),('x','@y2'),('spectrum','@spectrum'),('flag','@flag')]
## END OF INITIALIZE FUNCTION
#########################################################################################################################################################################
## LOAD_VAR FUNCTION
def load_var(site_file_list,site):
	'''
	Function called by set_site() and set_var() to load the variables matching site,variable, and date inputs
	'''

	global broken_text, broken, cache_dic, all_var

	date_val = date_input.value # can be of the form 'firstdate'-'lastdate' or just 'firstdate'
	site = site_input.value

	dum_text.value = str(time.time()) # click the timer button to start the loading countdown in the 'load_div' widget

	# check if there already is cached data that matches the inputs
	no_cached_data = False
	try:
		type(cache_dic[date_val][site][var_input.value])
	 	if mode == 'comp':
	 		type(cache_dic[date_val][site][var_input2.value])
	except KeyError:
		no_cached_data = True

	if no_cached_data:
		# loop over the TCCON files for the selected site
		broken = False # boolean check that will be set to True before breaking out of the loop below
		for filenum,site_file in enumerate(site_file_list):
			if netcdf:
				f = netCDF4.Dataset(os.path.join(tccon_path,site_file),'r') # netcdf file reader
				all_var = [var for var in f.variables]
			else:
				df = pd.read_csv(os.path.join(tccon_path,site_file),header=2) # read the .eof.csv file
				all_var = list(df)

			# setup some initializations if it is the first file
			if filenum==0:
				initialize(all_var) # fills variable inputs with options,resets the data source, and setup hovertool tooltips

			if var_input.value != '':	# only tries to read variables if a 'var_input' option has been selected
				if netcdf:
					nctime =  f.variables['time'][:] # fractional days since 1970
				else:			
					# for .eof.csv files convert year,day,hour in fractional days since 1970
					nctime = np.array([(datetime(int(df.year[i]),1,1)+timedelta(days=df.day[i]-1)-datetime(1970,1,1)).days+df.hour[i]/24 for i in range(len(df))])

				# use the value of the 'date_input' widget to determine the range of dates over which data should be fetched
				try: # for the minimum of the range
					mindate = date_input.value.split('-')[0] # if the 'date_input' widget is empty, this will raise an exception
					mindate = calendar.timegm(datetime(int(mindate[:4]),int(mindate[4:6]),int(mindate[6:8])).timetuple())/24/3600 # if the date is entered wrong, this will raise an exception
				except: # catch all exceptions
					mindate = calendar.timegm(datetime(1970,1,1).timetuple())/24/3600 # if an exception has been caught, just use a very early date

				try: # for the maximum of the range 
					maxdate = date_input.value.split('-')[1] # if the 'date_input' widget is empty, or if only the minimum date is given, this will raise an exception
					maxdate = calendar.timegm(datetime(int(maxdate[:4]),int(maxdate[4:6]),int(maxdate[6:8])).timetuple())/24/3600 # if the date is entered wrong, this will raise an exception
				except: # catch all exceptions
					maxdate = calendar.timegm(datetime(2050,1,1).timetuple())/24/3600 # if an exception has been caught, just use a very late date

				# perform a quick check on date ranges based on file names to avoid looping over files for nothing
				if filenum==0:
					all_file_min = site_file_list[0][2:10]		# minimum YYYYMMDD of all files
					all_file_max = site_file_list[-1][11:19]	# maximum YYYYMMDD of all files
					all_file_min_date = calendar.timegm(datetime(int(all_file_min[:4]),int(all_file_min[4:6]),int(all_file_min[6:8])).timetuple())/24/3600
					all_file_max_date = calendar.timegm(datetime(int(all_file_max[:4]),int(all_file_max[4:6]),int(all_file_max[6:8])).timetuple())/24/3600
					# break out of the for loop if min or max dates from the file names are not compatible with the date input
					if (mindate>all_file_max_date) or (maxdate<all_file_min_date):
						broken_text = site +' date range: '+all_file_min+'-'+all_file_max
						broken = True
					if mindate>maxdate:
						broken_text = 'Wrong date input'
						broken = True
					if broken:
						if netcdf:
							f.close() # close the netcdf reader
						break

				# for each file, fastforward to next iteration of the for loop if the dates in the file name are not compatible with the date input
				file_min = site_file[2:10]
				file_max = site_file[11:19]
				file_min_date = calendar.timegm(datetime(int(file_min[:4]),int(file_min[4:6]),int(file_min[6:8])).timetuple())/24/3600
				file_max_date = calendar.timegm(datetime(int(file_max[:4]),int(file_max[4:6]),int(file_max[6:8])).timetuple())/24/3600
				if (mindate>file_max_date) or (maxdate<file_min_date):
					if netcdf:
						f.close() # close the netcdf reader
					continue

				newtime = nctime[(nctime>=mindate) & (nctime<maxdate)] # list of times that satisfy the 'date_input' value (still fractional days since 1970)

				# check that there is actual data in the time range, if not, go to the next (next iteration in for loop)
				if len(newtime)==0:
					continue

				start_id = np.where(nctime==newtime[0])[0][0] # ID of the starting time in the full time list of the file
				end_id = np.where(nctime==newtime[-1])[0][0]+1 # ID of the end time in the full time list of the file

				try: # attempt to fetch data using start_id and end_id
					add_x = np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in nctime[start_id:end_id]])
					if netcdf:
						add_y1 = f.variables[var_input.value][start_id:end_id]
						if mode == 'comp':
							add_y2 = f.variables[var_input2.value][start_id:end_id]
					else:
						add_y1 = np.array(df[var_input.value][start_id:end_id])
						if mode == 'comp':
							add_y2 = np.array(df[var_input2.value][start_id:end_id])
					if public:
						add_colo = np.array([main_color]*len(add_x))
					else:
						if netcdf:
							add_colo = np.array([main_color if int(elem)==0 else flagged_color for elem in f.variables['flag'][start_id:end_id]]) # data with flag=0 is red; data with flag!=0 is grey
							add_spectrum = np.array([''.join(elem) for elem in f.variables['spectrum'][start_id:end_id]])# name of spectra for the HoverTool
							add_flag = f.variables['flag'][start_id:end_id].astype(int) # flags for the HoverTool
						else:
							add_colo = np.array([main_color if int(elem)==0 else flagged_color for elem in df['flag'][start_id:end_id]])
							add_spectrum = np.array([elem for elem in df['spectrum'][start_id:end_id]])
							add_flag = np.array([' '.join([str(flag),all_var[flag],'=',str(df[all_var[flag]][start_id:end_id][ID])]) for ID,flag in enumerate(df['flag'][start_id:end_id])])
				except:
					pass # do nothing if any exception is raised (potentially index errors if the files have columns of inconsistent lengths)
				else: # if the try didnt raise any exceptions, update the 'source' of the plots with new data
					if public and mode=='simple':
						source.data.update({	'x' : np.append(source.data['x'],add_x),
												'y1' : np.append(source.data['y1'],add_y1),
												'colo' : np.append(source.data['colo'],add_colo),})
					elif not public and mode=='simple':
						source.data.update({	'x' : np.append(source.data['x'],add_x),
												'y1' : np.append(source.data['y1'],add_y1),
												'colo' : np.append(source.data['colo'],add_colo),
												'flag' : np.append(source.data['flag'],add_flag),
												'spectrum': np.append(source.data['spectrum'],add_spectrum),})
					elif public and mode=='comp':
						source.data.update({	'x' : np.append(source.data['x'],add_x),
												'y1' : np.append(source.data['y1'],add_y1),
												'y2' : np.append(source.data['y2'],add_y2),
												'colo' : np.append(source.data['colo'],add_colo),})
					elif not public and mode=='comp':
						source.data.update({	'x' : np.append(source.data['x'],add_x),
												'y1' : np.append(source.data['y1'],add_y1),
												'y2' : np.append(source.data['y2'],add_y2),
												'colo' : np.append(source.data['colo'],add_colo),
												'flag' : np.append(source.data['flag'],add_flag),
												'spectrum': np.append(source.data['spectrum'],add_spectrum),})
					print site_file,'data streamed'
			if netcdf:
				f.close() # close the netcdf reader		
		else: # else clause of the for loop, this will execute if no 'break' was encountered
			if len(source.data['x'])==0: # no data
				if date_val=='': # date_val empty mean the whole time range of the site has been evaluated
					load_div.text = site+' has no data' # if you see this one then there is a problem with the netcdf file ...
					if mode == 'comp':
						load_div.text = 'Select a variable for Figure 2'
				elif len(date_val)==8: # if only 'firstdate' is given
					load_div.text = site+' has no data after '+date_val
				elif len(date_val)>8: # if 'firstdate-lastdate' is given
					load_div.text = site+' has no data for '+date_val
		if broken:
			load_div.text = broken_text

		if len(source.data['y1'])!=0:
			add_cache(date_val,site,source.data)
	
	else: #else clause of 'if no_cached_data'
		initialize(all_var) # fills variable inputs with options,resets the data source, and setup hovertool tooltips
		print 'Using cached data for',date_val,site
		source.data.update({	'x' : cache_dic[date_val][site]['x'],
								'y1' : cache_dic[date_val][site][var_input.value]['value'],
								'y2' : cache_dic[date_val][site][var_input2.value]['value'],
								'colo' : cache_dic[date_val][site]['colo']['value'],
								'flag' : cache_dic[date_val][site]['flag']['value'],
								'spectrum': cache_dic[date_val][site]['spectrum']['value'],})
	dum_text.value = str(time.time()) # click the timer button again to end the loading countdown
	load_div.text = 'Data loaded'
	print 'load_var() DONE: if all the data is not showing quickly, you should use the date_input widget to select a smaller subset of data'
## END OF LOAD_VAR FUNCTION
#########################################################################################################################################################################
## SET_SITE FUNCTION
def set_site(attr,old,new):
	'''
	callback for the 'site_input' dropdown widget.
	If no variable has been selected before, this only fills the 'var_input' dropdown widget with options.
	Otherwise it reads and displays data based on the values of all the input widgets
	'''

	global cache_dic

	site = site_input.value # the selected TCCON site
	prefix = [key for key in T_FULL if T_FULL[key]==site][0] # TCCON 2 letters abbreviation of the site
	site_file_list = [i for i in tccon_file_list if prefix in i] # a list of the associated files
	fig.yaxis[0].axis_label = var_input.value
	if mode == 'comp':
		fig2.yaxis[0].axis_label = var_input2.value
		fig3.yaxis[0].axis_label = var_input.value
		fig3.xaxis[0].axis_label = var_input2.value 
	fig.title.text = T_FULL[prefix]+', '+T_LOC[prefix] # site name + site location
	notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix]+"'>"+site+"</a>") # updated the information widget with a link to the site's webpage

	# now check if there are options in the 'var_input' widget
	no_select_var = False
	if len(var_input.options)==0:
		no_select_var = True

	load_var(site_file_list,site)

	if no_select_var:
		load_div.text = 'Select a variable for Figure 1'
## END OF SET_SITE FUNCTION
#########################################################################################################################################################################
## SET_VAR FUNCTION
def set_var(attr,old,new):
	'''
	callback for the 'var_input' dropdown widget.
	It reads and displays data based on the values of all the input widgets
	'''

	global broken, cache_dic

	date_val = date_input.value
	if date_val not in cache_dic:
		cache_dic[date_val] = {}

	try: # attempt to get a variable
		site = site_input.value # the selected TCCON site
		prefix = [key for key in T_FULL if T_FULL[key]==site][0] # TCCON 2 letters abbreviation of the site
		site_file_list = [i for i in tccon_file_list if prefix in i] # a list of the associated files
	except: 
		pass # do nothing if an exception occurs (not sure if there can even be one here; maybe an empty 'site_input')
	else: # if no exception occured ...
		fig.yaxis[0].axis_label = var_input.value
		if mode == 'comp':
				fig2.yaxis[0].axis_label = var_input2.value
				fig3.yaxis[0].axis_label = var_input.value
				fig3.xaxis[0].axis_label = var_input2.value
		
		load_var(site_file_list,site)
## END OF SET_VAR FUNCTION
#########################################################################################################################################################################
## INPUT WIDGETS SECTION
# assign the python callbacks to the input widgets
site_input.on_change('value',set_site)
var_input.on_change('value',set_var)
date_input.on_change('value',set_var)
if mode == 'comp':
	var_input2.on_change('value',set_var)
	
	# widgets specific to the comparison mode
	table_source = ColumnDataSource( data = {'N':[0],'R':[0]} ) # the data source of the table
	data_table = DataTable(source=table_source, columns=[ TableColumn(field='N',title='N'),TableColumn(field='R',title='R'),], width=200, height=55)
	select_text = Div(text='',width = 450) # text div that will be updated with the selected range of date within the BoxSelect tool
	
	# assign JS callbacks to the source and box selection tool
	source.js_on_change('data', CustomJS(args={'source':source,'dt':data_table},code=correlation_code)) # use the dummy widget because source.js_on_change('data') does not work
	source.js_on_change('selected', CustomJS(args={'source':source,'dt':data_table},code=correlation_code))
	fig.select_one(BoxSelectTool).callback = CustomJS(args={'txt':select_text},code = box_select_code)
	fig2.select_one(BoxSelectTool).callback = CustomJS(args={'txt':select_text},code = box_select_code)
	fig3.select_one(BoxSelectTool).callback = CustomJS(code = """
	setTimeout(function(){
		button_list = document.getElementsByTagName('button');
	
		for(i=0;i<button_list.length;i++){
			if(button_list[i].textContent.includes("Scale")){button_list[i].click()}
		}
	}, 1000); // I need the timeout because it sometimes misses without it""")
## END OF INPUT WIDGETS SECTION
#########################################################################################################################################################################
## CENTER FUNCTION
def center():
	"""
	callback to update the y axis range based on the data.
	the range is by default auto adjusted to show all the data even huge outliers.
	This function takes the flag = 0 data (or all the data if it's from public files) and scales the y axis based on the mean of that data subset

	In short it will zoom in by disregarding outliers.

	This callback is triggered by clicks on the 'center_button'.
	At the end of the BoxSelectTool callback, the 'center_button' is automatically clicked.
	"""

	data_list = [np.array(source.data['y1'])[np.array(source.data['colo'])==main_color].astype(np.float)] # data in red is all the data from public files and flag=0 data from private files
	if mode == 'comp':
		data_list += [source.data['y2'][source.data['colo']==main_color].astype(np.float)]

	for ID,main_data in enumerate(data_list):
		negatives = main_data[main_data<0]
		positives = main_data[main_data>0]

		mean_main_data = np.mean(main_data)

		try:
			if len(negatives)==0:
				main_data = main_data[(main_data<1.5*mean_main_data) & (main_data>0.5*mean_main_data)] # resample the data to get within +/- 50% of the mean
			elif len(positives)==0:
				main_data = main_data[(main_data>1.5*mean_main_data) & (main_data<0.5*mean_main_data)] # resample the data to get within +/- 50% of the mean
			elif 0.5<(len(positives)/len(negatives))<1.5:
				main_data = main_data[(main_data<2*np.mean(positives)) & (main_data>2*np.mean(negatives))] # resample the data to get within +20% of the positive mean and +20% of the negative mean

			if list(main_data).count(main_data[0])!=len(main_data): # the variable is not a constant

				min_y = min(main_data)
				max_y = max(main_data)

				ampli = (max_y-min_y)/10

				min_y = min_y - ampli
				max_y = max_y + ampli

				if ID == 0:
					fig.y_range.start = min_y
					fig.y_range.end = max_y
					if mode == 'comp':
						fig3.y_range.start = min_y
						fig3.y_range.end = max_y
				if ID == 1: # only happens when mode==comp
					fig2.y_range.start = min_y
					fig2.y_range.end = max_y
					fig3.x_range.start = min_y
					fig3.x_range.end = max_y				

			else: # if the variable is a constant, update the 'load_div' to let the user know why nothing happened
				load_div.text = 'Variable is constant'
		except IndexError:
			load_div.text = 'IndexError: Variable is constant'
# END OF CENTER FUNCTION

center_button = Button(label='Scale without extrema',name='test',width=220) # button to 'center' the plot on the 'good' data
center_button.on_click(center) # assign the callback function to the button

# put the figure by itself in a grid layout (I can better control where the toolbar will show if i do that)
if mode == 'comp':
	figrid = gridplot([[fig]], toolbar_location = 'left')
	figrid2 = gridplot([[fig2]], toolbar_location = 'left')
	figrid3 = gridplot([[fig3]], toolbar_location = 'left')
	figrid.children[0].logo = None
	figrid2.children[0].logo = None
	figrid3.children[0].logo = None

	for cur_grid in [figrid,figrid2,figrid3]:
		cur_grid.children[0].tools = [i for i in cur_grid.children[0].tools if type(i)==bokeh.models.tools.HoverTool]

	if public:
		figroup = gridplot([[figrid],[figrid2],[figrid3,widgetbox(select_text,data_table,hover_button,name='bottom box')]],toolbar_location='above')
	else:
		figroup = gridplot([[figrid],[figrid2],[figrid3,widgetbox(select_text,data_table,center_button,hover_button,name='bottom box')]],toolbar_location='above')
elif mode == 'simple':
	figrid = gridplot([[fig]], toolbar_location = 'above')

# callback to hide the dummy widget box after the timeout callback
dum_hide_code = """
widgetbox_list = document.getElementsByClassName("bk-widget-box");

for(i=0;i<widgetbox_list.length;i++){
	if(widgetbox_list[i].textContent.includes('hidden')){widgetbox_list[i].style.display="none"}
}

toolbar_list = document.getElementsByClassName("bk-toolbar-box");
if(toolbar_list.length>1){
	for(i=0;i<toolbar_list.length;i++){
		if(toolbar_list[i].textContent.includes("Hover")){toolbar_list[i].style.display="none"}
	}
}
"""
dum_hide = TextInput() # dummy text input widget; it will be used in a 'timeout' callback after the page load in order to make the dummy widgets invisible
dum_hide.js_on_change('value',CustomJS(code=dum_hide_code))

if public:	
	# put the figure and all the widget in a single grid layout object
	if mode == 'simple':
		#final grid for the plot
		grid = gridplot([[figrid,widgetbox(site_input,var_input,dumdiv,date_input,load_div,notes_div,width=600,name='side box')],[widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=600,name='dummy box')]], toolbar_location = None)
	elif mode == 'comp':
		grid = gridplot([[figroup,widgetbox(site_input,var_input,var_input2,dumdiv,date_input,load_div,notes_div,width=600,name='side box')],[widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=600,name='dummy box')]], toolbar_location = None)

else: # for the 'private mode' only the index of the html divs containing the dummy widgets changes because of the addtion of the 'center_button' in the layout
	if mode == 'simple':
		grid = gridplot([[figrid,widgetbox(site_input,var_input,dumdiv,date_input,load_div,notes_div,width=600,name='side box')],[center_button],[widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=600,name='dummy box')]], toolbar_location = None)
	elif mode == 'comp':
		grid = gridplot([[figroup,widgetbox(site_input,var_input,var_input2,dumdiv,date_input,load_div,notes_div,width=600,name='side box')],[widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=600,name='dummy box')]], toolbar_location = None)

figrid.children[0].merge_tools = False # need to do that due to a bug in bokeh 0.12.6 that prevent gridplot to display the HoverTool in the toolbar
if mode == 'comp':
	figrid2.children[0].merge_tools = False
	figrid3.children[0].merge_tools = False

def hide_dummy():
	'''
	timeout callback that will be executed 1 second after the page is open.
	changes the value of the 'dum_hide' dummy TextInput widget, thise triggers the widget's javascript callback that makes all the dummy widgets invisible
	'''
	dum_hide.value = 'hide'

curdoc().title='TCCON' # this changes the title of the internet tab in which the document will be displayed
curdoc().add_root(grid) # this add the grid layout to the document
curdoc().add_timeout_callback(hide_dummy,1000) # this schedules the 'hide_dummy' callback to be triggered 1000 milliseconds after page load; if this fails to trigger too often, increase the number
