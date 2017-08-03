####################
# Code description #
####################
"""
This code will produce an interactive plots reading into TCCON netcdf files to display time series of the data.

In the same directory as this code, make a 'TCCON' folder with netcdf files from TCCON (http://tccon.ornl.gov/ ~566MB for the public files as of 2017-08)
You can also edit the 'tccon_path' variable if you already have tccon files in a different location.

This code needs to be run by a bokeh server like this:

bokeh serve --show TCCON_archive.py

If the selected site's data is in different files it will take longer to load new variables.
"""

import sys
import os
import netCDF4
from datetime import datetime
import time
import calendar
import numpy as np
from functools import partial

from bokeh.io import show,curdoc,save
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, TextInput, Div, CustomJS, Button, TextInput, Select, HoverTool
from bokeh.layouts import gridplot, widgetbox, LayoutDOM

tccon_path = os.path.join(os.getcwd(),'full_archive') ## this is the only line that may need editing; full path to the folder containing the tccon netcdf files

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
		 }

T_site = {key:'https://tccon-wiki.caltech.edu/Sites/' for key in T_FULL} # dictionary mapping each site prefixe to its webpage
for key in T_FULL:
	if any(char.isdigit() for char in T_FULL[key]):
		T_site[key] += '_'.join(T_FULL[key].split()[:-1])
	else:
		T_site[key] += '_'.join(T_FULL[key].split())

tccon_file_list = [i for i in os.listdir(tccon_path) if '.nc' in i] # list of the files in the 'TCCON' folder
prefix_list = list(set([i[:2] for i in tccon_file_list])) # list of TCCON 2 letters abbreviations from the files in the 'TCCON' folder doing list(set(a)) prevents repeated elements in the final list

# determine if the files are from the public or private archive
public = True
f = netCDF4.Dataset(os.path.join(tccon_path,tccon_file_list[0]),'r')
if 'flag' in [v for v in f.variables]:
	public = False
f.close()

source = ColumnDataSource(data={'x':[],'y':[],'colo':[]}) # this is the object that stores the data to be plotted; it will be filled during callbacks

TOOLS = "box_zoom,wheel_zoom,pan,redo,undo,hover,reset" # the tools that will be available in the figure's toolbar
fig = figure(output_backend="webgl",plot_width=600,plot_height=450,x_axis_type='datetime',tools=TOOLS,active_inspect=[],active_drag="box_zoom")
fig.xaxis.axis_label = 'Time'

fig.scatter(x='x',y='y',color='colo',alpha=0.7,source=source) # this is the plot, it will read data from the 'source' object

site_input = Select(title='Site:',options = sorted([T_FULL[i] for i in prefix_list]),width=220) # dropdown widget to select the TCCON site
var_input = Select(title='Variable to plot:',width=220) # dropdown widget to select the variable to plot 
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

notes_div = Div(text=notes,width=600)
load_div = Div(text='Select a site',width=400) # the text of this widget will be updated with information on the state of callbacks

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
dum_button = Button(button_type='success') # the dummy button itself
dum_button.callback = CustomJS(args={'load_div':load_div},code=dum_button_code) # the callback of the button

dum_text = TextInput() # dummy text input widget; it will be used to trigger the dummy button callback when the dropdown widgets are used.
dum_hide = TextInput() # dummy text input widget; it will be used in a 'timeout' callback after the page load in order to make the dummy widgets invisible
dum_alert = Div(text='<b>OOPS !</b> the widgets below are not supposed to be visible. Try refreshing the page',width=700) # if the 'timeout' callback is set too early, the dummy widgets won't be hidden

def set_site(attr,old,new):
	'''
	callback for the 'site_input' dropdown widget.
	If no variable has been selected before, this only fills the 'var_input' dropdown widget with options.
	Otherwise it reads and displays data based on the values of all the input widgets
	'''
	site = site_input.value # the selected TCCON site
	prefix = [key for key in T_FULL if T_FULL[key]==site][0] # TCCON 2 letters abbreviation of the site
	site_file_list = [i for i in tccon_file_list if prefix in i] # a list of the associated netcdf files
	fig.yaxis.axis_label = var_input.value 
	fig.title.text = T_FULL[prefix]+', '+T_LOC[prefix] # site name + site location
	notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix]+"'>"+site+"</a>") # updated the information widget with a link to the site's webpage
	# now check if there are options in the 'var_input' widget
	var_check = True 
	if len(var_input.options)!=0:
		dum_text.value = str(time.time()) # click the timer button to start the loading countdown in the 'load_div' widget
	else:
		var_check = False
		load_div.text = 'Select a variable'
	
	# loop over the TCCON files for the selected site
	data_check = False
	broken = False
	for filenum,site_file in enumerate(site_file_list):
		f = netCDF4.Dataset(os.path.join(tccon_path,site_file),'r') # netcdf file reader

		# setup some initializations if it is the first file
		if filenum==0:
			var_input.options = sorted([v for v in f.variables if True not in [elem in v for elem in skip_list]]) # fill the 'var_input' dropdown with options
			# configure the source and HoverTool based on the type of file
			if public:
				source.data.update({'x':[],'y':[],'colo':[]})
				fig.select_one(HoverTool).tooltips = [('value','@y')]
			else:
				source.data.update({'x':[],'y':[],'colo':[],'flag':[],'spectrum':[]})
				fig.select_one(HoverTool).tooltips = [('value','@y'),('spectrum','@spectrum'),('flag','@flag')]

		if var_input.value != '':	# only tries to read variables if a'var_input' option has been selected 
			nctime =  f.variables['time'][:] # fractional days since 1970

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
					f.close() # close the netcdf reader
					break

			# for each file, fastforward to next iteration of the for loop if the dates in the file name are not compatible with the date input
			file_min = site_file[2:10]
			file_max = site_file[11:19]
			file_min_date = calendar.timegm(datetime(int(file_max[:4]),int(file_max[4:6]),int(file_max[6:8])).timetuple())/24/3600
			file_max_date = calendar.timegm(datetime(int(file_max[:4]),int(file_max[4:6]),int(file_max[6:8])).timetuple())/24/3600
			if (mindate>file_max_date) or (maxdate<file_min_date):
				f.close() # close the netcdf reader
				continue

			newtime = nctime[(nctime>=mindate) & (nctime<maxdate)] # list of times that satisfy the 'date_input' value (still fractional days since 1970)

			# check that there is actual data in the time range, if not, go to the next (next iteration in for loop)
			if len(newtime)!=0: 
				data_check = True
			else:
				continue

			start_id = np.where(nctime==newtime[0])[0][0] # ID of the starting time in the full time list of the file
			end_id = np.where(nctime==newtime[-1])[0][0]+1 # ID of the end time in the full time list of the file

			try: # attempt to fetch data using start_id and end_id
				add_x = np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in newtime])
				add_y = f.variables[var_input.value][start_id:end_id]
				if public:
					add_colo = np.array(['red']*len(add_x))
				else:
					add_colo = np.array(['red' if int(elem)==0 else 'grey' for elem in f.variables['flag'][start_id:end_id]]) # data with flag=0 is red; data with flag!=0 is grey
					add_spectrum = [''.join(elem) for elem in f.variables['spectrum'][start_id:end_id]] # name of spectra for the HoverTool
					add_flag = f.variables['flag'][start_id:end_id] # flags for the HoverTool
			except:
				pass # do nothing if any exception is raised (potentially index errors if the files have columns of inconsistent lengths)
			else: # if the try didnt raise any exceptions, update the 'source' of the plots with new data
				if public:
					source.data.update({	'x' : np.append(source.data['x'],add_x),
											'y' : np.append(source.data['y'],add_y),
											'colo' : np.append(source.data['colo'],add_colo)})
				else:
					source.data.update({	'x' : np.append(source.data['x'],add_x),
											'y' : np.append(source.data['y'],add_y),
											'colo' : np.append(source.data['colo'],add_colo),
											'flag' : np.append(source.data['flag'],add_flag),
											'spectrum':np.append(source.data['spectrum'],add_spectrum),})
		else: # if 'var_input' is empty
			data_check = True

		f.close() # close the netcdf reader
	else: # else clause of the for loop, this will execute if no 'break' was encountered		
		if data_check is not True: # no data
			date_val = date_input.value # can be of the form 'firstdate'-'lastdate' or just 'firstdate'
			if date_val=='': # date_val empty mean the whole time range of the site has been evaluated
				load_div.text = site+' has no data' # if you see this one then there is a problem with the netcdf file ...
			elif len(date_val)==8: # if only 'firstdate' is given
				load_div.text = site+' has no data after '+date_val
			else: # if 'firstdate-lastdate' is given
				load_div.text = site+' has no data for '+date_val

	if var_check: # 'var_input' is not empty
		dum_text.value = str(time.time()) # click the timer button again to end the loading countdown

	if broken:
		load_div.text = broken_text
		

def read_nc(attr,old,new):
	'''
	callback for the 'var_input' dropdown widget.
	It reads and displays data based on the values of all the input widgets
	'''
	try: # attempt to get a variable
		site = site_input.value # the selected TCCON site
		prefix = [key for key in T_FULL if T_FULL[key]==site][0] # TCCON 2 letters abbreviation of the site
		site_file_list = [i for i in tccon_file_list if prefix in i] # a list of the associated netcdf files
	except: 
		pass # do nothing if an exception occurs (not sure if there can even be one here; maybe an empty 'site_input')
	else: # if no exception occured ...
		fig.yaxis.axis_label = var_input.value
		dum_text.value = str(time.time())	# click the timer button to start the loading countdown
		
		# loop over the TCCON files for the selected site
		data_check = False
		broken = False
		for filenum,site_file in enumerate(site_file_list):
			f = netCDF4.Dataset(os.path.join(tccon_path,site_file),'r') # netcdf file reader

			if filenum==0: # if it is the first file, configure the source and HoverTool based on the type of file
				if public:
					source.data.update({'x':[],'y':[],'colo':[]})
					fig.select_one(HoverTool).tooltips = [('value','@y')]
				else:
					fig.select_one(HoverTool).tooltips = [('value','@y'),('spectrum','@spectrum'),('flag','@flag')]

					source.data.update({'x':[],'y':[],'colo':[],'flag':[],'spectrum':[]})

			nctime =  f.variables['time'][:] # fractional days since 1970

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
					f.close() # close the netcdf reader
					break

			# for each file, fastforward to next iteration of the for loop if the dates in the file name are not compatible with the date input
			file_min = site_file[2:10]
			file_max = site_file[11:19]
			file_min_date = calendar.timegm(datetime(int(file_max[:4]),int(file_max[4:6]),int(file_max[6:8])).timetuple())/24/3600
			file_max_date = calendar.timegm(datetime(int(file_max[:4]),int(file_max[4:6]),int(file_max[6:8])).timetuple())/24/3600
			if (mindate>file_max_date) or (maxdate<file_min_date):
				f.close() # close the netcdf reader
				continue

			newtime = nctime[(nctime>=mindate) & (nctime<maxdate)]  # list of times that satisfy the 'date_input' value (still fractional days since 1970)

			# check that there is actual data in the time range, if not, go to the next (next iteration in for loop)
			if len(newtime)!=0:
				data_check = True
			else:
				continue

			start_id = np.where(nctime==newtime[0])[0][0] # ID of the starting time in the full time list of the file
			end_id = np.where(nctime==newtime[-1])[0][0]+1 # ID of the end time in the full time list of the file

			try: # attempt to fetch data using start_id and end_id
				add_x = np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in newtime])
				add_y = f.variables[var_input.value][start_id:end_id]
				if public:
					add_colo = np.array(['red']*len(add_x))
				else:
					add_colo = np.array(['red' if int(elem)==0 else 'grey' for elem in f.variables['flag'][start_id:end_id]]) # data with flag=0 is red; data with flag!=0 is grey
					add_spectrum = [''.join(elem) for elem in f.variables['spectrum'][start_id:end_id]]# name of spectra for the HoverTool
					add_flag = f.variables['flag'][start_id:end_id]# flags for the HoverTool
			except:
				pass # do nothing if any exception is raised (potentially index errors if the files have columns of inconsistent lengths)
			else: # if the try didnt raise any exceptions, update the 'source' of the plots with new data
				if public:
					source.data.update({	'x' : np.append(source.data['x'],add_x),
											'y' : np.append(source.data['y'],add_y),
											'colo' : np.append(source.data['colo'],add_colo)})
				else:
					source.data.update({	'x' : np.append(source.data['x'],add_x),
											'y' : np.append(source.data['y'],add_y),
											'colo' : np.append(source.data['colo'],add_colo),
											'flag' : np.append(source.data['flag'],add_flag),
											'spectrum':np.append(source.data['spectrum'],add_spectrum),})
			
			f.close() # close the netcdf reader		
		else: # else clause of the for loop, this will execute if no 'break' was encountered		
			if data_check is not True: # no data
				date_val = date_input.value # can be of the form 'firstdate'-'lastdate' or just 'firstdate'
				if date_val=='': # date_val empty mean the whole time range of the site has been evaluated
					load_div.text = site+' has no data' # if you see this one then there is a problem with the netcdf file ...
				elif len(date_val)==8: # if only 'firstdate' is given
					load_div.text = site+' has no data after '+date_val
				else: # if 'firstdate-lastdate' is given
					load_div.text = site+' has no data for '+date_val
			broken = False

		dum_text.value = str(time.time()) # click the timer button again to end the loading countdown

		if broken:
			load_div.text = broken_text

# assign the python callbacks to the input widgets
site_input.on_change('value',set_site)
var_input.on_change('value',read_nc)
date_input.on_change('value',read_nc)

def center():
	"""
	callback to update the y axis range based on the data.
	the range is by default auto adjusted to show all the data even huge outliers.
	This function takes the flag = 0 data (or all the data if it's from public files) and scales the y axis based on the mean of that data subset

	In short it will zoom in by disregarding outliers.
	"""

	main_data = source.data['y'][source.data['colo']=='red'].astype(np.float) # data in red is all the data from public files and flag=0 data from private files
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

			fig.y_range.start = min_y
			fig.y_range.end = max_y

		else: # if the variable is a constant, update the 'load_div' to let the user know why nothing happened
			load_div.text = 'Variable is constant'
	except IndexError:
		load_div.text = 'Variable is constant'


center_button = Button(label='Scale without extremas',name='test',width=220) # button to 'center' the plot on the 'good' data
center_button.on_click(center) # assign the callback function to the button

figrid = gridplot([[fig]] , toolbar_location = 'above') # put the figure by itself in a grid layout (I can better control where the toolbar will show if i do that)

if public:
	# callback of the dummy 'dum_hide' TextInput widget to hide all the dummy widgets, it will be triggered in a 'timeout' callback 1 second after the page is opened
	dum_hide.js_on_change('value',CustomJS(code="""document.getElementsByClassName("bk-widget-box")[1].style.display="none";""")) 	
	# callback of the dummy 'dum_text' TextInput widget to trigger a button click on the dummy button 'dum_button'
	dum_text.js_on_change('value',CustomJS(args={},code="document.getElementsByTagName('button')[0].click();"))
	# put the figure and all the widget in a single grid layout object
	grid = gridplot([[figrid,widgetbox(site_input,var_input,dumdiv,date_input,load_div,notes_div,width=600)],[widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=0,height=0)]], toolbar_location = None)

else: # for the 'private mode' only the index of the html divs containing the dummy widgets changes because of the addtion of the 'center_button' in the layout
	dum_hide.js_on_change('value',CustomJS(code="""document.getElementsByClassName("bk-widget-box")[2].style.display="none";"""))
	dum_text.js_on_change('value',CustomJS(args={},code="document.getElementsByTagName('button')[1].click();"))
	grid = gridplot([[figrid,widgetbox(site_input,var_input,dumdiv,date_input,load_div,notes_div,width=600)],[center_button],[widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=0,height=0)]], toolbar_location = None)

figrid.children[0].merge_tools = False # need to do that due to a bug in bokeh 0.12.6 that prevent gridplot to display the HoverTool in the toolbar

def hide_dummy():
	'''
	timeout callback that will be executed 1 second after the page is open.
	changes the value of the 'dum_hide' dummy TextInput widget, thise triggers the widget's javascript callback that makes all the dummy widgets invisible
	'''
	dum_hide.value = 'hide'

curdoc().title='TCCON' # this changes the title of the internet tab in which the document will be displayed
curdoc().add_root(grid) # this add the grid layout to the document
curdoc().add_timeout_callback(hide_dummy,1000) # this schedules the 'hide_dummy' callback to be triggered 1 second after page load
