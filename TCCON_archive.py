####################
# Code description #
####################
"""
This code will produce interactive plots  that will read into TCCON netcdf files to display time series of the data.

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
- 'comp' : there will be four plots, two figures for time series, and two 'correlation' figures that will plot the y axis of the time series figures against each other
	 	   In the time series figures, the left axis is for the first site variable while the right axis is for the second site variable

If the selected site's data is in different files it will take longer to load new variables.

The program will create a cache_dic.npy file in which it will save full time series of variables that correspond to previous inputs.
The size of this cache file will be kept under 200 MB. If you like you can increase that size by changing the 'cache_max_size' value (in bytes)

e.g.	Reading a new variable in 'comp' mode for Lamont for the entire time series can take several minutes when done for the first time.
		The same data read from the cache file will load in 20-50 seconds.
		For Lamont in 'comp' mode, a set of two variables (+time,flag,color, and spectrum arrays) is ~ 50 MB for the full time series
"""

#############
# Libraries #
#############

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
from bokeh.models import ColumnDataSource, TextInput, Div, CustomJS, Button, TextInput, Select, HoverTool, BoxSelectTool, DataTable, TableColumn, LinearAxis, DataRange1d, Tool
from bokeh.layouts import gridplot, widgetbox

# to ignore warnings
import warnings
warnings.filterwarnings('ignore')

#############
#############

layout_mode = sys.argv[1]
if layout_mode not in ['simple','comp']:
	print 'bokeh serve --show TCCON_archive.py --args A'
	print 'A is either "simple" or "comp"'
	sys.exit()

#########################################################################################################################################################################
## MODIFIABLE SECTION
tccon_path = os.path.join(os.getcwd(),'TCCON') ## this is the only line that may need editing; full path to the folder containing the tccon netcdf files

cache_max_size = 2E8 # maximum size of the cache file (in bytes), any new cached data after that will remove the oldest data following certain rules (see add_cache function)

# if you modify the plotting colors, you will need to remove your cache_dic.npy file
main_color = '#9ACD32' # this will be the color used for the flag=0 data; I use css 'YellowGreen' (#9ACD32) by default
flag_color = 'grey' # this will be the color used for the flag!=0 data
hover_color = 'red' # this will be the color used for data hovered by the mouse

main_color2 = '#DDA0DD' # the main color for data from the second site in 'comp' mode; I use css 'Plum' (#DDA0DD) by default

# this is the mode of the hovertool for the Figure 1 and Figure 2 in 'comp' mode
# 'width' to have all y data select over a x range (you can only pan the box selection left-right)
# 'height' can only pan the box selection up-down and it selects all the x data
# 'both' you can select a specific rectangle of data
boxselecttool_dimensions = 'both'

skip_list = ['_Version','ak_','prio','checksum','graw','spectrum','year','ada','aicf','adcf','time'] # variables including those keywords won't be shown in the 'var_input' dropdowns

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

# list of sites that will be available in the site_input selection widget; I add an empty string at the beginning so that the first site can be selected right away
ordered_site_list = ['']+sorted([T_FULL[i] for i in prefix_list]) 

#################################################################################
## FIGURES SETUP

site_input = Select(title='Site:',options = ordered_site_list,width=220) # dropdown widget to select the TCCON site
date_input = TextInput(title='start-end (yyyymmdd):',width=220) # text input widget to specify dates between which data should be fetched
flag_input = TextInput(title='Flag (an integer):',value='',width=130) # text input widget to specify the flag of the data to show
load_button = Button(label='Load Data',width=100) # this button will be used to start updating the plots according to the user inputs

TOOLS = "box_zoom,wheel_zoom,pan,box_select,redo,undo,hover" # the tools that will be available in the figure's toolbar

plot_width = 700
if layout_mode == 'simple':
	plot_height = 400
elif layout_mode == 'comp':
	plot_height = 200

fig = figure(output_backend="webgl",plot_width=plot_width,plot_height=plot_height,x_axis_type='datetime',tools=TOOLS,active_inspect=[],active_drag="box_zoom") # Figure 1
fig.min_border_left = 90
	
if layout_mode == 'comp':
	fig.min_border_right = 80

	source = ColumnDataSource(data={'x':[],'y1':[],'y2':[],'colo':[]}) # this is the object that stores the data to be plotted for the first site; it will be filled during callbacks
	source2 = ColumnDataSource(data={'x':[],'y1':[],'y2':[],'colo':[]}) # this is the object that stores the data to be plotted for the second site; it will be filled during callbacks

	# second figure to display time series of another variable
	fig2 = figure(output_backend="webgl",plot_width=plot_width,plot_height=plot_height,x_axis_type='datetime',x_range=fig.x_range,tools=TOOLS,active_inspect=[],active_drag="box_zoom") # Figure 2
	fig2.xaxis[0].axis_label = 'Time'
	fig2.scatter(x='x',y='y2',color='colo',hover_color=hover_color,alpha=0.7,source=source) # this is plotted in the second figure, it will read data from the 'source' object
	
	fig2.extra_y_ranges = {'second_var': DataRange1d()}
	new_fig_axis2= LinearAxis(y_range_name='second_var')
	fig2.add_layout(new_fig_axis2,'right')
	fig2.scatter(x='x',y='y2',color='colo',hover_color=hover_color,alpha=0.7,y_range_name='second_var',source=source2) # this is plotted in the second figure, it will read data from the 'source2' object	
	
	fig2.select_one(BoxSelectTool).dimensions = boxselecttool_dimensions

	# Figure 3, will plot y axis of Figure 1 vs y axis of Figure 2 using data from 'source' (the first site)
	fig3 = figure(output_backend="webgl",plot_width=350,plot_height=350,tools=TOOLS,active_inspect=[],active_drag="box_zoom")
	fig3.min_border_right = 80
	fig3.min_border_bottom = 170
	fig3.scatter(x='y2',y='y1',color='colo',hover_color=hover_color,alpha=0.7,source=source)

	# Figure 4, will plot y axis of Figure 1 vs y axis of Figure 2 using data from 'source2' (the second site)
	fig4 = figure(output_backend="webgl",plot_width=350,plot_height=350,tools=TOOLS,active_inspect=[],active_drag="box_zoom")
	fig4.min_border_left = 90
	fig4.scatter(x='y2',y='y1',color='colo',hover_color=hover_color,alpha=0.7,source=source2)

	# input widgets
	site_input2 = Select(title='Site:',options = ordered_site_list,width=220) # dropdown widget to select the second site

	var_input = Select(title='Figure 1 variable:',width=220) # dropdown widget to select the variable to plot in Figure 1
	var_input2 = Select(title='Figure 2 variable:',width=220) # dropdown widget to select the variable to plot in Figure 2
	var_input3 = Select(title='Figure 1 variable:',width=220) # dropdown widget to select the variable to plot in Figure 1
	var_input4 = Select(title='Figure 2 variable:',width=220) # dropdown widget to select the variable to plot in Figure 2
elif layout_mode == 'simple':
	source = ColumnDataSource(data={'x':[],'y1':[],'colo':[]}) # this is the object that stores the data to be plotted; it will be filled during callbacks
	fig.xaxis[0].axis_label = 'Time'
	var_input = Select(title='Variable to plot:',width=220) # dropdown widget to select the variable to plot in Figure 1

fig.scatter(x='x',y='y1',color='colo',hover_color=hover_color,alpha=0.7,source=source) # this is the plot of the first figure, it will read data from the 'source' object
if layout_mode == 'comp':
	fig.extra_y_ranges = {'first_var': DataRange1d()}
	new_fig_axis = LinearAxis(y_range_name='first_var')
	fig.add_layout(new_fig_axis,'right')
	fig.scatter(x='x',y='y1',color='colo',hover_color=hover_color,alpha=0.7,y_range_name='first_var',source=source2) # this is plotted in the first figure, it will read data from the 'source2' object	

fig.select_one(BoxSelectTool).dimensions = boxselecttool_dimensions

## END FIGURES SETUP
#################################################################################
## SIDE WIDGETS SETUP

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

notes_div = Div(text=notes,width=400) # this will display the 'notes' above
status_div = Div(text='Select a site',width=400) # the text of this widget will be updated with information on the state of callbacks
select_text = Div(text='',width = 400) # text div that will be updated with the selected range of date within the BoxSelect tool
dumdiv = Div(text='',height=10) # dummy empty Div widget for spacing

## END OF SIDE WIDGETS SETUP
#################################################################################
## TOOLS SETUP

# this code will trigger when the source data is changed in order to update the data table that shows correlations of selected variables
correlation_code = """
var inds = cb_obj.selected['1d'].indices;
var data = cb_obj.data;
var tab = dt.source.data;

tab['Site'][0] = site_input.value;
tab['Site'][1] = site_input2.value;

var ym1 = 0;
var ym2 = 0;

var T1 = 0;
var T2 = 0;
var T3 = 0;

tab['N'][ROWID] = inds.length;

if (inds.length == 0) {
	tab['R'][ROWID] = 0;
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

tab['R'][ROWID] = (T1/Math.sqrt(T2*T3)).toFixed(3);

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
# the boxselect code of the correlation figure is different since it will not have time on any of its axes !, it just clicks the rescale button after selection
corfig_box_select_code = """
setTimeout(function(){
	button_list = document.getElementsByTagName('button');

	for(i=0;i<button_list.length;i++){
		if(button_list[i].textContent.includes("Scale")){button_list[i].click()}
	}
}, 1000); // I need the timeout because it sometimes misses without it"""

# in 'comp' mode, each of the 3 plots has a different hover tool, this button can click them all at once for convenience.
hover_button_code="""
if(cb_obj.button_type.includes("success")){
cb_obj.button_type = 'danger';
cb_obj.label = 'Disable hover tools'
} else {
	
cb_obj.button_type = 'success';
cb_obj.label= 'Enable hover tools';
}

var toolbar_list = document.getElementsByClassName("bk-toolbar-box");
if(toolbar_list.length>1){
	for(i=0;i<toolbar_list.length;i++){
		if(toolbar_list[i].textContent.includes("Hover")){}
	}
}

var toolbar_button_list = document.getElementsByClassName("bk-toolbar-button");
if(toolbar_button_list.length>1){
	for(i=0;i<toolbar_button_list.length;i++){
		if(toolbar_button_list[i].textContent.includes("Hover")){toolbar_button_list[i].click()}
	}
}
"""
hover_button = Button(label='Enable hover',button_type='success',width=140)
hover_button.callback = CustomJS(code=hover_button_code)

#
reset_tool_code = """
import * as p from "core/properties"
import {ResetTool,ResetToolView} from "models/tools/actions/reset_tool"

export class NewResetToolView extends ResetToolView
	
	NewResetToolView::doit = ->
	  @plot_view.clear_state()
	  @plot_view.reset_range()
	  @plot_view.reset_selection()
	  custom_event()
	  if @model.reset_size
	    return @plot_view.reset_dimensions()
	  return

export class NewResetTool extends ResetTool

	default_view : NewResetToolView
	type : "NewResetTool"
	tool_name : "NewReset"
	icon : "bk-tool-icon-reset"
"""

class NewResetTool(Tool):
	__implementation__ = reset_tool_code

## END OF TOOLS SETUP
#################################################################################
## DUMMY WIDGETS SETUP

# javascript code for a dummy (invisible) button, it starts and stops a timer that will be displayed in the 'status_div' when data is loading
dum_button_code = """
if (cb_obj.button_type.includes('success')){
var start = new Date();	
var intervalID = setInterval(function(){var current = new Date(); var diff=((current-start)/1000.0).toFixed(1); status_div.text='Loading data: '+diff.toString()+' s';	}, 100)
cb_obj.button_type = 'danger';
} else {
var noIntervals = setInterval(function(){});
for (var i = 0; i<noIntervals; i++) { window.clearInterval(i);}
status_div.text='Data loaded';
cb_obj.button_type = 'success';
}
"""
dum_button = Button(label='dummy',button_type='success',width=200) # the dummy button itself
dum_button.callback = CustomJS(args={'status_div':status_div},code=dum_button_code) # the callback of the button

dum_text = TextInput() # dummy text input widget; it will be used to trigger the dummy button callback when the dropdown widgets are used.
# callback of the dummy 'dum_text' TextInput widget to trigger a button click on the dummy button 'dum_button'
dum_text_code = """
button_list = document.getElementsByTagName('button');

for(i=0;i<button_list.length;i++){
	if(button_list[i].textContent.includes("dummy")){button_list[i].click()}
}
"""
dum_text.js_on_change('value',CustomJS(code=dum_text_code))

# callback to hide the dummy widget box after the timeout callback
# also changes the css of the 'site_input' widgets to display them in bold and with the apprioriate color
# this callback will also be triggered at the start of the set_site function (otherwise the css goes back to the default css when an option is selected)
dum_hide_code = """
console.log('DUM_HIDE_CODE');
if (document.getElementById('mystyle')===null){
	var widgetbox_list = document.getElementsByClassName("bk-widget-box");

	for(i=0;i<widgetbox_list.length;i++){
		if(widgetbox_list[i].textContent.includes('hidden')){widgetbox_list[i].style.display="none"}
	}

	var toolbar_list = document.getElementsByClassName("bk-toolbar-box");
	if(toolbar_list.length>1){
		for(i=0;i<toolbar_list.length;i++){
			if(toolbar_list[i].textContent.includes("Hover")){toolbar_list[i].style.display="none"}
		}
	}

	var custom_css = document.createElement("style");
	custom_css.id = 'mystyle';
	custom_css.type = "text/css";
	custom_css.innerHTML = ".bk-tooltip>div:not(:first-child) { display: none }";
	document.body.appendChild(custom_css);

	var custom_script = document.createElement("script");
	custom_script.id = 'myscript';
	custom_script.type = "text/javascript";
	custom_script.innerHTML = `	function custom_event(){
		// make the background color of the site_input widgets correspond to the color of data in the plot
		setTimeout( function(){
			var first_site = document.getElementsByTagName('label')[0];
			var second_site = document.getElementsByTagName('label')[1];
			first_site.style["font-weight"] = "bold";
			second_site.style["font-weight"] = "bold";
			first_site.style["color"] = "%s";
			second_site.style["color"] = "%s";
		}, 10)

		// make the TextInput widgets smaller
		setTimeout(function(){
			var textinputs = document.getElementsByClassName('bk-widget-form-group bk-layout-fixed');
			for(i=0;i<textinputs.length;i++){
				if(textinputs[i].textContent==='start-end (yyyymmdd):'){textinputs[i].children[1].style.width='160px'}
				if(textinputs[i].textContent==='Flag (an integer):'){textinputs[i].children[1].style.width='20px'}
			}
		}, 10)
	}`;
	document.body.appendChild(custom_script);
}

custom_event();

""" % (main_color,main_color2)

dum_hide = TextInput() # dummy text input widget; it will be used in a 'timeout' callback after the page load in order to make the dummy widgets invisible
dum_hide.js_on_change('value',CustomJS(code=dum_hide_code))

# if the 'timeout' callback is set too early, the dummy widgets won't be hidden.
# but if it is set too late, the user will have time to see the dummy widgets before they disappear.
dum_alert = Div(text='<b>OOPS !</b> the widgets below are supposed to be hidden. Try refreshing the page',width=700)

dum_box = widgetbox(dum_alert,dum_button,dum_text,dum_hide,width=600,name='dummy box')

def hide_dummy():
	'''
	timeout callback that will be executed 1 second after the page is open.
	changes the value of the 'dum_hide' dummy TextInput widget, this triggers the widget's javascript callback that makes all the 'dum_box' invisible
	'''
	dum_hide.value = 'hide'
## END OF DUMMY WIDGETS SETUP
#################################################################################
## CACHE SETUP
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
## END OF CACHE SETUP
#################################################################################

## END OF SETUP SECTION
#########################################################################################################################################################################
## ADD_DATA FUNCTION
def add_data(data,x=None,y1=None,y2=None,colo=None,flag=None,spectrum=None):
	'''
	function to update a data dictionary based on the layout mode and data type
	'''
	if public and layout_mode=='simple':
		return {'x' : np.append(data['x'],x),
				'y1' : np.append(data['y1'],y1),
				'colo' : np.append(data['colo'],colo),}
	elif not public and layout_mode=='simple':
		return {'x' : np.append(data['x'],x),
				'y1' : np.append(data['y1'],y1),
				'colo' : np.append(data['colo'],colo),
				'flag' : np.append(data['flag'],flag),
				'spectrum': np.append(data['spectrum'],spectrum),}
	elif public and layout_mode=='comp':
		return {'x' : np.append(data['x'],x),
				'y1' : np.append(data['y1'],y1),
				'y2' : np.append(data['y2'],y2),
				'colo' : np.append(data['colo'],colo),}
	elif not public and layout_mode=='comp':
		return {'x' : np.append(data['x'],x),
				'y1' : np.append(data['y1'],y1),
				'y2' : np.append(data['y2'],y2),
				'colo' : np.append(data['colo'],colo),
				'flag' : np.append(data['flag'],flag),
				'spectrum': np.append(data['spectrum'],spectrum),}
## END OF ADD_DATA FUNCTION
#########################################################################################################################################################################
## ADD_CACHE FUNCTION
def add_cache(date_val,site,source_data,max_size=cache_max_size,first_var='',second_var=''):
	'''
	This function will add data to the cache_dic dictionnary in order to make load time shorter.
	It will keep the size of the cache file below 'max_size' (in bytes)
	Whenever that size is exceeded, it will remove the oldest cached data following certain priority rules.
	'''

	global cache_dic

	try:
		cache_dic = np.load(cache_file).item()
	except IOError:
		print cache_file+' not found'
		print 'Creating new empty '+cache_file
		cache_dic = {}
		np.save(cache_file,cache_dic)

	stamp = datetime.now()

	keep_loop = True
	rewritten = False
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
				print 'Caching new time data:',date_val,site,len(source_data['x']),'values'
				cache_dic[date_val][site]['x'] = source_data['x']

			for var in source_data:
				if var not in ['x','colo','']:
					nice_var = var
					if var=='y1':
						nice_var = first_var
					elif var=='y2':
						nice_var = second_var

					try:
						type(cache_dic[date_val][site][nice_var])
					except KeyError:
						print 'Caching new variable:',date_val,site,nice_var,len(source_data[var]),'values'
						cache_dic[date_val][site][nice_var] = {}
						cache_dic[date_val][site][nice_var]['stamp'] = stamp
						cache_dic[date_val][site][nice_var]['value'] = source_data[var]
		else:
			rewritten = True
			while os.path.getsize(cache_file)>max_size:
				site_list = list(set([[site for site in cache_dic[date_val] if site!='stamp'] for date_val in cache_dic][0])) # list(set()) removes duplicates
				var_list = list(set([[[var for var in cache_dic[date_val][site] if var not in ['x','flag','spectrum','colo','stamp']] for site in cache_dic[date_val] if site!='stamp'] for date_val in cache_dic][0][0]))

				# if there are more than 5 date_val, remove the oldest date_val
				if len(cache_dic)>3:
					del_date = min([(cache_dic[date_val]['stamp'],date_val) for date_val in cache_dic if date_val!=''])[1]
					del cache_dic[del_date]
					print 'Removing',del_date,'from cache'
				# if there are more than 5 sites, remove the oldest site
				elif len(site_list)>2:
					del_site = min([[(cache_dic[date_val][site]['stamp'],(date_val,site)) for site in cache_dic[date_val] if site!='stamp'] for date_val in cache_dic][0])[1]
					del cache_dic[del_site[0]][del_site[1]]
					print 'Removing',del_site[0],del_site[1],'from cache'
				# if there are more than 20 variables, remove the oldest variable (except for times, flags, colors, and spectrum names; those will only be removed when a site or date_val is removed)
				elif len(var_list)>10:
					del_var = min([[[(cache_dic[date_val][site][var]['stamp'],(date_val,site,cache_dic[date_val][site][var]['value'])) for var in cache_dic[date_val][site] if var not in ['x','flag','spectrum','colo','stamp']] for site in cache_dic[date_val] if site !='stamp'] for date_val in cache_dic][0][0])
					del cache_dic[del_var[0]][del_var[1]][del_var[2]]
					print 'Removing',del_var[0],del_var[1],del_var[2],'from cache'
				print 'Saving new cache file ...'
				np.save(cache_file,cache_dic)
				print cache_file,'size:',os.path.getsize(cache_file)/1E6,'MB'

		if os.path.getsize(cache_file)<max_size:
			keep_loop = False
		else:
			print 'Cache is over',max_size/1E6,'MB after caching iteration',it
			print 'If this does not stop iterating it means that the current selection alone is over the max_size'
			it+=1
	if not rewritten:
		print 'Saving new cache file ...'
		np.save(cache_file,cache_dic)
		print cache_file,'size:',os.path.getsize(cache_file)/1E6,'MB'
## END OF ADD_CACHE FUNCTION
#########################################################################################################################################################################
## INITIALIZE FUNCTION
def initialize(var_list,site_source,site_ID,reset=True):
	'''
	Function that fills the var_input widgets with options and resets the data source
	Also sets up the hovertool tooltips according to the file type and mode
	'''
	if site_ID == 1:
		site_var_input = var_input
		if layout_mode == 'comp':
			corfig = fig3
			site_var_input2 = var_input2
	elif site_ID == 2:
		corfig = fig4
		site_var_input = var_input3
		site_var_input2 = var_input4 

	if len(site_var_input.options)==0:
		var_list = ['']+sorted([var for var in var_list if True not in [elem in var for elem in skip_list]])
		site_var_input.options = var_list # fill the 'var_input' dropdown with options
		site_var_input.value = ''
		if layout_mode == 'comp':
			site_var_input2.options = var_list
			site_var_input2.value = ''
	elif var_list==[]:
		site_var_input.options = var_list
		site_var_input.value = ''
		if layout_mode == 'comp':
			site_var_input2.options = var_list
			site_var_input2.value = ''

	# configure the source and HoverTool based on the type of file
	if reset:
		if public and layout_mode=='simple':
			site_source.data.update({'x':[],'y1':[],'colo':[]})
			fig.select_one(HoverTool).tooltips = [('y','@y1'),('date','@x{%F %T}')]
			fig.select_one(HoverTool).formatters = {'x':'datetime'}
		elif not public and layout_mode=='simple':
			site_source.data.update({'x':[],'y1':[],'colo':[],'flag':[],'spectrum':[]})
			fig.select_one(HoverTool).tooltips = [('value','@y1'),('spectrum','@spectrum'),('flag','@flag'),('date','@x{%F %T}')]
			fig.select_one(HoverTool).formatters = {'x':'datetime'}
		elif public and layout_mode=='comp':
			site_source.data.update({'x':[],'y1':[],'y2':[],'colo':[]})
			fig.select_one(HoverTool).tooltips = [('y','@y1'),('Fig2 y','@y2'),('date','@x{%F %T}')]
			fig2.select_one(HoverTool).tooltips = [('y','@y2'),('Fig1 y','@y1'),('date','@x{%F %T}')]
			corfig.select_one(HoverTool).tooltips = [('y','@y1'),('x','@y2'),('date','@x{%F %T}')]
		elif not public and layout_mode=='comp':
			site_source.data.update({'x':[],'y1':[],'y2':[],'colo':[],'flag':[],'spectrum':[]})
			fig.select_one(HoverTool).tooltips = [('y','@y1'),('spectrum','@spectrum'),('flag','@flag'),('Fig2 y','@y2'),('date','@x{%F %T}')]
			fig2.select_one(HoverTool).tooltips = [('y','@y2'),('spectrum','@spectrum'),('flag','@flag'),('Fig1 y','@y1'),('date','@x{%F %T}')]
			corfig.select_one(HoverTool).tooltips = [('y','@y1'),('x','@y2'),('spectrum','@spectrum'),('flag','@flag'),('date','@x{%F %T}')]

		if layout_mode == 'comp':
			for curfig in [fig,fig2,corfig]:
				curfig.select_one(HoverTool).formatters = {'x':'datetime'}

	dum_hide.value = str(time.time()) # update the css of labels
## END OF INITIALIZE FUNCTION
#########################################################################################################################################################################
## SET_SITE FUNCTION
def set_site(site='',site_ID=0):
	'''
	callback for the 'site_input' dropdown widget.
	This will only fill the 'var_input' dropdown widgets with options.
	'''

	dum_hide.value = str(time.time()) # update the css of site_input label

	if site_ID==1:
		site_source = source
	elif site_ID==2:
		site_source = source2

	if site == '':
		initialize([],site_source,site_ID)
		if layout_mode == 'simple':
			status_div.text = 'Select a site'
		elif layout_mode =='comp':
			if site_input.value==site_input2.value=='':
				status_div.text = 'Select a site'
			elif site_input.value==''!=site_input2.value:
				if '' in [var_input3.value,var_input4.value]:
					status_div.text = site_input2.value+' still has an empty variable input'
				else:
					status_div.text = 'Data ready to load'
			elif site_input2.value==''!=site_input.value:
				if '' in [var_input.value,var_input2.value]:
					status_div.text = site_input.value+' still has an empty variable input'
				else:
					status_div.text = 'Data ready to load'
		return
	
	prefix = [key for key in T_FULL if T_FULL[key]==site][0] # TCCON 2 letters abbreviation of the site
	site_file_list = [i for i in tccon_file_list if prefix in i] # a list of the associated files

	if layout_mode=='comp' and site_ID==1:
		corfig = fig3
	elif site_ID==2:
		corfig = fig4
	
	# update titles, labels, and notes
	fig.yaxis[0].axis_label = var_input.value
	if layout_mode == 'comp':
		fig.yaxis[1].axis_label = var_input3.value
		fig2.yaxis[0].axis_label = var_input2.value
		fig2.yaxis[1].axis_label = var_input4.value

		fig3.yaxis[0].axis_label = var_input.value
		fig3.xaxis[0].axis_label = var_input2.value

		fig3.yaxis[0].axis_label = var_input.value
		fig3.xaxis[0].axis_label = var_input2.value

		fig4.yaxis[0].axis_label = var_input3.value
		fig4.xaxis[0].axis_label = var_input4.value

		corfig.title.text =  T_FULL[prefix]+', '+T_LOC[prefix] # site name + site location

		if site_input.value==site_input2.value and ('' not in [site_input.value,site_input2.value]):
			notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix]+"'>"+site+"</a>") # updated the information widget with a link to the site's webpage
		elif '' not in [site_input.value,site_input2.value]:
			prefix1 = [key for key in T_FULL if T_FULL[key]==site_input.value][0]
			prefix2 = [key for key in T_FULL if T_FULL[key]==site_input2.value][0]
			notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix1]+"'>"+site_input.value+"</a>, <a href='"+T_site[prefix2]+"'>"+site_input2.value+"</a>") # updated the information widget with a link to the site's webpage
		elif site_input2.value=='' and site_input.value=='':
			notes_div.text = notes
		elif site_input.value=='':
			prefix2 = [key for key in T_FULL if T_FULL[key]==site_input2.value][0]
			notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix2]+"'>"+site_input2.value+"</a>") # updated the information widget with a link to the site's webpage
		elif site_input2.value=='':
			prefix1 = [key for key in T_FULL if T_FULL[key]==site_input.value][0]
			notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix1]+"'>"+site_input.value+"</a>") # updated the information widget with a link to the site's webpage

	elif layout_mode == 'simple':
		fig.title.text =  T_FULL[prefix]+', '+T_LOC[prefix] # site name + site location

		notes_div.text = notes.replace("a></font>","a></font>, <a href='"+T_site[prefix]+"'>"+site+"</a>") # updated the information widget with a link to the site's webpage

	load_var(site_file_list,site,site_source,site_ID,mode="set_site")
## END OF SET_SITE FUNCTION
#########################################################################################################################################################################
## SET_VAR FUNCTION
def set_var(site='',site_ID=0):
	'''
	This function reads and displays data based on the values of the input widgets of a given site
	'''

	if site_ID==1:
		site_source = source
	elif site_ID==2:
		site_source = source2

	if site == '':
		initialize([],site_source,site_ID)
		return

	dum_hide.value = str(time.time()) # update the css of site_input label

	prefix = [key for key in T_FULL if T_FULL[key]==site][0] # TCCON 2 letters abbreviation of the site
	site_file_list = [i for i in tccon_file_list if prefix in i] # a list of the associated files

	fig.yaxis[0].axis_label = var_input.value
	if layout_mode == 'comp':
			fig.yaxis[1].axis_label = var_input3.value
			fig2.yaxis[0].axis_label = var_input2.value
			fig2.yaxis[1].axis_label = var_input4.value

			fig3.yaxis[0].axis_label = var_input.value
			fig3.xaxis[0].axis_label = var_input2.value

			fig3.yaxis[0].axis_label = var_input.value
			fig3.xaxis[0].axis_label = var_input2.value

			fig4.yaxis[0].axis_label = var_input3.value
			fig4.xaxis[0].axis_label = var_input4.value
	
	load_var(site_file_list,site,site_source,site_ID,mode="set_var")
## END OF SET_VAR FUNCTION
#########################################################################################################################################################################
## LOAD_DATA FUNCTION
if layout_mode == 'simple':
	save_inputs = ['']*5
elif layout_mode =='comp':
	save_inputs = ['']*8

def load_data():
	'''
	callback for the load_button, it will load data based on the values of the input widgets for both sites
	'''
	global save_inputs

	if layout_mode=='simple' and site_input.value=='':
		status_div.text = 'Select a site'
		return	
	elif layout_mode=='comp':
		if site_input.value==site_input2.value=='':
			status_div.text = 'Select a site'
			return

	if layout_mode == 'simple':
		current_inputs = [site_input.value,var_input.value,var_input2.value,date_input.value,flag_input.value]
	elif layout_mode == 'comp':
		current_inputs = [site_input.value,site_input2.value,var_input.value,var_input2.value,var_input3.value,var_input4.value,date_input.value,flag_input.value]

	no_new_inputs = save_inputs==current_inputs
	save_inputs = current_inputs
	
	if '' not in [site_input.value,var_input.value,var_input2.value]:
		if no_new_inputs:
			status_div.text = 'Data already loaded'
		else:
			set_var(site=site_input.value,site_ID=1)
	elif site_input.value!='' and ('' in [var_input.value,var_input2.value]):
		status_div.text = site_input.value+' still has an empty variable input'
		return
		
	if layout_mode == 'comp':
		if '' not in [site_input2.value,var_input3.value,var_input4.value]:
			if no_new_inputs:
				status_div.text = 'Data already loaded'
			else:
				set_var(site=site_input2.value,site_ID=2)
		elif site_input2.value!='' and ('' in [var_input3.value,var_input4.value]):
			status_div.text = site_input2.value+' still has an empty variable input'
## END OF LOAD_DATA FUNCTION
#########################################################################################################################################################################
## LOAD_VAR FUNCTION
def load_var(site_file_list,site,site_source,site_ID,mode=""):
	'''
	Function called by set_site() and set_var() to load the variables matching site,variable, and date inputs
	'''

	global cache_dic, all_var

	if site_ID==1:
		no_flag_color = main_color
		first_var = var_input.value
		if layout_mode == 'comp':
			second_var = var_input2.value
	elif site_ID==2:
		no_flag_color = main_color2
		first_var = var_input3.value
		second_var = var_input4.value

	date_val = date_input.value # can be of the form 'firstdate'-'lastdate' or just 'firstdate'

	dum_text.value = str(time.time()) # click the timer button to start the loading countdown in the 'status_div' widget

	filled_site_inputs = False
	if layout_mode == 'simple':
		if first_var!='':
			filled_site_inputs = True
	elif layout_mode == 'comp':
		if '' not in [first_var,second_var]:
			filled_site_inputs = True

	# check if there already is cached data that matches the inputs
	no_cached_data = False
	try:
		type(cache_dic[date_val][site][first_var])
	 	if layout_mode == 'comp':
	 		type(cache_dic[date_val][site][second_var])
	except KeyError:
		no_cached_data = True

	if not public:
		flag_mode = flag_input.value
		if flag_mode != '':
			if not float(flag_mode).is_integer():
				dum_text.value = str(time.time()+2) # click the timer button to start the loading countdown in the 'status_div' widget
				time.sleep(0.1)
				status_div.text = "The flag should be an integer"
				return

	if no_cached_data or mode=='set_site':
		# loop over the TCCON files for the selected site
		for filenum,site_file in enumerate(site_file_list):
			if netcdf:
				f = netCDF4.Dataset(os.path.join(tccon_path,site_file),'r') # netcdf file reader
				all_var = [var for var in f.variables if 'run' not in var]
			else:
				df = pd.read_csv(os.path.join(tccon_path,site_file),header=2) # read the .eof.csv file
				all_var = [var for var in list(df) if 'run' not in var]

			# setup some initializations if it is the first file
			if filenum==0:		
				if not filled_site_inputs:
					initialize(all_var,site_source,site_ID,reset=False) # fills variable inputs with options, does not reset the data source
					dum_text.value = str(time.time()+3) # click the timer button again to stop the loading countdown in the 'status_div' widget
					time.sleep(0.1)
					status_div.text = site+' still has an empty variable input'
					return
				if mode == 'set_site':
					initialize(all_var,site_source,site_ID,reset=False) # fills variable inputs with options, does not reset the data source
					dum_text.value = str(time.time()+7) # click the timer button again to stop the loading countdown in the 'status_div' widget
					time.sleep(0.1)
					status_div.text = 'Data ready to load'
					return
				
				initialize(all_var,site_source,site_ID) # fills variable inputs with options,resets the data source, and setup hovertool tooltips
				save_data = {key:[] for key in site_source.data}
				new_data = {key:[] for key in site_source.data}					
	
			if netcdf:
				nctime =  f.variables['time'][:] # fractional days since 1970
			else:			
				# for .eof.csv files convert year,day,hour in fractional days since 1970
				nctime = np.array([(datetime(int(df.year[i]),1,1)+timedelta(days=df.day[i]-1)-datetime(1970,1,1)).days+df.hour[i]/24 for i in range(len(df))])

			# use the value of the 'date_input' widget to determine the range of dates over which data should be fetched
			try: # for the minimum of the range
				mindate = date_val.split('-')[0] # if the 'date_input' widget is empty, this will raise an exception
				mindate = calendar.timegm(datetime(int(mindate[:4]),int(mindate[4:6]),int(mindate[6:8])).timetuple())/24/3600 # if the date is entered wrong, this will raise an exception
			except: # catch all exceptions
				mindate = calendar.timegm(datetime(1970,1,1).timetuple())/24/3600 # if an exception has been caught, just use a very early date

			try: # for the maximum of the range 
				maxdate = date_val.split('-')[1] # if the 'date_input' widget is empty, or if only the minimum date is given, this will raise an exception
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
					dum_text.value = str(time.time()+3) # click the timer button to start the loading countdown in the 'status_div' widget
					time.sleep(0.1)
					status_div.text = site +' date range: '+all_file_min+'-'+all_file_max
					if netcdf:
						f.close() # close the netcdf reader
					return
				if mindate>maxdate:
					dum_text.value = str(time.time()+3) # click the timer button to start the loading countdown in the 'status_div' widget
					time.sleep(0.1)
					status_div.text = 'Wrong date input'
					if netcdf:
						f.close() # close the netcdf reader
					return

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
					add_y1 = f.variables[first_var][start_id:end_id]
					add_y2 = ''
					if layout_mode == 'comp':
						add_y2 = f.variables[second_var][start_id:end_id]

				else:
					add_y1 = np.array(df[first_var][start_id:end_id])
					add_y2 = ''
					if layout_mode == 'comp':
						add_y2 = np.array(df[second_var][start_id:end_id])

				if not public:
					if netcdf:
						add_spectrum = np.array([''.join(elem) for elem in f.variables['spectrum'][start_id:end_id]])# name of spectra for the HoverTool
						add_flag = f.variables['flag'][start_id:end_id].astype(int) # flags for the HoverTool
					else:
						add_spectrum = np.array([elem for elem in df['spectrum'][start_id:end_id]])
						add_flag = np.array([' '.join([str(flag),all_var[flag],'=',str(df[all_var[flag]][start_id:end_id][ID])]) for ID,flag in enumerate(df['flag'][start_id:end_id])])
				else:
					add_spectrum = ''
					add_flag = ''
			except KeyError:
				print 'KeyError'
			else: # if the try didnt raise any exceptions, update the 'source' of the plots with new data
				if not public:
					if netcdf:
						add_colo = np.array([no_flag_color if int(elem)==0 else flag_color for elem in add_flag])
					else:
						add_colo = np.array([no_flag_color if int(elem.split()[0])==0 else flag_color for elem in add_flag])

					save_data = add_data(save_data,x=add_x,y1=add_y1,y2=add_y2,colo=add_colo,flag=add_flag,spectrum=add_spectrum)

					if flag_mode != '':
						if netcdf:
							inds = add_flag==int(flag_mode)
						else:
							inds = [flag_mode==elem.split()[0] for elem in add_flag]
						add_x = add_x[inds]
						add_y1 = add_y1[inds]
						add_colo = add_colo[inds]
						add_spectrum = add_spectrum[inds]
						add_flag = add_flag[inds]
						if layout_mode == 'comp':
							add_y2 = add_y2[inds]

						if True not in inds:
							continue

				elif public:
					add_colo = np.array([no_flag_color]*len(add_x))
					save_data = add_data(save_data,x=add_x,y1=add_y1,y2=add_y2,colo=add_colo,flag=add_flag,spectrum=add_spectrum)

				new_data = add_data(new_data,x=add_x,y1=add_y1,y2=add_y2,colo=add_colo,flag=add_flag,spectrum=add_spectrum)					

				if filenum == 0:
					print 'load_var() ...'
				print site_file,'data read:',len(add_x),'new values'
			if netcdf:
				f.close() # close the netcdf reader		
		else: # else clause of the for loop, this will execute if no 'break' was encountered
			print 'Updating',site,'data source ...'
			site_source.data.update(new_data)
			del new_data			
			if len(site_source.data['x'])==0: # no data
				add_flag_message = ''
				if flag_mode != '':
					add_flag_message = ' with flag='+flag_mode
				dum_text.value = str(time.time()+8) # click the timer button again to end the loading countdown
				time.sleep(0.1)
				if date_val=='': # date_val empty mean the whole time range of the site has been evaluated
					status_div.text = site+' has no data'+add_flag_message # if you see this one then there is a problem with the netcdf file ...
				elif len(date_val)==8: # if only 'firstdate' is given
					status_div.text = site+' has no data'+add_flag_message+' after '+date_val
				elif len(date_val)>8: # if 'firstdate-lastdate' is given
					status_div.text = site+' has no data'+add_flag_message+' for '+date_val
				return
			print site,'data source updated:',len(site_source.data['x']),'values'
		if len(site_source.data['y1'])!=0:
			if layout_mode == 'simple':
				add_cache(date_val,site,save_data,first_var=first_var)
			elif layout_mode == 'comp':
				add_cache(date_val,site,save_data,first_var=first_var,second_var=second_var)

	else: #else clause of 'if no_cached_data', this will execute if there is cached data corresponding to the inputs
		initialize(all_var,site_source,site_ID) # fills variable inputs with options,resets the data site_source, and setup hovertool tooltips
		if not filled_site_inputs:
			dum_text.value = str(time.time()+1) # click the timer button to start the loading countdown in the 'status_div' widget
			time.sleep(0.1)
			status_div.text = site+' still has an empty variable input'
			return

		if not public:
			if netcdf:
				new_colo = np.array([no_flag_color if int(elem)==0 else flag_color for elem in cache_dic[date_val][site]['flag']['value']])
			else:
				new_colo = np.array([no_flag_color if int(elem.split()[0])==0 else flag_color for elem in cache_dic[date_val][site]['flag']['value']])
		
			if flag_mode != '':
				if netcdf:
					inds = cache_dic[date_val][site]['flag']['value']==int(flag_mode)
				else:
					inds = [flag_mode==elem.split()[0] for elem in cache_dic[date_val][site]['flag']['value']]

				if True not in inds:
					dum_text.value = str(time.time()+5) # click the timer button again to end the loading countdown
					time.sleep(0.1)
					add_flag_message = ' with flag='+flag_mode
					if date_val=='': # date_val empty mean the whole time range of the site has been evaluated
						status_div.text = site+' has no data'+add_flag_message # if you see this one then there is a problem with the netcdf file ...
					elif len(date_val)==8: # if only 'firstdate' is given
						status_div.text = site+' has no data'+add_flag_message+' after '+date_val
					elif len(date_val)>8: # if 'firstdate-lastdate' is given
						status_div.text = site+' has no data'+add_flag_message+' for '+date_val
					return
			else:
				inds = [True for elem in cache_dic[date_val][site]['flag']['value']]

			add_x = cache_dic[date_val][site]['x'][inds]
			add_y1 = cache_dic[date_val][site][first_var]['value'][inds]
			add_colo = new_colo[inds]
			add_flag = cache_dic[date_val][site]['flag']['value'][inds]
			add_spectrum = cache_dic[date_val][site]['spectrum']['value'][inds]
			add_y2 = ''
			if layout_mode == 'comp':
				add_y2 = cache_dic[date_val][site][second_var]['value'][inds]
		
		elif public:
			add_x = cache_dic[date_val][site]['x']
			add_y1 = cache_dic[date_val][site][first_var]['value']
			add_colo = np.array([no_flag_color]*len(cache_dic[date_val][site]['x']))
			add_flag = ''
			add_spectrum = ''
			add_y2 = ''
			if layout_mode == 'comp':
				add_y2 = cache_dic[date_val][site][second_var]['value']

		print 'load_var() ...'
		print 'Using cached data for',date_val,site
		site_source.data.update(add_data(site_source.data,x=add_x,y1=add_y1,y2=add_y2,colo=add_colo,flag=add_flag,spectrum=add_spectrum))

	if len(site_source.data['x'])!=0:
		dum_text.value = str(time.time()+6) # click the timer button again to end the loading countdown
		print 'If all the data is not showing quickly, you should use the date_input widget to select a smaller subset of data'
	print 'load_var() DONE'
## END OF LOAD_VAR FUNCTION
#########################################################################################################################################################################
## INPUT WIDGETS CALLBACKS SECTION
# assign the python callbacks to the input widgets
site_input.on_change('value',lambda attr,old,new: set_site(site=site_input.value,site_ID=1))
var_input.on_change('value',lambda attr,old,new: set_site(site=site_input.value,site_ID=1))
var_input2.on_change('value',lambda attr,old,new: set_site(site=site_input.value,site_ID=1))
date_input.js_on_change('value',CustomJS(code=dum_hide_code)) # conserve the size of the input widget after using it

load_button.on_click(load_data)

if layout_mode == 'comp':
	site_input2.on_change('value',lambda attr,old,new: set_site(site=site_input2.value,site_ID=2))
	var_input3.on_change('value',lambda attr,old,new: set_site(site=site_input2.value,site_ID=2))
	var_input4.on_change('value',lambda attr,old,new: set_site(site=site_input2.value,site_ID=2))

	flag_input.js_on_change('value',CustomJS(code=dum_hide_code)) # conserve the size of the input widget after using it
	
	# widgets specific to the comparison layout_mode
	table_source = ColumnDataSource( data = {'Site':['',''],'N':[0,0],'R':[0,0]} ) # the data source of the table
	data_table = DataTable(source=table_source, reorderable=False, columns=[ TableColumn(field='Site',title='Site'),TableColumn(field='N',title='N'),TableColumn(field='R',title='R'),], width=200, height=75)
		
	# assign JS callbacks to the source and box selection tool
	source.js_on_change('data', CustomJS(args={'site_input':site_input,'site_input2':site_input2,'dt':data_table},code=correlation_code.replace('ROWID','0')))
	source.js_on_change('selected', CustomJS(args={'site_input':site_input,'site_input2':site_input2,'dt':data_table},code=correlation_code.replace('ROWID','0')))
	source2.js_on_change('data', CustomJS(args={'site_input':site_input,'site_input2':site_input2,'dt':data_table},code=correlation_code.replace('ROWID','1')))
	source2.js_on_change('selected', CustomJS(args={'site_input':site_input,'site_input2':site_input2,'dt':data_table},code=correlation_code.replace('ROWID','1')))

	fig.select_one(BoxSelectTool).callback = CustomJS(args={'txt':select_text},code = box_select_code)
	fig2.select_one(BoxSelectTool).callback = CustomJS(args={'txt':select_text},code = box_select_code)
	fig3.select_one(BoxSelectTool).callback = CustomJS(code = corfig_box_select_code)
	fig4.select_one(BoxSelectTool).callback = CustomJS(code = corfig_box_select_code)

## END OF INPUT WIDGETS CALLBACKS SECTION
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

	if layout_mode == 'comp':
		var_list = [i.value for i in [var_input,var_input2,var_input3,var_input4]]
	elif layout_mode =='simple':
		var_list = [var_input.value]

	# get all the data sources in a list
	if not public:
		data_list = [np.array(source.data['y1'])[np.array(source.data['colo'])!=flag_color].astype(np.float)] # data in red is all the data from public files and flag=0 data from private files
		if layout_mode == 'comp':
			data_list += [	np.array(source.data['y2'])[np.array(source.data['colo'])!=flag_color].astype(np.float),
							np.array(source2.data['y1'])[np.array(source2.data['colo'])!=flag_color].astype(np.float),
							np.array(source2.data['y2'])[np.array(source2.data['colo'])!=flag_color].astype(np.float),]
	else:
		data_list = [np.array(source.data['y1'])]
		if layout_mode == 'comp':
			data_list += [	np.array(source.data['y2']),
							np.array(source2.data['y1']),
							np.array(source2.data['y2']),]

	for ID,current_data in enumerate(data_list): # loop over the data sources
		# ID = 0 ; left Y axis of the first figure
		# ID = 1 ; left Y axis of the second figure
		# ID = 2 ; right Y axis of the first figure
		# ID = 3 ; right Y axis of the second figure

		if var_list[ID]=='': # if no variable is selected, go to next iteration of the for loop
			continue

		# in 'comp' mode, if the variable is the same on both y axes of a figure, add the two sites data sources before doing the scaling calculations
		if layout_mode == 'comp':
			if ID==0 and var_input.value==var_input3.value: # if the y axes of the first figure show the same variable
				current_data = np.append(current_data,data_list[2])
			elif ID==1 and var_input2.value==var_input4.value: # if the y axes of the second figure show the same variable
				current_data = np.append(current_data,data_list[3])

			# if the two sites data sources have been added together, the scaling has already been done with the first y axis, so make the second y axis range equal to the first
			elif ID==2 and var_input.value==var_input3.value: # if the y axes of the first figure show the same variable
				min_y = fig.y_range.start
				max_y = fig.y_range.end
				fig.extra_y_ranges['first_var'].start = min_y
				fig.extra_y_ranges['first_var'].end = max_y
				fig4.y_range.start = min_y
				fig4.y_range.end = max_y
				continue
			elif ID==3 and var_input2.value==var_input4.value: # if the y axes of the second figure show the same variable
				min_y = fig2.y_range.start
				max_y = fig2.y_range.end
				fig2.extra_y_ranges['second_var'].start = min_y
				fig2.extra_y_ranges['second_var'].end = max_y
				fig4.x_range.start = min_y
				fig4.x_range.end = max_y
				continue

		negatives = current_data[current_data<0] # negative values in the data
		positives = current_data[current_data>0] # positives values in the data
		zeroes = current_data[current_data==0]

		mean_current_data = np.mean(current_data)

		try:
			# I tried to define some rules to obtain a good scaling in different situations
			if len(zeroes)==len(current_data): # if all the data is exactly 0; the constant variable exception will handle it
 				pass
 			elif len(zeroes)>0.1*len(current_data): # if more than 10% of the data is exactly 0, it might be a dummy value, so I don't include those
 				current_data = current_data[current_data!=0]
			elif len(negatives)==0:
				current_data = current_data[(current_data<1.5*mean_current_data) & (current_data>0.5*mean_current_data)] # resample the data to get within +/- 50% of the mean
			elif len(positives)==0:
				current_data = current_data[(current_data>1.5*mean_current_data) & (current_data<0.5*mean_current_data)] # resample the data to get within +/- 50% of the mean
			elif 0.5<(len(positives)/len(negatives))<1.5:
				current_data = current_data[(current_data<2*np.mean(positives)) & (current_data>2*np.mean(negatives))] # resample the data to get within +20% of the positive mean and +20% of the negative mean

			if list(current_data).count(current_data[0])!=len(current_data): # if the variable is not a constant, start calculations

				min_y = min(current_data)
				max_y = max(current_data)

				ampli = (max_y-min_y)/10

				min_y = min_y - ampli
				max_y = max_y + ampli

				if ID == 0:
					fig.y_range.start = min_y
					fig.y_range.end = max_y
					if layout_mode == 'comp':
						fig3.y_range.start = min_y
						fig3.y_range.end = max_y

				elif ID == 1: # only happens when mode==comp
					fig2.y_range.start = min_y
					fig2.y_range.end = max_y				
					fig3.x_range.start = min_y
					fig3.x_range.end = max_y

				elif ID ==2: # only happens when mode==comp
					fig.extra_y_ranges['first_var'].start = min_y
					fig.extra_y_ranges['first_var'].end = max_y
					fig4.y_range.start = min_y
					fig4.y_range.end = max_y
				
				elif ID ==3: # only happens when mode==comp
					fig2.extra_y_ranges['second_var'].start = min_y
					fig2.extra_y_ranges['second_var'].end = max_y
					fig4.x_range.start = min_y
					fig4.x_range.end = max_y												
			else: # if the variable is a constant, update the 'status_div' to let the user know why nothing happened
				status_div.text = var_list[ID]+' is constant: {:3.2E}'.format(current_data[0])
		except IndexError:
			status_div.text = 'Cannot scale: missing data'

		dum_hide.value = str(time.time()) # update the css of labels
## END OF CENTER FUNCTION
#########################################################################################################################################################################
## LAYOUT PLOT ELEMENTS
center_button = Button(label='Scale without extrema',name='test',width=160) # button to 'center' the plot on the 'good' data
center_button.on_click(center) # assign the callback function to the button

# some blue lines that will be used to visually separate groups of widgets
linediv = Div(text='<hr width="100%" color="lightblue">',width=430)
linediv2 = Div(text='<hr width="100%" color="lightblue">',width=430)
linediv3 = Div(text='<hr width="100%" color="lightblue">',width=430)

fig.toolbar.tools += [NewResetTool()]

# put the figure by itself in a grid layout (I can better control where the toolbar will show if i do that)
if layout_mode == 'comp':

	# add the custom reset tool to figures
	fig2.toolbar.tools += [NewResetTool()]
	fig3.toolbar.tools += [NewResetTool()]
	fig4.toolbar.tools += [NewResetTool()]

	figrid = gridplot([[fig],[fig2]], toolbar_location = 'above')
	figrid2 = gridplot([[fig3,fig4]], toolbar_location = 'above')

	for curgrid in [figrid,figrid2]:
		curgrid.children[0].merge_tools = False # need to do that due to a bug in bokeh 0.12.6 that prevent gridplot to display the HoverTool in the toolbar
		curgrid.children[0].tools = [i for i in curgrid.children[0].tools if type(i)==bokeh.models.tools.HoverTool]

	dumdiv2 = Div(text='',width=50) # dummy div for spacing
	if not public:
		side_box = gridplot([[site_input,site_input2],[var_input,var_input3],[var_input2,var_input4],[linediv],[date_input,flag_input],[linediv2],[load_button],[status_div],[linediv3],[center_button,dumdiv2,hover_button],[select_text],[data_table],[notes_div]],toolbar_location=None)
	else:
		side_box = gridplot([[site_input,site_input2],[var_input,var_input3],[var_input2,var_input4],[linediv],[date_input],[linediv2],[load_button],[status_div],[linediv3],[center_button,dumdiv2,hover_button],[select_text],[data_table],[notes_div]],toolbar_location=None)

	figroup = gridplot([[figrid],[figrid2]],toolbar_location='left')
elif layout_mode == 'simple':
	figroup = gridplot([[fig]], toolbar_location = 'left')
	figroup.children[0].merge_tools = False
	if not public:
		side_box = gridplot([[site_input],[var_input],[linediv],[date_input,flag_input],[linediv2],[load_button],[status_div],[linediv3],[center_button],[notes_div]],toolbar_location=None)
	else:
		side_box = gridplot([[site_input],[var_input],[linediv],[date_input],[linediv2],[load_button],[status_div],[linediv3],[center_button],[notes_div]],toolbar_location=None)

# final layout
grid = gridplot([[figroup,side_box],[dum_box]], toolbar_location = None)

## END OF LAYOUT PLOT ELEMENTS
#########################################################################################################################################################################
	
curdoc().title='TCCON' # this changes the title of the internet tab in which the document will be displayed
curdoc().add_root(grid) # this add the grid layout to the document
curdoc().add_timeout_callback(hide_dummy,2000) # this schedules the 'hide_dummy' callback to be triggered 1000 milliseconds after page load; if this fails to trigger too often, increase the number