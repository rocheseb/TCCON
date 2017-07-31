#!/var/lib/py27_sroche/bin/python
 # -*- coding: utf-8 -*-

from __future__ import print_function # allows the use of Python 3.x print(function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# Code description #
####################

'''
Read TCCON data from .eof, .eof.csv, or .ncdf files in a given folder. 
Files must be from the same site !  
File names must be such that they appear in chronological order if sorted by name.

Produce .html plots of the data in /path/to/folder/SAVE

How to use:

Can be used with or without commandline arguments. If no arguments are given, the user will be asked to provide the necessary input after the code starts.

python TCCON_eof.py /path/to/folder flag

- the first argument is the path to the folder containing the TCCON files, (.eof, .eof.csv, or .nc)
- the second argument is the flag that can be a number or 'all'. If the flag is 'all', all the data will be read. Otherwise only the data with a specific flag will be read.

the last 3 lines of the code should be modified to customize the .html file name or the chain of character that appears in the internet tab when you open the plot.

How to modify the plot structure:

Two python objects can be edited by the user to customize the plot, they are in the section "modify here"

- colors_dict : dictionary that associate a color to a key, all variables that contain a given key will be ploted with the associated color

- bok_struct :  this is the crucial part of the code. The structure of this ordered dictionary controls the structure of the plot ( Number of panels, figures in panels, and lines plotted in figures ).

- you need to know the exact variables names of the variables you want to plot. Be careful as some variables don't have the exact same name in .eof and .nc files ...


TIPS:

if you think setting up the bok_struct dictionary is too confusing/complicated, you should try just using 'Key' panels. 
Key panels include 'Key' in the panel name, and their attribute 'lines' is a tuple containing 1 or several keywords.
All variables that include the keywords in their name will be read and available for plotting in two figures, the y axis from the two figures will be plotted against each other in a third figure.

For example, to have a one panel plot with all variables that include 'co2_6220' or 'x':

bok_struct =  OrderedDict([
			('Key_test',OrderedDict([
					('custom',{
								'lines':('co2_6220','x',),
								'plot_height':250,
								'plot_width':1000,
								}),	
					])),
				])

You should probably create a file to save different bok_struct objects that you may find useful
'''

#############
# Libraries #
#############

#general
import os
import sys

#special arrays with special functions
import numpy as np

#time handling
import time
from datetime import datetime,timedelta
import calendar

#netcdf reader
import netCDF4

#round up
from math import ceil

# html plots
from bokeh.plotting import figure
from bokeh.models import Panel, Tabs, CustomJS, ColumnDataSource, RadioGroup, VBox, Div, BoxSelectTool, DataTable, TableColumn, AutocompleteInput, Select
from bokeh.layouts import gridplot, widgetbox
from bokeh.resources import CDN
from bokeh.embed import file_html

# dictionaries with sorted keys
import collections
from collections import OrderedDict

# to ignore warnings (like for the save_source which is just meant to "pass" python variables to a callback but will return warnings because of columns of different sizes)
import warnings
warnings.filterwarnings('ignore')

###############
# Modify here #
###############

# the name of the html file that this code produces, make sure it finishes with '.html'
save_name = 'test.html'

# tab_name will appear in the internet tab when you open the page
tab_name = 'TCCON'

# associate a keyword with a color. Variables including the keyword will be plotted with that color. (be careful with co, co2 and o2 ! )
# note: if you use 'all' for the flag, only flag 0 data will use the colors, flag != 0 data will be grey
colors_dict = {
				'co2':'red',
				'n2o':'orange',
				'ch4':'green',
				'co_':'black',
				'hdo':'cyan',
				'h2o':'blue',
				'hf':'purple',
				'air':'pink',
				'_o2':'firebrick',
				'hcl':'goldenrod',
				'sza':'fuchsia',
				}

# dictionary with a structure that will define the final plot layout.
'''
-> first key => Panel
	-> subkey 1 => figure in panel
		-> 'lines' => variables to plot in figure (specific lines/scatter lists)
'''
			# this indent is for panels
					# this indent is for figures under a panel (the key will be used as figure title)
								# this indent if for lines in a figure, given in a list named 'lines' (the keys will be used as Y axis label and must be exact variable names from the EOF file header or from the .nc file)

# to add an error plot below regular plots, set the 'errlines' parameter to True.
# this only works for variables that have a '_error' extension (e.g. xco2_ppm, xco2_ppm_error)								
bok_struct = OrderedDict([
			('Key_co2',OrderedDict([
					('custom',{
								'lines':('x','vsf','sza','column','S_G','fvsi','p_out','h_out','t_out'),
								'plot_height':250,
								'plot_width':800,
								}),	
					])), # end of 'Key_co2_6220' panel
			])
# below is an example of a more complex bok_struct object containing different kinds of tabs
"""
bok_struct = OrderedDict([
			('xgas',OrderedDict([	
					('xair',{
								'lines':['xair'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),
					('xhf_ppt',{
								'lines':['xhf_ppt'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),
					('xh2o_ppm',{
								'lines':['xh2o_ppm'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),
					('xhdo_ppm',{
								'lines':['xhdo_ppm'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),
					('xco_ppb',{
								'lines':['xco_ppb'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),
					('xch4_ppm',{
								'lines':['xch4_ppm'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),																																								
					('xco2_ppm',{
								'lines':['xco2_ppm'],
								'errlines':False,
								'plot_height':125,
								'plot_width':1000,
								}),
					])), # end of 'xgas' Panel

			# the 'Custom' panel is a special case, it can only have one figure with several lines in it.
			('Custom_xgas',OrderedDict([
					('custom',{
								'lines':['xair','xhf_ppt','xh2o_ppm','xhdo_ppm','xco_ppb','xn2o_ppb','xch4_ppm','xco2_ppm'],
								'errlines':True,
								'plot_height':500,
								'plot_width':1000,
								}),	
					])), # end of Custom_xgas panel
			('Custom_col',OrderedDict([
					('custom',{
								'lines':['column_hf','column_h2o','column_hdo','column_co','column_n2o','column_ch4','column_co2'],
								'errlines':True,
								'plot_height':500,
								'plot_width':1000,
								}),	
					])), # end of Custom_col panel
			('Custom_VSF',OrderedDict([
					('custom',{
								'lines':['vsf_air','vsf_hf','vsf_h2o','vsf_hdo','vsf_co','vsf_n2o','vsf_ch4','vsf_co2'],
								'errlines':True,
								'plot_height':500,
								'plot_width':1000,
								}),	
					])), # end of 'Custom_VSF' panel

			# panels that include 'Key' in their name are another special case, the 'lines' attribute must be a tuple of keywords
			# all variables that include each keyword will be availabe (e.g. ('co2_6220') will get 25 variables !)
			# 'Key' panels do not have an 'errlines' parameter
			# the search for the given keys in variables names is CASE SENSITIVE
			('Key_co2',OrderedDict([
					('custom',{
								'lines':('co2',),
								'plot_height':250,
								'plot_width':1000,
								}),	
					])), # end of 'Key_co2_6220' panel
			('other',OrderedDict([
					('fvsi_%',{
								#'lines':['fvsi_%'], #eof
								'lines':['fvsi'], #netcdf
								'plot_height':150,
								'plot_width':1000,
								'errlines':False,
								}),			

					('hout_%RH',{
								#'lines':['hout_%RH'], #eof
								'lines':['hout_RH'], #netcdf
								'plot_height':150,
								'plot_width':1000,
								'errlines':False,
								}),														
					('wspd_m/s',{
								#'lines':['wspd_m/s'], #eof
								'lines':['wspd_m_s'], #netcdf
								'plot_height':150,
								'plot_width':1000,
								'errlines':False,
								}),	
					('wdir_deg',{
								'lines':['wdir_deg'],
								'plot_height':150,
								'plot_width':1000,
								'errlines':False,
								}),
					('pout_hPa',{
								'lines':['pout_hPa'],
								'plot_height':150,
								'plot_width':1000,
								'errlines':False,
								}),	
					])), # end of 'other' panel
			])
"""
#############
# Functions #
#############

def merged_tccon_data(path,diag_var=[],diag_key=[],flag=''):
	'''
	reads 
	'''

	file_list=sorted([i for i in os.listdir(path) if '.' in i]) # list of file_list in the given path

	# determine the file type for the read_tccon function
	if '.nc' in file_list[0]:
		mode = 'netcdf'
	if '.eof' in file_list[0]:
		mode = 'eof'
	if ('.nc' not in file_list[0]) and ('.netcdf' not in file_list[0]) and ('.eof' not in file_list[0]):
		print('You need netCDF or .eof tccon file_list to run this program')
		sys.exit()

	# loop over file_list and merge all the data in a dictionary, this will skip any data that is not properly time sorted
	print('\nGetting data:')
	all_files_data = {}
	for tccon_file in file_list:

		print('\t-',tccon_file)

		file_path = os.path.join(path,tccon_file)

		file_data = read_tccon(file_path,mode=mode,variables=diag_var,key_variables=diag_key,flag=flag)

		if len(all_files_data.keys())==0:
			for key,value in file_data.iteritems():
				all_files_data[key] = np.array(value)
		else:
			for key in ['xtime']+[i for i in file_data if i!='xtime']:
				for time_id,new_time in enumerate(file_data['xtime']):
					if new_time>all_files_data['xtime'][len(file_data)+time_id-1]:
						all_files_data[key] = np.append(all_files_data[key],file_data[key][time_id])
					else:
						print('Time overlap:',new_time,'<',all_files_data['xtime'][-1])

	return {key:value for key,value in all_files_data.iteritems()}

def read_tccon(path,mode='eof',variables=[],key_variables=[],flag='all'):
	'''
	 can read tccon .eof, .eof.csv, or .nc file_list with specified variables of two types:
	 "variables" is a list of exact variables names to be read
	 "key_variables" is a list of keywords, every variable that includes one of the keywords in their name will be read
	 if those lists are not specified all variables (>1200) will be read
	 flag is "all" by default and data with any flag will be read, you can give the appropriate integer if you want to only read data with a specific flag
	 the mode can be set to 'eof', 'csv', or 'netcdf'
	'''

	DATA = {}

	if mode.lower() not in ['eof','netcdf']:
		print('in function read_tccon() : The file format was not recognized')
		return DATA

	if mode.lower() == 'eof':
		infile = open(path,'r')
		content = infile.readlines()
		infile.close()

		if 'csv' in path:
			header = content[2].split(',')
		else:
			header = content[2].split()

		if variables == []:
			variables = header

		# add all the variables that includes keys in key_variables, and make sure there are no repeat variables
		add_var = []
		if key_variables != []:
			for key in key_variables:
				add_var += [var for var in header if key in var]
			for var in add_var:
				if var not in variables:
					variables.append(var)

		if flag == 'all':
			if 'csv' in path:
				content_T = np.array([line.split(',') for line in content[3:]]).T
			else:
				content_T = np.array([line.split() for line in content[3:]]).T
		else:
			if 'csv' in path:
				content_T = np.array([line.split(',') for line in content[3:] if line.split(',')[0]==flag]).T
			else:
				content_T = np.array([line.split() for line in content[3:] if line.split()[0]==flag]).T

		for var in variables:
			if var not in header:
				print('in function read_tccon(): wrong variable input:',var,'not in',path)
				return {}				
			try:
				DATA[var] = np.array( content_T[header.index(var)], dtype=np.float )
			except ValueError:
				print('Skipping string variable:',var)
				
	if mode.lower() == 'netcdf':
		f = netCDF4.Dataset(path,'r')

		header = [v for v in f.variables]

		if variables == []:
			variables = header

		# add all the variables that includes keys in key_variables, and make sure there are no repeat variables
		add_var = []
		if key_variables != []:
			for key in key_variables:
				add_var += [var for var in header if ((key in var) and ('prio' not in var) and ('ak_' not in var))]
			for var in add_var:
				if var not in variables:
					variables.append(var)

		for var in variables:
			try:
				if flag == 'all':
					DATA[var] = np.array( f.variables[var][:], dtype=np.float )
				else:
					DATA[var] = np.array( [f.variables[var][ID] for ID,elem in enumerate(list(f.variables['flag'][:])) if int(elem)==int(flag)] ,dtype=np.float ) # this is kind of slow ...
			except ValueError:
				'Skipping string variable:',var

		f.close()

	DATA['xtime'] = [datetime(int(DATA['year'][i]),1,1)+timedelta(days=int(DATA['day'][i])-1,hours=float(DATA['hour'][i])) for i in range(len(DATA['hour']))]

	del DATA['year']
	del DATA['day']
	del DATA['hour']

	return DATA

def descend_values(dic):
	'''
	gets all the values in a dictionnary and put them in a 1d list
	'''
	for key in sorted(dic.keys()):
		if type(dic[key]) in [dict,collections.OrderedDict]:
			for x in descend_values(dic[key]):
				yield x
		else:
			yield dic[key]

def flatten(obj,keep=[]):
	'''
	gets all the items in a nested iterable oject and put them in a 1d list
	keep is a list of object types that should not be flattened ( in python 3 strings will be "exploded" if str is not in the keep list)
	'''
	if (type(obj) in [dict,collections.OrderedDict]) and ((dict not in keep) and (collections.OrderedDict not in keep)):
		for x in descend_values(obj):
			yield x
	else:	
		for item in obj:
			if hasattr(item,'__iter__') and (type(item) not in keep):
				for x in flatten(item,keep=keep):
					yield x
			else:
				yield item

#########
# Setup #
#########

######## /!\ important time handling to make sure times don't get shifted from UTC due to computer environment variables when using datetime objects ###########
os.environ['TZ'] = 'UTC'
time.tzset()
# This will not change your system timezone settings
################################################################################################################################################################

argu = sys.argv

path='1'
while os.path.isdir(path)==False:
	if len(argu)>1:
		path = argu[1]
	else:
		print('\n\nPut in a folder all the .eof files you want to process, give the path to that folder')
		path=raw_input('Give the path to your folder /YOUR/PATH/TO/FOLDER  :\n')
	if os.path.isdir(path)==False:
		print('/!\\ You gave a wrong path /!\\\n')
		if len(argu)>1:
			sys.exit()

save_path=os.path.join(path,'SAVE')
if not os.path.isdir(save_path):
	os.makedirs(save_path)

print('\nPlots will be saved in',save_path,'\n')

flag_check = True
while flag_check:
	flag = 'all'
	if len(argu)>2:
		flag = argu[2]
	else:
		flag = raw_input('Select data flag to be used (type "all" to use all data)\n')
	if flag != 'all':
		try:
			test = int(flag)
			flag_check = False
		except ValueError:
			print('The flag must be a number (or blank "")')
			flag_check = True
			if len(argu)>2:
				sys.exit()
	else:
		flag_check = False

# callback javascript codes
err_code = """
all = S_all.data;
main = S_main.data;
err = S_err.data;
varlist = S_save.data["varlist"];
colors = S_save.data["colors"][0];

var val = cb_obj.active;
var vartoplot = varlist[val];

main_laby.axis_label = vartoplot;

for (var key in colors) {
if(vartoplot.includes(key)){
var colo = colors[key]
}
}

var y = all[vartoplot];
var min = Math.min.apply(null,y);
var max = Math.max.apply(null,y);
var ampli = max - min;
mainy.start = min - 0.1*ampli;
mainy.end = max + 0.1*ampli;

var yer = all[varlist[val]+'_error'];
var min = Math.min.apply(null,yer);
var max = Math.max.apply(null,yer);
var ampli = max - min;
erry.start = min - 0.1*ampli;
erry.end = max + 0.1*ampli;

main["y"] = [];
main["colo"] = [];
err["y"] = [];

for (i=0;i<y.length;i++) {
main["y"].push(all[varlist[val]][i]);
err["y"].push(all[varlist[val]+'_error'][i]);

if (all["flag"][i]=="0") {main["colo"].push(colo);} else {main["colo"].push("grey");}
	
}

S_main.change.emit();
S_err.change.emit();
"""

code = """
all = S_all.data;
main = S_main.data;
varlist = S_save.data["varlist"];
colors = S_save.data["colors"][0];

var val = cb_obj.active;
var vartoplot = varlist[val];

main_laby.axis_label = vartoplot;

for (var key in colors) {
if(vartoplot.includes(key)){
var colo = colors[key]
}
}

var y = all[vartoplot];
var min = Math.min.apply(null,y);
var max = Math.max.apply(null,y);
var ampli = max - min;
mainy.start = min - 0.1*ampli;
mainy.end = max + 0.1*ampli;

main["y"] = [];
main["colo"] = [];

for (i=0;i<y.length;i++) {
main["y"].push(all[vartoplot][i]);

if (all["flag"][i]=="0") {main["colo"].push(colo);} else {main["colo"].push("grey");}
	
}

S_main.change.emit();
"""

input_code= """
var all = S_all.data;
var main = S_main.data;
var fill = S_fill.data;
var varlist = S_save.data["varlist"];
var colors = S_save.data["colors"][0];
var colo = '';

var vartoplot = cb_obj.value;

console.log(vartoplot);

main_laby.axis_label = vartoplot;
fill_lab.axis_label = vartoplot;

for (var key in colors) {
if(vartoplot.includes(key)){
colo = colors[key]
} else {colo = 'red'}
}

var y = all[vartoplot];
var min = Math.min.apply(null,y);
var max = Math.max.apply(null,y);
var ampli = max - min;
mainy.start = min - 0.1*ampli;
mainy.end = max + 0.1*ampli;

main["y"] = [];
main["colo"] = [];

fill["colo"] = [];
if (cb_obj.title.includes("1")) {
fill["y"]=[];
}
if (cb_obj.title.includes("2")) {
fill["x"]=[];
}

for (i=0;i<y.length;i++) {
main["y"].push(all[vartoplot][i]);

if (all["flag"][i]=="0") {main["colo"].push(colo);fill["colo"].push("red");} else {main["colo"].push("grey");fill["colo"].push("grey");}

if (cb_obj.title.includes("1")) {fill["y"].push(all[vartoplot][i]);}
if (cb_obj.title.includes("2")) {fill["x"].push(all[vartoplot][i]);}
	
}

console.log(fill["y"]);
console.log(fill["x"]);

S_main.change.emit();
S_fill.change.emit();
"""

key_source_code = """
var inds = cb_obj.selected['1d'].indices;
var d1 = cb_obj.data;
var d2 = s2.data;
var tab = dt.source.data;
var dcor = S_fill.data;

s2.selected['1d'].indices = inds;
S_fill.selected['1d'].indices = inds;

S_fill.change.emit();
s2.change.emit();

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
	ym1 += d1['y'][inds[i]];
	ym2 += d2['y'][inds[i]];
}

ym1 /= inds.length;
ym2 /= inds.length;

for (i=0; i < inds.length; i++){
	T1 += (d1['y'][inds[i]] - ym1)*(d2['y'][inds[i]] - ym2);
	T2 += Math.pow(d1['y'][inds[i]] - ym1,2);
	T3 += Math.pow(d2['y'][inds[i]] - ym2,2);
}

tab['R'][0] = (T1/Math.sqrt(T2*T3)).toFixed(3);

dt.change.emit();
"""

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
"""

key_notes = """
<font size=4><b>Notes:</b></font><font size=2></br>
</br>
Type in the input widgets on the right to select the variable to display in Figure 1 and Figure 2</br>
The input widgets use auto completion to suggest available variables, you can double click on suggestions to fill in the entry faster</br>
</br>
Use the "Box Select" tool to select data in Figure 1, number of points (N) and correlations (R) will be shown in the table</br>
</br>
Data in grey has a flag != 0</br></br>If the axis labels do not update, use the "Reset" tool</br>
</br>
<img src="https://upload.wikimedia.org/wikipedia/commons/2/24/Warning_icon.svg" height="20" width="20"> The "Save" tool from the toolbar will save each figure to a different .png</font>
"""

# flatten the bok_struct dictionary to get all the different variable names that should be read in the eof/netcdf files
bok_struct_values = list(descend_values(bok_struct))
all_var = list(flatten([i for i in bok_struct_values if type(i)==list])) # all the variables names that should be read

all_key = list(flatten([i for i in bok_struct_values if type(i)==tuple])) # all the variables key, if a key is 'co2', all variables including 'co2' in their name will be added

# if a figure has 'errlines' set to true, get the corresponding error variables, this only work for variables that have and '_error' extension !
for val_id,val in enumerate(bok_struct_values):
	if val is True:
		all_var += [var+'_error' for var in bok_struct_values[val_id+1]]

# get rid of repeat variables
diag_var = []
for var in all_var:
	if var not in diag_var:
		diag_var.append(var)

# get rid of repeat keys
diag_key = []
for var in all_key:
	if var not in diag_key:
		diag_key.append(var)

diag_var += ['year','day','hour','flag'] # append some default variables to be read. So we can create datetime objects, and also the quality flags

# tools for the plotting
TOOLS = "pan,wheel_zoom,box_zoom,undo,redo,reset,save"

# special bokh object to store data inside the HTML page
all_source = ColumnDataSource(data=merged_tccon_data(path=path,diag_var=diag_var,diag_key=diag_key,flag=flag), id='all_source')

main_source_list = [] # this source will be empty and filled from all_source via callbacks
err_source_list = [] # this source will be empty and filled from all_source via callbacks
save_source_list = [] # this will just contain python objects I want to pass to the javascript callbacks

#############
# Main code #
#############

# making the plot based on the bok_struct dictionary

tabs = [] # the different panels in the final bokeh plot
for panel_key in bok_struct:
	fig_list = [] # the different figures in the panel
	if True not in [elem in panel_key for elem in ['Custom','Key','Diag']]: # general case
		for fig_key in bok_struct[panel_key]:

			width = bok_struct[panel_key][fig_key]['plot_width']
			height = bok_struct[panel_key][fig_key]['plot_height']
			
			if panel_key==bok_struct.keys()[0] and len(fig_list)==0:
				fig_list.append( figure(output_backend = "webgl", title=fig_key,plot_width=width,plot_height=height,x_axis_type='datetime',tools=TOOLS) )
				save_fig = fig_list[0]
			else: 
				fig_list.append( figure(output_backend = "webgl", title=fig_key,plot_width=width,plot_height=height,x_axis_type='datetime',x_range=save_fig.x_range,tools=TOOLS) ) 

			if bok_struct[panel_key][fig_key]['errlines'] is True:
				fig_list.append( figure(output_backend = "webgl", plot_width=width,plot_height=int(height/2),x_axis_type='datetime',x_range=save_fig.x_range,tools=TOOLS ) )

			for plot_var in bok_struct[panel_key][fig_key]['lines']:
				colo = [colors_dict[key] for key in colors_dict if key in plot_var]

				if bok_struct[panel_key][fig_key]['errlines'] is True:
					if len(colo)==1:
						colo = colo[0]
						fig_list[-2].scatter(x='xtime',y=plot_var,color=colo,source=all_source)
					else:
						fig_list[-2].scatter(x='xtime',y=plot_var,source=all_source)
				
					fig_list[-1].scatter(x='xtime',y=plot_var+'_error',color='black',source=all_source)
					fig_list[-1].yaxis.axis_label = 'Error'

					fig_list[-2].yaxis.axis_label = plot_var
				else:
					if len(colo)==1:
						colo = colo[0]
						fig_list[-1].scatter(x='xtime',y=plot_var,color=colo,source=all_source)
					else:
						fig_list[-1].scatter(x='xtime',y=plot_var,source=all_source)

					fig_list[-1].yaxis.axis_label = plot_var

			grid = gridplot( [[fig] for fig in fig_list],toolbar_location='left',tools=TOOLS )

	if 'Custom' in panel_key: # first special case, custom plots

		main_source_list.append( ColumnDataSource(data={"x":all_source.data['xtime'],"y":[],"colo":[]}) )
		if True in list(descend_values(bok_struct[panel_key])):
			err_source_list.append( ColumnDataSource(data={"x":all_source.data['xtime'],"y":[]}) )
		save_source_list.append( ColumnDataSource(data={"colors":[colors_dict]}) )

		for fig_key in bok_struct[panel_key]:
			width = bok_struct[panel_key][fig_key]['plot_width']
			height = bok_struct[panel_key][fig_key]['plot_height']

			# the figure in the custom panel is not linked to other panels figures

			fig_list.append( figure(output_backend = "webgl", plot_width=width,plot_height=height,x_axis_type='datetime',tools=TOOLS,y_axis_label='x_gas') )

			if bok_struct[panel_key][fig_key]['errlines'] is True:
				fig_list.append( figure(output_backend = "webgl", plot_width=width,plot_height=int(height/5),x_axis_type='datetime',x_range=fig_list[0].x_range,tools=TOOLS) )

			var_list = bok_struct[panel_key][fig_key]['lines']

			save_source_list[-1].data['varlist'] = var_list

			for plot_var in var_list:

				fig_list[0].scatter(x="x",y="y",color='colo',source=main_source_list[-1])

				if bok_struct[panel_key][fig_key]['errlines'] is True:
					fig_list[1].scatter(x="x",y="y",color='black',source=err_source_list[-1])				
					fig_list[1].yaxis.axis_label='Error'

			fig_list[0].xaxis.axis_label='Time'
			fig_list[0].yaxis.axis_label='x_gas'
			
			if bok_struct[panel_key][fig_key]['errlines'] is True:
				callback=CustomJS(	args=dict(
											S_all=all_source,
											S_main=main_source_list[-1],
											S_err=err_source_list[-1],
											S_save=save_source_list[-1],
											mainy=fig_list[0].y_range,
											erry=fig_list[1].y_range,
											main_laby=fig_list[0].yaxis[0],
											),  
									code=err_code)
			else:
				callback=CustomJS(	args=dict(
											S_all=all_source,
											S_main=main_source_list[-1],
											S_save=save_source_list[-1],
											mainy=fig_list[0].y_range,
											main_laby=fig_list[0].yaxis[0],
											),  
									code=code)

			radiogroup = RadioGroup(labels=var_list,active=0,callback=callback)
			radiobox = VBox(radiogroup,width=100)

			if bok_struct[panel_key][fig_key]['errlines'] is True:
				grid = gridplot([[fig_list[1]],[fig_list[0],radiobox]],toolbar_location='left',tools=TOOLS)
			else:
				grid = gridplot([[fig_list[0],radiobox]],toolbar_location='left',tools=TOOLS)

	if 'Key' in panel_key: # second special case, key plots

		fig_key = [i for i in bok_struct[panel_key]][0] # there is only one key in "Key" panels
		key_list = bok_struct[panel_key][fig_key]['lines']
		var_list = list(flatten([[var for var in sorted(all_source.data.keys()) if key in var] for key in key_list]))
		key_tools = "box_zoom,wheel_zoom,box_select,pan,undo,redo,reset,save"

		# sources
		main_source_list.append( ColumnDataSource(data={"x":np.array([]),"y":np.array([]),"colo":[]}) )
		main_source_list.append( ColumnDataSource(data={"x":all_source.data['xtime'],"y":np.array([]),"colo":np.array([])}) )
		main_source_list.append( ColumnDataSource(data={"x":all_source.data['xtime'],"y":np.array([]),"colo":np.array([])}) )
		save_source_list.append( ColumnDataSource(data={"colors":[colors_dict],"varlist":var_list}) )
		table_source = ColumnDataSource( data = {'N':[0],'R':[0]} ) # the data source of the table

		# figures
		width = bok_struct[panel_key][fig_key]['plot_width']
		height = bok_struct[panel_key][fig_key]['plot_height']
		fig_list.append( figure(output_backend = "webgl", title='Figure 1',plot_width=width,plot_height=height,x_axis_type='datetime',tools=key_tools) )
		fig_list[-1].select_one(BoxSelectTool).dimensions = 'width'
		fig_list.append( figure(output_backend = "webgl", title='Figure 2',plot_width=width,plot_height=height,x_axis_type='datetime',x_range=fig_list[0].x_range,tools=key_tools) )
		fig_list[-1].select_one(BoxSelectTool).dimensions = 'width'
		fig_list.append( figure(output_backend = "webgl", title='Fig1 y VS Fig2 y',plot_width=400,plot_height=400) )

		# glyphrenderers
		fig_list[0].scatter(x="x",y="y",color='colo',source=main_source_list[-1])
		fig_list[1].scatter(x="x",y="y",color='colo',source=main_source_list[-2])
		fig_list[2].scatter(x="x",y="y",color='colo',source=main_source_list[-3])

		# widgets
		menu = [(plot_var,plot_var) for plot_var in var_list]
		#input_0 = AutocompleteInput(title="Figure 1:", value=None,completions=sorted(all_source.data.keys()),width=200)
		#input_1 = AutocompleteInput(title="Figure 2:", value=None,completions=sorted(all_source.data.keys()),width=200)
		input_0 = Select(title="Figure 1:", value=None,options=sorted(all_source.data.keys()),width=200)
		input_1 = Select(title="Figure 2:", value=None,options=sorted(all_source.data.keys()),width=200)			
		data_table = DataTable(source=table_source, columns=[ TableColumn(field='N',title='N'),TableColumn(field='R',title='R'),], width=200, height=55)
		select_text = Div(text='',width = 450) # text div that will be updated with the selected range of date within the BoxSelect tool
		notes = Div(text=key_notes,width=600)
		dum = Div(text='',height=70) # dummy div widget for widget spacing in the layout
		dum2 = Div(text='',height=70) # dummy div widget for widget spacing in the layout

		# callbacks
		main_source_list[-1].callback = CustomJS(args = dict(s2=main_source_list[-2],dt=data_table,S_fill=main_source_list[-3]), code=key_source_code) # update selection in other plots and fill the table after using the box_select tool
		main_source_list[-2].callback = CustomJS(args = dict(s2=main_source_list[-1],dt=data_table,S_fill=main_source_list[-3]), code=key_source_code) # update selection in other plots and fill the table after using the box_select tool

		fig_list[0].select_one(BoxSelectTool).callback = CustomJS(args=dict(txt=select_text),code=box_select_code) # make the BoxSelect tool update the 'txt' Div widget with the currently selected range of dates.
		fig_list[1].select_one(BoxSelectTool).callback = CustomJS(args=dict(txt=select_text),code=box_select_code) # make the BoxSelect tool update the 'txt' Div widget with the currently selected range of dates.

		input_0.callback=CustomJS(	args=dict(
									S_all=all_source,
									S_main=main_source_list[-1],
									s2=main_source_list[-2],
									S_fill=main_source_list[-3],
									S_save=save_source_list[-1],
									mainy=fig_list[0].y_range,
									main_laby=fig_list[0].yaxis[0],
									fill_lab=fig_list[2].yaxis[0],
									dt = data_table,
									),  
							code=input_code+key_source_code.replace('cb_obj','S_main'))

		input_1.callback=CustomJS(	args=dict(
									S_all=all_source,
									S_main=main_source_list[-2],
									s2=main_source_list[-1],
									S_fill=main_source_list[-3],
									S_save=save_source_list[-1],
									mainy=fig_list[1].y_range,
									main_laby=fig_list[1].yaxis[0],
									fill_lab=fig_list[2].xaxis[0],
									dt = data_table,
									),  
							code=input_code+key_source_code.replace('cb_obj','S_main'))

		# layout the final grid
		notebox = widgetbox(select_text,data_table,notes,width=650)
		dropbox_0 = widgetbox(dum,input_0,width=210) # I use the dummy div widget to have the input button ~aligned with the center of the figure 
		dropbox_1 = widgetbox(dum2,input_1,width=210)

		grid = gridplot([[fig_list[0],dropbox_0],[fig_list[1],dropbox_1],[fig_list[2],notebox]],toolbar_location='left')

	tabs.append( Panel(child=grid,title=panel_key) )

if len(tabs) == 1:
	final = tabs[0].child
else:
	final=Tabs(tabs=tabs)

print('\nWritting',save_name,'...')
outfile=open(os.path.join(save_path,save_name),'w')
outfile.write(file_html(final,CDN,tab_name))
outfile.close()