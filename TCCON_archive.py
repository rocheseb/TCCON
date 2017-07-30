####################
# Code description #
####################
"""
This code will produce an interactive plots reading into TCCON netcdf files to display time series of the data.

In the same directory as this code, make a 'TCCON' folder with netcdf files from the TCCON public archive (http://tccon.ornl.gov/ 566MB for all the files)

This code needs to be run by a bokeh server like this:

bokeh serve --show TCCON_archive.py
"""

import os
import netCDF4
from datetime import datetime
import time
import calendar
import numpy as np

from bokeh.io import show,curdoc,save
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, TextInput, Div, CustomJS, AutocompleteInput, Button, TextInput, Select
from bokeh.layouts import gridplot, widgetbox

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

tccon_path = os.path.join(os.getcwd(),'TCCON')
tccon_file_list = os.listdir(tccon_path)

txt = Div(text='',width=200)

source = ColumnDataSource(data={'x':[],'y':[]})

fig = figure(output_backend="webgl",plot_width=600,plot_height=400,x_axis_type='datetime')
fig.xaxis.axis_label = 'Time'

fig.scatter(x='x',y='y',source=source)

#site_input = AutocompleteInput(title='Site:',completions = T_FULL.values())
#var_input = AutocompleteInput(title='Variable to plot:')
site_input = Select(title='Site:',options = sorted(T_FULL.values()))
var_input = Select(title='Variable to plot:')
date_input = TextInput(title='start-end yyyymmdd:')

def set_site(attr,old,new):
	site = site_input.value
	prefix = [key for key in T_FULL if T_FULL[key]==site][0]
	site_file = [i for i in tccon_file_list if prefix in i][0]
	f = netCDF4.Dataset(os.path.join(tccon_path,site_file),'r')
	#var_input.completions = [v for v in f.variables]
	var_input.options = sorted([v for v in f.variables])
	if var_input.value != '':
		nctime =  f.variables['time'][:] # fractional days since 1970
		try:
			mindate = date_input.value.split('-')[0]
			mindate = calendar.timegm(datetime(int(mindate[:4]),int(mindate[4:6]),int(mindate[6:8])).timetuple())/24/3600
		except:
			mindate = calendar.timegm(datetime(1970,1,1).timetuple())/24/3600

		try:
			maxdate = date_input.value.split('-')[1]
			maxdate = calendar.timegm(datetime(int(maxdate[:4]),int(maxdate[4:6]),int(maxdate[6:8])).timetuple())/24/3600
		except:
			maxdate = calendar.timegm(datetime(2050,1,1).timetuple())/24/3600

		newtime = nctime[(nctime>=mindate) & (nctime<maxdate)]

		start_id = np.where(nctime==newtime[0])[0][0]
		end_id = np.where(nctime==newtime[-1])[0][0]+1

		try:
			source.data = {'x':np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in newtime]),'y':f.variables[var_input.value][start_id:end_id]}
			fig.yaxis.axis_label = var_input.value
		except:
			pass
	f.close()
	fig.title.text = T_FULL[prefix]+', '+T_LOC[prefix]

def read_nc(attr,old,new):
	try:
		site = site_input.value
		prefix = [key for key in T_FULL if T_FULL[key]==site][0]
		site_file = [i for i in tccon_file_list if prefix in i][0]
		f = netCDF4.Dataset(os.path.join(tccon_path,site_file),'r')
	except NameError:
		pass
	else:
		nctime =  f.variables['time'][:] # fractional days since 1970
		try:
			mindate = date_input.value.split('-')[0]
			mindate = calendar.timegm(datetime(int(mindate[:4]),int(mindate[4:6]),int(mindate[6:8])).timetuple())/24/3600
		except:
			mindate = calendar.timegm(datetime(1970,1,1).timetuple())/24/3600

		try:
			maxdate = date_input.value.split('-')[1]
			maxdate = calendar.timegm(datetime(int(maxdate[:4]),int(maxdate[4:6]),int(maxdate[6:8])).timetuple())/24/3600
		except:
			maxdate = calendar.timegm(datetime(2050,1,1).timetuple())/24/3600

		newtime = nctime[(nctime>=mindate) & (nctime<maxdate)]

		start_id = np.where(nctime==newtime[0])[0][0]
		end_id = np.where(nctime==newtime[-1])[0][0]+1

		try:
			source.data = {'x':np.array([datetime(*time.gmtime(i*24*3600)[:6]) for i in newtime]),'y':f.variables[var_input.value][start_id:end_id]}
			fig.yaxis.axis_label = var_input.value
		except:
			pass

		f.close()

site_input.on_change('value',set_site)
var_input.on_change('value',read_nc)
date_input.on_change('value',read_nc)

grid = gridplot([[fig,widgetbox(site_input,var_input,date_input,width=300)]], toolbar_location = 'left')

curdoc().add_root(grid)
