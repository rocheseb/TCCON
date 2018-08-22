#!/usr/bin/env python2.7
 # -*- coding: utf-8 -*-

from __future__ import print_function # allows the use of Python 3.x print function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# Code Description #
####################
'''
# README #

Linefit is an instrument line shape-calculating program written by Frank Hase at the Institute for Meteorology und Climate Research (IMK) in Karlsruhe. 

	Frank Hase, Thomas Blumenstock, and Clare Paton-Walsh, Analysis of the Instrumental Line Shape of High-Resolution Fourier Transform IR Spectrometers with Gas Cell Measurements and New Retrieval Software, Applied Optics, Vol. 38, Issue 15, pp. 3417-3422 (1999)

This app can be used to run linefit 14.5 and display its outputs.

This app requires python 2.7.x (not tested with python 3.x) with bokeh installed.

	Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html

### Python ###

I suggest downloading python from https://www.anaconda.com/download/
Choose your operating system and get Python2.7

To install bokeh use the command (with windows you need to run the terminal as administator):

	conda install -c bokeh bokeh

To install a package:

	conda install PackageName

If the "conda" command does not find the package, you can also use pip:

	pip install PackageName

If you encounter error messages related to bokeh when running the app, you can try to revert to an earlier version of the package with:

	conda install bokeh=0.12.10

### How to use this app ###

This can read both OPUS and .dpt files.

Using .dpt files requires extra steps.

- Put the lft_app folder in the linefit/lft145/ directory
- Spectrum file names need to follow this naming convention: YYMMDD_site_cell_MOPD_X_num.ext
	- YYMMDD year month day
	- site: two letter site abbreviation (e.g. eu for Eureka, oc for Lamont)
	- cell: one of 'hbr', 'n2o', 'hcl'
	- X: 'v' for vented instrument, 'e' for evacuated
	- MOPD: the maximum optical path difference in cm
	- num: an index number for the cell test (there might be more than one per day)
	- ext: either 'dpt' or a number
	
		e.g. 180308_eu_HCl_45_e_0.dpt or 180308_eu_HCl_45_e_0.0 
		For the first HCl cell test with 45 MOPD in an evacuated instrument at Eureka on March 8 2018
		For opus files the extension number does not matter

- For several tests in one day : 161122_eu_HCl_45_e_0.dpt, 161122_eu_HCl_45_e_1.dpt etc.
- Spectra must be placed in lft_app/spectra/cut
- For '.DPT' files:
	- .dpt (data point table) files can be generated in OPUS via the pop-up window generated with "Save as" 
	- they must have no headers and be cut e.g. between ~5200-5900 wavenumbers for HCl cells
	- In lft_app/spectra/cut/temp.dat write the spectrum filename, scanner temperature, and aperture size.
	
			spectrumfilename1,temperature1,apt_size1
			spectrumfilename2,temperature2,apt_size2
			etc.

	- in lft_app/lft_setup.py add the focal length of collimator of your instrument
- For opus files the parameters are read directly from the file
	
- In lft_app/lft_setup.py, add your cell information (follow the template)

- To run the app, navigate to the linefit/lft145/ directory in your terminal and use the command

	bokeh serve --show lft_app --args light

The --show option will pop up the browser.

By default the spectrum itself will be plotted in the browser and also saved in the data dictionary. This can lead to very large files and more loading time.

The --args light option will not display or save the whole spectrum.

While the server is running, the app will be available in the browser at localhost:5006/lft_app

- Python dictionaries of the data are saved in lft_app/saved_sessions/
- PDF documents with the plots are saved in lft_app/pdf/

There are two example spectra from Eureka in lft_app/spectra/cut/ so the app can be tested directly.

### Rationg of spectra ###

Spectra should be ratioed to ~1 to be used with the linefit extended mode:

- HCl cells: 
	- no background
	- In the code I fit a 2nd order polynomial to the spectrum without the lines and use that to ratio the spectrum to normalize it to ~1 (seems more consistent than using a fixed numbers)

- HBr cells:
	- background
	- the background file should be cut the same way as the spectrum if .dpt files are used, have the same file name but starting with 'ref_' (e.g. ref_180308_eu_HBr_180_e_0.dpt)
	- put the HBr background files in lft_app/spectra/background/
	- the spectra are ratioed with the background
	- the resulting ratioed spectrum is ratioed with its own average to normalize it to ~1

- N2O cells:
	- background, but different resolution from the spectrum
	- the rationg of spectrum with background shoul be done in OPUS
	- the resulting spectrum should be placed in lft_app/spectra/cut
	- it will be ratioed with its own average to normalize it to ~ 1

### Other info ###

N2O and HBr cell spectra are processed in a loop until the cell pressure converges; this usually take 2-3 linefit runs.

The python dictionaries saved in lft_app/saved_sessions/ can be merged with a utility program lft_app/utils/merge_sessions.py

The merged file can then be loaded from the browser.

### DISCLAIMER ###

If any warning or error message is given by linefit, this app will hang, you should then run linefit from the terminal to figure out what the problem is.

The app may hang if there is any convergence problem.

There will be more detailed outputs in the terminal than in the browser.

If a significant spectral detuning is detected, the app will use it to update the input file and re-run linefit.

### Contact ###

sebastien.roche@mail.utoronto.ca

'''
####################
# Import libraries #
####################

# class to read opus files
from opus import *

# generic libraries
import os
import sys
import platform
import subprocess
from functools import partial

import collections
from collections import OrderedDict

import time
from datetime import datetime, timedelta

# interactive plotting library
import bokeh
from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, CustomJS, Button, Div, TextInput, Select, Panel, Tabs, Legend, DataRange1d, RadioButtonGroup, Legend
from bokeh.layouts import gridplot, widgetbox, Column, Row
from bokeh.palettes import viridis

# special arrays with special functions for easier vectorized operations
import numpy as np

import lft_setup # where the site instrument and cell information is kept

# for generating the pdf report with plots
import pylab as pl
import matplotlib.backends.backend_pdf as mpl_pdf

# string parsing libraries
import parse
import re

# template for modifying input file
from jinja2 import Template

# color palettes for plots
from bokeh.palettes import Viridis256,viridis

#####################
# General Functions #
#####################

def correct_Pressure(P,T):
	Tin = float(T)
	Pref = float(P)
	return Pref*Tin/296.0

def compute_pressure(col,T,l=0.02):
	R = 8.314	# gas constant
	Na = 6.02e23	# avogadro number

	return 0.01 * col * R * T / (Na * l) # 0.01 to convert from Pa to hPa

#cell column = 7.243e24 * p[mbar] * l[m] / T[K])
# I have noticed that using that formula often gives worst column scale factors than not using it
def correct_column(P,T,l=0.1):
	Tin = float(T)
	Pin = float(P)
	return 7.243e24*Pin*l/Tin

def linediv(color='lightblue',width=400):
	"""
	Function to generate a Div widget with a straight line.
	This is because a same model cannot be used several times in a same document.
	"""
	return Div(text='<hr width="100%" color="{}">'.format(color),width=width)

#########
# Setup #
#########

argu = sys.argv
ignore_spec = False
if 'light' in argu:
	# bokeh serve lftp_app --args light 
	ignore_spec = True # the spectrum will not be plotted or saved in the numpy file

system = platform.system()
if system == 'Windows':
	lft_command = ['lft145.exe']
elif system == 'Linux':
	lft_command = ['./lft145_ifort']
elif system == 'Darwin':
	lft_command = ['./lft145_gfortran']

specname_fmt = '{}_{}_{}_{:d}_{}_{:d}.{}' # formatted string for spectrum file names YYMMDD_site_cell_X_MOPD_num

fmt = '{:.2f},.false.,{:.4e},{:.3f},.false.,{:.3f},0.0075\n' # formatted string to read and edit the input file

reg_dpt = re.compile('.*[.]dpt$',re.IGNORECASE) # regular expression that will be used to select .dpt files
reg_opus = re.compile('.*[.][\d]+$',re.IGNORECASE) # regular expression that will be used to select opus files
reg_npy = re.compile('.*[.]npy$',re.IGNORECASE) # regular expression that will be used to select .npy files

site_data = lft_setup.site_data()

cell_map = {'hcl':site_data['hcl'],'hbr':site_data['hbr'],'n2o':site_data['n2o']} # read cell data

# This is a color list of 20 "high contrast" colors.
kelly_colors = [ 	'#F3C300','#875692', '#F38400', '#A1CAF1','#BE0032', '#C2B280', '#848482','#008856', '#E68FAC', '#0067A5',
				 	'#F99379', '#604E97', '#F6A600','#B3446C', '#DCD300', '#882D17','#8DB600', '#654522', '#E25822','#2B3D26',		]

# You can add any number of colors to that list, for example do something like below
#kelly_colors += ['blue','red','green','chartreuse','magenta','firebrick']
# but it is not really possible to have only high contrast colors for more than 20 colors

# if more than 20 lines are plotted the color palette will be switched to viridis(n) with n the number of lines
# you can check the viridis palette here: https://bokeh.pydata.org/en/latest/docs/reference/palettes.html
# above 20 lines, each time a line is added, all line colors will change because viridis(n) does a linear mapping with the 256 colors of Viridis256
# this will maximize contrast within the color palette

window_dict = {}

# MIR microwindows for HCl
window_dict['hcl_test'] = [
				'(2843.1,2844.1)', # 1
				#'(2596.8,2597.8)',
				#'(2625.2,2626.2)',
				'(2862.5,2863.5)', # 2
				'(2864.6,2865.6)', # 3
				'(2903.6,2904.6)', # 4
				'(2905.7,2906.7)', # 5
				'(2923.2,2924.2)', # 6
				'(2923.2,2924.2)', # 7
				'(2923.2,2924.2)', # 8
				'(2925.4,2926.4)', # 9
				'(2960.6,2961.6)', # 10
				#'(2962.8,2963.8)', # 11
				'(3011.6,3012.6)', # 12
				'(3013.9,3014.9)', # 13
				'(3027.3,3028.3)', # 14
				#'(3029.6,3030.6)', # 15
				]

# NIR microwindows for HCL
window_dict['hcl'] = [	
				'(5683.0,5684.0)', # 1 
				'(5687.1,5688.1)', # 2
				'(5701.5,5702.5)', # 3
				'(5705.6,5706.6)', # 4
				'(5718.7,5719.7)', # 5 
				'(5734.6,5735.6)', # 6 
				'(5738.8,5739.8)', # 7
				'(5749.3,5750.3)', # 8
				'(5753.5,5754.5)', # 9
				'(5762.7,5763.7)', # 10
				'(5766.9,5767.9)', # 11
				'(5774.8,5775.8)', # 12
				'(5779.0,5780.0)', # 13
				]

# MIR microwindows for N2O
window_dict['n2o'] = [
				'(2167.03,2185.25)', # 1
				'(2222.825,2223.019)', # 2	
				'(2224.457,2224.715)', # 3
				]

# MIR microwindows for HBr
window_dict['hbr'] = [
				'(2590.32, 2590.72)', # 1
				'(2590.71, 2591.11)', # 2
				'(2605.60, 2606.00)', # 3
				'(2606.00, 2606.40)', # 4
				'(2620.39, 2620.79)', # 5
				'(2620.80, 2621.20)', # 6
				'(2634.70, 2635.10)', # 7
				'(2635.10, 2635.50)', # 8
				'(2648.50, 2648.90)', # 9
				'(2648.90, 2649.30)', # 10
				'(2661.76, 2662.16)', # 11
				'(2662.18, 2662.58)', # 12
				'(2674.52, 2674.92)', # 13
				'(2674.94, 2675.34)', # 14
				]

app_path = os.path.dirname(__file__) # the app should be in ... /linefit/lft145/lft_app
spec_path = os.path.join(app_path,'spectra','cut') # ... /linefit/lft145/lft_app/spectra/cut
ref_path = os.path.join(app_path,'spectra','background') # ... /linefit/lft145/lft_app/spectra/background
save_path = os.path.join(app_path,'saved_sessions') # ... /linefit/lft145/lft_app/saved_sessions
static_path = os.path.join(app_path,'static') # ... /linefit/lft145/lft_app/static

wdir = os.sep.join(app_path.split(os.sep)[:-1]) # get the working directory; path to ... /linefit/lft145
erg_path = os.path.join(wdir,'ergs') # ... /linefit/lft145/ergs

print('app_path:',app_path)
print('spec_path:',spec_path)
print('ref_path:',ref_path)
print('save_path:',save_path)
print('static_path:',static_path)
print('wdir:',wdir)
print('erg_path:',erg_path)

for path in [app_path,spec_path,ref_path,save_path,static_path,erg_path]:
	if not os.path.exists(path):
		print('WARNING: Missing path',path)

spec_input_code = """
if (cb_obj.value===""){
	status_div.text = 'Select a Spectrum';
} else {
	status_div.text='Use the button to run linefit with the selected spectrum';
}
"""

TOOLS = "box_zoom,wheel_zoom,pan,redo,undo,reset,save" # tools that will be available to interact with the plots

all_data = {'ID':0} # this dictionary will store all the data for plots; 'ID' will store the ID of the appriopriate color from 'kelly_colors'

# Input file templates:
template_map = {}
for cell in ['hcl','n2o','hbr']:
	with open(os.path.join(app_path,'lft14_{}_template.inp'.format(cell)),'r') as infile:
		template_map[cell] = Template(infile.read())

template_inputs = {
	'maxopd':None,			# maximum optical path difference (cm)
	'maxir':None,			# maximum inclination of rays in interferometer (apterture radius)/(focal length of colimator)
	'temperature':None,		# scanner temperature (K)
	'column':None,			# column density of first gas
	'pressure':None,		# pressure of first gas
	'regularisation':None,	# regularisation factor (for both modulation and phase)
	'spectrum': None,		# full path to the spectrum
	'window_list':None,		# list of tuples for microwindow boundaries
}

##################
# Main functions #
##################

def execute(cmd,cwd=os.getcwd()):
	'''
	Function to execute a prompt command and print the output
	'''
	popen = subprocess.Popen(cmd,stdout=subprocess.PIPE,stdin=subprocess.PIPE,cwd=cwd)
	output,err = popen.communicate(input="0\n")
	return_code = popen.wait()
	if return_code:
		raise subprocess.CalledProcessError(return_code, cmd)
	print(output)

	output = output.splitlines()
	if "stop program" in output[-1]:
		spectral_detuning = float(output[-2])
		status_div = curdoc().select_one({"name":"status_div"})
		
		status_div.text+='<br>- Spectral detuning detected: {:9.2e}'.format(spectral_detuning)
		status_div.text+='<br>- Rerun linefit with spectral detuning corrected'

		input_file = os.path.join(wdir,'lft14.inp')
		with open(input_file,'r') as infile:
			content = infile.readlines()

		for i,line in enumerate(content):
			if ".true.,1.0" in line:
				content[i] = ".true.,1.0,{:.2E}\n".format(spectral_detuning)

		with open(input_file,'w') as outfile:
			outfile.writelines(content)

		execute(cmd,cwd=cwd)

def busy(func):
	'''
	Decorator function to display a loading animation when the program is working on something
	'''
	def wrapper():
		# show the loading animation for time consuming tasks
		with open(os.path.join(static_path,'loader.html'),'r') as infile:
			loader_css = infile.read().replace('\n',"")
		curdoc().select_one({'name':'loader'}).text = loader_css

		func() # run the decorated function

		# hide the loading animation
		curdoc().select_one({'name':'loader'}).text = ""
	return wrapper

def get_inputs(spectrum,mode):
	'''
	retrieves info from the spectrum name
	
	Inputs:
		- spectrum: spectrum file name
		- mode: 'dpt' or 'opus'
	Outputs:
		- site: two letter site abbreviation
		- cell: 'hcl', 'n2o', or 'hbr'
		- MOPD: maximum optical path difference (cm)
		- APT: aperture size (mm)
		- temperature: scanner temperature (K)
		- window_list: list of microwindows that corresponds to the cell
	'''

	date,site,cell,MOPD,ev,num,ext = parse.parse(specname_fmt,spectrum)

	cell = cell.lower()
	site = site.lower()
	ev = ev.lower()

	window_list = window_dict[cell]

	for site in cell_map[cell]: # check the cell data for the appriopriate site
		if '_{}_'.format(site) in spectrum:
			curdoc().select_one({"name":"status_div"}).text += '<br>- {} {} cell specs found'.format(cell_map[cell][site]['location'],cell)
			break
	else: # if the site is not specified, stop the program)
		print("\nSite not recognized: filename must be YYMMDD_site_cell_OPD_ev_num.dpt , with 'site' the two letter site abbreviation")
		print("\ne.g.170315_eu_HCl_45_e_0.dpt")
		print("\nCheck that the cell information for the site is entered correctly in lft_setup.py")
		curdoc().select_one({"name":"status_div"}).text += '<br>- Site not recognized'
		sys.exit()
	print(site)

	if mode == 'dpt':
		#get the temperature and aperture size
		with open(os.path.join(spec_path,'temp.dat'),'r') as infile:
			content = infile.read().splitlines()

		for line in content:
			if spectrum.lower() in line.lower():
				break

		line = line.split(',')
		temperature = float(line[1]) # scanner temperature in Kelvin
		APT = float(line[2]) # aperture diameter in millimeters
	elif mode == 'opus':
		opus_file = Opus(os.path.join(spec_path,spectrum))
		opus_file.get_data(request='p') # only get parameters
		parameters = opus_file.param[0]
		temperature = parameters['TSC']+273.15
		APT = float(parameters['APT'].split()[0])
		site_data['FOC'][site] = parameters['FOC'] # if it existed, the value from lft_setup will be overwritten

	return site,cell,str(MOPD),APT,temperature,window_list

def modify_input_file(spectrum,site,cell,MOPD,APT,temperature,window_list):
	'''
	Update the linefit input file to correspond to the selected spectrum and regularisation factor.
	For each site, the cell information must be added to the cell_data.py file

	Spectrum file names need to follow this naming convention: YYMMDD_site_cell_X_MOPD_num.dpt
	YYMMDD year month day
	site: two letter site abbreviation (e.g. eu for Eureka, oc for Lamont)
	cell: one of 'hbr', 'n2o', 'hcl'
	X: 'v' for vented instrument, 'e' for evacuated
	MOPD: the maximum optical path difference in cm
	num: an index number for the cell test (there might be more than one per day)

	e.g. 180308_eu_HCl_45_e_0
	for the first HCl cell test with 45 MOPD in an evacuated instrument at Eureka on March 8 2018
	'''

	N_windows = len(window_list)

	reg = curdoc().select_one({'name':'reg_input'}).value

	maxir = '{:>8.6f}'.format(APT/2.0/site_data['FOC'][site]) # maximum inclination of rays in the interferometer (aperture radius / focal length of collimator)

	template_inputs.update({
		'maxopd':MOPD,											# maximum optical path difference (cm)
		'maxir':maxir,											# maximum inclination of rays in interferometer (apterture radius)/(focal length of colimator)
		'window_list':window_list,								# list of tuples for microwindow boundaries
		'N_windows':N_windows,
		'MW_list':['MW'+str(i+1) for i in range(N_windows)],		
		'regularisation':reg,									# regularisation factor (for both modulation and phase)
		'spectrum': os.path.join('lft_app','spectra',spectrum),	# full path to the spectrum
		'temperature':'{:.2f}'.format(temperature),				# scanner temperature (K)
		})

	if cell == 'hcl':
		# HCl35
		newP = correct_Pressure(cell_map[cell][site]['effp_h35cl_296k'],temperature) # update the cell pressure
		newcol = correct_column(newP,temperature,l=0.1) # update the cell column as described in the input file
		print(newP,newcol)
		#newP = cell_map[cell][site]['effp_h35cl_296k']  # uncomment to use the effective pressure from the TCCON wiki
		newcol = cell_map[cell][site]['h35cl_column'] # uncomment to use the column from the TCCON wiki
		print(newP,newcol)

		template_inputs.update({
				'column':'{:.4e}'.format(newcol),						# column density of first gas
				'pressure':'{:.3f}'.format(newP),						# pressure of first gas
			})

		# HCl37
		newP = correct_Pressure(cell_map[cell][site]['effp_h37cl_296k'],temperature) # update the cell pressure
		newcol = correct_column(newP,temperature,l=0.1) # update the cell column as described in the input file
		print(newP,newcol)
		#newP = cell_map[cell][site]['effp_h37cl_296k'] # uncomment to use the effective pressure from the TCCON wiki
		newcol = cell_map[cell][site]['h37cl_column'] # uncomment to use the column from the TCCON wiki
		print(newP,newcol)

		template_inputs.update({
			'column_2':'{:.4e}'.format(newcol),	# column density of second gas
			'pressure_2':'{:.3f}'.format(newP),	# pressure of second gas
			})

	if cell in ['hbr','n2o']:
		#newP = correct_Pressure(cell_map[cell][site]['pressure'],temperature) # update the cell pressure
		#newcol = correct_column(newP,temperature,l=0.02) # update the cell column as described in the input file
		newP = cell_map[cell][site]['pressure'] # uncomment to use the initial cell pressure
		newcol = cell_map[cell][site]['column'] # uncomment to use the initial cell column

		template_inputs.update({
			'column':'{:.4e}'.format(newcol),
			'pressure':'{:.3f}'.format(newP),
			})

	with open(os.path.join(wdir,'lft14.inp'),'w') as outfile: #rewrite input file
		outfile.writelines(template_map[cell].render(**template_inputs).replace("\\n","\n"))

	curdoc().select_one({"name":"status_div"}).text+='<br>- Input file updated'
	print('\n\t- Input file updated')

@busy
def setup_linefit():
	'''
	setup a linefit run for the spectrum selected in spec_input
	'''

	global all_data

	status_div = curdoc().select_one({"name":"status_div"})

	dum_fig = curdoc().select_one({"name":"dum_fig"})

	if len(dum_fig.legend[0].items)>0:
		all_data['ID'] += 1

	curdoc().select_one({"name":"spec_input"}).js_on_change('value', CustomJS(args={'status_div':curdoc().select_one({"name":"status_div"})},code=spec_input_code))

	spectrum = curdoc().select_one({"name":"spec_input"}).value

	if reg_opus.match(spectrum):
		mode = 'opus'
	elif reg_dpt.match(spectrum):
		mode = 'dpt'

	if spectrum=='':
		status_div.text = "Select a spectrum"
		return

	reg = curdoc().select_one({"name":"reg_input"}).value

	try:
		float(reg)
	except:
		status_div.text = "Regularisation factor must be a number"
		return

	dum_leg_labels = [elem.label['value'] for elem in dum_fig.legend[0].items]
	already_done = [elem for elem in dum_leg_labels if ((spectrum.split('.')[0] in elem) and ('reg={}'.format(reg) in elem))]!=[]
	if already_done:
		status_div.text = "{} already analysed with reg={}".format(spectrum,reg)
		return
	
	status_div.text = "<b>Now doing:</b> <br>{}<br>reg= {}".format(spectrum,reg)
	print('\nNow doing',spectrum,'with reg=',reg)

	if mode == 'dpt':
		# preliminary check on the temp file to make sure it has the spectrum
		# I write my own inputfile called 'temp.dat' that has lines with 'SpectrumName,Scannertemperature,ApertureSize'
		with open(os.path.join(spec_path,'temp.dat'),'r') as infile:
			content = infile.readlines()

		speclist = [line.split(',')[0] for line in content[1:]] # all the SpectrumName in the file

		#if the spectrum is not listed in the temp file, go to next spectrum
		if spectrum.lower() not in [spec.lower() for spec in speclist]:
			status_div.text = spectrum+':</br>scanner temperature not listed in the temp file'
			print(spectrum,'scanner temperature not listed in the temp file')
			all_data['ID'] += -1
			return

	site,cell,MOPD,APT,temperature,window_list = get_inputs(spectrum,mode)
	
	colo = check_colors(add_one=True)
	
	spectrum_path = os.path.join(spec_path,spectrum)
	
	if mode == 'dpt':
		# check that the spectral range is ordered (dpt files are written with decreasing wavenumbers and linefit wants increasing wavenumbers)
		# if it is not ordered, orders it.
		
		check_spectrum(spectrum_path,spectrum)
		# also check the background spectrum for MIR cells
		if cell == 'hbr':
			ref_path = os.path.join(ref_path,'ref_'+spectrum)
			check_spectrum(ref_path,'ref_'+spectrum)

	# comment out to not ratio the spectrum, the temp file still needs to be in lft_app/spectra/cut and the spectrum will need to be directly in lft_app/spectra
	ratio_spectrum(spectrum_path,spectrum,cell,mode) 

	# update the input file; make sure that it modifies everything that you need !
	# the regularisation factors are updated from the browser
	modify_input_file(spectrum,site,cell,MOPD,APT,temperature,window_list)

	run_linefit(cell)

	# store results in all_data and update plots
	linefit_results(spectrum,colo)

	curdoc().select_one({'name':'status_div'}).text += "<br><b>DONE</b>"

def check_colors(add_one=False):
	'''
	If there are more tests to be displayed than colors in the kelly_colors list, change the kelly_colors color palette for a linear mapping of Viridis256 colors
	'''

	global all_data

	add = 0 # when used in update_doc()
	if add_one:
		add = 1	# when used in setup_linefit()

	test_list = sorted([i for i in all_data.keys() if 'reg' in i])

	N_tests = len(test_list)

	if N_tests >= len(kelly_colors):
		print('\n\t - Too many lines for kelly_colors, using large palette')
		new_palette = viridis(N_tests)[::-1]	# linear mapping from 256 Viridis colors to the N differents tests
		for i,test in enumerate(test_list):
			all_data[test]['color'] = new_palette[i]
		if add_one:
			update_colors()
		return new_palette[-1]
	else:
		if add_one:
			return kelly_colors[all_data['ID']]
		else:
			for i,test in enumerate(test_list):
				all_data[test]['color'] = kelly_colors[i]			

def run_linefit(cell):
	'''
	Run linefit and printst he output in the terminal as it is running.

	Once for HCl cells
	In a loop for HBr and N2O cells, until the pressure converges
	'''

	status_div = curdoc().select_one({"name":"status_div"})

	status_div.text+='<br>- Running linefit ...'
	print('\n\t- Running linefit ...')
	print('\t executing command {} in {}'.format(lft_command[0],wdir))
	execute(lft_command,cwd=wdir)

	if cell in ['n2o','hbr']:
		iteration = 1
		conv = False
		while not conv:
			# open the input file
			with open(os.path.join(wdir,'lft14.inp'),'r') as infile:
				content = infile.readlines()

			# read the column and pressure
			for i,line in enumerate(content):
				if 'species parameters:' in line:
					temperature,column,pressure,pressure = parse.parse(fmt,content[i+8])
					break

			# read the column scale factor
			with open(os.path.join(erg_path,'colparms.dat'),'r') as infile:
				col_content = infile.readlines()

			scale_factor = np.mean([float(x) for x in col_content[1:]])

			# compute the scaled column and the new pressure
			new_column = scale_factor * column
			new_pressure = compute_pressure(new_column,temperature)

			print('\n\t- pres,new_pres,dif:',pressure,new_pressure,abs(pressure-new_pressure))
			# check for convergence
			if abs(pressure-new_pressure)<0.001:
				conv = True
				break

			# force stop if a certain number of iterations is reached
			iteration +=1
			if iteration>5:
				break

			# replace the pressure in the input file
			content[i+8] = fmt.format(temperature,column,new_pressure,new_pressure)

			with open(os.path.join(wdir,'lft14.inp'),'w') as outfile:
				content = outfile.writelines(content)

			# run a new iteration of linefit
			status_div.text+='<br>- Running linefit iteration {}'.format(iteration)
			print('\n\t- Running linefit iteration',iteration)
			execute(lft_command,cwd=wdir)
				
		if conv:
			print('\n\t- convergence after',iteration,'iterations')
			status_div.text+='<br>- pressure converged'
		else:
			print('\n\t- no convergence')
			status_div.text+='<br>- pressure <b>did not<b> converge'

def check_spectrum(spectrum_path,spectrum):
	'''
	check that the spectral range is ordered (dpt files are written with decreasing wavenumbers and lienfit wants increasing wavenumbers)
	if it is not ordered, orders it.
	'''

	status_div = curdoc().select_one({"name":"status_div"})

	with open(spectrum_path,'r') as infile:
		content = infile.readlines()

	specrange = [line.split()[0] for line in content]

	if specrange!=sorted(specrange):
		status_div.text += '<br>- Reordering {}'.format(spectrum)
		print('\n\t- Reordering',spectrum)

		with open(spectrum_path,'w') as outfile:
			outfile.writelines(content[::-1])

def ratio_spectrum(spectrum_path,spectrum,cell,mode):
	'''
	For HCl cells, fit a second order polynomial to a spectrum in order to ratio it to ~1

	For N2O and HBr cells, use a background spectrum to do the ratio
	'''

	if mode == 'dpt':
		x,y = np.loadtxt(spectrum_path,unpack=True)	
	elif mode == 'opus':
		opus_file = Opus(spectrum_path)
		opus_file.get_data()
		x = opus_file.xdata[1] # need to check if this always corresponds to the correct data block
		y = opus_file.ydata[1]
		# cut the spectrum
		if cell == 'hcl':
			minwn,maxwn = (5200,5900)
		elif cell == 'n2o':
			minwn,maxwn = (2100,2300)
		elif cell == 'hbr':
			minwn,maxwn = (2500,2700)

		y = y[(x>=minwn) & (x<=maxwn)]
		x = x[(x>=minwn) & (x<=maxwn)]
	
	if cell == 'hcl':
		if 0.7<np.mean(y)<1.3: # spectrum was already ratioed
			return

		# resample_y_max takes the max y values for each interval of 14 wavenumbers, it avoids getting points in the lines themselves
		resample_y_max = [max(y[(x<(start+15)) & (x>start)]) for start in range(int(x[0]),int(x[-1]),14)]
		resample_x = np.array(range(int(x[0]),int(x[-1]),14))

		fit = np.polyfit(resample_x,resample_y_max,2) # second order polynomial fit

		# I substract the 0.0004 because the way I fit will always be a bit too high as it will include values above the base line
		base_y = fit[0]*x**2+fit[1]*x+fit[2]-0.0004 # use the fit to get the fit line for each x

	elif cell == 'hbr':
		ref_path = os.path.join(ref_path,'ref_'+spectrum)
		if mode == 'dpt':	
			xref,yref = np.loadtxt(ref_path,unpack=True)
		elif mode == 'opus':
			ref_opus_file = Opus(ref_path)
			ref_opus_file.get_data()
			xref = ref_opus_file.xdata[1]
			yref = ref_opus_file.ydata[1]
			yref = yref[(xref>=minwn) & (xref<=maxwn)]
			xref = xref[(xref>=minwn) & (xref<=maxwn)]
			
		if len(x)!=len(xref):
			print('WARNING: background and cell spectra have different number of points !')
			sys.exit() # the correct cutting of the spectra should be handled outside this code.

		# ratio the cell spectrum and background spectrum
		y = y/yref
		base_y = np.mean(y)

	elif cell == 'n2o': # n2o spectra should be ratioed with their background in OPUS, then the resulting spectrum is saved in a .dpt file
		# just take the average intensity to do the ratio
		base_y = np.mean(y)

	new_y = y/base_y # ratio to ~1

	np.savetxt(os.path.join('lft_app','spectra',spectrum),np.transpose([x,new_y]),fmt='%10.5f\t%.5f') # write the ratioed spectrum in lft_app/spectra

	curdoc().select_one({"name":"status_div"}).text += '<br>- Spectrum ratioed to ~{:.4f}'.format(np.mean(new_y))


def linefit_results(spectrum,colo):
	'''
	After linefit has finished running, this stores the results in all_data and updates the plots
	'''

	global all_data, ignore_spec

	# string containing info on the current spectrum + current regularisation factor
	test = '{} reg={}'.format(spectrum.split('.')[0],curdoc().select_one({'name':'reg_input'}).value)

	# Add a new entry in the all_data dictionary
	all_data[test] = {'statevec':0,'ILS':0,'AK':0,'mw1':0,'mw2':0,'mw3':0,'mw4':0,'mw5':0,'mw6':0,'mw7':0,'mw8':0,'mw9':0,'mw10':0,'mw11':0,'mw12':0,'mw13':0,}

	# Update the title in the ils_fits_panel and ak_panel
	curdoc().select_one({"name":"cur_spec_div"}).text = "<font size=3 color='teal'><b>{}</b></font>".format(test)
	curdoc().select_one({"name":"cur_spec_div2"}).text = "<font size=3 color='teal'><b>{}</b></font>".format(test)
	# Update status
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding ME and PE plot'
	
	###############################
	############################### Modulation Efficiency and Phase Error
	print('\n\t- ME and PE')

	with open(os.path.join(erg_path,'modulat.dat'),'r') as infile:
		content = infile.readlines()

	content = [line.split() for line in content[1:]]
	content = np.array([[float(elem) for elem in row] for row in content]).T

	all_data[test]['color'] = colo
	all_data[test]['ME'] = {'x':content[0],'y':content[1]}
	all_data[test]['PE'] = {'x':content[0],'y':content[2]}

	cur_date_string = test[:6]
	cur_date_time_struct = time.strptime(cur_date_string,'%y%m%d')
	cur_date = datetime(*cur_date_time_struct[:6])
	all_data[test]['series'] = {'x':[cur_date],'y':[content[1][-1]],'name':[test]}

	ME_source = ColumnDataSource(data=all_data[test]['ME'])
	PE_source = ColumnDataSource(data=all_data[test]['PE'])
	series_source = ColumnDataSource(data=all_data[test]['series'])

	ME_line = curdoc().select_one({"name":"ME_fig"}).line(x='x',y='y',color=colo,line_width=2,source=ME_source,name='{} ME line'.format(test),legend=test)
	ME_line.on_change('visible',partial(show_hide,test=test))
	update_legend(test)
	curdoc().select_one({"name":"PE_fig"}).line(x='x',y='y',color=colo,line_width=2,source=PE_source,name='{} PE line'.format(test))
	curdoc().select_one({"name":"series_fig"}).scatter(x='x',y='y',color=colo,size=5,source=series_source,name='{} series scatter'.format(test))

	###############################
	############################### Column
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding column plot'
	print('\n\t- column')

	with open(os.path.join(erg_path,'colparms.dat'),'r') as infile:
		content = infile.readlines()

	if 'hcl' in spectrum.lower():
		secondID = [i for i,v in enumerate(content) if 'Species:   2' in v][0] #hcl37, I don't plot it though.
		second_spec = [float(i) for i in content[secondID+1:]]
		first_spec = [float(i) for i in content[1:secondID]]
	else:
		first_spec = [float(i) for i in content[1:]]

	mw = range(1,len(first_spec)+1) # just the microwindow numbers

	all_data[test]['COL'] = {'x':mw,'y':first_spec}

	COL_source = ColumnDataSource(data=all_data[test]['COL'])

	curdoc().select_one({"name":"column_fig"}).scatter(x='x',y='y',color=colo,source=COL_source,name='{} column scatter'.format(test))
	curdoc().select_one({"name":"column_fig"}).line(x='x',y='y',color=colo,line_width=2,source=COL_source,name='{} column line'.format(test))
	
	###############################
	############################### ILS
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding ILS'
	print('\n\t- Adding ILS')

	with open(os.path.join(erg_path,'ilsre.dat'),'r') as infile:
		content = infile.readlines()

	content = np.array([[float(elem) for elem in line.split()] for line in content]).T
 
	all_data[test]['ILS'] = {'x':content[0],'y':content[1]}

	curdoc().select_one({"name":"ILS_line"}).data_source.data.update(all_data[test]['ILS'])

	###############################
	if not ignore_spec:
		############################### Spectrum
		curdoc().select_one({"name":"status_div"}).text+='<br>- Adding spectrum'
		print('\n\t- Adding spectrum')

		x,y = np.loadtxt(os.path.join('lft_app','spectra','cut',spectrum),unpack=True)
	 
		all_data[test]['spec'] = {'x':x,'y':y}

		curdoc().select_one({"name":"spec_line"}).data_source.data.update(all_data[test]['spec'])
		spec_fig = curdoc().select_one({'name':'spec_fig'})
		spec_fig.x_range.start = np.min(all_data[test]['spec']['x'])
		spec_fig.x_range.end = np.max(all_data[test]['spec']['x'])
		
		###############################
	############################### Microwindows
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding spectral fits with residuals'
	print('\n\t- Adding spectral fits with residuals')

	mwfiles = [i for i in os.listdir(erg_path) if 'specre' in i]

	for it,mwfile in enumerate(mwfiles):
		with open(os.path.join(erg_path,mwfile),'r') as infile:
			content = infile.readlines()

		content = np.array([[float(elem) for elem in line.split()] for line in content]).T

		resid = 100*(content[1]-content[2])/content[2] # (measured-calculated)/calculated

		all_data[test]['mw{}'.format(it+1)] = {'x':content[0],'meas':content[1],'calc':content[2],'resid':resid}
		all_data[test]['rms_resid_mw{}'.format(it+1)] = "{:.5f}".format(np.sqrt(np.mean(resid**2)))

	mw_fig = curdoc().select_one({"name":"mw_fig"})
	mw_fig.title.text = "Microwindow 1"

	resid_fig = curdoc().select_one({"name":"resid_fig"})
	resid_fig.title.text = 'RMS = '+all_data[test]['rms_resid_mw1']

	curdoc().select_one({"name":"meas_line"}).data_source.data.update(all_data[test]['mw1'])
	curdoc().select_one({"name":"calc_line"}).data_source.data.update(all_data[test]['mw1'])
	curdoc().select_one({"name":"resid_line"}).data_source.data.update(all_data[test]['mw1'])

	# setup mw_buttons
	MW_buttons = curdoc().select_one({'name':'MW_buttons'}).labels = ['MW {}'.format(i+1) for i in range(len(mwfiles))]

	###############################
	############################### Averaging Kernels
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding averaging kernels'
	print('\n\t- Adding averaging kernels')

	AKapo_fig = curdoc().select_one({"name":"AKapo_fig"})
	AKphase_fig = curdoc().select_one({"name":"AKphase_fig"})

	with open(os.path.join(erg_path,'actparms.dat'),'r') as infile:
		content = infile.read().splitlines()

	sumID = 0
	all_data[test]['statevec'] = OrderedDict([('iterations',{})])
	for line in content:
		line = line.strip().split()
		try:
			it = int(line[0])
		except ValueError:
			try:
				all_data[test]['statevec'][line[0]] = int(line[1])
			except ValueError:
				all_data[test]['statevec'][" ".join(line[:-1])] = int(line[-1])
			if 'apo' in line:
				AKapoID = sumID
			elif 'phase' in line:
				AKphaseID = sumID
			sumID += int(line[-1])			
		else:
			all_data[test]['statevec']['iterations'][it] = np.array(line[1:]).astype(np.float)

	with open(os.path.join(erg_path,'kernel.dat'),'r') as infile:
		content = infile.read().splitlines()

	content = np.array([line.split() for line in content]).astype(np.float)

	points = AKphaseID-AKapoID

	all_data[test]['AK'] = {
							'x':[all_data[test]['ME']['x'] for i in range(points)],
							'AKapo':content[AKapoID:AKapoID+points,AKapoID:AKapoID+points].tolist(),
							'AKphase':content[AKphaseID:AKphaseID+points,AKphaseID:AKphaseID+points].tolist(),
							'opd':all_data[test]['ME']['x'],
							'fix':[-0.7 for i in range(points)],
							'color':viridis(points),
							}

	for elem in ["AKapo_line","AKphase_line","AKapo_scatter","AKphase_scatter"]:
		curdoc().select_one({"name":elem}).data_source.data.update(all_data[test]['AK'])

	AKapo_fig.title.text = 'AK apo; DOFS = {:5.3f}'.format(np.sum(np.diag(all_data[test]['AK']['AKapo'])))
	AKphase_fig.title.text = 'AK phase; DOFS = {:5.3f}'.format(np.sum(np.diag(all_data[test]['AK']['AKphase'])))

	###############################
	###############################

	add_button(test)

	print(spectrum,'DONE')

def add_button(test):
	'''
	add a button to the 'button_box' corresponding to a spectrum + regularisation factor
	Used at the end of the linefit_results function and in the update_doc function
	'''
	global all_data

	button = Button(label=test,width=180)
	remove_button = Button(label='X',width=30,tags=[test],css_classes=["remove_button"])

	button_box = curdoc().select_one({"name":"button_box"})
	button_box.children += [Row(children=[Column(children=[remove_button]),Column(children=[button])])]

	all_data['cur_clicks'] = [0 for i in range(len(button_box.children)-1)]
	all_data['prev_clicks'] = [0 for i in range(len(button_box.children)-1)]

	button.on_click(change_spectrum)
	remove_button.on_click(remove_test)

def remove_test():
	'''
	remove a test from all_data and remake the document
	'''

	global all_data

	button_box = curdoc().select_one({"name":"button_box"})

	for row in button_box.children[1:]:
		elem = row.children[0].children[0].children[0]
		if elem.clicks==1:
			break

	test = elem.tags[0]

	del all_data[test]

	doc_maker()

@busy
def update_doc():
	'''
	Use the current all_data to fill the plots in the document
	'''
	global all_data

	status_div = curdoc().select_one({"name":"status_div"})

	status_div.text = 'Updating plots ...'

	test_list = sorted([key for key in all_data if 'reg' in key])

	colo = check_colors()

	# add the buttons corresponding to the current all_data dictionary keys
	for test in test_list:

		add_button(test)

		colo = all_data[test]['color']

		# Fill data sources
		ME_source = ColumnDataSource(data=all_data[test]['ME'])
		PE_source = ColumnDataSource(data=all_data[test]['PE'])
		COL_source = ColumnDataSource(data=all_data[test]['COL'])
		series_source = ColumnDataSource(data=all_data[test]['series'])
		# Update lines
		ME_line = curdoc().select_one({"name":"ME_fig"}).line(x='x',y='y',color=colo,line_width=2,source=ME_source,name='{} ME line'.format(test),legend=test,tags=[])
		ME_line.on_change('visible',partial(show_hide,test=test))
		update_legend(test)
		curdoc().select_one({"name":"PE_fig"}).line(x='x',y='y',color=colo,line_width=2,source=PE_source,name='{} PE line'.format(test))
		curdoc().select_one({"name":"column_fig"}).line(x='x',y='y',color=colo,line_width=2,source=COL_source,name='{} column line'.format(test))
		curdoc().select_one({"name":"column_fig"}).scatter(x='x',y='y',color=colo,source=COL_source,name='{} column scatter'.format(test))
		curdoc().select_one({"name":"series_fig"}).scatter(x='x',y='y',color=colo,size=5,source=series_source,name='{} series scatter'.format(test))

	status_div.text = "Ready"

def show_hide(attr,old,new,test):
	'''
	link visible attributes of lines
	'''
	ME_line = curdoc().select_one({'name':'{} ME line'.format(test)})

	PE_line = curdoc().select_one({'name':'{} PE line'.format(test)})
	column_line = curdoc().select_one({'name':'{} column line'.format(test)})
	column_scatter = curdoc().select_one({'name':'{} column scatter'.format(test)})
	series_scatter = curdoc().select_one({'name':'{} series scatter'.format(test)})

	for rend in [PE_line,column_line,column_scatter,series_scatter]:
		rend.visible = ME_line.visible

def update_legend(test):
	'''
	When a new line is added, move it over to dum_fig and remove it from the ME_fig
	'''

	ME_fig = curdoc().select_one({'name':'ME_fig'})
	dum_fig = curdoc().select_one({'name':'dum_fig'})

	glyph = curdoc().select_one({'name':'{} ME line'.format(test)})

	legend = ME_fig.legend[0]

	dum_fig.legend[0].items += legend.items
	dum_fig.renderers += [glyph]
	
	legend.visible = False
	ME_fig.renderers = [i for i in ME_fig.renderers if type(i)!=bokeh.models.annotations.Legend]

def update_colors():
	'''
	Use info in all_data to update the line colors
	'''
	global all_data

	test_list = sorted([key for key in all_data if 'reg' in key])

	for test in test_list:

		colo = all_data[test]['color']

		curdoc().select_one({'name':'{} ME line'.format(test)}).glyph.line_color = colo
		curdoc().select_one({'name':'{} PE line'.format(test)}).glyph.line_color = colo
		curdoc().select_one({'name':'{} column line'.format(test)}).glyph.line_color = colo
		curdoc().select_one({'name':'{} column scatter'.format(test)}).glyph.fill_color = colo
		curdoc().select_one({'name':'{} series scatter'.format(test)}).glyph.fill_color = colo

def change_spectrum():
	'''
	callback for the spectrum buttons in 'button_box'
	update the plots in 'ils_fits_grid' to correspond to the desired spectrum
	'''
	global all_data, ignore_spec

	button_box = curdoc().select_one({"name":"button_box"})

	# list of current number of clicks for each spectrum button (including the one that just got clicked)
	all_data['cur_clicks'] = [row.children[1].children[0].children[0].clicks for row in button_box.children[1:]]

	# compare cur_clicks to prev_clicks to know which button has just been clicked
	for i in range(len(all_data['cur_clicks'])):
		if all_data['cur_clicks'][i]!=all_data['prev_clicks'][i]:
			test = button_box.children[1:][i].children[1].children[0].children[0].label
	# set prev_clicks equal to cur_clicks 
	all_data['prev_clicks'] = [i for i in all_data['cur_clicks']]

	curdoc().select_one({"name":"cur_spec_div"}).text = "<font size=3 color='teal'><b>{}</b></font>".format(test)
	curdoc().select_one({"name":"cur_spec_div2"}).text = "<font size=3 color='teal'><b>{}</b></font>".format(test)
	
	# Select microwindow and residuals figures
	mw_fig = curdoc().select_one({"name":"mw_fig"})
	resid_fig = curdoc().select_one({"name":"resid_fig"})
	# Get current microwindow
	cur_MW = mw_fig.title.text.split()[1]	
	# Update titles
	resid_fig.title.text = 'RMS = '+all_data[test]['rms_resid_mw'+cur_MW]
	# Get ILS, microwindow, and averaging kernel data for the new spectrum
	new_mw_data =  all_data[test]['mw'+cur_MW]
	new_ILS_data =  all_data[test]['ILS']
	new_ak_data = all_data[test]['AK']

	# Update lines
	curdoc().select_one({"name":"meas_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"calc_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"resid_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"ILS_line"}).data_source.data.update(new_ILS_data)
	curdoc().select_one({"name":"AKapo_line"}).data_source.data.update(new_ak_data)
	curdoc().select_one({"name":"AKphase_line"}).data_source.data.update(new_ak_data)

	# Update AK figure titles with DOFs
	AKapo_fig = curdoc().select_one({"name":"AKapo_fig"})
	AKphase_fig = curdoc().select_one({"name":"AKphase_fig"})
	AKapo_fig.title.text = 'AK apo; DOFS = {:5.3f}'.format(np.sum(np.diag(new_ak_data['AKapo'])))
	AKphase_fig.title.text = 'AK phase; DOFS = {:5.3f}'.format(np.sum(np.diag(new_ak_data['AKphase'])))

	if not ignore_spec:
		curdoc().select_one({"name":"spec_line"}).data_source.data.update(new_spec_data)
		# update spectrum range
		new_spec_data = all_data[test]['spec']
		spec_fig = curdoc().select_one({'name':'spec_fig'})
		spec_fig.x_range.start = np.min(all_data[test]['spec']['x'])
		spec_fig.x_range.end = np.max(all_data[test]['spec']['x'])

	cell = test.split('_')[2].lower()

	# update the microwindow buttons
	MW_buttons = curdoc().select_one({'name':'MW_buttons'})
	MW_buttons.labels = ['MW {}'.format(i+1) for i in range(len(window_dict[cell]))]

def change_microwindow(attr,old,new):
	'''
	callback for the microwindow buttons
	update the plots in 'mw_grid' so that they correspond to the desired microwindow
	'''
	global all_data

	new_MW = str(new+1)

	# Get current spectrum
	cur_spec_div = curdoc().select_one({"name":"cur_spec_div"})
	if 'Spectrum' in cur_spec_div.text:
		curdoc().select_one({"name":"status_div"}).text = "No spectrum selected"
		return	

	# Get the corresponding key for the all_data dictionary
	test = cur_spec_div.text.split('<b>')[1].split('</b>')[0]
	
	# Get the new microwindow data
	new_mw_data = all_data[test]['mw'+new_MW]

	# Select microwindow and residuals figures
	mw_fig = curdoc().select_one({"name":"mw_fig"})
	resid_fig = curdoc().select_one({"name":"resid_fig"})
	# Update titles
	mw_fig.title.text = "Microwindow "+new_MW
	resid_fig.title.text = 'RMS = '+all_data[test]['rms_resid_mw'+new_MW]
	# Update lines
	curdoc().select_one({"name":"meas_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"calc_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"resid_line"}).data_source.data.update(new_mw_data)

	# update the x axis range
	mw_fig.x_range.start = new_mw_data['x'][0]
	mw_fig.x_range.end = new_mw_data['x'][-1]


def pdf_report(save_name):
	'''
	save all the plots in a multipage pdf document
	''' 
	pdf = mpl_pdf.PdfPages(os.path.join(app_path,'pdf',save_name+'.pdf'))
	mepefig,ax = pl.subplots(4,1)
	mepefig.set_size_inches(10,8)
	mepelines = []
	mepefig_indiv,ax_indiv = pl.subplots(3,1)
	mepefig_indiv.set_size_inches(10,8)

	mw_fig, mw = pl.subplots(2,1)
	ilsfig = pl.figure()
	mw_fig.set_size_inches(8,7)
	ilsfig.set_size_inches(8,8)
	mindate = datetime(2050,1,1)
	maxdate = datetime(1900,1,1)
	for test in sorted(all_data.keys()):

		if 'reg' in test:

			cell = test.split('_')[2].lower()
			
			# ILS
			ILS_source = all_data[test]['ILS']		
			pl.figure(ilsfig.number)
			pl.plot(ILS_source['x'],ILS_source['y'])
			pl.xlabel('Wavenumber (cm-1)')
			pl.ylabel('Response')
			pl.title(test)
			pl.tight_layout()
			pdf.savefig(ilsfig,bbox_inches='tight')
			pl.clf()
			
			# MEPECOL
			ME_source = all_data[test]['ME']
			PE_source = all_data[test]['PE']
			col_source = all_data[test]['COL']
			series_source = all_data[test]['series']
			if mindate>min(series_source['x']):
				mindate = min(series_source['x'])
			if maxdate<max(series_source['x']):
				maxdate = max(series_source['x'])
			
			# Individual MEPECOL
			pl.figure(mepefig_indiv.number)
			ax_indiv[0].plot(ME_source['x'],ME_source['y'],color=all_data[test]['color'])
			ax_indiv[1].plot(PE_source['x'],PE_source['y'],color=all_data[test]['color'])
			ax_indiv[2].plot(col_source['x'],col_source['y'],color=all_data[test]['color'])
			ax_indiv[0].set_title(test)
			ax_indiv[0].set_ylabel("ME")
			ax_indiv[1].set_ylabel("Phase Error")
			ax_indiv[1].set_xlabel("OPD (cm)")
			ax_indiv[2].set_ylabel("Column sf")
			ax_indiv[2].set_xlabel("Microwindow number")
			pl.tight_layout()
			pdf.savefig(mepefig_indiv,bbox_inches='tight')
			for elem in ax_indiv:
				elem.lines.pop(0)
			pl.cla()
			
			# Summary MEPECOL
			newline = ax[0].plot(ME_source['x'],ME_source['y'],color=all_data[test]['color'])
			mepelines += [newline[0]]
			ax[1].plot(PE_source['x'],PE_source['y'],color=all_data[test]['color'])
			ax[2].plot(col_source['x'],col_source['y'],color=all_data[test]['color'])
			ax[3].scatter(series_source['x'],series_source['y'],color=all_data[test]['color'])
			
			# microwindows
			pl.figure(mw_fig.number)
			for i in range(1,len(window_dict[cell])+1):
				mw_source = all_data[test]['mw{}'.format(i)]
				meas = mw[0].plot(mw_source['x'],mw_source['meas'],color='blue',label='measured')
				calc = mw[0].plot(mw_source['x'],mw_source['calc'],color='red',label='calculated')
				mw[0].set_title(test+': Microwindow {}'.format(i))
				resid = mw[1].plot(mw_source['x'],mw_source['resid'],color='black')
				mw[0].set_ylabel('Transmission')
				mw[1].set_ylabel('% Residuals')
				mw[1].set_xlabel('Wavenumber (cm-1)')
				mw[1].set_title('RMS='+all_data[test]['rms_resid_mw{}'.format(i)])
				mw[0].legend(loc=4)
				mw[0].relim()
				mw[0].autoscale()
				mw[0].get_xaxis().get_major_formatter().set_scientific(False)
				mw[1].get_xaxis().get_major_formatter().set_scientific(False)
				mw[0].get_xaxis().get_major_formatter().set_useOffset(False)
				mw[1].get_xaxis().get_major_formatter().set_useOffset(False)
				pl.subplots_adjust(hspace=0.2)
				pdf.savefig(mw_fig,bbox_inches='tight')
				mw[0].lines.remove(meas[0])
				mw[0].lines.remove(calc[0])
				mw[1].lines.remove(resid[0])
				pl.cla()
			
	pl.figure(mepefig.number)
	for elem in ax:
		box = elem.get_position()
		elem.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	lgd = pl.legend(mepelines,[test for test in all_data if 'reg' in test],bbox_to_anchor=(1.02, 4),loc='center left', borderaxespad=0)
	ax[0].set_ylabel("ME")
	ax[1].set_ylabel("Phase Error")
	ax[1].set_xlabel("OPD (cm)")
	ax[2].set_ylabel("Column sf")
	ax[2].set_xlabel("Microwindow number")
	ax[3].set_ylabel("ME at MOPD")
	ax[3].set_xlabel("Date")
	ax[3].set_xlim(mindate-timedelta(days=1),maxdate+timedelta(days=1))
	pl.tight_layout()
	pdf.savefig(mepefig,bbox_inches='tight',bbox_extra_artist=[lgd])

	pl.close('all')

	pdf.close()

@busy
def save_session():
	'''
	save the current all_data
	'''
	global all_data

	status_div = curdoc().select_one({"name":"status_div"})
	session_input = curdoc().select_one({"name":"session_input"})
	save_input = curdoc().select_one({"name":"save_input"})

	status_div.text = "Saving current session ..."

	save_name = save_input.value

	status_div.text += "<br> -saving the source dictionary"
	np.save(os.path.join(save_path,save_name),all_data) # save the all_data dictionary in a .npy file

	status_div.text += "<br> -writting PDF report"
	pdf_report(save_name) # save a pdf document with all the plots

	print('\nCurrent session data saved in :','lft_app/saved_sessions/{}.npy'.format(save_name))
	print('\nPDF report saved in :','lft_app/pdf/{}.pdf'.format(save_name))

	# update the status_div
	status_div.text = 'Current session data saved in:<br><b>lft_app/saved_sessions/{}.npy</b><br>PDF report saved in:<b><br>lft_app/pdf/{}.pdf</b>'.format(save_name,save_name)
	
	# now that a new session has been saved, update the options of the session_input widget
	session_input.options = ['']+[i for i in os.listdir(save_path) if reg_npy.match(i)]

def load_session():
	'''
	load a previously saved all_data and remake the whole document
	'''

	global all_data, ignore_spec

	curdoc().select_one({"name":"status_div"}).text = "Loading new session ..."
	print('\n\t- Loading new session ...')

	saved_session = os.path.join(save_path,curdoc().select_one({"name":"session_input"}).value)
	
	all_data = np.load(saved_session).item()

	test_list = [i for i in all_data.keys() if 'reg' in i]
	if not ignore_spec:
		# check that all tests in the saved data have a stored spectrum.
		# if they don't all have a spectrum, switch to the ignore_spec mode.
		for test in test_list:
			try:
				dum = all_data[test]['spec']['x'][0]
			except KeyError:
				ignore_spec = True
				break
	if ignore_spec:
		# set the 'spec' of each test to None to free some memory.
		for test in test_list:
			all_data[test]['spec'] = None

	doc_maker() # rebuild the entire document using the new all_data dictionary

	curdoc().select_one({"name":"spec_input"}).value = ""

	curdoc().select_one({"name":"status_div"}).text = "New session loaded"
	print('\n\t- New session loaded')

def update_dropdowns():
	'''
	update the dropdown lists
	'''

	session_input = curdoc().select_one({"name":"session_input"})
	spec_input = curdoc().select_one({"name":"spec_input"})

	# update the dropdown of saved sessions
	session_input.options = ['']+[i for i in os.listdir(save_path) if reg_npy.match(i)]

	#update the dropdown of spectra
	spec_input.options = ['']+[i for i in os.listdir(spec_path) if reg_dpt.match(i) or reg_opus.match(i)]	

def doc_maker():
	'''
	make the whole document
	'''

	global all_data, ignore_spec

	curdoc().clear() # removes everything in the current document

	## WIDGETS
	# Inputs
	spec_input = Select(title='Spectrum:',options = ['']+[i for i in os.listdir(spec_path) if reg_dpt.match(i) or reg_opus.match(i)],width=150,css_classes=["spec_input"],name="spec_input")
	reg_input = TextInput(value='1.8',title='Regularisation factor:',width=150,css_classes=["small_input"],name="reg_input")
	session_input = Select(title='Previous sessions:',width=150,options=['']+[i for i in os.listdir(save_path) if reg_npy.match(i)],css_classes=["spec_input"],name="session_input")
	save_input = TextInput(title='Save name',value="_".join(str(datetime.now())[:-7].split()).replace(':','-'),css_classes=["save_input"],name="save_input")
	loop_input = TextInput(title='Loop key',value="HCl_45",width=100,css_classes=["small_input"],name="loop_input")
	# BUTTONS
	lft_button = Button(label='Run linefit', width=80, css_classes=["custom_button"],name="lft_button")
	save_button = Button(label='Save Session', width=90, css_classes=["custom_button"],name="save_button")
	load_button = Button(label='Load Session', width=90, css_classes=["custom_button"],name="load_button")
	loop_button = Button(label='loop',width=90, css_classes=["custom_button"],name="loop_button")
	refresh_button = Button(label='Update dropdowns',width=105,css_classes=["custom_button"],name="refresh_button")
	MW_buttons = RadioButtonGroup(labels=[''],active=0,width=850,name='MW_buttons') # buttons that will switch between the different microwindows; start empty, will be updated later
	# Button callbacks
	lft_button.on_click(setup_linefit)
	save_button.on_click(save_session)
	load_button.on_click(load_session)
	loop_button.on_click(linefit_loop)
	refresh_button.on_click(update_dropdowns)
	MW_buttons.on_change('active',change_microwindow)
	
	# Text
	status_text = Div(text='<font size=2 color="teal"><b>Status:</b></font>',name="status_text")
	status_div = Div(text='Select a spectrum',width=300,name="status_div") # will display information on app status
	cur_spec_div = Div(text="<font size=3 color='teal'><b>Spectrum</b></font>",width=400,name="cur_spec_div") # will display the current spectrum
	cur_spec_div2 = Div(text="<font size=3 color='teal'><b>Spectrum</b></font>",width=400,name="cur_spec_div2") # duplicate widget for the averaging kernels panel
	suptitle = Div(text='<font size=5 color="teal"><b>Linefit 14.5</b></font>',width=150,name='suptitle') # big title displayed at the top of the webpage
	# Loader gif
	loader = Div(text="",width=40,height=40,name="loader")
	# Spacing DIVs
	space_div = Div(text='',height=30,name="dum_div")
	space_div2 = Div(text='',height=15,name="dum_div")
	dum_div = Div(text='',height=10,name="dum_div")
	dum_div2 = Div(text='',height=10,name="dum_div2")
	dum_div3 = Div(text='',height=15,name="dum_div3")
	# Separation lines
	line_div,line_div2,line_div3,line_div4 = [linediv(width=220) for i in range(4)]

	## FIGURES
	# Modulation efficiency
	ME_fig = figure(plot_width=600,plot_height=175,tools=TOOLS,active_drag="box_zoom",min_border_left=100,min_border_bottom=40,name="ME_fig")
	ME_fig.yaxis.axis_label = 'Modulation Efficiency'
	ME_fig.xaxis.axis_label = 'OPD (cm)'
	# Phase Error
	PE_fig = figure(plot_width=600,plot_height=175,tools=TOOLS,active_drag="box_zoom",min_border_left=100,min_border_bottom=40,x_range=ME_fig.x_range,name="PE_fig")
	PE_fig.yaxis.axis_label = 'Phase Error (rad)'
	PE_fig.xaxis.axis_label = 'OPD (cm)'
	# Column
	column_fig = figure(plot_width=600,plot_height=175,tools=TOOLS,active_drag="box_zoom",min_border_left=100,min_border_bottom=40,name="column_fig")
	column_fig.yaxis.axis_label = 'column scale factor'
	column_fig.xaxis.axis_label = 'Microwindow #'
	# ME at MOPD time series
	series_fig = figure(plot_width=600,plot_height=175,tools=TOOLS,active_drag="box_zoom",min_border_left=100,x_axis_type='datetime',min_border_bottom=40,name="series_fig")
	series_fig.xaxis.axis_label = 'Date'
	series_fig.yaxis.axis_label = 'ME at MOPD'
	# ILS
	ILS_fig = figure(title='ILS',plot_width=350,plot_height=360,min_border_left=80,min_border_bottom=50,min_border_right=30,y_range=DataRange1d(start=-25, end=100),x_range=DataRange1d(start=-1.2,end=1.2),tools=TOOLS,active_drag="box_zoom",name="ILS_fig")
	ILS_fig.yaxis.axis_label = 'Response'
	ILS_fig.xaxis.axis_label = 'Wavenumber (cm-1)'
	# Microwindows and residuals
	mw_fig = figure(plot_width=450,plot_height=200,tools=TOOLS,active_drag="box_zoom",min_border_bottom=30,title='Microwindow 1',name="mw_fig")
	resid_fig = figure(x_range=mw_fig.x_range,plot_width=450,plot_height=170,min_border_bottom=50,y_range=DataRange1d(start=-1,end=1),tools=TOOLS,active_drag="box_zoom",name="resid_fig")
	resid_fig.yaxis.axis_label = '% Residuals'
	resid_fig.xaxis.axis_label = 'Wavenumber (cm-1)'
	resid_fig.title.text = 'RMS = '
	# Averaging kernels
	AKapo_fig = figure(title='AK apo',plot_width=450,plot_height=400,min_border_left=80,min_border_bottom=50,min_border_right=30,x_range=DataRange1d(start=-0.8,end=1.1),tools="box_select,tap,pan,box_zoom,wheel_zoom,redo,undo,reset,save",active_drag="box_select",name="AKapo_fig")
	AKphase_fig = figure(title='AK phase',plot_width=450,plot_height=400,min_border_left=80,min_border_bottom=50,min_border_right=30,x_range=AKapo_fig.x_range,y_range=AKapo_fig.y_range,tools="box_select,tap,pan,box_zoom,wheel_zoom,redo,undo,reset,save",active_drag="box_select",name="AKphase_fig")
	for elem in [AKapo_fig,AKphase_fig]:
		elem.xaxis.axis_label = 'AK'
		elem.yaxis.axis_label = 'OPD (cm)'

	if not ignore_spec:
		# Spectrum
		spec_fig = figure(title='Spectrum,',plot_width=800,plot_height=360,min_border_left=80,min_border_bottom=50,min_border_right=30,tools=TOOLS,active_drag="box_zoom",name="spec_fig")
		spec_fig.yaxis.axis_label = 'Intensity (??)'
		spec_fig.xaxis.axis_label = 'Wavenumber (cm-1)'
		spec_source = ColumnDataSource(data={'x':[],'y':[]})
		spec_fig.line(x='x',y='y',source=spec_source,name="spec_line")
	
	## SOURCES
	ILS_source = ColumnDataSource(data={'x':[],'y':[]},name="ILS_source")
	mw_source = ColumnDataSource(data={'x':[],'meas':[],'calc':[],'resid':[]},name="mw_source")
	ak_source = ColumnDataSource(data={'x':[],'AKapo':[],'AKphase':[],'opd':[],'fix':[],'color':[]},name="ak_source")

	## LINES
	# Microwindows
	mw_fig.line(x='x',y='meas',color='blue',legend='measured',source=mw_source,name="meas_line")
	mw_fig.line(x='x',y='calc',color='red',legend='calculated',source=mw_source,name="calc_line")
	mw_fig.legend.location="bottom_right"
	mw_fig.legend.click_policy = 'hide'
	# Residuals
	resid_fig.line(x='x',y='resid',color='black',source=mw_source,name="resid_line")
	# ILS
	ILS_fig.line(x='x',y='y',source=ILS_source,name="ILS_line")
	# AK
	AKapo_fig.multi_line(xs='AKapo',ys='x',color='color',source=ak_source,name="AKapo_line")
	AKphase_fig.multi_line(xs='AKphase',ys='x',color='color',source=ak_source,name="AKphase_line")
	AKapo_fig.scatter(x='fix',y='opd',color='color',source=ak_source,name="AKapo_scatter")
	AKphase_fig.scatter(x='fix',y='opd',color='color',source=ak_source,name="AKphase_scatter")
	
	## LEGEND
	dum_fig = figure(plot_width=265,plot_height=850,outline_line_alpha=0,toolbar_location=None,name='dum_fig')
	dum_fig.x_range.end = 1005
	dum_fig.x_range.start = 1000
	dum_fig.renderers = [Legend(click_policy='hide',location='top_left',border_line_alpha=0)]

	## Laying out plot objects
	# Grid for modulation efficiency, phase error, column scale factor, and ME at MOPD time series
	MEPECOL_grid = gridplot([[space_div2],[ME_fig],[PE_fig],[column_fig],[series_fig]],toolbar_location='left',name="MEPECOL_grid")
	# Panel for the MEPECOL grid and the legend
	MEPECOL_panel = Panel(child=gridplot([[MEPECOL_grid,dum_fig]],toolbar_location=None),title='Summary',name="MEPECOL_panel")
	
	# Subgrid with the microwindow and residuals figures
	mw_grid = gridplot([[mw_fig],[resid_fig]],toolbar_location=None)
	# Subgrid2 with the ILS figure and the mw_grid subgrid 
	if ignore_spec:
		ils_fits_grid = gridplot([[ILS_fig,mw_grid]],toolbar_location='left')
	else:
		ils_fits_grid = gridplot([[ILS_fig,mw_grid],[spec_fig]],toolbar_location='left')
	# Grid with the 'cur_spec_div', the buttons for microwindows and the 'ils_fits_grid'
	ils_fits_grid = gridplot([[cur_spec_div],[MW_buttons],[ils_fits_grid]],toolbar_location=None)
	# Panel for the ils_fits_grid
	ils_fits_panel = Panel(child=ils_fits_grid,title='ILS and fits',name="ils_fits_panel")

	# Grid for the ME and phase averaging kernels
	ak_grid = gridplot([[AKapo_fig,AKphase_fig]],toolbar_location='left')
	ak_grid = gridplot([[cur_spec_div2],[ak_grid]],toolbar_location=None)

	# Panel of the ak_grid
	ak_panel = Panel(child=ak_grid,title='Averaging kernels')

	# put the ils_fits_panel and MEPECOL_panel in a Tabs() object
	final = Tabs(tabs=[MEPECOL_panel,ils_fits_panel,ak_panel],width=920,name='final')

	# put all the widgets in a widget box
	widget_box = widgetbox(space_div,refresh_button,session_input,load_button,line_div,dum_div,spec_input,dum_div2,reg_input,line_div2,lft_button,line_div4,save_input,save_button,line_div3,loop_input,loop_button,dum_div3,loader,status_text,status_div,css_classes=['side_widgets'],name="widget_box")

	# empty widget box. After linefit is run, it will be filled with buttons that select the spectrum to be displayed in the ils_fits_panel
	button_box = Column(children=[widgetbox(width=255)],name='button_box')

	# put the widget_box in a grid
	side_box = gridplot([[widget_box]],toolbar_location=None,name="side_box")

	# put the page title and the final Tabs() in a grid
	sub_grid = gridplot([[suptitle],[final]],toolbar_location=None,name="sub_grid")

	# put 'sub_grid', the button_box, and 'side_box' in a grid
	grid = gridplot([[sub_grid,button_box,side_box]],toolbar_location=None,name="grid")

	# add that grid to the document
	curdoc().add_root(grid)

	# use the all_data dictionary to fill all the plots with lines (when a previous session is loaded)
	update_doc()

def linefit_loop():
	'''
	run linefit for all the spectra that include in their name the keyword given in the 'loop_input'
	'''

	keyword = curdoc().select_one({"name":"loop_input"}).value

	reg_key = re.compile(keyword.replace('*','.*'),re.IGNORECASE)

	spec_input = curdoc().select_one({"name":"spec_input"})
	status_div = curdoc().select_one({"name":"status_div"})

	# get list of spectra in lft_app/spectra that include the keyword in their name
	select_spectra = [elem for elem in spec_input.options if reg_key.match(elem)]

	# loop over those spectra and run linefit for each
	for i,spectrum in enumerate(select_spectra):
		spec_input.value = spectrum
		setup_linefit()
		status_div.text = 'Spectrum {}/{} done'.format(i+1,len(select_spectra))
		time.sleep(3) # I put as small delay to let the lines render before the next iteration

	save_session() # save the current session when the loop is finished

	status_div.text = "{} loop finished, {} spectra analysed".format(keyword,len(select_spectra))

# this is displayed in the browser tab
curdoc().title = 'LINEFIT 14.5'

# fill the document
doc_maker()