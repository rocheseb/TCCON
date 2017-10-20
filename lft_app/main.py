#!/usr/bin/env python2.7
 # -*- coding: ascii -*-

from __future__ import print_function # allows the use of Python 3.x print function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# Code Description #
####################
'''
The folder 'lft_app' should be in the same directory as lft145.exe

- reads all the .DPT ( not OPUS format !!! ) spectra in the 'spectra' folder of 'lft_app'
- those spectra must have been cut and ratioed to ~1 beforehand
- reads the scanner temperature for each spectrum in the 'temp' file of the 'spectra' folder
- modify linefit input file and runs linefit
- plot Modulation efficiency and phase error vs OPD
- plot column scale factor vs microwindow
- plot ILS and fits in each microwindow

You must create a file named 'temp' (with no extensions) in /lft_app/spectra
In it you should list the spectrum file names with associated temperatures like this:

spectrumfilename1,temperature1
spectrumfilename2,temperature2

etc.

DISCLAIMER: if any warning or error message is given by linefit, this app will hang, you should then run linefit from the terminal to figure out what the problem is
The app may hang if there is any convergence problem, or if a significant spectral detuning is detected.
If a significant spectral detuning is detected. Run linefit from the terminal, notice the value of spectral residuals given after the warning, and add it to the input file for each microwindow

e.g. for a spectral detuning of -2.67E-06

species parameters:
    for each species:
        gas T,fit of gas T (F/T),column of species [m-2], ptot [mbar], fit of total pressure (F/T), ppart[mbar], default gamma
            (cell column = 7.243e24 * p[mbar] * l[m] / T[K])
            (first-guess values in case of retrieval)
        for each MW: take species into account(T/F),column scaling factor, spectral scaling factor of species - 1

$
293.15,.false.,1.3310e+22,4.62,.false.,4.62,0.0075
.true.,1.0,-2.67E-06
'''
####################
# Import libraries #
####################

import os
import sys
import subprocess

import collections
from collections import OrderedDict

import time
from datetime import datetime, timedelta

import bokeh
from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, CustomJS, Button, Div, TextInput, Select, Panel, Tabs, Legend, DataRange1d
from bokeh.layouts import gridplot, widgetbox, Column

import numpy as np

from cell_data import hcl_cells, hbr_cells

import pylab as pl
import matplotlib.backends.backend_pdf as mpl_pdf

#############
# Functions #
#############

def corP(P,T):
	Tin = float(T)
	Pref = float(P)
	return str(Pref*Tin/296.0)

#cell column = 7.243e24 * p[mbar] * l[m] / T[K])
# I have noticed that using that formula often gives worst column scale factors than not using it
def corcol(P,T,l=0.1):
	Tin = float(T)
	Pin = float(P)
	return str(7.243e24*Pin*l/Tin)

#########
# Setup #
#########

hcl_cell_data = hcl_cells()
hbr_cell_data = hbr_cells()

# add colors to the kelly_colors list if you want to have more than 20 lines available
kelly_colors = [ 	'#F3C300','#875692', '#F38400', '#A1CAF1','#BE0032', '#C2B280', '#848482','#008856', '#E68FAC', '#0067A5',
				 	'#F99379', '#604E97', '#F6A600','#B3446C', '#DCD300', '#882D17','#8DB600', '#654522', '#E25822','#2B3D26',		]

# for example do something like below
kelly_colors += ['blue','red','green','chartreuse','magenta','firebrick']

WINDOWS = {}

# MIR mircrowindows for HBr
WINDOWS['HBr_test'] = [
				(2590.32, 2590.72), # 1
				(2590.71, 2591.11), # 2
				(2605.60, 2606.00), # 3
				(2606.00, 2606.40), # 4
				(2620.39, 2620.79), # 5
				(2620.80, 2621.20), # 6
				(2634.70, 2635.10), # 7
				(2635.10, 2635.50), # 8
				(2648.50, 2648.90), # 9
				(2648.90, 2649.30), # 10
				(2661.76, 2662.16), # 11
				(2662.18, 2662.58), # 12
				(2674.52, 2674.92), # 13
				]

# MIR microwindows for HCl
WINDOWS['HCl_test'] = [
				(2843.1,2844.1), # 1
				#(2596.8,2597.8),
				#(2625.2,2626.2),
				(2862.5,2863.5), # 2
				(2864.6,2865.6), # 3
				(2903.6,2904.6), # 4
				(2905.7,2906.7), # 5
				(2923.2,2924.2), # 6
				(2923.2,2924.2), # 7
				(2923.2,2924.2), # 8
				(2925.4,2926.4), # 9
				(2960.6,2961.6), # 10
				#(2962.8,2963.8), # 11
				(3011.6,3012.6), # 12
				(3013.9,3014.9), # 13
				(3027.3,3028.3), # 14
				#(3029.6,3030.6), # 15
				]

# NIR microwindows for HCL
WINDOWS['HCl'] = [	
				(5683.0,5684.0), # 1 
				(5687.1,5688.1), # 2
				(5701.5,5702.5), # 3
				(5705.6,5706.6), # 4
				(5718.7,5719.7), # 5 
				(5734.6,5735.6), # 6 
				(5738.8,5739.8), # 7
				(5749.3,5750.3), # 8
				(5753.5,5754.5), # 9
				(5762.7,5763.7), # 10
				(5766.9,5767.9), # 11
				(5774.8,5775.8), # 12
				(5779.0,5780.0), # 13
				]

app_path = os.path.dirname(__file__) # the app should be in /linefit/lft145/lft_app
specpath = os.path.join(app_path,'spectra') # /linefit/lft145/lft_app/spectra
save_path = os.path.join(app_path,'saved_sessions') # /linefit/lft145/lft_app/saved_sessions

wdir = os.sep.join(app_path.split(os.sep)[:-1]) #get the working directory; path to /linefit/lft145

ergpath = os.path.join(wdir,'ergs') # /linefit/lft145/ergs

spec_input_code = """
if (cb_obj.value===""){
	status_div.text = 'Select a Spectrum';
} else {
	status_div.text='Use the button to run linefit with the selected spectrum';
}
"""

TOOLS = "box_zoom,wheel_zoom,pan,redo,undo,reset,save" # tools that will be available to interact with the plots

source_dict = {'ID':0} # this dictionary will store all the data for plots; 'ID' will store the ID of the appriopriate color from 'kelly_colors'

def dumfig(width=600,height=600,legend={}):
	'''
	Need to make a dummy figure to get the legend somewhere by itself ....
	
	legend is a dict/ordereddict with strcuture:
		- {'name':'color'}
	or  - {'name':{'legend':'name','line_dash':'...','color':'...'}}
	or  - {'name':{'legend':'name','marker':'...','color':'...'}}
	'''
	dumx=range(10)
	dumfig=figure(outline_line_alpha=0,plot_height=height,plot_width=width,toolbar_location=None)
	for key in legend:
		if type(legend[key]) in [dict,collections.OrderedDict]:
			try:
				dumfig.line(x=dumx,y=dumx,visible=False,color=legend[key]['color'],line_width=2,line_dash=legend[key]['line_dash'],legend=legend[key]['legend'])
			except KeyError:
				dumfig.scatter(x=dumx,y=dumx,visible=False,color=legend[key]['color'],marker=legend[key]['marker'],legend=legend[key]['legend'])
		else:
			dumfig.line(x=dumx,y=dumx,visible=False,color=legend[key],line_width=2,legend=key)
	dumfig.renderers=[rend for rend in dumfig.renderers if (type(rend)==bokeh.models.renderers.GlyphRenderer or type(rend)==bokeh.models.annotations.Legend)]
	dumfig.renderers[0].border_line_alpha=0
	dumfig.renderers[0].spacing=6
	dumfig.renderers[0].location='top_left'
	for rend in dumfig.renderers:
		if type(rend)==bokeh.models.renderers.GlyphRenderer:
			rend.visible = False

	return dumfig

def modify_input_file(spectrum,Temperature,mwindows):
	'''
	Update the linefit input file to correspond to the selected spectrum and regularisation factor.
	For each site, the cell information must be added to the cell_data.py file
	.dpt spectra must be on lft_app/spectra with the following naming convention

	HCl_MAXOPD_specnum_site_yymmdd.dpt

	e.g. HCl_45_0_eu_170315.dpt for a HCl spectrum from Eureka on March 15 2017, with 45cm MOPD, and spectrum number = 0

	If the spectrum has significant spectral detuning (can happen for spectra in a vented instrument), put the number before HCl
	e.g. 0_HCl_45_eu_170315.dpt
	'''

	reg = curdoc().select_one({'name':'reg_input'}).value
	reg_phase_mod = (reg,reg)

	if '_45_' in spectrum:
		OPDmax = '45.0'
	elif '_65_' in spectrum:
		OPDmax = '65.0'
	elif '_75_' in spectrum:
		OPDmax = '75.0'

	for site in hcl_cell_data: # check the cell data for the appriopriate site
		if '_'+site+'_' in spectrum:
			curdoc().select_one({"name":"status_div"}).text += '<br>- '+hcl_cell_data[site]['location']+' cell specs found'
			break
	else: # if the site is not specified, default to Eureka (I need to change that)
		print("\nSite not recognized: filename must be HCl_OPD_site_YYMMDD.dat , with 'site' the two letter site abbreviation")
		print("\nCheck that the cell information for the site is entered correctly in cell_data.py")
		curdoc().select_one({"name":"status_div"}).text += '<br>- Site not recognized'
		return
	print(site)

	if 'hcl' in spectrum.lower():
		infile = open(os.path.join(app_path,'lft14_hcl.inp'),'r')
		content = infile.readlines()
		infile.close()

		for i in range(len(content)):
			if 'spectra\\' in content[i]:
				content[i] = 'lft_app\\spectra\\'+spectrum+'\n' #path to spectrum for each microwindow
			if 'number of microwindows' in content[i]:
				for wid,window in enumerate(mwindows):
					content[i+5+wid] = str(window)+'\n'
			if '1.2909e22' in content[i]:
				newP = corP(hcl_cell_data[site]['effp_h35cl_296k'],Temperature) # update the cell pressure
				newcol = corcol(newP,Temperature) # update the cell column as described in the input file
				print(newP,newcol)
				#newP = hcl_cell_data[site]['effp_h35cl_296k']  # use the effective pressure from the TCCON wiki
				newcol = hcl_cell_data[site]['h35cl_column'] # use the column from the TCCON wiki
				print(newP,newcol)
				content[i] = Temperature+',.false.,'+newcol+','+newP+',.false.,'+newP+',0.0075\n' #change retrieval parameters for 1st species (HCl35 for hcl cell)
				if ('_hcl' in spectrum.lower()) and ('ratio' not in spectrum.lower()):
					curdoc().select_one({"name":"status_div"}).text += "<br>- spectral detuning: -2.67E-06"
					for ite in range(1,14):
						content[i+ite] = '.true.,1.0,-2.67E-06\n'			
			if '1.2836e22' in content[i]:
				newP = corP(hcl_cell_data[site]['effp_h37cl_296k'],Temperature) # update the cell pressure
				newcol = corcol(newP,Temperature) # update the cell column as described in the input file
				print(newP,newcol)
				#newP = hcl_cell_data[site]['effp_h37cl_296k'] # use the effective pressure from the TCCON wiki
				newcol = hcl_cell_data[site]['h37cl_column'] # use the column from the TCCON wiki
				print(newP,newcol)
				content[i] = Temperature+',.false.,'+newcol+','+newP+',.false.,'+newP+',0.0075\n'  #change retrieval parameters for 2nd species (HCl37 for hcl cell)
				if ('_hcl' in spectrum.lower()) and ('ratio' not in spectrum.lower()):
					for ite in range(1,14):
						content[i+ite] = '.true.,1.0,-2.67E-06\n'
			if 'gas cell parameters' in content[i]:
				content[i+13] = Temperature+'\n'
			if '1.196e-3' in content[i]:
				content[i-2] = OPDmax+'\n'	# change OPDmax value
			if 'Reg Modulation' in content[i]:
				content[i+4] = '1.0e4,%s,0.0,%s,0.0\n' % reg_phase_mod

	if 'hbr' in spectrum.lower():
		# the hbr section is not complete
		infile = open(os.path.join(app_path,'lft14_hbr.inp'),'r')
		content = infile.readlines()
		infile.close()

		for i in range(len(content)):
			if 'spectra\\' in content[i]:
				content[i] = 'spectra\\EUREKA-HCL\\'+spectrum+'\n' #path to spectrum for each microwindow
			if 'number of microwindows' in content[i]:
				for wid,window in enumerate(mwindows):
					content[i+5+wid] = str(window)+'\n'
			if '1.2909e22' in content[i]:
				newP = corP(hbr_cell_data[site]['pressure'],Temperature) # update the cell pressure
				newcol = corcol(newP,Temperature) # update the cell column as described in the input file
				# newP = hbr_cell_data[site]['pressure']
				newcol = hbr_cell_data[site]['column'] # cell column
				content[i] = Temperature+',.false.,'+newcol+','+newP+',.false.,'+newP+',0.0075\n' #change retrieval parameters for HBr cell
			if 'gas cell parameters' in content[i]:
				content[i+13] = Temperature+'\n'
			if '1.196e-3' in content[i]:
				content[i-2] = OPDmax+'\n'	# change OPDmax value		
			if 'Reg Modulation' in content[i]:
				content[i+4] = '1.0e4,%s,0.0,%s,0.0\n' % reg_phase_mod

	outfile = open(os.path.join(wdir,'lft14.inp'),'w') #rewrite input file
	outfile.writelines(content)
	outfile.close()

	curdoc().select_one({"name":"status_div"}).text+='<br>- Input file updated'
	print('\n\t- Input file updated')

def run_linefit():
	'''
	run linefit for the spectrum selected in spec_input
	'''

	global source_dict

	dum_leg = curdoc().select_one({"name":"dum_leg"})

	if len(dum_leg.legend[0].items)>1:
		source_dict['ID'] += 1

	curdoc().select_one({"name":"spec_input"}).js_on_change('value', CustomJS(args={'status_div':curdoc().select_one({"name":"status_div"})},code=spec_input_code))

	spectrum = curdoc().select_one({"name":"spec_input"}).value

	if spectrum=='':
		status_div.text = "Select a spectrum"
		return

	reg = curdoc().select_one({"name":"reg_input"}).value

	status_div = curdoc().select_one({"name":"status_div"})

	try:
		float(reg)
	except:
		status_div.text = "Regularisation factor must be a number"
		return

	dum_leg_labels = [elem.label['value'] for elem in dum_leg.legend[0].items]
	already_done = [elem for elem in dum_leg_labels if ((spectrum.split('.')[0] in elem) and ('reg='+reg in elem))]!=[]
	if already_done:
		status_div.text = spectrum+" already analysed with reg="+reg 
		return
	
	status_div.text = "<b>Now doing:</b> <br>"+spectrum+'<br>reg='+reg
	print('\nNow doing',spectrum,'with reg=',reg)

	if 'hcl_45' in spectrum.lower():
		mwindows = WINDOWS['HCl']
	if 'hcl_65' in spectrum.lower():
		mwindows = WINDOWS['HCl']

	colo = kelly_colors[source_dict['ID']]

	# get the scanner temperature
	# I write my own inputfile called 'temp' that has lines with 'SpectrumName,ScannerTemperature'
	infile = open(os.path.join(specpath,'temp'),'r')
	content = infile.readlines()
	infile.close()

	speclist = [line.split(',')[0] for line in content] # all the SpectrumName in the file

	#if the spectrum is not listed in the temp file, go to next spectrum
	if spectrum.lower() not in [spec.lower() for spec in speclist]:
		status_div.text = spectrum+':</br>scanner temperature not listed in the temp file'
		print(spectrum,'scanner temperature not listed in the temp file')
		return

	# get the temperature
	Temperature = [line.split(',')[1] for line in content if spectrum.lower() in line.lower()][0].split('\n')[0] #second element in comma separated line minus the line return

	# check that the spectral range is ordered (dpt files are written with decreasing wavenumbers and lienfit wants increasing wavenumbers)
	# if it is not ordered, orders it.
	spec_path = os.path.join(specpath,spectrum)
	infile = open(spec_path,'r')
	content = infile.readlines()
	infile.close()

	specrange = [line.split()[0] for line in content]

	if specrange!=sorted(specrange):
		status_div.text += '<br>- Reordering '+spectrum
		print('\n\t- Reordering',spectrum)

		outfile = open(spec_path,'w')
		outfile.writelines(content[::-1])
		outfile.close()

	ratio_spectrum(spec_path)

	modify_input_file(spectrum,Temperature,mwindows)

	status_div.text+='<br>- Running linefit ...'
	print('\n\t- Running linefit ...')
	proc = subprocess.Popen(['lft145.exe'], stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) # run linefit without shell output
	#give_input = proc.communicate(input="1\n")[0]
	tmp = proc.stdout.read() # this saves the shell output of the linefit code, I don't use it later but it's there if you need to check the messages.

	linefit_results(spectrum,colo)

	curdoc().select_one({'name':'status_div'}).text += "<br><b>DONE</b>"

def ratio_spectrum(spec_path):
	'''
	Fit a second order polynomial to a spectrum in order to ratio it to ~1
	'''
	x,y = np.loadtxt(spec_path,unpack=True)

	if 0.7<np.mean(y)<1.3: # spectrum was already ratioed
		return

	# resample_y_max takes the max y values for each interval of 14 wavenumbers, it avoids getting points in the lines themselves
	resample_y_max = [max(y[(x<(start+15)) & (x>start)]) for start in range(int(x[0]),int(x[-1]),14)]
	resample_x = np.array(range(int(x[0]),int(x[-1]),14))

	fit = np.polyfit(resample_x,resample_y_max,2) # second order polynomial fit

	# I substract the 0.0004 because the way I fit will always be a bit too high as it will include values above the base line
	base_y = fit[0]*x**2+fit[1]*x+fit[2]-0.0004 # use the fit to get the fit line for each x

	new_y = y/base_y # ratio to ~1

	np.savetxt(spec_path,np.transpose([x,new_y])) # overwrite the spectrum with the ratioed one in lft_app/spectra

	curdoc().select_one({"name":"status_div"}).text += '<br>- Spectrum ratioed to ~{:.2f}'.format(np.mean(new_y))


def linefit_results(spectrum,colo):
	'''
	after linefit has finished running, this store the results in source_dict and updates the plots
	'''

	global source_dict

	# string containing info on the current spectrum + current regularisation factor
	cur_name = spectrum.split('.')[0]+' reg='+curdoc().select_one({'name':'reg_input'}).value

	# Add a new entry in the source_dict dictionnary
	source_dict[cur_name] = {'ILS':0,'mw1':0,'mw2':0,'mw3':0,'mw4':0,'mw5':0,'mw6':0,'mw7':0,'mw8':0,'mw9':0,'mw10':0,'mw11':0,'mw12':0,'mw13':0,}

	# Update the title in the diag_panel
	curdoc().select_one({"name":"cur_spec_div"}).text = "<font size=3 color='teal'><b>"+cur_name+"</b></font>"
	# Update legend
	curdoc().select_one({"name":"dum_leg"}).line(range(2),range(2),color=colo,line_width=2,visible=False,legend=cur_name)
	# Update status
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding ME and PE plot'
	
	###############################
	############################### Modulation Efficiency and Phase Error
	print('\n\t- ME and PE')

	infile = open(os.path.join(ergpath,'modulat.dat'),'r')
	content = infile.readlines()
	infile.close()

	content = [line.split() for line in content[1:]]
	content = np.array([[float(elem) for elem in row] for row in content]).T

	source_dict[cur_name]['color'] = colo
	source_dict[cur_name]['ME'] = {'x':content[0],'y':content[1]}
	source_dict[cur_name]['PE'] = {'x':content[0],'y':content[2]}

	cur_date_string = cur_name.split('_')[-1][:6]
	cur_date_time_struct = time.strptime(cur_date_string,'%y%m%d')
	cur_date = datetime(*cur_date_time_struct[:6])
	source_dict[cur_name]['series'] = {'x':[cur_date],'y':[content[1][-1]]}

	ME_source = ColumnDataSource(data=source_dict[cur_name]['ME'])
	PE_source = ColumnDataSource(data=source_dict[cur_name]['PE'])
	series_source = ColumnDataSource(data=source_dict[cur_name]['series'])

	curdoc().select_one({"name":"ME_fig"}).line(x='x',y='y',color=colo,line_width=2,source=ME_source)
	curdoc().select_one({"name":"PE_fig"}).line(x='x',y='y',color=colo,line_width=2,source=PE_source)
	curdoc().select_one({"name":"series_fig"}).scatter(x='x',y='y',color=colo,size=5,source=series_source)

	###############################
	############################### Column
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding column plot'
	print('\n\t- column')

	infile = open(os.path.join(ergpath,'colparms.dat'),'r')
	content = infile.readlines()
	infile.close()

	if 'hcl' in spectrum.lower():
		secondID = [i for i,v in enumerate(content) if 'Species:   2' in v][0] #hcl37, I don't plot it though.
		second_spec = [float(i) for i in content[secondID+1:]]
		first_spec = [float(i) for i in content[1:secondID]]
	else:
		first_spec = [float(i) for i in content[1:]]

	mw = range(1,len(first_spec)+1) # just the microwindow numbers

	source_dict[cur_name]['COL'] = {'x':mw,'y':first_spec}

	COL_source = ColumnDataSource(data=source_dict[cur_name]['COL'])

	curdoc().select_one({"name":"column_fig"}).line(x='x',y='y',color=colo,line_width=2,source=COL_source)
	
	###############################
	############################### ILS
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding ILS'
	print('\n\t- Adding ILS')

	infile = open(os.path.join(ergpath,'ilsre.dat'),'r')
	content = infile.readlines()
	infile.close()

	content = np.array([[float(elem) for elem in line.split()] for line in content]).T
 
	source_dict[cur_name]['ILS'] = {'x':content[0],'y':content[1]}

	ILS_fig = curdoc().select_one({"name":"ILS_fig"})

	curdoc().select_one({"name":"ILS_line"}).data_source.data.update(source_dict[cur_name]['ILS'])
	
	###############################
	############################### Microwindows
	curdoc().select_one({"name":"status_div"}).text+='<br>- Adding spectral fits with residuals'
	print('\n\t- Adding spectral fits with residuals')

	mwfiles = [i for i in os.listdir(ergpath) if 'specre' in i]

	for it,mwfile in enumerate(mwfiles):
		infile = open(os.path.join(ergpath,mwfile),'r')
		content = infile.readlines()
		infile.close()

		content = np.array([[float(elem) for elem in line.split()] for line in content]).T

		resid = 100*(content[1]-content[2])/content[2]

		source_dict[cur_name]['mw'+str(it+1)] = {'x':content[0],'meas':content[1],'calc':content[2],'resid':resid}
		source_dict[cur_name]['rms_resid_mw'+str(it+1)] = "{:.5f}".format(np.sqrt(np.mean(resid**2)))

	mw_fig = curdoc().select_one({"name":"mw_fig"})
	mw_fig.title.text = "Microwindow 1"

	resid_fig = curdoc().select_one({"name":"resid_fig"})
	resid_fig.title.text = 'RMS = '+source_dict[cur_name]['rms_resid_mw1']

	curdoc().select_one({"name":"meas_line"}).data_source.data.update(source_dict[cur_name]['mw1'])
	curdoc().select_one({"name":"calc_line"}).data_source.data.update(source_dict[cur_name]['mw1'])
	curdoc().select_one({"name":"resid_line"}).data_source.data.update(source_dict[cur_name]['mw1'])

	###############################
	###############################

	add_button(cur_name)

	print(spectrum,'DONE')

def add_button(cur_name):
	'''
	add a button to the 'button_box' corresponding to a spectrum + regularisation factor
	Used at the end of the linefit_results function and in the update_doc function
	'''
	global source_dict

	new_button = Button(label=cur_name,width=180)

	button_box = curdoc().select_one({"name":"button_box"})
	button_box.children += [new_button]

	source_dict['cur_name_clicks'] = [0 for i in range(len(button_box.children))]
	source_dict['prev_name_clicks'] = [0 for i in range(len(button_box.children))]

	new_button.on_click(change_spectrum)

def update_doc():
	'''
	Use the current source_dict to fill the plots in the document
	'''
	global source_dict

	# add the buttons corresponding to the current source_dict dictionnary keys
	for cur_name in sorted([key for key in source_dict if 'reg' in key]):

		add_button(cur_name)

		colo = source_dict[cur_name]['color']

		# Fill data sources
		ME_source = ColumnDataSource(data=source_dict[cur_name]['ME'])
		PE_source = ColumnDataSource(data=source_dict[cur_name]['PE'])
		COL_source = ColumnDataSource(data=source_dict[cur_name]['COL'])
		series_source = ColumnDataSource(data=source_dict[cur_name]['series'])
		# Update lines
		curdoc().select_one({"name":"ME_fig"}).line(x='x',y='y',color=colo,line_width=2,source=ME_source)
		curdoc().select_one({"name":"PE_fig"}).line(x='x',y='y',color=colo,line_width=2,source=PE_source)
		curdoc().select_one({"name":"column_fig"}).line(x='x',y='y',color=colo,line_width=2,source=COL_source)
		curdoc().select_one({"name":"series_fig"}).scatter(x='x',y='y',color=colo,size=5,source=series_source)
		# Update legend
		curdoc().select_one({"name":"dum_leg"}).line(range(2),range(2),color=colo,line_width=2,visible=False,legend=cur_name)

def change_spectrum():
	'''
	callback for the spectrum buttons in 'button_box'
	update the plots in 'spec_grid' to correspond to the desired spectrum
	'''
	global source_dict

	button_box = curdoc().select_one({"name":"button_box"})

	# list of current number of clicks for each spectrum button (including the one that just got clicked)
	source_dict['cur_name_clicks'] = [i.clicks for i in button_box.children]
	# compare cur_clicks to prev_clicks to know which button has just been clicked
	for i in range(len(source_dict['cur_name_clicks'])):
		if source_dict['cur_name_clicks'][i]!=source_dict['prev_name_clicks'][i]:
			cur_name = button_box.children[i].label
	# set prev_clicks equal to cur_clicks 
	source_dict['prev_name_clicks'] = [i for i in source_dict['cur_name_clicks']]

	curdoc().select_one({"name":"cur_spec_div"}).text = "<font size=3 color='teal'><b>"+cur_name+"</b></font>"
	
	# Select microwindow and residuals figures
	mw_fig = curdoc().select_one({"name":"mw_fig"})
	resid_fig = curdoc().select_one({"name":"resid_fig"})
	# Get current microwindow
	cur_MW = mw_fig.title.text.split()[1]	
	# Update titles
	resid_fig.title.text = 'RMS = '+source_dict[cur_name]['rms_resid_mw'+cur_MW]
	# Get ILS and microwindow data for the new spectrum
	new_mw_data =  source_dict[cur_name]['mw'+cur_MW]
	new_ILS_data =  source_dict[cur_name]['ILS']
	# Update lines
	curdoc().select_one({"name":"meas_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"calc_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"resid_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"ILS_line"}).data_source.data.update(new_ILS_data)

def change_microwindow():
	'''
	callback for the microwindow buttons
	update the plots in 'mw_grid' so that they correspond to the desired microwindow
	'''
	global source_dict

	mw_buttons = [curdoc().select_one({"name":"MW"+str(it+1)+"_button"}) for it in range(13)]

	# list of current number of clicks for each microwindow button (including the one that just got clicked)
	source_dict['cur_clicks'] = [i.clicks for i in mw_buttons]
	# compare cur_clicks to prev_clicks to know which button has just been clicked
	for i in range(len(source_dict['cur_clicks'])):
		if source_dict['cur_clicks'][i]!=source_dict['prev_clicks'][i]:
			new_MW = str(i+1)
	# set prev_clicks equal to cur_clicks 
	source_dict['prev_clicks'] = [i for i in source_dict['cur_clicks']]

	# Get current spectrum
	cur_spec_div = curdoc().select_one({"name":"cur_spec_div"})
	if 'Spectrum' in cur_spec_div.text:
		curdoc().select_one({"name":"status_div"}).text = "No spectrum selected"
		return	

	# Get the corresponding key for the source_dict dictionnary
	cur_name = cur_spec_div.text.split('<b>')[1].split('</b>')[0]
	
	# Get the new microwindow data
	new_mw_data = source_dict[cur_name]['mw'+new_MW]

	# Select microwindow and residuals figures
	mw_fig = curdoc().select_one({"name":"mw_fig"})
	resid_fig = curdoc().select_one({"name":"resid_fig"})
	# Update titles
	mw_fig.title.text = "Microwindow "+new_MW
	resid_fig.title.text = 'RMS = '+source_dict[cur_name]['rms_resid_mw'+new_MW]
	# Update lines
	curdoc().select_one({"name":"meas_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"calc_line"}).data_source.data.update(new_mw_data)
	curdoc().select_one({"name":"resid_line"}).data_source.data.update(new_mw_data)

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
	for cur_name in sorted(source_dict.keys()):

		if 'reg' in cur_name:

			# ILS
			ILS_source = source_dict[cur_name]['ILS']		
			pl.figure(ilsfig.number)
			pl.plot(ILS_source['x'],ILS_source['y'])
			pl.xlabel('Wavenumber (cm-1)')
			pl.ylabel('Response')
			pl.title(cur_name)
			pl.tight_layout()
			pdf.savefig(ilsfig,bbox_inches='tight')
			pl.clf()

			# MEPECOL
			ME_source = source_dict[cur_name]['ME']
			PE_source = source_dict[cur_name]['PE']
			col_source = source_dict[cur_name]['COL']
			series_source = source_dict[cur_name]['series']
			if mindate>min(series_source['x']):
				mindate = min(series_source['x'])
			if maxdate<max(series_source['x']):
				maxdate = max(series_source['x'])

			# Individual MEPECOL
			pl.figure(mepefig_indiv.number)
			ax_indiv[0].plot(ME_source['x'],ME_source['y'],color=source_dict[cur_name]['color'])
			ax_indiv[1].plot(PE_source['x'],PE_source['y'],color=source_dict[cur_name]['color'])
			ax_indiv[2].plot(col_source['x'],col_source['y'],color=source_dict[cur_name]['color'])
			ax_indiv[0].set_title(cur_name)
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
			newline = ax[0].plot(ME_source['x'],ME_source['y'],color=source_dict[cur_name]['color'])
			mepelines += [newline[0]]
			ax[1].plot(PE_source['x'],PE_source['y'],color=source_dict[cur_name]['color'])
			ax[2].plot(col_source['x'],col_source['y'],color=source_dict[cur_name]['color'])
			ax[3].scatter(series_source['x'],series_source['y'],color=source_dict[cur_name]['color'])

			# microwindows
			pl.figure(mw_fig.number)
			for i in range(1,14):
				mw_source = source_dict[cur_name]['mw'+str(i)]
				meas = mw[0].plot(mw_source['x'],mw_source['meas'],color='blue',label='measured')
				calc = mw[0].plot(mw_source['x'],mw_source['calc'],color='red',label='calculated')
				mw[0].set_title(cur_name+': Microwindow '+str(i))
				resid = mw[1].plot(mw_source['x'],mw_source['resid'],color='black')
				mw[0].set_ylabel('Transmission')
				mw[1].set_ylabel('% Residuals')
				mw[1].set_xlabel('Wavenumber (cm-1)')
				mw[1].set_title('RMS='+source_dict[cur_name]['rms_resid_mw'+str(i)])
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
	lgd = pl.legend(mepelines,[cur_name for cur_name in source_dict if 'reg' in cur_name],bbox_to_anchor=(1.02, 4),loc='center left', borderaxespad=0)
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

def save_session():
	'''
	save the current source_dict
	'''
	global source_dict

	curdoc().select_one({"name":"status_div"}).text = "Saving current session ..."

	save_name = curdoc().select_one({"name":"save_input"}).value

	np.save(os.path.join(save_path,save_name),source_dict) # save the source_dict dictionnary in a .npy file

	pdf_report(save_name) # save a pdf document with all the plots

	print('\nCurrent session data saved in :','lft_app/saved_sessions/'+save_name+'.npy')
	print('\nPDF report saved in :','lft_app/pdf/'+save_name+'.pdf')

	# update the status_div
	curdoc().select_one({"name":"status_div"}).text = 'Current session data saved in:<br>lft_app/saved_sessions/'+save_name+'.npy<br>PDF report saved in:<br>lft_app/pdf/'+save_name+'.pdf'
	
	# now that a new session has been saved, update the options of the session_input widget
	curdoc().select_one({"name":"session_input"}).options = ['']+[i for i in os.listdir(save_path) if '.npy' in i]

def load_session():
	'''
	load a previously saved source_dict and remake the whole document
	'''

	global source_dict

	curdoc().select_one({"name":"status_div"}).text = "Loading new session ..."

	saved_session = os.path.join(save_path,curdoc().select_one({"name":"session_input"}).value)
	
	source_dict = np.load(saved_session).item()

	doc_maker() # rebuild the entire document using the new source_dict

	curdoc().select_one({"name":"spec_input"}).value = ""

	curdoc().select_one({"name":"status_div"}).text = "New session loaded"

def doc_maker():
	'''
	make the whole document
	'''

	global source_dict

	curdoc().clear() # removes everything in the current document

	## WIDGETS
	# Inputs
	spec_input = Select(title='Spectrum:',options = ['']+[i for i in os.listdir(specpath) if (('.dpt' in i) or ('.DPT' in i))],width=150,css_classes=["spec_input"],name="spec_input")
	reg_input = TextInput(value='1.8',title='Regularisation factor:',width=150,css_classes=["small_input"],name="reg_input")
	session_input = Select(title='Previous sessions:',width=150,options=['']+[i for i in os.listdir(save_path) if '.npy' in i],css_classes=["spec_input"],name="session_input")
	save_input = TextInput(title='Save name',value="_".join(str(datetime.now())[:-7].split()),css_classes=["save_input"],name="save_input")
	loop_input = TextInput(title='Loop key',value="HCl_45",width=100,css_classes=["small_input"],name="loop_input")
	# BUTTONS
	lft_button = Button(label='Run linefit', width=80, css_classes=["custom_button"],name="lft_button")
	save_button = Button(label='Save Session', width=90, css_classes=["custom_button"],name="save_button")
	load_button = Button(label='Load Session', width=90, css_classes=["custom_button"],name="load_button")
	loop_button = Button(label='loop',width=90, css_classes=["custom_button"],name="loop_button")
	# Button callbacks
	lft_button.on_click(run_linefit)
	save_button.on_click(save_session)
	load_button.on_click(load_session)
	loop_button.on_click(linefit_loop)
	## Create the buttons that will switch between the different microwindows
	# list of buttons spaced with DIV objects
	MW_buttons = [elem for part in [[Button(label='MW '+str(it+1),width=40,name="MW"+str(it+1)+"_button"),Div(text='',width=27)] for it in range(13)] for elem in part]
	# prev_clicks and cur_clicks are lists of number of clicks for each button, used to know which button was clicked last
	source_dict['prev_clicks'] = [0 for i in range(13)]
	source_dict['cur_clicks'] = [0 for i in range(13)]
	# mw_buttons callbacks
	for elem in MW_buttons:
		# try except to get the button and not the DIVs
		try:
			elem.label
		except AttributeError:
			pass
		else:
			elem.on_click(change_microwindow)
	# Text
	status_text = Div(text='<font size=2 color="teal"><b>Status:</b></font>',name="status_text")
	status_div = Div(text='Select a spectrum',width=300,name="status_div") # will display information on app status
	cur_spec_div = Div(text="<font size=3 color='teal'><b>Spectrum</b></font>",width=400,name="cur_spec_div") # will display the current spectrum
	suptitle = Div(text='<font size=5 color="teal"><b>Linefit 14.5</b></font>',width=300,name='suptitle') # big title displayed at the top of the webpage
	# Spacing DIVs
	space_div = Div(text='',height=30,name="dum_div")
	dum_div = Div(text='',height=10,name="dum_div")
	dum_div2 = Div(text='',height=10,name="dum_div2")
	dum_div3 = Div(text='',height=15,name="dum_div3")
	# Separation lines
	line_div = Div(text='<hr width="100%" color="lightblue">',width= 185,name="line_div")
	line_div2 = Div(text='<hr width="100%" color="lightblue">',width= 185,name="line_div2")
	line_div3 = Div(text='<hr width="100%" color="lightblue">',width= 185,name="line_div3")
	line_div4 = Div(text='<hr width="100%" color="lightblue">',width= 185,name="line_div4")

	## FIGURES
	# Modulation efficiency
	ME_fig = figure(plot_width=600,plot_height=160,tools=TOOLS,min_border_left=100,min_border_bottom=40,name="ME_fig")
	ME_fig.yaxis.axis_label = 'Modulation Efficiency'
	ME_fig.xaxis.axis_label = 'OPD (cm)'
	# Phase Error
	PE_fig = figure(plot_width=600,plot_height=160,tools=TOOLS,min_border_left=100,min_border_bottom=40,x_range=ME_fig.x_range,name="PE_fig")
	PE_fig.yaxis.axis_label = 'Phase Error (rad)'
	PE_fig.xaxis.axis_label = 'OPD (cm)'
	# Column
	column_fig = figure(plot_width=600,plot_height=160,tools=TOOLS,min_border_left=100,min_border_bottom=40,name="column_fig")
	column_fig.yaxis.axis_label = 'column scale factor'
	column_fig.xaxis.axis_label = 'Microwindow #'
	# ME at MOPD time series
	series_fig = figure(plot_width=600,plot_height=160,tools=TOOLS,min_border_left=100,x_axis_type='datetime',min_border_bottom=40,name="series_fig")
	series_fig.xaxis.axis_label = 'Date'
	series_fig.yaxis.axis_label = 'ME at MOPD'
	# ILS
	ILS_fig = figure(title='ILS',plot_width=350,plot_height=360,min_border_left=80,min_border_bottom=50,min_border_right=30,y_range=DataRange1d(start=-25, end=100),x_range=DataRange1d(start=-1.2,end=1.2),tools=TOOLS,name="ILS_fig")
	ILS_fig.yaxis.axis_label = 'Response'
	ILS_fig.xaxis.axis_label = 'Wavenumber (cm-1)'
	# Microwindows and residuals
	mw_fig = figure(plot_width=450,plot_height=200,tools=TOOLS,min_border_bottom=30,title='Microwindow 1',name="mw_fig")
	resid_fig = figure(x_range=mw_fig.x_range,plot_width=450,plot_height=170,min_border_bottom=50,y_range=DataRange1d(start=-1,end=1),tools=TOOLS,name="resid_fig")
	resid_fig.yaxis.axis_label = '% Residuals'
	resid_fig.xaxis.axis_label = 'Wavenumber (cm-1)'
	resid_fig.title.text = 'RMS = '

	## SOURCES
	ILS_source = ColumnDataSource(data={'x':[],'y':[]})
	mw_source = ColumnDataSource(data={'x':[],'meas':[],'calc':[],'resid':[]})

	## LINES
	# Microwindows
	mw_fig.line(x='x',y='meas',color='blue',legend='measured',source=mw_source,name="meas_line")
	mw_fig.line(x='x',y='calc',color='red',legend='calculated',source=mw_source,name="calc_line")
	mw_fig.legend.location="bottom_right"
	# Residuals
	resid_fig.line(x='x',y='resid',color='black',source=mw_source,name="resid_line")
	# ILS
	ILS_fig.line(x='x',y='y',source=ILS_source,name="ILS_line")

	## LEGEND
	dum_leg = dumfig(width=250,height=850,legend={'lft145':'black'})
	dum_leg.name = 'dum_leg'

	## Laying out plot objects
	# Grid for modulation efficiency, phase error, column scale factor, and ME at MOPD time series
	MEPECOL_grid = gridplot([[ME_fig],[PE_fig],[column_fig],[series_fig]],toolbar_location='left',name="MEPECOL_grid")
	# Panel for the MEPECOL grid and the legend
	MEPECOL_panel = Panel(child=gridplot([[MEPECOL_grid,dum_leg]],toolbar_location=None),title='Summary',name="MEPECOL_panel")
	
	# Subgrid with the microwindow and residuals figures
	mw_grid = gridplot([[mw_fig],[resid_fig]],toolbar_location=None)
	# Subgrid2 with the ILS figure and the mw_grid subgrid 
	spec_grid = gridplot([[ILS_fig,mw_grid]],toolbar_location='left')
	# Grid with the 'cur_spec_div', the buttons for microwindows and the 'spec_grid'
	diag_grid = gridplot([[cur_spec_div],[Column(elem) for elem in MW_buttons],[spec_grid]],toolbar_location=None)
	# Panel for the diag_grid
	diag_panel = Panel(child=diag_grid,title='ILS and fits',name="diag_panel")

	# put the diag_panel and MEPECOL_panel in a Tabs() object
	final = Tabs(tabs=[MEPECOL_panel,diag_panel],width=900,name='final')

	# put all the widgets in a widget box
	widget_box = widgetbox(space_div,session_input,load_button,line_div,dum_div,spec_input,dum_div2,reg_input,line_div2,lft_button,line_div4,save_input,save_button,line_div3,loop_input,loop_button,dum_div3,status_text,status_div,css_classes=['side_widgets'],name="widget_box")

	# empty widget box. After linefit is run, it will be filled with buttons that select the spectrum to be displayed in the diag_panel
	button_box = widgetbox(width=210,name="button_box")

	# put the widget_box in a grid
	side_box = gridplot([[widget_box]],toolbar_location=None,name="side_box")

	# put the page title and the final Tabs() in a grid
	sub_grid = gridplot([[suptitle],[final]],toolbar_location=None,name="sub_grid")

	# put 'sub_grid', the button_box, and 'side_box' in a grid
	grid = gridplot([[sub_grid,button_box,side_box]],toolbar_location=None,name="grid")

	# add that grid to the document
	curdoc().add_root(grid)

	# use the source_dict dictionnary to fill all the plots with lines (when a previous session is loaded)
	update_doc()

def linefit_loop():
	'''
	run linefit for all the spectra that include in their name the keyword given in the 'loop_input'
	'''

	keyword = curdoc().select_one({"name":"loop_input"}).value

	spec_input = curdoc().select_one({"name":"spec_input"})

	# get list of spectra in lft_app/spectra that include the keyword in their name
	select_spectra = [elem for elem in spec_input.options if keyword.lower() in elem.lower()]

	# loop over those spectra and run linefit for each
	for spectrum in select_spectra:
		spec_input.value = spectrum
		run_linefit()
		time.sleep(30) # give linefit time to finish running and the plots to finish rendering

	save_session() # save the current session when the loop is finished

	curdoc().select_one({"name":"status_div"}).text = keyword+" loop finished"

# this is displayed in the browser tab
curdoc().title = 'LINEFIT 14.5'

# fill the document
doc_maker()