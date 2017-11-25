#!/var/lib/py27_sroche/bin/python
 # -*- coding: utf-8 -*-

from __future__ import print_function # allows the use of Python 3.x print(function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# Code Description #
####################

'''
Produces interactive .html spectra using files in a given folder.
Can be run with or without arguments. To run with arguments:
python TCCON_spectra.py arg1 arg2
arg1 is the /path/to/spectra
arg2 is a common chain of characters in the name of the spectra to be plotted
if arg2 = n it will make plots for all the spectra in /path/to/spectra
Plots will be saved in /path/to/spectra/SAVE

Simple modifications:
- The color of the species are determined based on the dictionary "colors". If you retrieve a window that include species that are not specified in that dictionary, you need to add them with an associated color or the code will give a KeyError
- If you prefer scatter plots instead of line plots, just search and replace ".line" with ".scatter" 
'''

####################
# import libraries #
####################

# manipulate paths
import os.path

# special arrays with special functions
import numpy as np

# prompt interactions
import sys

# interactive html plots with bokeh
from bokeh.plotting import figure
from bokeh.models import Legend, CustomJS, ColumnDataSource, HoverTool, CheckboxGroup, Button, Range1d
from bokeh.layouts import gridplot,widgetbox
from bokeh.resources import CDN
from bokeh.embed import file_html

# to declare dictionaries with sorted keys
from collections import OrderedDict

#############
# Functions #
#############

# loadbar to be displayed in the prompt to keep track of a loop's iterations
def progress(i,tot,bar_length=20,word=''):
	if tot==0:
		tot=1
	percent=float(i+1)/tot
	hashes='#' * int(round(percent*bar_length))
	spaces=' ' * (bar_length - len(hashes))
	sys.stdout.write("\rPercent:[{0}] {1}%".format(hashes + spaces, int(round(percent * 100)))+"    "+str(i+1)+"/"+str(tot)+'  Now plotting: '+word)
	sys.stdout.flush()

# spt files are the spectrum files output by GFIT/GFIT2
def read_spt(path):

	DATA = {}

	infile=open(path,'r')
	content=infile.readlines()
	infile.close()

	head = content[2].split()

	DATA['header'] = head

	content_T = np.array([[elem for elem in line.split()] for line in content[3:]],dtype=np.float).T # transpose of the file content after the header, so content_T[i] is the ith column

	DATA['columns'] = {}
	for var in head:
		DATA['columns'][var] = content_T[head.index(var)]

	DATA['sza'] = float(content[1].split()[3])
	DATA['zobs'] = float(content[1].split()[4])

	resid = 100.0*(DATA['columns']['Tm']-DATA['columns']['Tc']) #the % residuals, tm and tc are transmittances so we just need to multiply by 100 to get %
	rms_resid = np.sqrt(np.mean(np.square(resid)))  #rms of residuals

	DATA['resid'] = resid
	DATA['rms_resid'] = rms_resid

	DATA['params'] = content[1].split()

	return DATA

#########
# Setup #
#########

# hardcode colors of elements
# this should include all the standard tccon species
# it is ok for different species to have the same color if they are not retrieved in the same window
colors = {
		'co2':'red',
		'lco2':'red',
		'wco2':'red',
		'2co2':'olive',
		'3co2':'hotpink',
		'4co2':'indigo',
		'0co2':'green',
		'ch4':'green',
		'co':'darkred',
		'th2o':'blue',
		'h2o':'blue',
		'hdo':'cyan',
		'hcl':'magenta',
		'hf':'pink',
		'n2o':'darkorange',
		'o2':'purple',
		'ao2':'purple',
		'bo2':'purple',
		'0o2':'green',
		'solar':'goldenrod',
		'other':'salmon',
		}

argu = sys.argv # commandline arguments

path='not_a_directory_1e7fwf8wrf78wf' #this should just be a non-existent directory name
while os.path.isdir(path)==False:
	if len(argu)>1:
		path = argu[1] # first commandline argument is the path to the spectra
	else:
		# if no argument is given, ask for the path
		print('Please create a folder and put your processed TCCON spectra in it\n')
		path=raw_input('Give the path to your folder /YOUR/PATH/TO/SPECTRA   :\n')
	if os.path.isdir(path)==False:
		print('/!\\ You gave a wrong path /!\\\n') #error message if the given path doesn't exist)
		if len(argu)>1:
			sys.exit()

#path='/home/sroche/TESTPATH' #if you want to hardcode the path, comment out the while loop above and just change this line.
save_path=os.path.join(path,'SAVE')
if not os.path.isdir(save_path):
	os.makedirs(save_path)

spectra = os.listdir(path) # list of everything in the given directory

# if A = N or n, select all spectra in the given directory (all files with a '.', so the directory must only have folders and spectra files !).
# if A = Y or y, ask the question on selection after the program starts.
# otherwise select all spectra that include the given keyword in their name.
keyword = 'N'
if len(argu)>2:
	if argu[2] in ['N','n']:
		keyword = argu[2]
	else:
		keyword = 'Y'
		select_spectra = sorted([i for i in spectra if (argu[2] in i)])
else:
	keyword = raw_input('Do you want to select specific spectra? (Y/N)\n')
	if keyword in ['Y','y']:
		select=raw_input('Type the string that must be contained in the spectrum file name:\n')
		select_spectra=sorted([i for i in spectra if (select in i)])

if keyword in ['N','n']:
	select_spectra=sorted([i for i in spectra if ('.' in i)])

#############
# Main code #
#############

# loop over the selected spectra
for spectrum in select_spectra:

	progress(select_spectra.index(spectrum),len(select_spectra),word=spectrum) #fancy loadbar

	#read the spectrum file
	spt_data = read_spt(os.path.join(path,spectrum))
	header = spt_data['header']
	
	species = np.array([spt_data['columns'][var] for var in header])
	SZA = str(spt_data['sza'])
	zobs = str(spt_data['zobs'])

	freq=species[0] # the frequency list
	tm=species[1]	# measured transmittance list
	tc=species[2]	# calculated transmittance list
	residuals = spt_data['resid'] # 100*(calculated - measured)
	sigma_rms = spt_data['rms_resid'] # sqrt(mean(residuals**2))

	## start bokeh plot
	TOOLS = "box_zoom,wheel_zoom,pan,undo,redo,reset,crosshair,save" #tools for bokeh figures

	# spectrum figure 
	fig = figure(title=spectrum+'; SZA='+SZA+'°; zobs='+zobs+'km; %resid=100*(Measured-Calculated); RMSresid='+('%.4f' % sigma_rms)+'%',plot_width = 1000,plot_height=400,tools=TOOLS,y_range=Range1d(-0.04,1.04),outline_line_alpha=0)
	# residual figure
	fig_resid = figure(plot_width=1000,plot_height=150,x_range=fig.x_range,tools=TOOLS,y_range=Range1d(-3,3),outline_line_alpha=0)

	# axes labels
	fig_resid.xaxis.axis_label = 'Wavenumber (cm-1)'
	fig_resid.yaxis.axis_label = '% Residuals'
	fig.yaxis.axis_label = 'Transmittance'
	
	N_plots = range(len(species)-1) # a range list from 0 to the number of plots, used by the checkbox group
	
	# group of checkboxes that will be used to toggle line and HoverTool visibility
	checkbox = CheckboxGroup(labels=[header[j+3] for j in range(len(species)-3)]+['Measured','Calculated'],active=N_plots,width=200)
	
	# plotting species lines
	plots = []
	for j in range(len(species)-3):
		try:
			plots.append(fig.line(x=freq,y=species[j+3],color=colors[header[j+3]],line_width=2,name=header[j+3]))
		except KeyError:
			print('KeyError:',header[j+3],'is not specified in the "colors" dictionary, you need to add it with an associated color')
			sys.exit()
		# each line has a associated hovertool with a callback that looks at the checkboxes status for the tool visibility.
		fig.add_tools( HoverTool(mode='vline',line_policy='prev',renderers=[plots[j]],names=[header[j+3]],tooltips=OrderedDict( [('name',header[j+3]),('index','$index'),('(x;y)','(@x{0.00} ; @y{0.000})')] ) ) )

	# adding the measured spectrum
	plots.append(fig.line(x=freq,y=tm,color='black',line_width=2,name='Tm'))
	fig.add_tools( HoverTool(mode='vline',line_policy='prev',renderers=[plots[j+1]],names=['Tm'],tooltips=OrderedDict( [('name','Measured'),('index','$index'),('(x;y)','(@x{0.00} ; @y{0.000})')] ) ) )
	
	# adding the calculated spectrum
	plots.append(fig.line(x=freq,y=tc,color='chartreuse',line_width=2,name='Tc'))
	fig.add_tools( HoverTool(mode='vline',line_policy='prev',renderers=[plots[j+2]],names=['Tc'],tooltips=OrderedDict( [('name','Calculated'),('index','$index'),('(x;y)','(@x{0.00} ; @y{0.000})')] ) ) )

	# legend outside of the figure
	fig_legend=Legend(items=[(header[j+3],[plots[j]]) for j in range(len(species)-3)]+[('Measured',[plots[-2]]),('Calculated',[plots[-1]])],location=(0,0),border_line_alpha=0)
	fig.add_layout(fig_legend,'right')
	fig.legend.click_policy = "hide"
	fig.legend.inactive_fill_alpha = 0.6

	# now the residual figure
	fig_resid.line(x=freq,y=residuals,color='black',name='residuals')
	fig_resid.add_tools(HoverTool(mode='vline',line_policy='prev',names=['residuals'],tooltips={'index':'$index','(x;y)':'($x{0.00} ; $y{0.000})'}))

	# set up a dummy legend for the residual figure so that it aligns with the spectrum figure
	dummy = fig_resid.line(x=freq,y=[0 for i in range(len(freq))],color='white',visible=False,alpha=0)
	fig_resid_legend=Legend(items=[('               ',[dummy])],location=(0,0),border_line_alpha=0)
	fig_resid.add_layout(fig_resid_legend,'right')
	
	# checkbox group callback
	checkbox_iterable = [('p'+str(i),plots[i]) for i in N_plots]+[('checkbox',checkbox)]
	checkbox_code = ''.join(['p'+str(i)+'.visible = checkbox.active.includes('+str(i)+');' for i in N_plots])
	checkbox.callback = CustomJS(args={key: value for key,value in checkbox_iterable}, code=checkbox_code)

	# button to uncheck all checkboxes
	clear_button = Button(label='Hide all lines',width=200)
	clear_button_code = """checkbox.active=[];"""+checkbox_code
	clear_button.callback = CustomJS(args={key: value for key,value in checkbox_iterable}, code=clear_button_code)

	# button to check all checkboxes
	check_button = Button(label='Show all lines',width=200)
	check_button_code = """checkbox.active="""+str(N_plots)+""";"""+checkbox_code
	check_button.callback = CustomJS(args={key: value for key,value in checkbox_iterable}, code=check_button_code)

	# put all the widgets in a box
	group=widgetbox(clear_button,check_button,width=120)

	# define the grid with the figures and widget box
	grid = gridplot([[fig,group],[fig_resid]],toolbar_location='left')

	# write the HTML file
	outfile=open(os.path.join(save_path,spectrum+'.html'),'w')
	outfile.write(file_html(grid,CDN,spectrum[:12]+spectrum[-3:]))
	outfile.close()
	## end bokeh plot
print('\n')

sys.exit() # to make sure the program doesn't hang after it's finished
