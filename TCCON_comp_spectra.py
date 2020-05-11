#!/var/lib/py27_sroche/bin/python
 # -*- coding: utf-8 -*-

from __future__ import print_function # allows the use of Python 3.x print(function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

####################
# Code Description #
####################

'''
Produces interactive .html spectra using files in a given folder.

All the spectra in a given folder will be in one html documents, each spectrum will be available in a separate panel.

Spectra sharing the same starting letter will have shared axes (for shared panning/zooming)
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
from bokeh.plotting import figure, output_file
from bokeh.models import Legend, CustomJS, ColumnDataSource, HoverTool, CheckboxGroup, Button, PreText, Range1d, Panel, Tabs, CrosshairTool
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

def read_spt(path):
	"""
	spt files are the spectrum files output by GFIT/GFIT2
	"""

	DATA = {}

	with open(path,'r') as infile:
		content=infile.readlines()

	head = content[2].split()

	DATA['header'] = head

	content_T = np.array([[elem for elem in line.split()] for line in content[3:]],dtype=np.float).T # transpose of the file content after the header, so content_T[i] is the ith column

	DATA['columns'] = {}
	for var in head:
		DATA['columns'][var] = content_T[head.index(var)]

	DATA['sza'] = float(content[1].split()[4])
	DATA['zobs'] = float(content[1].split()[5])

	resid = 100.0*(DATA['columns']['Tm']-DATA['columns']['Tc']) #the % residuals, tm and tc are transmittances so we just need to multiply by 100 to get %
	rms_resid = np.sqrt(np.mean(np.square(resid)))  #rms of residuals

	DATA['resid'] = resid
	DATA['rms_resid'] = rms_resid

	DATA['params'] = content[1].split()

	return DATA

def add_vlinked_crosshairs(fig1, fig2):
	js_move = '''if(cb_obj.x >= fig.x_range.start && cb_obj.x <= fig.x_range.end && cb_obj.y >= fig.y_range.start && cb_obj.y <= fig.y_range.end)
					{ cross.spans.height.computed_location = cb_obj.sx }
				else 
					{ cross.spans.height.computed_location = null }'''
	js_leave = 'cross.spans.height.computed_location = null'

	cross1 = CrosshairTool()
	cross2 = CrosshairTool()
	fig1.add_tools(cross1)
	fig2.add_tools(cross2)
	args = {'cross': cross2, 'fig': fig1}
	fig1.js_on_event('mousemove', CustomJS(args = args, code = js_move))
	fig1.js_on_event('mouseleave', CustomJS(args = args, code = js_leave))
	args = {'cross': cross1, 'fig': fig2}
	fig2.js_on_event('mousemove', CustomJS(args = args, code = js_move))
	fig2.js_on_event('mouseleave', CustomJS(args = args, code = js_leave))

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
		'2ch4':'olive',
		'3ch4':'hotpink',
		'co':'darkred',
		'2co':'darkgray',
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
		'n2':'red',
		'o3':'red',
		'no':'pink',
		'nh3':'yellow',
		'hcn':'darkgray',
		'solar':'goldenrod',
		'other':'salmon',
		}

argu = sys.argv # commandline arguments

path='not_a_directory_1e7fwf8wrf78wf' #this should just be a non-existent directory name

if len(argu)>1:
	path = argu[1] # first commandline argument is the path to the spectra
else:
	print('You must give your /PATH/TO/SPECTRA as argument\n')
	sys.exit()

if os.path.isdir(path)==False:
	print('/!\\ You gave a wrong path /!\\\n') #error message if the given path doesn't exist
	if len(argu)>1:
		sys.exit()

save_path=os.path.join(path,'SAVE')
if not os.path.isdir(save_path):
	os.makedirs(save_path)

spectra = [i for i in os.listdir(path) if 'sfddaa' in i] # list of everything in the given directory

window_dic = {'y':'6339','w':'6073','z':'6220','s':'4852','k':'6500'} # dictionary to associate a spectrum prefix letter to a window center wavenumber.

TOOLS = "box_zoom,wheel_zoom,pan,undo,redo,reset,crosshair,save" #tools for bokeh figures, the first tool in the list will be active by default

select_spectra=sorted([i for i in spectra if ('.' in i)])

#############
# Main code #
#############

source_list = {}

window_tabs = []

tabs = []

save_figs = {}

# loop over the selected spectra
for spec_ID,spectrum in enumerate(select_spectra):

	if spectrum[0] not in save_figs.keys():
		first_check = True

	progress(select_spectra.index(spectrum),len(select_spectra),word=spectrum) #fancy loadbar

	#read the spectrum file
	spt_data = read_spt(os.path.join(path,spectrum))

	header = spt_data['header']
	SZA = str(spt_data['sza'])
	zobs = str(spt_data['zobs'])

	species = np.array([spt_data['columns'][var] for var in header])
	freq=species[0] # the frequency list
	tm=species[1]	# measured transmittance list
	tc=species[2]	# calculated transmittance list
	cont = species[3] # continuum
	not_gas = 4 # number of column that are not retrieved species
	residuals = spt_data['resid'] # 100*(calculated - measured)
	sigma_rms = spt_data['rms_resid'] # sqrt(mean(residuals**2))

	N_plots = list(range(len(species)-2)) # a range list from 0 to the number of plots, used by the checkbox group

	source_list[spectrum] = ColumnDataSource(data={var:spt_data['columns'][var] for var in header})
	source_list[spectrum].data['resid'] = residuals
	source_list[spectrum].data['const'] = np.zeros(len(residuals))

	if first_check:
		# spectrum figure 
		fig = figure(title=spectrum+'; SZA='+SZA+'°; zobs='+zobs+'km; %resid=100*(Measured-Calculated); RMSresid='+('%.4f' % sigma_rms)+'%',plot_width = 1000,plot_height=400,tools=TOOLS,toolbar_location=None,y_range=Range1d(-0.04,1.04),outline_line_alpha=0,active_inspect="crosshair",active_drag="box_zoom")
		# residual figure
		fig_resid = figure(plot_width=1000,plot_height=150,x_range=fig.x_range,tools=TOOLS,toolbar_location=None,y_range=Range1d(-1,1),active_inspect="crosshair",active_drag="box_zoom")
		
		save_figs[spectrum[0]] = [fig,fig_resid]
	else:
		fig = figure(title=spectrum+'; SZA='+SZA+'°; zobs='+zobs+'km; %resid=100*(Measured-Calculated); RMSresid='+('%.4f' % sigma_rms)+'%',plot_width = 1000,plot_height=400,tools=TOOLS,toolbar_location=None,y_range=save_figs[spectrum[0]][0].y_range,x_range=save_figs[spectrum[0]][0].x_range,outline_line_alpha=0,active_inspect="crosshair",active_drag="box_zoom")
		fig_resid = figure(plot_width=1000,plot_height=150,tools=TOOLS,toolbar_location=None,y_range=save_figs[spectrum[0]][1].y_range,x_range=save_figs[spectrum[0]][1].x_range,active_inspect="crosshair",active_drag="box_zoom")
	# axes labels
	fig_resid.xaxis.axis_label = u'Wavenumber (cm\u207B\u00B9)'
	fig_resid.yaxis.axis_label = '% Residuals'
	fig_resid.xaxis.axis_label_text_font_size = "12pt"
	fig_resid.xaxis.axis_label_text_font_size = "12pt"
	fig_resid.xaxis.major_label_text_font_size = "12pt"
	fig_resid.yaxis.axis_label_text_font_size = "12pt"
	fig_resid.yaxis.major_label_text_font_size = "12pt"
	fig.yaxis.axis_label = 'Transmittance'
	fig.yaxis.axis_label_text_font_size = "12pt"
	fig.yaxis.major_label_text_font_size = "12pt"
	fig.xaxis.major_label_text_font_size = "12pt"

	# group of checkboxes that will be used to toggle line and HoverTool visibility
	checkbox = CheckboxGroup(labels=header[not_gas:]+['Measured','Calculated'],active=N_plots,width=200)
	
	# plotting species lines
	plots = []
	for j,gas in enumerate(header[not_gas:]):
		try:
			plots.append(fig.line(x='Freq',y=gas,color=colors[gas],line_width=2,name=gas,source=source_list[spectrum]))
		except KeyError:
			print('KeyError:',gas,'is not specified in the "colors" dictionary, you need to add it with an associated color')
			sys.exit()
		# each line has a associated hovertool with a callback that looks at the checkboxes status for the tool visibility.
		#fig.add_tools( HoverTool(mode='vline',line_policy='prev',renderers=[plots[-1]],names=[gas],tooltips=OrderedDict( [('name',gas),('index','$index'),('(x;y)','($~x{0.00} ; @'+gas+'{0.000})')] ) ) )

	# adding the measured spectrum
	plots.append(fig.line(x='Freq',y='Tm',color='black',line_width=2,name='Tm',source=source_list[spectrum]))
	#fig.add_tools( HoverTool(mode='vline',line_policy='prev',renderers=[plots[-1]],names=['Tm'],tooltips=OrderedDict( [('name','Measured'),('index','$index'),('(x;y)','($~x{0.00} ; @Tm{0.000})')] )) )
	
	# adding the calculated spectrum
	plots.append(fig.line(x='Freq',y='Tc',color='chartreuse',line_width=2,name='Tc',source=source_list[spectrum]))
	#fig.add_tools( HoverTool(mode='vline',line_policy='prev',renderers=[plots[-1]],names=['Tc'],tooltips=OrderedDict( [('name','Calculated'),('index','$index'),('(x;y)','($~x{0.00} ; @Tc{0.000})')] )) )

	# legend outside of the figure
	fig_legend=Legend(items=[(header[j+not_gas],[plots[j]]) for j in range(len(species)-not_gas)]+[('Measured',[plots[-2]]),('Calculated',[plots[-1]])],location=(0,0),border_line_alpha=0)
	fig.add_layout(fig_legend,'right')
	fig.legend.click_policy = "hide"
	fig.legend.inactive_fill_alpha = 0.6

	# now the residual figure
	fig_resid.line(x='Freq',y='resid',color='black',name='residuals',source=source_list[spectrum])
	fig_resid.line(x='Freq',y='const',color='red',source=source_list[spectrum])
	#fig_resid.add_tools(HoverTool(mode='vline',line_policy='prev',names=['residuals'],tooltips={'index':'$index','(x;y)':'($~x{0.00} ; @resid{0.000})'}))

	# set up a dummy legend for the residual figure so that it aligns with the spectrum figure
	dummy = fig_resid.line(x=freq,y=[0 for i in range(len(residuals))],color='white',visible=False,alpha=0)
	fig_resid_legend=Legend(items=[('                 ',[dummy])],location=(0,0),border_line_alpha=0)
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

	# add vertical crosshair linked between the spectrum and residual figures
	add_vlinked_crosshairs(fig,fig_resid)

	# define the grid with the figures and widget box
	grid = gridplot([[fig,group],[fig_resid]],tools=TOOLS,toolbar_location='left')

	tabs.append( Panel(child=grid,title=spectrum.split('.')[-1]) )

	# put spectra in a window specific Tab if the next spectrum has a new prefix, or if the current spectrum is the last one
	try:
		if select_spectra[spec_ID+1][0]!=spectrum[0]:
			window_tabs.append( Panel(child=Tabs(tabs=tabs),title=window_dic[spectrum[0]]) )
			tabs = [] # reset the tabs list, will be filled by spectra from the next window
	except IndexError:
		window_tabs.append( Panel(child=Tabs(tabs=tabs),title=window_dic[spectrum[0]]) )
		
	first_check = False

final=Tabs(tabs=window_tabs,width=1100)

# write the HTML file
outfile=open(os.path.join(save_path,'comp_spectra_syn2.html'),'w')
outfile.write(file_html(final,CDN,'GFIT2 spectra'))
outfile.close()

print('\n')

sys.exit() # to make sure the program doesn't hang after it's finished
