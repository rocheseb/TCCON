import os, sys

import readoutputs as grd

from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.io import show

GGGPATH = os.environ['GGGPATH']

path1 = os.path.join(GGGPATH,'models','gnd','ncep')
path2 = os.path.join(GGGPATH,'models','gnd','ncepy')

file_list = os.listdir(path1)

fig_dict = {}
for mod_file in file_list:

	data1 = grd.read_mod(os.path.join(path1,mod_file))
	data2 = grd.read_mod(os.path.join(path2,mod_file))

	for var in ['H2O','RH','Temperature','Height']:
		fig_dict[var] = []

		dif = (data2[var]-data1[var])/data1[var]

		fig = figure(plot_width=125,plot_height=175)

		fig.line(dif,data1['Pressure'])
		fig.scatter(dif,data1['Pressure'])

		fig_dict[var] += [fig]

for var in fig_dict:
	grid = gridplot(fig_dict[var],ncols=7)
	show(grid)
	a=raw_input()