# README #

https://bitbucket.org/rocheseb/tccon

Various python codes to plot data from TCCON files.

All codes are written to run with Python 2.7 and bokeh 0.12.10

libraries:

	- netCDF4: http://unidata.github.io/netcdf4-python/
	
	- bokeh: http://bokeh.pydata.org/en/latest/docs/installation.html
	
	- NumPy: http://www.numpy.org/


### What is this repository for? ###

This repository will contain some programs that can be used to visualize some of the outputs of GGG, the retrieval algorithm of TCCON (http://www.tccon.caltech.edu/)

### How do I get set up? ###

Each code contains a description section at the top that explains how to use it.

### tccon_app ###

install bokeh

put tccon .eof or .nc files in tccon_app/data

To run the app, in the same directory as tccon_app, run this command from a terminal

	bokeh serve --show tccon_app

Your browser should pop up with the plots. The app will be available at http://localhost:5006/tccon_app

### Who do I talk to? ###

* Repo owner or admin
sebastien.roche@mail.utoronto.ca
