### README ###

This app requires python 2.7.x (not tested with python 3.x) and bokeh 0.12.10 and netCDF4 installed.

	https://bokeh.pydata.org/en/latest/docs/installation.html

Public TCCON files can be downloaded from http://tccon.ornl.gov/

## How to use this app:

	- Put TCCON .eof.csv or .nc files in the 'data' folder of this app, only one type of file should be in the 'data' folder

	- Run the app with the command "bokeh serve --show tccon_app"

	- While the server is running, the app will be available in the browser at localhost:5006/tccon_app

## init.py ##

For private TCCON files there are 1200+ variables. I only make a subset of those available for selection in the variable dropdown widgets.
You can access more variables by editing the 'skip_list' list in the init.py program

This code can read from any .eof.csv files.

You just need to add the appropriate key:value pairs to the T_FULL and T_LOC dictionaries in init.py

# cache #

By default the program will cache data corresponding to specific sets of inputs to make loading faster.

This is not useful when using netcdf files, to disable the caching just set cache_max_size to 0 in init.py


## Contact ##

sebastien.roche@mail.utoronto.ca