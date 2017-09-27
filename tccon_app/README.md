This app requires python 2.7.x (not tested with python 3.x) with bokeh 0.12.9 and netCDF4 installed.

	- Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html
	- netCDF4: http://unidata.github.io/netcdf4-python/

Public TCCON files can be downloaded from http://tccon.ornl.gov/

How to use this app:

	- Put TCCON .eof.csv or .nc files in the 'data' folder, only one type of file should be in the 'data' folder

	- Run the app with the command "bokeh serve --show tccon_app"

	- While the server is running, the app will be available in the browser at localhost:5006/tccon_app

For private TCCON files there are 1200+ variables. I only make a subset of those available for selection in the variable dropdown widgets.
You can access more variables by editing the 'skip_list' list in init.py

This code can read from any .eof.csv files. So it can also plot data from an EM27. You just need to add the appropriate key:value pairs to the T_FULL and T_LOC dictionaries in init.py
