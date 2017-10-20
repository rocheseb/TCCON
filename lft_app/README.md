# README #

This app requires python 2.7.x (not tested with python 3.x) with bokeh 0.12.10 and netCDF4 installed.

	Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html
	netCDF4: http://unidata.github.io/netcdf4-python/

### How to use this app ###

- Put the lft_app folder in the linefit/lft145/ directory
- Spectrum name e.g. HCl_45_eu_161122.dpt for a HCl test from Eureka on November 22nd 2016, with 45 cm max OPD
- For several tests in one day : HCl_45_1_eu_161122.dpt, HCl_45_2_eu_161122.dpt etc.
- Spectra in .DPT format, with no headers, must be placed in lft_app/spectra/ (cut the spectra between ~5200-5900 wavenumbers)
- In lft_app/spectra/temp write the spectrum filename and scanner temperature

- Run the app from the linefit/lft145/ directory with the command

	bokeh serve --show lft_app

While the server is running, the app will be available in the browser at localhost:5006/lft_app

- Python dictionaries of the data are saved in lft_app/saved_sessions/
- PDF document with all the plots are saved in lft_app/

### DISCLAIMER ###

Always run linefit itself first, this app hides the linefit terminal output and will not work if linefit outputs error messages or requires user input.
It should be used as a visualization tool after linefit has been run from the terminal without errors.

lft_app/lft14_hcl.inp is the template used to rewrite the actual linefit input file

Some things still need to be edited manually in the template e.g.

	- maximum inclination of rays in the interferometer when changing the aperture size

Check the 'modify_input_file' function in 'main.py' to see what is edited programatically and what is not

If the spectra in lft_app/spectra/ are not ratioed to ~1, they will be ratioed using a 2nd order polynomial fit (see ratio_spectrum function in 'main.py')

### Contact ###

sebastien.roche@mail.utoronto.ca
