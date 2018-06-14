# README #

This app requires python 2.7.x (not tested with python 3.x) with bokeh and netCDF4 installed.
It also uses the parse package.

	Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html
	netCDF4: http://unidata.github.io/netcdf4-python/

### Python ###

I suggest downloading python from https://www.anaconda.com/download/
Choose your operating system and get Python2.7

To install bokeh use the command (with windows you need to run the terminal as administator):

	conda install -c bokeh bokeh

To install netCDF4:

	conda install netCDF4

If the "conda" command does not find the package, you can also use pip:

	pip install PackageName

If you encounter error messages related to bokeh when running the app, you can try to revert to an earlier version of the package with:

	conda install bokeh=0.12.10

### How to use this app ###

This does not yet read OPUS files.

.dpt (data point table) files can be generated in OPUS via the pop-up window generated with "Save as" 

- Put the lft_app folder in the linefit/lft145/ directory
- Spectrum file names need to follow this naming convention: YYMMDD_site_cell_X_MOPD_num.dpt
	- YYMMDD year month day
	- site: two letter site abbreviation (e.g. eu for Eureka, oc for Lamont)
	- cell: one of 'hbr', 'n2o', 'hcl'
	- X: 'v' for vented instrument, 'e' for evacuated
	- MOPD: the maximum optical path difference in cm
	- num: an index number for the cell test (there might be more than one per day)
	
		e.g. 180308_eu_HCl_45_e_0 for the first HCl cell test with 45 MOPD in an evacuated instrument at Eureka on March 8 2018

- For several tests in one day : 161122_eu_HCl_45_e_0.dpt, 161122_eu_HCl_45_e_1.dpt etc.
- Spectra in .DPT format, with no headers, must be placed in lft_app/spectra/cut (cut the spectra, e.g. between ~5200-5900 wavenumbers for HCl cells)
- In lft_app/spectra/cut/temp.dat write the spectrum filename, scanner temperature, and aperture size. A 4th optional parameter can be added if the spectrum has a significant spectral detuning (see DISCLAIMER below).
	
		spectrumfilename1,temperature1,apt_size1
		spectrumfilename2,temperature2,apt_size2,spectral_detuning2
		etc.
	
- In lft_app/lft_setup.py, add your cell information and the Focal Length of Collimator of your instrument (follow the template)

- To run the app, navigate to the linefit/lft145/ directory in your terminal and use the command

	bokeh serve --show lft_app

The --show option will pop up the browser.

While the server is running, the app will be available in the browser at localhost:5006/lft_app

- Python dictionaries of the data are saved in lft_app/saved_sessions/
- PDF documents with the plots are saved in lft_app/pdf/

By default the spectrum itself will be plotted in the browser and also saved in the data dictionary.
This can lead to very large files and more loading time. To avoid this the app can be run in light mode with:

	bokeh serve --show lft_app --args light

There are two example spectra from Eureka in lft_app/spectra/cut/

### Rationg of spectra ###

Spectra should be ratioed to ~1 to be used with the linefit extended mode:

- HCl cells: 
	- no background
	- I fit a 2nd order polynomial to the spectrum without the lines and use that to ratio the spectrum to normalize it to ~1 (seems more consistent than using a fixed numbers)

- HBr cells:
	- background
	- the background file should be cut the same way as the spectrum, have the same file name but starting with 'ref_' (e.g. ref_180308_eu_HBr_180_e_0.dpt)
	- put the HBr background files in lft_app/spectra/background/
	- the spectra are ratioed with the background
	- the resulting ratioed spectrum is ratioed with its own average to normalize it to ~1

- N2O cells:
	- background, but different resolution from the spectrum
	- the rationg of spectrum with background is done in OPUS
	- the resulting spectrum should be placed in lft_app/spectra/cut
	- it will be ratioed with its own average to normalize it to ~ 1

### Other info ###

N2O and HBr cell spectra are processed in a loop until the cell pressure converges; this usually take 2-3 linefit runs.

The python dictionaries saved in lft_app/saved_sessions/ can be merged with a utility program lft_app/utils/merge_sessions.py

The merged file can then be loaded from the browser.

### DISCLAIMER ###

DISCLAIMER: if any warning or error message is given by linefit, this app will hang, you should then run linefit from the terminal to figure out what the problem is

The app may hang if there is any convergence problem, or if a significant spectral detuning is detected.

There will be more detailed outputs in the terminal than in the browser.

If a significant spectral detuning is detected. Run linefit from the terminal, notice the value of spectral residuals given after the warning, and add it to the temp.dat line of the spectrum:

	spectrumfilename1,temperature1,apt_size1,spectral_detuning1

### Contact ###

sebastien.roche@mail.utoronto.ca
