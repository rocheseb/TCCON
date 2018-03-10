# README #

This app requires python 2.7.x (not tested with python 3.x) with bokeh 0.12.10 (not yet tested with latest versions) and netCDF4 installed.
It also uses the parse package.

	Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html
	netCDF4: http://unidata.github.io/netcdf4-python/

### How to use this app ###

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
- In lft_app/spectra/cut/temp write the spectrum filename, scanner temperature, and aperture size
	
		spectrumfilename1,temperature1,apt_size1
		spectrumfilename2,temperature2,apt_size2
		etc.
	
- In lft_app/cell_data.py, add your cell information (follow the template)

- Run the app from the linefit/lft145/ directory with the command

	bokeh serve --show lft_app

While the server is running, the app will be available in the browser at localhost:5006/lft_app

- Python dictionaries of the data are saved in lft_app/saved_sessions/
- PDF documents with all the plots are saved in lft_app/pdf/

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

### Contact ###

sebastien.roche@mail.utoronto.ca
