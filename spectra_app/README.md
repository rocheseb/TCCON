# README #

Put the GGG output spectra (from the GGGPATH/spt folder) in the spectra_app/spectra folder

Run the app from a terminal, in the same directory as spectra_app run:

	bokeh serve --show spectra_app

A browser window will pop up with the app.

If you close the browser window, the server is still running and the app still available at http://localhost:5006/spectra_app

To shut the server down, use ctrl+C in the terminal

Each time a spectrum is loaded, a static .html file with the plot will be saved in spectra_app/save/
