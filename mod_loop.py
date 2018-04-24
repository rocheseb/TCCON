"""
small code to run mod_maker multiple times
"""

from __future__ import print_function # allows the use of Python 3.x print function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

import subprocess

def execute(cmd):
	'''
	function to execute a prompt command and print the output as it is produced
	'''
	popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
	for stdout_line in iter(popen.stdout.readline, ""):
		yield stdout_line
	popen.stdout.close()
	return_code = popen.wait()
	if return_code:
		raise subprocess.CalledProcessError(return_code, cmd)

# dictionary mapping TCCON site abbreviations to their lat-lon-alt data, and full names
site_dict = {
			'pa':{'name': 'Park Falls','loc':'Wisconsin, USA','lat':45.945,'lon':269.727,'alt':442},
			'oc':{'name': 'Lamont','loc':'Oklahoma, USA','lat':36.604,'lon':262.514,'alt':320},
			'wg':{'name': 'Wollongong','loc':'Australia','lat':-34.406,'lon':150.879,'alt':30},
			'db':{'name': 'Darwin','loc':'Australia',},
			'or':{'name': 'Orleans','loc':'France','lat':47.97,'lon':2.113,'alt':130},
			'bi':{'name': 'Bialystok','loc':'Poland','lat':53.23,'lon':23.025,'alt':180},
			'br':{'name': 'Bremen','loc':'Germany','lat':53.1037,'lon':8.849517,'alt':30},
			'jc':{'name': 'JPL 01','loc':'California, USA','lat':34.202,'lon':241.825,'alt':390},
			'jf':{'name': 'JPL 02','loc':'California, USA','lat':34.202,'lon':241.825,'alt':390},
			'ra':{'name': 'Reunion Island','loc':'France','lat':-20.901,'lon':55.485,'alt':87},
			'gm':{'name': 'Garmisch','loc':'Germany','lat':47.476,'lon':11.063,'alt':743},
			'lh':{'name': 'Lauder 01','loc':'New Zealand','lat':-45.038,'lon':169.684,'alt':370},
			'll':{'name': 'Lauder 02','loc':'New Zealand','lat':-45.038,'lon':169.684,'alt':370},
			'tk':{'name': 'Tsukuba 02','loc':'Japan','lat':63.0513,'lon':140.1215,'alt':31},
			'ka':{'name': 'Karlsruhe','loc':'Germany','lat':49.1002,'lon':8.4385,'alt':119},
			'ae':{'name': 'Ascenssion Island','loc':'United Kingdom','lat':-7.933333,'lon':345.583333,'alt':0},
			'eu':{'name': 'Eureka','loc':'Canada','lat':80.05,'lon':273.58,'alt':610},
			'so':{'name': 'Sodankyla','loc':'Finland','lat':67.3668,'lon':26.6310,'alt':188},
			'iz':{'name': 'Izana','loc':'Spain','lat':28,'lon':344,'alt':0},
			'if':{'name': 'Idianapolis','loc':'Indiana, USA','lat':39.861389,'lon':273.996389,'alt':270},
			'df':{'name': 'Dryden','loc':'California, USA','lat':34959917,'lon':242.118931,'alt':700},
			'js':{'name': 'Saga','loc':'Japan','lat':33.240962,'lon':130.288239,'alt':7},
			'fc':{'name': 'Four Corners','loc':'USA','lat':36.79749,'lon':251.51991,'alt':1643},
			#'ci':{'name': 'Pasadena','loc':'California, USA','lat':34.13623,'lon':241.873103,'alt':230},
			'ci':{'name': 'Pasadena','loc':'California, USA','lat':34.136,'lon':241.873,'alt':230},
			'rj':{'name': 'Rikubetsu','loc':'Japan','lat':43.4567,'lon':143.7661,'alt':380},
			'pr':{'name': 'Paris','loc':'France','lat':48.846,'lon':2.356,'alt':60},
			'ma':{'name': 'Manaus','loc':'Brazil','lat':-3.2133,'lon':299.4017,'alt':50},
			'sp':{'name': 'Ny-Alesund','loc':'Norway','lat':78.92324,'lon':11.92298,'alt':0},
			'et':{'name': 'East Trout Lake','loc':'Canada','lat':54.353738,'lon':255.013333,'alt':501.8},
			'an':{'name': 'Anmyeondo','loc':'Korea','lat':36.5382,'lon':126.331,'alt':30},
			'bu':{'name': 'Burgos','loc':'Philippines','lat':18.5325,'lon':120.6496,'alt':35},
			'we':{'name': 'Jena','loc':'Austria','lat':50.91,'lon':11.57,'alt':211.6},
			}

for site in site_dict:
	print('\n\nNOW DOING:',site_dict[site]['name'])
	for mode in ['ncep','merraglob','fpglob','fpitglob']:
		for line in execute(['python','mod_maker.py',site,'20171210-20171217',mode]):
			print(line,end="")