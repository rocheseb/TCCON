def hcl_cells():

	cell = {}

	template = {
				'ID':None,
				'length':None,
				'batch':None,
				'owner':None,
				'text_on_cell_body':None,
				'effp_h35cl_296k':None,
				'effp_h37cl_296k':None,
				'calibration_run':None,
				'h35cl_column':None,
				'h37cl_column':None,
				'location':None,
				}
	"""
	# copy paste this template to add a cell
	# site
	cell['site'] = {key:template[key] for key in template}
	cell['site']['location'] =
	cell['site']['owner'] = 
	cell['site']['ID'] =
	cell['site']['length'] =
	cell['site']['batch'] =
	cell['site']['text_on_cell_body'] =
	cell['site']['effp_h35cl_296k'] = 
	cell['site']['effp_h37cl_296k'] = 
	cell['site']['h35cl_column'] = 
	cell['site']['h37cl_column'] = 
	cell['site']['calibration_run'] = 
	"""
	
	# Eureka
	cell['eu'] = {key:template[key] for key in template}
	cell['eu']['location'] = 'Eureka'
	cell['eu']['owner'] = 'University of Toronto'
	cell['eu']['ID'] = 30
	cell['eu']['length'] = 100
	cell['eu']['batch'] = 'CT1303'
	cell['eu']['text_on_cell_body'] = "cell #28"
	cell['eu']['effp_h35cl_296k'] = '4.78'
	cell['eu']['effp_h37cl_296k'] = '4.817'
	cell['eu']['h35cl_column'] = '1.2909e+22'
	cell['eu']['h37cl_column'] = '1.2836e+22'
	cell['eu']['calibration_run'] = 1312

	# East Trout Lake
	cell['et'] = {key:template[key] for key in template}
	cell['et']['location'] = 'East Trout Lake'
	cell['et']['owner'] = 'University of Toronto'
	cell['et']['ID'] = 62
	cell['et']['length'] = 100
	cell['et']['batch'] = 'CT1510'
	cell['et']['text_on_cell_body'] = 'cell #56'
	cell['et']['effp_h35cl_296k'] = '4.673'
	cell['et']['effp_h37cl_296k'] = '4.683'
	cell['et']['h35cl_column'] = '1.3310e+22'
	cell['et']['h37cl_column'] = '1.3260e+22'
	cell['et']['calibration_run'] = 1510

	return cell

def hbr_cells():

	cell = {}

	template = {
				'ID':None,
				'length':None,
				'batch':None,
				'owner':None,
				'text_on_cell_body':None,
				'calibration_run':None,
				'location':None,
				'pressure':None,
				'column':None,
				}

	# Eureka
	cell['eu'] = template
	cell['eu']['location'] = 'Eureka'
	cell['eu']['pressure'] = '1.52'
	cell['eu']['column'] = '7.60e20'

	return cell