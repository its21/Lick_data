import re

def string2hours(string):
	""" 
	Converts a Sexagesimal string to a decimal hour
	
	Parameters
	----------
	string : str
		A string representing a sexagismal time, in the form HH:MM:SS.SS.
		
	Returns
	-------
	decHours : float
		The converted sexagismal time to decimal time.
	
	"""
	string = string.strip() # trim leading/trailing whitespace

	div = '[|:|\s]' # allowable field delimiters "|", ":", whitespace
	
	# (optional +/-) degrees + min + decimal seconds
	sdec = '([+-]?)(\d{1,3})' + div + '(\d{1,2})' + div + '(\d{1,2}\.?\d+?)'
	
	co_re= re.compile(sdec)
	co_search= co_re.search(string)
	
	if co_search is None:
		raise ValueError("Invalid input string: '%s'" % string)
	elems = co_search.groups()

	plusminus = elems[0]
	hours = float(elems[1])
	minutes = float(elems[2])
	seconds = float(elems[3])
	
	# Check for nonsense values
	if hours > 24.0:
		raise ValueError("Hour value must be < 24.")
	if minutes >= 60.0:
		raise ValueError("Minute value must be < 60.")
	if seconds >= 60.0:
		raise ValueError("Second value must be < 60.")
		
	# Convert dec
	decHours = hours + minutes/60.0 + seconds/3600.0
	
	if plusminus is "-":
		decHours = -1.0 * decHours
	
	return decHours

def string2deg(string):
	""" 
	Converts a Sexagesimal string to a decimal degrees
	
	Parameters
	----------
	string : str
		A string representing a sexagismal time, in the form [+|-]DD:MM:SS.SS
		
	Returns
	-------
	decHours : float
		The converted sexagismal degrees to decimal degrees
	
	"""
	string = string.strip() # trim leading/trailing whitespace

	div = '[|:|\s]' # allowable field delimiters "|", ":", whitespace
	
	# (optional +/-) degrees + min + decimal seconds
	sdec = '([+-]?)(\d{1,3})' + div + '(\d{1,2})' + div + '(\d{1,2}\.?\d+?)'
	
	co_re= re.compile(sdec)
	co_search= co_re.search(string)
	
	if co_search is None:
		raise ValueError("Invalid input string: %s" % string)
	elems = co_search.groups()

	plusminus = elems[0]
	degrees = float(elems[1])
	arcminutes = float(elems[2])
	arcseconds = float(elems[3])
	
	# Check for nonsense values
	if degrees > 90.0:
		raise ValueError("Degree value must be <= 90.")
	if arcminutes >= 60.0:
		raise ValueError("Arcminute value must be < 60.")
	if arcseconds >= 60.0:
		raise ValueError("Arcsecond value must be < 60 (was %f)." % arcseconds)
		
	# Convert dec
	decDegrees = degrees + arcminutes/60.0 + arcseconds/3600.0
	
	if plusminus is "-":
		decDegrees = -1.0 * decDegrees
	
	return decDegrees

	