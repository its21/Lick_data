import sys, urllib
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import pylab 
from astropy.io.votable import parse_single_table
import pywcs
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import angles	# pip install angles
import matplotlib.cm as cm
from scipy import ndimage

rrl = np.loadtxt('RRL_targets_phases.txt')
No_obj = len(rrl)
name_index = np.arange(0, No_obj, 1)

for i in range(79, No_obj):
	#clean all existing windows
	plt.clf()

	ra=rrl[i,1]
	dec=rrl[i,2]
	name= 'HAC'+str(i)

	# Construct URL to download DSS image cutout, and save to tmp.fits
	#change x = 3 and y = 3 (in &x=3&y=3&) if you want to change the field of view! 3x3 arcmin
	image_url = 'http://archive.eso.org/dss/dss/image?ra=%.5f&dec=%.5f&x=5&y=5&Sky-Survey=DSS1&mime-type=download-fits' % ((ra), (dec))
	print "Downloading DSS image..."
	urllib.urlretrieve(image_url, 'tmp.fits')

	# Download USNO-B1 catalog for the position
	#catalog_url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-B1&RA=%.5f&DEC=%.5f&SR=0.017&VERB=1' % ((ra), (dec))
	catalog_url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-B1&RA=%.5f&DEC=%.5f&SR=0.05&VERB=1' % ((ra), (dec))
	print "Downloading USNO-B1 catalog..."
	urllib.urlretrieve(catalog_url, 'tmp.cat')

	# Read RA, Dec and magnitude from XML format USNO catalog
	catalog = parse_single_table("tmp.cat")
	cat_ra = catalog.array['RA'] 
	cat_dec = catalog.array['DEC']
	cat_R1mag = catalog.array['R1']  
	cat_R2mag = catalog.array['R2']  

	# Make 4xN array, containing RA, Dec and mag for ref stars
	np_catalog = np.column_stack((cat_ra,cat_dec,cat_R1mag,cat_R2mag))
	#np.set_printoptions(threshold='nan')

	# Iterate over ref star catalog, rejecting stars with null photometry values.
	# Average R1 amd R2 mags
	# Only keep reference stars between 13th and 18th mag
	"""
	ref_catalog = []

	for row in np_catalog:
		if row[3] != 0 and row[2] != 0:
			mag = np.average((row[2],row[3]))
			if mag > 13 and mag < 16:		
				ref_catalog.extend((row[0], row[1], np.average((row[2],row[3]))))

	print "!!!!!!!!!!!"
	# Convert ref_catalog to numpy array, reshape it to proper dimensions
	np_ref_catalog = np.array(ref_catalog)
	np_ref_catalog = np.reshape(np_ref_catalog, (-1,3))
    
	print np_ref_catalog
	"""	
	# Shuffle ref_catalog and pick top three stars as references
	# Or change to sort on brightest stars
	#np.random.shuffle(np_ref_catalog)	
	#ref1= np_ref_catalog[0]
	#ref2= np_ref_catalog[1]
	#ref3= np_ref_catalog[2]

	dss_image = pyfits.open("tmp.fits")

	# Get pixel coordinates of SN, reference stars in DSS image
	wcs = pywcs.WCS(dss_image[0].header)
	#ref1_pix = wcs.wcs_sky2pix(np.array([ref1[0:2]], np.float_), 1)
	#ref2_pix = wcs.wcs_sky2pix(np.array([ref2[0:2]], np.float_), 1)
	#ref3_pix = wcs.wcs_sky2pix(np.array([ref3[0:2]], np.float_), 1)
	target_pix = wcs.wcs_sky2pix([(np.array([ra,dec], np.float_))], 1)

	print target_pix[0,0]
	#print ref1_pix
	#print ref2_pix
	#print ref3_pix

	# Plot finder chart
	#plt.set_cmap('gray_r')
	#if you need to rotate the finding chart by 18deg to the east
	dss_image_rot = ndimage.rotate(dss_image[0].data, 18)
	deg_field = 5/60.
	extent = [ra - deg_field/2. , ra + deg_field/2., dec-deg_field/2., dec+ deg_field/2. ]
	plt.imshow(np.sqrt(dss_image_rot),origin='lower', extent = extent, aspect = 'auto', cmap = cm.binary_r, interpolation = 'nearest') #, vmax=0.9*np.max(dss_image_rot)

	# Mark target
	#plt.scatter(ra,dec, marker = 'D',color = 'red')
	plt.plot([ra+deg_field/10.,ra+deg_field/40.],[dec,dec], 'k-', lw=2, color = 'red')
	plt.plot([ra,ra],[dec+deg_field/40.,dec+deg_field/10.], 'k-', lw=2, color = 'red')
	#plt.annotate(name, xy=(target_pix[0,0], target_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points')

	#my_polygon = [target_pix[0,0],(target_pix[0,0])],[target_pix[0,1],(target_pix[0,1])]
	#my_polygon2 = rotatePolygon(my_polygon,18)
	#plt.scatter(my_polygon2[0], marker = 'D',color = 'green')

	#plt.scatter([target_pix[0,0],(target_pix[0,0])],[target_pix[0,1],(target_pix[0,1])], marker = 'D',color = 'red')
	#plt.plot([target_pix[0,0]+10,(target_pix[0,0]+3)],[target_pix[0,1],(target_pix[0,1])], 'k-', lw=2, color = 'red')
	#plt.plot([target_pix[0,0],(target_pix[0,0])],[target_pix[0,1]+3,(target_pix[0,1])+10], 'k-', lw=2, color = 'red')
	#plt.annotate(name, xy=(target_pix[0,0], target_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points')

	# Mark and label reference stars
	"""
	plt.plot([ref1_pix[0,0]+10,(ref1_pix[0,0]+3)],[ref1_pix[0,1],(ref1_pix[0,1])], 'k-', lw=2, color = 'blue')
	plt.plot([ref1_pix[0,0],(ref1_pix[0,0])],[ref1_pix[0,1]+3,(ref1_pix[0,1])+10], 'k-', lw=2, color = 'blue')
	plt.annotate("R1", xy=(ref1_pix[0,0], ref1_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points', color = 'blue')
	plt.plot([ref2_pix[0,0]+10,(ref2_pix[0,0]+3)],[ref2_pix[0,1],(ref2_pix[0,1])], 'k-', lw=2, color = 'blue')
	plt.plot([ref2_pix[0,0],(ref2_pix[0,0])],[ref2_pix[0,1]+3,(ref2_pix[0,1])+10], 'k-', lw=2, color = 'blue')
	plt.annotate("R2", xy=(ref2_pix[0,0], ref2_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points', color = 'blue')
	#plt.plot([ref3_pix[0,0]+10,(ref3_pix[0,0]+3)],[ref3_pix[0,1],(ref3_pix[0,1])], 'k-', lw=2)
	#plt.plot([ref3_pix[0,0],(ref3_pix[0,0])],[ref3_pix[0,1]+3,(ref3_pix[0,1])+10], 'k-', lw=2)
	#plt.annotate("R3", xy=(ref3_pix[0,0], ref3_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points')
	"""
	# Set limits to size of DSS image
	pylab.xlim([extent[0]+1/60.,extent[1]-1/60.])
	pylab.ylim([extent[2]+1/60.,extent[3]-1/60.])

	# Plot compass
	plt.plot([(dss_image[0].data.shape[0])-10,(dss_image[0].data.shape[0]-40)],[10,10], 'k-', lw=2)
	plt.plot([(dss_image[0].data.shape[0])-10,(dss_image[0].data.shape[0])-10],[10,40], 'k-', lw=2)
	plt.annotate("N", xy=((dss_image[0].data.shape[0])-10, 40),  xycoords='data',xytext=(-4,5), textcoords='offset points')
	plt.annotate("E", xy=((dss_image[0].data.shape[0])-40, 10),  xycoords='data',xytext=(-12,-5), textcoords='offset points')

	# Set axis tics (not implemented correctly yet)
	plt.tick_params(labelbottom='off')
	plt.tick_params(labelleft='off')
	plt.axes().xaxis.set_major_locator(MaxNLocator(3))
	plt.axes().yaxis.set_major_locator(MaxNLocator(3))
	plt.axes().set_xlabel('3\'')
	plt.axes().set_ylabel('3\'')

	# Set size of window (leaving space to right for ref star coords)
	plt.subplots_adjust(right=0.8,left=0.05)

	# List coords, mag of references etc
	plt.text(1.02, 0.80,str(name), transform=plt.axes().transAxes)
	plt.text(1.02, 0.75,str(str(ra)+" "+str(dec)),transform=plt.axes().transAxes)
	"""
	plt.text(1.02, 0.65,str('R1, mag='+str(ref1[2])), transform=plt.axes().transAxes)
	plt.text(1.02, 0.6,str(str(ref1[0])+" "+str(ref1[1])),transform=plt.axes().transAxes)
	plt.text(1.02, 0.5,str('R2, mag='+str(ref2[2])), transform=plt.axes().transAxes)
	plt.text(1.02, 0.45,str(str(ref2[0])+" "+str(ref2[1])),transform=plt.axes().transAxes)
	"""
	# Save to pdf
	pylab.savefig(str(name+'_finder.png'))
	#plt.show()
