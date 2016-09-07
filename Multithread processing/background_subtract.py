def load_im(image,array):
	jj=0
	hdulist = pyfits.open(image)
	for hdu in hdulist:
	  if hdu.header['NAXIS'] != 0:
	    im_dim=(4094, 2096) #hdu.shape
	    im_ndim=len(im_dim)
	    if im_ndim==2:
		  array[0,jj,:,:]=hdu.data
		  jj+=1
	hdulist.close()
	return array

def weight_im(weight_image,image):
	for chip_num in range(im_nchip):
		image[0][chip_num] = image[0][chip_num]*weight_image[0][chip_num]/np.max(weight_image[0][chip_num])
	return image

def find_nearest(array,nvals,val):
	sorted = np.sort(np.abs(array))
	keep_vals = sorted[0:nvals]
	inds = []
	for value in keep_vals:
		inds.append(list(np.abs(array)).index(value))
	inds = np.array(inds)
	return inds

def median_subtract(back_cube,med_data,image,chip_num,q):
	print "Calculating medians..."
	for ii in range(im_size[0]):
		for jj in range(im_size[1]):
			med_data[ii][jj] = stats.nanmedian([back_cube[kk][chip_num][ii][jj] for kk in range(len(back_cube))])

	print "Subtracting background..."
	image = image-med_data
	q.put([image, med_data])

def median_global_subtract(back_cube,med_data,image,chip_num,q):
	med_data[:][:] = stats.nanmedian([back_cube[kk][chip_num] for kk in range(len(back_cube))])

	bv_mask = np.isnan(med_data)
	chip_disp = stats.nanstd(med_data,axis=None)
	chip_med = stats.nanmedian(med_data,axis=None)
	blah = chip_med+np.random.randn(im_size[0],im_size[1])*chip_disp
	
	blah[~bv_mask] = 0.
	med_data[bv_mask] = blah[bv_mask]

	image = image-med_data
	q.put([image, med_data, chip_num])
	
def tps_subtract(back_cube,image,chip_num):
	back_time = time.time()
	zz = np.zeros( (im_size[0], im_size[1]) )
	med_size=320

	z_chk = []
	for ii in range(len(back_im_cube)):
		z_chk.append([])

	nx, ny = 10, 10
	xx1, yy1 = np.meshgrid(np.arange(0, im_size[0], nx), 
						 np.arange(im_size[1], 0, -1*ny))

	xx2, yy2 = np.meshgrid(np.arange(0, im_size[0], 1), 
						np.arange(im_size[1], 0, -1))	

	cnt = 0
	for kk in range(len(back_im_cube)):
		xx = 0
		x_chk = []
		y_chk = []
		while xx < im_size[0]:
			x_chk_size = med_size
			if xx+x_chk_size > im_size[0]: x_chk_size = im_size[0]-xx
			yy = 0
			while yy < im_size[1]:
				y_chk_size=med_size
				if yy+y_chk_size > im_size[1]: y_chk_size = im_size[1]-yy
				x_chk.append(int(xx+x_chk_size/2))
				y_chk.append(int(yy+y_chk_size/2))
				z_chk[cnt].append(stats.nanmedian(back_im_cube[kk][chip_num][xx:xx+x_chk_size,yy:yy+y_chk_size],axis=None))
				yy+=y_chk_size
			xx+=x_chk_size
		cnt+=1

	x_chk = np.array(x_chk)
	y_chk = np.array(y_chk)
	z_chk = np.array(z_chk)

	z_fit = []
	for ii in range(len(z_chk[0])):#range(len(z_chk[0])):
# 		print np.median([z_chk[jj][ii] for jj in range(len(z_chk))])
		z_fit.append(np.median([z_chk[jj][ii] for jj in range(len(z_chk))]))

	z_fit = np.array(z_fit)

# 	nx, ny = 10, 10
# 	xx, yy = np.meshgrid(np.arange(0, im_size[0], nx), 
# 						 np.arange(im_size[1], 0, ny))

	interp_type = 'thin_plate'
	print "Interpolating..."
	tps = interpolate.Rbf(x_chk,y_chk,z_fit,function=interp_type)#'thin_plate')
	print "Interpolation done..."

	zz_i = tps(xx1,yy1)

# 	xx, yy = np.meshgrid(np.arange(0, im_size[0], 1), 
# 						 np.arange(im_size[1], 0, 1))

#	print tps_data
#  	for ii in range(im_size[0]):
#  		for jj in range(im_size[1]):
#  			zz[ii][jj] = tps(ii,jj).T
	zz = tps(xx2,yy2).T

	print "Plotting..."
	plt.figure()
	plt.subplot(211)
 	plt.imshow(zz_i, extent=(0, im_size[0], 0, im_size[1]))
# 	plt.imshow(tps_data, extent=(0, im_size[0], 0, im_size[1]))
	plt.colorbar()
	plt.subplot(212)
	plt.scatter(x_chk, y_chk, c=z_fit)
	plt.colorbar()
	plt.xlim(0,im_size[0])
	plt.ylim(0,im_size[1])
	plt.savefig('/Users/matt/Desktop/backs/chip%s_back%i_d%i_fit.pdf' % (chip_names[chip_num],med_size,int(do_dither)))

	plt.figure()
 	plt.imshow(zz, extent=(0, im_size[0], 0, im_size[1]))
	plt.colorbar()
# 	plt.xlim(0,im_size[0])
# 	plt.ylim(0,im_size[1])
	plt.savefig('/Users/matt/Desktop/backs/chip%s_back%i_d%i_%s.pdf' % (chip_names[chip_num],med_size,int(do_dither),interp_type))

	image = image - zz
	print "Time for Chip#%i background subtraction = %f seconds." % (chip_num+1, time.time()-back_time)

	return image, zz
# 	q.put([image, tps_data])

def infill(back_im_cube):
	for ii in range(len(back_im_cube)):
#		print "Filling NaNs, background frame %i/%i" % (ii+1,len(back_im_cube))
		mask_im_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
		mask_im_data = load_im(mask_ims[ii],mask_im_temp)
		for jj in range(im_nchip):
			print "Filling NaNs, background frame %i/%i, Chip#%i/%i" % (ii+1,len(back_im_cube),jj+1,im_nchip)
			bv_mask = (mask_im_data[0][jj] == 1.)
#			bv_mask = (mask_im_data[0][jj] > 0.)

			chip_disp = stats.nanstd(back_im_cube[ii][jj],axis=None)

			blah = np.random.rand(im_size[0],im_size[1])
			blah = (blah-0.5)*2.
			blah = stats.nanmedian(back_im_cube[ii][jj],axis=None)+blah*chip_disp

			blah[~bv_mask] = 0.
			back_im_cube[ii][jj][bv_mask] = blah[bv_mask]

	return back_im_cube

import pyfits
import numpy as np
import glob
from multiprocessing import Process, Queue
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import interpolate, stats
import itertools
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os,sys,subprocess
import time
#from random import randint
import random
import inpaint

start_time = time.time()

if len(sys.argv) < 2 :
	print "ERROR - missing subtraction method"
	exit()
# if sys.argv[1] != "median":
# 	if sys.argv[1] != "median global":
# 		if sys.argv[1] != "bispline":
# 			print "ERROR - background subtraction method must be one of [median, median global, bispline]"
# 			exit()

if sys.argv[1] not in ["median", "median global", "bispline", "dynamic tps"]:
	print "ERROR - background subtraction method must be one of [median, median global, bispline, tps]"
	exit()

im_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
#im_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
do_tile = '1'
do_dither = '1'
do_filt = 'i'
n_backs = 3

# im_dir = sys.argv[2]
# do_tile = sys.argv[3]
# do_dither = sys.argv[4]
# do_filt = sys.argv[5]
# n_backs = int(sys.argv[6])
# ss_im_out = sys.argv[7]
# sky_im_out = sys.argv[8]

#back_dir = '/Volumes/MyPassport/masks/'
#mask_dir = '/Volumes/MyPassport/masks/'
back_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
mask_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
#back_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
#mask_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'

# back_dir = im_dir
# mask_dir = im_dir

test_image = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.fits'
test_weight = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.WEIGHT.fits'
test_mask = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.MASK.fits'

im_h = pyfits.open(test_image)[0].header
chip_names = []

#Open image file, and scale it by the weight map. ???SHOULD THE IMAGE BE SCALED????
print "\nProcessing image %s for background subtraction..." % test_image
print "\nLoading image..."
hdulist = pyfits.open(test_image)

im_nchip=0
for hdu in hdulist:
	if hdu.header['NAXIS'] != 0:
		im_dim = np.array(hdu.data).shape#(4096,2048)#(4094, 2046) #hdu.shape
		im_ndim=len(im_dim)
		chip_names.append(hdu.header['DETPOS'])
		if im_ndim==2:
			if im_nchip==0: im_size=im_dim
			im_nchip+=1
hdulist.close()

im_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
#w_im_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
im_mask_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
im_masked = np.zeros((1,im_nchip,im_size[0],im_size[1]))
im_data = load_im(test_image,im_data)
im_masked[:] = im_data[:]
#im_masked = load_im(test_image,im_masked)
#w_im_data =load_im(test_weight,w_im_data)
im_mask_data = load_im(test_mask,im_mask_data)

#Weight scale the original image
# im_data = weight_im(w_im_data,im_data)

#Create masked image data
print "Creating masked image..."
for ii in range(len(im_data[0])):
 	bv_mask = (im_mask_data[0][ii] == 1.)
	im_masked[0][ii][bv_mask] = np.nan

#del w_im_data

#hdulist_back_out = hdulist
#for ii in range(im_nchip):
#      hdulist_back_out[ii+1].data = back_im_cube[0][ii]
#hdulist_back_out.writeto('/Users/matt/Desktop/deleteme_image.fits',clobber=True)
#hdulist_back_out.close()
#exit()

####################################################################
####Open background and mask files and scale them by weight maps####
####################################################################

#Determine the background files closest in time to the image
# temp_back_ims = glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.fits')
# temp_back_weights = glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.WEIGHT.fits')
# temp_mask_ims = glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.MASK.fits')
if do_dither == '1':
	temp_back_ims = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.fits')))
	temp_back_weights = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.WEIGHT.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.WEIGHT.fits')))
 	temp_mask_ims = np.concatenate((glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.MASK.fits'),glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.MASK.fits')))
#	temp_mask_ims = np.concatenate((glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.SEGMENTATION.fits'),glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.SEGMENTATION.fits')))
if do_dither != '1' and int(do_dither) <= 5:
	temp_back_ims = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.fits'),glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.fits')))
	temp_back_weights = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.WEIGHT.fits'),glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.WEIGHT.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.WEIGHT.fits')))
 	temp_mask_ims = np.concatenate((glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.MASK.fits'),glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.MASK.fits'),glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.MASK.fits')))
#	temp_mask_ims = np.concatenate((glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.SEGMENTATION.fits'),glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.SEGMENTATION.fits'),glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.SEGMENTATION.fits')))
if int(do_dither) > 5:
	temp_back_ims = glob.glob(back_dir+'survey_t*_*d*_*'+do_filt+'_*short.fits')
	temp_back_weights = glob.glob(back_dir+'survey_t*_*d*_*'+do_filt+'_*short.WEIGHT.fits')
 	temp_mask_ims = glob.glob(mask_dir+'survey_t*_*d*_*'+do_filt+'_*short.MASK.fits')
#	temp_mask_ims = glob.glob(mask_dir+'survey_t*_*d*_*'+do_filt+'_*short.SEGMENTATION.fits')

# print temp_back_ims
# print temp_back_weights
print temp_mask_ims

im_mjd = im_h['MJD-OBS']
print "\nImage observation date = %f" % im_mjd

#Delete tiles with bright, extended objects in them
#SCABS: tile1 (CenA), tile11 (wCen), tile24(wCen)
if do_dither == '1':
	gv_dithers = [do_dither,str(int(do_dither)+1)]
if do_dither != '1' and int(do_dither) <= 5:
	gv_dithers = [str(int(do_dither)-1),do_dither,str(int(do_dither)+1)]
if int(do_dither) > 5:
	gv_dithers = [str(ii+1) for ii in range(15)]
	
#print gv_dithers
bv_tiles = ['1','11','24']
for gv_d in gv_dithers:
	for bv in bv_tiles:
		if back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.fits' in temp_back_ims:
			temp_back_ims = np.delete(temp_back_ims, list(temp_back_ims).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.fits'),0)
			temp_back_weights = np.delete(temp_back_weights, list(temp_back_weights).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.WEIGHT.fits'), 0)
 			temp_mask_ims = np.delete(temp_mask_ims, list(temp_mask_ims).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.MASK.fits'), 0)
#			temp_mask_ims = np.delete(temp_mask_ims, list(temp_mask_ims).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.SEGMENTATION.fits'), 0)

for imfile in temp_back_ims:
 	if imfile.replace('short.fits','short.MASK.fits') not in temp_mask_ims:
#	if imfile.replace('short.fits','short.SEGMENTATION.fits') not in temp_mask_ims:
		temp_back_ims = np.delete(temp_back_ims, list(temp_back_ims).index(imfile), 0)		
for wfile in temp_back_weights:
 	if wfile.replace('short.WEIGHT.fits','short.MASK.fits') not in temp_mask_ims:
#	if wfile.replace('short.WEIGHT.fits','short.SEGMENTATION.fits') not in temp_mask_ims:
		temp_back_weights = np.delete(temp_back_weights, list(temp_back_weights).index(wfile), 0)		
		
# print len(temp_back_ims), len(temp_back_weights), len(temp_mask_ims)
# print temp_back_ims
# print temp_back_weights
# print temp_mask_ims

back_mjds = []
for ii in range(len(temp_back_ims)):
	back_mjds.append(float(subprocess.Popen("dfits "+temp_back_ims[ii]+" | fitsort MJD-OBS | awk 'NR>1 {print $2}'", shell=True, stdout=subprocess.PIPE).stdout.read().replace('\n','')))
back_mjds = np.array(back_mjds)

time_diffs = back_mjds - im_mjd
#print len(temp_back_ims), len(temp_back_weights), len(temp_mask_ims), len(time_diffs)
#print time_diffs

back_ims = []
back_weights = []
mask_ims = []

# print len(temp_back_ims), len(temp_back_weights), len(temp_mask_ims)

for idx in find_nearest(time_diffs,n_backs,0):
	back_ims.append(temp_back_ims[idx])
	back_weights.append(temp_back_weights[idx])
	mask_ims.append(temp_mask_ims[idx])

back_ims = np.array(back_ims)
back_weights = np.array(back_weights)
mask_ims = np.array(mask_ims)

# print back_ims
# print back_weights
# print mask_ims

print "\nUsing the following frames for the background subtraction: " 
print back_ims

back_mjds = []
for ii in range(n_backs):
	back_mjds.append(float(subprocess.Popen("dfits "+temp_back_ims[ii]+" | fitsort MJD-OBS | awk 'NR>1 {print $2}'", shell=True, stdout=subprocess.PIPE).stdout.read().replace('\n','')))

#Load background images and weight maps
back_im_cube = np.zeros( (len(back_ims), im_nchip, im_size[0], im_size[1]) )
for ii in range(len(back_ims)):
	print "\nLoading background image %i/%i into data cube..." % (ii+1,len(back_ims))
	back_im_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
	back_im_cube[ii] = load_im(back_ims[ii],back_im_temp)

#print back_im_cube
#print len(back_im_cube), len(back_im_cube[0]), len(back_im_cube[0][0]), len(back_im_cube[0][0][0])

#back_w_cube = np.zeros( (len(back_weights), im_nchip, im_size[0], im_size[1]) )

#Weight scale the backgroundi mages
# for ii in range(len(back_ims)):
# 	print "\nWeighting background image %i/%i..." % (ii+1,len(back_weights))
# 	back_weight_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
# 	back_w_data = load_im(back_weights[ii],back_weight_temp)
# 	back_im_cube[ii] = weight_im(back_w_data[0],back_im_cube[ii])
# 
# del back_w_data

#Load mask images and mask the background images
for ii in range(len(mask_ims)):
	print "\nMasking background image %i/%i..." % (ii+1,len(mask_ims))
	mask_im_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
	mask_im_data = load_im(mask_ims[ii],mask_im_temp)

	for jj in range(len(back_im_cube[0])):
 		bv_mask = (mask_im_data[0][jj] == 1.)
#		bv_mask = (mask_im_data[0][jj] != 0.)
		back_im_cube[ii][jj][bv_mask] = np.nan

for ii in range(len(back_im_cube)):
	print "Normalizing background image %i/%i..." % (ii+1,len(back_im_cube))
	for jj in range(im_nchip):
		back_im_cube[ii][jj] = stats.nanmedian(im_masked[0][jj]/back_im_cube[ii][jj])*back_im_cube[ii][jj]

# plt.figure()
# plt.imshow(back_im_cube[0][0])#,vmin=np.min(back_im_cube[0][0]),vmax=500)
# plt.colorbar()

# print "Saving background data..."
# hdulist_back_out = hdulist
# for ii in range(im_nchip):
# 	hdulist_back_out[ii+1].data = back_im_cube[0][ii]
# hdulist_back_out.writeto('/Users/matt/Desktop/deleteme_nofill.fits',clobber=True)
# # hdulist_back_out.writeto('%s' % sky_im_out,clobber=True)
# hdulist_back_out.close()

if sys.argv[1] == "dynamic tps":
	back_im_cube = infill(back_im_cube)


#print "Saving background data..."
#for jj in range(len(back_im_cube)):
#	hdulist_back_out = hdulist
#	for ii in range(im_nchip):
#		hdulist_back_out[ii+1].data = back_im_cube[jj][ii]
#	hdulist_back_out.writeto('/Users/matt/Desktop/back%i_filled.fits' % (jj+1),clobber=True)

# hdulist_back_out.writeto('%s' % sky_im_out,clobber=True)
#hdulist_back_out.close()

del mask_im_data

########################################################################
####Calculate the background level from all of the background images####
########################################################################

# for ii in range(10):
#  	print "Chip #%i before subtraction:" % (ii+1)
# # 	print back_im_cube[0][ii][0][0], back_im_cube[1][ii][0][0], back_im_cube[2][ii][0][0]
# # 	print np.median([back_im_cube[0][ii][0][0], back_im_cube[1][ii][0][0], back_im_cube[2][ii][0][0]])
# # 	print np.median([back_im_cube[kk][ii][:][:] for kk in range(len(back_im_cube))])
#  	print im_data[0][ii][0][0]

if sys.argv[1] == "median" or sys.argv[1] == "median global":# or sys.argv[1] == "dynamic tps":
	q = Queue()
	results = np.zeros( (im_nchip, 2, im_size[0], im_size[1]) )
	results = list(results)
	results.append(0)
	results = np.array(results)
	med_back_data = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )

	n_procs_max = 10
	chip_start = 0
	max_chips = 60
	procs = []

	while chip_start < max_chips:#im_nchips:
		print "Processing Chip #%i-%i..." % (chip_start+1, chip_start+n_procs_max)
		for chip_num in range(n_procs_max):
	#  	  print "Processing Chip #%i" % (chip_num+1+chip_start)
		  if sys.argv[1] == "median":
			  procs.append(Process(target=median_subtract, args=(back_im_cube,med_back_data[0][chip_num+chip_start],im_data[0][chip_num+chip_start],chip_num+chip_start,q)))
		  if sys.argv[1] == "median global":
			  procs.append(Process(target=median_global_subtract, args=(back_im_cube,med_back_data[0][chip_num+chip_start],im_data[0][chip_num+chip_start],chip_num+chip_start,q)))
		  #print "Proc %i start" % (chip_num+1+chip_start)
		  if sys.argv[1] == "bispline":
			  procs.append(Process(target=bspline_subtract, args=(back_im_cube,med_back_data[0][chip_num+chip_start],im_data[0][chip_num+chip_start],chip_num+chip_start,q)))
# 	 	  if sys.argv[1] == "dynamic tps":
# 	 	  	  procs.append(Process(target=tps_subtract, args=(back_im_cube,med_back_data[0][chip_num+chip_start],im_data[0][chip_num+chip_start],chip_num+chip_start,q)))

		  procs[chip_num+chip_start].start()
		  time.sleep(3.0)
		for chip_num in range(n_procs_max):
		  #print "Queue %i start" % (chip_num+1+chip_start)
		  results[chip_num+chip_start] = q.get()#result1 = q.get()
		  #im_data[0][chip_num+chip_start] = results[chip_num+chip_start][0] #result1[0]
		  #med_back_data[0][chip_num+chip_start] = results[chip_num+chip_start][1] #result1[1]
		  im_data[0][results[chip_num+chip_start][2]] = results[chip_num+chip_start][0] #result1[0]
		  med_back_data[0][results[chip_num+chip_start][2]] = results[chip_num+chip_start][1] #result1[1]
		chip_start+=n_procs_max

#Spline surface background subtraction
#Fit a spline to the surface of each chip, for each background image

if sys.argv[1] == "dynamic tps":
	med_back_data = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
	for chip_num in range(1):#range(im_nchip):
		print "Processing Chip #%i..." % (chip_num+1)
		im_data[0][chip_num], med_back_data[0][chip_num] = tps_subtract(back_im_cube,im_data[0][chip_num],chip_num)

# for ii in range(10):
#  	print "Chip #%i after subtraction:" % (ii+1)
# # 	print back_im_cube[0][ii][0][0], back_im_cube[1][ii][0][0], back_im_cube[2][ii][0][0]
# # 	print np.median([back_im_cube[0][ii][0][0], back_im_cube[1][ii][0][0], back_im_cube[2][ii][0][0]])
# # 	print np.median([back_im_cube[kk][ii][:][:] for kk in range(len(back_im_cube))])
#  	print im_data[0][ii][0][0], med_back_data[0][ii][0][0]

#exit()

#Save the background subtracted image and the background image
hdulist = pyfits.open(test_image)

print "Saving image data..."
hdulist_out = hdulist
for ii in range(im_nchip):
	hdulist_out[ii+1].data = im_data[0][ii]
hdulist_out.writeto('/Users/matt/Desktop/deleteme_image.fits',clobber=True)
# hdulist_out.writeto('%s' % ss_im_out,clobber=True)
hdulist_out.close()

print "Saving background data..."
hdulist_back_out = hdulist
for ii in range(im_nchip):
	hdulist_back_out[ii+1].data = med_back_data[0][ii]
hdulist_back_out.writeto('/Users/matt/Desktop/deleteme_back.fits',clobber=True)
# hdulist_back_out.writeto('%s' % sky_im_out,clobber=True)
hdulist_back_out.close()

print "Done in %5.2f seconds." % (time.time() - start_time)
