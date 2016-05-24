#! /usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import sys,os
import subprocess
import numpy as np
import random
import time
import cv2 as cv
import pyfits
from pyfits import getheader
import multiprocessing, Queue
import ctypes
#import agpy
#from agpy import gaussfit

class Worker(multiprocessing.Process):
 
	def __init__(self, work_queue, result_queue):
 
		# base class initialization
		multiprocessing.Process.__init__(self)
 
		# job management stuff
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.kill_received = False
 
	def run(self):
		while not self.kill_received:
 
			# get a task
			try:
				i_range, psf_file = self.work_queue.get_nowait()
			except Queue.Empty:
				break
 
			# the actual processing
			print "Adding artificial stars - index range=", i_range

			radius=16
			x_c,y_c=( (psf_size[1]-1)/2, (psf_size[2]-1)/2 )
			x,y=np.meshgrid(np.arange(psf_size[1])-x_c,np.arange(psf_size[2])-y_c)
			distance = np.sqrt(x**2 + y**2)

			for i in range(i_range[0],i_range[1]):
				psf_xy=np.zeros(psf_size[1:3], dtype=float)
				j=0
				for i_order in range(psf_order+1):
					j_order=0
					while (i_order+j_order < psf_order+1):
						psf_xy += psf_data[j,:,:] * ((mock_y[i]-psf_offset[1])/psf_scale[1])**i_order * ((mock_x[i]-psf_offset[0])/psf_scale[0])**j_order
						j_order+=1
						j+=1

#				psf_xy = np.ma.masked_where(distance > radius, psf_xy)
#				psf_xy = np.ma.masked_less(psf_xy, 0.)
#				psf_xy = np.ma.filled(psf_xy, 0.)

				psf_factor=10.**( (30.-mock_mag[i])/2.5)/np.sum(psf_xy)
				psf_xy *= psf_factor

				npsf_xy=cv.resize(psf_xy,(npsf_size[0],npsf_size[1]),interpolation=cv.INTER_LANCZOS4)
				npsf_factor=10.**( (30.-mock_mag[i])/2.5)/np.sum(npsf_xy)
				npsf_xy *= npsf_factor

				im_rangex=[max(mock_x[i]-npsf_size[1]/2,0), min(mock_x[i]-npsf_size[1]/2+npsf_size[1], im_size[1])]
				im_rangey=[max(mock_y[i]-npsf_size[0]/2,0), min(mock_y[i]-npsf_size[0]/2+npsf_size[0], im_size[0])]
				npsf_rangex=[max(-1*(mock_x[i]-npsf_size[1]/2),0), min(-1*(mock_x[i]-npsf_size[1]/2-im_size[1]),npsf_size[1])]
				npsf_rangey=[max(-1*(mock_y[i]-npsf_size[0]/2),0), min(-1*(mock_y[i]-npsf_size[0]/2-im_size[0]),npsf_size[0])]

				im_data[im_rangey[0]:im_rangey[1], im_rangex[0]:im_rangex[1]] += npsf_xy[npsf_rangey[0]:npsf_rangey[1], npsf_rangex[0]:npsf_rangex[1]]
#				im_data[im_rangey[0]:im_rangey[1], im_rangex[0]:im_rangex[1]] = 10.

				#print "Before i=", i, " sum= ", np.sum(im_data[im_rangey[0]:im_rangey[1], im_rangex[0]:im_rangex[1]])
				#print "After  i=", i, " sum= ", np.sum(im_data[im_rangey[0]:im_rangey[1], im_rangex[0]:im_rangex[1]])

			print 'Done'
			self.result_queue.put(id)
			# store the result


class Worker_sex(multiprocessing.Process):
 
	def __init__(self, work_queue, result_queue):
 
		# base class initialization
		multiprocessing.Process.__init__(self)
 
		# job management stuff
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.kill_received = False

	def run(self):
		while not self.kill_received:
 
			# get a task
			try:
				i_thread, i_range = self.work_queue.get_nowait()
			except Queue.Empty:
				break

			fwhm=1.0
			pixel_scale=0.263
			weight_type='MAP_WEIGHT'
			checkimage_type='NONE'
			checkimage_file='NONE'
			satur_level=4.3e5
			analysis_thresh=2.0
			detect_minarea=3
			detect_thresh=1.4
			phot_apertures=",".join(["%.2f" % x for x in 2*np.array((0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.))*fwhm/pixel_scale])
			filter_name='sex_config/gauss_5.0_9x9_t1_z.conv'
			xml_name='scabs_tile2_z_sex.xml'
		
			# the actual processing
			log_file="decam_completeness_sex_thread%d.log" % i_thread
			for i in range(i_range[0],i_range[1]):
				mock_mag_cat_file=mock_im_file[i].replace('mock_im','sex_cat')
				mock_mag_cat_file = mock_mag_cat_file.replace('fits','ldac')
# 				print mock_mag_cat_file
#/Volumes/Q6/matt/stacks/completeness/mock/sex_cat_t1_i_1_0001.ldac	
# 				command = "sex %s -c sex_config/cfht_wircam.sex -PARAMETERS_NAME sex_config/cfht_wircam_psfmodel.param -CATALOG_TYPE FITS_LDAC -CATALOG_NAME %s -SEEING_FWHM %.2f -WEIGHT_TYPE %s -WEIGHT_THRESH 0. -WEIGHT_IMAGE %s -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s -SATUR_LEVEL %d -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 30 -BACK_SIZE 250 -BACK_FILTERSIZE 3 -MASK_TYPE CORRECT -ANALYSIS_THRESH %.2f -DETECT_MINAREA %d -DETECT_THRESH %.2f -DEBLEND_MINCONT 0.0000001 -INTERP_TYPE ALL -INTERP_MAXXLAG 1 -INTERP_MAXYLAG 1 -FLAG_TYPE OR -FLAG_IMAGE %s -PHOT_AUTOPARAMS 2.3,4.0 -PHOT_FLUXFRAC 0.5 -PHOT_APERTURES %s -PIXEL_SCALE %.4f -FILTER Y -FILTER_NAME %s -WRITE_XML Y -XML_NAME %s -PSF_NAME %s -PSF_NMAX 1" % (mock_im_file[i], sex_cat_file[i], fwhm, weight_type, weight_file[i], checkimage_type, checkimage_file, satur_level, analysis_thresh, detect_minarea, detect_thresh, flag_file[i], phot_apertures, pixel_scale, filter_name, xml_name, psf_file[i] )
#				command = "sex %s -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_TYPE FITS_LDAC -CATALOG_NAME %s -SEEING_FWHM %.2f -WEIGHT_TYPE %s -WEIGHT_THRESH 0. -WEIGHT_IMAGE %s -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s -SATUR_LEVEL %d -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 30 -BACK_SIZE 250 -BACK_FILTERSIZE 3 -MASK_TYPE CORRECT -ANALYSIS_THRESH %.2f -DETECT_MINAREA %d -DETECT_THRESH %.2f -DEBLEND_MINCONT 0.0000001 -INTERP_TYPE ALL -INTERP_MAXXLAG 1 -INTERP_MAXYLAG 1 -PHOT_AUTOPARAMS 2.3,4.0 -PHOT_FLUXFRAC 0.5 -PHOT_APERTURES %s -PIXEL_SCALE %.4f -FILTER Y -FILTER_NAME %s -WRITE_XML Y -XML_NAME %s -PSF_NAME %s -PSF_NMAX 1" % (mock_im_file[i], sex_cat_file[i], fwhm, weight_type, weight_file[i], checkimage_type, checkimage_file, satur_level, analysis_thresh, detect_minarea, detect_thresh, phot_apertures, pixel_scale, filter_name, xml_name, psf_file[i] )
				command = "sex %s -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_TYPE FITS_LDAC -CATALOG_NAME %s -SEEING_FWHM %.2f -WEIGHT_TYPE %s -WEIGHT_THRESH 0. -WEIGHT_IMAGE %s -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s -SATUR_LEVEL %d -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 30 -BACK_SIZE 250 -BACK_FILTERSIZE 3 -MASK_TYPE CORRECT -ANALYSIS_THRESH %.2f -DETECT_MINAREA %d -DETECT_THRESH %.2f -DEBLEND_MINCONT 0.0000001 -INTERP_TYPE ALL -INTERP_MAXXLAG 1 -INTERP_MAXYLAG 1 -PHOT_AUTOPARAMS 2.3,4.0 -PHOT_FLUXFRAC 0.5 -PHOT_APERTURES %s -PIXEL_SCALE %.4f -FILTER Y -FILTER_NAME %s -WRITE_XML Y -XML_NAME %s -PSF_NAME %s -PSF_NMAX 1" % (mock_im_file[i], mock_mag_cat_file, fwhm, weight_type, weight_file[i], checkimage_type, checkimage_file, satur_level, analysis_thresh, detect_minarea, detect_thresh, phot_apertures, pixel_scale, filter_name, xml_name, psf_file[i] )
				print command
				with open(log_file, "a") as log:
					result=subprocess.call(command, stderr=log, stdout=log, shell=True)
#				log=subprocess.Popen(command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

				log.close()
				print "SExtractor thread: %d - iteration: %d is done!" % (i_thread, i)

			self.result_queue.put(id)


#		fig = plt.figure()
#		ax = fig.add_subplot(111, projection='3d')
#		x,y=np.meshgrid(np.arange(npsf_size[0])-12,np.arange(npsf_size[1])-12)
#		surf = ax.plot_surface(x, y, npsf_xy, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=True)
#		plt.show()
#		if i >= 0: break

if __name__ == "__main__":

#		global im_data, im_size, mock_x, mock_y, mock_mag
#		global psf_data, psf_size, npsf_size, psf_order, psf_offset, psf_scale, psf_pixstep

		n_cpu=2
		n_core=6
		n_processes=n_cpu*n_core*1

# 		input_mock_file=sys.argv[1]

#		if os.path.exists(output_im_file):
#			os.remove(output_im_file)	
# 		input_data_dtype=np.dtype({'names':['im_file','weight_file','psf_file','mock_mag_file','mock_im_file','sex_cat_file'],'formats':['S200','S200','S200','S200','S200','S200']})
# 		input_data=np.loadtxt(input_mock_file, skiprows=1, dtype=input_data_dtype)
# 		im_file=input_data['im_file']
# 		weight_file=input_data['weight_file']
# #		flag_file=input_data['flag_file']
# 		psf_file=input_data['psf_file']
# 		mock_mag_file=input_data['mock_mag_file']
# 		mock_im_file=input_data['mock_im_file']
# 		sex_cat_file=input_data['sex_cat_file']
# 		input_n=im_file.size

# 		input_data_dtype=np.dtype({'names':['im_file','weight_file','psf_file','mock_mag_file','mock_im_file','sex_cat_file'],'formats':['S200','S200','S200','S200','S200','S200']})
# 		input_data=np.loadtxt(input_mock_file, skiprows=1, dtype=input_data_dtype)
		im_file=[sys.argv[1]]
		weight_file=[sys.argv[2]]
#		flag_file=input_data['flag_file']
		psf_file=[sys.argv[3]]
		mock_mag_file=[sys.argv[4]]
		mock_im_file=[sys.argv[5]]
		sex_cat_file=[sys.argv[6]]
		input_n=1#im_file.size

		#n_processes=np.minimum(n_processes,input_n)
		print n_processes
		
		for i in range(input_n):
			print mock_im_file[i]
			if os.path.exists(mock_im_file[i]):
				print "Removing file ", mock_im_file[i]
				os.remove(mock_im_file[i])	
			if os.path.exists(sex_cat_file[i]):
				print "Removing file ", sex_cat_file[i]
				os.remove(sex_cat_file[i])	

# First, add artificial stars
		for i in range(0,input_n):
			hdulist = pyfits.open(psf_file[i])
			psf_h = hdulist[1].header
			psf_data = (hdulist[1].data)[0][0]
			hdulist.close()

			psf_order=psf_h['POLDEG1']
			psf_offset=[psf_h['POLZERO1'],psf_h['POLZERO2']]
			psf_scale=[psf_h['POLSCAL1'],psf_h['POLSCAL2']]
			psf_pixstep=psf_h['PSF_SAMP']
			psf_size=psf_data.shape
			npsf_size=(np.array(psf_size[1:3])*psf_pixstep).astype(int)

			mock_data=np.loadtxt(mock_mag_file[i], skiprows=1)
			mock_n=mock_data[:,0].size
#			mock_random=random.sample(xrange(mock_n),mock_n)
			mock_sort=np.argsort(mock_data[:,1])

			mock_x=mock_data[mock_sort,0]
			mock_y=mock_data[mock_sort,1]
			mock_mag=mock_data[mock_sort,2]

			print "Reading file ", im_file[i]
			hdu=pyfits.open(im_file[i])
			data=hdu[0].data
			im_size=data.shape

			im_data_base = multiprocessing.Array(ctypes.c_float, im_size[0]*im_size[1])
			im_data = np.ctypeslib.as_array(im_data_base.get_obj())
			im_data = im_data.reshape(im_size[0], im_size[1])
			im_data[:] = data
#			im_data[:] = 0.
			data=0
# 			print im_data.base.base, im_data_base.get_obj()
# 			assert im_data.base.base is im_data_base.get_obj()

		# run
		# load up work queue
			tic=time.time()
			j_step=np.int(np.ceil( mock_n*1./n_processes ))
			j_range=range(0,mock_n,j_step)
			j_range.append(mock_n)

			work_queue = multiprocessing.Queue()
			for j in range(np.size(j_range)-1):
				if work_queue.full():
					print "Oh no! Queue is full after only %d iterations" % j
				work_queue.put( (j_range[j:j+2], psf_file[i]) )
 
			# create a queue to pass to workers to store the results
			result_queue = multiprocessing.Queue()
			procs=[]
 
			# spawn workers
			for j in range(n_processes):
				worker = Worker(work_queue, result_queue)
				procs.append(worker)
				worker.start()
 
			# collect the results off the queue
#			while not work_queue.empty():
#				result_queue.get()
			for j in range(n_processes):
				result_queue.get()

			for p in procs:
				p.join()

			print 'Final Done'
			print "Writing file ", mock_im_file[i]
			hdu[0].data=im_data
			hdu.writeto(mock_im_file[i])
			print "%f s for parallel computation." % (time.time() - tic)


# Second, run Sextractor
		n_processes=n_cpu*n_core
		n_processes=np.minimum(n_processes,input_n)

		tic=time.time()
		j_step=np.int(np.ceil( input_n*1./n_processes ))
		j_range=range(0,input_n,j_step)
		j_range.append(input_n)

		work_queue = multiprocessing.Queue()
		for j in range(np.size(j_range)-1):
			if work_queue.full():
				print "Oh no! Queue is full after only %d iterations" % j
			work_queue.put( (j+1, j_range[j:j+2]) )

		# create a queue to pass to workers to store the results
		result_queue = multiprocessing.Queue()
		procs=[]

		# spawn workers
		for j in range(n_processes):
			worker = Worker_sex(work_queue, result_queue)
			procs.append(worker)
			worker.start()
			time.sleep(30)
 
		# collect the results off the queue
#		while not work_queue.empty():
#			result_queue.get()

		for j in range(n_processes):
			result_queue.get()

		for p in procs:
			p.join()

