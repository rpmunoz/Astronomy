# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Astromatic tools are being developed by Roberto Pablo Munoz, PhD

import os
import sys
import shlex
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits

def create_sex_cat(im_file, weight_file='', cat_file='', zp=30., saturation=65536., back_size=64, seeing=0.6, pixel_scale=0.186, phot_radius=1.8):

	# seeing: FWHM of the stellar sources in the image in units of arcsec
	# pixel_scale: Pixel scale of the image in units of arcsec/pixel
	# phot_radius: Radius of the circular aperture photometry in units of arcsec

	sex_im_file=im_file
	sex_weight_file=('NONE' if weight_file=='' else weight_file)
	sex_cat_file=(im_file.replace('.fits','_cat.ldac') if cat_file=='' else cat_file)
	
	sex_weight_type=('NONE' if weight_file=='' else 'MAP_WEIGHT')
	
	sex_conf_file='astromatic_config/sex_default.sex'
	sex_param_file='astromatic_config/sex_default.param'
	sex_conv_file='astromatic_config/sex_default.conv'
	sex_detect_thresh='2.'
	sex_analysis_thresh='2.'
	sex_detect_minarea='3'
	sex_deblend_mincont='0.001'
	sex_zp='{:.2f}'.format(zp)
	sex_seeing='{:.2f}'.format(seeing)
	sex_pixel_scale='{:.3f}'.format(pixel_scale)
	sex_phot_apertures='{:.3f}'.format(2.*phot_radius/pixel_scale)
	sex_satur_level='{:.2f}'.format(saturation)
	sex_mem_pixstack='1000000'
	sex_cat_type='FITS_LDAC' #'ASCII_HEAD'
	sex_back_size='{:d}'.format(back_size)
	sex_back_filtersize='3'
	sex_check_type='NONE'
	sex_check_file='NONE'
	
	cmd='sex -c '+sex_conf_file+' -DETECT_THRESH '+sex_detect_thresh+' -DETECT_MINAREA '+sex_detect_minarea+' -ANALYSIS_THRESH '+sex_analysis_thresh+\
		' -DEBLEND_MINCONT '+sex_deblend_mincont+' -BACK_SIZE '+sex_back_size+' -BACK_FILTERSIZE '+sex_back_filtersize+\
		' -PARAMETERS_NAME '+sex_param_file+' -FILTER_NAME '+sex_conv_file+' -CATALOG_TYPE '+sex_cat_type+' -WEIGHT_TYPE '+sex_weight_type+\
		' -MAG_ZEROPOINT '+sex_zp+' -SATUR_LEVEL '+sex_satur_level+' -MEMORY_PIXSTACK '+sex_mem_pixstack+' -CHECKIMAGE_TYPE '+sex_check_type+\
		' -CHECKIMAGE_NAME '+sex_check_file+' -SEEING_FWHM '+sex_seeing+' -PIXEL_SCALE '+sex_pixel_scale+' -PHOT_APERTURES '+sex_phot_apertures+\
		' -CATALOG_NAME '+sex_cat_file+' -WEIGHT_IMAGE '+sex_weight_file+' '+ sex_im_file

	cmd_args=shlex.split(cmd)
	print cmd
	proc = subprocess.Popen(cmd_args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	print stderr[0:300]


#cat_ref_file='data/2MASS_0338-3531_r70.cat'
#cat_file='data/XCS0083000101_r_d1.ldac'
#def convert_sextractor_to_scamp(cat_file)

#hdulist=fits.open(cat_ref_file)
#hdulist.info()
#cat_h0, cat_data0 = (hdulist[0].header, hdulist[0].data)
#cat_h1, cat_data1 = (hdulist[1].header, hdulist[1].data)
#cat_h2, cat_data2 = (hdulist[2].header, hdulist[2].data)
#hdulist.close()

#print cat_data.columns
#cat_data['NUMBER']

#fits_read, cat_fcb, cat_data0, cat_h0, exten=0
#fits_read, cat_fcb, cat_data1, cat_h1, exten=1
#fits_read, cat_fcb, cat_data2, cat_h2, exten=2
#fits_close, cat_fcb