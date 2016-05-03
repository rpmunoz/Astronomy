import os
import sys
import shlex
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits

def create_sextractor_cat(im_file, wim_file='', zp=30., saturation=65536.)

	sex_conf_file='sex/default.sex'
	sex_param_file='sex/default.param'
	sex_conv_file='sex/default.conv'
	sex_detect_thresh='2.'
	sex_analysis_thresh='2.'
	sex_detect_minarea='3'
	sex_deblend_mincont='0.001'
	sex_zp='{:.2f}'.format(zp)
	sex_satur_level='{:.2f}'.format(saturation)
	sex_mem_pixstack='1000000'
	sex_cat_type='ASCII_HEAD'
	sex_weight_type='MAP_WEIGHT'
	sex_back_size='64'
	sex_back_filtersize='3'
	
	cmd='sex -c '+sex_conf_file+' -DETECT_THRESH '+sex_detect_thresh+' -DETECT_MINAREA '+sex_detect_minarea+' -ANALYSIS_THRESH '
	+sex_analysis_thresh+' -DEBLEND_MINCONT '+sex_deblend_mincont+' -BACK_SIZE '+sex_back_size+' -BACK_FILTERSIZE '
	+sex_back_filtersize+' -PARAMETERS_NAME '+sex_param_file+' -FILTER_NAME '+sex_conv_file+' -CATALOG_TYPE '
	+sex_cat_type+' -WEIGHT_TYPE '+sex_weight_type+' -MAG_ZEROPOINT '+sex_zp+' -SATUR_LEVEL '+sex_satur_level+' -MEMORY_PIXSTACK '
	sex_mem_pixstack+' -CHECKIMAGE_TYPE '+sex_check_type+' -CHECKIMAGE_NAME '+sex_check_name+' -CATALOG_NAME '
	+sex_cat_name+' -WEIGHT_IMAGE '+wim_file+' '+ im_file
cmd_args=shlex.split(cmd)
print cmd
proc = subprocess.Popen(cmd_args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = proc.communicate()
print stderr[0:200]


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