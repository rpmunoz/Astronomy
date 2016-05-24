def draw_rect(xvals, yvals, col, alph):
	xy = zip(xvals,yvals)
	poly = Polygon(xy, facecolor=col, alpha=alph)
	return poly

def draw_arrow(x,y,dx,dy):
	arrow = Arrow(x,y,dx,dy,width=500.0,color='k')
	return arrow

def comp_mag_calc(pct,nx,ny,dm):
	nbins_x = nx
	nbins_y = ny
	dmag = dm
	xmax = np.max([np.max(mock_in_data['X']),np.max(mock_out_data['X'])])
	ymax = np.max([np.max(mock_in_data['Y']),np.max(mock_out_data['Y'])])
	bs_x = int(xmax/nbins_x)
	bs_y = int(ymax/nbins_y)
	xchk = 0
	ychk = 0

	mag_out = []
	for ii in range(nbins_x):
		temp_mags_out = []
		if xchk+bs_x < xmax:
			xlims = [xchk,xchk+bs_x]
		else:
			xlims = [xchk,xmax]

		ychk = 0

		for jj in range(nbins_y):
			if ychk+bs_y < ymax:
				ylims = [ychk,ychk+bs_y]
			else:
				ylims = [ychk,ymax]

			mask = np.logical_and(mock_in_data['X'] > xlims[0], mock_in_data['X'] <= xlims[1])
			temp_mock_in_data = mock_in_data[mask]
			mask = np.logical_and(temp_mock_in_data['Y'] > ylims[0], temp_mock_in_data['Y'] <= ylims[1])
			temp_mock_in_data = temp_mock_in_data[mask]

			mask = np.logical_and(mock_out_data['X'] > xlims[0], mock_out_data['X'] <= xlims[1])
			temp_mock_out_data = mock_out_data[mask]
			mask = np.logical_and(temp_mock_out_data['Y'] > ylims[0], temp_mock_out_data['Y'] <= ylims[1])
			temp_mock_out_data = temp_mock_out_data[mask]
	
			if len(temp_mock_in_data) > 0 and len(temp_mock_out_data) > 0:
				mag_chk = np.max([np.min(temp_mock_in_data['MAG']),np.min(temp_mock_out_data['MAG_AUTO'])])
				mag_max = np.min([np.max(temp_mock_in_data['MAG']),np.max(temp_mock_out_data['MAG_AUTO'])])
				mag_outs = []
				comps = []
				while mag_chk < mag_max:
					mask = np.logical_and(temp_mock_in_data['MAG'] >= mag_chk, temp_mock_in_data['MAG'] < mag_chk+dmag)
					temp_mag_in = temp_mock_in_data['MAG'][mask]
					mask = np.logical_and(temp_mock_out_data['MAG'] >= mag_chk, temp_mock_out_data['MAG'] < mag_chk+dmag)
					temp_mag_out = temp_mock_out_data[mask]

					mag_outs.append(mag_chk+dmag/2.)
					if len(temp_mag_in) > 0:
						comps.append(float(len(temp_mag_out['MAG'])/float(len(temp_mag_in))))
					else:
						comps.append(0.) 
					if mag_chk + dmag < mag_max:
						mag_chk += dmag
					else:
						mag_chk += mag_max-mag_chk

				mag_outs = np.array(mag_outs)
				comps = np.array(comps)

				spl = UnivariateSpline(mag_outs,comps,k=3,s=0)
				xx = np.linspace(np.min(mag_outs)-0.1,np.max(mag_outs)+0.1,1e4)

				comp_mags = spl(xx)
				comp = comp_mags[(np.abs(comp_mags-pct).argmin())]
				mag = xx[(np.abs(comp_mags-pct).argmin())]
	
				temp_mags_out.append(mag)

			else:
				temp_mags_out.append(0.)

			ychk += bs_y

		mag_out.append(temp_mags_out)

		xchk += bs_x

	return np.array(mag_out)

def plot_comp_mags(mags_out,pct,nx,ny,mean,std):
	xx = np.linspace(0,xmax,nx)
	yy = np.linspace(0,ymax,ny)		

	mags_out = np.array(mags_out)
	mags_out[(mags_out == 0.0)] = np.nan

	xgrid, ygrid = np.meshgrid(xx, yy)

# 	gs = gridspec.GridSpec(2,1, height_ratios=[0.05,1], width_ratios=[1])
# 	gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)

	fig = plt.figure(figsize=(12,12))
	ax = fig.add_subplot(111)
# 	ax = fig.add_subplot(gs[1,0])
# 	im = plt.imshow(mags_out.T, interpolation='lanczos', origin='lower',
# 					cmap=cm, vmin=xlim2[0], vmax=xlim2[1],
# 					extent=(np.max(ygrid), np.min(ygrid), np.max(xgrid),np.min(xgrid)))
	im = plt.imshow(mags_out.T, interpolation='lanczos', origin='lower',
					cmap=cm, vmin=xlim2[0], vmax=xlim2[1],
					extent=(np.max(ygrid), np.min(ygrid), np.max(xgrid),np.min(xgrid)))
	# levels = [15.,17.,19.,21.,23.]
	# CS = plt.contour(mag50_out.T, levels,
	#                  origin='lower',
	#                  linewidths=2,
	#                  extent=(np.max(ygrid), np.min(ygrid), np.max(xgrid),np.min(xgrid)))
	# plt.clabel(CS, levels[1::2],  # label every second level
	#            inline=1,
	#            fmt='%1.1f',
	#            fontsize=14)
		   
	ns = draw_arrow(25000,30000,0.,-2000)
	ew = draw_arrow(25000,30000,2000,0.)

	ax.add_patch(ns)
	ax.add_patch(ew)

	plt.setp(ax.get_xticklabels(), fontsize=14)
	ax.xaxis.set_tick_params(length=8,width=2)
	ax.yaxis.set_tick_params(length=8,width=2)
	plt.setp(ax.get_yticklabels(), fontsize=14)

	plt.annotate(r"$\mu_{%s',%s\%%}=%.2f$" % (do_filter,pct,mean),xy=(27600,1250),color='k',fontsize=20)
	plt.annotate(r"$\sigma_{%s',%s\%%}=%.2f$" % (do_filter,pct,std),xy=(27600,2250),color='k',fontsize=20)

	plt.annotate(r"N",xy=(24900,27900),color='k',fontsize=16)
	plt.annotate(r"E",xy=(27600,29800),color='k',fontsize=16)
	plt.xlabel(r'Y Pixel Value',fontsize=20)
	plt.ylabel(r'X Pixel Value',fontsize=20)
# 
# 	if do_tile == '1':
# 		plt.title(r"Tile %s" % do_tile, fontsize=18)
# 	if do_tile != '1':
# 		plt.title(r"Outer Ring", fontsize=18)

# 	cbax = fig.add_subplot(gs[0,0])
# 	cb = Colorbar(ax=cbax, mappable=im, orientation='horizontal', ticklocation='top')
# 	cb.set_label(r"%s%% $%s'$-band Completeness Magnitude" % (pct,do_filter),fontsize=18)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("top", size="5%", pad=0.05)

	CB = plt.colorbar(im,cax=cax,orientation='horizontal',ticklocation='top')#orientation='horizontal',ticklocation='top')
	CB.set_label(r"%s%% $%s'$-band Completeness Magnitude" % (pct,do_filter),fontsize=18)


	plt.savefig("t%s_%s_mag%s_compmap.pdf" % (do_tile,do_filter,pct),bbox_inches='tight')

def fit(xs,ys):
	est = smooth.NonParamRegression(xs, ys, method=npr_methods.LocalPolynomialKernel(q=2))
	est.fit()
	return est

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
from scipy.spatial.distance import cdist
from astropy.coordinates import SkyCoord
from matplotlib.ticker import NullFormatter
from astropy.stats import sigma_clip
from astropy.table import Table, Column
from scipy.interpolate import UnivariateSpline, interp2d, Rbf
from matplotlib.patches import Polygon, Arrow# mpatches, Polygon
import sys
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
import pyqt_fit.bootstrap as bs
import numpy.ma as ma
import time
from astroML.resample import bootstrap

# a = [[0.,1.],[2.,3.]]
# mask = (np.array(a) > 0.)
# print mask
# exit()

prefix = '/Users/matt/Projects/CNTAC2014A/CN2014A-92/completeness/'

do_tile = sys.argv[1]
do_filter = sys.argv[2]

if do_filter == 'i':
	cat_file = '%ssex/survey_tile%s_%s_psf.ldac' % (prefix,do_tile,do_filter)
else:
	cat_file = '%ssex/survey_tile%s_%s_psf_ALIGNi.ldac' % (prefix,do_tile,do_filter)
	
mock_in_files = glob.glob('%smock/mock_mag_t%s_%s_*.dat' % (prefix,do_tile,do_filter))
mock_out_files = glob.glob('%smock/matched_t%s_%s*.ldac' % (prefix,do_tile,do_filter))

print mock_in_files
print mock_out_files
print cat_file

x_coords_in = []
y_coords_in = []
x_coords_out = []
y_coords_out = []
mags_in_tot = []
mags_in_out = []
mags_out = []

for ii in range(len(mock_out_files)):
   print "Processing file: %s" % mock_out_files[ii]
   temp_cat_hdu = fits.open(cat_file)
   temp_mock_out_hdu = fits.open(mock_out_files[ii])

   temp_cat_data = temp_cat_hdu[2].data
   temp_mock_in_data = ascii.read(mock_in_files[ii])
   temp_mock_out_data = temp_mock_out_hdu[1].data

   mask = (temp_cat_data['MAG_AUTO'] < 99.0)
   temp_cat_data = temp_cat_data[mask]
   # mask = (mock_in_data['MAG'][:,0] < 99.0)
   # mock_in_data = mock_in_data[mask]
   mask = (temp_mock_out_data['MAG_AUTO'] < 99.0)
   temp_mock_out_data = temp_mock_out_data[mask]

   x_coords_in = np.concatenate((np.array(x_coords_in),np.array(temp_mock_in_data['X'])))
   y_coords_in = np.concatenate((np.array(y_coords_in),np.array(temp_mock_in_data['Y'])))
   x_coords_out = np.concatenate((np.array(x_coords_out),np.array(temp_mock_out_data['X_IMAGE'])))
   y_coords_out = np.concatenate((np.array(y_coords_out),np.array(temp_mock_out_data['Y_IMAGE'])))
   mags_in_tot = np.concatenate((np.array(mags_in_tot),np.array(temp_mock_in_data['MAG'])))
   mags_in_out = np.concatenate((np.array(mags_in_out),np.array(temp_mock_out_data['MAG'])))
   mags_out = np.concatenate((np.array(mags_out),np.array(temp_mock_out_data['MAG_AUTO'])))

mock_in_data = {'X': x_coords_in, 'Y': y_coords_in, 'MAG': mags_in_tot}  
mock_in_data = Table(mock_in_data)
mock_out_data = {'X': x_coords_out, 'Y': y_coords_out, 'MAG': mags_in_out, 'MAG_AUTO': mags_out}
mock_out_data = Table(mock_out_data)

if do_tile == '1':
  if do_filter == 'u':
	  xlim1 = [19.,26.]
  if do_filter == 'g':
	  xlim1 = [16.,24.5]
  if do_filter == 'r':
	  xlim1 = [16.,24.]
  if do_filter == 'i':
	  xlim1 = [15.,25.]
  if do_filter == 'z':
	  xlim1 = [15.,24.]
if do_tile == '2':
  if do_filter == 'u':
	  xlim1 = [19.,25.5]
  if do_filter == 'g':
	  xlim1 = [16.,24.5]
  if do_filter == 'r':
	  xlim1 = [16.,24.]
  if do_filter == 'i':
	  xlim1 = [16.,23.5]
  if do_filter == 'z':
	  xlim1 = [16.5,23.]

#Non-parametric regression with bootstrap
grid = np.r_[xlim1[0]:xlim1[1]:512j]
xs = mock_out_data['MAG']
ys = mock_out_data['MAG_AUTO']-mock_out_data['MAG']


# plt.figure()
# plt.plot(xs, ys, 'b,')

ys = sigma_clip(ys,sigma=10.0,iters=None)
mask = (ma.getmaskarray(ys) == False)	
xs = xs[mask]
ys = ys[mask]

# plt.plot(xs,ys,'r,')
# plt.ylim(-1.0,1.0)
# plt.show()
# exit()

# lin_fit = smooth.NonParamRegression(xs, ys, method=npr_methods.LocalPolynomialKernel(q=1))
# cub_fit = smooth.NonParamRegression(xs, ys, method=npr_methods.LocalPolynomialKernel(q=3))
quad_fit = fit(xs,ys)#smooth.NonParamRegression(xs, ys, method=npr_methods.LocalPolynomialKernel(q=2))
# plt.figure()
# plt.plot(grid, quad_fit(grid), 'r', lw=2)
# plt.show()
# exit()
# print quad_fit(25.5)
# exit()
# quad_fit.fit()
# cub_fit.fit()
# lin_fit.fit()

# print "Bootstrapping..."
# time = time.time()
# boots = bs.bootstrap(fit, xs, ys, eval_points = grid, CI = (95,99))
# print "Bootstrapping took %.2f seconds." % (time.time()-time)

# plt.figure(figsize=(14,6))
# plt.plot(xs, ys, ',', label='Data')
# plt.plot(grid, quad_fit(grid), 'k', label='quadratic', linewidth=2)
# plt.plot(grid, boots.CIs[0][0,0], 'g--', label='95% CI', linewidth=2)
# plt.plot(grid, boots.CIs[0][0,1], 'g--', linewidth=2)
# plt.plot(grid, boots.CIs[1][0,0], 'r--', label='95% CI', linewidth=2)
# plt.plot(grid, boots.CIs[1][0,1], 'r--', linewidth=2)
# 
# # plt.fill_between(grid, result.CIs[0][0,0], result.CIs[0][0,1], color='g', alpha=0.25)
# plt.xlim(19,26)
# plt.ylim(-0.5,0.5)
# plt.show()
# exit()

# blah_mag = mock_in_data['MAG'][(mock_in_data['MAG'] > 27.)]
# blah_auto = mock_out_data['MAG_AUTO'][(mock_out_data['MAG_AUTO'] > 27.)]
# 
# print blah_mag
# print blah_auto
# exit()


print len(mock_in_data['X']), len(mock_out_data['X'])

dmag = 0.2
mag_chk = np.max([np.min(mock_in_data['MAG']),np.min(mock_out_data['MAG_AUTO'])])
#If tile = 2, and filter = i, set mag_max = 25.5
#If tile = 1, and filter = z, set mag_max = 26.5
#If tile = 2, and filter = z, set mag_max = 25.5
mag_max = 25.5#np.min([np.max(mock_in_data['MAG']),np.max(mock_out_data['MAG_AUTO'])])
mag_outs = []
comps = []
diffs = []
diff_errs = []
boots_out = []
error_fileout = open('tile%s_filter%s_errors.txt' % (do_tile,do_filter),'w')
print >> error_fileout, "mag1	mag2	1sig_boot"
while mag_chk < mag_max:
	print mag_chk, mag_max
	mask = np.logical_and(mock_in_data['MAG'] >= mag_chk, mock_in_data['MAG'] < mag_chk+dmag)
	temp_mag_in = mock_in_data['MAG'][mask]
	mask = np.logical_and(mock_out_data['MAG'] >= mag_chk, mock_out_data['MAG'] < mag_chk+dmag)
	temp_mag_out = mock_out_data[mask]

# 	plt.figure()
# 	plt.plot(temp_mag_out['MAG'],temp_mag_out['MAG_AUTO']-temp_mag_out['MAG'],',')

 	temp_diffs = sigma_clip(temp_mag_out['MAG_AUTO']-temp_mag_out['MAG'],sigma=2.3,iters=None)
	mask = (ma.getmaskarray(temp_diffs) == False)	
# 	temp_mag_out = temp_mag_out[mask]
	temp_diffs = temp_diffs[mask]

# 	plt.plot(temp_mag_out['MAG'],temp_diffs,'r,')

	print len(temp_diffs)
	if len(temp_diffs) > 1:
		temp_boots = bootstrap(temp_diffs,10000,np.std)
	else:
		temp_boots = 0.0
	print >> error_fileout, "%.2f	%.2f	%.2e" %	(mag_chk, mag_chk+dmag, temp_boots)
	boots_out.append(temp_boots)
# 	temp_x = mag_chk+dmag/2.
# 	plt.plot([mag_chk,mag_chk+dmag],[quad_fit(mag_chk),quad_fit(mag_chk+dmag)],'g',lw=2)
# 	plt.plot([temp_x,temp_x],[quad_fit(temp_x)-3.*boots,quad_fit(temp_x)+3.*boots],'k',lw=2)
# 	plt.plot([temp_x,temp_x],[quad_fit(temp_x)-boots,quad_fit(temp_x)+boots],'g',lw=2)
# 	print 3.*boots
# 
# 	plt.show()
# 	exit()

	mag_outs.append(mag_chk+dmag/2.)
	comps.append(float(len(temp_mag_out['MAG'])/float(len(temp_mag_in))))
	diffs.append(quad_fit(mag_chk+dmag/2.))

# 	if len(temp_diffs) > 2:
# 		diffs.append(quad_fit(mag_chk+dmag/2.))
# 	else:
# 		diffs.append(0.)
# 	diffs.append(np.mean(sigma_clip(temp_mag_out['MAG_AUTO']-temp_mag_out['MAG'],sigma=2.3,iters=None)))
# 	diff_errs.append(np.std(sigma_clip(temp_mag_out['MAG_AUTO']-temp_mag_out['MAG'],sigma=2.3,iters=None)))
 
	if mag_chk + dmag < mag_max:
		mag_chk += dmag
	else:
		mag_chk += mag_max-mag_chk

error_fileout.close()
mag_outs = np.array(mag_outs)
comps = np.array(comps)
boots_out = np.array(boots_out)
diffs = np.array(diffs)
# diff_errs = np.array(diff_errs)

nullfmt = NullFormatter() 

left, width = 0.1, 0.8
bottom, height = 0.1, 0.6
bottom_h = 0.7#left_h = left + width + 0.02

rect1 = [left, bottom, width, height]
rect2 = [left, bottom_h, width, 0.2]

spl = UnivariateSpline(mag_outs,comps,k=3,s=0)
xx = np.linspace(np.min(mag_outs)-0.1,np.max(mag_outs)+0.1,1e4)

comp_mags = spl(xx)
# print (np.abs(comp_mags-0.9).argmin())
comp100 = comp_mags[(np.abs(comp_mags-0.995).argmin())]
mag100 = xx[(np.abs(comp_mags-0.995).argmin())]
comp90 = comp_mags[(np.abs(comp_mags-0.9).argmin())]
mag90 = xx[(np.abs(comp_mags-0.9).argmin())]
comp50 = comp_mags[(np.abs(comp_mags-0.5).argmin())]
mag50 = xx[(np.abs(comp_mags-0.5).argmin())]

print mag100, mag90, mag50	
		
plt.figure(figsize=(14,10))

ax1 = plt.axes(rect1)
ax2 = plt.axes(rect2)

ax2.xaxis.set_major_formatter(nullfmt)
if do_tile == '1':
  if do_filter == 'u':
	  xlim1 = [19.,26.]
  if do_filter == 'g':
	  xlim1 = [16.,24.5]
  if do_filter == 'r':
	  xlim1 = [16.,24.]
  if do_filter == 'i':
	  xlim1 = [15.,25.]
  if do_filter == 'z':
	  xlim1 = [15.,24.]
if do_tile == '2':
  if do_filter == 'u':
	  xlim1 = [19.,25.5]
  if do_filter == 'g':
	  xlim1 = [16.,24.5]
  if do_filter == 'r':
	  xlim1 = [16.,24.]
  if do_filter == 'i':
	  xlim1 = [16.,23.5]
  if do_filter == 'z':
	  xlim1 = [16.5,23.]


ylim1 = [-0.05,1.05]

rect50 = draw_rect([xlim1[0],xlim1[0],mag50,mag50],[0.,comp50,comp50,0.], 'dodgerblue', 1.0)
rect90 = draw_rect([xlim1[0],xlim1[0],mag90,mag90],[0.,comp90,comp90,0.], 'blue', 1.0)
rect100 = draw_rect([xlim1[0],xlim1[0],mag100,mag100],[0.,1.,1.,0.], 'midnightblue', 1.)
# patches.append(rect50)
# collection = PatchCollection
ax1.add_patch(rect50)
ax1.add_patch(rect90)
ax1.add_patch(rect100)
# plt.show()

ax1.plot(xlim1,[1.0,1.0],'k-',lw=2)
ax1.plot([mag100,mag100],[ylim1[0],comp100],'k-',lw=2)
ax1.plot([np.min(mag_outs)-0.1,np.max(mag_outs)+0.1],[0.9,0.9],'k--',lw=2)
ax1.plot([mag90,mag90],[ylim1[0],comp90],'k--',lw=2)
ax1.plot([np.min(mag_outs)-0.1,np.max(mag_outs)+0.1],[0.5,0.5],'k--',lw=2)
ax1.plot([mag50,mag50],[ylim1[0],comp50],'k--',lw=2)
ax1.plot([np.min(mag_outs)-0.1,np.max(mag_outs)+0.1],[0.,0.],'k-',lw=2)
ax1.plot(mag_outs,comps,'b--')
ax1.plot(xx,spl(xx),'r-',lw=2)
ax1.set_xlim(xlim1)
ax1.set_ylim(ylim1)
ax1.set_xlabel(r"$%s'_{\rm input}$ [mag]" % do_filter,fontsize=20)
ax1.set_ylabel(r"Completeness", fontsize=20)
plt.setp(ax1.get_xticklabels(), fontsize=14)
ax1.xaxis.set_tick_params(length=8,width=2)
ax1.yaxis.set_tick_params(length=8,width=2)
plt.setp(ax1.get_yticklabels(), fontsize=14)

ax1.annotate(r"$m_{%s', 90\%%} = %.2f$ mag" % (do_filter,mag90),xy=(mag90,comp90+0.01),fontsize=16)
ax1.annotate(r"$m_{%s', 50\%%} = %.2f$ mag" % (do_filter,mag50),xy=(mag50,comp50+0.01),fontsize=16)

ax2.plot([np.min(mag_outs)-0.1,np.max(mag_outs)+0.1],[0.,0.],'k--')
ax2.plot(mock_out_data['MAG'],mock_out_data['MAG_AUTO']-mock_out_data['MAG'],',')
# ax2.plot(mag_outs,diffs,'r',lw=2)
ax2.plot(grid, quad_fit(grid), 'r', lw=2)
for ii in range(len(mag_outs)):
# 	ax2.plot([mag_outs[ii]],[diffs[ii]],'ro',mec='r',ms=2)
# 	ax2.plot([mag_outs[ii],mag_outs[ii]],[diffs[ii]-diff_errs[ii],diffs[ii]+diff_errs[ii]],'r-',lw=2)
	ax2.plot([mag_outs[ii],mag_outs[ii]],quad_fit(mag_outs[ii])+[-1.*boots_out[ii],boots_out[ii]],'r-',lw=2)
ax2.set_xlim(xlim1)
ax2.set_ylim(-0.5,0.5)
ax2.set_ylabel(r"$%s'_{\rm output}-%s'_{\rm input}$" % (do_filter,do_filter), fontsize=20)
plt.setp(ax2.get_xticklabels(), fontsize=14)
ax2.xaxis.set_tick_params(length=8,width=2)
ax2.yaxis.set_tick_params(length=8,width=2)
plt.setp(ax2.get_yticklabels(), fontsize=14)

if do_tile == '1':
	plt.title(r"Tile %s" % do_tile, fontsize=18)
if do_tile != '1':
	plt.title(r"Outer Ring", fontsize=18)
plt.savefig("tile%s_filter%s_complete.png" % (do_tile, do_filter),bbox_inches='tight')

###Split the area into bins, and calculate the 50 and 90% completeness magnitudes in each bin.
#Interpolate between the points to draw the completeness map

nbins_x_50 = 25
nbins_y_50 = 25
nbins_x_90 = 20
nbins_y_90 = 20

mag50_out = comp_mag_calc(.5,nbins_x_50,nbins_y_50,0.25)
mag90_out = comp_mag_calc(.9,nbins_x_90,nbins_y_90,0.5)

xmax = np.max([np.max(mock_in_data['X']),np.max(mock_out_data['X'])])
ymax = np.max([np.max(mock_in_data['Y']),np.max(mock_out_data['Y'])])

if do_tile == '1':
  if do_filter == 'u':
	  xlim2 = [16.,25.5]
  if do_filter == 'g':
	  xlim2 = [15.,24.]
  if do_filter == 'r':
	  xlim2 = [15.,23.5]
  if do_filter == 'i':
	  xlim2 = [15.,24.]
  if do_filter == 'z':
	  xlim2 = [15.,23.5]
if do_tile == '2':
  if do_filter == 'u':
	  xlim2 = [16.,25.5]
  if do_filter == 'g':
	  xlim2 = [15.,24.]
  if do_filter == 'r':
	  xlim2 = [15.,24.]
  if do_filter == 'i':
	  xlim2 = [16.,23.5]
  if do_filter == 'z':
	  xlim2 = [16.,23.]

x_chk = [xlim1[1]-4.,xlim1[1]]
# mask = np.logical_and(mag50_out >= xlim1[0], mag50_out <= xlim1[1])
mask = np.logical_and(mag50_out >= x_chk[0], mag50_out <= x_chk[1])
mag50_mean = np.median(mag50_out[mask])
mag50_std = np.std(mag50_out[mask])
# mask = np.logical_and(mag90_out >= xlim1[0], mag90_out <= xlim1[1])
mask = np.logical_and(mag90_out >= x_chk[0], mag90_out <= x_chk[1])
mag90_mean = np.median(mag90_out[mask])
mag90_std = np.std(mag90_out[mask])

cm = plt.cm.get_cmap('coolwarm')
cm.set_bad('grey',1.)

plot_comp_mags(mag50_out,'50',nbins_x_50,nbins_y_50,mag50_mean,mag50_std)
plot_comp_mags(mag90_out,'90',nbins_x_90,nbins_y_90,mag90_mean,mag90_std)
plt.show()
