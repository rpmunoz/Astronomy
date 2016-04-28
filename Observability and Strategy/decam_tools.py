# Licensed under a 3-clause BSD style license - see LICENSE.rst
# DECam tools are being developed by Roberto Pablo Munoz, PhD

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from sklearn.cluster import AffinityPropagation, MeanShift
from scipy.spatial.distance import pdist
from astroquery.simbad import Simbad

def plot_decam_footprint(cat_data, ra=0., dec=0., target='', filter='', title='DECam'):

	if cat_data.__class__.__name__ == 'Table':
		if filter!='':
			gv=(cat_data['filter']==filter)
			if np.sum(gv)==0:
				print 'There are no observations in the filter ', filter
				return
		else:
			gv=np.arange(len(cat_data))

		ra=np.asarray(cat_data['ra'][gv])
		dec=np.asarray(cat_data['dec'][gv])
		program=cat_data['dtpi'][gv]+' '+cat_data['dtpropid'][gv]
		exptime=np.asarray(cat_data['exposure'][gv])
		title=target+' '+filter
	else:
		ra=np.asarray(ra)
		dec=np.asarray(dec)
		program=''

	decam_corners=np.array([[1652.0,-2994.2],[1652.0,-3521.1],[585.1,-3521.1], #shape 1
						[585.1,-2994.2],[-583.4,-2991.8], #Missing chip South
						[-583.4,-3521.6],[-1651.3,-3521.6],[-1651.3,-2991.8], #shape 2
						[-2210.7,-2930.6],[-2210.7,-2400.3], #shape 6
						[-2769.5,-2339.3],[-2769.5,-1808.5], #shape 11
						[-3322.1,-1746.9],[-3322.1,-1217.0], # shape 17
						[-3329.6,-1156.4],[-3329.6,-626.1], # shape 23
						[-3887.0,-564.0],[-3887.0,-33.4], #shape 30
						[-3886.2,28.5],[-3886.2,557.9], #shape 37
						[-3328.0,621.2],[-3328.0,1150.2], #shape 43
						[-3328.0,1212.0],[-3328.0,1741.0], #shape 49
						[-2769.6,1805.1],[-2769.6,2333.2], #shape 54
						[-2210.4,2397.3],[-2210.4,2925.0], #shape 58
						[-1651.7,2989.0],[-1651.7,3516.3],[-583.1,3516.3], #shape 60
						[-583.1,2989.0],[583.8,2987.6], # Missing chip North
						[583.8,3517.3],[1652.8,3517.3],[1652.8,2987.6], #shape 59
						[2211.8,2926.8],[2211.8,2396.0], #shape 55
						[2770.8,2335.2],[2770.8,1803.6], #shape 50
						[3330.1,1743.0],[3330.1,1212.0], #shape 44
						[3329.7,1151.8],[3329.7,620.6], #shape 38
						[3887.7,558.8],[3887.7,29.4], #shape 31
						[3887.1,-32.9],[3887.1,-562.9], #shape 24
						[3328.7,-626.1],[3328.7,-1155.1], #shape 18
						[3328.0,-1217.6],[3328.0,-1745.7], #shape 12
						[2769.3,-1810.1],[2769.3,-2336.9], #shape 7
						[2210.3,-2401.9],[2210.3,-2928.9], #shape 3
						[1652.0,-2994.2],[1652.0,-3521.1] #shape 1
						])
	decam_corners[:,0]=decam_corners[:,0]*-1
	for i in range(decam_corners.shape[0]-1):
		if np.abs(decam_corners[i,1]-decam_corners[i+1,1])<100:
			decam_corners[i:i+2,1]=np.mean(decam_corners[i:i+2,1])
			i=i+1

	
	program_uniq=np.unique(program)
#	color_uniq= ['red','green','blue','yellow','cyan','magenta']
	color_uniq= plt.cm.jet(np.linspace(0,1,len(program_uniq)))

#	print 'List of PI and Observing programs'
#	print program_uniq

	fig, ax = plt.subplots(figsize=(8,8))

	program_legend=[]
	program_label=[]
	for i in range(len(program_uniq)):
		gv=(program==program_uniq[i])
		x=ra[gv]
		y=dec[gv]
		t=exptime[gv]
		label=program_uniq[i]
		color=color_uniq[i]
#		text_coo=(np.median(ra[gv]),np.median(dec[gv]))
#		ax.text(text_coo[0], text_coo[1], text, horizontalalignment='left')

		patches=[]
		for j in range(len(x)):
			corners=decam_corners/3600. # Convert from arcseconds to degrees
			corners[:,0]=x[j]+corners[:,0]/np.cos(np.radians(y[j]))
			corners[:,1]=y[j]+corners[:,1]
			polygon=Polygon(corners,closed=True)
			patches.append(polygon)
	
		p = PatchCollection(patches, alpha=1./len(x), color=color)
		ax.add_collection(p)

		program_legend.append(Line2D([0], [0], linestyle="none", marker="o", alpha=0.4, markersize=10, markerfacecolor=color))
		program_label.append(label)

		coo=np.concatenate((np.transpose([x]), np.transpose([y])), axis=1)
#		print label, coo, len(coo)
		if len(coo)>1:
			init = np.max(pdist(coo))
#			print 'Affinity init value ', init 
#			af = AffinityPropagation(affinity="euclidean", damping=0.95).fit(coo)
#			af = AffinityPropagation(preference=-1*init).fit(coo)
#			cluster_centers_indices = af.cluster_centers_indices_
#			labels = af.labels_

			ms = MeanShift(bandwidth=0.5, bin_seeding=True)
			ms.fit(coo)
			ms_labels = ms.labels_
			n_tiles = len(np.unique(ms_labels))

			print 'Program: ', label, ' - N of tiles: ', n_tiles
			for i in range(n_tiles):
				class_members = (ms_labels == i)
				tile_ra=np.median(coo[class_members,0])
				tile_dec=np.median(coo[class_members,1])
				tile_exptime=np.sum(t[class_members])
				tile_n_images=np.sum(class_members)
				print("Tile {:d} - Coo: {:.4f}  {:.4f} - N images: {:d} - Total Exposure: {:.1f} s".format(*(i+1,tile_ra,tile_dec, tile_n_images, tile_exptime)))
		else:
			tile_ra=x[0]
			tile_dec=y[0]
			tile_exptime=t[0]
#			print tile_ra, tile_dec, tile_exptime
			print 'Program: ', label, ' - N of tiles: ', 1
			print("Tile {:d} - Coo: {:.4f}  {:.4f} - Exposure: {:.1f} s".format(*(1,tile_ra,tile_dec,tile_exptime)))
		print '---'

	ax.legend(program_legend, program_label, numpoints=1, loc="center left", bbox_to_anchor=(1., 0.5))

	delta_x=np.abs( max(ra)-min(ra) + 2*1.2) # /np.cos(np.radians(np.mean(dec))) )
	delta_y=np.abs( max(dec)-min(dec) + 2*1.2 )
	delta_xy = np.max([delta_x,delta_y])
	plot_xrange=list( np.mean([min(ra),max(ra)]) + delta_xy/(2*np.cos(np.radians(np.mean(dec))))*np.array([1,-1]) )
	plot_yrange=list( np.mean([min(dec),max(dec)]) + delta_xy/2*np.array([-1,1]) )

#	plot_xrange=(min(ra)-1.2/np.cos(np.radians(np.mean(dec))),max(ra)+1.2/np.cos(np.radians(np.mean(dec))))
#	plot_yrange=(min(dec)-1.2, max(dec)+1.2)	

	ax.set_xlim(plot_xrange)
	ax.set_ylim(plot_yrange)
	ax.set_title(title)
	ax.set_xlabel('R.A. (deg)')
	ax.set_ylabel('Dec (deg)')


	if target!='':
		customSimbad = Simbad()
		customSimbad.add_votable_fields('ra(d)','dec(d)','flux(V)','flux(K)')
		result = customSimbad.query_object(target)
		ax.scatter(result['RA_d'][0], result['DEC_d'][0], s=100, c='yellow', marker="*")
		ax.text(result['RA_d'][0], result['DEC_d'][0]-np.abs(plot_yrange[1]-plot_yrange[0])/50, target, horizontalalignment='center', verticalalignment='top')

	plt.show()
	
