template_cat_file='reference_catalog/2MASS_0338-3531_r70.cat'

input_cat_file='reference_catalog/NGVS-2-5.I2_ascii.cat'
output_cat_file='reference_catalog/NGVS-2-5.I2_scamp.cat'

cgdelete, /all
readcol, input_cat_file, temp_ra, temp_dec, temp_mag_auto, temp_magerr_auto, temp_erra_world, temp_errb_world, temp_flux_radius, temp_flags, FORMAT='X,X,X,F,F,F,F,X,X,X,X,X,F,X,F,X,X,X,X,F,F'

gv=where(temp_mag_auto LT 99. AND temp_flags LE 3 AND temp_flux_radius GT 0, n_gv)
plot_xrange=[-0.5,4.*median(temp_flux_radius[gv])]
plot_yrange=[max(temp_mag_auto[gv]), min(temp_mag_auto[gv])]

;cgdisplay, 800, 600, wid=0
cgwindow, wxsize=800, wysize=600, wobject=win1
cgplot, temp_flux_radius[gv], temp_mag_auto[gv], psym=cgsymcat('filled circle'), symsize=0.4, xrange=plot_xrange, yrange=plot_yrange, /xstyle, /ystyle, /window 

n_gv=n_elements(gv)
gv_sort=gv[sort(temp_mag_auto[gv])]
gv_mag=gv_sort[n_gv*0.001:n_gv*0.5]

n_gv_mag=n_elements(gv_mag)
gv_sort=gv_mag[sort(temp_flux_radius[gv_mag])]
gv_mag_radius=gv_sort[n_gv_mag*0.01:n_gv_mag*0.5]

cgplot, temp_flux_radius[gv_mag_radius], temp_mag_auto[gv_mag_radius], psym=cgsymcat('filled circle'), symsize=0.4, color='red', /addcmd, /over

temp_mean=biweight_mean(temp_flux_radius[gv_mag_radius], temp_sigma)
temp_bin= round( (3.5 * temp_sigma) / n_Elements(gv_mag_radius)^0.3333 * 100)/100.
temp_yhist=cghistogram(temp_flux_radius[gv_mag_radius], binsize=temp_bin, locations=temp_xhist)
temp_xhist+=temp_bin/2.
temp_max=max(temp_yhist, gv_max)
flux_radius_mean=temp_xhist[gv_max]

cgwindow, wxsize=800, wysize=600, wobject=win2
cghistoplot, temp_flux_radius[gv_mag_radius], binsize=temp_bin, color='red', xtitle='Flux radius', /window
cgplot, flux_radius_mean*[1,1], [0,1e6], line=2, color='blue', /addcmd, /over

cgset, win1
cgplot, flux_radius_mean*[1,1], [0,1e6], line=2, color='blue', /addcmd, /over
gv_low=where( temp_flux_radius[gv_mag_radius] LE flux_radius_mean, n_gv_low)
gv_sort=gv_mag_radius[gv_low[sort(temp_flux_radius[gv_mag_radius[gv_low]])]]
gv_mag_radius_low=gv_sort[n_gv_low*0.05:n_gv_low-1]
flux_radius_range=flux_radius_mean+abs(flux_radius_mean - min(temp_flux_radius[gv_mag_radius_low]))*[-1,1]
cgplot, flux_radius_range[0]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over
cgplot, flux_radius_range[1]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over

gv_low_high=where(temp_flux_radius[gv_mag_radius] GE flux_radius_range[0] AND temp_flux_radius[gv_mag_radius] LE flux_radius_range[1], n_gv_low_high)
flux_radius_mean=biweight_mean(temp_flux_radius[gv_mag_radius[gv_low_high]], flux_radius_sigma)

flux_radius_range= flux_radius_mean + 2*flux_radius_sigma*[-1,1]
gv_stars=gv_mag_radius[where(temp_flux_radius[gv_mag_radius] GE flux_radius_range[0] AND temp_flux_radius[gv_mag_radius] LE flux_radius_range[1], n_gv_stars)]
cgplot, temp_flux_radius[gv_stars], temp_mag_auto[gv_stars], psym=cgsymcat('filled circle'), symsize=0.4, color='green', /addcmd, /over

create_struct, cat_data, '', ['X_WORLD','Y_WORLD','ERRA_WORLD','ERRB_WORLD','MAG','MAGERR','OBSDATE'], 'D,D,E,E,E,E,D', dim=n_gv_stars
cat_data.x_world=temp_ra[gv_stars]
cat_data.y_world=temp_dec[gv_stars]
cat_data.erra_world=temp_erra_world[gv_stars]
cat_data.errb_world=temp_errb_world[gv_stars]
cat_data.mag=temp_mag_auto[gv_stars]
cat_data.magerr=temp_magerr_auto[gv_stars]
cat_data.obsdate=replicate(2000., n_elements(cat_data))

;mwrfits, input_cat, 'test.fits'
;fits_open, 'test.fits', cat_fcb
;fits_read, cat_fcb, cat_data0, cat_h0, exten=0
;fits_read, cat_fcb, cat_data1, cat_h1, exten=1

fits_open, template_cat_file, cat_fcb
fits_read, cat_fcb, cat_data0, cat_h0, exten=0
fits_read, cat_fcb, cat_data1, cat_h1, exten=1
fits_read, cat_fcb, cat_data2, cat_h2, exten=2
fits_close, cat_fcb

fxaddpar, cat_h2, 'NAXIS2', n_elements(cat_data)
cat_data1=reform(cat_data1, [n_elements(cat_data1),1])

fits_open, output_cat_file, fcb_out, /write
fits_write, fcb_out, cat_data0, cat_h0
fits_write, fcb_out, cat_data1, cat_h1
fits_close, fcb_out

fxbcreate, table_u, output_cat_file, cat_h2
for j=0L, n_tags(cat_data)-1 do begin
	fxbwritm, table_u, j+1, cat_data.(j)
endfor
fxbfinish, table_u

END
