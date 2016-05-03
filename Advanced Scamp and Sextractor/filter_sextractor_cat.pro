pro filter_sextractor_cat, input_dir, input_file, FLUX_KEY=flux_key, FWHM_KEY=fwhm_key, FLAGS_KEY=flags_key, N_CHIP=n_chip

	n_cat_sample=20L
	n_mjd_sample=5L
	flux_key = (n_elements(flux_key) GT 0) ? flux_key : 'FLUX_AUTO'
	fwhm_key = (n_elements(fwhm_key) GT 0) ? fwhm_key : 'FWHM_IMAGE'
	flags_key = (n_elements(flags_key) GT 0) ? flags_key : 'FLAGS'
	n_chip = (n_elements(n_chip) GT 0) ? fix(n_chip) : 4

	readcol, input_file, cat_file, FORMAT='A'
	n_cat_file=n_elements(cat_file)

	create_struct, cat_data, '', ['file_in','file_out','mjd','fwhm','mag_range','fwhm_range'], 'A,A,D,F,F(2),F(2)', dim=n_cat_file
	cat_data.file_in=input_dir+'/'+cat_file
	cat_file_ext=(strsplit(cat_file[0], '.', /extract))[-1]
	cat_data.file_out=input_dir+'/'+repstr(cat_file, '.'+cat_file_ext, '_stars.'+cat_file_ext)

	for i=0, n_elements(cat_data)-1 do begin
		im_h=readfits(cat_data[i].file_in, exten=1)
		im_h=string(reform(im_h, 80, n_elements(im_h)/80))
		cat_data[i].mjd=fxpar(im_h, 'MJD')
	endfor

	mjd_floor=floor(cat_data.mjd)
	mjd_uniq=mjd_floor[uniq(mjd_floor, sort(mjd_floor))]

	cgdelete, /all
	cgwindow, wxsize=800, wysize=600, wxpos=0, wypos=0, wobject=win1
	cgwindow, wxsize=800, wysize=600, wxpos=600, wypos=0, wobject=win2
	for i=0L, n_elements(mjd_uniq)-1 do begin

		print, 'Processing MJD '+string(mjd_uniq[i], FORMAT='(D0.4)')

		gv_mjd_uniq=where( floor(cat_data.mjd) EQ mjd_uniq[i], n_gv_mjd_uniq)
		n_mjd_sample=n_gv_mjd_uniq/10

		mjd_range=floor(min(cat_data[gv_mjd_uniq].mjd)*1d3)/1d3 + (ceil(max(cat_data[gv_mjd_uniq].mjd)*1d3)/1d3-floor(min(cat_data[gv_mjd_uniq].mjd)*1d3)/1d3)*dindgen(n_mjd_sample)/(n_mjd_sample-1)

;		gv_file=(sort(randomu(seed, n_gv_mjd)))[0:n_cat_sample-1]

		for j=0, n_elements(mjd_range)-2 do begin
			gv_mjd=where(cat_data.mjd GE mjd_range[j] AND cat_data.mjd LT mjd_range[j+1], n_gv_mjd)
			print, FORMAT='("Processing MJD range: ", F0.4," ",F0.4)', mjd_range[j], mjd_range[j+1]
			print, 'Number of catalogs in MJD range ', n_gv_mjd

			for k=0L, n_gv_mjd-1 do begin
				cat_sex=mrdfits( cat_data[gv_mjd[k]].file_in, 2, cat_sex_h, COLUMNS=[flux_key, fwhm_key, flags_key], /silent)
				gv_sex=where( cat_sex.(0) GT 0. AND cat_sex.(1) GT 0., n_gv_sex)

				create_struct, temp_data, '', ['mag','fwhm','flags'], 'F,F,I', dim=n_gv_sex
				temp_data.mag=30.-2.5*alog10( (cat_sex.(0))[gv_sex] )
				temp_data.fwhm=(cat_sex.(1))[gv_sex]
				temp_data.flags=(cat_sex.(2))[gv_sex]

				if k EQ 0 then cat_data_mjd=temp_data else cat_data_mjd=[cat_data_mjd,temp_data]
			endfor

			gv_flags=where(cat_data_mjd.flags LE 3, n_gv_flags)
			gv_sort=gv_flags[sort(cat_data_mjd[gv_flags].mag)]
			gv_mag=gv_sort[n_gv_flags*0.001:n_gv_flags*0.5]

			n_gv_mag=n_elements(gv_mag)
			gv_sort=gv_mag[sort(cat_data_mjd[gv_mag].fwhm)]
			gv_mag_fwhm=gv_sort[n_gv_mag*0.01:n_gv_mag*0.4]

			temp_mean=biweight_mean(cat_data_mjd[gv_mag_fwhm].fwhm, temp_sigma)
			temp_bin= round( (3.5 * temp_sigma) / n_elements(gv_mag_fwhm)^0.3333 * 100)/100.
			temp_yhist=cghistogram(cat_data_mjd[gv_mag_fwhm].fwhm, binsize=temp_bin, locations=temp_xhist)
			temp_xhist+=temp_bin/2.
			temp_max=max(temp_yhist, gv_max)
			fwhm_mean=temp_xhist[gv_max]

			cgset, win1
			cghistoplot, cat_data_mjd[gv_mag_fwhm].fwhm, binsize=temp_bin, color='red', xtitle='FWHM (pixels)', /window
			cgplot, fwhm_mean*[1,1], [0,1e6], line=2, color='blue', /addcmd, /over

			gv_low=where( cat_data_mjd[gv_mag_fwhm].fwhm LE fwhm_mean, n_gv_low)
			gv_sort=gv_mag_fwhm[gv_low[sort(cat_data_mjd[gv_mag_fwhm[gv_low]].fwhm)]]
			gv_mag_fwhm_low=gv_sort[n_gv_low*0.1:n_gv_low-1]
			fwhm_range=fwhm_mean+abs(fwhm_mean - min(cat_data_mjd[gv_mag_fwhm_low].fwhm))*[-1,1]

			gv_low_high=where(cat_data_mjd[gv_mag_fwhm].fwhm GE fwhm_range[0] AND cat_data_mjd[gv_mag_fwhm].fwhm LE fwhm_range[1], n_gv_low_high)
			fwhm_mean=biweight_mean(cat_data_mjd[gv_mag_fwhm[gv_low_high]].fwhm, fwhm_sigma)
			fwhm_sigma >= 0.2 
			fwhm_range= fwhm_mean + 2*fwhm_sigma*[-1,1]

			gv_stars=gv_mag_fwhm[where(cat_data_mjd[gv_mag_fwhm].fwhm GE fwhm_range[0] AND cat_data_mjd[gv_mag_fwhm].fwhm LE fwhm_range[1], n_gv_stars)]
			mag_range=[min(cat_data_mjd[gv_stars].mag), max(cat_data_mjd[gv_stars].mag)]

			plot_xrange=[-0.5,4.*median(cat_data_mjd[gv_flags].fwhm)]
			plot_yrange=[max(cat_data_mjd[gv_flags].mag), min(cat_data_mjd[gv_flags].mag)-1]

			cgset, win2
			cgplot, cat_data_mjd.fwhm, cat_data_mjd.mag, psym=cgsymcat('filled circle'), color='black', symsize=0.5, xrange=plot_xrange, yrange=plot_yrange, /xstyle, /ystyle, /window, xtitle='FWHM (pixels)', ytitle='Mag', title=string(mjd_range[j], FORMAT='(F0.3)')
			cgplot, cat_data_mjd[gv_flags].fwhm, cat_data_mjd[gv_flags].mag, psym=cgsymcat('filled circle'), color='blue', symsize=0.5, /over, /addcmd 
			cgplot, cat_data_mjd[gv_mag_fwhm].fwhm, cat_data_mjd[gv_mag_fwhm].mag, psym=cgsymcat('filled circle'), color='orange', symsize=0.5, /over, /addcmd 
			cgplot, cat_data_mjd[gv_stars].fwhm, cat_data_mjd[gv_stars].mag, psym=cgsymcat('filled circle'), color='red', symsize=0.5, /over, /addcmd 
			cgplot, fwhm_mean*[1,1], [0,1e6], line=2, color='blue', /addcmd, /over
			cgplot, fwhm_range[0]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over
			cgplot, fwhm_range[1]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over

			cat_data[gv_mjd].fwhm=replicate(fwhm_mean, n_gv_mjd)
			cat_data[gv_mjd].fwhm_range=rebin(fwhm_range, [2,n_gv_mjd])
			cat_data[gv_mjd].mag_range=rebin(mag_range, [2,n_gv_mjd])
			wait, 1

		endfor

	endfor

	for i=0L, n_elements(cat_data)-1 do begin
	
		fits_open, cat_data[i].file_in, fcb_in  ; Lee la tabla fits sin intervenerla           
		fits_read, fcb_in, cat_data0, cat_h0, exten=0
		writefits, cat_data[i].file_out, cat_data0, cat_h0

		for j=0L, n_chip-1 do begin
			temp_data=mrdfits(cat_data[i].file_in, 2*(j+1), cat_sex_h, COLUMNS=[flux_key, fwhm_key, flags_key], /silent)
			create_struct, cat_sex, '', ['mag','fwhm','flags'], 'F,F,I', dim=n_elements(temp_data)
			cat_sex.mag=30.-2.5*alog10( temp_data.(0) )
			cat_sex.fwhm=temp_data.(1)
			cat_sex.flags=temp_data.(2)
				
			if j EQ 0 then cat_sex_all=cat_sex else cat_sex_all=[cat_sex_all,cat_sex]
		endfor

		gv_stars=where(cat_sex_all.mag GT cat_data[i].mag_range[0] AND cat_sex_all.mag LT cat_data[i].mag_range[1] AND cat_sex_all.fwhm GT cat_data[i].fwhm_range[0] AND cat_sex_all.fwhm LT cat_data[i].fwhm_range[1] AND cat_sex_all.flags LE 3, n_gv_stars)
		fwhm_mean=biweight_mean(cat_sex_all[gv_stars].fwhm, fwhm_sigma)
		cat_data[i].fwhm_range = fwhm_mean + abs(cat_data[i].fwhm_range[1]-cat_data[i].fwhm_range[0])/2.*[-1,1]
;		gv_stars=where(cat_sex_all.mag GT cat_data[i].mag_range[0] AND cat_sex_all.mag LT cat_data[i].mag_range[1] AND cat_sex_all.fwhm GT cat_data[i].fwhm_range[0] AND cat_sex_all.fwhm LT cat_data[i].fwhm_range[1] AND cat_sex_all.flags LE 3, n_gv_stars)

		for j=0L, n_chip-1 do begin
			temp_data=mrdfits(cat_data[i].file_in, 2*(j+1), cat_sex_h, COLUMNS=[flux_key, fwhm_key, flags_key], /silent)
			create_struct, cat_sex, '', ['mag','fwhm','flags'], 'F,F,I', dim=n_elements(temp_data)
			cat_sex.mag=30.-2.5*alog10( temp_data.(0) )
			cat_sex.fwhm=temp_data.(1)
			cat_sex.flags=temp_data.(2)

			gv_stars=where(cat_sex.mag GT cat_data[i].mag_range[0] AND cat_sex.mag LT cat_data[i].mag_range[1] AND cat_sex.fwhm GT cat_data[i].fwhm_range[0] AND cat_sex.fwhm LT cat_data[i].fwhm_range[1] AND cat_sex.flags LE 3, n_gv_stars)
;			fwhm_mean=biweight_mean(cat_sex[gv_stars].fwhm, fwhm_sigma)
;			cat_data[i].fwhm_range = fwhm_mean + abs(cat_data[i].fwhm_range[1]-cat_data[i].fwhm_range[0])/2.*[-1,1]
;			gv_stars=where(cat_sex.mag GT cat_data[i].mag_range[0] AND cat_sex.mag LT cat_data[i].mag_range[1] AND cat_sex.fwhm GT cat_data[i].fwhm_range[0] AND cat_sex.fwhm LT cat_data[i].fwhm_range[1] AND cat_sex.flags LE 3, n_gv_stars)
;			gv_stars=where(cat_sex.flags LE 3, n_gv_stars)

			cgset, win1
			cgplot, cat_sex.fwhm, cat_sex.mag, psym=cgsymcat('filled circle'), color='black', symsize=0.5, xrange=cat_data[i].fwhm_range*[0.5,4.], yrange=reverse(cat_data[i].mag_range)+[2.,-1.], /xstyle, /ystyle, /window, xtitle='FWHM (pixels)', ytitle='Mag'
			cgplot, cat_sex[gv_stars].fwhm, cat_sex[gv_stars].mag, psym=cgsymcat('filled circle'), color='red', /addcmd, /over
			cgplot, cat_data[i].fwhm_range[0]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over
			cgplot, cat_data[i].fwhm_range[1]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over
			wait, 0.5

      fits_read, fcb_in, cat_data1, cat_h1, exten=2*j+1
      fits_read, fcb_in, cat_data2, cat_h2, exten=2*j+2

      cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
      cat_data2=cat_data2[*,[gv_stars]]
      fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]

      fits_open, cat_data[i].file_out, fcb_out, /append
      fits_write, fcb_out, cat_data1, cat_h1
      fits_write, fcb_out, cat_data2, cat_h2
      fits_close, fcb_out

		endfor
		fits_close, fcb_in

	endfor

end
