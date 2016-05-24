pro decam_completeness

cd, '/Volumes/Q6/matt/2014A-0610/completeness' ;'/Users/rmunoz/Research/ngvs/completeness/code'
im_dir='/Volumes/Q6/matt/stacks'
sex_dir='/Volumes/Q6/matt/stacks/completeness/results/sex' ;'/Volumes/RAID/Research/NGVS/completeness/results/sex'
psfex_dir='/Volumes/Q6/matt/stacks';completeness/results/psfex' ;'/Volumes/RAID/Research/NGVS/completeness/results/psfex'
mock_dir='/Volumes/Q6/matt/stacks/completeness/mock'
do_debug=0

if file_test(sex_dir, /dir) EQ 0 then file_mkdir, sex_dir, /noexpand
if file_test(mock_dir, /dir) EQ 0 then file_mkdir, mock_dir, /noexpand

do_tile = '2'
tile_suffix=[do_tile]
do_filter = 'z'

if do_filter eq 'i' then begin
	psf_file=psfex_dir+'/survey_tile'+do_tile+'_'+do_filter+'_psfex.psf' ;Use the previously generated catalogues
	im_file=im_dir+'/survey_tile'+do_tile+'_'+do_filter+'_short.fits'
	im_wgt_file=im_dir+'/survey_tile'+do_tile+'_'+do_filter+'_short.WEIGHT.fits'
endif else begin
	psf_file=psfex_dir+'/survey_tile'+do_tile+'_'+do_filter+'_psfex_ALIGNi.psf' ;Use the previously generated catalogues
	im_file=im_dir+'/survey_tile'+do_tile+'_'+do_filter+'_short_ALIGNi.fits'
	im_wgt_file=im_dir+'/survey_tile'+do_tile+'_'+do_filter+'_short_ALIGNi.WEIGHT.fits'
endelse
;im_file=im_dir+'/survey_tile'+do_tile+'_'+do_filter+'_short.fits'
;im_wgt_file=im_dir+'/survey_tile'+do_tile+'_'+do_filter+'_short.WEIGHT.fits'
;im_flag_file='/Volumes/RAID/Data/NGVS/WIRCAM/v7/mosaics/'+'ngvs_'+tile_suffix+'_Ks_image_v7.MEGACAM_MEGACAM_BEST.SIGWEIGHTED_LANCZOS2.FLAG.fits'

im_cat_file=sex_dir+'/survey_tile'+do_tile+'_'+do_filter+'_psf.ldac' ;'/ngvs_cat_sex_'+tile_suffix+'.ldac' ;Used to avoid crowding issues
im_check_file=sex_dir+'/check/survey_tile'+do_tile+'_'+do_filter+'.CHECK_SEGMENTATION.fits'   ;'/ngvs_check_seg_sex_'+tile_suffix+'.fits'

;psf_mock_mag_file=psfex_dir+'/survey_tile'+do_tile+'_'+do_filter+'_psf.ldac' ;
;mock_mags = readfits(psf_mock_mag_file, h, /EXTEN)
;fxbopen, mock_mag_tab, psf_mock_mag_file, 2, mag_h
;fxbread, mock_mag_tab, mock_mag_data, 'MAG_PSF'
;mock_mag_limits=[min(mock_mag_data),max(mock_mag_data[where(mock_mag_data NE 99.,n_gv)])]

mock_n=10L ;number of realizations
;mock_mag_limits=[min(mock_mag_data),max(mock_mag_data[where(mock_mag_data NE 99.,n_gv)])] ; magnitude range
mock_mag_limits=[15.,30.] ; magnitude range, adjust the bright range if it's complete
mock_nstar=4L ; number of stars per cell (mesh)
mock_mesh=[200.,200.] ; Mesh size in pixels
;ab_vega=1.827
;
mock_coo_file=mock_dir+'/mock_coo_'+do_tile+'_'+do_filter+'.dat';coordinate file for mock sources
a=byte(tile_suffix)
a_size=size(a,/dim)
b=byte(repstr(string(indgen(mock_n)+1,FORMAT='(I4)'),' ','0'))
b_size=size(b,/dim)
;print, a, a_size, b, b_size
;stop
mock_mag_file=reform(mock_dir+'/mock_mag_t'+do_tile+'_'+do_filter+'_'+string(rebin(reform(a, [a_size,1]),[a_size,b_size[1]]))+'_'+string(rebin(reform(b, [b_size[0],1,b_size[1]]),[b_size[0],a_size[0],b_size[1]]))+'.dat', [n_elements(im_file),mock_n])
mock_im_file=reform(mock_dir+'/mock_im_t'+do_tile+'_'+do_filter+'_'+string(rebin(reform(a, [a_size,1]),[a_size,b_size[1]]))+'_'+string(rebin(reform(b, [b_size[0],1,b_size[1]]),[b_size[0],a_size[0],b_size[1]]))+'.fits', [n_elements(im_file),mock_n])
sex_cat_file=reform(mock_dir+'/sex_cat_t'+do_tile+'_'+do_filter+'_'+string(rebin(reform(a, [a_size,1]),[a_size,b_size[1]]))+'_'+string(rebin(reform(b, [b_size[0],1,b_size[1]]),[b_size[0],a_size[0],b_size[1]]))+'.ldac', [n_elements(im_file),mock_n])

mock_completeness_file = mock_dir+'/mock_completeness_t'+do_tile+'_'+do_filter+'.dat'

im_mask_region=[{tile:'t+0+0', xrange:[15535L,17759L], yrange:[16623L,19383L]}]
im_blank_region=[{tile:'t-1+0', xrange:[9000L,-1], yrange:[0L,5850L]}]
;
for i=0L, n_elements(im_file)-1 do begin
;  goto, mock_plots
;  goto, mock_psf
;  goto, mock_coo

  fxbopen, table_u, psf_file[i], 1, psf_h
  fxbread, table_u, psf_data, 1
  free_lun, table_u

  psf_order=long(fxpar(psf_h, 'POLDEG1'))
  psf_offset=double( [fxpar(psf_h, 'POLZERO1'),fxpar(psf_h, 'POLZERO2')] )
  psf_scale=double( [fxpar(psf_h, 'POLSCAL1'),fxpar(psf_h, 'POLSCAL2')] )
  psf_pixstep=double(fxpar(psf_h, 'PSF_SAMP'))
  psf_size=size(psf_data, /dim)

;  command='sex '+im_file[i]+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+im_cat_file[i]+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+im_wgt_file[i]+' -WRITE_XML N -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+ im_check_file[i] +' -SEEING_FWHM 1.00 -BACK_SIZE 256 -BACK_FILTERSIZE 3 -DETECT_THRESH 2. -ANALYSIS_THRESH 3.' ;-FLAG_TYPE OR -FLAG_IMAGE '+im_flag_file[i]
  command='sex '+im_file[i]+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+im_cat_file[i]+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+im_wgt_file[i]+' -WRITE_XML N -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+ im_check_file[i] +' -SEEING_FWHM 1.00 -BACK_SIZE 256 -BACK_FILTERSIZE 3 -DETECT_THRESH 2. -ANALYSIS_THRESH 3.' ;-FLAG_TYPE OR -FLAG_IMAGE '+im_flag_file[i]
  print, command
  spawn, command
;
  mock_coo:
;Generates valid coordinates (voids)
  im_check_data=(readfits_big(im_check_file[i], im_h) EQ 0)
  im_wgt_data=(readfits_big(im_wgt_file[i], im_h) GT 0)
  
  im_h=headfits(im_file[i])
  extast, im_h, im_ast
  im_size=im_ast.naxis

  mock_xrange=(0.+findgen(ceil(im_size[0]/mock_mesh[0]))*mock_mesh[0])<im_size[0]
  mock_yrange=(0.+findgen(ceil(im_size[1]/mock_mesh[1]))*mock_mesh[1])<im_size[1]
  
  mock_xy=lonarr(2,n_elements(mock_xrange)*n_elements(mock_yrange)*mock_nstar)
  mock_xy[*]=-1
  mock_id=0L
	mock_psf_radius=22
	mock_check_radius=10
	mock_star_radius=16
  for j=0L, n_elements(mock_xrange)-1 do begin
 		print, 'Processing coo iteration ', strn(j), ' of ', strn(n_elements(mock_xrange))

    for k=0L, n_elements(mock_yrange)-1 do begin
      x=(mock_xrange[j] + randomu(seed, mock_nstar*10)*(mock_mesh[0]-1)) < (im_size[0]-mock_psf_radius-1) > mock_psf_radius
      y=(mock_yrange[k] + randomu(seed, mock_nstar*10)*(mock_mesh[1]-1)) < (im_size[1]-mock_psf_radius-1) > mock_psf_radius

      l=0L
      temp_n=0L
      repeat begin
        mock_range=[(mock_id-(n_elements(mock_yrange)+2)*mock_nstar)>0L,(mock_id-1)>0]
        do_star = ( total(im_check_data[x[l]-mock_check_radius:x[l]+mock_check_radius,y[l]-mock_check_radius:y[l]+mock_check_radius] EQ 0) EQ 0 AND total(im_wgt_data[x[l]-mock_check_radius:x[l]+mock_check_radius,y[l]-mock_check_radius:y[l]+mock_check_radius] EQ 0) EQ 0 AND total(sqrt(total((mock_xy[*,mock_range[0]:mock_range[1]]-[x[l],y[l]]#make_array(mock_range[1]-mock_range[0]+1, value=1., /float))^2,1)) LT 2.*mock_star_radius) EQ 0 )
        if do_star then begin
          mock_xy[*,mock_id]=round([x[l],y[l]])
          mock_id++
          temp_n++
        endif
        l++
      endrep until temp_n GE mock_nstar OR l GE mock_nstar*10  
    endfor
  endfor
  
  gv=where(mock_xy[0,*] GT 0, n_gv)
  if n_gv GT 0 then begin
    mock_xy=mock_xy[*,gv]

    xy2ad, mock_xy[0,*], mock_xy[1,*], im_ast, temp_ra, temp_dec
    openw, lun, mock_coo_file[i], /get_lun
    printf, lun, '#    X      Y     ALPHA_J2000     DELTA_J2000'
    for j=0L, n_gv-1 do begin
      printf, lun, mock_xy[0,j], mock_xy[1,j], temp_ra[j], temp_dec[j], FORMAT='(I8,4X,I8,4X,F10.6,4X,F10.6)'
    endfor
    free_lun, lun
  endif $
  else begin
		print, 'ERROR'
		exit
	endelse

  mock_mag:
;Generate magnitudes at the coordinates generated above.
  readcol, mock_coo_file[i], mock_x, mock_y, FORMAT='I,I,F,F'
  n_mock_x=n_elements(mock_x)
  
  for j=0L, mock_n-1 do begin
 		print, 'Processing mag iteration ', strn(j), ' of ', strn(mock_n)

	  mock_mag=mock_mag_limits[0]+randomu(seed, n_mock_x)*(mock_mag_limits[1]-mock_mag_limits[0])
    openw, lun, mock_mag_file[i,j], /get_lun
    printf, lun, '#    X      Y     MAG'
		printf, lun, string(transpose(mock_x),FORMAT='(I8)')+'  '+string(transpose(mock_y),FORMAT='(I8)')+'  '+string(transpose(mock_mag),FORMAT='(F8.4)'), FORMAT='(A)'
;    for k=0L, n_mock_x-1 do begin
;      printf, lun, mock_x[k], mock_y[k], mock_mag[k], FORMAT='(I8,4X,I8,4X,F8.4)'
;    endfor
    free_lun, lun
  endfor 

 
;	continue
  mock_psf:
;Generate sources based on the psf catalogues, call the python script for multiprocessing
  forprint, replicate(im_file[i], mock_n), replicate(im_wgt_file[i], mock_n), replicate(psf_file[i], mock_n), mock_mag_file[i,*], mock_im_file[i,*], sex_cat_file[i,*], FORMAT='(A,4X,A,4X,A,4X,A,4X,A,4X,A,4X,A)', textout=mock_completeness_file[i], COMMENT='#    im_file   weight_file   flag_file   psf_file   mock_mag_file   mock_im_file   sex_cat_file', /SILENT

;	gv=90+indgen(10)
;	n_gv=n_elements(gv)
;  forprint, replicate(im_file[i], n_gv), replicate(im_wgt_file[i], n_gv), replicate(im_flag_file[i], n_gv), replicate(psf_file[i], n_gv), mock_mag_file[i,gv], mock_im_file[i,gv], sex_cat_file[i,gv], FORMAT='(A,4X,A,4X,A,4X,A,4X,A,4X,A,4X,A)', textout=mock_completeness_file[i], COMMENT='#    im_file   weight_file   flag_file   psf_file   mock_mag_file   mock_im_file   sex_cat_file', /SILENT

;	command='python decam_completeness.py '+mock_completeness_file[i]
;	print, command
;	spawn, command

;stop
;	continue

   	for j=0L, mock_n-1 do begin
		print, 'Processing mock_mag_file number '+strn(j+1)
		command='python decam_completeness.py '+im_file[i]+' '+im_wgt_file[i]+' '+psf_file[i]+' '+mock_mag_file[i,j]+' '+mock_im_file[i,j]+' '+im_cat_file[i]
		print, command
		spawn, command
	endfor
	stop

	mock_plots:

	mock_in_list_x=list()
	mock_in_list_y=list()
	mock_in_list_mag=list()
	mock_in_mag_match=list()

	mock_out_list_x=list()
	mock_out_list_y=list()
	mock_out_list_mag=list()
 
	for j=0L, mock_n-1 do begin
		print, 'Processing files: ', mock_mag_file[i,j]
		print, sex_cat_file[i,j]

    readcol, mock_mag_file[i,j], temp_x, temp_y, temp_mag, FORMAT='I,I,F', /silent
		temp_n=n_elements(temp_x)
		create_struct, temp_cat, '', ['x','y','mag'], 'F,F,F', dim=temp_n
		temp_cat.x=temp_x+0.5
		temp_cat.y=temp_y+0.5
		temp_cat.mag=temp_mag
		gv=sort(temp_cat.x)
		mock_in_cat=temp_cat[gv]
		n_mock_in_cat=n_elements(mock_in_cat)
		temp_cat=0 

		fxbopen, table_u, sex_cat_file[i,j], 2, cat_h
;		fxbreadm, table_u, ['X_IMAGE','Y_IMAGE','MAG_APER','MAG_AUTO','FLAGS','IMAFLAGS_ISO'], temp_x, temp_y, temp_mag_aper, temp_mag_auto, temp_flags, temp_imaflags
		fxbreadm, table_u, ['X_IMAGE','Y_IMAGE','MAG_APER','MAG_AUTO'], temp_x, temp_y, temp_mag_aper, temp_mag_auto
		free_lun, table_u
		temp_n=n_elements(temp_x)
		create_struct, temp_cat, '', ['x','y','mag_aper','mag_auto'], 'F,F,F(7),F', dim=temp_n
		temp_cat.x=temp_x
		temp_cat.y=temp_y
		temp_cat.mag_aper=temp_mag_aper
		temp_cat.mag_auto=temp_mag_auto
;		temp_cat.flags=temp_flags
;		temp_cat.imaflags=temp_imaflags
		gv=sort(temp_cat.x)
		mock_out_cat=temp_cat[gv]
		n_mock_out_cat=n_elements(mock_out_cat)
		temp_cat=0 

		if n_mock_out_cat LT n_mock_in_cat*0.1 then begin
			print, 'WARNING! The number of sources in the output catalog is less than 10% of input catalog'
			stop
		endif

		temp_x=0
		temp_y=0
		temp_mag=0
		temp_mag_aper=0
		temp_mag_auto=0
		temp_flags=0
		temp_imaflags=0

		match_n=100L
		mock_in_range=[transpose(floor(findgen(match_n)/match_n*(n_mock_in_cat-1))>0L) , transpose(floor((findgen(match_n)+1)/match_n*(n_mock_in_cat-1))-1)]
		mock_in_range[1,match_n-1]=n_mock_in_cat-1

		mock_out_range=mock_in_range
		for k=0L, match_n-1 do begin
			temp=min(abs(mock_out_cat.x - mock_in_cat[mock_in_range[0,k]].x), gv)
			mock_out_range[0,k]=gv-500
			temp=min(abs(mock_out_cat.x - mock_in_cat[mock_in_range[1,k]].x), gv)
			mock_out_range[1,k]=gv+500
;			if k LT match_n-1 then mock_out_range[0,k+1]=gv-2000
		endfor
		mock_out_range[1,match_n-1]=n_mock_out_cat-1
		mock_out_range= (mock_out_range > 0L < (n_mock_out_cat-1))

		;gv_match= where( (temp=sqrt(min( (mock_in_cat[mock_in_range[0,k]:mock_in_range[1,k]].x#make_array(n_mock_out_cat,value=1.,/double) - make_array(n_mock_in_cat,value=1.,/double)#mock_out_cat.x)^2 + (mock_in_cat.y#make_array(n_mock_out_cat,value=1.,/double) - make_array(n_mock_in_cat,value=1.,/double)#mock_out_cat.y)^2, id_match, dim=2))) LT 1., n_gv_match)
		gv_match_mock_in=list()
		gv_match_mock_out=list()
		for k=0L, match_n-1 do begin
;			print, 'Matching catalogs k='+strn(k)
			n_mock_in_range=mock_in_range[1,k]-mock_in_range[0,k]+1
			n_mock_out_range=mock_out_range[1,k]-mock_out_range[0,k]+1
			
			gv_match= where( (temp=min( (mock_in_cat[mock_in_range[0,k]:mock_in_range[1,k]].x#make_array(n_mock_out_range,value=1.,/float) - make_array(n_mock_in_range,value=1.,/float)#mock_out_cat[mock_out_range[0,k]:mock_out_range[1,k]].x)^2 + (mock_in_cat[mock_in_range[0,k]:mock_in_range[1,k]].y#make_array(n_mock_out_range,value=1.,/float) - make_array(n_mock_in_range,value=1.,/float)#mock_out_cat[mock_out_range[0,k]:mock_out_range[1,k]].y)^2, id_match, dim=2)) LT 4., n_gv_match)
			gv_match_mock_in.add, mock_in_range[0,k] + (id_match[gv_match] mod n_mock_in_range), /extract
			gv_match_mock_out.add,  mock_out_range[0,k] + id_match[gv_match]/n_mock_in_range, /extract
		endfor
		gv_match_mock_in=gv_match_mock_in.toarray(type='long')
		gv_match_mock_out=gv_match_mock_out.toarray(type='long')
		
;		plotsym, 0, 0.4, /fill
;		plot, mock_in_cat[gv_match_mock_in].mag, reform(mock_out_cat[gv_match_mock_out].mag_aper[2,*])-mock_in_cat[gv_match_mock_in].mag, psym=8, yrange=[-0.1,0.1]
;		if j EQ 0 then begin
;			gv=where(mock_in_cat[gv_match_mock_in].mag LT 19., n_gv)
;			mag_offset=median( reform(mock_out_cat[gv_match_mock_out[gv]].mag_aper[2,*])-mock_in_cat[gv_match_mock_in[gv]].mag)
;			print, mag_offset
;		endif

		if do_debug then begin
			grid_xrange=[2706.80, 3045.15]
			grid_yrange=[0.00000, 356.217]

		    gv_in=where( mock_in_cat.x GE grid_xrange[0] AND mock_in_cat.x LT grid_xrange[1] AND mock_in_cat.y GE grid_yrange[0] AND mock_in_cat.y LT grid_yrange[1], n_gv_in)
  		    gv_out=where( mock_in_cat[gv_match_mock_in].x GE grid_xrange[0] AND mock_in_cat[gv_match_mock_in].x LT grid_xrange[1] AND mock_in_cat[gv_match_mock_in].y GE grid_yrange[0] AND mock_in_cat[gv_match_mock_in].y LT grid_yrange[1], n_gv_out)

			xhist_bin=0.4
	    	yhist_in=histogram(mock_in_cat[gv_in].mag, bin=xhist_bin, locations=xhist_in, min=18., max=21.999, /nan)
	    	xhist_in += xhist_bin/2.
 	   		yhist_out=histogram(mock_in_cat[gv_match_mock_in[gv_out]].mag, bin=xhist_bin, locations=xhist_out, min=18., max=21.999, /nan)
	    	xhist_out += xhist_bin/2

	    	temp=yhist_out*1./yhist_in
	 	  	plot, xhist_out+ab_vega, temp, xtitle='K (AB)', ytitle='mock_out/mock_in', xrange=[20.,24.], yrange=[0.,1.2], ystyle=1

			gv1=sort(mock_in_cat[gv_in].mag)
			forprint, mock_in_cat[gv_in[gv1]].x, mock_in_cat[gv_in[gv1]].y, mock_in_cat[gv_in[gv1]].mag, text=1
			gv2=sort(mock_in_cat[gv_match_mock_in[gv_out]].mag)
			print, mock_in_cat[gv_match_mock_in[gv_out[gv2]]].mag
			stop
		endif

		mock_in_list_x.add, round(mock_in_cat.x), /extract
		mock_in_list_y.add, round(mock_in_cat.y), /extract
		mock_in_list_mag.add, mock_in_cat.mag, /extract
		mock_in_mag_match.add, mock_in_cat[gv_match_mock_in].mag, /extract

		mock_out_list_x.add, round(mock_in_cat[gv_match_mock_in].x), /extract
		mock_out_list_y.add, round(mock_in_cat[gv_match_mock_in].y), /extract
		mock_out_list_mag.add, reform(mock_in_cat[gv_match_mock_in].mag), /extract ;mock_in_cat[gv_match_mock_in].mag, /extract
;		mock_out_x.add, round(mock_out_cat[gv_match_mock_out].x), /extract
;		mock_out_y.add, round(mock_out_cat[gv_match_mock_out].y), /extract
;		mock_out_mag.add, reform(mock_out_cat[gv_match_mock_out].mag_aper[2,*]), /extract ;mock_in_cat[gv_match_mock_in].mag, /extract

		mock_in_cat=0
		mock_out_cat=0
		gv_match_mock_in=0
		if (j EQ 0) then begin ;OR ((j+1) mod 50 EQ 0) then begin
			print, 'Plotting iteration j='+strn(j+1)
			xhist_bin=0.1
			yhist_in=histogram(mock_in_list_mag.toarray(type='float'), bin=xhist_bin, locations=xhist_in, min=mock_mag_limits[0], max=mock_mag_limits[1]-1e-4, /nan)
			xhist_in += xhist_bin/2.

			yhist_out=histogram(mock_out_list_mag.toarray(type='float'), bin=xhist_bin, locations=xhist_out, min=mock_mag_limits[0], max=mock_mag_limits[1]-1e-4, /nan)
			xhist_out += xhist_bin/2.

			plot, xhist_out+1.82, yhist_out*1./yhist_in, xrange=[20.,28.], yrange=[0.,1.2]
			oplot, [0,100], 0.9*[1,1], line=2  
			oplot, [0,100], 0.5*[1,1], line=2  
		endif

	endfor

	mock_in_x=mock_in_list_x.toarray(type='long')
	mock_in_list_x.remove, /ALL
	mock_in_y=mock_in_list_y.toarray(type='long')
	mock_in_list_y.remove, /ALL
	mock_in_mag=mock_in_list_mag.toarray(type='float')
	mock_in_list_mag.remove, /ALL
	mock_out_x=mock_out_list_x.toarray(type='long')
	mock_out_list_x.remove, /ALL
	mock_out_y=mock_out_list_y.toarray(type='long')
	mock_out_list_y.remove, /ALL
	mock_out_mag=mock_out_list_mag.toarray(type='float')
	mock_out_list_mag.remove, /ALL
	mock_in_mag_match=mock_in_mag_match.toarray(type='float')

;	save, mock_n, mock_nstar, mock_mesh, mock_in_x, mock_in_y, mock_out_x, mock_out_y, mock_in_mag, mock_out_mag, mock_mag_limits, ab_vega, file='ngvsir_completeness_'+tile_suffix[i]+'_partial.sav'
	save, mock_n, mock_nstar, mock_mesh, mock_in_x, mock_in_y, mock_out_x, mock_out_y, mock_in_mag, mock_out_mag, mock_mag_limits, file='decam_completeness_'+tile_suffix[i]+'_partial.sav'
;	mock_in_x=0
;	mock_in_y=0
;	mock_in_mag=0
;	mock_in_mag_match=0
;	mock_out_x=0
;	mock_out_y=0
;	mock_out_mag=0
;	continue

	gv_in=where(mock_in_x GT 2000L AND mock_in_x LT 18500 AND mock_in_y GT 2000L AND mock_in_y LT 15000L, n_gv_in)
	gv_out=where(mock_out_x GT 2000L AND mock_out_x LT 18500 AND mock_out_y GT 2000L AND mock_out_y LT 15000L, n_gv_out)

	xhist_bin=0.1
	yhist_in=histogram(mock_in_mag[gv_in], bin=xhist_bin, locations=xhist_in, min=mock_mag_limits[0], max=mock_mag_limits[1]-1e-4, /nan)
	xhist_in += xhist_bin/2.
	yhist_out=histogram(mock_out_mag[gv_out], bin=xhist_bin, locations=xhist_out, min=mock_mag_limits[0], max=mock_mag_limits[1]-1e-4, /nan)
	xhist_out += xhist_bin/2.

	temp_x=xhist_out;+ab_vega
	temp_y=yhist_out*1./yhist_in
	gv=where(temp_x GT mock_mag_limits[0]);+ab_vega+0.2)

	w=window(dim=[800,600], margin=[0.2,0.2,0.2,0.2])
;	p=plot(temp_x[gv], temp_y[gv]/max(temp_y[gv]), xtitle='K!Doutput!N (AB)', ytitle='Completeness', xrange=mock_mag_limits, yrange=[0.,1.2], xstyle=1, ystyle=1, /current)
	p=plot(temp_x[gv], temp_y[gv], xtitle='K!Dinput!N (AB)', ytitle='Completeness', xrange=mock_mag_limits, yrange=[0.,1.2], xstyle=1, ystyle=1, /current)
	p=plot([0,100], 0.9*[1,1], linestyle=2, color="red", /over)
	p=plot([0,100], 0.5*[1,1], linestyle=2, color="red", /over)
	w.save, 'decam_completeness_'+do_filter+'_'+tile_suffix[i]+'.png', BORDER=10, RESOLUTION=300
	w.close

	w=window(dim=[800,600])
	p=plot(mock_in_mag_match[gv_out[0:20000]], (mock_out_mag-mock_in_mag_match)[gv_out[0:20000]], xtitle='K!Dinput!N (AB)', ytitle='K!Doutput!N - K!Dinput!N (AB)', xrange=mock_mag_limits, yrange=[-0.25,0.25], xstyle=1, ystyle=1, symbol="circle", sym_filled=1, sym_size=0.2, sym_transparency=70, linestyle="none", /current)
	p=plot([0,100], 0.*[1,1], linestyle=2, color="red", /overplot)
	w.save, 'decam_completeness_'+do_filter+'_mag_comparison_'+tile_suffix[i]+'.png', BORDER=10, RESOLUTION=600
	w.close

;	mock_in_x=0
;	mock_in_y=0
;	mock_in_mag=0
;	mock_in_mag_match=0
;	mock_out_x=0
;	mock_out_y=0
;	mock_out_mag=0
;	continue

	im_h=headfits(im_file[i])
	im_size=[sxpar(im_h,'NAXIS1'),sxpar(im_h,'NAXIS2')]
	grid_nx=100
	grid_ny=100
	grid_xrange=(findgen(grid_nx+1)/grid_nx*(im_size[0]-1))>0L
	grid_yrange=(findgen(grid_ny+1)/grid_ny*(im_size[1]-1))>0L
	grid_data=fltarr(grid_nx,grid_ny)

  	grid_x=((grid_xrange+shift(grid_xrange,-1))/2.)[0:-2]
  	grid_x=grid_x#make_array(grid_ny, /float, value=1.)
  	grid_y=((grid_yrange+shift(grid_yrange,-1))/2.)[0:-2]
  	grid_y=make_array(grid_nx, /float, value=1.)#grid_y

	xhist_bin=0.4
	for ii=0L, grid_nx-1 do begin
		print, 'Running iteration ii='+strn(ii)
		for jj=0L, grid_ny-1 do begin
			gv_in=where( mock_in_x GE grid_xrange[ii] AND mock_in_x LT grid_xrange[ii+1] AND mock_in_y GE grid_yrange[jj] AND mock_in_y LT grid_yrange[jj+1], n_gv_in)
			gv_out=where( mock_out_x GE grid_xrange[ii] AND mock_out_x LT grid_xrange[ii+1] AND mock_out_y GE grid_yrange[jj] AND mock_out_y LT grid_yrange[jj+1], n_gv_out)

			if n_gv_in GT 1000L AND n_gv_out GT 100L then begin
				yhist_in=histogram(mock_in_mag[gv_in], bin=xhist_bin, locations=xhist_in, min=19., max=21.999, /nan)
				xhist_in += xhist_bin/2.
				yhist_out=histogram(mock_out_mag[gv_out], bin=xhist_bin, locations=xhist_out, min=19., max=21.999, /nan)
				xhist_out += xhist_bin/2.

				temp=yhist_out*1./yhist_in
				gv=where(temp GT 0., n_gv)
				if jj EQ 5 then	begin
					plot, xhist_out+1.82, temp, xtitle='K (AB)', ytitle='mock_out/mock_in', yrange=[0.,1.2], ystyle=1
					print, max(yhist_out*1./yhist_in)
					print, 'mag_50 ', interpol(xhist_out, temp, 0.5)+1.82 
				endif

				mag_completeness=interpol(xhist_out, temp, 0.5) 
				grid_data[ii,jj]=mag_completeness+1.82
			endif else $
				grid_data[ii,jj]=!values.f_nan

		endfor
	endfor
	save, /variables, file='decam_completeness_'+tile_suffix[i]+'.sav'
;	mock_in_x=0
;	mock_in_y=0
;	mock_in_mag=0
;	mock_in_mag_match=0
;	mock_out_x=0
;	mock_out_y=0
;	mock_out_mag=0

	mag_limits=mock_mag_limits;[22.4,24.0]
	mag_step=0.2
	mag_levels=(mag_limits[1]-mag_limits[0])/mag_step+1

	ngrid_data=grid_data
;	print, im_mask_region
;	stop
	gv=where(im_mask_region.tile EQ tile_suffix[i], n_gv)
	if n_gv EQ 1 then begin
		gv_mask=where(grid_x GE im_mask_region[gv].xrange[0] AND grid_x LE im_mask_region[gv].xrange[1] AND grid_y GE im_mask_region[gv].yrange[0] AND grid_y LE im_mask_region[gv].yrange[1], n_gv_mask)
		if n_gv_mask GT 0 then ngrid_data[gv_mask]=!values.f_nan
	endif

	gv=where(ngrid_data GE mag_limits[0]-0.2 AND ngrid_data LE mag_limits[1]+0.2, n_gv, COMPLEMENT=bv, NCOMPLEMENT=n_bv)
	if n_bv GT 0 then ngrid_data[bv]=!VALUES.F_NAN
		ngrid_data=smooth(ngrid_data, 2, /EDGE, /NAN)
;	
		gv=where(finite(ngrid_data))
		ngrid_nx=grid_nx*2
		ngrid_ny=grid_ny*2
		ngrid_x=((findgen(ngrid_nx)+0.5)*(im_size[0]-1.)/ngrid_nx)#make_array(ngrid_ny, /float, value=1.)
		ngrid_y=make_array(ngrid_nx, /float, value=1.)#((findgen(ngrid_ny)+0.5)*(im_size[1]-1.)/ngrid_ny)
;		ngrid_data=tri_surf(ngrid_data[gv], grid_x[gv], grid_y[gv], BOUNDS=[ngrid_x[0,0],ngrid_y[0,0],ngrid_x[-1,-1],ngrid_y[-1,-1]], NX=ngrid_nx, NY=ngrid_ny )
;;	ngrid_data=min_curve_surf(ngrid_data[gv], grid_x[gv], grid_y[gv], XPOUT=ngrid_x, YPOUT=ngrid_y)
;;  ngrid_data=griddata(grid_x[gv], grid_y[gv], grid_data[gv], /MIN_CURV, START=[ngrid_x[0],ngrid_y[0]], DELTA=(im_size-1.)/[ngrid_nx,ngrid_ny], DIM=[ngrid_nx, ngrid_ny])

  		im_h=headfits(im_file[i])
		extast, im_h, im_ast
		xy2ad, [im_size[0]/2.,im_size[0]/2.,0,im_size[0]], [0,im_size[1],im_size[1]/2.,im_size[1]/2.], im_ast, temp_ra_c, temp_dec_c
		temp_ra=round(im_ast.crval[0]*6)/6. + 10./60*(findgen(21)-10)
		temp_dec=round(im_ast.crval[1]*6)/6. + 10./60*(findgen(21)-10)
		ad2xy, temp_ra, replicate(temp_dec_c[0],21), im_ast, temp_x_ra, temp_y_ra
		ad2xy, replicate(temp_ra_c[2],21), temp_dec, im_ast, temp_x_dec, temp_y_dec
		gv_ra=where(temp_x_ra GE 0. AND temp_x_ra LT im_size[0], n_gv_ra)
		gv_dec=where(temp_y_dec GE 0. AND temp_y_dec LT im_size[1], n_gv_dec)

		ngrid_data=(ngrid_data > (mag_limits[0]+0.02) < (mag_limits[1]-0.02))

;		gv=where(im_blank_region.tile EQ tile_suffix[i], n_gv)
;		if n_gv EQ 1 then begin
;			gv_xrange=where(im_blank_region[gv].xrange EQ -1, n_gv_xrange)
;			gv_yrange=where(im_blank_region[gv].yrange EQ -1, n_gv_yrange)
;			if n_gv_xrange GT 0 then im_blank_region[gv].xrange[gv_xrange]=im_size[0]-1
;		if n_gv_yrange GT 0 then im_blank_region[gv].yrange[gv_yrange]=im_size[1]-1

;		gv_blank=where(ngrid_x GE im_blank_region[gv].xrange[0] AND ngrid_x LE im_blank_region[gv].xrange[1] AND ngrid_y GE im_blank_region[gv].yrange[0] AND ngrid_y LE im_blank_region[gv].yrange[1], n_gv_blank)
;		if n_gv_blank GT 0 then ngrid_data[gv_blank]=!values.f_nan
;	endif


	
  loadct, 25, RGB_TABLE=color_table
  color_table=transpose(color_table)
  color_id=round(indgen(mag_levels)*255./(mag_levels-1))
	w=window(dim=[800,600])
	im=contour(ngrid_data, ngrid_x, ngrid_y, /FILL, position=[0.35,0.1,0.95,0.9], title='Tile '+strmid(tile_suffix[i],1,10), c_value=mag_limits[0]+(mag_limits[1]-mag_limits[0])*findgen(mag_levels)/(mag_levels-1),  RGB_TABLE=25, RGB_IND=color_id, axis_style=0, font_size=12, /current, xtitle='X (pixels)', ytitle='Y (pixels)', xrange=[0,im_size[0]-1], yrange=[0,im_size[1]-1])
	radec, temp_ra[gv_ra], temp_dec[gv_dec], temp_ra1, temp_ra2, temp_ra3, temp_dec1, temp_dec2, temp_dec3
	tickname=string(floor(temp_ra[gv_ra]),FORMAT='(I0)')+':'+repstr(string( round((temp_ra[gv_ra]-floor(temp_ra[gv_ra]))*60.),FORMAT='(I2)'),' ','0')
	xaxis=axis('X', location=[0.,0.], TITLE='RA (degrees)', TICKVALUES=reverse(temp_x_ra[gv_ra]), TICKNAME=reverse(tickname), MINOR=4)
	xaxis=axis('X', location=[0.,im_size[1]], TICKVALUES=reverse(temp_x_ra[gv_ra]), TICKNAME=replicate(' ',n_gv_ra), MINOR=4, tickdir=1)
	tickname=string(floor(temp_dec[gv_dec]),FORMAT='(I0)')+':'+repstr(string( round((temp_dec[gv_dec]-floor(temp_dec[gv_dec]))*60.),FORMAT='(I2)'),' ','0')
	yaxis=axis('Y', location=[0.,0.], TITLE='Dec (degress)', TICKVALUES=temp_y_dec[gv_dec], TICKNAME=tickname, MINOR=4)
	xaxis=axis('Y', location=[im_size[0],0.], TICKVALUES=temp_y_dec[gv_dec], TICKNAME=replicate(' ',n_gv_dec), MINOR=4, tickdir=1)
	cbar=colorbar(target=im, orientation=1, position=[0.15,0.1,0.20,0.9], title='50% completeness mag', major=(mag_limits[1]-mag_limits[0])/0.4+1, tickname=string(mag_limits[0]+findgen((mag_limits[1]-mag_limits[0])/0.4+1)*0.4, FORMAT='(F0.1)')) 

	w.save, 'ngvsir_completeness_'+tile_suffix[i]+'.png', BORDER=10, RESOLUTION=300
	w.close

endfor

end
