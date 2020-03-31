pro batch_dem_jp2_v1_1, t_start,cadence,n_obs,fits_base_dir,jp2_base_dir,get_fits=get_fits,min_snr=min_snr,sat_lvl=sat_lvl,serr_per=serr_per,gauss_stdev=gauss_stdev,gauss_mean=gauss_mean,fe_min=fe_min
  ;creates a batch of dem from synoptic fits data
  ;t_start is the date for the first observation
  ;cadence is time between observations (in seconds) this MUST be a multiple of 120.
  ;n_obs is number of observations
  ;get_fits, set this flag to download data, leave it if data is predownloaded
  ; Only calculate the DEM for pixels which has this minimum SNR
  if (n_elements(min_snr) lt 1) then min_snr=3.0
  ; Ignore data with values above this level
  ; This is the default AIA saturation of >15,000 DN/px
  if (n_elements(sat_lvl) lt 1) then sat_lvl=1.5e4

  ;Systematic uncertainty (in % terms) to add to the data
  if (n_elements(serr_per) lt 1) then serr_per=0.0
  
  if (n_elements(fe_min) lt 1) then fe_min=0.0

  ;set up base directories for fits & output jpegs


  t_start_s = anytim(t_start,/utime)
  ;work out the start time
  ;and we loop over all sets of observations

  for i=0,n_obs-1 do begin
    skip=0
    ;fetch date & time
    obs_time=anytim(t_start_s)+double(i)*double(cadence)
    ;stringify
    ystr=strmid(anytim(obs_time,/ccsds),0,4)
    mstr=strmid(anytim(obs_time,/ccsds),5,2)
    dstr=strmid(anytim(obs_time,/ccsds),8,2)
    hour=strmid(anytim(obs_time,/ccsds),11,2)
    mins=strmid(anytim(obs_time,/ccsds),14,2)
    fits_dir=fits_base_dir+ystr+'/'+mstr+'/'+dstr+'/'
    jp2_dir=jp2_base_dir+ystr+'/'+mstr+'/'+dstr+'/'
    ;check the fits directory exists
    ftest=file_test(fits_dir,/directory)
    if (ftest ne 1) then file_mkdir, fits_dir
    ;6 coronal channels
    wavenum=['94','131','171','193','211','335']
    ;synoptic filenames
    ffs='AIA'+ystr+mstr+dstr+'_'+hour+mins+'_0'+['094','131','171','193','211','335']+'.fits'
    if (n_elements(get_fits) ge 1) then begin
      ;here we downloads synoptic data
      ;synoptic folder url
      durl='http://jsoc.stanford.edu/data/aia/synoptic/'+ystr+'/'+mstr+'/'+dstr+'/H'+hour+'00/'
      ;total url
      furls=durl+ffs
      ;wget
      wget1=ssw_wget_mirror(furls,fits_dir,/spawn)
    endif
    ;check for missing data, go to next iteration if not all files there.
    ftest=file_test(fits_dir+ffs)
    for ll=0,5 do begin
      if ((ftest[ll] eq 0)) then begin
        skip=1
      endif
    endfor
    if (skip eq 1) then begin
      print, 'data missing for:'
      print, anytim(obs_time, /ccsds)
      print, 'skipping'
      continue
    endif

    ;read the fits

    read_sdo,fits_dir+ffs,ind,data,/uncomp_delete

    ; Just double check in the correct order of 94, 131, 171, 193, 211, 335
    wv_srt=sort(ind.wavelnth)
    ind=ind[wv_srt]
    data=data[*,*,wv_srt]
    ;calculate the warm component of AIA 94 channel.
    a94_fe18=make_array(1024,1024)
    a94_fe18[*,*]=data[*,*,0]-data[*,*,4]/120d0-data[*,*,2]/450d0
    ;set any pixels where the warm compoent is
    a94_fe18[where(a94_fe18 lt fe_min)]=0

    nx=n_elements(data[*,0,0])
    ny=n_elements(data[0,*,0])
    nf=n_elements(data[0,0,*])

    xc=ind[0].xcen
    yc=ind[0].ycen
    dx=ind[0].cdelt1
    dy=ind[0].cdelt2

    ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ; Filter out the "bad" pixels
    ; Only need to set the actual bad ones to 0 (not all filters at that location)
    ; as the DEM code will only perform the calc if all filters are non-zero at that pixel location

    ; If the value is above sat_lvl get rid of it
    ; assuming bad from saturation
    id=where(data ge sat_lvl,nid)
    if (nid gt 1) then data[id]=0.0

    ; Anything less than 0 just set to 0
    id=where(data le 0,nid)
    if (nid gt 1) then data[id]=0.0

    ; Set the edge few pixels to 0 to minimise problems - from prepd data that has been shifted (?)
    edg0=10
    data[0:edg0-1,*,*]=0.0
    data[*,0:edg0-1,*]=0.0
    data[nx-edg0:nx-1,*,*]=0.0
    data[*,nx-edg0:nx-1,*]=0.0



    ;The synoptic data is in DN/px

    ; Proper way but a bit slow?
    ;  for i=0,nf-1 do edata[*,*,i]=aia_bp_estimate_error(reform(data[*,*,i]),replicate(wavenum[i],nx,ny),n_sample=16)
    ; So instead quick approx based on aia_bp_estimate_error.pro
    npix=4096.^2/(nx*ny)
    edata=fltarr(nx,ny,nf)
    gains=[18.3,17.6,17.7,18.3,18.3,17.6]
    dn2ph=gains*[94,131,171,193,211,335]/3397.0
    rdnse=1.15*sqrt(npix)/npix
    drknse=0.17
    qntnse=0.288819*sqrt(npix)/npix
      ; error in DN/px
    data(data le 0)=0
    for j=0, nf-1 do begin
      etemp=sqrt(rdnse^2.+drknse^2.+qntnse^2.+(dn2ph[j]*abs(data[*,*,j]))/(npix*dn2ph[j]^2))
      esys=serr_per*data[*,*,j]/100.
      edata[*,*,j]=sqrt(etemp^2. + esys^2.)
    endfor

    ; Get rid of data with too large an uncertaintity
    id=where(data[*,*,0:5]/edata[*,*,0:5] le min_snr,nid);
    if (nid gt 1) then data[id]=0.0

    ; For the DEM code need the data in DN/px/s;
    durs=ind.exptime
    for k=0, nf-1 do data[*,*,k]=data[*,*,k]/durs[k]
    for k=0, nf-1 do edata[*,*,k]=edata[*,*,k]/durs[k];

    ;  What temperature binning do you want for the DEM?
    ; temps variable are the bin edges
    ; These are the bin edges
    ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
    ;  temps=[0.5,1,2,4,6,8,11,14,19]*1d6
    ;  logtemps=alog10(temps)

    ;  ;or more bins and then rebin at the end?
    ;  0.1 binsizing
    ;logtemps=5.7+findgen(17)*0.1
    ;  0.05 binsizing
    logtemps=5.7+findgen(17)*0.1

    temps7=10d^logtemps
    temps=temps7(0:14)
    ; This is is the temperature bin mid-points
    mlogt=get_edges(logtemps,/mean)
    nt=n_elements(mlogt)

    ; Given that the responses are time dependent, use one on a monthly basis
    ; Calculate and save if not already done so, then load in
    respfile='resp/aia_resp_'+mstr+ystr+'.dat'
    ; Need to make the response functions in advance and save them?
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm,timedepend_date=ystr+'/'+mstr+'/'+dstr)

    ; Only want the coronal ones without 304A
    idc=[0,1,2,3,4,6]

    tr_logt=tresp.logte

    ; Don't need the response outside of the T range we want for the DEM
    gdt=where(tr_logt ge min(logtemps) and tr_logt le max(logtemps),ngd)
    tr_logt=tr_logt[gdt]
    TRmatrix=tresp.all[*,idc]
    TRmatrix=TRmatrix[gdt,*]
    TRfe18= (TRmatrix[*,0]-TRmatrix[*,4]/120.-TRmatrix[*,2]/450.) >0.
    ; remove low peak
    TRfe18[where(tr_logt le 6.4)]=0.
    
    TR7=dblarr(n_elements(TRmatrix[*,0]),7)
    TR7[*,0:5]=TRmatrix
    TR7[*,6]=TRfe18
    
    data7=dblarr(nx,ny,7)
    data7[*,*,0:5]=data 
    data7[*,*,6]=a94_fe18
    
    edata7=dblarr(nx,ny,7)
    edata7[*,*,0:5]=edata
    ;10% error?
    edata7[*,*,6]=a94_fe18*0.1+1
    ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ; Use the initial weighting to the DEM calc or not (default yes)
    dem_norm0=dblarr(nx,ny,nt-2)
    dem_norm07=dblarr(nx,ny,nt)
    dem_norm_temp=[0.01,0.0355,0.1614,0.4915,1.0024,1.3689,1.2518,0.7664,0.3142,0.0863,0.0159,0.0020,0.0002,0.0001]
    dem_norm_temp7=[0.01,0.0355,0.1614,0.4915,1.0024,1.3689,1.2518,0.7664,0.3142,0.0863,0.0249,0.0040,0.001,0.0005,0.0002,0.0001]
    dem_norm_temp7=[0.0303,0.0836,0.1888,0.3487,0.5269,0.6511,0.6580,0.5438,0.3676,0.2032,0.0919,0.0340,0.0103,0.0025,0.0005,0.0001]
    dem_norm_temp7=[0.0103,0.0336,0.0888,0.3487,0.5269,0.6511,0.6580,0.5438,0.3676,0.2032,0.0919,0.0340,0.0203,0.0065,0.0015,0.0005]
    for xx=0,nx-1 do begin
      for yy=0,ny-1 do begin
        dem_norm0[xx,yy,*]=dem_norm_temp
        dem_norm07[xx,yy,*]=dem_norm_temp7
      endfor
    endfor
    ; Do DEM calculation


    dn2dem_pos_nb, data7,edata7,TR7,tr_logt,temps7,dem7,edem7,elogt,chisq,dn_reg,/timed,dem_norm0=dem_norm07
    dn2dem_pos_nb, data,edata,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed,dem_norm0=dem_norm0
    ;any negative data gets wiped
    ;debug purposes we print
    dem_neg=n_elements(where(dem lt 0))
    print,'number of < 0 elements =',dem_neg-1,'  fraction of dem computed =',double(n_elements(dem gt 0))/16./1024./1024.
    for xx=0,nx-1 do begin
      for yy =0,ny-1 do begin
        for kk=0,nt-1 do begin
          if (dem7[xx,yy,kk] le 0) then begin
            dem7[xx,yy,*]=0
            continue
          endif
        endfor
        for kk=0,nt-3 do begin
          if (dem[xx,yy,kk] le 0) then begin
            dem[xx,yy,*]=0
            continue
          endif
        endfor
      endfor
    endfor
dem_combi=dblarr(nx,ny,nt)
dem_combi[*,*,0:nt-3]=dem
mask=array_indices(a94_fe18,where(a94_fe18 ge 3))
sze=size(mask)
for ii=0,sze[2]-1 do dem_combi[mask[0,ii],mask[1,ii],*]=dem7[mask[0,ii],mask[1,ii],*]
stop
    ;rebin the dem
    dem_rebin=rebin(dem,[1024,1024,nt/2])
    ;sanity check there are no negative elements
    print,'number of < 0 elements in final DEM bins =',n_elements(where(dem_rebin lt 0))-1


    ;check jp2 directory
    ftest=file_test(jp2_dir,/directory)
    if (ftest ne 1) then file_mkdir, jp2_dir
    ;loop over all temp bins
    nimage=7
    for l=0,nimage-1 do begin

      ;set up structure containing relevant parameters for input to the header
      params={tbin1:logtemps(l*2),$
        tbin2:logtemps((l+1)*2),$
        jp2_dir:jp2_dir,$
        hh:hour,$
        dd:dstr,$
        yy:ystr,$
        mm:mstr,$
        mmm:mins,$
        nt:nt,$
        gauss_std:gauss_stdev,$
        sat_lvl:sat_lvl,$
        min_snr:min_snr,$
        serr_per:serr_per,$
        img:l+1,$
        nimage:nimage}

      ;call dem_details for generating structure containing header information
      details=dem_details(params)

      ;make jp2000
      hv_dem2jp2,dem_rebin[*,*,l],details,ind(0)
    endfor


  endfor

end
