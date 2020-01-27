pro hv_dem2jp2,dem,details, aia_ind
; small number for logs
  d_s=1e-20
  
  ;fetch general jp2gen info, deprecate soon
  ;g=hvs_gen()
  
  ;add some tags to the aia header
;  aia_ind = add_tag(aia_ind,info.observatory,'hv_observatory')
;  aia_ind = add_tag(aia_ind,info.instrument,'hv_instrument')
;  aia_ind = add_tag(aia_ind,info.detector,'hv_detector')
;  aia_ind = add_tag(aia_ind,measurement,'hv_measurement')
;  aia_ind = add_tag(aia_ind,0.0,'hv_rotation')
;  aia_ind = add_tag(aia_ind,progname,'hv_source_program')
  hv_rotation=0
  temp1=strcompress(string(details.tbin1,format='(F4.2)'),/remove_all)
  temp2=strcompress(string(details.tbin2,format='(F4.2)'),/remove_all)
  img = bytscl(alog10(dem+d_s),min=alog10(details.dmin),max=alog10(details.dmax))
  
  loc=details.outdir
  
  date =details. yy + '_' +  details.mm + '_' +  details.dd
  time =  details.hh + '_' +  details.mmm + '_' +   '00' + '_' +  '000'
  jp2_filename = date + '__' + time + '__' +  'dem_reginv_T_'+ temp1+'-'+temp2+'.jp2'
  
  dminString = trim(details.dmin)
  dmaxString = trim(details.dmax)
  hv_img_function = 'Two-dimensional image data IMG'
  
  ;header details
  ;line feed
  lf=string(10b)
  
  ;header details from fits
  xh = ''
  ntags = n_tags(aia_ind)
  tagnames = tag_names(aia_ind)
  tagnames = HV_XML_COMPLIANCE(tagnames)
  jcomm = where(tagnames eq 'COMMENT')
  jhist = where(tagnames eq 'HISTORY')
  jhv = where(strupcase(strmid(tagnames[*],0,3)) eq 'HV_')
  jhva = where(strupcase(strmid(tagnames[*],0,4)) eq 'HVA_')
  indf1=where(tagnames eq 'TIME_D$OBS',ni1)
  if ni1 eq 1 then tagnames[indf1]='TIME-OBS'
  indf2=where(tagnames eq 'DATE_D$OBS',ni2)
  if ni2 eq 1 then tagnames[indf2]='DATE-OBS'
  xh='<?xml version="1.0" encoding="UTF-8"?>'+lf
  
  ;
  xh+='<meta>'+lf
  ;
  ; FITS keywords
  ;
  xh+='<fits>'+lf
  for j=0,ntags-1 do begin
    if ( (where(j eq jcomm) eq -1) and $
      (where(j eq jhist) eq -1) and $
      (where(j eq jhv) eq -1)   and $
      (where(j eq jhva) eq -1) )then begin
      ;            xh+='<'+tagnames[j]+' descr="">'+strtrim(string(header.(j)),2)+'</'+tagnames[j]+'>'+lf
      value = HV_XML_COMPLIANCE(strtrim(string(aia_ind.(j)),2))
      xh+='<'+tagnames[j]+'>'+value+'</'+tagnames[j]+'>'+lf
    endif
  endfor
  ;

  xh+='</fits>'+lf
  ;
  ; Explicitly encode the allowed Helioviewer JP2 tags
  ;
  xh+='<helioviewer>'+lf
  ;
  ; Original rotation state
  ;
;  xh+='<HV_ROTATION>'+HV_XML_COMPLIANCE(strtrim(string(hv_rotation),2))+'</HV_ROTATION>'+lf
;  ;
;  ; JP2GEN version
;  ;
;  xh+='<HV_JP2GEN_VERSION>'+HV_XML_COMPLIANCE(trim(g.source.jp2gen_version))+'</HV_JP2GEN_VERSION>'+lf
;  ;
;  ; JP2GEN branch revision
;  ;
;  xh+='<HV_JP2GEN_BRANCH_REVISION>'+HV_XML_COMPLIANCE(trim(g.source.jp2gen_branch_revision))+'</HV_JP2GEN_BRANCH_REVISION>'+lf
;  ;
  ; HVS setup file
  ;
  ;xh+='<HV_HVS_DETAILS_FILENAME>'+HV_XML_COMPLIANCE(trim(info.hvs_details_filename))+'</HV_HVS_DETAILS_FILENAME>'+lf
  ;
  ; HVS setup file version
  ;
  ;xh+='<HV_HVS_DETAILS_FILENAME_VERSION>'+HV_XML_COMPLIANCE(trim(info.hvs_details_filename_version))+'</HV_HVS_DETAILS_FILENAME_VERSION>'+lf
  ;
  ; JP2 comments
  ;
  ;xh+='<HV_COMMENT>'+hv_comment+'</HV_COMMENT>'+lf
  ;
  ; Explicit support from the Helioviewer Project
  ;
  ;xh+='<HV_SUPPORTED>TRUE</HV_SUPPORTED>'+lf
  ;
  ; Clipping values and scaling function
  ;
  xh+='<derived_data>'+lf
  xh+='<DEM_JP2GEN_VERSION>'+details.demjp2_version+'</DEM_JP2GEN_VERSION>'+lf
  xh+='<CONTACT_EMAIL>'+'alasdair.wilson@glasgow.ac.uk'+'</CONTACT_EMAIL>'+lf
  xh+='<GITHUB>'+'https://github.com/alasdairwilson/demjp2'+'</GITHUB>'+lf
  xh+='<DERIVED_QUANTITY>DEM</DERIVED_QUANTITY>'+lf
  xh+='<DERIVATION_METHOD>Regularized Inversion - Hannah and Kontar (2012)</DERIVATION_METHOD>'+lf
  xh+='<DATA_PRODUCT>Temperature</DATA_PRODUCT>'+lf
  xh+='<DEM_JP2GEN_VERSION>0.9</DEM_JP2GEN_VERSION>'+lf
  xh+='<TEMP_RANGE_LOGT>'+HV_XML_COMPLIANCE(temp1)+'-'+HV_XML_COMPLIANCE(temp2)+'</TEMP_RANGE_LOGT>'+lf
  xh+='<TEMP_RANGE_K>'+HV_XML_COMPLIANCE(string(10^details.tbin1,FORMAT='(E8.2)'))+'-'+HV_XML_COMPLIANCE(string(10^details.tbin2,FORMAT='(E8.2)'))+'</TEMP_RANGE_K>'+lf
  xh+='<DEM_UNIT>'+HV_XML_COMPLIANCE(details.units)+'</DEM_UNIT>'+lf
  xh+='<SYS_ERR_PER>'+HV_XML_COMPLIANCE(strcompress(string(details.serr_per),/remove_all))+'</SYS_ERR_PER>'+lf
  xh+='<SAT_LVL>'+HV_XML_COMPLIANCE(strcompress(string(details.sat_lvl),/remove_all))+'</SAT_LVL>'+lf
  xh+='<MIN_SNR>'+HV_XML_COMPLIANCE(strcompress(string(details.min_snr),/remove_all))+'</MIN_SNR>'+lf
  xh+='<N_TEMP_BINS>'+HV_XML_COMPLIANCE(strcompress(string(details.nt),/remove_all))+'</N_TEMP_BINS>'+lf
  xh+='<IMG_SCALE>log 10</IMG_SCALE>'+lf
  xh+='<IMG_DMIN>'+HV_XML_COMPLIANCE(dminString)+'</IMG_DMIN>'+lf
  xh+='<IMG_DMAX>'+HV_XML_COMPLIANCE(dmaxString)+'</IMG_DMAX>'+lf
  xh+='<N_HV_IMGS>'+'image:'+HV_XML_COMPLIANCE(strcompress(string(details.img),/remove_all))+' of '+HV_XML_COMPLIANCE(strcompress(string(details.nimage),/remove_all))+'</N_HV_IMGS>'+lf
  xh+='</derived_data>'+lf
  
  ;xh+='<HV_IMG_FUNCTION>'+hv_img_function+'</HV_IMG_FUNCTION>'+lf
  ;
  ; Close the Helioviewer information
  ;
  xh+='</helioviewer>'+lf
  ;
  ; Enclose all the XML elements in their own container
  ;
  xh+='</meta>'+lf
  
  oJP2 = OBJ_NEW('IDLffJPEG2000',loc + jp2_filename,/WRITE,$
    PROGRESSION = 'RPCL',xml=xh)
  oJP2->SetData,img
  OBJ_DESTROY, oJP2
  print,' '
  print,' created ' + loc + jp2_filename

  end
  