function dem_details, dem_param
  ;sorts some infromation regarding the DEM into a structure for inclusion in the header file for jp2k for helioviewer
  ;most of these are extracted from either the corresponding AIA header or from a structure containing parameters used in the dem setup.
  
    
    



    
;date
;  tobs=hv_parse_ccsds(index.date_obs)
  det={yy:dem_param.yy,$
  mm:dem_param.mm,$
  dd:dem_param.dd,$
  hh:dem_param.hh,$
  mmm:dem_param.mmm,$
  dmin:1e19,$ ; min and max for bytescale   
  dmax:1e23,$
  outdir:dem_param.jp2_dir,$
  tbin1:dem_param.tbin1,$
  tbin2:dem_param.tbin2,$
  nt:dem_param.nt,$
  sat_lvl:dem_param.sat_lvl,$
  min_snr:dem_param.min_snr,$
  serr_per:dem_param.serr_per,$
  img:dem_param.img,$
  nimage:dem_param.nimage,$
  units:'cm−5 K−1',$
  demjp2_version:'1.0'};files
  
 
 return, det
end