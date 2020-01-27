pro dem_compare_aia, dem,data,temps,trmatrix,fname

;data is aia observations from aia_read/prep
;dem is the dem
;temps is the temp range the DEM covers
;trmatrix is the temp response matrix for the AIA channels
;fname is the prefix for filenames
;set our data scales
;94
datamin=make_array(6)
datamax=make_array(6)
datamin[0]=alog10(0.1)
datamax[0]=alog10(30)

;131
datamin[1]=alog10(0.7)
datamax[1]=alog10(500)

;171
datamin[2]=sqrt(0.1)
datamax[2]=sqrt(12000)

;193
datamin[3]=alog10(20)
datamax[3]=alog10(2500)

;211
datamin[4]=alog10(7)
datamax[4]=alog10(1500)

;335
datamin[5]=alog10(0.4)
datamax[5]=alog10(80)
wavenum=['94','131','171','193','211','335']

;convert our dem  per temp bin to dn per aia_channel
dem[where(dem le 0)]=0 ;fix negatives just in case
dn_rec=make_array(1024,1024,6)
temp_spacing=temps(1:14)-temps(0:13) ;form converting dem to EM
aia_recon=trmatrix(1:28:2,*) ;5.7 to 7.1 in 0.1 steps
for i=0,1023 do for j=0,1023 do for k=0,5 do dn_rec(i,j,k)=total(aia_recon(*,k)*dem(i,j,*)*temp_spacing)
dn_rec[where(dn_rec le 0, /null)]=0 ; fix zeros in dn_rec (not necessary)
data[where(dn_rec le 0, /null)]=0 ; for comparison we wipe the data where the dem failed



device, /encapsulated, /color, /isolatin1, /inches, set_font='Helvetica',/TT_FONT, $
bits=8, xsize=14, ysize=20,file=fname+'01.eps'

!p.multi=[0,2,3]
for i =0,1 do begin
  aia_lct,r,g,b,wavelnth=wavenum[i],/load
  plot_image,alog10(data(*,*,i)),max=(datamax[i]),min=(datamin[i])
  plot_image,alog10(dn_rec(*,*,i)),max=(datamax[i]),min=(datamin[i])
endfor
  aia_lct,r,g,b,wavelnth=wavenum[2],/load
  plot_image,sqrt(data(*,*,i)),max=(datamax[2]),min=(datamin[2])
  plot_image,sqrt(dn_rec(*,*,i)),max=(datamax[2]),min=(datamin[2])
device,/close
device, /encapsulated, /color, /isolatin1, /inches, set_font='Helvetica',/TT_FONT, $
  bits=8, xsize=14, ysize=20,file=fname+'02.eps'
!p.multi=[0,2,3]  
for i=3,5 do begin
  aia_lct,r,g,b,wavelnth=wavenum[i],/load
  plot_image,alog10(data(*,*,i)),max=(datamax[i]),min=(datamin[i])
  plot_image,alog10(dn_rec(*,*,i)),max=(datamax[i]),min=(datamin[i])
endfor
device, /close


end