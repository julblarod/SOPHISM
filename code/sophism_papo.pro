PRO sophism_papo

; ==================================================================
; DESCRIPTION
;   For the case of Pupil Apodisation in SOPHISM.
;   Convolution with OTF, and etalon transmission are performed here.   
;
; CALLED FROM
;    sophism
;
; SUBROUTINES
;    sophism_apo2d
;
; MAIN INPUTS
;   otfs          Array of wavelength dependant OTFs (from OTF module)
;   fffav         Etalon transmission curve for each selected
;                    wavelength, averaged over the pupil (from etalon module)
;   int_range     Wavelength range for integration around the selected
;                    position [Angstroms]
;   wavesind      Wavelength positions selected for simulation
;
; OUTPUT
;   Data convolved with OTF, spectral
;
; LIMITATIONS
;   Data must have gone through modulation (Polarization
;   module). Filtergraph and Optics modules also must have been run.
;
; VERSION HISTORY
;   20-Apr-11 J. Hirzberger, MPS
;   J. Blanco. May 2012. v1_0. Adapted to be included in
;      sophism. Integration range selectable.
;   J. Blanco. Sep 2012. v2_0. Changes of etalon transmission
;   J. Blanco. Nov 2012. v2_1. Neighbouring wavelengths weighting
;      included when no sampled data at selected wavelengths
;   J. Blanco. Dec 2012. v2_2. Prints for talk mode
;   J. Blanco. Jan 2013. v2_3. Corrected bug for number of files (to
;      not take into account already-discarded frames)
;   J. Blanco. Feb 2013. v2_4. Corrected format for screen print of
;      sample numbers. 
;   J. Blanco. Jul 2014. v2_5. Corrected ntim according to new cycrep scheme
;   J. Blanco. Jun 2016. v2_6. Included start at some given sample
;      option (nstart)
;   J. Blanco. Apr 2017. v2_7. Adapted for no time series case
;
; ==================================================================

;--------------------------------------------------------------------- 
; Initialize common blocks
;--------------------------------------------------------------------- 

progind=5
print,'sophism_papo'

restore,'../settings/settings_sophism.sav'
;load outputs from etalon module
restore,info.saves+info.files(3)+'.sav'
;load outputs from OTF module
restore,info.saves+info.files(4)+'.sav'

;ntim=info.modnbuff*info.nstate*info.cycrep*info.etalpos
ntim=info.modnbuff*info.nstate*info.etalpos
if (info.dataser eq 0) then ntim=1

;--------------------------------------------------------------------- 
; Convolution with monochromatic OTFs
;--------------------------------------------------------------------- 

; scan positions in 'pixel' units (not wavelength)
wsteps=wavesind

npos = info.etalpos
nwave=n_elements(index)

;integration interval
int_range=fix(round(info.intrange/info.wavesamp))

;limit the int_range in the low part
int_range=wsteps < int_range
;careful!!! in this way, some integration ranges may be smaller than
;others. It should be 0 most of the 'other' parts, but still....

;and limit in the upper part
wrange=where(wsteps+int_range gt info.wavepoint-1)
if (min(wrange) ne -1) then int_range[wrange]=info.wavepoint-1-wsteps[wrange]
;again, careful!!! integration ranges may be different in diff. positions

; normalization factor (with the averaged transmission for the pupil)
norm=fltarr(nwave)
for pos=0,nwave-1 do norm[pos]=total(fffav[wsteps[pos]-int_range[pos]:wsteps[pos]+int_range[pos],pos],1)

;??????????
if (info.startmed eq 1 AND min(where(info.routines eq 1)) eq progind) then nstart=info.nstart else nstart=0
for nf=nstart,ntim-1 do begin
;for nf=0,ntim-1 do begin
;???????????
; define array of degraded Stokes profiles (only 1 in Stokes dim because it
; has already been modulated)
   ima=readfits(info.saves+info.files(2)+'_'+strtrim(nf,2)+'.fits*',head,/compr,/sil)
   szp=size(ima)
;   lcc = fltarr(nwave,1,info.sz,info.sz) 
   lcc = fltarr(nwave,szp(2),info.sz,info.sz) 

;loop over modulation states, in case of no data series
   for ppp=0,szp(2)-1 do begin
; loop over scan positions
      FOR j=0,nwave-1 DO BEGIN 
; integrated image
         imc = fltarr(info.sz,info.sz) 
; loop over integration range
         FOR i=-int_range[j],int_range[j] DO BEGIN 
            if (info.talk eq 1) then print,format='(52(%"\b"),"Integrating sample ",I5,", wavelength ",I3,", step ",I5,$)',nf+1,j+1,i+1
;convolution of image with OTF for the wavelength position and
;integration step
;         im = reform(ima[i+wsteps[j],0,*,*])
            im = reform(ima[i+wsteps[j],ppp,*,*])
; apodization
            if (info.sscene_apod gt 0) then sophism_apo2d,im,info.sscene_apod
;convolution
            fim = fft(im,-1)         
            otf = reform(otfs[*,*,wsteps[j]+i,j])
            fim = fim*otf
            imr = fft(fim,1)

;integral and etalon transmission
            imc=imc+fffav[wsteps[j]+i,j]*sqrt(float(imr*conj(imr)))

         ENDFOR ;int_range
; normalization
;         lcc[j,0,*,*] = imc/norm[j]
         lcc[j,ppp,*,*] = imc/norm[j]

      ENDFOR ;scan pos
   endfor ;ppp

;now 'interpolate' to real selected wavelength position (if it was not
;sampled in input data)
   lccpre=lcc
;   lcc=fltarr(npos,1,info.sz,info.sz) ;lcc*0.
   lcc=fltarr(npos,szp(2),info.sz,info.sz) ;lcc*0.
   mind=max(index)
   for ii=0,mind do begin
      indy=where(index eq ii,nindy)
;      for bb=0,nindy-1 do lcc[ii,0,*,*]=lcc[ii,0,*,*]+(wg[indy(bb)]*lccpre[indy(bb),*,*,*])
      for bb=0,nindy-1 do lcc[ii,*,*,*]=lcc[ii,*,*,*]+(wg[indy(bb)]*lccpre[indy(bb),*,*,*])
   endfor 

   sxaddpar,head,"HISTORY",'Pupil Apodization Mode'
   sxaddpar,head,"HISTORY",'Etalon Module'
   sxaddpar,head,"FRATIO",info.frat
   for ii=0,info.etalpos-1 do sxaddpar,head,"WVLTRAN"+strtrim(ii,2),format='(F8.3)',llo(ii),'Transmitted wavelength (A)'
   for ii=0,info.etalpos-1 do sxaddpar,head,"FWHMREAL"+strtrim(ii,2),format='(F6.3)',fwhm_real(ii),'FWHM (A)'
   sxaddpar,head,"HISTORY",'Optics Module'
   sxaddpar,head,"TYPE",info.otftype,'Diffraction or Zemax'
   if (info.sscene_apod gt 0) then sxaddpar,head,"APODIZAT",format='(F5.2)',info.sscene_apod,'%'
   zerhead=where(info.zerco ne 0)
   if (max(zerhead) ne -1) then for cc=0,n_elements(zerhead) do sxaddpar,head,"ZERCO"+strtrim(zerhead+1,2),info.zerco(cc),'Zernike coefficients included'
   
   sxaddpar,head,"HISTORY",'Pupil Apod Module'

   if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',lcc,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',lcc,head


endfor ;ntim
print,''

;show example of data convolved
if (info.talk eq 1) then begin
   window,progind,title='Convolved image example',/free
   tvframe,lcc[0,0,*,*],/asp
endif


;#######################################################################
;dual-beam part
if (info.dualbeam eq 1) then begin
   print,''
   print,'Starting papo dual-beam part'

   for nf=0,ntim-1 do begin
; define array of degraded Stokes profiles (only 1 in Stokes dim because it
; has already been modulated)
      ima=readfits(info.saves+info.files(2)+'_2_'+strtrim(nf,2)+'.fits*',head,/compr,/sil)
      lcc = fltarr(nwave,szp(2),info.sz,info.sz) 

      for ppp=0,szp(2)-1 do begin
; loop over scan positions
         FOR j=0,nwave-1 DO BEGIN
; integrated image
            imc = fltarr(info.sz,info.sz)
; loop over integration range
            FOR i=-int_range[j],int_range[j] DO BEGIN
               if (info.talk eq 1) then print,format='(57(%"\b"),"DUAL Integrating sample ",I5,", wavelength ",I3,", step ",I5,$)',nf+1,j+1,i+1
;convolution of image with OTF for the wavelength position and
;integration step
;               im = reform(ima[i+wsteps[j],0,*,*])
               im = reform(ima[i+wsteps[j],ppp,*,*])
; apodization 
               if (info.sscene_apod gt 0) then sophism_apo2d,im,info.sscene_apod
               fim = fft(im,-1)          
               otf = reform(otfs[*,*,wsteps[j]+i,j])
               fim = fim*otf 
               imr = fft(fim,1)
;integral and etalon transmission
               imc=imc+fffav[wsteps[j]+i,j]*sqrt(float(imr*conj(imr)))

            ENDFOR ;int_range
 ; normalization
;         lcc[j,0,*,*] = imc/norm[j]
            lcc[j,ppp,*,*] = imc/norm[j]

         ENDFOR ;scan pos
      endfor 

;now 'interpolate' to real selected position, if it was not sampled in
;input data
;no need to do it for llo and fwhm, it was done on single-beam part
      lccpre=lcc
      lcc=fltarr(npos,1,info.sz,info.sz) ;lcc*0.
      mind=max(index)
      for ii=0,mind do begin
         indy=where(index eq ii,nindy)
         for bb=0,nindy-1 do lcc[ii,0,*,*]=lcc[ii,0,*,*]+(wg[indy(bb)]*lccpre[indy(bb),*,*,*])
      endfor


      sxaddpar,head,"HISTORY",'Pupil Apodization Mode'
      sxaddpar,head,"HISTORY",'Etalon Module'
      sxaddpar,head,"FRATIO",info.frat
      for ii=0,info.etalpos-1 do sxaddpar,head,"WVLTRAN"+strtrim(ii,2),format='(F8.3)',llo(ii),'Transmitted wavelength (A)'
      for ii=0,info.etalpos-1 do sxaddpar,head,"FWHMREAL"+strtrim(ii,2),format='(F6.3)',fwhm_real(ii),'FWHM (A)'
      sxaddpar,head,"HISTORY",'Optics Module'
      sxaddpar,head,"TYPE",info.otftype,'Diffraction or Zemax'
      if (info.sscene_apod gt 0) then sxaddpar,head,"APODIZAT",format='(F5.2)',info.sscene_apod,'%'
      zerhead=where(info.zerco ne 0)
      if (max(zerhead) ne -1) then for cc=0,n_elements(zerhead) do sxaddpar,head,"ZERCO"+strtrim(zerhead+1,2),info.zerco(cc),'Zernike coefficients included'
      sxaddpar,head,"HISTORY",'Pupil Apod Module'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',lcc,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',lcc,head

   endfor ;ntim
   print,''

endif ;dualbeam

;########################################################################




END
