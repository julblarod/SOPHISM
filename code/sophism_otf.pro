pro sophism_otf

; ==============================================================================
;
; DESCRIPTION
;   Define OTF from diffraction or loading wavefronts, adding
;   aberrations if selected.
; 
; CALLED FROM 
;    sophism
;
; SUBROUTINES
;    radius_aper2
;    check_rpupil
;    zernike2
;
; MAIN INPUTS
;    zerco          Coefficients for the zernike polynomials
;                      composition of a wavefront
;    zercofile      ASCII file with zernike coefficients to generate a
;                      wavefront
;    frontrms       Normalization of wavefront
;    fff            Transmission curve of the etalon, for amplitudes
;                      of pupils in pupil apodization case
;    phase_et       Wave phases of the etalon, for pupil apodization case
;
; OUTPUT
;    Optical Transfer Functions for each wavelength in case of pupil
;       apodisation. In normal case, the OTF and the convolved images
;
; VERSION HISTORY 
;   Alex Feller   v0.1    2010-09-21 
;   J. Blanco. 2011. v1_0. Fit into sophism
;   J. Blanco. 2012. v2_0. Included zernikes treatment for
;      aberrations. Loading option for zernikes collection. 
;   J. Blanco. April 2012. v3_0. Pupil apodization part included
;   J. Blanco. Sept 2012. v3_1. Wavelength, pupil apodization and aliasing
;      corrections
;   J. Blanco. Dec 2012. v3_2. Add printings for talk mode
;   J. Blanco. Jan 2013. v3_3. Corrected bug for number of files (to
;      not take into account already-discarded frames)
;   J. Blanco. May 2013. v3_4. Added the wavefront to the save file
;   J. Blanco. Oct 2013. v3_5. Modified read_ascii for zercofile to
;      consider comments when starting with ';'
;   J. Blanco. Jul 2014. v3_6. Corrected ntim according to new cycrep scheme
;   J. Blanco. Apr 2017. v3_61. Added ntim=1 for the case of no data series
;
; ==============================================================================

;module index for restoring previous data in the module chain
progind=4

restore,'../settings/settings_sophism.sav'

;ntim=info.modnbuff*info.nstate*info.cycrep*info.etalpos
ntim=info.modnbuff*info.nstate*info.etalpos

print,''
print,'sophism_otf'
print,''

; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------

;calculate radius of pupil

;rpup (cutfreq) is function of lambda. Not big changes at this range but...
;have to load the wavelengths from sophism_linbo?
;for now, take only the line center, shifted by doppler
llo_f=info.wl
dop=info.scvel/299792.5d0
ldop=llo_f*dop
llo_f=llo_f+ldop

cutfreq=(info.telap/(llo_f*1.e-10))*(1./(!radeg*3600.))
deltanu=1./(info.sz*info.sssamp)
rpup=cutfreq/(2.*deltanu)
;now check for the case of aliasing
tam=check_rpupil(info.sz/2,rpup)
;borders to cut the extended pupil into the 'real' one
low=(2*tam-info.sz)/2
high=low+info.sz-1
;define variable to read zemax wavefront into, or just leave empty for
;later
front=fltarr(tam,tam)

; ------------------------------------------------------------------------------
; ZEMAX OTF
; ------------------------------------------------------------------------------

;start with Zemax wavefront. Not implement yet
if info.otftype eq 'zemax' then begin

   print, 'ZEMAX'
   print,'not ready yet'
   return
;careful! probably not yet correctly implemented.
;careful! Fixed wavefront filename. Change
;   header=strarr(16)
;   front=fltarr(64,64)
;;front_tot=fltarr(64,64,9)
;   openr,5,'../PHI_wavefront_0_0.txt'
;   readf,5,header                           
;   readf,5,front
;;front_tot(*,*,ff)=front                               
;   close,5
; stop
;??????
;  if (info.frontrms ne 1.) then begin
;      frontmed=total(front*zernis(*,*,0))/total(zernis(*,*,0))
;      front=(front-frontmed)*(zernis(*,*,0) ne 0)
;      front=double(front)
;      frontrms=sqrt((moment(front))(1))
;      front=front/frontrms*(info.frontrms*2*!pi)
;      front=float(front)
;   endif 

endif


; ------------------------------------------------------------------------------
; PRE-CALCULATED OTF
; ------------------------------------------------------------------------------
;maybe useful in the future...
;if info.otfspar.type eq 'load' then begin
;if (info.talk eq 1) then print,'Pre-calculated'
;restore,'polychromatic_psf_mtf.sav'
;endif


;+++++++++++++++++++++++++++++++++++++++++++++++++
; Aberrations
;+++++++++++++++++++++++++++++++++++++++++++++++++

;Calculate Zernikes and contributions

;Define the zernike polinomials
   zernis=zernike2(tam/2.,rpup,45)

   phi=zernis(*,*,0)*0.

;Add zernike coefficients from ascii file (e.g. entrance window
;aberrations)
   if (info.zercoload eq 1) then begin
      zercofile=read_ascii(info.zercofile,comment=';')
      zerconums=reform(zercofile.field1(1,*))
;remove defocus, since it will be compensated... or not?
;      zerconums(4)=0.
;convert from waves to radians
      zerconums=zerconums*2.*!pi
;add contributions of the different zernikes
      for cc=0,n_elements(zerconums)-1 do phi=phi+zerconums(cc)*zernis(*,*,cc)
   endif

;Add all the different aberrations from input
   for cc=0,n_elements(info.zerco)-1 do phi=phi+info.zerco(cc)*zernis(*,*,cc+2)
;Add phases from Zemax wavefront (first resize to actual pupil size)
   front=congrid(front,tam,tam,/inter,cub=-0.5,/center)
   phi=front+phi
   
;wavefront normalization
   if (info.frontnorm eq 1 AND max(abs(phi)) gt 0) then begin
      phim=total(phi*zernis(*,*,0))/total(zernis(*,*,0))
      phi=phi-phim
      phirms=sqrt(total(phi^2.*zernis(*,*,0))/total(zernis(*,*,0)))
      phi=(phi/phirms)*(2.*!pi*info.frontrms)
      if (info.talk eq 1) then print,format='("Wavefront rms [lambda]  ",F5.3)',sqrt(total(phi^2.*zernis(*,*,0))/total(zernis(*,*,0)))/2./!pi
   endif 
   if (info.talk eq 1) then print,format='("Wavefront rms [lambda]  ",F5.3)',sqrt(total(phi^2.*zernis(*,*,0))/total(zernis(*,*,0)))/2./!pi

;estudiar mejor porque habria que hacer opcion para simultaneo y para seguido
;Add PD defocus
;   phi=phi+info.pddefoc*zernis(*,*,2+2)

;stop
;*************************************************************************
; this first part below, without pupil apodization

   if (info.routines(5) ne 1) then begin

      otfs=fltarr(tam,tam)
;calculate the OTF from the above wavefront and a plain pupil
      hh=zernis(*,*,0)*complex(cos(phi),sin(phi))
      otf=otfunc(hh,norm)
;cut down to 'real' size the OTF (aliasing case)
      otf=otf[low:high,low:high]
      otf=shift(otf,-info.sz/2.,-info.sz/2.)
      otfs=otf

      if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
      if (progma eq -1) then progma=0
      if (info.talk eq 1) then print,'Applying OTF'
      head=headfits(info.saves+info.files(progma)+'_0.fits*')
      headmod=sxpar(head,'HISTORY')
      etal=where(headmod eq 'Etalon Module')
;in case the etalon module hasn't been run, the ntim may have
;to include the observation sets and so, so use the ntim
;calculated when starting simulation
      if (etal eq -1) then ntim=info.ntim-fix(total(info.ndeaths))
      if (info.dataser eq 0) then ntim=1

prg0=0.
prg=prg0   
      for nf=0,ntim-1 do begin
         ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
         sizy=size(ima)
         nlam=sizy(1)
         npol=sizy(2)

         for ii = 0,npol-1 do begin
            for jj=0,nlam-1 do begin 
               a = reform(ima[jj,ii,*,*])
;apply apodisation if enabled
               if (info.sscene_apod gt 0) then sophism_apo2d,a,info.sscene_apod 
;convolve data with OTF
               ima[jj,ii,*,*] = fft(fft(a, -1) * otfs, 1)
;calculate progress of process to show on screen
               prg = round(100.*float(jj+1)*float(ii+1)*float(nf+1)/(npol*nlam*ntim))
               if (info.talk eq 1 AND prg gt prg0) then begin
                  print,format='(16(%"\b"),"  Progress ",I3," %",$)', prg
                  prg0=prg
               endif 
            endfor ;jj
         endfor ;ii

         sxaddpar,head,"HISTORY",'Optic Module'
         sxaddpar,head,"TYPE",info.otftype,'Diffraction or Zemax'
         if (info.zercoload eq 1) then sxaddpar,head,"WFE_FILE",info.zercofile
         if (info.sscene_apod gt 0) then sxaddpar,head,"APODIZAT",format='(F5.2)',info.sscene_apod,'%'
         zerhead=where(info.zerco ne 0)
         if (max(zerhead) ne -1) then for cc=0,n_elements(zerhead)-1 do sxaddpar,head,"ZERCO"+strtrim(zerhead(cc)+1,2),info.zerco(zerhead(cc)),'Zernike coefficients included'

         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',ima,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',ima,head
      endfor ;nf
      print,''

;#####################################################################
; the dual-beam part

      if (info.dualbeam eq 1) then begin
         print,'Starting OTF dual-beam part'
         if (info.talk eq 1) then print,'Applying OTF dual-beam'

prg0=0.
prg=prg0
         for nf=0,ntim-1 do begin 
            ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
            sizy=size(ima)
            nlam=sizy(1)
            npol=sizy(2)

            for ii = 0,npol-1 do begin
               for jj=0,nlam-1 do begin 
                  a = reform(ima[jj,ii,*,*])
                  if (info.sscene_apod gt 0) then sophism_apo2d,a,info.sscene_apod
                  ima[jj,ii,*,*] = fft(fft(a, -1) * otfs, 1)

                  prg = round(100.*float(jj+1)*float(ii+1)*float(nf+1)/(npol*nlam*ntim))
                     if (info.talk eq 1 AND prg gt prg0) then begin
                        print,format='(16(%"\b"),"  Progress ",I3,"%",$)', prg
                        prg0=prg
                     endif 
                  endfor ;jj
            endfor ;ii

            sxaddpar,head,"HISTORY",'Optic Module'
            sxaddpar,head,"TYPE",info.otftype,'Diffraction or Zemax'
            if (info.zercoload eq 1) then sxaddpar,head,"WFE_FILE",info.zercofile
            if (info.sscene_apod gt 0) then sxaddpar,head,"APODIZAT",format='(F5.2)',info.sscene_apod,'%'
            zerhead=where(info.zerco ne 0)
            if (max(zerhead) ne -1) then for cc=0,n_elements(zerhead) do sxaddpar,head,"ZERCO"+strtrim(zerhead+1,2),info.zerco(cc),'Zernike coefficients included'

            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',ima,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',ima,head
         endfor ;nf
         print,''

      endif ;dualbeam
;#####################################################################


   endif else begin

;**************************************************************************
;this process below, for the case of Pupil Apodization

;load etalon information
      restore,info.saves+info.files(3)+'.sav'
;define final OTFs with proper dimensions
      otfs=fltarr(info.sz,info.sz,info.wavepoint,n_elements(index))

;phase shift across pupil (form etalon module). Create an empty array
;of appropriate size in case it is disabled
      phase_et=fltarr(tam,tam,info.wavepoint,n_elements(index))
      if (info.phases eq 1) then begin
         restore,info.saves+info.files(progind-1)+'_phase.sav'
      endif
    
;calculate the OTFs for every wavelength
      for lamb=0,n_elements(index)-1 do begin
         for step=0,info.wavepoint-1 do begin
;preparing the amplitude of the pupils
            pupigran=sqrt(reform(fff[*,*,step,lamb]))

            hh=pupigran*complex(cos(phi+phase_et[*,*,step,lamb]),sin(phi+phase_et[*,*,step,lamb]))
            otf=otfunc(hh,norm)
;cutting borders of extra size in case of aliasing
            otf=otf[low:high,low:high]
            otf=shift(otf,-info.sz/2.,-info.sz/2.)

            otfs[*,*,step,lamb]=otf
         endfor ;step
      endfor ;lamb
   endelse 

; ------------------------------------------------------------------------------
; Finalize
; ------------------------------------------------------------------------------

save, phi,otfs, filename=info.saves+info.files(progind)+'.sav'


end

