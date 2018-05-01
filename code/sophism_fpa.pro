pro sophism_fpa

; ==============================================================================
;
; DESCRIPTION
;    Simulate detector readout incl. resampling, noise, ...
;
; CALLED FROM
;    sophism
;
; SUBROUTINES
;    headfits
;
; MAIN INPUTS
;
;
; OUTPUT
;    Detector images with noises (darks, flats,...) added.
;    gain_table      Gain table (flat field)
;    dc_im           Dark current
;
; VERSION HISTORY 
;   Alex Feller   v0.2    2010-12-01
;   J. Piqueras           2011-10-01
;   J. Blanco. 
;   J. Blanco. 2012. v2_0. Added gain table, conversion to photons,
;      dual-beam mode. 
;   J. Blanco. Dec 2012. v2_1. Display images in verbose mode.   
;   J. Blanco. Jan 2013. v2_2. Flat generation and fringes added.
;   V. Martinez, J. Blanco. Feb 2013. v2_3. Corrections to photon and
;      readout noise calculations. Inclusion of transmittance. Option
;      to not include photon noise and dark.
;   J. Blanco. Jan 2014. v2_4. Changed flat generation, simplified and
;      more variation in the FOV. Hinted options to aproximate max or
;      min of the final flat (normalized units)
;   J. Blanco. May 2014. v2_5. Changed fwhm_real for its
;      mean. Re-order the noises generation.
;   J. Blanco. Jun 2014. v2_6. Modified flat generation, option to
;      give a max-min range of final flat. Added Loading/Saving photon
;      noise option.
;   J. Blanco. Jul 2014. v2_7. Corrected ntim according to new cycrep
;      scheme. Transmittance changed, in case of dual-beam, to 50% in
;      each detector. Corrected bug in fringes generation. Loading
;      option for gain table, modified variable name.
;   J. Blanco. Sep 2014. v2_8. Removed the setting of negative values to
;      0 in the case of no modulation. Added double keyword to the
;      photon noise generation
;   J. Blanco. Nov 2014. v2_9. Introduced dead/hot pixel
;      generation. Added cosmic rays generator 
;   J. Blanco. May 2015. v2_95. Introduced factor for flux conversion
;      for input data, when it is not in erg/s/cm2/...
;   J. Blanco. Feb 2016. v2_96. Corrected gain_table generation for
;      'None' case
;   J. Blanco. Apr 2016. v2_97. New fringes generation, based on
;      theory, lambda-dependent
;   J. Blanco. Nov 2016. v2_98. Correction of bugs when loading noises
;      in dual-beam
;   J. Blanco. Apr 2017. v2_99. Added ndet=1 for the case of no data series
;
; ==============================================================================


; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load settings
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print,''
print,'sophism_fpa'
print,''

restore,'../settings/settings_sophism.sav'

progind=6
if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if (progma eq -1) then progma=0

;to check on things regarding the modules the data have passed
head=headfits(info.saves+info.files(progma)+'_0.fits*')
headmod=sxpar(head,'HISTORY')
etal=where(headmod eq 'Etalon Module')
poli=where(headmod eq 'Polarization Module')

; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Module configuration
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Transmittance of the system, to convert to 50% in case of dual-beam
if (info.dualbeam eq 1) then transsys=info.transsys*0.5 else transsys=info.transsys

;  Shutter type (snapshot or rolling)
shutter_type = info.fpashutt

;careful!! We are not considering interpolation case. The solar
;change during shutter sohuld be small enough?
; Rolling shutter interpolation (on or off)
roll_interpol = 'off'

;temperature
fpa_temp = info.fpa_temp

;  [krad] Total Ionization Dose
;careful!! Right now, fixed parameter
fpa_tid = 0

; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Module derived parameters
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; number of samples per detector exposure
nexp = fix(round(info.etim*1e-3/info.tsamp))

; number of samples per detector frame.
; It may be different from nexp, because this is the number of samples
; corresponding to frame time, which is the largest of exposure time
; vs readout time.
nframe = info.nframe

;tdeaths is no longer an issue, since those frames have been already
;removed during modulation or there is no modulation at all, no need
;to wait for LCVR states
;nmcycl=info.cycrep

;  Total number of frames read from the detector
ndet=info.modnbuff*info.nacc*info.etalpos ;info.modnbuff*info.nacc*nmcycl*info.etalpos  
;if (etal eq -1 AND poli eq -1) then ndet=ndet*info.obsset
if (etal eq -1) then ndet=ndet*info.obsset
if (info.dataser eq 0) then ndet=1

; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Module specific data: Definition of the FPA
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;careful!! A lot of this, fixed parameters. Should any of them be input?

; [DN/e-] FPA conversion gain 
fpa_gain = 0.023
; Dark Current
;  [deg C] DC vs T, temperature vector
dc_temp = [-20,-10,0,10,20,30,40,50,60]
;  [e-/s] DC vs T, DC values
dc_val = [100,100,100,105,180,185,800,1100,4100] 
 ;  [e-/s/krad] DC TID factor
dc_tid = 20 
; Quantum Efficiency and Fill Factor
;  [Adim] QE x FF @ 617 nm
qe_val = info.fpaquan
; Readout (temporal) noise
;  [e- rms] Average temporal (dark) noise [Note that rms=std because mu=0, substracted before computing]
rout_sigma = info.noisigma

; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load input data
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;DN to photon conversion, based on solar energy input, aperture area, etc.
hh=6.62606957d-27                           ; erg x s
cc=2.9979d10                                ; cm/s
telap=info.telap*1e2                        ; cm
fpaplaterad=info.fpaplate/3600.*(!pi/180.)  ; strad

;if etalon module was run, use the spectral fwhm calculated there and
;wavelength array
;if (max(strmatch(sxpar(head,'HISTORY'),'Etalon Module')) eq 1) then begin
if (etal ne -1) then begin
   restore,info.saves+info.files(3)+'.sav' 
   fwhm=mean(fwhm_real)*1d-8                      ; cm -- fwhm_real is in AA
   lambd=llo*1d-8                           ; cm
endif else begin
   fwhm=85.*1d-11                           ; cm
;   lambd=info.wl*1d-8+(findgen(info.wavepoint)*info.wavesamp*1d-11)-(info.wavepoint/2.*info.wavesamp*1d-11)                 ; cm
   lambd=info.lll*1d-8
endelse

;considering the MHD data with proper units,i.e. erg/s/cm2/strad/cm
dnphconv=(!pi*(telap/2.)^2.)*(fpaplaterad^2.)*fwhm/(hh*cc/lambd)
dnphconv=dnphconv*info.etim*1e-3

;include in the conversion to photons the factor for the flux, in case
;input data was normalized or so. It should be around 3e14 erg/s/cm2/...
if (info.fpafluxconv eq 1) then dnphconv=dnphconv*info.fpafluxfactor

;Dark Current
dc_im=fltarr(ndet,info.sz,info.sz)
;load dark if selected
if (info.fpadarksel eq 'Load') then begin
   if (info.talk eq 1) then print,'Loading dark'
   restore,info.saves+info.fpadarkname
endif 
;start dark generation if selected
if (info.fpadarksel eq 'Generate') then begin
   if (info.talk eq 1) then print,'Generating dark'
   dc = dc_val[where(min(abs(dc_temp-fpa_temp)) eq abs(dc_temp-fpa_temp))] ;  [e-/s]
   dc = dc + fpa_tid * dc_tid   ;  [e-/s] 
   dc_mean = (info.etim*1e-3)*dc ;  [e-]
   seed = systime(/s)
   seed = long((seed-long(seed))*1.e6)
 ; Note: /POISSON could be Gaussian if dc_mean is not very low
;   dc_im = randomn(seed, sizy[3], sizy[4], POISSON=dc_mean[0])
   dc_im = randomn(seed, ndet, info.sz,info.sz, POISSON=dc_mean[0])
endif ;generate dark

; gain table detector 1
gain_table=fltarr(info.sz,info.sz)+1.
if (info.fpagainsel eq 'Generate') then gain_table=info.fpagain*randomn(seed,info.sz,info.sz)+1.
if (info.fpagainsel eq 'Load') then begin
   if (info.talk eq 1) then print,'Loading gain table'
   restore,info.saves+info.fpagainname
endif

; FPA flat
flat=fltarr(1,1,info.sz,info.sz)+1.
;load flat if selected
if (info.fpaflatsel eq 'Load') then begin
   if (info.talk eq 1) then print,'Loading flat'
   restore,info.saves+info.fpaflatname
endif
;start flat generation if selected
if (info.fpaflatsel eq 'Generate') then begin
   if (info.talk eq 1) then print,'Generating flat'
   nflam=1
   nfpol=1
;number of flats to be generated if selected for each wavelength
;and/or polarization state
   if (info.fpaflatmore(0) eq 1) then nflam=info.etalpos
   if (info.fpaflatmore(1) eq 1) then nfpol=info.modnbuff
;for the case of flats changing in time. Probably some relation
;between them needed, not completely random?
         ;if (info.fpaflatmore(2) eq 1) then nfl=ndet
;define the flat
   flat=fltarr(nflam,nfpol,info.sz,info.sz)
;;define arrays for horizontal and vertical for the surface coefficients
;   xaxis=fltarr(info.sz,info.sz)
;   yaxis=xaxis
;   for ll=0,info.sz-1 do xaxis[*,ll]=findgen(info.sz)
;   for ll=0,info.sz-1 do yaxis[ll,*]=findgen(info.sz)

   for ff=0,nflam-1 do begin
      for fff=0,nfpol-1 do begin
;generate random image
         prefl=randomu(seed,info.sz,info.sz)
;fit a surface ;;;to obtain polynomial coefficients
         prefl=sfit(prefl/mean(prefl),3) ;info.fpaflorder,kx=kxx)
;;'amplify' coefficients to get higher effects
;         prefl=fltarr(info.sz,info.sz)
;         for ii=0,info.fpaflorder do begin
;            for jj=0,info.fpaflorder do prefl=prefl+kxx(ii,jj)*xaxis^jj*yaxis^ii*fix(randomu(seed)*10.)
;         endfor 
;to increase the range over the FOV, take the power 3 of the surface,
;scaling with the size of the image array respect to 100 pixels
;        if (info.fpaflmax eq 0 AND info.fpaflmin eq 0) then $
         prefl=prefl^(3.*round(info.sz/100.))
;possibility to set a maximum or minimum of the range for the final
;flat
;check if normalizaton afterwards would strongly modify the
;value. Probably a bit but...
;        if (info.fpaflmax ne 0) then prefl=prefl^(alog(info.fpaflmax)/alog(max(prefl)))
;        if (info.fpaflmin ne 0 AND info.fpaflmax eq 0) then prefl=prefl^(alog(info.fpaflmax)/alog(min(prefl)))         
;to stablish a max-min range for the flat
         if (info.fpaflrang ne 0) then begin
            xx=findgen(100)/10.
            rangn=max(prefl)^xx-min(prefl)^xx
            minrang=min(abs(rangn-info.fpaflrang),minrangind)
            prefl=prefl^xx(minrangind)
;maybe could try to reach a compromise between target range and target
;minimum (or maximum)? But could be far from both
;           xx=findgen(100)/10.
;           rangn=max(prefl)^xx-min(prefl)^xx
;           minin=min(prefl)^xx
;           difrang=rangn-info.fpaflrang
;           difmin=minin-info.fpaflmin
;           mintot=min(abs(difrang-difmin),mintotind)
;           prefl=prefl^xx(mintotind)
         endif
         prefl=prefl/mean(prefl)
         flat(ff,fff,*,*)=prefl
      endfor  
   endfor 
endif ;generate flat
sizfl=size(flat)


;add fringes
;check the phase change to different etalon and polarization. And
;fringes constant or always will change? 4dimensions like flat?
nfdet=1
nflam=1
nfpol=1
;fffffffffffffffffffffffffffffffffffffffffffff
;fring=0
fring=fltarr(nfdet,nflam,nfpol,info.sz,info.sz)
;;fring=flat*0.

;if (info.fpafring eq 1) then begin 
;;   nfdet=1
;;   nflam=1
;;   nfpol=1
;   fringdirmod=(info.fpafringdir mod 180)
;   if (fringdirmod le 45) then fringdir=fringdirmod/45.
;   if (fringdirmod gt 45 AND fringdirmod lt 135) then fringdir=(2.-fringdirmod/45.)
;   if (fringdirmod ge 135) then fringdir=(fringdirmod/45.-4.)
;;   if (info.fpafringtime eq 1) then begin 
;;      nflam=info.etalpos
;;      nfpol=info.modnbuff
;;      fring=fltarr(nflam,nfpol,info.sz,info.sz)
;;   endif  
;   if (total(info.fpafringtime) ge 1) then begin
;depending on fringes variation, create different ones for frames,
;lambdas and/or polarizations
;In this case, use etalpos and modnbuff, because fringes should only
;change when etalon/LCVRs change. So if no etalon/modulation module is
;selected, apply same fringes to whole wavelength/stokes dimension 
;      if (info.fpafringtime[0] eq 1) then nfdet=ndet
;      if (info.fpafringtime[1] eq 1) then nflam=info.etalpos
;      if (info.fpafringtime[2] eq 1) then nfpol=info.modnbuff
;   endif 
;   fring=fltarr(nfdet,nflam,nfpol,info.sz,info.sz)
;   for nff=0,nfdet-1 do begin
;y cambiar por ej. la direccion segun sean franjas de etalon o de LCVRs?
;      for ff=0,nflam-1 do begin
;         for fff=0,nfpol-1 do begin
;            pha=randomu(seed)       
;            if (fringdirmod gt 45 AND fringdirmod lt 135) then begin 
;               for j=0,info.sz-1 do fring[nff,ff,fff,j,*]=sin((findgen(info.sz)+j*fringdir)/(!pi*info.fpafringwidth)+!pi*pha)*info.fpafringamp/100. 
;            endif else begin 
;               for j=0,info.sz-1 do fring[nff,ff,fff,*,j]=sin((findgen(info.sz)+j*fringdir)/(!pi*info.fpafringwidth)+!pi*pha)*info.fpafringamp/100.
;            endelse 
;         endfor 
;      endfor 
;   endfor 
;endif
;fring=fring+1.

if (info.fpafring eq 1) then begin
;model for spatial variation of the fringing-element
   szlambda=(size(lambd))(1)
   L=1e-8
   a=70/!radeg
;   x=fltarr(200,200)
;   for jj=0,199 do x[*,jj]=findgen(200)
   x=fltarr(info.sz,info.sz)
   for jj=0,info.sz-1 do x[*,jj]=findgen(info.sz)
   y=transpose(x)
   xc=-info.sz/2
   yc=info.sz/2
   ;d=L*(x*cos(a)-y*sin(a))+d
   ddd=L*((x-xc)^2.-(y-yc)^2.)+info.fpafrinsizd
;stop
;fringes equations
;   phi=(4*!pi*d*n*cos(theta))/lambd
;   frin=(r^2-1)^2/(r^4-2*r^2*cos(2*phi)+1)
;right now, only variation with lambda, not polarization or time
   phi=fltarr(1,szlambda,1,info.sz,info.sz)
   for uu=0,szlambda-1 do phi[0,uu,0,*,*]=(4*!pi*ddd*info.fpafrinn*cos(info.fpafrintheta))/lambd[uu]
   fring=(info.fpafrinref^2-1)^2/(info.fpafrinref^4-2*info.fpafrinref^2*cos(2*phi)+1)
endif else fring=fring+1.
;fffffffffffffffffffffffffffffffffffffffffffff
sizfr=size(fring)
ipol_old=0
polis=0

;make dead pixels as 0
maskdead=fltarr(info.sz,info.sz)+1
if (info.fpadeadpix gt 0) then begin
   for pixi=0,info.fpadeadpix-1 do maskdead[randomu(seed)*(info.sz-1),randomu(seed)*(info.sz-1)]=0.
endif 

;for the hot pixels, provide saturation level as setting? Guess so
maskhot=fltarr(info.sz,info.sz)
if (info.fpahotpix gt 0) then begin
   for pixi=0,info.fpahotpix-1 do maskhot[randomu(seed)*(info.sz-1),randomu(seed)*(info.sz-1)]=1
   wwhot=where(maskhot gt 0)
endif 


; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Detector output
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (info.talk eq 1) then begin
   print,'FPA detector simulation'
   print,'Join time samples in exposure; adding noises'
endif 


;///////
;defining array to hold photon noise for saving/loading
;should take all lambdas and pols, in case some of it has not been run
if (info.fpaphotnoi eq 'Generate') then begin
   lamax=fix(sxpar(head,'NAXIS1'))
   polax=fix(sxpar(head,'NAXIS2'))
   photnoise=fltarr(info.sz,info.sz,lamax,polax,ndet)
endif
if (info.fpaphotnoi eq 'Load') then begin
   restore,info.saves+info.fpaphotnoiname
;reform the photon noise array in case of various observation sets,
;because that option changes the wavelength dimension
   ima=readfits(info.saves+info.files(progma)+'_0.fits*',head,/comp,/sil)
   sizy=size(ima)
   szph=size(photnoise)
   if (etal ne -1 AND info.obsset gt 1 AND szph(3) ne sizy(1)) then begin
      if (szph(3)*info.obsset ne sizy(1)) then begin
         print,''
         print,'Error. Photon noise loaded not compatible with given simulation, different spectral positiones selected.'
         print,''
         stop
      endif else begin
         photpre=photnoise
         photnoise=fltarr(sizy(3),sizy(4),sizy(1),sizy(2),ndet)
         for ooo=0,info.obsset-1 do photnoise(*,*,ooo*szph(3):szph(3)*(ooo+1)-1,*,*)=photpre
      endelse 
   endif 
    
endif

; Loop for every detector frame
for f=0,ndet-1 do begin
  
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Exposure and shutter implementation
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   case shutter_type of
      'snapshot': begin
;case snapshot shutter, just add all the time samples forming a frame,
;i.e. falling into exposure time (or readout time, the largest)
         for i=0,nexp-1 do begin 
;only take into account nexp samples. If nframe is larger than nexp,
;those are lost
            ima=readfits(info.saves+info.files(progma)+'_'+strtrim(fix(f*nframe)+i,2)+'.fits*',head,/comp,/sil)
            sizy=size(ima)
            if (f eq 0 AND i eq 0) then imout=ima else imout+=ima
         endfor
      end
      'rolling': begin
         case roll_interpol of
;careful!! Since, in principle, this first part with Interpolation ON 
;is not used, it hasn't been corrected for a long time
            'on': begin
               print,'Not really functional now. Sorry.'
               stop
               tmp_index = [floor(indgen(nframe)*(im_sz[4]-1)/(nframe-1))]
               for i=0,im_sz[3]-1 do begin
                  for j=0,im_sz[4]-1 do begin
                     tmp_inputs = [tseries[f*nframe:(f+1)*nframe-1].i[i,j]]
                     imout[f].i[i,j] = total(interpol(tmp_inputs, tmp_index, j+indgen(nexp)))
                  endfor
               endfor
            end

            'off': begin
; number of rows (of the final detector image) per input image (time sample)
               factor_im = info.sz/nframe 
               for i=0,nframe-1 do begin
                  ima=readfits(info.saves+info.files(progma)+'_'+strtrim(fix(f*nframe)+i,2)+'.fits*',head,/comp,/sil)
                  if (f eq 0 AND i eq 0) then begin
                     sizy=size(ima)
;define a variable like the data but with extra temporal
;dimension for time samples of a frame time. Those will be added to
;form a single exposure of exposure time
                     tseries=fltarr(nframe,sizy(1),sizy(2),sizy(3),sizy(4))
                  endif
                  tseries[i,*,*,*,*]=ima 
               endfor 
               imout=ima*0.
;loop over rows. After factor_im, changes to next time sample within
;the frame. This represents a discrete evolution of the rolling shutter.
;Around the border of factor_im, both previous and following time
;sample are added to create a mixing
               for j=0,sizy[4]-1 do begin
                  tmp_index = ([floor((j+indgen(nexp))/factor_im)])+f*nframe < ((f+1)*nframe-1)
                  imout[*,*,*,j]=total(tseries[tmp_index,*,*,*,j],1)
               endfor
            end
            else: message, 'Invalid interpolation option'
         endcase
      end                       ;rolling
      else: message, 'Invalid shutter type'
   endcase

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Photon conversion
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;conversion from MHD units to photons, including transmittance
;since nexp samples were added, it has to be normalized by that
  for lam=0,sizy(1)-1 do imout[lam,*,*,*] = transsys * imout[lam,*,*,*]*dnphconv[lam]/nexp  ;  [phot]
;stop
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; FPA flat & fringes
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ffr=(f mod sizfr[1])
;to use flats and/or fringes of different modulation states with the
;corresponding data
  ipol=floor(f/info.nacc)
;to reset to 0 the modulation flat/fringes after the 4 (or 2 in
;lineal) states, use mod
;if no modulation, it will only use one flat/fringes due to
;modulation, because LCVRs are not changing
  if (ipol gt ipol_old) then begin
     polis=(polis+1) mod info.modnbuff
     ipol_old=ipol
  endif
;basically, in case flats or fringes have no variation with modulation
  polfr=(polis mod sizfr[3])
  polfl=(polis mod sizfl[2])
;need to loop in lambda because ends up as 2D arrays. And in
;polarization in case no modulation was done (if it was, that
;dimension is 1 element)
  for lam=0,sizy[1]-1 do begin 
;basically, in case flats or fringes have no variation with lambda
     lamfl=(lam mod sizfl[1])
     lamfr=(lam mod sizfr[2])
     for pol=0,sizy[2]-1 do begin
;        polfl=(pol mod sizfl[2])
;        polfr=(pol mod sizfr[3])
;stop
        imout[lam,pol,*,*]=reform(imout[lam,pol,*,*])*reform(flat[lamfl,polfl,*,*])*reform(fring[ffr,lamfr,polfr,*,*])
     endfor 
  endfor
;stop

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; FPA gain table
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  for lam=0,sizy[1]-1 do begin 
     for pol=0,sizy[2]-1 do imout[lam,pol,*,*]=imout[lam,pol,*,*]*gain_table
  endfor

;stop
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; QExFF (ph -> e-)
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;conversion to electrons
  imout=imout*qe_val  ; [e-]

;stop
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Photon (electron) Noise
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (info.fpaphotnoi eq 'Generate') then begin
;original image, to obtain later only the photon noise and save
    imoutpre=imout
     for lam=0,sizy[1]-1 do begin 
        for pol=0,sizy[2]-1 do begin
           if (poli eq -1) then begin
;if there are pixels with value </= 0 (FDT out of disk,...), poisson doesn't
;work. Give those pixels 1e-6 value (poisson gives 0). Then remove again?
              iii=reform(imout[lam,0,*,*])
              zero=where(iii le 0)
              if (max(zero) ne -1) then iii[zero]=1e-6
;in case data is not modulated, use only I.
;              for xx=0,sizy[3]-1 do for yy=0,sizy[4]-1 do imout[lam,pol,xx,yy]=imout[lam,pol,xx,yy]+(randomn(seed,1,poisson=reform(imout[lam,0,xx,yy]))-imout[lam,0,xx,yy])
if (lam eq 0 AND pol eq 0) then print,'Careful! Check this makes sense!'
              for xx=0,sizy[3]-1 do for yy=0,sizy[4]-1 do imout[lam,pol,xx,yy]=imout[lam,pol,xx,yy]+(randomn(seed,1,poisson=reform(iii[xx,yy]),/double)-iii[xx,yy])
           endif else begin
;if there are pixels with value </= 0 (FDT out of disk,...), poisson doesn't
;work. Give those pixels 1e-6 value (poisson gives 0). Then remove again?
              zero=where(imout le 0)
              if (max(zero) ne -1) then imout[zero]=1e-6
              for xx=0,sizy[3]-1 do for yy=0,sizy[4]-1 do imout[lam,pol,xx,yy]=randomn(seed,1,poisson=reform(imout[lam,pol,xx,yy]),/double)
           endelse 
;storing the noise in a variable
;should take all lambdas and pols, in case some of it has not been run
           photnoise[*,*,lam,pol,f]=imout[lam,pol,*,*]-imoutpre[lam,pol,*,*]
        endfor ;pol
     endfor  ;lam
  endif 
;loadphot:
  if (info.fpaphotnoi eq 'Load') then begin
     for lam=0,sizy[1]-1 do begin
        for pol=0,sizy[2]-1 do begin
           imout[lam,pol,*,*]=imout[lam,pol,*,*]+reform(photnoise[*,*,lam,pol,f])
        endfor 
     endfor
  endif
;stop
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read (dark) noise + Dark current
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  rnoise_im = randomn(seed, sizy[3], sizy[4], /NORMAL) * rout_sigma  ;  [e-]
  for lam=0,sizy[1]-1 do begin 
     for pol=0,sizy[2]-1 do imout[lam,pol,*,*]=imout[lam,pol,*,*]+dc_im[f,*,*]+randomn(seed, sizy[3], sizy[4], /NORMAL) * rout_sigma
  endfor

;in case of no modulation, don't set to 0 the negative values, of course
  if (poli ne -1) then begin
     tmp = where(imout[*,*,*,*] lt 0)
;  Negative values are set to zero
     if(tmp[0] ne -1) then imout[tmp] = 0
  endif

; At this point, the bias correction would be already done, because the average value of the bias image has not been added, only the std (rms)

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; DSNU, PRNU, FPN...
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; FPA output chain (e- (--> V) --> DN)
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  imout[*,*,*,*] = fpa_gain * imout[*,*,*,*]

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Cosmic Rays
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (info.cosmic eq 1) then begin
;mostly derived from 'Image simulations for Heliospheric Imager
;for the SECCHI/STEREO project' by Davis & Harrison 2003
;number of cosmic rays per frame. Calculate the total hits of cosmic
;rays for the whole observation from the hits/s and the total time,
;and divide by the number of frames to get an average cosmic hits per
;frame. This includes total frames, also the ones discarded since they
;are part of the observation although not stored.
;May be misleading to use the same for snapshot/rolling shutter, since
;the former only works during texp, the latter for the whole frtim, but..
     ncosm=round(info.cosmhit*info.tobs/(info.ntim/info.nframe))
;generate length of events as random exponential distribution (more
;probable short traces), to a maximum of 25 pixels length, rounding
;all up so no 0 length appears
;the maximum length should be assessed, it's just a made up number of
;25 pixels now (it should depend on the detector material and
;thickness?). And it's made so that there is always some 25
;pixel event, should re-check.
;Change to a 5 percent lenght of the image size for cosmic ray length
;     cosmlen=abs(randomn(seed,ncosm))
     cosmlen=randomu(seed,ncosm,gamma=1)
     cosmlen=ceil(cosmlen/max(cosmlen)*info.sz*0.05);25.)

     for ncc=0,ncosm-1 do begin 
;generate the cosmic ray traces, with random angles of incidence
;place them in random positions in the FOV
;Limit to 90 degrees because 180 could get into negative coordinates
;and it is probably not so important the 90-180 direction
;        cosmray=findgen(cosmlen[ncc])*cos(randomu(seed)*180./!radeg)
        ang=randomu(seed)*90./!radeg
        xcosm=findgen(cosmlen[ncc])*cos(ang)+randomu(seed)*info.sz
        ycosm=findgen(cosmlen[ncc])*sin(ang)+randomu(seed)*info.sz
        if (max(xcosm) gt info.sz) then begin
           ww=max(where(xcosm le info.sz))
           xcosm=xcosm(0:ww)
           ycosm=ycosm(0:ww)
           cosmlen(ncc)=(size(xcosm))(1)
        endif 
        if (max(ycosm) gt info.sz) then begin
           ww=max(where(ycosm le info.sz))
           xcosm=xcosm(0:ww)
           ycosm=ycosm(0:ww)
           cosmlen(ncc)=(size(ycosm))(1)
        endif 

;set the cosmic traces intensity to saturation, although they may vary
;in reality.
;place them in random positions in the FOV and in all
;wavelengths/modulations.
;Wavelengths will be selected in the accumulation module so no
;problem. If the etalon was not run, then there is no temporal
;selection of wavelengths, so consider that all are affected by cosmic
;ray.
;If there is modulation, only one element in that dimension, so apply
;cosmic to it. If there is no modulation, that dimension is the actual
;Stokes vector, so apply to all.
;        imout[*,*,findgen(cosmlen[ncc])+randomu(seed)*info.sz,cosmray+randomu(seed)*info.sz]=info.fpasatval
        for nccl=0,cosmlen[ncc]-1 do imout[*,*,xcosm[nccl],ycosm[nccl]]=info.fpasatval
;stop
     endfor ;ncc

  endif ;cosmic

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Dead/Hot pixels
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;done like this, hot and dead pixels are at the same positions for all
;frames and accumulations. But then hot pixels are added when
;accumulating, so they go off the roof. Divide by number of
;accumulations? Not, no?
  for ll=0,sizy[1]-1 do begin
     for pp=0,sizy[2]-1 do begin
        imhot=reform(imout[ll,pp,*,*])
        if (info.fpahotpix gt 0) then imhot[wwhot]=info.fpasatval
        imout[ll,pp,*,*]=imhot
        imout[ll,pp,*,*]=imout[ll,pp,*,*]*maskdead
     endfor
  endfor

  sxaddpar,head,"HISTORY",'FPA Module'
  sxaddpar,head,'TRANSSYS',format='(F6.3)',transsys,'Transmission of the system'
  sxaddpar,head,'SHUTTER',shutter_type,'Type of Shutter'
  sxaddpar,head,'QUANEFF',format='(F8.3)',qe_val,'Quantum Efficiency'
  sxaddpar,head,'ROUTSIG',format='(F8.3)',rout_sigma,'Readout Noise Sigma'
;stop
  if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(f,2)+'.fits',imout,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(f,2)+'.fits',imout,head
  imout=imout*0.

  if (info.talk eq 1) then print,format='(21(%"\b"),"Detector frames ", I3," %",$)',(f+1)/float(ndet)*100.

endfor ;ndet
print,''

save,gain_table,dc_im,filename=info.saves+info.files(progind)+'.sav'
save,flat,fring,filename=info.saves+info.files(progind)+'_flat.sav'
;saving photon noise
save,photnoise,filenam=info.saves+info.files(progind)+'_phot.sav'


;display images of dark and gain table if selected
if (info.talk eq 1) then begin
   window,xs=600,ys=600,title='Dark and Gain table',/free
   if (info.dualbeam eq 0) then !p.multi=[0,3,1] else !p.multi=[0,3,2]
   tvframe,dc_im[0,*,*],/bar,/aspect,btit='Dark'
   tvframe,gain_table,/bar,/aspect,btit='Gain Table'
   tvframe,flat[0,0,*,*]*fring[0,0,0,*,*],/bar,/aspect,btit='Flat'
endif



; ///////////////////////////////
; Dual-beam
;////////////////////////////////

if (info.dualbeam eq 1) then begin
   if (info.talk eq 1) then print,'FPA dual-beam detector simulation'

;Dark Current
   dc_im=fltarr(ndet,info.sz,info.sz)
;load dark if selected
   if (info.fpadarksel eq 'Load') then begin
      if (info.talk eq 1) then print,'Loading dark'
      posi=strpos(info.fpadarkname,'.sav')
      if (posi eq -1) then begin
         print,'File not found'
         stop
      endif else restore,info.saves+strmid(info.fpadarkname,0,posi)+'_2'+strmid(info.fpadarkname,posi)
   endif 
;start dark generation if selected
   if (info.fpadarksel eq 'Generate') then begin
      if (info.talk eq 1) then print,'Generating dark'
      dc = dc_val[where(min(abs(dc_temp-fpa_temp)) eq abs(dc_temp-fpa_temp))] ;  [e-/s]
      dc = dc + fpa_tid * dc_tid ;  [e-/s] 
      dc_mean = (info.etim*1e-3)*dc ;  [e-]
      seed = systime(/s)
      seed = long((seed-long(seed))*1.e6)
 ; Note: /POISSON could be Gaussian if dc_mean is not very low
;   dc_im = randomn(seed, sizy[3], sizy[4], POISSON=dc_mean[0])
      dc_im = randomn(seed, ndet, info.sz,info.sz, POISSON=dc_mean[0])
   endif ;generate dark

;gain table
   gain_table=fltarr(info.sz,info.sz)+1.
   if (info.fpagainsel eq 'Generate') then gain_table=info.fpagain*randomn(seed,info.sz,info.sz)+1.
   if (info.fpagainsel eq 'Load') then begin
      posi=strpos(info.fpagainname,'.sav')
      if (posi eq -1) then begin
         print,'File not found'
         stop
      endif else restore,info.saves+strmid(info.fpagainname,0,posi)+'_2'+strmid(info.fpagainname,posi)
   endif 

; FPA flat
   flat=fltarr(1,1,info.sz,info.sz)+1.
;load flat if selected
   if (info.fpaflatsel eq 'Load') then begin
      posi=strpos(info.fpaflatname,'.sav')
      if (posi eq -1) then begin
         print,'File not found'
         stop
      endif else restore,info.saves+strmid(info.fpaflatname,0,posi)+'_2'+strmid(info.fpaflatname,posi)
   endif
;start flat generation
   if (info.fpaflatsel eq 'Generate') then begin
;   nfl=1
      nflam=1
      nfpol=1
;number of flats to be generated if selected for each wavelength
;and/or polarization state
      if (info.fpaflatmore(0) eq 1) then nflam=info.etalpos
      if (info.fpaflatmore(1) eq 1) then nfpol=info.modnbuff
;for the case of flats changing in time. Probably some relation
;between them needed, not completely random?
         ;if (info.fpaflatmore(2) eq 1) then nfl=ndet
;define the flat
      flat=fltarr(nflam,nfpol,info.sz,info.sz)
;define arrays for horizontal and vertical for the surface coefficients
;      xaxis=fltarr(info.sz,info.sz)
;      yaxis=xaxis
;      for ll=0,info.sz-1 do xaxis[*,ll]=findgen(info.sz)
;      for ll=0,info.sz-1 do yaxis[ll,*]=findgen(info.sz)
;   for ff=0,nfl-1 do begin
      for ff=0,nflam-1 do begin
         for fff=0,nfpol-1 do begin
;generate random image
            prefl=randomu(seed,info.sz,info.sz)
;fit a surface to obtain polynomial coefficients
            prefl=sfit(prefl/mean(prefl),3) ;,info.fpaflorder,kx=kxx)
;'amplify' coefficients to get higher effects
;            prefl=fltarr(info.sz,info.sz)
;            for ii=0,info.fpaflorder do begin
;               for jj=0,info.fpaflorder do prefl=prefl+kxx(ii,jj)*xaxis^jj*yaxis^ii*fix(randomu(seed)*10.)
;            endfor 
            prefl=prefl^(3.*round(info.sz/100.))
            if (info.fpaflrang ne 0) then begin
               xx=findgen(100)/10.
               rangn=max(prefl)^xx-min(prefl)^xx
               minrang=min(abs(rangn-info.fpaflrang),minrangind)
               prefl=prefl^xx(minrangind)
            endif
            prefl=prefl/mean(prefl)
            flat(ff,fff,*,*)=prefl
         endfor 
      endfor 
   endif ;generate flat
   sizfl=size(flat)

;add fringes
;check the phase change to different etalon and polarization. And
;fringes constant or always will change? 4dimensions like flat?
   nfdet=1
   nflam=1
   nfpol=1
   fring=fltarr(nfdet,nflam,nfpol,info.sz,info.sz)
;fffffffffffffffffffffffffffffffffffffffffffff
;;fring=flat*0.
;   if (info.fpafring eq 1) then begin 
;;   nflam=1
;;   nfpol=1
;      fringdirmod=(info.fpafringdir mod 180)
;      if (fringdirmod le 45) then fringdir=fringdirmod/45.
;      if (fringdirmod gt 45 AND fringdir lt 135) then fringdir=(2.-fringdirmod/45.)
;      if (fringdirmod ge 135) then fringdir=(fringdirmod/45.-4.)
;;   if (info.fpafringtime eq 1) then begin 
;;      nflam=info.etalpos
;;      nfpol=info.modnbuff
;;      fring=fltarr(nflam,nfpol,info.sz,info.sz)
;;   endif  
;      if (total(info.fpafringtime) ge 1) then begin
;         if (info.fpafringtime[0] eq 1) then nfdet=ndet
;         if (info.fpafringtime[1] eq 1) then nflam=info.etalpos
;         if (info.fpafringtime[2] eq 1) then nfpol=info.modnbuff
;      endif 
;      fring=fltarr(nfdet,nflam,nfpol,info.sz,info.sz)
;      for nff=0,nfdet-1 do begin
;         for ff=0,nflam-1 do begin
;            for fff=0,nfpol-1 do begin
;               pha=randomu(seed)       
;;         for j=0,info.sz-1 do fring[ff,fff,*,j]=sin((findgen(info.sz)+j*info.fpafringdir/45.)/(!pi*info.fpafringwidth)+!pi*pha)*info.fpafringamp/100.
;               if (fringdirmod gt 45 AND fringdirmod lt 135) then begin 
;                  for j=0,info.sz-1 do fring[ff,fff,j,*]=sin((findgen(info.sz)+j*fringdir)/(!pi*info.fpafringwidth)+!pi*pha)*info.fpafringamp/100. 
;               endif else begin 
;                  for j=0,info.sz-1 do fring[ff,fff,*,j]=sin((findgen(info.sz)+j*fringdir)/(!pi*info.fpafringwidth)+!pi*pha)*info.fpafringamp/100.
;               endelse 
;            endfor 
;         endfor 
;      endfor 
;   endif
;   fring=fring+1.

   if (info.fpafring eq 1) then begin
;model for spatial variation of the fringing-element
      szlambda=(size(lambd))(1)
      L=1e-8
      a=70/!radeg
;   x=fltarr(200,200)
;   for jj=0,199 do x[*,jj]=findgen(200)
      x=fltarr(info.sz,info.sz)
      for jj=0,info.sz-1 do x[*,jj]=findgen(info.sz)
      y=transpose(x)
      xc=-info.sz/2
      yc=info.sz/2
   ;d=L*(x*cos(a)-y*sin(a))+d
      ddd=L*((x-xc)^2.-(y-yc)^2.)+info.fpafrinsizd
;stop
;fringes equations
;   phi=(4*!pi*d*n*cos(theta))/lambd
;   frin=(r^2-1)^2/(r^4-2*r^2*cos(2*phi)+1)
;right now, only variation with lambda, not polarization or time
      phi=fltarr(1,szlambda,1,info.sz,info.sz)
      for uu=0,szlambda-1 do phi[0,uu,0,*,*]=(4*!pi*ddd*info.fpafrinn*cos(info.fpafrintheta))/lambd[uu]
      fring=(info.fpafrinref^2-1)^2/(info.fpafrinref^4-2*info.fpafrinref^2*cos(2*phi)+1)
   endif else fring=fring+1.
;fffffffffffffffffffffffffffffffffffffffffffff
   sizfr=size(fring)
   ipol_old=0
   polis=0


;make dead pixels as 0
   maskdead=fltarr(info.sz,info.sz)+1
   if (info.fpadeadpix gt 0) then begin
      for pixi=0,info.fpadeadpix-1 do maskdead[randomu(seed)*(info.sz-1),randomu(seed)*(info.sz-1)]=0.
   endif 

;for the hot pixels, provide saturation level as setting? Guess so
   maskhot=fltarr(info.sz,info.sz)
   if (info.fpahotpix gt 0) then begin
      for pixi=0,info.fpahotpix-1 do maskhot[randomu(seed)*(info.sz-1),randomu(seed)*(info.sz-1)]=1
      wwhot=where(maskhot gt 0)
   endif 


; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Detector output
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if (info.fpaphotnoi eq 'Generate') then begin
      lamax=fix(sxpar(head,'NAXIS1'))
      polax=fix(sxpar(head,'NAXIS2'))
      photnoise=fltarr(info.sz,info.sz,lamax,polax,ndet)
   endif
   if (info.fpaphotnoi eq 'Load') then begin
      posi=strpos(info.fpaphotnoiname,'.sav')
      if (posi eq -1) then begin
         print,'File not found'
         stop
      endif else restore,info.saves+strmid(info.fpaphotnoiname,0,posi)+'_2'+strmid(info.fpaphotnoiname,posi)
   endif 

   ; Loop for every detector frame
   for f=0,ndet-1 do begin

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Exposure and shutter implementation
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      case shutter_type of
         'snapshot': begin      ; Exposure (etim) + Readout (itim-etim or etim)
            for i=0,nexp-1 do begin 
               ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(fix(f*nframe)+i,2)+'.fits*',head,/comp,/sil)
               sizy=size(ima)
               if (f eq 0 AND i eq 0) then imout=ima else imout+=ima
            endfor
         end
         'rolling': begin       ; Readout time
            case roll_interpol of
               'on': begin      ;  Interpolation is ON              
                  tmp_index = [floor(indgen(nframe)*(im_sz[4]-1)/(nframe-1))]
                  for i=0,im_sz[3]-1 do begin
                     for j=0,im_sz[4]-1 do begin
                        tmp_inputs = [tseries[f*nframe:(f+1)*nframe-1].i[i,j]]
                        imout[f].i[i,j] = total(interpol(tmp_inputs, tmp_index, j+indgen(nexp)))
                     endfor
                  endfor
               end
               'off': begin          ;  Interpolation is OFF
                  factor_im = info.sz/nframe ;  # rows per input image
                  for i=0,nframe-1 do begin
                     ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(fix(f*nframe)+i,2)+'.fits*',head,/comp,/sil)
                     if (f eq 0 AND i eq 0) then begin
                        sizy=size(ima)
                        tseries=fltarr(nframe,sizy(1),sizy(2),sizy(3),sizy(4))
                        imout=ima*0.
                     endif
                     tseries[i,*,*,*,*]=ima 
                  endfor 
                  imout=ima*0.
                  for j=0,sizy[4]-1 do begin
                     tmp_index = ([floor((j+indgen(nexp))/factor_im)])+f*nframe < ((f+1)*nframe-1)
                     imout[*,*,*,j]=total(tseries[tmp_index,*,*,*,j],1)
                  endfor 
               end
               else: message, 'Invalid interpolation option'
            endcase
         end
         else: message, 'Invalid shutter type'
      endcase

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Photon conversion
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      for nlam=0,sizy(1)-1 do imout[nlam,*,*,*] = transsys * imout[nlam,*,*,*]*dnphconv[nlam]/nexp ;  [phot]

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; FPA flat
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ffr=(f mod sizfr[1])
      ipol=floor(f/info.nacc)
      if (ipol gt ipol_old) then begin
         polis=(polis+1) mod info.modnbuff
         ipol_old=ipol
      endif
      polfr=(polis mod sizfr[3])
      polfl=(polis mod sizfl[2])

      for lam=0,sizy[1]-1 do begin 
         lamfl=(lam mod sizfl[1])
         lamfr=(lam mod sizfr[2])
         for pol=0,sizy[2]-1 do begin
;            polfl=(pol mod sizfl[2])
;            polfr=(pol mod sizfr[2])
            imout[lam,pol,*,*]=reform(imout[lam,pol,*,*])*reform(flat[lamfl,polfl,*,*])*reform(fring[ffr,lamfr,polfr,*,*])
         endfor 
      endfor

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; FPA gain table
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      for lam=0,sizy[1]-1 do begin 
         for pol=0,sizy[2]-1 do imout[lam,pol,*,*]=imout[lam,pol,*,*]*gain_table
      endfor

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; QExFF (ph -> e-)
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      imout=qe_val * imout

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Photon (electron) Noise
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      if (info.fpaphotnoi eq 'Generate') then begin
         imoutpre=imout
 ;        if (info.fpaphotnoi eq 1) then begin
         for lam=0,sizy[1]-1 do begin 
            for pol=0,sizy[2]-1 do begin
               if (poli eq -1) then begin
                  iii=reform(imout[lam,0,*,*])
                  zero=where(iii le 0)
                  if (max(zero) ne -1) then iii[zero]=1e-6
                  for xx=0,sizy[3]-1 do for yy=0,sizy[4]-1 do imout[lam,pol,xx,yy]=imout[lam,pol,xx,yy]+(randomn(seed,1,poisson=reform(imout[lam,0,xx,yy]))-imout[lam,0,xx,yy])
               endif else begin
                  zero=where(imout le 0)
                  if (max(zero) ne -1) then imout[zero]=1e-6
                  for xx=0,sizy[3]-1 do for yy=0,sizy[4]-1 do imout[lam,pol,xx,yy]=randomn(seed,1,poisson=reform(imout[lam,pol,xx,yy]))
               endelse 
               photnoise[*,*,lam,pol,f]=imout[lam,pol,*,*]-imoutpre[lam,pol,*,*]
            endfor ;pol 
         endfor ;lam
      endif  
 
      if (info.fpaphotnoi eq 'Load') then begin
         for lam=0,sizy[1]-1 do begin
            for pol=0,sizy[2]-1 do begin
               imout[lam,pol,*,*]=imout[lam,pol,*,*]+reform(photnoise[*,*,lam,pol,f])
            endfor 
         endfor
      endif
     
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read (dark) noise + Dark Current
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      rnoise_im = randomn(seed, sizy[3], sizy[4], /NORMAL) * rout_sigma ;  [e-]
      
      for lam=0,sizy[1]-1 do begin 
         for pol=0,sizy[2]-1 do imout[lam,pol,*,*]=imout[lam,pol,*,*]+dc_im[f,*,*]+randomn(seed, sizy[3], sizy[4], /NORMAL) * rout_sigma
      endfor

      if (poli ne -1) then begin
         tmp = where(imout[*,*,*,*] lt 0)
         if(tmp[0] ne -1) then imout[tmp] = 0
      endif

  ; At this point, the bias correction would be already done, because the average value of the bias image has not been added, only the std (rms)

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; DSNU, PRNU, FPN...
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; FPA output chain (e- (--> V) --> DN)
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      imout[*,*,*,*] = fpa_gain * imout[*,*,*,*]

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Cosmic Rays
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      if (info.cosmic eq 1) then begin
         ncosm=round(info.cosmhit*info.tobs/(info.ntim/info.nframe))
;     cosmlen=abs(randomn(seed,ncosm))
         cosmlen=randomu(seed,ncosm,gamma=1)
         cosmlen=ceil(cosmlen/max(cosmlen)*info.sz*0.05) ;25.)

         for ncc=0,ncosm-1 do begin 
;generate the cosmic ray traces, with random angles of incidence
;place them in random positions in the FOV
;        cosmray=findgen(cosmlen[ncc])*cos(randomu(seed)*180./!radeg)
            ang=randomu(seed)*180./!radeg
            xcosm=findgen(cosmlen[ncc])*cos(ang)+randomu(seed)*info.sz
            ycosm=findgen(cosmlen[ncc])*sin(ang)+randomu(seed)*info.sz
            if (max(xcosm) gt info.sz) then begin
               ww=max(where(xcosm le info.sz))
               xcosm=xcosm(0:ww)
               ycosm=ycosm(0:ww)
               cosmlen(ncc)=(size(xcosm))(1)
            endif 
            if (max(ycosm) gt info.sz) then begin
               ww=max(where(ycosm le info.sz))
               xcosm=xcosm(0:ww)
               ycosm=ycosm(0:ww)
               cosmlen(ncc)=(size(ycosm))(1)
            endif 
            
;        imout[*,*,findgen(cosmlen[ncc])+randomu(seed)*info.sz,cosmray+randomu(seed)*info.sz]=info.fpasatval
            for nccl=0,cosmlen[ncc]-1 do imout[*,*,xcosm[nccl],ycosm[nccl]]=info.fpasatval
         endfor ;ncc

      endif ;cosmic


  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Dead/Hot pixels
  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      for ll=0,sizy[1]-1 do begin
         for pp=0,sizy[2]-1 do begin
            imhot=reform(imout[ll,pp,*,*])
            if (info.fpahotpix gt 0) then imhot[wwhot]=info.fpasatval
            imout[ll,pp,*,*]=imhot
            imout[ll,pp,*,*]=imout[ll,pp,*,*]*maskdead
         endfor
      endfor



      sxaddpar,head,"HISTORY",'FPA Module'
      sxaddpar,head,'SHUTTER',shutter_type,'Type of Shutter'
      sxaddpar,head,'QUANEFF',qe_val,'Quantum Efficiency'
      sxaddpar,head,'ROUTSIG',rout_sigma,'Readout Noise Sigma'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(f,2)+'.fits',imout,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(f,2)+'.fits',imout,head
      imout=imout*0.

      if (info.talk eq 1) then print,format='(26(%"\b"),"Dual-detector frames ", I3," %",$)',(f+1)/float(ndet)*100.

   endfor ;ndet
   print,''

   save,gain_table,dc_im,filename=info.saves+info.files(progind)+'_2.sav'
   save,flat,fring,filename=info.saves+info.files(progind)+'_flat_2.sav'
   save,photnoise,filenam=info.saves+info.files(progind)+'_phot_2.sav'

   if (info.talk eq 1) then begin
      tvframe,dc_im[0,*,*],/bar,/aspect,btit='Dark'
      tvframe,gain_table,/bar,/aspect,btit='Gain Table'
      tvframe,flat[0,0,*,*]*fring[0,0,0,*,*],/bar,/aspect,btit='Flat'
   endif


endif ;dual-beam

!p.multi=0

;display an example of detector image
if (info.talk eq 1) then begin
   window,xs=500,ys=500,title='Detector frame example',/free
   imout=readfits(info.saves+info.files(progind)+'_'+strtrim(f-1,2)+'.fits*',/comp,/sil)
   tvframe,reform(imout[0,0,*,*])/mean(imout[0,0,*,*]),/asp,/bar,btit='I!DC!N'
endif 


end
   
