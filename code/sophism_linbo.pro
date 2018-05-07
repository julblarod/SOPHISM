PRO sophism_linbo

; ==================================================================
; DESCRIPTION
;   Simulates spectral PSF of PHI with 1 LiNbO3 etalon. Program
;   includes telecentric angle variation, dispersion, voltage
;   tuning, temperature tuning and mean tilt angle over etalon
;   tuning.
;   Depending on Pupil Apodisation module, it creates 2D arrays of
;   incident angles, etc, or makes an average of the effect if
;   it's not enabled.
;
; CALLED FROM 
;   sophism
;
; SUBROUTINES
;   sophism_linbo_lambint
;   sophism_linbo_geometry
;   headfits
;
; MAIN INPUTS
;   Hardwired = tuning constants (from TCE117)
;   temp = Etalon (oven) temperature in celsius
;   volt/wavearr = arrays of voltages or wavelength positions for
;     etalon tuning
;   frat = F-number
;   refl = Reflectivity of the etalon
;   fabfin = Fabrication finesse (over the incident area, not the
;     average over the whole etalon)
;   sizd = Width of the etalon. Constant now; 2D variable when available.
;   tilt = Inclination of etalon with respect to 
;   scvel = Velocity of S/C (for doppler shift)
;   fwhm = Prefilter FWHM
;   nc = Number of cavities of prefilter
;
; OUTPUT
;   Outputs are wavelengths (ll) and PSF (fff)
;
; VERSION HISTORY
;   vmp@iac.es. December 2010.
;   J. Blanco. Sept 2011. v1_0. Added parameters for etalon tuning by
;      wavelength input instead of voltage. Integration in SOPHISM code 
;   J. Blanco. May 2012. v2_0. Pupil apodisation part and related subroutines 
;   J. Blanco. Aug 2012. v2_1. Corrections for pupil apodisation
;   J. Blanco. Nov 2012. v2_2. Tuning jittering included. Neighbouring
;      wavelengths when no sampled data exists at selected wavelengths
;   J. Blanco. Dec 2012. v2_3. Printing and plotting for talk mode
;   J. Blanco. Dec 2012. v2_4. Correction of average calc. for fffav
;   J. Blanco. Jan 2013. v2_5. Corrected bug for number of files (to
;      not take into account already-discarded frames). Added
;      possibility of loading prefilter curve from idl-save
;      file. Still incomplete (depends on how curve is provided)
;   J. Blanco. May 2014. v2_6. Changed to save the llo and fwhm_real
;      after 'interpolating' neighbouring positions, not before as
;      previously. Brought calculation from papo for pupil apod
;      mode. Minor changes to plots.
;   J. Blanco. Jul 2014. v2_7. Corrected ntim according to new cycrep
;      scheme. Modified non-pupil apodization loop over data for
;      discardings when polarization module has not been used.
;   J. Blanco. Aug 2014. v2_8. Corrected ntim in case of not
;      polarization module for the discards. Corrected progress
;      on-screen display for dual-beam.
;   H. Waller. Feb 2017. v3_0. Included prefilter variation in several
;      modifications 
;   J. Blanco. Mar 2017. v3_1. Corrected prefilter file reading part,
;      now considering ASCII input file
;
; ==================================================================

;program index in SOPHISM code
progind=3

print,''
print,'sophism_linbo'
print,''

;loading settings file
restore,'../settings/settings_sophism.sav'

;ntim=info.modnbuff*info.nstate*info.cycrep*info.etalpos
ntim=info.modnbuff*info.nstate*info.etalpos
if (info.dataser eq 0) then ntim=1

;renaming variables from settings file. Mostly useless, should be removed
sizd=info.spectsizd

refl=info.refl
fabfin=info.fabfin
volt=info.voltarr
frat=info.frat
temp=info.etal_temp
tilt=info.tilt
scvel=info.scvel
fwhm=info.pfhwb
nc=info.ncav1
plotfig=0
;voltage mode (wavevolt=0) or wavelength mode (wavevolt=1)
if (info.wavevolt eq 1) then wavearr=info.wavearr else wavearr=info.voltarr

; lambda Fe, corrected with doppler
llo_f=info.wl
dop=scvel/299792.5d0
ldop=llo_f*dop
llo_f=llo_f+ldop

;use same wavelength array as data, no interpolation of the etalon curve
;needed. Already doppler-shifted
lll=info.lll

;define variables like usual for not-FOV prefilter case
lll_f=lll
lll_f_dop = lll_f-lll_f * dop
;FOV wavelength shift of pre-filter
if (info.pffov eq 1) then begin
;extend the prefilter wavelength domain so the shifting later does not
;create 'border' problems
;extend it with a bit more than the maximum shift
   lllext=2*ceil(abs(info.pffovmax)/info.wavesamp)+2
   lll_f=dindgen(info.wavepoint+lllext)*info.wavesamp*1d-3
   lll_f=lll_f+info.wl-lll_f(info.wavepoint+lllext-1)/2.
   lll_f=lll_f+lll_f*(info.scvel/299792.5d0)
   lll_f_dop = lll_f-lll_f * dop

;;for x_fil = 0, info.sz-1 do begin
;;   for y_fil = 0, info.sz-1 do begin
;;filter = 1./(1.+((lll_dop-info.wl)/(fwhm/2.))^(2.*nc))
;;filter = interpol(filter, lll_dop, lll)
;filter = 1./(1.+((lll_f_dop-info.wl)/(fwhm/2.))^(2.*nc))
;filter = interpol(filter, lll_f_dop, lll_f)
;;      filter[*, x_fil, y_fil] = filterCur; nc cavity prefilter
;;   endfor
;;endfor
endif 

;prefilter options
; prefilter (careful if lll is shifted by doppler)
;lll_dop=lll-lll*dop
;filter=1./(1.+((lll_dop-info.wl)/(fwhm/2.))^(2.*nc)) ; nc cavity prefilter
filter = 1./(1.+((lll_f_dop-info.wl)/(fwhm/2.))^(2.*nc))
;loading idl-save file with prefilter curve. But what will be in
;there?
;actually, ascii file with prefilter curve and lambdas
if (info.pfload eq 1) then begin 
   filterpre=read_ascii(info.pfcurve) ;restore,info.pfcurve
   lll_dopcur=reform(filterpre.field1[0,*])
   filtercur=reform(filterpre.field1[1,*])
;if the wavelength range of the loaded prefilter is smaller than of
;the data, the interpolation will send it to negative values
;make use of starting and ending points of lorentzian prefilter and
;interpolate then
   nlll=n_elements(lll_f_dop)
;   if (min(lll_dopcur) gt min(lll_dop)) then begin
   if (min(lll_dopcur) gt min(lll_f_dop)) then begin
      lll_dopcur=[lll_f_dop[0:1],lll_dopcur]
      filtercur=[filter[0:1],filtercur]
   endif
;   if (max(lll_dopcur) lt max(lll_dop)) then begin
   if (max(lll_dopcur) lt max(lll_f_dop)) then begin
      lll_dopcur=[lll_dopcur,lll_f_dop[nlll-2:nlll-1]]
      filtercur=[filtercur,filter[nlll-2:nlll-1]]
   endif 
   filter=filtercur
   lll_dop=lll_dopcur
   lll_f_dop=lll_dop
endif 
;pero deberia ponerlo correctamente, porque luego multiplico uno por
;otro...
;dejo el interpol o mejor 'corto' sin mas, para no meter cosas raras?
filter=interpol(filter,lll_f_dop,lll_f)
;the interpolation could create negative values if the range is
;smaller originally than lll and the curve was falling
;Prevent those nonsense negatives
wfn=where(filter lt 0,nwfn)
if (nwfn gt 0) then filter[wfn]=0.


; first compute etalon at 35 degress and 0 volt
; refractive index with dispersion from CSIRO.
temper=35.+273.15d0

n01=4.913d0
n02=0.1173d0+1.65e-8*(temper^2)
n03=(lll/1.d4)^2-(0.212+2.7d-8*(temper^2))^2
n04=2.78d-2*(lll/1.d4)^2
nn_dis=sqrt(n01+n02/n03-n04)

; uncomment this line to see the effect of dispersion
;nn_dis[*]=mean(nn_dis)

dd0=sizd ;
dd0=dd0[0]

opd_dis=nn_dis*dd0

; temperature & Voltage sensitivity
cll=[5250.d0,6328.d0] ; CSIRO's manual
;cvv=[2.82d-4,3.44d-4]; Angstroms/Volt
;ctt=[0.0252,0.0293] ; Angstrom/degree
cvv=[2.82d-5,3.44d-5]; mu/Volt ;nm/Volt realmente
ctt=[2.52d-3,2.93d-3] ; mu/degree ;nm/Volt realmente

a=interpol(cvv,cll,llo_f)
if (info.talk eq 1) then print,'OPD voltage tuning coefficient is (microns/Volt)=',a
b=interpol(ctt,cll,llo_f)
if (info.talk eq 1) then print,'OPD thermal tuning coefficient is (microns/Degree)=',b

; find wavelength transmitted at maximum
; or voltage for a given wavelength

case (info.wavevolt) of

   0: begin
;start case of voltage mode
      nwave=info.etalpos
;define arrays for number of voltage positions
      llo=dblarr(nwave)
      mm0=intarr(nwave)
      opd_disarr=dblarr(n_elements(lll),nwave)
      opd0=dblarr(nwave)
;add 'tuning jittering'
      if (info.voltjit eq 1) then begin
         voltrand=randomn(seed,nwave)
         if (nwave eq 1) then voltrand=randomn(seed,6)
         voltrand=voltrand-mean(voltrand)
         voltrand=voltrand/stddev(voltrand)*info.voltjitrms
         if (nwave eq 1) then voltrand=voltrand(0)
         volt=volt+voltrand
      endif 

      for ww=0,nwave-1 do begin
; tune wavelength for T and Volt (for each voltage)
         opd_disarr[*,ww]=nn_dis*dd0+(temp-35.)*b+volt[ww]*a
; nn at Fe lambda
         opd0[ww]=interpol(opd_disarr[*,ww],lll,llo_f)
         opd_check=opd0[ww]
; transmission order
         mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llo_f))
; transmitted lambda
         llo[ww]=1.d4*2.*opd0[ww]/float(mm0[ww]) 
; opd at transmitted lambda
         opd0[ww]=interpol(opd_disarr[*,ww],lll,llo[ww]) 
; need to iterate the referactive index (here opd) for new lambda
         while (abs(opd0[ww]-opd_check) gt 1.d-7) do begin
            opd_check=opd0[ww]
            mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llo[ww]))
; transmitted lambda
            llo[ww]=1.d4*2.*opd0[ww]/float(mm0[ww])
; opd at transmitted lambda
            opd0[ww]=interpol(opd_disarr[*,ww],lll,llo[ww]) 
         endwhile
      endfor
      opd_dis=opd_disarr
   
;this would be incorrect, since there may be no input data sampling at
;the transmitted lambda for the voltage.
;Should select the two closest wavelengths and continue through steps
;like in option 1 below.
      sophism_linbo_lambint,llo,index,wg,waves,wavesind
      
;now, obtain opd's for these wavelengths, repeating the process
;start with the calculations for the new lambdas a la 'lambda mode'
      opd_dis=nn_dis*dd0+(temp-35.)*b
      nwave=n_elements(waves)
      volt=fltarr(nwave)
      llp=dblarr(nwave)
      mm0=intarr(nwave)
      opd_disarr=dblarr(n_elements(lll),nwave)
      opd0=dblarr(nwave)
      for ww=0,nwave-1 do begin
;;careful!!! incluir tilt?
         llp[ww]=waves[ww]
         opd0[ww]=interpol(opd_dis,lll,llp[ww])
         mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llp[ww]))
         llo=1.d4*2.*opd0[ww]/float(mm0[ww])
         volt[ww]=(llp[ww]-llo)*mm0[ww]/2./1.d4/a
         opd_disarr[*,ww]=nn_dis*dd0+(temp-35.)*b+volt[ww]*a
         opd0[ww]=interpol(opd_disarr[*,ww],lll,llp[ww])
      endfor
      llo=waves
      opd_dis=opd_disarr

   end

   1: begin
   ;start case of wavelength mode

;two ways depending on tuning jittering. Without it, much simpler
      if (info.voltjit eq 0) then begin
;array of wavelength positions
         llo=llo_f+wavearr*1d-3
;check if given positions are sampled in input data or should be
;'interpolated'
         sophism_linbo_lambint,llo,index,wg,waves,wavesind
; tune wavelength for T and Volt (with Volt=0)
         opd_dis=opd_dis+(temp-35.)*b
;define arrays for number of wavelength positions
         nwave=n_elements(waves)
         volt=fltarr(nwave)
         llp=dblarr(nwave)
         mm0=intarr(nwave)
         opd_disarr=dblarr(n_elements(lll),nwave)
         opd0=dblarr(nwave)
         for ww=0,nwave-1 do begin
;careful!!! incluir tilt?
            llp[ww]=waves[ww]
            opd0[ww]=interpol(opd_dis,lll,llp[ww])
            mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llp[ww]))
            llo=1.d4*2.*opd0[ww]/float(mm0[ww])
            volt[ww]=(llp[ww]-llo)*mm0[ww]/2./1.d4/a
            opd_disarr[*,ww]=nn_dis*dd0+(temp-35.)*b+volt[ww]*a
            opd0[ww]=interpol(opd_disarr[*,ww],lll,llp[ww])
         endfor

         llo=llp
         opd_dis=opd_disarr

      endif else begin

;tuning jittering case   
; tune wavelength for T and Volt (with Volt=0)
         opd_dis=opd_dis+(temp-35.)*b
;define arrays for number of wavelength positions
         nwave=info.etalpos
         volt=fltarr(nwave)
         llp=dblarr(nwave)
         mm0=intarr(nwave)
         opd_disarr=dblarr(n_elements(lll),nwave)
         opd0=dblarr(nwave)

         for ww=0,nwave-1 do begin
;careful!!! incluir tilt?
            llp[ww]=llo_f+wavearr[ww]*1d-3
            opd0[ww]=interpol(opd_dis,lll,llp[ww])
            mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llp[ww]))
            llo=1.d4*2.*opd0[ww]/float(mm0[ww])
            volt[ww]=(llp[ww]-llo)*mm0[ww]/2./1.d4/a
            opd_disarr[*,ww]=nn_dis*dd0+(temp-35.)*b+volt[ww]*a
            opd0[ww]=interpol(opd_disarr[*,ww],lll,llp[ww])
         endfor
      
;add 'tuning jittering'
         voltrand=randomn(seed,nwave)
         if (nwave eq 1) then voltrand=randomn(seed,6)
         voltrand=voltrand-mean(voltrand)
         voltrand=voltrand/stddev(voltrand)*info.voltjitrms
         if (nwave eq 1) then voltrand=voltrand(0)
         volt=volt+voltrand
;find wavelengths and opds for the new voltages a la 'voltage mode'         
         llo=dblarr(nwave)
         mm0=intarr(nwave)
         opd0=dblarr(nwave)
         opd_disarr=dblarr(n_elements(lll),nwave)
         for ww=0,nwave-1 do begin
            opd_disarr[*,ww]=nn_dis*dd0+(temp-35.)*b+volt[ww]*a
            opd0[ww]=interpol(opd_disarr[*,ww],lll,llo_f)
            opd_check=opd0[ww]
            mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llo_f))
            llo[ww]=1.d4*2.*opd0[ww]/float(mm0[ww]) 
            opd0[ww]=interpol(opd_disarr[*,ww],lll,llo[ww]) 
            while (abs(opd0[ww]-opd_check) gt 1.d-7) do begin
               opd_check=opd0[ww]
               mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llo[ww]))
               llo[ww]=1.d4*2.*opd0[ww]/float(mm0[ww])
               opd0[ww]=interpol(opd_disarr[*,ww],lll,llo[ww]) 
            endwhile
         endfor
         opd_dis=opd_disarr
;check if wavelengths obtained are sampled in input data
        sophism_linbo_lambint,llo,index,wg,waves,wavesind

;repeat now the lambda mode calculations to obtain the opd's at
;these new wavelengths
         opd_dis=nn_dis*dd0+(temp-35.)*b
         nwave=n_elements(waves)
         volt=fltarr(nwave)
         llp=dblarr(nwave)
         mm0=intarr(nwave)
         opd_disarr=dblarr(n_elements(lll),nwave)
         opd0=dblarr(nwave)
         for ww=0,nwave-1 do begin
;;careful!!! incluir tilt?
            llp[ww]=waves[ww]
            opd0[ww]=interpol(opd_dis,lll,llp[ww])
            mm0[ww]=fix(round(2.*opd0[ww]*1.d4/llp[ww]))
            llo=1.d4*2.*opd0[ww]/float(mm0[ww])
            volt[ww]=(llp[ww]-llo)*mm0[ww]/2./1.d4/a
            opd_disarr[*,ww]=nn_dis*dd0+(temp-35.)*b+volt[ww]*a
            opd0[ww]=interpol(opd_disarr[*,ww],lll,llp[ww])
         endfor
         opd_dis=opd_disarr
         llo=waves
      endelse 
   end 

endcase

nn0=opd0/dd0 ; this is not perfect as it does
	     ; not use the relation dd0(V,T). It is
	     ; only used for angle and FSR
             ; probably OK

; compute some etalon parameters
; reflectivity
rr0=refl 
; reflecttivity finesse
ffr=!pi*sqrt(rr0)/(1.-rr0)

; telecentric configuration
alpha=1./(2.d0*frat)
;etendue finesse 
ffe=(llo*1.d-4)*nn0/(dd0*(sin(alpha))^2) 

; total finesse
ff0_par=1./sqrt((1./ffr)^2+(1./fabfin)^2) ; etendue finesse should not be included here ; it is included later 
cte0=2.d0*ff0_par/!pi ; for later use
; total finesse including etendue
ff0=1./sqrt((1./ffr)^2+(1./fabfin)^2+(1./ffe)^2) 
 
if (info.talk eq 1) then begin
   print,'Applied Voltage: ',volt
   print,'Order: ',mm0
   print,'Transmitted lambda (A): ',llo
   print,'Normal incidence OPD (microns): ',opd0
   print,'Refractive index for nominal dd0: ',nn0
   print,'Maximum Incidence Angle (deg): ',alpha*180./!pi
   print,'Reflective Finesse: ', ffr
   print,'Telecentric Finesse: ',ffe
   print,'Fabrication Finesse: ', fabfin
   print,'Total Finesse: ',ff0 
endif

;*************************************************************************
; here first, the calculations for pupil apodization mode

if (info.routines(5) eq 1) then begin 

;calculate the pupil radius for every wavelength
;necessary? The difference is not much (0.02 pix for info.sz=64) and it
;complicates develop and computer time
   deltanu=1./(info.sz*info.sssamp)
;cutfreq=(info.telap/(lll*1.e-10))*(1./(!radeg*3600.))
   cutfreq=(info.telap/(llo_f*1.e-10))*(1./(!radeg*3600.))
   rpup=cutfreq/(2.*deltanu)
   tam=check_rpupil(info.sz/2,rpup)
   rd=radius_aper2(tam/2.,rpup,maskpup)

;define final variables to store results for all wavel. points
   fsr0=dblarr(nwave)
   fff=dblarr(tam,tam,info.wavepoint,nwave)
   fwhm_real=dblarr(nwave)
   phase_et=dblarr(tam,tam,info.wavepoint,nwave)

; dispersion 
   for ww=0,nwave-1 do begin
      delta0=0.
      fp0_arr=0.
      delta=0.

      disper=deriv(lll/1.d4,opd_dis[*,ww]/dd0)*(lll/1.d4)
      fac_disper=interpol(disper,lll,llo[ww])
      fsr0[ww]=llo[ww]^2*1e-4/(2.*(nn0[ww]-fac_disper)*dd0)
      fwhm0=fsr0[ww]/ff0_par

;calculate 2D array of incident angles
; introducing opd_dis/dd0 for the internal angles (snellius)
; calculation. Not completely correct but enough?
      sophism_linbo_geometry,tam,rpup,llo[ww],opd_dis[*,ww]/dd0,anglesarr_int,intmask

      llo_tilt=(1.d4*2.*opd0[ww]/float(mm0[ww]))*cos(anglesarr_int)*intmask
      ; opd0[ww]=interpol(opd_dis[*,ww],lll,llo_tilt) ; nn at transmited lambda
      ; nn0[ww]=opd0[ww]/dd0

      delta0=2.d0*!pi*opd_dis[0,ww]*cos(anglesarr_int[*,*,0])/(lll[0]/1.d4)
      for jj=1,n_elements(lll)-1 do begin
         delta=2.d0*!pi*opd_dis[jj,ww]*cos(anglesarr_int[*,*,jj])/(lll[jj]/1.d4)
         delta0=[[[delta0]],[[delta]]]
      endfor

      fp0_arr=1./(1+(cte0^2)*((sin(delta0))^2))
      for jj=0,n_elements(lll)-1 do fp0_arr[*,*,jj]=fp0_arr[*,*,jj]*intmask

;apply prefilter to etalon transmission
      for xx=0,tam-1 do begin 
         for yy=0,tam-1 do begin
            fffpre=reform(fp0_arr[xx,yy,*])*filter
            fff[xx,yy,*,ww]=fffpre
         endfor 
      endfor 
;mejor no interpolar sino tomar los correspondientes, ahora que ya
;existen, correcto?

;obtain averaged transmission for the whole pupil, used in later routines
      fffav=total(total(fff,1),1)/float(total(maskpup))

      tmparr=abs(lll-llo[ww])
      nn=where(tmparr eq min(tmparr))
      nn=nn(0)
      blue=interpol(lll(0:nn),fffav(0:nn,ww),max(fffav[*,ww])/2.)
      blue=blue(0)
      red=interpol(lll(nn:*),fffav(nn:*,ww),max(fffav[*,ww])/2.)
      red=red(0)
      fwhm_real[ww]=double(red)-blue
      if (info.talk eq 1) then begin 
         print,'FSR (A)=', fsr0[ww]
         print,'Estimated intrinsic (no telecentric cone) FWHM (A)=', fwhm0
         print,'Real FWHM (A)=', fwhm_real[ww]
      endif 

;transmitted phases
      phase_et[*,*,*,ww]=atan((refl*sin(2*delta0))/(1-refl*cos(2*delta0)))

   endfor ;ww

;interpolate the wavelengths in case of neighourin to get the selected
;positions
   llopre=llo
   llo=fltarr(info.etalpos)
   fwhmpre=fwhm_real
   fwhm_real=fltarr(info.etalpos)
   mind=max(index)
   for ii=0,mind do begin
      indy=where(index eq ii,nindy)
      for bb=0,nindy-1 do begin
         llo[ii]=llo[ii]+double(wg[indy(bb)]*llopre[indy(bb)])
         fwhm_real[ii]=fwhm_real[ii]+double(wg[indy(bb)]*fwhmpre[indy(bb)])
      endfor
   endfor

;save lots of stuff
   save,opd_dis,lll,llo,wg,index,wavesind,fff,fffav,filter,nn0,fsr0,ff0,volt,fwhm_real, filename=info.saves+info.files(progind)+'.sav'

; save the phases in a different file
   save,phase_et,filename=info.saves+info.files(progind)+'_phase.sav'

;display examples of pupil intensities and plot etalon curve
if (info.talk eq 1) then begin
   fwhmpos=round(fwhm_real(0)/(info.wavesamp*1d-3))
   pupint=fltarr(tam*5,tam)
   pupint[2*tam:3*tam-1,*]=fff[*,*,wavesind(0),0]
   pupint[0:tam-1,*]=fff[*,*,wavesind(0)-fwhmpos,0]
   pupint[tam:2*tam-1,*]=fff[*,*,wavesind(0)-fwhmpos/2.,0]
   pupint[3*tam:4*tam-1,*]=fff[*,*,wavesind(0)+fwhmpos/2.,0]
   pupint[4*tam:5*tam-1,*]=fff[*,*,wavesind(0)+fwhmpos,0]

   window,title='Pupil intensities at -fwhm, -fwhm/2, 0, fhwm/2, fwhm',xsiz=800,ysiz=200,/free
   loadct,39
   tvframe,pupint/max(pupint),/asp,/bar
   loadct,0

   restore,'../data/fts6173.save'
   iii=iii*1.d-4
   iii=fft_shift(iii,ldop/mean(deriv(xxx)))
   iip=interpol(iii,xxx,lll)
   window,title='Etalon+Prefilter Transmission (aver. over pupil)',/free
   plot,lll,fffav[*,0],xrange=[fix(llo_f-4.),fix(llo_f+4.)],yrange=[0,1],xtitle=xtitle, ytitle=ytitle
   oplot,lll,iip
   oplot,lll,filter,lines=1
   plots,[llo(0),llo(0)],[0,1],lines=2
   plots,[llo(0)-fsr0(0),llo(0)-fsr0(0)],[0,1],lines=2
   plots,[llo(0)+fsr0(0),llo(0)+fsr0(0)],[0,1],lines=2

   if (info.etalpos gt 1) then begin
      oplot,lll,fffav[*,nwave-1],lines=3
      plots,[llopre(nwave-1),llopre(nwave-1)],[0,1],lines=4
      plots,[llopre(nwave-1)-fsr0(nwave-1),llopre(nwave-1)-fsr0(nwave-1)],[0,1],lines=4
      plots,[llopre(nwave-1)+fsr0(nwave-1),llopre(nwave-1)+fsr0(nwave-1)],[0,1],lines=4
   endif

endif ;talk


;***********************************************************************
; below, the classical way, no pupil-apod

endif else begin

; total angles to be computed
   nn_tilt=81                
   nns=(nn_tilt-1)/2.d0

;define final variables to store results for all wavel. points
   fsr0=dblarr(nwave)
; transmitted lambdas for each incidence angle
   llo_arr=dblarr(nn_tilt,nwave) 
   fp0=dblarr(n_elements(lll),nwave)
;   fff=dblarr(n_elements(lll),nwave)
   if (info.pffov eq 1) then fff=dblarr(info.sz,info.sz,n_elements(lll),nwave) else fff=dblarr(n_elements(lll),nwave)
   fwhm_real=dblarr(nwave)
   telcon=dblarr(nwave)
   phase_et=dblarr(n_elements(lll),nwave)

   for ww=0,nwave-1 do begin

      disper=deriv(lll/1.d4,opd_dis[*,ww]/dd0)*(lll/1.d4)
      fac_disper=interpol(disper,lll,llo[ww])
      fsr0[ww]=llo[ww]^2*1e-4/(2.*(nn0[ww]-fac_disper)*dd0)
      fwhm0=fsr0[ww]/ff0_par

      ; loop over telecentric incidence cone angles 

      angle_arr=(alpha*(findgen(nn_tilt)-nns)/nns)*180.d0/!pi+tilt ;degress
      fp0_arr=dblarr(n_elements(lll),nn_tilt)
; angles inside etalon
      theta0_arr=dblarr(nn_tilt)
; refractive indexes 
      nn0_arr=dblarr(nn_tilt)

      nn0_ref=nn0[ww]
      llo_ref=llo[ww]

      for ind=0,nn_tilt-1 do begin	
         nn_check=nn0_ref
; angle inside etalon
         theta0=180.d0*asin(sin(double(angle_arr[ind])*!pi/180.d0)/nn0[ww])/!pi 
         llo_tilt=(1.d4*2.*opd0[ww]/float(mm0[ww]))*cos(theta0*!pi/180.d0)
         ; llo_tilt=llo_ref*cos(theta0*!pi/180.d0) ; this would be wrong
; nn at transmited lambda
         opd0[ww]=interpol(opd_dis[*,ww],lll,llo_tilt)
         nn0[ww]=opd0[ww]/dd0
         
         while (abs(nn0[ww]-nn_check) gt 1.d-7) do begin
            nn_check=nn0[ww]
; angle inside etalon 
            theta0=180.d0*asin(sin(double(angle_arr[ind])*!pi/180.d0)/nn0[ww])/!pi
            llo_tilt=(1.d4*2.*opd0[ww]/float(mm0[ww]))*cos(theta0*!pi/180.d0)
            ; llo_tilt=llo_ref*cos(theta0*!pi/180.d0) ; this would be wrong
; nn at transmited lambda
            opd0[ww]=interpol(opd_dis[*,ww],lll,llo_tilt)
            nn0[ww]=opd0[ww]/dd0
         endwhile
         
         if (ind eq fix(nns) AND info.talk eq 1) then begin
            print,'Angle inside etalon (deg)=',theta0
            print,'Transmitted lambda at this tilt (A)=',llo_tilt
         endif
         
         theta0_arr[ind]=theta0
         llo_arr[ind,ww]=llo_tilt
         nn0_arr[ind]=nn0[ww]
         
         delta0=2.d0*!pi*opd_dis[*,ww]*cos(theta0*!pi/180.)/(lll/1.d4)
         fp0_arr[*,ind]=1./(1+(cte0^2)*((sin(delta0))^2))

      endfor

;telecentric cone
      telcon[ww]=max(llo_arr[*,ww])-min(llo_arr[*,ww]) 

; The reason why high inclination get so broad and
; assymetric is a combination of greater separations for
; higher tilts and the non-linear effect

      fp0[*,ww]=total(fp0_arr,2)/float(nn_tilt)
      dd=where(lll gt min(llo_arr[*,ww]) and lll lt max(llo_arr[*,ww]))
      if (dd[0] ne -1) then begin
         intmax=max(fp0[dd,ww])
      endif else begin
         dd=where(abs(lll-(min(llo_arr[*,ww])*.5+max(llo_arr[*,ww])*.5)) eq $
                  min(abs(lll-(min(llo_arr[*,ww])*.5+max(llo_arr[*,ww])*.5))))
         intmax=max(fp0[dd,ww])
      endelse
      llo2=lll[dd[where(fp0[dd,ww] eq intmax)]]
      llo[ww]=llo2[0]
     
; etalon times filter 

;hhhhhhhhhhhhhhhhhhhhhh
      if (info.pffov eq 1) then begin
         para_fff = fltarr(info.sz, info.sz)
         for xx = 0, info.sz-1 do begin 
            for yy = 0, info.sz-1 do begin 
               para_fff[xx, yy] = (xx-info.sz/2)^2.+(yy-info.sz/2)^2.
            endfor
         endfor
         para_fff = (para_fff / max(para_fff)) * info.pffovmax 
;jjjjjjjjjjjj
;these 3 next lines just to find the 'unperturbed' prefilter position
;and calculate fwhm there afterwards
         wwfilt=where(para_fff eq 0,nwwfilt)
         if (nwwfilt eq -1) then wwfilt=where(abs(para_fff) eq min(abs(para_fff)))
         wwfilt=wwfilt[0]
         ywwfilt=wwfilt/info.sz
         xwwfilt=wwfilt-(ywwfilt*info.sz)
;         filterpre=filter
         filterpara = make_array(info.wavepoint, info.sz, info.sz)
;         filter = make_array(info.wavepoint, info.sz, info.sz)
;jjjjjjjjjjjj
;      largefilter = make_array(info.wavepoint+110, value=min(filter))
;      largefilter[54:53+info.wavepoint] = filter
;Just FYI filter does not need to be 3D can be 1D
         for x_fff = 0, info.sz-1 do begin
            for y_fff = 0, info.sz-1 do begin
;               filtshift = fft_shift(filter, para_fff[x_fff, y_fff] / info.wavesamp)	
;               filterpara[*, x_fff, y_fff] = filtshift[lllext/2:lllext/2+info.wavepoint-1]
;               fff[x_fff, y_fff, *, ww] = fp0[*, ww] * filterpara[*, x_fff, y_fff] 
;probably could do this in just one go without the extra filtshift
;variable. Maybe useful for saving memory?
               filtshift = fft_shift(filter, para_fff[x_fff, y_fff] / info.wavesamp)	
               filterpara[*, x_fff, y_fff] = filtshift[lllext/2:lllext/2+info.wavepoint-1]
               fff[x_fff, y_fff, *, ww] = fp0[*, ww] * filterpara[*, x_fff, y_fff] 
;               fff[x_fff, y_fff, *, ww] = fp0[*, ww] * filter;uniform etaloncase
            endfor
         endfor
      endif else begin    
         fff[*,ww]=fp0[*,ww]*filter
      endelse
;hhhhhhhhhhhhhhhhhhhhhh

      tmparr=abs(lll-llo[ww])

      nn=where(tmparr eq min(tmparr))
      nn=nn(0)
;jjjjjjjjjjjjjjjjjjjjjjj
      if (info.pffov eq 1) then begin
         blue=interpol(lll(0:nn),fff[xwwfilt,ywwfilt,0:nn,ww],max(fff[xwwfilt,ywwfilt,*,ww])/2.)
         blue=blue(0)
         red=interpol(lll(nn:*),fff[xwwfilt,ywwfilt,nn:*,ww],max(fff[xwwfilt,ywwfilt,*,ww])/2.)
         red=red(0)
      endif else begin 
         blue=interpol(lll(0:nn),fff(0:nn,ww),max(fff[*,ww])/2.)
         blue=blue(0)
         red=interpol(lll(nn:*),fff(nn:*,ww),max(fff[*,ww])/2.)
         red=red(0)
      endelse 
      fwhm_real[ww]=double(red)-blue
;jjjjjjjjjjjjjjjjjjjjjjj

   endfor ;ww nwave

;print some info, plot the etalon transmission, filter, and atlas
;spectral curve fts I
   if (info.talk eq 1) then begin

      print,'FSR (A): ', fsr0
      print,'Estimated intrinsic (no telecentric cone) FWHM (A): ', fwhm0
      print,'Maximum transmitted intensity (prefilter not included): ',intmax
      print,'Effective transmitted lambda (A): ',llo
      print,'Telecentric blueshift (A): ',llo_ref-llo
      print,'Telecentric cone total lambda shift (A): ',telcon
      print,'Real FWHM (A): ', fwhm_real

   ;restoring profile from FTS atlas
      restore,'../data/fts6173.save'
      iii=iii*1.d-4
      
      print,'Distance to Fe I line (A): ',llo-llo_f   
      iii=fft_shift(iii,ldop/mean(deriv(xxx)))
      iip=interpol(iii,xxx,lll)
      window,title='Etalon+Prefilter Transmission',/free
;por ahora hacer un case para elegir entre FOV prefilter o no (o con
;roughness del etalon)
;jjjjjjjjjjjjjjjjjjjjjjj
      case info.pffov of
         0: begin
            plot,lll,fff[*,0],xrange=[fix(llo_f-4.),fix(llo_f+4.)],yrange=[0,1],$
                 title=string(fix(volt[0]))+' Volts    '+ $ 
                 strmid(strcompress(string(temp),/remove_all),0,6)+' !uo!nC      '+ $
                 strmid(strcompress(string(tilt),/remove_all),0,5)+' degrees     '+ $
                 strmid(strcompress(string(scvel),/remove_all),0,3)+' km/s  ' 
            oplot,lll,iip
            plots,[llo_f,llo_f],[0,1],lines=1
            plots,[llo(0),llo(0)],[0,1],lines=1
            plots,[llo(0)-fsr0(0),llo(0)-fsr0(0)],[0,1],lines=2
            plots,[llo(0)+fsr0(0),llo(0)+fsr0(0)],[0,1],lines=2
            oplot,lll,filter,lines=1
            if (info.etalpos gt 1) then begin
               oplot,lll,fff[*,nwave-1],lines=3
               plots,[llo(nwave-1),llo(nwave-1)],[0,1],lines=3
               plots,[llo(nwave-1)-fsr0(nwave-1),llo(nwave-1)-fsr0(nwave-1)],[0,1],lines=3
               plots,[llo(nwave-1)+fsr0(nwave-1),llo(nwave-1)+fsr0(nwave-1)],[0,1],lines=3
            endif

            for ww=0,nwave-1 do begin
               print,'Secondary Blue side maximum=',interpol(fff[*,ww],lll,llo[ww]-fsr0[ww])
               print,'Secondary Red side maximum=',interpol(fff[*,ww],lll,llo[ww]+fsr0[ww])

; stray-light
               ddb=where(lll ge llo[ww]-fwhm_real[ww]*2 and lll le llo[ww]+fwhm_real[ww]*2)
               ddl=where(lll lt llo[ww]-fwhm_real[ww]*2)         
               ddg=where(lll gt llo[ww]+fwhm_real[ww]*2)        
;jjjjjjjjjjjjjjjj
;correction in case no indexes are found above and/or below
               if (min(ddl) lt 0) then ddl=0
               if (min(ddg) lt 0) then ddg=n_elements(lll)-1
;jjjjjjjjjjjjjjjj
               print,'Spectral Stray-light (%)=',100*total([fff(ddl,ww),fff(ddg,ww)])/total(fff(ddb,ww))
               print,'Spectral Stray-light weighted with I (%)=',100*total([fff(ddl,ww)*iip(ddl),fff(ddg,ww)*iip(ddg)])/total(fff(ddb,ww)*iip(ddb))
            endfor ;ww nwave (II)
            end  
         
         1: begin
            loadct,39
            ;first in white, the 'unperturbed' case
            plot,lll,fff[xwwfilt,ywwfilt,*,0],xrange=[fix(llo_f-4.),fix(llo_f+4.)],yrange=[0,1],$
                 title=string(fix(volt[0]))+' Volts    '+ $ 
                 strmid(strcompress(string(temp),/remove_all),0,6)+' !uo!nC      '+ $
                 strmid(strcompress(string(tilt),/remove_all),0,5)+' degrees     '+ $
                 strmid(strcompress(string(scvel),/remove_all),0,3)+' km/s  ' 
            oplot,lll,iip
            plots,[llo_f,llo_f],[0,1],lines=1
            plots,[llo(0),llo(0)],[0,1],lines=1
            plots,[llo(0)-fsr0(0),llo(0)-fsr0(0)],[0,1],lines=2
            plots,[llo(0)+fsr0(0),llo(0)+fsr0(0)],[0,1],lines=2
            oplot,lll,filterpara[*,xwwfilt,ywwfilt],lines=1
            oplot,lll,fff[0,0,*,0],col=250,thick=2.5
            oplot,lll,filterpara[*,0,0],lines=1,col=250,thick=2.5
            if (info.etalpos gt 1) then begin
               oplot,lll,fff[xwwfilt,ywwfilt,*,nwave-1],lines=3
               plots,[llo(nwave-1),llo(nwave-1)],[0,1],lines=3
               plots,[llo(nwave-1)-fsr0(nwave-1),llo(nwave-1)-fsr0(nwave-1)],[0,1],lines=3
               plots,[llo(nwave-1)+fsr0(nwave-1),llo(nwave-1)+fsr0(nwave-1)],[0,1],lines=3
            endif

            for ww=0,nwave-1 do begin
               print,'Unperturbed Secondary Blue side maximum=',interpol(fff[xwwfilt,ywwfilt,*,ww],lll,llo[ww]-fsr0[ww])
               print,'Unperturbed Secondary Red side maximum=',interpol(fff[xwwfilt,ywwfilt,*,ww],lll,llo[ww]+fsr0[ww])
               print,'Secondary Blue side maximum=',interpol(fff[0,0,*,ww],lll,llo[ww]-fsr0[ww])
               print,'Secondary Red side maximum=',interpol(fff[0,0,*,ww],lll,llo[ww]+fsr0[ww])

; stray-light
               ddb=where(lll ge llo[ww]-fwhm_real[ww]*2 and lll le llo[ww]+fwhm_real[ww]*2)
               ddl=where(lll lt llo[ww]-fwhm_real[ww]*2)         
               ddg=where(lll gt llo[ww]+fwhm_real[ww]*2)        
;jjjjjjjjjjjjjjjj
;correction in case no indexes are found above and/or below
               if (min(ddl) lt 0) then ddl=0
               if (min(ddg) lt 0) then ddg=n_elements(lll)-1
;jjjjjjjjjjjjjjjj
               print,'Unperturbed Spectral Stray-light (%)=',100*total([reform(fff(xwwfilt,ywwfilt,ddl,ww)),reform(fff(xwwfilt,ywwfilt,ddg,ww))])/total(fff(xwwfilt,ywwfilt,ddb,ww))
               print,'Unperturbed Spectral Stray-light weighted with I (%)=',100*total([reform(fff(xwwfilt,ywwfilt,ddl,ww))*iip(ddl),reform(fff(xwwfilt,ywwfilt,ddg,ww))*iip(ddg)])/total(fff(xwwfilt,ywwfilt,ddb,ww)*iip(ddb))
               print,'Spectral Stray-light (%)=',100*total([reform(fff(0,0,ddl,ww)),reform(fff(0,0,ddg,ww))])/total(fff(0,0,ddb,ww))
               print,'Spectral Stray-light weighted with I (%)=',100*total([reform(fff(0,0,ddl,ww))*iip(ddl),reform(fff(0,0,ddg,ww))*iip(ddg)])/total(fff(0,0,ddb,ww)*iip(ddb))
            endfor ;ww nwave (II)
            loadct,0
            end 
      endcase 
;jjjjjjjjjjjjjjjjjjjjjj
   endif ;talk

; saving in a .sav the prefilter & fpi curve
;      save,wg,index,wavesind,llo,opd_dis,lll,fff,filter,nn0,fp0,fsr0,ff0,fwhm_real,volt, filename=info.saves+info.files(progind)+'.sav'

   if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
   if (progma eq -1) then progma=0
   head=headfits(info.saves+info.files(progma)+'_0.fits*')
   headmod=sxpar(head,'HISTORY')
   poli=where(headmod eq 'Polarization Module')
;in case the polarization module hasn't been run, the ntim may
;have to include samples for etalon discard, so use the ntim
;calculated when starting simulation
   if (poli eq -1) then ntim=info.ntim

;'interpolate' neighbouring info to real selected positions, if it was
;not sampled in input data 
   llopre=llo
   llo=fltarr(info.etalpos)
   fwhmpre=fwhm_real
   fwhm_real=fltarr(info.etalpos)
   mind=max(index)
   for ii=0,mind do begin
      indy=where(index eq ii,nindy)
      for bb=0,nindy-1 do begin
         llo[ii]=llo[ii]+double(wg[indy(bb)]*llopre[indy(bb)])
         fwhm_real[ii]=fwhm_real[ii]+double(wg[indy(bb)]*fwhmpre[indy(bb)])
      endfor 
   endfor

; saving in a .sav the prefilter & fpi curve
      save,wg,index,wavesind,llo,opd_dis,lll,fff,filter,nn0,fp0,fsr0,ff0,fwhm_real,volt, filename=info.saves+info.files(progind)+'.sav'

   if (info.talk eq 1) then print,'Data-etalon convolution'

prg0=0.
prg=prg0
nst=0.
lamind=0
nout=0
   for nf=0,ntim-1 do begin
;in case the discaretal is selected but polarization module has not
;been run, discard here the appropriate data
;for dual-beam should not be done because only makes sense if
;polarizing, right?
;a 0 in the first index of ndeaths because there is no polarization
      if (poli eq -1 AND nf eq nst+info.nstate+info.ndeaths(0,lamind)) then begin
         nst=nst+info.nstate+info.ndeaths(0,lamind)
         nout=nout+info.ndeaths(0,lamind)
         lamind=lamind+1
      endif 

;if polarization has been run, all remaining frames are considered
      if ((poli eq -1 AND nf lt nst+info.nstate) OR (poli ne -1)) then begin
         ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/compr,/sil)
         sizima=size(ima)
         imafpi=fltarr(nwave,sizima(2),sizima(3),sizima(4))
         for ww=0,nwave-1 do begin
;normalization to 'conserve' the input level of intensity
;correct? Wouldn't we be integrating all the entering light? We
;should have higher level of intensity, no?
;jjjjjjjjjjjjjjjjjjj
;in case of FOV prefilter, need to calculate normalization for each
;pixel since the contribution will/may be different 
            if (info.pffov eq 1) then norm=total(fff[*,*,*,ww],3) else norm=total(fff[*,ww],1)
;jjjjjjjjjjjjjjjjjjj
            imafpipre=fltarr(sizima(2),sizima(3),sizima(4))
;hhhhhhhhhhhhh
;jjj change dimension order of fff to be similar as ima? May speed
;some things up or less if-loops. Think about it
            if (info.pffov eq 1) then begin
               for sto = 0, sizima(2)-1 do begin
                  for step = 0, info.wavepoint-1 do begin      
                     imafpipre[sto, *, *] = imafpipre[sto, *, *] + reform(ima[step, sto, *, *])*fff[*, *, step, ww] ;changed 14/12/16
                  endfor
                  imafpipre[sto,*,*]=reform(imafpipre[sto,*,*])/norm
               endfor
               imafpi[ww,*,*,*]=imafpipre   
            endif else begin
               for step=0,info.wavepoint-1 do imafpipre = imafpipre + reform(ima[step,*,*,*])*fff[step,ww]
               imafpi[ww,*,*,*]=imafpipre/norm
            endelse 
;hhhhhhhhhhhhhh
            prg=round(100.*float(ww+1)*float(nf+1))/(float(ntim)*float(nwave))
            if (info.talk eq 1 AND prg gt prg0) then begin
               print,format='(16(%"\b"),"  Progress ",I3," %",$)',prg
               prg0=prg
            endif 
         endfor ;ww nwave (III)
;print,nf
      
;now 'interpolate' data to real selected position, if it was not
;sampled in input data
         imafpipre=imafpi
         imafpi=fltarr(info.etalpos,sizima(2),sizima(3),sizima(4))
         mind=max(index)
         for ii=0,mind do begin
            indy=where(index eq ii,nindy)
            for bb=0,nindy-1 do imafpi[ii,*,*,*]=imafpi[ii,*,*,*]+(wg[indy(bb)]*imafpipre[indy(bb),*,*,*])
         endfor

; Header
         sxaddpar,head,"HISTORY",'Etalon Module'
         sxaddpar,head,"FRATIO",info.frat
         for ii=0,info.etalpos-1 do sxaddpar,head,"WVLTRAN"+strtrim(ii,2),format='(F8.3)',llo(ii),'Transmitted wavelen. (A) [incl. telec. bluesh]'
         for ii=0,info.etalpos-1 do sxaddpar,head,"FWHMREAL"+strtrim(ii,2),format='(F6.3)',fwhm_real(ii),'FWHM (A)'
         if (info.discaretal eq 1) then for net=0,n_elements(info.ndeathetal)-1 do sxaddpar,head,"TDEATET"+strtrim(net,2),info.ndeathetal(net)*info.tsamp,'Time for Etalon change'

         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(fix(nf-nout),2)+'.fits',imafpi,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(fix(nf-nout),2)+'.fits',imafpi,head
      endif 

   endfor ;nf
   print,''


;##########################################################################
; below, the dual-beam mode part. Repeating last part from above
; (i.e. the data part)

   if (info.dualbeam eq 1) then begin
      prg0=0.
      prg=prg0
      nst=0.
      lamind=0
      if (info.talk eq 1) then print,'Starting etalon dual-beam part'
      for nf=0,ntim-1 do begin
         if (poli eq -1 AND nf eq nst+info.nstate+info.ndeaths(0,lamind)) then begin
            nst=nst+info.nstate+info.ndeaths(0,lamind)
            nout=nout+info.ndeaths(0,lamind)
            lamind=lamind+1 
         endif 
         if ((poli eq -1 AND nf lt nst+info.nstate) OR (poli ne -1)) then begin
            ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(nf,2)+'.fits*',head,/compr,/sil)
            sizima=size(ima)
            imafpi=fltarr(nwave,sizima(2),sizima(3),sizima(4))
            for ww=0,nwave-1 do begin
;jjjjjjjjjjjjjj
;               norm=total(fff[*,ww],1)
               if (info.pffov eq 1) then norm=total(fff[*,*,*,ww],3) else norm=total(fff[*,ww],1)
;jjjjjjjjjjjjjj
               imafpipre=fltarr(sizima(2),sizima(3),sizima(4))
;hhhhhhhhhhhhhh
               if (info.pffov eq 1) then begin
                  for sto = 0, sizima(2)-1 do begin
                     for step = 0, info.wavepoint-1 do begin      
                        imafpipre[sto, *, *] = imafpipre[sto, *, *] + reform(ima[step, sto, *, *])*fff[*, *, step, ww] ;changed 14/12/16
                     endfor
                     imafpipre[sto,*,*]=reform(imafpipre[sto,*,*])/norm
                  endfor
                  imafpi[ww,*,*,*]=imafpipre   
               endif else begin
                  for step=0,info.wavepoint-1 do imafpipre = imafpipre + reform(ima[step,*,*,*])*fff[step,ww]
                  imafpi[ww,*,*,*]=imafpipre/norm
               endelse 
;hhhhhhhhhhhhhhh
               prg=round(100.*float(ww+1)*float(nf+1))/(float(ntim)*float(nwave))
               if (info.talk eq 1 AND prg gt prg0) then begin
                  print,format='(21(%"\b"),"  DUAL Progress ",I3," %",$)',prg
                  prg0=prg
               endif 
            endfor ;ww nwave (III)

;now 'interpolate' to real selected position, if it was not sampled in
;input data
;no need to interpolate llo and fwhm since it was done on the
;single-beam part just above and it carries
            imafpipre=imafpi
            imafpi=fltarr(info.etalpos,sizima(2),sizima(3),sizima(4))
            mind=max(index)
            for ii=0,mind do begin
               indy=where(index eq ii,nindy)
               for bb=0,nindy-1 do imafpi[ii,*,*,*]=imafpi[ii,*,*,*]+(wg[indy(bb)]*imafpipre[indy(bb),*,*,*])
            endfor

; Header
            sxaddpar,head,"HISTORY",'Etalon Module'
            sxaddpar,head,"FRATIO",format='(F8.3)',info.frat
            for ii=0,info.etalpos-1 do sxaddpar,head,"WVLTRAN"+strtrim(ii,2),format='(F8.3)',llo(ii),'Transmitted wavelen. (A) [incl. telec. bluesh]'
            for ii=0,info.etalpos-1 do sxaddpar,head,"FWHMREAL"+strtrim(ii,2),format='(F6.3)',fwhm_real(ii),'FWHM (A)'
            if (info.discaretal eq 1) then for net=0,n_elements(info.ndeathetal)-1 do sxaddpar,head,"TDEATET"+strtrim(net,2),info.ndeathetal(net)*info.tsamp,'Time for Etalon change'

            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf-nout,2)+'.fits',imafpi,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf-nout,2)+'.fits',imafpi,head
         endif 
      endfor ;nf
      print,''

;###################################################################


   endif 

endelse


end
