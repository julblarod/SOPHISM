PRO sophism_loadascii

; ==============================================================================
;
; DESCRIPTION
;    Loads SOPHISM settings specified in an ASCII file. Modifies and
;    adds some variables for the simulation run.
;
; CALLED FROM
;    sophism
;
; SUBROUTINES
;    headfits
;
; MAIN INPUTS
;    sophism_ascii.pro      ASCII file with the variables and values
;                              to be used in the simulation
; OUTPUT
;    Settings file ready to be used for simulation.
;
; VERSION HISTORY 
;    J. Blanco. 2011. v0_1. First steps
;    J. Blanco. 2012. v1_0. Major changes to whole routine
;    J. Blanco. Dec 2012. v2_0. Changes in some variables to follow
;       main routine. Included etalon voltage error, report option
;       to be run and modules for selection
;    J. Blanco. Jan 2013. v2_1. Corrected bug 'info.routines'
;    J. Blanco. Feb 2013. v2_2. Corrected ntim. Included check for
;       output folder existence. Corrected progctrl_ref
;    J. Blanco. May 2013. v2_3. Added variable obsset to replicate the
;       observation run.
;    J. Blanco. June 2014. v2_4. Added compression module, etalon
;       discarding, frames or samples discard scheme. Added check for
;       etalon/polarization module in data header, for ntim
;       purposes. Consider etalon speed for calculating discarded
;       frames as a user setting, not fixed.
;    J. Blanco. Jul 2014. v2_5. Cycle repetition correction. Keep only
;       one element of ndeathlcvr if not running polarization module. 
;    J. Blanco. Sep 2014. v2_6. Change of cycle repetition scheme to
;       associate it to array's Stokes dimension. In case of no
;       polarization, set cycles to 1.
;    J. Blanco. Nov 2014. v2_7. Corrected bug in defocus conversion 
;       from mm to rad, adding focal length setting.
;    J. Blanco. Jul 2015. v2_8. Corrected bug when calculating file
;       number, in case of starting simulation from the beginning
;
; ==============================================================================

@../settings/sophism_ascii.pro

;only for reference. Don't touch
;progctrl_ref=['sophism_input','sophism_jitter','sophism_polmeas','sophism_linbo_2d','sophism_otf_2d_b','sophism_papo_2d_b','sophism_fpa','sophism_accu','sophism_demod','sophism_inversion_2d','sophism_compression','sophism_report'] 
progctrl_ref=['sophism_input','sophism_jitter','sophism_polmeas','sophism_linbo','sophism_otf','sophism_papo','sophism_fpa','sophism_accu','sophism_demod','sophism_inversion','sophism_compression','sophism_report'] 
    
routines=fltarr(12)
;find the modules selected to be run
for tt=0,n_elements(progctrl)-1 do routines(where(progctrl(tt) eq progctrl_ref))=1
;Data files list
dataf=file_search(datafilo,count=numfil) 

;check if saves variable ends with / and add if not
if (strmid(saves,strlen(saves)-1,1) ne '/') then saves=saves+'/'
;check for saves folder. Create it if not existing
testsaves=file_search(saves,count=numsaves,/test_directory)
if (numsaves eq 0) then spawn,'mkdir '+saves

;'Replicated' spatial dimension (pix)
sz=sz0*sscene_mag
;Spatial sampling at Sun-S/C distance (arcsec)
sssamp=(sssamp0*1.e3)*(1./(scdis*1.495e11))*!radeg*3600. 

;choose if select a narrower range for the wavelength dimension of the
;input data (because some have +- 7 AA in 992 wavelength points)
case subset of
   'Complete': wavepoint=wavepointor
   'Subset':   wavepoint=wavepointor/2
   'Minimal':  wavepoint=wavepointor/4
endcase
;this could create problems for simulations loading previous data
;(progma ne -1), if the user doesn't enter the appropriate
;wavepointor (the original one, that is) or uses Subset instead of
;Complete for the new wavepoint array

;initialize to 0 the number of files, in case the if loop below is not entered
numprogmas=0
;update some variables calculated above when loading data obtained
;from a previous simulation run
etal=-1
poli=-1
if (max(progma) gt -1) then begin
   testprogma=file_search(saves+files(progma)+'_*.fits*',count=numprogmas)
;in case pupil apodization and loading linbo or otf, it won't
;find any fits but there should be a series, so force it. And read
;header from polarization files
   if (numprogmas eq 0) then begin
      if (routines(5) eq 1 AND (progma eq 3 OR progma eq 4)) then begin
         dataser=1 
         testprogma=file_search(saves+files(2)+'_*.fits*',count=numprogmas)
      endif else begin
         print,'No files found'
         retall
      endelse  
   endif 
   hh=headfits(testprogma(0))
   sz=sxpar(hh,'NAXIS3')
   if (fparesamp eq 1) then sssamp=fpaplate
   histo=sxpar(hh,'HISTORY')
   etal=where(histo eq 'Etalon Module')
   poli=where(histo eq 'Polarization Module')
endif

;adding the two first zeroes, of tip&tilt
zerco=[0,0,zerco] 
;convert focus displacement to zernike coefficient if needed
if (defocopt eq 1) then zerco(2)=!pi*zerco(2)/8./sqrt(3.)/wl/((flen/telap)^2.)

;replicate observations in case observation sets is higher than one
if (obsset gt 1) then begin 
   if (acqsche eq 'Fast Polarization') then begin 
      wavearr=sophism_repeat(wavearr,obsset)
      voltarr=sophism_repeat(voltarr,obsset)
   endif else cycrep=obsset
endif

;to create cycle repetitions
;depending on etalon mode (volts/wavelength) positions
if (wavevolt eq 1) then begin
;   wavearr=fix(strsplit(wavearrstr,',',/extract))
;cycle repetitions
;   if ((routines(2) eq 1 OR poli ne -1) AND cycrep gt 1) then begin
;      wavearrpre=wavearr
;      wavearr=0
;      for ll=0,n_elements(wavearrpre)-1 do wavearr=[wavearr,replicate(wavearrpre(ll),cycrep)]
;      wavearr=wavearr[1:*]
;   endif
   etalpos=n_elements(wavearr)
   lambdapos=wavearr
endif else begin
;   voltarr=fix(strsplit(voltarrstr,',',/extract))
;cycle repetitions
;   if ((routines(2) eq 1 OR poli ne -1) AND cycrep gt 1) then begin
;      voltarrpre=voltarr
;      voltarr=0
;      for ll=0,n_elements(voltarrpre)-1 do voltarr=[voltarr,replicate(voltarrpre(ll),cycrep)]
;      voltarr=voltarr[1:*]
;   endif 
   etalpos=n_elements(voltarr)
   lambdapos=voltarr
endelse 
 
;depending on etalon mode (volts/wavelength) positions, prepare variables
;if (wavevolt eq 1) then begin
;   etalpos=n_elements(wavearr) 
;   lambdapos=wavearr
;endif else begin 
;   etalpos=n_elements(voltarr)
;   lambdapos=voltarr
;endelse

;convert central wavelength to double
wl=double(wl)
;create the array of wavelengths for the data and to use everywhere else
lll=dindgen(wavepoint)*wavesamp*1d-3
lll=lll+wl-lll(wavepoint-1)/2.
;calibrate with doppler because of s/c velocity
lll=lll+lll*(scvel/299792.5d0)

;the longest time between exposure and readout, wins
frtim=etim > routtim
;number of time samples in a frame
nframe=fix(round(frtim*1e-3/tsamp))
;number of time samples in a modulation state
nstate=nacc*nframe


;total observation time is the state/frame time, times x modulation
;states, times wavel/volt positions, times cycle repetitions
;plus 'latence' times (times of LCVRs AND etalon state changes)
;Depending on discarding scheme, the number of images in ndeaths is
;different
if (discar eq 0) then begin 
;samples discarding,i.e.,discard only time samples not the whole frame 
   tdeathlcvr=[tdeath1,tdeath2,tdeath3,tdeath4]
   ndeathlcvr=fix(round(tdeathlcvr*1e-3/tsamp))
endif else begin
;frame discarding,i.e.,discard whole frame with all the corresponding samples
   tdeathlcvr=[tdeath1,tdeath2,tdeath3,tdeath4]
;as soon as tdeath is not 0, a frame is discarded. Compare tdeaths
;with frtim to see is more than one frame is discarded
   frdeaths=ceil(tdeathlcvr/frtim)
;calculate how many samples are discarded
   ndeathlcvr=fix(round(frdeaths*frtim*1e-3/tsamp))
;Possible future options: discard frame only when a certain % of frtim
;is 'lost', use camera on/off scheme instead of discarding (will need
;times of that... and synchronization in reality)
endelse

;cycles
if (cycrep gt 1) then begin
   ndeathlcvr=sophism_repeat(ndeathlcvr,cycrep)
   modnbuff=modnbuff*cycrep
endif

;calculate the voltage distance between spectral positions to obtain
;the time needed for each change and the discarded frames
difetalpos=abs(shift(lambdapos,-1)-lambdapos)
;tvolt=1500.    ;Volt/s
;if spectral positions given in wavelength, must change to voltage
;can it be left like here? From differences directly. Or should
;calculate voltages of the lambdas and then the difference?
if (wavevolt eq 1) then begin
   cll=[5250.d0,6328.d0]     ; CSIRO's manual
   cvv=[2.82d-5,3.44d-5]     ; nm/Volt
   llo_f=wl+wl*scvel/299792.5d0
   aet=interpol(cvv,cll,llo_f)*10. ;AA/Volt
   difetalpos=difetalpos*1e-3/aet
endif
if (discar eq 0) then tdeathetal=(difetalpos/tvolt)*1e3 else tdeathetal=ceil((difetalpos/tvolt)*1e3/frtim)*frtim
ndeathetal=fix(round(tdeathetal*1e-3/tsamp))

if (discaretal eq 0 OR routines(3) eq 0) then ndeathetal=ndeathetal*0.

;read a data header to check if it has gone through modulation and/or
;etalon not to put ones in modnbuff or etalpos below
;heady=headfits(saves+files(progma)+'_0.fits*')

;check for below if data have gone through Etalon and/or Polarization
;module, not to make etalpos and modnbuff = 1
;headdat=headfits(testprogma)
;headmod=sxpar(headdat,'HISTORY')
;etalhead=where(headmod eq 'Etalon Module')
;polihead=where(headmod eq 'Polarization Module')

;number of time samples in a complete modulation period
;if polarization module is not used, only 'one state' considered and
;no time expended in LCVRs changes
if (poli eq -1 AND routines(2) eq 0) then begin
   modnbuff=1
   cycrep=1
   ndeathlcvr=ndeathlcvr*0.
   ndeathlcvr=ndeathlcvr(0)
endif 

;number of time samples in complete observation
;if etalon module (and pupil apodization) is not used, only 'one
;wavelength' considered
if (etal eq -1 AND routines(3) eq 0 AND routines(5) eq 0) then begin
   etalpos=1
   ndeathetal=ndeathetal(0)
endif 

;prepare ndeaths matrix without repeating last lcvr discard and etalon
;discard
;always in the same way, whatever choosing for modulation and etalon
;discard or modules run 
ndeaths=fltarr(n_elements(ndeathlcvr),n_elements(ndeathetal))
for nlam=0,n_elements(ndeathetal)-1 do begin
   ndelcet=ndeathlcvr[n_elements(ndeathlcvr)-1] > ndeathetal[nlam]
   if (n_elements(ndeathlcvr) gt 1) then ndeaths[*,nlam]=[ndeathlcvr[0:n_elements(ndeathlcvr)-2],ndelcet] else ndeaths[*,nlam]=[ndelcet]
endfor 


;if (cycrep gt 1) then begin
;   modnbuff=modnbuff*cycrep
;   ndeathsmid=ndeaths
;   ndeaths=fltarr(modnbuff,etalpos*cycrep)
;   for ll=0,etalpos-1 do ndeaths[*,ll*cycrep:ll*cycrep+cycrep-1]=sophism_repeat(ndeathsmid[*,ll],cycrep)
;   for ll=0,etalpos-1 do ndeaths[*,ll*cycrep:ll*cycrep+cycrep-2]=ndeathlcvr(3)
;endif


;adding the discarded frames in the 'modulation dimension'. nperiod
;ends up as as an array of etalpos size 
nperiod=modnbuff*nstate+total(ndeaths,1)
;ntim=total(nperiod*cycrep)
ntim=total(nperiod)
;addendum to the ntim calculation in case of no etalon and polariz modules
;if (etal eq -1 AND poli eq -1 AND routines(2) eq 0 AND routines(3) eq 0 AND routines(5) eq 0) then ntim=ntim*obsset
if (etal eq -1 AND routines(3) eq 0 AND obsset gt 1) then begin
   ntim=ntim*obsset
   ndd=ndeaths
   for oo=1,obsset-1 do ndd=[[ndd],[ndeaths]]
   ndeaths=ndd
endif
;total time of observation
tobs=ntim*tsamp

stop
if (numfil gt 1 OR dataser eq 1 OR numprogmas gt 1) then ntim=round(tobs/tsamp) else ntim=1.

;modules to be included in the report
repomodules=intarr(11)
;find the modules selected to be run. Re-use the progctrl_ref since
;modules are the same (except report)
for tt=0,n_elements(repoctrl)-1 do repomodules(where(repoctrl(tt) eq progctrl_ref))=1



;estan mal pero vamos, para ir pensando
namvar=scope_varname(level=0)

info=create_struct(namvar(0),scope_varfetch(namvar(0)))

for tt=1,n_elements(namvar)-1 do begin
;intentar no guardar namvar, null, progctrl_ref?
   null=execute('if (n_elements(namvar(tt)) ne 0 AND namvar(tt) ne ''INFO'' AND namvar(tt) ne ''NAMVAR'' AND namvar(tt) ne ''NULL'' AND namvar(tt) ne ''PROGCTRL_REF'' AND namvar(tt) ne ''HEADDAT'' AND namvar(tt) ne ''HEADMOD'' AND namvar(tt) ne ''ETALHEAD'' AND namvar(tt) ne ''POLIHEAD'' AND namvar(tt) ne ''FRDEATHS'' AND namvar(tt) ne ''TDEATHLCVR'' AND namvar(tt) ne ''LAMBDAPOS'' AND namvar(tt) ne ''DIFETALPOS'' AND namvar(tt) ne ''CLL'' AND namvar(tt) ne ''CVV'' AND namvar(tt) ne ''AET'' AND namvar(tt) ne ''LLO_F'' AND namvar(tt) ne ''LL'' AND namvar(tt) ne ''WAVEARRPRE'' AND namvar(tt) ne ''VOLTARRPRE'' AND namvar(tt) ne ''TDEATHETAL'' AND namvar(tt) ne ''NDELCET'' AND namvar(tt) ne ''NDEATHSMID'' AND namvar(tt) ne ''NDD'' AND namvar(tt) ne ''OO'') then info=create_struct(info,namvar(tt),scope_varfetch(namvar(tt)))')
endfor

if (numfil eq 0 AND max(progma) eq -1) then begin
   print,'No files found'
   return
endif

save,filenam='../settings/settings_sophism.sav',info

END
