@sophism_init_aux
@sophism_aux

; ==============================================================================
;
; DESCRIPTION
;    Main routine for the SOPHISM simulator. 
;    Sets default values for the variables. Opens a widget for changing
;    those values or selecting previously saved settings.
;    Selection of modules to be run on the simulation, input data
;    files, filenames of outputs.
;
; USAGE INSTRUCTIONS
;    Compile and run. Change necessary settings, select modules and 'Start'
; 
; SUBROUTINES
;    sophism_input
;    sophism_jitter
;    sophism_polmeas
;    sophism_linbo
;    sophism_otf
;    sophism_papo
;    sophism_fpa
;    sophism_accu
;    sophism_demod
;    sophism_inversion
;    sophism_compression
;    sophism_report
;
; OUTPUT
;    settings_sophism.sav (and similar file with selected name)
;       Contains the structure-variable info with all the settings for the
;       simulation run.
;
; VERSION HISTORY 
;    J. Blanco. 2011
;    J. Blanco. 2012. 
;    J. Blanco. Nov 2012. Added Inversion and Report
;       modules. Added 'tuning jittering' in etalon.
;    J. Blanco. Dec 2012. Added verbose option. Added options to
;       report module to force report in modules not selected in the
;       run, in case they were run before and data are present.
;    J. Blanco. Dec 2012. v2_3. Added input fields for LCVR errors values
;    J. Blanco. Jan 2013. v2_4. Check for existing output
;       directory. Demodulation parameters split to mark there the
;       adhoc option. Included file selection for input data, fringes,
;       prefilter loading 
;    J. Blanco. Feb 2013. v2_5. Transmittance ready, filcuth for low
;       frequency cutoff Hinode filter, issunder for ISS
;       underperformance. Check if info.saves variable ends in /. When
;       recalculating ntim=tobs/tsamp, make round or may cause
;       problems. Added nogui option. Allow no dark and no photon
;       noise.
;    J. Blanco. Apr 2013. v2_6. Corrected bug when preparing
;       voltarr. Also with numprogmas. Also in defocus conversion from
;       mm to zernike coeff. in rad (missing 1e-7)
;    J. Blanco. May 2013. v2_7. Added variable obsset to replicate the
;       observation run.
;    J. Blanco. Oct 2013. v2_8. Corrected bug from obsset.
;    J. Blanco. Jan 2014. v3_0. Option for discarding whole frames
;       (not just the samples of the tdeath). Option to consider
;       discarding when changing etalon position. Option to enable or
;       disable the rearrange of Stokes dimension (I,V,Q,U--I,Q,U,V)
;       for modulation. Corrected bug for obsset.
;    J. Blanco. May 2014. v4_0. Added compression module to list and
;       settings. Added inversion continuum setting. Readout (frame)
;       time as modifiable setting.
;    J. Blanco. Jun 2014. v4_1. Added etalon speed for calculating
;       discarded frames. Split FPA settings in two. Photon noise
;       selection of generation, loading or nothing. Cycle repetition
;       correction. Added flat intensity range. Keep only one element
;       of ndeathlcvr if not running polarization module.
;    J. Blanco. Jul 2014. v4_2. Selection for generating/loading gain
;       table, modified variable name.
;    J. Blanco. Aug 2014. v4_3. Re-organization of polarization cycles
;       and observation sets variables and cases.
;    J. Blanco. Sep 2014. v4_4. Change of cycle repetition scheme to
;       associate it to array's Stokes dimension. In case of no
;       polarization, set cycles to 1.
;    J. Blanco. Oct 2014. v4_5. Update default LCVRs'
;       retardances. Corrected bug in defocus conversion 
;       from mm to rad, adding focal length setting. Change variable
;       name 'Wavelength dimension' to 'Wavelength extension' because
;       it was repeated.
;    J. Blanco. Nov 2014. v4_6. Added hot/dead pixels
;       generation. Added cosmic rays generator. Start option for
;       beginning the simulation in a given sample, not from the first
;       one (for times when simulation was cut in the middle of a
;       module run).
;    J. Blanco. May 2015. v4_7. Introduced factor for flux conversion
;       for input data, when it is not in erg/s/cm2/... Aesthetic
;       changes in GUI, making some fields only modifiable when
;       marking the corresponding options (e.g. start sample)
;    J. Blanco. Feb 2016. v4_8. Introduced birefringence settings,
;       including distance EW-Entrance in Global. 
;    J. Blanco. Apr 2016. v4_9. Modified fringes settings because of
;       new way of producing them
;    J. Blanco. Mar 2017. v5_0. Included settings for C-Milos
;       inversion, limited inversion FOV, and for FOV-prefilter
;
; ==============================================================================


PRO help,ev
dummy
end


FUNCTION save_info, ev 
  
;Description: This function implements the changes being made through
;the graphical interface to the settings.


  ; Get the 'stash' structure. 
  WIDGET_CONTROL, ev.TOP, GET_UVALUE=info 
  ; Get the value and user value from the widget that 
  ; generated the event. 
  WIDGET_CONTROL, ev.ID, GET_VALUE=val, GET_UVALUE=uval 
  ; Set the value of the correct field in the 'retStruct' 
  ; structure, using the field's index number (stored in 
  ; the widget's user value). 

case strmid(uval,0,3) of
 'ZER': begin
          uval=strmid(uval,3,3)
          case uval of
             'TIP': info.zerco(0:1)=float(strsplit(val,',',/ext))
             'DEF': info.zerco(2)=val
             'AST': info.zerco(3:4)=float(strsplit(val,',',/ext))
             'COM': info.zerco(5:6)=float(strsplit(val,',',/ext))
             'SPH': info.zerco(9)=val
             'COL': begin
                indiinfo=where(tag_names(info) eq 'ZERCOLOAD')
                typeinfo=size(info.(indiinfo),/type)
                max=n_elements(info.(indiinfo))
                info.(indiinfo)=(fix(val,type=typeinfo))(0:max-1)
             end 
             'COF': begin
                indiinfo=where(tag_names(info) eq 'ZERCOFILE')
                typeinfo=size(info.(indiinfo),/type)
                max=n_elements(info.(indiinfo))
                info.(indiinfo)=(fix(val,type=typeinfo))(0:max-1)               
             end
          endcase 
       end

 'BIR': begin
          if (uval eq 'BIREFSEL') then val=ev.value
          indiinfo=where(tag_names(info) eq uval)
          typeinfo=size(info.(indiinfo),/type)
          max=n_elements(info.(indiinfo))
          info.(indiinfo)=(fix(val,type=typeinfo))(0:max-1)
       end 

 'MUE': begin
          if (uval eq 'MUELLERRORS') then begin
             vali=['muellerramp1','muellerrpha1','muellerrang1']
             widget_control,widget_info(ev.top,find_by_uname=vali(ev.value)),sensitive=ev.select
             vali=['muellerramp2','muellerrpha2','muellerrang2']
             widget_control,widget_info(ev.top,find_by_uname=vali(ev.value)),sensitive=ev.select
          endif
          indiinfo=where(tag_names(info) eq uval)
          typeinfo=size(info.(indiinfo),/type)
          max=n_elements(info.(indiinfo))
          info.(indiinfo)=(fix(val,type=typeinfo))(0:max-1)
        end 

 'RET': begin
          uval2=strmid(uval,3,3)
          case uval2 of
             'DEG': val=float(strsplit(val,',',/ext))
             else: begin
                if (uval eq 'RETERRORS') then begin
                   vali=['reterramp1','reterrpha1','reterrang1']
                   widget_control,widget_info(ev.top,find_by_uname=vali(ev.value)),sensitive=ev.select
                   vali=['reterramp2','reterrpha2','reterrang2']
                   widget_control,widget_info(ev.top,find_by_uname=vali(ev.value)),sensitive=ev.select
                endif 
             end 
          endcase
          indiinfo=where(tag_names(info) eq uval)
          typeinfo=size(info.(indiinfo),/type)
          max=n_elements(info.(indiinfo))
          info.(indiinfo)=(fix(val,type=typeinfo))(0:max-1)
        end 

 'FPA': begin
          uval2=strmid(uval,3,7)
          if (uval2 eq 'FLATSEL') then val=ev.value
          indiinfo=where(tag_names(info) eq uval)
          typeinfo=size(info.(indiinfo),/type)
          max=n_elements(info.(indiinfo))
          info.(indiinfo)=(fix(val,type=typeinfo))(0:max-1)
       end 

 'FSS': info.files(0)=val
 'FJI': info.files(1)=val
 'FPO': info.files(2)=val
 'FSP': info.files(3)=val
 'FOT': info.files(4)=val
 'FPU': info.files(5)=val
 'FFP': info.files(6)=val
 'FAC': info.files(7)=val
 'FDE': info.files(8)=val
 else:  begin
          indiinfo=where(tag_names(info) eq uval)
          typeinfo=size(info.(indiinfo),/type)
          max=n_elements(info.(indiinfo))
          if (typeinfo eq 2 AND max gt 1 AND uval ne 'ROUTINES' AND uval ne 'REPOMODULES') then info.(indiinfo)=fix(strsplit(val,',',/extract)) else info.(indiinfo) = (fix(val,type=typeinfo))(0:max-1)
        end
endcase
  ; Reset the top-level widget's user value to the updated 
  ; 'stash' structure. 
  WIDGET_CONTROL, ev.TOP, SET_UVALUE=info 
END 


;***********************************************************


FUNCTION save_info_str,ev,vali=vall

;Description: This function implements the changes being made through
;the graphical interface to the settings.
;Just like the function above but for some special cases.

  WIDGET_CONTROL, ev.TOP, GET_UVALUE=info
  WIDGET_CONTROL, ev.ID, GET_VALUE=val, GET_UVALUE=uval

        indiinfo=where(tag_names(info) eq uval.tag)
        typeinfo=size(info.(indiinfo),/type)
        max=n_elements(info.(indiinfo))
        if (size(val,/type) eq 8) then begin
         val=vall
         widget_control,ev.id,set_uvalue=uval.tag
        endif
        info.(indiinfo)=(fix(val(ev.index),type=typeinfo))(0:max-1)

  WIDGET_CONTROL,ev.TOP,SET_UVALUE=info                             
END


;***********************************************************


PRO sophism_event,ev

;Description: interprets and performs changes in the widget, like
;blocking/un-blocking settings depending on selections and starts
;the simulation. In this last case, some changes to the settings are
;also performed before starting.

  common mapp,idmapped

  widget_control,ev.id,get_uvalue=uvalue


  IF (TAG_NAMES(ev, /STRUCTURE_NAME) EQ 'WIDGET_CONTEXT') $
    THEN BEGIN
 ; Obtain the widget ID of the context menu base.                            
    contextBase = WIDGET_INFO(ev.ID, FIND_BY_UNAME = 'contextMenu')
  ; Display the context menu and send its events to the                          ; other event handler routines                             
    WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X,ev.Y, contextBase
  ENDIF



  case uvalue(0) of
     'onoff':  begin
                 if (ev.value eq 'start') then begin
                    print,'Comienza la sim con' 
                    widget_control,ev.top,get_uvalue=info

;case of loading the settings through an ASCII file.
                    if (info.settloadasc eq 1) then begin
                       if (info.fsettasc ne '../settings/sophism_ascii.pro') then spawn,'cp '+info.fsettasc+' ../settings/sophism_ascii.pro'
                       sophism_loadascii
                       restore,'../settings/settings_sophism.sav'
                       progctrl=[info.progctrl]
                       for tt=0,n_elements(progctrl)-1 do call_procedure,progctrl(tt)
                       retall
                    endif

;'Replicated' spatial dimension (pix)
                    info.sz=info.sz0*info.sscene_mag
;Spatial sampling at Sun-S/C distance (arcsec)
                    info.sssamp=(info.sssamp0*1e3)*(1./(info.scdis*1.495e11))*!radeg*3600.

;create the variables for the 'levels' option of the jittering
;filter with the given groups of parameters
                    if (info.filtype eq 'levels') then begin
                       info=create_struct(info,'fillevel',float(strsplit(info.fillevelpre,',',/extract)))
                       info=create_struct(info,'filrange',float(strsplit(info.filrangepre,',',/extract)))
                    endif

;convert focus displacement to zernike coefficient if needed
                    if (info.defocopt eq 1) then info.zerco(2)=-!pi*info.zerco(2)/8./sqrt(3.)/(info.wl*1e-7)/((info.flen/info.telap)^2.)

;convert focus displacement from PD option to zernike coefficient if needed
;                    if (info.pdset eq 1) then info.pddefoc=-!pi*info.pddefoc/8./sqrt(3.)/(info.wl*1e-7)/(info.frat^2.)

;replicate observations in case observation sets is higher than one
                    wavearrstr=info.wavearrstr
                    voltarrstr=info.voltarrstr
                    if (info.obsset gt 1) then begin
                       if (info.acqsche eq 'Fast Polarization') then begin 
;                          wavearrstr=[wavearrstr,sophism_repeat(','+info.wavearrstr,info.obsset-1)]
                          for rep=1,info.obsset-1 do wavearrstr=wavearrstr+','+info.wavearrstr
;                             voltarrstr=[voltarrstr,sophism_repeat(','+info.voltarrstr,info.obsset-1)]
                          for rep=1,info.obsset-1 do voltarrstr=voltarrstr+','+info.voltarrstr
                       endif else info.cycrep=info.obsset
                    endif 

;to create cycle repetitions
;apply just in case of modulation enabled. No sense "repeating"
;modualtion cycles without modulation
;depending on etalon mode (volts/wavelength) positions
                    if (info.wavevolt eq 1) then begin
                       voltarr=fix(strsplit(voltarrstr,',',/extract))
                       wavearr=fix(strsplit(wavearrstr,',',/extract))
;cycle repetitions
;                       if (info.routines(2) eq 1 AND info.cycrep gt 1) then begin
;                          wavearrpre=wavearr
;                          wavearr=0
;                          for ll=0,n_elements(wavearrpre)-1 do wavearr=[wavearr,replicate(wavearrpre(ll),info.cycrep)]
;                          wavearr=wavearr[1:*]
;                       endif
                       info=create_struct(info,'voltarr',voltarr)
                       info=create_struct(info,'wavearr',wavearr)
                       info=create_struct(info,'etalpos',n_elements(info.wavearr))
                       lambdapos=info.wavearr
                    endif else begin
                       voltarr=fix(strsplit(voltarrstr,',',/extract))
                       wavearr=fix(strsplit(wavearrstr,',',/extract))
;cycle repetitions
;                       if (info.routines(2) eq 1 AND info.cycrep gt 1) then begin
;                          voltarrpre=voltarr
;                          voltarr=0
;                          for ll=0,n_elements(voltarrpre)-1 do voltarr=[voltarr,replicate(voltarrpre(ll),info.cycrep)]
;                          voltarr=voltarr[1:*]
;                       endif 
                       info=create_struct(info,'voltarr',voltarr)
                       info=create_struct(info,'wavearr',wavearr)
                       info=create_struct(info,'etalpos',n_elements(info.voltarr))
                       lambdapos=info.voltarr
                    endelse 

;depending on etalon mode (volts/wavelength) positions, prepare variables
;                    info=create_struct(info,'wavearr',fix(strsplit(wavearrstr,',',/extract))) 
;                    info=create_struct(info,'voltarr',fix(strsplit(voltarrstr,',',/extract)))
;                    if (info.wavevolt eq 1) then begin
;                       info=create_struct(info,'etalpos',n_elements(info.wavearr))
;                       lambdapos=info.wavearr
;                    endif else begin
;                       info=create_struct(info,'etalpos',n_elements(info.voltarr))
;                       lambdapos=info.voltarr
;                    endelse 
;the longest time between exposure and readout, wins
                    frtim=info.etim > info.routtim
;number of time samples in a frame
                    nframe=fix(round(frtim*1e-3/info.tsamp))
                    info=create_struct(info,'nframe',nframe)
;number of time samples in a modulation state
                    nstate=info.nacc*nframe
                    info=create_struct(info,'nstate',nstate)
;total observation time is the state/frame time, times x modulation
;states, times wavel/volt positions, times cycle repetitions
;plus 'latence' times (times of LCVRs AND etalon state changes)
;Depending on discarding scheme, the number of images in ndeaths is
;different
                    if (info.discar eq 0) then begin 
;sample discarding,i.e.,discard only time samples not the whole frame 
                       tdeathlcvr=[info.tdeath1,info.tdeath2,info.tdeath3,info.tdeath4]
                       ndeathlcvr=fix(round(tdeathlcvr*1e-3/info.tsamp))
                    endif else begin
;frame discarding,i.e.,discard whole frame with all the corresponding samples
                       tdeathlcvr=[info.tdeath1,info.tdeath2,info.tdeath3,info.tdeath4]
;as soon as tdeath is not 0, a frame is discarded. Compare tdeaths
;with frtim to see is more than one frame is discarded
                       frdeaths=ceil(tdeathlcvr/frtim)
;calculate how many samples are discarded
                       ndeathlcvr=fix(round(frdeaths*frtim*1e-3/info.tsamp))
;Possible future options: discard frame only when a certain % of frtim
;is 'lost', use camera on/off scheme instead of discarding (will need
;times of that... and synchronization in reality)
                    endelse

;cycles
                    if (info.cycrep gt 1) then begin
                       ndeathlcvr=sophism_repeat(ndeathlcvr,info.cycrep)
                       info.modnbuff=info.modnbuff*info.cycrep
                    endif

;calculate the voltage distance between spectral positions to obtain
;the time needed for each change and the discarded frames
                    difetalpos=abs(shift(lambdapos,-1)-lambdapos)
;                    tvolt=1500. ;Volt/s
;if spectral positions given in wavelength, must change to voltage
;can it be left like here? From differences directly. Or should
;calculate voltages of the lambdas and then the difference?
                    if (info.wavevolt eq 1) then begin
                       cll=[5250.d0,6328.d0] ; CSIRO's manual
                       cvv=[2.82d-5,3.44d-5]; nm/Volt
                       llo_f=info.wl
                       dop=info.scvel/299792.5d0
                       ldop=llo_f*dop
                       llo_f=llo_f+ldop
                       aet=interpol(cvv,cll,llo_f)*10. ;AA/Volt
                       difetalpos=difetalpos*1e-3/aet
                    endif
                    if (info.discar eq 0) then tdeathetal=(difetalpos/info.tvolt)*1e3 else tdeathetal=ceil((difetalpos/info.tvolt)*1e3/frtim)*frtim
                    ndeathetal=fix(round(tdeathetal*1e-3/info.tsamp))

                    if (info.discaretal eq 0 OR info.routines(3) eq 0) then ndeathetal=ndeathetal*0.

;number of time samples in a complete modulation period
;if polarization module is not used, only 'one state' considered and
;no time expended in LCVRs changes
                    if (info.routines(2) eq 0) then begin
                       info.modnbuff=1
                       info.cycrep=1
                       ndeathlcvr=ndeathlcvr*0.
                       ndeathlcvr=ndeathlcvr(0)
                    endif 

;number of time samples in complete observation
;if etalon module is not used, only 'one wavelength' considered
                    if (info.routines(3) eq 0) then begin
                       info.etalpos=1
                       ndeathetal=ndeathetal(0)
                    endif 

;maybe not really needed? Probably yes, to check ndeathetal if
;Polarization module has not been enabled, and so on
                    info=create_struct(info,'ndeathetal',ndeathetal)
                    info=create_struct(info,'ndeathlcvr',ndeathlcvr)

;prepare ndeaths matrix always in the same way, whatever choosing
;for modulation and etalon discard or modules run
                    ndeaths=fltarr(n_elements(ndeathlcvr),n_elements(ndeathetal))
                    for nlam=0,n_elements(ndeathetal)-1 do begin
                       ndelcet=ndeathlcvr[n_elements(ndeathlcvr)-1] > ndeathetal[nlam]
                       if (n_elements(ndeathlcvr) gt 1) then ndeaths[*,nlam]=[ndeathlcvr[0:n_elements(ndeathlcvr)-2],ndelcet] else ndeaths[*,nlam]=[ndelcet]
                    endfor 

;prepare ndeaths matrix, size variable depending on modules and discards
;                    ndeaths=fltarr(info.modnbuff,info.etalpos)
;                    for nlam=0,info.etalpos-1 do begin
;                       ndelcet=ndeathlcvr[info.modnbuff-1] > ndeathetal[nlam]
;                       if (info.modnbuff gt 1) then ndeaths[*,nlam]=[ndeathlcvr[0:info.modnbuff-2],ndelcet] else ndeaths[*,nlam]=[ndelcet]
;                    endfor

;                    info=create_struct(info,'ndeaths',ndeaths)

;                    nperiod=info.modnbuff*nstate+total(ndeaths)
;adding the discarded frames in the 'modulation dimension'. nperiod
;ends up as as an array of etalpos size 
                    nperiod=info.modnbuff*nstate+total(ndeaths,1)
;                    info=create_struct(info,'nperiod',nperiod)

;                    ntim=nperiod*info.cycrep*info.etalpos
;                    ntim=total(nperiod*info.cycrep)
                    ntim=total(nperiod)
;addendum to the ntim calculation in case of no etalon and polariz modules
;                    if (info.routines(2) eq 0 AND info.routines(3) eq 0 AND info.routines(5) eq 0) then ntim=ntim*info.obsset

                    if (info.routines(3) eq 0 AND info.obsset gt 1) then begin
                       ntim=ntim*info.obsset
                       ndd=ndeaths
                       for oo=1,info.obsset-1 do ndd=[[ndd],[ndeaths]]
                       ndeaths=ndd
                       nperiod=info.modnbuff*nstate+total(ndeaths,1)
                    endif

                    info=create_struct(info,'ndeaths',ndeaths)
                    info=create_struct(info,'nperiod',nperiod)

;total time of observation
                    info.tobs=ntim*info.tsamp

;convert quantum well from percentage to fraction
                    info.fpawell=info.fpawell/100.
                    
;choose if select a narrower range for the wavelength dimension of the
;input data (because some have +- 7 AA in 992 wavelength points)
                    case info.subset of
                       'Complete': info.wavepoint=info.wavepointor
                       'Subset': info.wavepoint=info.wavepointor/2.
                       'Minimal': info.wavepoint=info.wavepointor/4.
                    endcase
;create the array of wavelengths for the data and to use everywhere else
                    lll=dindgen(info.wavepoint)*info.wavesamp*1d-3
                    lll=lll+info.wl-lll(info.wavepoint-1)/2.
;calibrate with doppler because of s/c velocity
                    lll=lll+lll*(info.scvel/299792.5d0)
                    info=create_struct(info,'lll',lll)

;when loading previous settings file, consider new selection of
;modules and 'starting data'
                    if (info.settloadon eq 1) then begin
                       rout=info.routines
                       progm=info.progma
                       startmed=info.startmed
                       nstart=info.nstart
                       repo=info.repomodules
                       verb=info.talk
                       restore,info.fsettl
                       info.routines=rout
                       info.progma=progm
                       info.startmed=startmed
                       info.nstart=nstart
                       info.repomodules=repo
                       info.talk=verb
                       if (info.numfil eq 0 AND max(info.progma) eq -1) then begin
                          print,''
                          print,'No input files selected!'
                          print,''
                          break
                       endif 
                    endif else begin
;find input files defined in the settings
;                       datas=file_search(info.datafilo,count=numfil)
;                       info=create_struct(info,'dataf',strarr(numfil))
;                       info.numfil=numfil
;                       info.dataf=datas
                       idd=widget_info(ev.top,find_by_uname='listfiles')
                       widget_control,idd,get_uvalue=dataf
                       if (n_elements(dataf) gt 1) then begin
                          dataf=dataf(1:*) 
                       endif else begin
                          if (max(info.progma) eq -1) then begin
                             print,''
                             print,'No input files selected!'
                             print,''
                             break
                          endif 
                       endelse
                       info=create_struct(info,'dataf',dataf)
                       info.numfil=n_elements(dataf)
                    endelse 

                    numprogmas=-1
                    if (max(info.progma) gt -1) then begin
                       testprogma=file_search(info.saves+info.files(info.progma)+'_*.fits*',count=numprogmas)
;in case pupil apodization and loading linbo or otf, it won't
;find any fits but there should be a series, so force it. And read
;header from polarization files
                       if (numprogmas eq 0) then begin
                          if (info.routines(5) eq 1 AND (info.progma eq 3 OR info.progma eq 4)) then begin
                             info.dataser=1 
                             testprogma=file_search(info.saves+info.files(2)+'_*.fits*',count=numprogmas)
                          endif else begin
                             print,'No files found'
                             return
                          endelse  
                       endif 
;this, when loading data (progma) but creating a new set of parameters
                       hh=headfits(testprogma(0))
                       info.sscene_mag=sxpar(hh,'MAGNIF')
                       info.sz=sxpar(hh,'NAXIS3')
                       if (info.fparesamp eq 1) then info.sssamp=info.fpaplate
                    endif 

;in case of producing the report without selecting any modules in its
;option tab, use the same modules selected for simulation
                    if (info.routines(11) eq 1 AND max(info.repomodules) eq 0) then info.repomodules=info.routines(0:10) 
stop
;save settings in default file and user's one, destroy widget 
                    save,filena=info.fsett,info
                    save,filena='../settings/settings_sophism.sav',info
                    widget_control,ev.top,/destroy
                    progctrl=['sophism_input','sophism_jitter','sophism_polmeas','sophism_linbo','sophism_otf','sophism_papo','sophism_fpa','sophism_accu','sophism_demod','sophism_inversion','sophism_compression','sophism_report']

                    indprog=where(info.routines eq 1)

                    if (info.numfil gt 1 OR info.dataser eq 1 OR numprogmas gt 1) then begin
                       info.ntim=round(info.tobs/info.tsamp)
                       save,filena=info.fsett,info
                       save,filena='../settings/settings_sophism.sav',info
                    endif
;check if no modules are selected for the simulation
                    if (indprog[0] eq -1) then begin
                       print, 'No Modules selected'
                       break
                    endif else begin
;check if saves variable ends with / and add if not
                       if (strmid(info.saves,strlen(info.saves)-1,1) ne '/') then info.saves=info.saves+'/'
;check for saves folder. Create it if not existing
                       testsaves=file_search(info.saves,count=numsaves,/test_directory)
                       if (numsaves eq 0) then spawn,'mkdir '+info.saves
;start the simulation with the selected modules.
                       progctrl=progctrl(indprog)
                       for tt=0,n_elements(progctrl)-1 do call_procedure,progctrl(tt)
                    endelse
                 endif else begin
                       widget_control,ev.top,/destroy
                       print,'Hasta la vista'
                 endelse
               end

;some 'widget magic' to allow/block some setting when chaging selections
     'TABI':   begin
                 if (ev.tab eq 1) then begin
                    idd=widget_info(ev.top,find_by_uname='groupfiles')
                    widget_control,idd,map=0
                    widget_control,idmapped,map=1
                 endif else begin
                    widget_control,idmapped,map=0
                    idd=widget_info(ev.top,find_by_uname='groupfiles')
                    widget_control,idd,map=1
                 endelse
               end

     'APARTE': break

     'SETTLOADON': begin
                     idd=widget_info(ev.top,find_by_uname='groupsettingl2')
                     widget_control,idd,sensitive=ev.select
                     aa=save_info(ev)
                   end

     'SETTLOADASC': begin
                     idd=widget_info(ev.top,find_by_uname='groupsettasc2')
                     widget_control,idd,sensitive=ev.select
                     aa=save_info(ev)
                   end

     'STARTMED': begin
                 idd=widget_info(ev.top,find_by_uname='startgr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

     'GLOBAL': begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupglob')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'INPU':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupinpu')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'DATAFILO': begin
                   idd=widget_info(ev.top,find_by_uname='listfiles')
                   widget_control,idd,get_uvalue=dataf
                   widget_control,ev.id,get_value=sele
                   if n_elements(sele) gt 1 then dataf=[dataf,sele(1:*)] else begin
                      if (max(strmatch(dataf,sele)) eq 1) then dataf=dataf[where(strmatch(dataf,sele) eq 0)] else dataf=[dataf,sele]
                   endelse 
                   if (n_elements(dataf) gt 1) then begin
                      dataf2=dataf(1:*)
                      for ff=0,n_elements(dataf)-2 do begin
                         barr=strpos(dataf(ff+1),'/',/reverse_search)+1
                         dataf2(ff)=strmid(dataf(ff+1),barr,100)
                      endfor 
                   endif else dataf2=''
                   widget_control,idd,set_value=dataf2,set_uvalue=dataf
                 end

     'LISTFILES': begin
                    widget_control,ev.id,get_uvalue=dataf
                    if (n_elements(dataf) gt 1) then dataf=dataf[where(strmatch(dataf,dataf(ev.index+1)) eq 0)]
                    if (n_elements(dataf) gt 1) then begin
                       dataf2=dataf(1:*)
                       for ff=0,n_elements(dataf)-2 do begin
                          barr=strpos(dataf(ff+1),'/',/reverse_search)+1
                          dataf2(ff)=strmid(dataf(ff+1),barr,100)
                       endfor 
                    endif else dataf2=''
                    widget_control,ev.id,set_value=dataf2,set_uvalue=dataf
                  end 

     'JITT':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupjitt')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'FILTYPE':begin
                 indi=ev.value
                 filtgrarr=['groupcutout','groupcutoff','groupcuthin','grouplevels']
                 val=['cutout','cutoff','Hinode','levels']
                 indoff=where(filtgrarr ne filtgrarr(indi))
                 for jj=0,2 do begin
                    iddoff=widget_info(ev.top,find_by_uname=filtgrarr(indoff(jj)))
                    widget_control,iddoff,sensitive=0
                 endfor
                 idd=widget_info(ev.top,find_by_uname=filtgrarr(indi))
                 widget_control,idd,sensitive=1
                 uval={tag:uvalue}
                 widget_control,ev.id,set_uvalue=uval
                 ev=create_struct(ev,'index',ev.value)
                 aa=save_info_str(ev,vali=val)
               end

     'ISS':    begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupiss')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'TRANS':  begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='grouptrans')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'OPTIC':  begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupoptic')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'ZERCOLOAD': begin
                 idd=widget_info(ev.top,find_by_uname='zercogr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

     'FRONTNORM': begin
                 idd=widget_info(ev.top,find_by_uname='optfrontgr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

     'PREFIL':  begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupprefil')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'PFLOAD': begin
                 idd=widget_info(ev.top,find_by_uname='pfcurvegr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end

     'PFFOV': begin
                idd=widget_info(ev.top,find_by_uname='pffovgr')
                widget_control,idd,sensitive=ev.select
                aa=save_info(ev)
              end

     'SPECT':  begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupspect')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'WAVEVOLT': begin
                   val=['voltarrstr','wavearrstr']
                   widget_control,widget_info(ev.top,find_by_uname=val(ev.value)),sensitive=1
                   widget_control,widget_info(ev.top,find_by_uname=val((ev.value+1) mod 2)),sensitive=0
                 end

     'VOLTJIT': begin
                 idd=widget_info(ev.top,find_by_uname='tunjitgr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

;     'SIZDOPT': begin
;                 idd=widget_info(ev.top,find_by_uname='sizdgr')
;                 widget_control,idd,sensitive=ev.select
;                 aa=save_info(ev)
;               end 

     'PUPAP':  begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='grouppupap')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'POLARI': begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='grouppolari')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'BIREFTREE': begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupbirefringence')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'BIREFSEL': begin
                 birefselection=['Load','Generate','None']
                 birefnames=['settpolari1533','settpolari1544']
                 if (ev.value eq 'None') then begin
                    idd=widget_info(ev.top,find_by_uname=birefnames(0))
                    widget_control,idd,sensitive=0
                    idd=widget_info(ev.top,find_by_uname=birefnames(1))
                    widget_control,idd,sensitive=0
                 endif else begin
                    idd=widget_info(ev.top,find_by_uname=birefnames(where(ev.value eq birefselection)))
                    widget_control,idd,sensitive=ev.select
                 endelse
                 if (ev.value eq 'Generate') then begin
                    idd=widget_info(ev.top,find_by_uname='settpolari166')
                    widget_control,idd,sensitive=ev.select
;                    idd=widget_info(ev.top,find_by_uname='settpolari1662')
;                    widget_control,idd,sensitive=0
                 endif
                 aa=save_info(ev)
              end 

     'BIREFLOADMODEL': begin
                 idd=widget_info(ev.top,find_by_uname='settpolari1662')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
              end 

     'FPA': begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupfpa')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'FPAFLUXCONV': begin
                 idd=widget_info(ev.top,find_by_uname='fpafluxgr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

     'FPAPHOTNOI': begin
                     photnoiselection=['Load','Generate','None']
                     photnoinames=['settfpa43']
                     if (ev.value eq 'Load') then begin
                        idd=widget_info(ev.top,find_by_uname=photnoinames(0))
                        widget_control,idd,sensitive=1
                     endif else begin 
                        idd=widget_info(ev.top,find_by_uname=photnoinames(0))
                        widget_control,idd,sensitive=0
                     endelse 
                     aa=save_info(ev)
                   end

     'COSMIC': begin
                 idd=widget_info(ev.top,find_by_uname='cosmicgr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

     'DARKFLAT': begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupdarkflat')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'FPADARKSEL': begin
                     darkselection=['Load','Generate','None']
                     darknames=['settfpa73']
                     if (ev.value eq 'Load') then begin
                        idd=widget_info(ev.top,find_by_uname=darknames(0))
                        widget_control,idd,sensitive=1
                     endif else begin 
                        idd=widget_info(ev.top,find_by_uname=darknames(0))
                        widget_control,idd,sensitive=0
                     endelse 
                     aa=save_info(ev)
                   end

     'FPAGAINSEL': begin
                     gainselection=['Load','Generate','None']
                     gainnames=['settfpa83','settfpa84']
                     if (ev.value eq 'None') then begin
                        idd=widget_info(ev.top,find_by_uname=gainnames(0))
                        widget_control,idd,sensitive=0
                        idd=widget_info(ev.top,find_by_uname=gainnames(1))
                        widget_control,idd,sensitive=0
                     endif else begin 
                        idd=widget_info(ev.top,find_by_uname=gainnames(where(ev.value eq gainselection)))
                        widget_control,idd,sensitive=ev.select
                     endelse 
                     aa=save_info(ev)
                   end

     'FPAFLATSEL': begin
                     flatselection=['Load','Generate','None']
                     fpanames=['settfpa103','settfpa104']
                     if (ev.value eq 'None') then begin
                        idd=widget_info(ev.top,find_by_uname=fpanames(0))
                        widget_control,idd,sensitive=0
                        idd=widget_info(ev.top,find_by_uname=fpanames(1))
                        widget_control,idd,sensitive=0
                     endif else begin 
                        idd=widget_info(ev.top,find_by_uname=fpanames(where(ev.value eq flatselection)))
                        widget_control,idd,sensitive=ev.select
                     endelse 
                     if (ev.value eq 'Generate') then begin
                        idd=widget_info(ev.top,find_by_uname='flatrange')
                        widget_control,idd,sensitive=ev.select
                     endif
                     aa=save_info(ev)
                   end

     'FPAFRING': begin
                 idd=widget_info(ev.top,find_by_uname='fringegr')
                 widget_control,idd,sensitive=ev.select
                 idd=widget_info(ev.top,find_by_uname='fringesettgr')
                 widget_control,idd,sensitive=ev.select
                 aa=save_info(ev)
               end 

     'ACCU':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupaccu')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'DEMO':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupdemo')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'INVER':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupinver')
                 widget_control,idd,map=1
                 idmapped=idd
               end
 
     'INVERSMALL': begin
                    idd=widget_info(ev.top,find_by_uname='fovinvergr')
                    widget_control,idd,sensitive=ev.select
                    aa=save_info(ev)
                   end 

     'COMPR':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='groupcompr')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     'REPO':   begin
                 widget_control,idmapped,map=0
                 idd=widget_info(ev.top,find_by_uname='grouprepo')
                 widget_control,idd,map=1
                 idmapped=idd
               end

     else: break
  endcase

END


;***********************************************************


PRO sophism,nogui=nogui,idlfile=idlfile,asciifile=asciifile,settfile

;  if keyword_set(nogui) then sophism_nogui,idlfile=idlfile,asciifile=asciifile,settfile
case keyword_set(nogui) of
   0: begin

      common mapp,idmapped

;general definition of widget
      main=widget_base(title='SOPHISM Main',xsize=700,ysize=630,/row)
      sele=widget_base(main,/align_left,/col,ysiz=630)
      sett=widget_base(main,/frame,xsi=500,ysiz=620,/align_right,uname='settbase')
      optis=widget_tab(sele,location=0,uvalue='TABI',ysiz=329)
      modopti=widget_base(optis,title='Modules',uvalue='palm',/row)
;button group for selecting which modules will be run
      modu=cw_bgroup(modopti,['Input','Jittering','Polarization','Etalon','Optics','Pupil Apod.','FPA','Accumulation','Demodulation','Inversion','Compression','Report'],/nonexclusive,row=12,uvalue='ROUTINES',/return_index,space=.3,event_func='save_info')
      progma=cw_bgroup(modopti,strarr(12),/exclusive,row=12,uvalue='PROGMA',/return_index,space=5.3,event_func='save_info')
      opti=widget_base(optis,title='Parameters',uvalue='palo')
      menu=widget_tree(opti,/align_top,uvalue='APARTE',ysiz=328)
      startsamp=widget_base(sele,/row)
      butts=widget_base(sele)


;-----------------------------------------------------------------
;tree of modules definition, to enter settings in each
      glob=widget_tree(menu,value='Global Parameters',uvalue='GLOBAL')
      inpu=widget_tree(menu,value='Input Data',uvalue='INPU')
      jitt=widget_tree(menu,value='Jittering',uvalue='JITT')
      iss=widget_tree(menu,value='ISS',uvalue='ISS')
      pola=widget_tree(menu,value='Polarization',uvalue='POLARI',/folder)
          bireftree=widget_tree(pola,value='Birefringence',uvalue='BIREFTREE')
      trans=widget_tree(menu,value='Transmittances',uvalue='TRANS')
      optic=widget_tree(menu,value='General Optics',uvalue='OPTIC')
      prefil=widget_tree(menu,value='Prefilter',uvalue='PREFIL')
      spect=widget_tree(menu,value='Etalon',uvalue='SPECT')
      pupapod=widget_tree(menu,value='Pupil Apodisation',uvalue='PUPAP')
      fpa=widget_tree(menu,value='Focal Plane',uvalue='FPA',/folder)
         darkflat=widget_tree(fpa,value='Darks Flats',uvalue='DARKFLAT')
      accu=widget_tree(menu,value='Accumulation',uvalue='ACCU')
      demo=widget_tree(menu,value='Demodulation',uvalue='DEMO')
      inver=widget_tree(menu,value='Inversion',uvalue='INVER')
      compr=widget_tree(menu,value='Compression',uvalue='COMPR')
      repo=widget_tree(menu,value='Report',uvalue='REPO')


;-----------------------------------------------------------------
;starting/exiting buttons

      startmed=cw_bgroup(startsamp,['Start sample: '],/nonexclusive,row=1,uvalue='STARTMED',set_value=0);,/return_index,event_func='save_info')
      startsamp2=widget_base(startsamp,/row,sensitive=0,map=1,uname='startgr') ;,xsiz=450,ysiz=40)
      nstart=cw_field(startsamp2,title='',uvalue='NSTART',value=0,/ret,xsiz=5,uname='nstart',event_func='save_info')

      talk=cw_bgroup(butts,['Verbose'],/nonexclusive,row=1,uvalue='TALK',set_value=0,/return_index,yoffs=125,event_func='save_info')
      compress=cw_bgroup(butts,['Compress files'],/nonexclusive,row=1,uvalue='COMPRESS',set_value=1,/return_index,yoffs=155,event_func='save_info')
      onoff=cw_bgroup(butts,['Start','Exit'],/return_name,col=2,button_uvalue=['start','exit'],yoffs=190,uvalue='onoff')


;-----------------------------------------------------------------
;filenames for outputs of modules
      groupfiles=widget_base(sett,/col,uname='groupfiles',map=1,/scroll,x_scroll_size=460,y_scroll_size=595)
      saves=cw_field(groupfiles,title='Saves Folder:  ',uvalue='SAVES',value='../data/tests/',/ret,xsiz=45,event_func='save_info')
      fsscene=cw_field(groupfiles,title='Fsscene:       ',uvalue='FSSCENE',value='sscene',/ret,xsiz=45,event_func='save_info')
      fjitter=cw_field(groupfiles,title='Fjitter:       ',uvalue='FJITTER',value='jitscene',/ret,xsiz=45,event_func='save_info')
      fotfs=cw_field(groupfiles,title='Fotfs:         ',uvalue='FOTFS',value='otfscene',/ret,xsiz=45,event_func='save_info')
      fspecfilt=cw_field(groupfiles,title='Fspecfilt:     ',uvalue='FSPECFILT',value='specscene',/ret,xsiz=45,event_func='save_info')
      fpupilap=cw_field(groupfiles,title='Fpupilap:      ',uvalue='FPUPILAP',value='pupapscene',/ret,xsiz=45,event_func='save_info')
      fpol=cw_field(groupfiles,title='Fpol:          ',uvalue='FPOL',value='polscene',/ret,xsiz=45,event_func='save_info')
      ffpa=cw_field(groupfiles,title='Ffpa:          ',uvalue='FFPA',value='fpascene',/ret,xsiz=45,event_func='save_info')
      faccu=cw_field(groupfiles,title='Faccu:         ',uvalue='FACCU',value='accuscene',/ret,xsiz=45,event_func='save_info')
      fdemod=cw_field(groupfiles,title='Fdemod:        ',uvalue='FDEMOD',value='demodscene',/ret,xsiz=45,event_func='save_info')
      finver=cw_field(groupfiles,title='Finver:        ',uvalue='FINVER',value='inverscene',/ret,xsiz=45,event_func='save_info')
      fcompr=cw_field(groupfiles,title='Fcompr:        ',uvalue='FCOMPR',value='comprscene',/ret,xsiz=45,event_func='save_info')
      freport=cw_field(groupfiles,title='Report File:   ',uvalue='FREPORT',value='../reports/sophism_report.pdf',/ret,xsiz=45,event_func='save_info')
      fsett=cw_field(groupfiles,title='Settings File: ',uvalue='FSETT',value='../settings/settings_sophism.sav',/ret,xsiz=45,event_func='save_info')
      groupsettingl=widget_base(groupfiles,/row,uname='groupsettingl',map=1)
      settbutt=widget_base(groupsettingl,uname='settbutt',map=1)
      settloadon=cw_bgroup(settbutt,[''],label_left='Load Settings:',/nonexclusive,row=1,uvalue='SETTLOADON')
      groupsettingl2=widget_base(groupsettingl,uname='groupsettingl2',map=1,sensitive=0)
      fsettl=cw_field(groupsettingl2,title='',uvalue='FSETTL',value='../settings/settings_sophism.sav',/ret,xsiz=40,event_func='save_info')
      groupsettasc=widget_base(groupfiles,/row,uname='groupsettasc',map=1)
      settbuttasc=widget_base(groupsettasc,uname='settbuttasc',map=1)
      settloadasc=cw_bgroup(settbuttasc,[''],label_left='Load ASCII:   ',/nonexclusive,row=1,uvalue='SETTLOADASC')
      groupsettasc2=widget_base(groupsettasc,uname='groupsettasc2',map=1,sensitive=0)
      fsettasc=cw_field(groupsettasc2,title='',uvalue='FSETTASC',value='../settings/sophism_ascii.pro',/ret,xsiz=40,event_func='save_info')


;-----------------------------------------------------------------
;global parameters
      globgroup=widget_base(sett,/col,uname='groupglob',map=0)
      settglob1=widget_base(globgroup,/row,xsiz=450,ysiz=35)
      scvel=cw_field(settglob1,title='S/C Velocity (km/s): ',uvalue='SCVEL',value=0,/ret,xsiz=5,event_func='save_info')
      scdis=cw_field(settglob1,title='S/C Distance (AU): ',uvalue='SCDIS',value='0.28',/ret,xsiz=5,event_func='save_info')
;  wTextglob = WIDGET_list(settglob1, VALUE=['?'],/CONTEXT_EVENTS,uvalue='contexito',xsize=0.5,ysize=2,/frame,xoffset=15)
;  contextBaseglob = WIDGET_BASE(wTextglob,/CONTEXT_MENU,UNAME="contextMenu",uvalue='contex')
;  cb11 = WIDGET_BUTTON(contextBaseglob, VALUE = 'Acquisition scheme allows for a scanning order of fixing a wavelength position and then changing polarization',EVENT_PRO = 'help',scr_xsiz=10)
      settglob2=widget_base(globgroup,/row,xsiz=450,ysiz=35)
      sclat=cw_field(settglob2,title='S/C Latitude (deg): ',uvalue='SCLAT',value=0,/ret,xsiz=5,event_func='save_info')
      temp=cw_field(settglob2,title='Temperature (C): ',uvalue='TEMP',value=20,/ret,xsiz=5,event_func='save_info')
      settglob3=widget_base(globgroup,/row,xsiz=450,ysiz=40)
      acqsche=cw_form(settglob3,'0,DROPLIST,Fast Polarization|Fast Wavelength|other,label_left=Acquisition Scheme: ,event=save_info_str,tag=ACQSCHE',uvalue='ACQSCHE')
      settglob4=widget_base(globgroup,/row,xsiz=450,ysiz=40)
      telap=cw_field(settglob4,title='Telescope Aperture (m): ',uvalue='TELAP',value='0.14',/ret,xsiz=5,event_func='save_info')
      settglob9=widget_base(globgroup,/row,xsiz=450,ysiz=40)
      dhrew=cw_field(settglob9,title='Distance EW-Entrance Pupil (mm): ',uvalue='DHREW',value='322.',/ret,xsiz=5,event_func='save_info')
      settglob5=widget_base(globgroup,/row,xsiz=450,ysiz=40)      
      flen=cw_field(settglob5,title='Effective focal length (at detector) (m): ',uvalue='FLEN',value='4.125',/ret,xsiz=7,event_func='save_info')
      settglob6=widget_base(globgroup,/row,xsiz=450,ysiz=40)
      obsset=cw_field(settglob6,title='Observation Sets: ',uvalue='OBSSET',value=1,/ret,xsiz=3,event_func='save_info')  
      settglob7=widget_base(globgroup,/row,xsiz=450,ysiz=70)
      discar=cw_bgroup(settglob7,['Samples','Frames'],/exclusive,row=2,label_left='Discarding scheme: ',uvalue='DISCAR',space=6,set_value=0,/return_index,event_func='save_info')
      settglob8=widget_base(globgroup,/row,xsiz=450,ysiz=40)
      dualbeam=cw_bgroup(settglob8,['Dual Beam'],/nonexclusive,row=1,uvalue='DUALBEAM',/return_index,event_func='save_info')
;      settglob7=widget_base(globgroup,/row,xsiz=450,ysiz=40)
;      pdset=cw_bgroup(settglob7,['Phase Div.'],/nonexclusive,row=1,uvalue='PDSET',/return_index,event_func='save_info')
;      pddefoc=cw_field(settglob7,title='Focus displacement (mm): ',uvalue='PDDEFOC',value='10.',/ret,xsiz=4,event_func='save_info')


;-----------------------------------------------------------------
;input data parameters
      inpugroup=widget_base(sett,/col,uname='groupinpu',map=0)
      settinput1=widget_base(inpugroup,/row,xsize=500,ysize=189)
      datafilo=cw_filesel(settinput1,filter=['.fits','.fits.gz'],/fix_filter,/multiple,uvalue='DATAFILO')
      listfiles=widget_list(settinput1,scr_xsiz=190,uname='listfiles',value='',uvalue='LISTFILES')
      settinput15=widget_base(inpugroup,/row,xsize=450,ysize=35)
      dataser=cw_bgroup(settinput15,['Series'],/nonexclusive,row=1,uvalue='DATASER',/return_index,event_func='save_info')
      settinput2=widget_base(inpugroup,/row,xsiz=450,ysiz=35)
      sscene_mag=cw_field(settinput2,title='Spatial Replication Factor: ',uvalue='SSCENE_MAG',value='1.',/ret,xsiz=5,event_func='save_info')
      fparesamp=cw_bgroup(settinput2,['Resample to CCD pixel size'],/nonexclusive,row=1,uvalue='FPARESAMP',set_value=1,/return_index,event_func='save_info')
      settinput3=widget_base(inpugroup,/row,xsiz=450,ysiz=35)
      sz0=cw_field(settinput3,title='Spatial Dimension (pix): ',uvalue='SZ0',value=288,/ret,xsiz=5,event_func='save_info')
      sssamp0=cw_field(settinput3,title='Spatial Sampling (km): ',uvalue='SSSAMP0',value='20.8',/ret,xsiz=5,event_func='save_info')
      settinput4=widget_base(inpugroup,/row,xsiz=450,ysiz=35)
      sscene_apod=cw_field(settinput4,title='Apodization (%): ',uvalue='SSCENE_APOD',value='0.',/ret,xsiz=5,event_func='save_info')
      settinput5=widget_base(inpugroup,/row,xsiz=450,ysiz=35)
      tsampor=cw_field(settinput5,title='Time Sampling Original (s): ',uvalue='TSAMPOR',value='1.',/ret,xsiz=5,event_func='save_info')
      tsamp=cw_field(settinput5,title='Time Sampling Interp (s): ',uvalue='TSAMP',value='0.020',/ret,xsiz=5,event_func='save_info')
      settinput6=widget_base(inpugroup,/row,xsiz=450,ysiz=35)
      wavepointor=cw_field(settinput6,title='Wavelength Dimension: ',uvalue='WAVEPOINTOR',value=992,/ret,xsiz=5,event_func='save_info')
      wavesamp=cw_field(settinput6,title='Wavelength Sampling (mA): ',uvalue='WAVESAMP',value='14.127',/ret,xsiz=7,event_func='save_info')
      settinput7=widget_base(inpugroup,/row,xsiz=450,ysiz=40)
      wl=cw_field(settinput7,title='Central wavelength (A): ',uvalue='WL',value='6173.341',/ret,xsiz=9,event_func='save_info')
      geff=cw_field(settinput7,title='Effective Lande factor: ',uvalue='GEFF',value='2.5',/ret,xsiz=5,event_func='save_info')
      settinput8=widget_base(inpugroup,/row,xsiz=450,ysiz=40)
      subset=cw_form(settinput8,'0,DROPLIST,Minimal|Subset|Complete,label_left=Wavelength extension: ,tag=SUBSET,set_value=1,event=save_info_str',uvalue='SUBSET')


;-----------------------------------------------------------------
;jittering parameters  
      jittgroup=widget_base(sett,/col,uname='groupjitt',map=0)
      settjitt1=widget_base(jittgroup,/row,xsize=450,ysize=35)
      isson=cw_bgroup(settjitt1,['ISS'],/nonexclusive,row=1,uvalue='ISSON',/return_index,event_func='save_info',xsiz=100)
      jitoffset=cw_field(settjitt1,title='Jitter offset: ',uvalue='JITOFFSET',value='0.,0.',/ret,xsiz=7,event_func='save_info')
      settjitt2=widget_base(jittgroup,/row,xsize=450,ysize=40)
      jitrms=cw_field(settjitt2,title='Jitter rms: ',uvalue='JITRMS',value='0.5',/ret,xsiz=5,event_func='save_info')
      jittype=cw_form(settjitt2,'0,DROPLIST,white noise|randomn,label_left=Jitter type: ,event=save_info_str,tag=JITTYPE',uvalue='JITTYPE')
      settjitt3=widget_base(jittgroup,/row,xsiz=450,ysiz=40)
      filtype=cw_form(settjitt3,'0,DROPLIST,cutout|cutoff|Hinode|levels,label_left=Jitter Frequency Filter Type: ,tag=FILTYPE,set_value=2',uvalue='FILTYPE')
      filrms=cw_bgroup(settjitt3,['Orig. rms filtered'],/nonexclusive,row=1,uvalue='FILRMS',set_value=1,/return_index,event_func='save_info',xsiz=200)
      settfilt1=widget_base(jittgroup,/row,xsiz=450,ysiz=40,uname='groupcuthin')
      filcuth=cw_field(settfilt1,title='Low cut out freq (Hz)',uvalue='FILCUTH',value='0.1',/ret,xsiz=5,event_func='save_info')
      settfilt2=widget_base(jittgroup,/row,xsiz=450,ysiz=40,uname='groupcutout',sensitive=0)
      filcut1=cw_field(settfilt2,title='Cut out freq 1 (Hz)',uvalue='FILCUT1',value='0.',/ret,xsiz=5,event_func='save_info')
      filcut2=cw_field(settfilt2,title='Cut out freq 2 (Hz)',uvalue='FILCUT2',value='10.',/ret,xsiz=5,event_func='save_info')
      settfilt3=widget_base(jittgroup,/row,xsiz=450,ysiz=40,yoffs=-40,uname='groupcutoff',sensitive=0)
      filcutoff=cw_field(settfilt3,title='Cut-off freq (Hz)',uvalue='FILCUTOFF',value='10.',/ret,xsiz=5,event_func='save_info')
      settfilt4=widget_base(jittgroup,/row,xsiz=450,ysiz=40,yoffs=-40,uname='grouplevels',sensitive=0)
      fillevelpre=cw_field(settfilt4,title='Levels: ',uvalue='FILLEVELPRE',value='2.',/ret,xsiz=10,event_func='save_info')
      filrangepre=cw_field(settfilt4,title='Range: ',uvalue='FILRANGEPRE',value='0.,10.',/ret,xsiz=18,event_func='save_info')
      filsmooth=cw_field(settfilt4,title='Smooth: ',uvalue='FILSMOOTH',value=0,/ret,xsiz=4,event_func='save_info')


;-----------------------------------------------------------------
;iss parameters
      issgroup=widget_base(sett,/col,uname='groupiss',map=0)
      settiss1=widget_base(issgroup,/row,xsiz=350,ysiz=40)
      issatten=cw_form(settiss1,'0,DROPLIST,kis|other,label_left=Attenuation type: ,event=save_info_str,tag=ISSATTEN',uvalue='ISSATTEN')
      settiss2=widget_base(issgroup,/row,xsiz=350,ysiz=40)
      issunder=cw_field(settiss2,title='ISS Performance (%): ',uvalue='ISSUNDER',value='100.',/ret,xsiz=5,event_func='save_info')


;-----------------------------------------------------------------
;transmittances
      transgroup=widget_base(sett,/col,uname='grouptrans',map=0)
      setttrans1=widget_base(transgroup,/row,xsiz=450,ysiz=40)
      transsys=cw_field(setttrans1,title='System transmittance: ',uvalue='TRANSSYS',value='0.075',/ret,xsiz=5,event_func='save_info')


;-----------------------------------------------------------------
;optics parameters
      opticgroup=widget_base(sett,/col,uname='groupoptic',map=0)
      settoptic1=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      otftype=cw_form(settoptic1,'0,DROPLIST,diffraction|zemax|other,label_left=OTF Type: ,event=save_info_str,tag=OTFTYPE',uvalue='OTFTYPE')
      settoptic2=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      zertip=cw_field(settoptic2,title='Tip-Tilt:    ',uvalue='ZERTIP',value='0.,0.',/ret,xsiz=5,event_func='save_info')
      settoptic3=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      defocopt=cw_bgroup(settoptic3,['Zernike Coefficient: ','Displacement (mm): '],/exclus,col=3,uvalue='DEFOCOPT',/return_index,space=.3,label_left='Defocus ',set_value=0)
      zerdefoc=cw_field(settoptic3,title='',uvalue='ZERDEFOC',value='0.',/ret,xsiz=5,event_func='save_info')
      settoptic4=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      zerastig=cw_field(settoptic4,title='Astigmatism: ',uvalue='ZERASTIG',value='0.,0.',/ret,xsiz=5,event_func='save_info')
      settoptic5=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      zercoma=cw_field(settoptic5,title='Coma:        ',uvalue='ZERCOMA',value='0.,0.',/ret,xsiz=5,event_func='save_info')
      settoptic6=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      zerspher=cw_field(settoptic6,title='Spherical:   ',uvalue='ZERSPHER',value='0.',/ret,xsiz=5,event_func='save_info')
      settoptic7=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      zercoload=cw_bgroup(settoptic7,['Load Zernikes file: '],/nonexclusive,row=1,uvalue='ZERCOLOAD');,/return_index,event_func='save_info')
      settoptic72=widget_base(settoptic7,/row,sensitive=0,map=1,uname='zercogr')
      zercofile=cw_field(settoptic72,title='',uvalue='ZERCOFILE',value='../data/hrew_zernikes.txt',/ret,xsiz=42,event_func='save_info')
      settoptic8=widget_base(opticgroup,/row,xsiz=450,ysiz=40)
      frontnorm=cw_bgroup(settoptic8,['Norm. rms wavefront (lambda): '],/nonexclusive,row=1,uvalue='FRONTNORM');,/return_index,event_func='save_info')
      settoptic82=widget_base(settoptic8,/row,sensitive=0,map=1,uname='optfrontgr') ;,xsiz=450,ysiz=40)
      frontrms=cw_field(settoptic82,title='',uvalue='FRONTRMS',value='0.1',/ret,xsiz=5,uname='frontrms',event_func='save_info')


;-----------------------------------------------------------------
;prefilter parameters
      prefilgroup=widget_base(sett,/col,uname='groupprefil',map=0)
      settpref1=widget_base(prefilgroup,/row,xsiz=450,ysiz=40)
      pfwl=cw_field(settpref1,title='Prefilter Peak Wavelen (A): ',uvalue='PFWL',value='6173.34',/ret,xsiz=9,event_func='save_info')
      pfhwb=cw_field(settpref1,title='Prefilter FWHM (A): ',uvalue='PFHWB',value='3.',/ret,xsiz=4,event_func='save_info')
      settpref2=widget_base(prefilgroup,/row,xsiz=450,ysiz=40)
      ncav1=cw_field(settpref2,title='Pref Number of Cavities: ',uvalue='NCAV1',value=2,/ret,xsiz=4,event_func='save_info')
      fixed_pf=cw_form(settpref2,'0,DROPLIST,fixed|tunable,label_left=Prefilter Type: ,event=save_info_str,tag=FIXED_PF',uvalue='FIXED_PF')
      settpref3=widget_base(prefilgroup,/row);,xsiz=450,ysiz=40)
      pfload=cw_bgroup(settpref3,['Pref. Curve file: '],/nonexclusive,row=1,uvalue='PFLOAD')
      settpref32=widget_base(settpref3,/row,sensitive=0,map=1,uname='pfcurvegr');,xsiz=450,ysiz=40,sensitive=0,map=1)
      pfcurve=cw_field(settpref32,title='',uvalue='PFCURVE',value='../data/prefilter.sav',/ret,xsiz=43,uname='pfcurve',event_func='save_info')
      settpref4=widget_base(prefilgroup,/row,xsiz=450,ysiz=40)
      pffov=cw_bgroup(settpref4,['FOV prefilter variaion'],/nonexclusive,row=1,uvalue='PFFOV',uname='pffov') ;,/return_index,event_func='save_info')
      settpref42=widget_base(settpref4,/row,sensitive=0,map=1,uname='pffovgr')
      pffovmax=cw_field(settpref42,title='Max prefilter shift (mA): ',uvalue='PFFOVMAX',value='-1600.',/ret,xsiz=5,event_func='save_info')


;-----------------------------------------------------------------
;etalon parameters
      spectgroup=widget_base(sett,/col,uname='groupspect',map=0)
      settspect1=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      etal_temp=cw_field(settspect1,title='Temperature (C): ',uvalue='ETAL_TEMP',value='20.',/ret,xsiz=5,event_func='save_info')
      wavevolt=cw_bgroup(settspect1,['Voltage','Wavelength'],/exclus,col=3,uvalue='WAVEVOLT',space=.3,label_left='Mode: ',set_value=1)
      settspect2=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      voltarrstr=cw_field(settspect2,title='App. Volt. (V): ',uvalue='VOLTARRSTR',value='-500.,500.',/ret,xsiz=15,event_func='save_info',uname='voltarrstr')
      wavearrstr=cw_field(settspect2,title='Wavel. (mA): ',uvalue='WAVEARRSTR',value='-40,40',/ret,xsiz=13,event_func='save_info',uname='wavearrstr')
      settspect3=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      voltjit=cw_bgroup(settspect3,['Tuning Jittering rms (V): '],/nonexclusive,row=1,uvalue='VOLTJIT');,/return_index,event_func='save_info')
      settspect32=widget_base(settspect3,/row,sensitive=0,map=1,uname='tunjitgr') ;,xsiz=450,ysiz=40)
      voltjitrms=cw_field(settspect32,title='',uvalue='VOLTJITRMS',value='1.5',/ret,xsiz=5,uname='voltjitrms',event_func='save_info')
      settspect4=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      spectsizd=cw_field(settspect4,title='Etalon Thickness (micron): ',uvalue='SPECTSIZD',value='257.',/ret,xsiz=6,event_func='save_info')
;      settspect5=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
;      sizdopt=cw_bgroup(settspect5,['Thickness variation map: '],/nonexclusive,row=1,uvalue='SIZDOPT')
;      settspect52=widget_base(settspect5,/row,xsiz=480,ysiz=40,uname='sizdgr',sensitive=0,map=1)
;      sizdfil=cw_field(settspect52,title='',uvalue='SIZDFIL',value='../data/RD_FE1-2final_ThicknessMap1.txt',/ret,xsiz=43,uname='sizdfil',event_func='save_info')
;      spectsizdrms=cw_field(settspect5,title='Etalon Thickness rms (wl/rms): ',uvalue='SPECTSIZDRMS',value='0.',/ret,xsiz=5,event_func='save_info')
      settspect6=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      fabfin=cw_field(settspect6,title='Etalon fabrication finesse: ',uvalue='FABFIN',value='200.',/ret,xsiz=5,event_func='save_info')
      settspect7=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      refl=cw_field(settspect7,title='Reflectivity: ',uvalue='REFL',value='0.925',/ret,xsiz=5,event_func='save_info')
      reflrms=cw_field(settspect7,title='Reflectivity rms (%): ',uvalue='REFLRMS',value='0.',/ret,xsiz=5,event_func='save_info')
      settspect8=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      frat=cw_field(settspect8,title='F Number: ',uvalue='FRAT',value='56.5',/ret,xsiz=5,event_func='save_info')
      tilt=cw_field(settspect8,title='Etalon Incidence Angle (deg): ',uvalue='TILT',value='0.',/ret,xsiz=5,event_func='save_info')
      settspect9=widget_base(spectgroup,/row,xsiz=450,ysiz=40)
      discaretal=cw_bgroup(settspect9,['Etalon Discarding Enabled'],/nonexclusive,row=1,uvalue='DISCARETAL',/return_index,event_func='save_info')
      tvolt=cw_field(settspect9,title='Etalon Speed (V/s): ',uvalue='TVOLT',value='1500.',/ret,xsiz=5,event_func='save_info')


;-----------------------------------------------------------------
;pupil apodisation parameters
      pupapgroup=widget_base(sett,/col,uname='grouppupap',map=0)
      settpupap1=widget_base(pupapgroup,/row,xsiz=450,ysiz=40)
      intrange=cw_field(settpupap1,title='Integration range (mA):  ',uvalue='INTRANGE',value='3200.',/ret,xsiz=5,event_func='save_info')
      settpupap2=widget_base(pupapgroup,/row,xsiz=450,ysiz=40)
      phases=cw_bgroup(settpupap2,['Phases'],/nonexclusive,row=1,uvalue='PHASES',/return_index,event_func='save_info')


;-----------------------------------------------------------------
;polarisation parameters
      polarigroup=widget_base(sett,/col,uname='grouppolari',map=0,space=-2)
      settpolari1=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      pol_temp=cw_field(settpolari1,title='Temperature (C): ',uvalue='POL_TEMP',value='20.',/ret,xsiz=5,event_func='save_info')
      settpolari2=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      modtype=cw_form(settpolari2,'0,DROPLIST,longit. ideal|longit. real|verctor. ideal|vector. real|other,label_left=Modulation Type: ,event=save_info_str,tag=MODTYPE,set_value=2',uvalue='MODTYPE')
;  moddmod_type=cw_form(settpolari1,'0,DROPLIST,basic|correl,label_left=Demodulation Type: ,event=save_info_str,tag=MODDMOD_TYPE',uvalue='MODDMOD_TYPE')
      modnbuff=cw_field(settpolari1,title='Demodulated Image Buffers: ',uvalue='MODNBUFF',value=4,/ret,xsiz=3,event_func='save_info')
;  settpolari3=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      modspar=cw_form(settpolari2,'0,DROPLIST,V|Q|U,label_left=Longitud. Stokes Param: ,event=save_info_str,tag=MODSPAR',uvalue='MODSPAR')
;  modclip=cw_field(settpolari3,title='Cross-corr. clip: ',uvalue='MODCLIP',value='0.0',/ret,xsiz=5,event_func='save_info')
      settpolari4=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      mueller=cw_form(settpolari4,'0,DROPLIST,Identity|Real|Other,label_left=Mueller matrix: ,event=save_info_str,tag=MUELLER',uvalue='MUELLER')
      settpolari5=widget_base(polarigroup,/row,xsiz=450,ysiz=117)
      muellerrors=cw_bgroup(settpolari5,['Amp Ratio','Phase','Angle'],row=3,/nonexc,/return_index,uvalue='MUELLERRORS',space=10,label_left='Mueller Errors: ',event_func='save_info')
      settpolari52=widget_base(settpolari5,/col,xsiz=95,ysiz=117,space=-5)
      muellerramp1=cw_field(settpolari52,title='Min err:',uname='muellerramp1',uvalue='MUELLERRAMP1',value='0.05',/ret,xsiz=6,event_func='save_info')
      muellerrpha1=cw_field(settpolari52,title='Min err:',uname='muellerrpha1',uvalue='MUELLERRPHA1',value='5',/ret,xsiz=3,event_func='save_info')
      muellerrang1=cw_field(settpolari52,title='Min err:',uname='muellerrang1',uvalue='MUELLERRANG1',value='10',/ret,xsiz=3,event_func='save_info')
      settpolari53=widget_base(settpolari5,/col,xsiz=95,ysiz=117,space=-5)
      muellerramp2=cw_field(settpolari53,title='Max err:',uname='muellerramp2',uvalue='MUELLERRAMP2',value='0.1',/ret,xsiz=6,event_func='save_info')
      muellerrpha2=cw_field(settpolari53,title='Max err:',uname='muellerrpha2',uvalue='MUELLERRPHA2',value='10',/ret,xsiz=3,event_func='save_info')
      muellerrang2=cw_field(settpolari53,title='Max err:',uname='muellerrang2',uvalue='MUELLERRANG2',value='20',/ret,xsiz=3,event_func='save_info')
      settpolari6=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      retang1=cw_field(settpolari6,title='LCVR1 orientation (deg): ',uvalue='RETANG1',value='0.',/ret,xsiz=5,event_func='save_info')
      retang2=cw_field(settpolari6,title='LCVR2 orientation (deg): ',uvalue='RETANG2',value='45.',/ret,xsiz=5,event_func='save_info')
      settpolari7=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      retdeg1=cw_field(settpolari7,title='Retardances LCVR1 (deg): ',uvalue='RETDEG1',value='225.,225.,315.,315.',/ret,xsiz=20,event_func='save_info')
      settpolari8=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      retdeg2=cw_field(settpolari8,title='Retardances LCVR2 (deg): ',uvalue='RETDEG2',value='235.,125.,55.,305.',/ret,xsiz=20,event_func='save_info')
      settpolari9=widget_base(polarigroup,/row,xsiz=450,ysiz=117)
      reterrors=cw_bgroup(settpolari9,['Amp Ratio','Phase','Angle'],row=3,/nonexc,/return_index,uvalue='RETERRORS',space=10,label_left='LCVR Errors: ',event_func='save_info')
      settpolari92=widget_base(settpolari9,/col,xsiz=95,ysiz=117,space=-5)
      reterramp1=cw_field(settpolari92,title='Min err:',uname='reterramp1',uvalue='RETERRAMP1',value='0.05',/ret,xsiz=6,event_func='save_info')
      reterrpha1=cw_field(settpolari92,title='Min err:',uname='reterrpha1',uvalue='RETERRPHA1',value='5',/ret,xsiz=3,event_func='save_info')
      reterrang1=cw_field(settpolari92,title='Min err:',uname='reterrang1',uvalue='RETERRANG1',value='10',/ret,xsiz=3,event_func='save_info')
      settpolari93=widget_base(settpolari9,/col,xsiz=95,ysiz=117,space=-5)
      reterramp2=cw_field(settpolari93,title='Max err:',uname='reterramp2',uvalue='RETERRAMP2',value='0.1',/ret,xsiz=6,event_func='save_info')
      reterrpha2=cw_field(settpolari93,title='Max err:',uname='reterrpha2',uvalue='RETERRPHA2',value='10',/ret,xsiz=3,event_func='save_info')
      reterrang2=cw_field(settpolari93,title='Max err:',uname='reterrang2',uvalue='RETERRANG2',value='20',/ret,xsiz=3,event_func='save_info')
      settpolari10=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      tdeath1=cw_field(settpolari10,title='Time Death 1 (ms): ',uvalue='TDEATH1',value=5,/ret,xsiz=5,event_func='save_info')
      tdeath2=cw_field(settpolari10,title='Time Death 2 (ms): ',uvalue='TDEATH2',value=94,/ret,xsiz=5,event_func='save_info')
      settpolari11=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      tdeath3=cw_field(settpolari11,title='Time Death 3 (ms): ',uvalue='TDEATH3',value=94,/ret,xsiz=5,event_func='save_info')
      tdeath4=cw_field(settpolari11,title='Time Death 4 (ms): ',uvalue='TDEATH4',value=80,/ret,xsiz=5,event_func='save_info')
      settpolari12=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      cycrep=cw_field(settpolari12,title='Cycle Repetitions: ',uvalue='CYCREP',value=1,/ret,xsiz=5,event_func='save_info')
      settpolari13=widget_base(polarigroup,/row,xsiz=450,ysiz=40)
      stokreorder=cw_bgroup(settpolari13,['Reorder Stokes dimension I,V,Q,U --> I,Q,U,V'],/nonexclusive,row=1,uvalue='STOKREORDER',set_value=1,/return_index,event_func='save_info')


      birefringroup=widget_base(sett,/col,uname='groupbirefringence',map=0)
      settpolari15=widget_base(birefringroup,/row,xsiz=450,ysiz=110,uname='settpolari15') ;,sensitive=0)
      birefsel=cw_bgroup(settpolari15,['Load','Generate','None'],row=3,/exc,/return_name,uvalue='BIREFSEL',uname='BIREFSEL',set_value=2,space=13,label_left='Birefringence: ')
      settpolari152=widget_base(settpolari15,/col,xsiz=450,ysiz=75)
      settpolari153=widget_base(settpolari152,/col,xsiz=450,ysiz=35,uname='settpolari1533',sensitive=0)
      birefloadname=cw_field(settpolari153,title='',uvalue='BIREFLOADNAME',value='../data/biref_mat.sav',/ret,xsiz=40,event_func='save_info')
      settpolari154=widget_base(settpolari152,/row,xsiz=450,ysiz=35,uname='settpolari1544',sensitive=0)
      birefptv=cw_field(settpolari154,title='Peak to Valley: ',uvalue='BIREFPTV',value='0.05',/ret,xsiz=5,event_func='save_info')
      birefscale=cw_field(settpolari154,title='Scale: ',uvalue='BIREFSCALE',value='32',/ret,xsiz=5,event_func='save_info')
      settpolari17=widget_base(birefringroup,/row,xsiz=450,ysiz=40,uname='settpolari176');,sensitive=0)
      birefvarmodel=cw_field(settpolari17,title='Variation model: ',uvalue='BIREFVARMODEL',value=0,/ret,xsiz=5,event_func='save_info')
      settpolari16=widget_base(birefringroup,/row,xsiz=450,ysiz=40,uname='settpolari166',sensitive=0)
      birefloadmodel=cw_bgroup(settpolari16,['Load Birefingence model: '],/nonexclusive,row=1,uvalue='BIREFLOADMODEL',uname='birefloadmodel');,/return_index,event_func='save_info')
      settpolari162=widget_base(settpolari16,/row,xsiz=300,ysiz=30,uname='settpolari1662',sensitive=0,map=1)
      birefloadmodelname=cw_field(settpolari162,title='',uvalue='BIREFLOADMODELNAME',value='../data/biref_model.sav',/ret,xsiz=40,uname='birefloadmodelname',event_func='save_info')

;-----------------------------------------------------------------
;fpa parameters
      fpagroup=widget_base(sett,/col,uname='groupfpa',map=0)
      settfpa1=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      fpa_temp=cw_field(settfpa1,title='Temperature (C): ',uvalue='FPA_TEMP',value='20.',/ret,xsiz=5,event_func='save_info')
      settfpa20=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
;      fpapixsiz=cw_field(settfpa20,title='Pixel Size (microns): ',uvalue='FPAPIXSIZ',value='10.',/ret,xsiz=5,event_func='save_info')
      fpaplate=cw_field(settfpa20,title='Plate Scale (arcsec/pix): ',uvalue='FPAPLATE',value='0.5',/ret,xsiz=5,event_func='save_info')
      settfpa2=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      fpashutt=cw_form(settfpa2,'0,DROPLIST,snapshot|rolling,label_left=Shutter type: ,set_value=1,event=save_info_str,tag=FPASHUTT',uvalue='FPASHUTT')
      settfpa3=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      fpasatval=cw_field(settfpa3,title='Saturation value: ',uvalue='FPASATVAL',value='65535.',/ret,xsiz=7,event_func='save_info')
      fpafluxconv=cw_bgroup(settfpa3,['Flux conv.: '],/nonexclusive,row=1,uvalue='FPAFLUXCONV');,/return_index,event_func='save_info')
      settfpa32=widget_base(settfpa3,/row,sensitive=0,map=1,uname='fpafluxgr');,xsiz=450,ysiz=40)
      fpafluxfactor=cw_field(settfpa32,title='',uvalue='FPAFLUXFACTOR',uname='fpafluxfactor',value='3e14',/ret,xsiz=7,event_func='save_info')
      settfpa4=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      etim=cw_field(settfpa4,title='Exposure Time (ms): ',uvalue='ETIM',value='20.',/ret,xsiz=5,event_func='save_info')
      routtim=cw_field(settfpa4,title='Readout Time (ms): ',uvalue='ROUTTIM',value='100.',/ret,xsiz=5,event_func='save_info')
      settfpa5=widget_base(fpagroup,/row,xsiz=450,ysiz=110)
;      fpaphotnoi=cw_bgroup(settfpa4,['Add Photon Noise'],/nonexclusive,row=1,uvalue='FPAPHOTNOI',set_value=1,/return_index,event_func='save_info')
      fpaphotnoi=cw_bgroup(settfpa5,['Load','Generate','None'],row=3,/exc,/return_name,uvalue='FPAPHOTNOI',uname='FPAPHOTNOI',set_value=1,space=13,label_left='Photon Noise: ')
      settfpa52=widget_base(settfpa5,/col,xsiz=450,ysiz=75)
      settfpa53=widget_base(settfpa52,/col,xsiz=450,ysiz=35,uname='settfpa43',sensitive=0)
      fpaphotnoiname=cw_field(settfpa53,title='',uvalue='FPAPHOTNOINAME',value='../data/fpascene_phot.save',/ret,xsiz=40,event_func='save_info')
      settfpa6=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      fpaquan=cw_field(settfpa6,title='Quantum Eff: ',uvalue='FPAQUAN',value='0.5',/ret,xsiz=5,event_func='save_info')
      fpawell=cw_field(settfpa6,title='Full Well (%): ',uvalue='FPAWELL',value='50.',/ret,xsiz=5,event_func='save_info')
      settfpa7=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      noitype=cw_form(settfpa7,'0,DROPLIST,normal|other,label_left=Readout Noise Sigma: ,event=save_info_str,tag=NOITYPE',uvalue='NOITYPE')
      noisigma=cw_field(settfpa7,title='Readout Noise Sigma: ',uvalue='NOISIGMA',value='70.',/ret,xsiz=7,event_func='save_info')
;      settfpa6=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
;      fpadark=cw_field(settfpa6,title='Dark Array: ',uvalue='FPADARK',value='../data/',/ret,xsiz=45,event_func='save_info')
      settfpa8=widget_base(fpagroup,/row,xsiz=450,ysiz=40)
      fpadeadpix=cw_field(settfpa8,title='Dead pixels number: ',uvalue='FPADEADPIX',value=0,/ret,xsiz=4,event_func='save_info')
      fpahotpix=cw_field(settfpa8,title='Hot pixels number: ',uvalue='FPAHOTPIX',value=0,/ret,xsiz=4,event_func='save_info')
      settfpa9=widget_base(fpagroup,/row);,xsiz=450,ysiz=40)
      cosmic=cw_bgroup(settfpa9,['Cosmic ray hits/sec: '],/nonexclusive,row=1,uvalue='COSMIC');,/return_index,event_func='save_info')
      settfpa92=widget_base(settfpa9,/row,sensitive=0,map=1,uname='cosmicgr');,xsiz=450,ysiz=40)
      cosmhit=cw_field(settfpa92,title='',uvalue='COSMHIT',value='16.8',/ret,xsiz=5,uname='cosmhit',event_func='save_info')

      darkflatgroup=widget_base(sett,/col,uname='groupdarkflat',map=0)
      settfpa10=widget_base(darkflatgroup,/row,xsiz=450,ysiz=110)
      fpadarksel=cw_bgroup(settfpa10,['Load','Generate','None'],row=3,/exc,/return_name,uvalue='FPADARKSEL',uname='FPADARKSEL',set_value=1,space=13,label_left='Dark: ')
      settfpa102=widget_base(settfpa10,/col,xsiz=450,ysiz=75)
      settfpa103=widget_base(settfpa102,/col,xsiz=450,ysiz=35,uname='settfpa73',sensitive=0)
      fpadarkname=cw_field(settfpa103,title='',uvalue='FPADARKNAME',value='../data/dark.save',/ret,xsiz=40,event_func='save_info')
      settfpa11=widget_base(darkflatgroup,/row,xsiz=450,ysiz=110)
      fpagainsel=cw_bgroup(settfpa11,['Load','Generate','None'],row=3,/exc,/return_name,uvalue='FPAGAINSEL',uname='FPAGAINSEL',set_value=1,space=13,label_left='Gain Table: ')
      settfpa112=widget_base(settfpa11,/col,xsiz=450,ysiz=75)
      settfpa113=widget_base(settfpa112,/col,xsiz=450,ysiz=35,uname='settfpa83',sensitive=0)
      fpagainname=cw_field(settfpa113,title='',uvalue='FPAGAINNAME',value='../data/gain.save',/ret,xsiz=40,event_func='save_info')
      settfpa114=widget_base(settfpa112,/col,xsiz=450,ysiz=35,uname='settfpa84',sensitive=1)      
      fpagain=cw_field(settfpa114,title='Gain Table rms: ',uvalue='FPAGAIN',value='0.01',/ret,xsiz=5,event_func='save_info')
      settfpa12=widget_base(darkflatgroup,/row,xsiz=450,ysiz=110)
      fpaflatsel=cw_bgroup(settfpa12,['Load','Generate','None'],row=3,/exc,/return_name,uvalue='FPAFLATSEL',uname='FPAFLATSEL',set_value=2,space=13,label_left='Flat: ')
      settfpa122=widget_base(settfpa12,/col,xsiz=450,ysiz=75)
      settfpa123=widget_base(settfpa122,/col,xsiz=450,ysiz=35,uname='settfpa103',sensitive=0)
      fpaflatname=cw_field(settfpa123,title='',uvalue='FPAFLATNAME',value='../data/flat.save',/ret,xsiz=40,event_func='save_info')
      settfpa124=widget_base(settfpa122,/col,xsiz=450,ysiz=40,uname='settfpa104',sensitive=0)
      fpaflatmore=cw_bgroup(settfpa124,['for each wavelength','for each pol. state'],/nonexclusive,row=1,uvalue='FPAFLATMORE',event_func='save_info')
      settfpa13=widget_base(darkflatgroup,/row,xsiz=450,ysiz=40,uname='flatrange',sensitive=0,map=1)
;      fpaflorder=cw_field(settfpa11,title='Generated flat order: ',uvalue='FPAFLORDER',value=1,/ret,xsiz=3,event_func='save_info')
      fpaflrang=cw_field(settfpa13,title='Generated flat range (norm): ',uvalue='FPAFLRANG',value=0.3,/ret,xsiz=3,event_func='save_info')
      settfpa14=widget_base(darkflatgroup,/row,xsiz=470,ysiz=40)
      fpafring=cw_bgroup(settfpa14,['Fringes'],/nonexclusive,row=1,uvalue='FPAFRING');,/return_index,event_func='save_info')
      settfpa142=widget_base(settfpa14,/row,sensitive=0,map=1,uname='fringegr');,xsiz=450,ysiz=40)
;      fpafringamp=cw_field(settfpa142,title='Amplitude (%) ',uvalue='FPAFRINGAMP',value='1.',/ret,xsiz=5,uname='fpafringamp',event_func='save_info')
;      fpafringdir=cw_field(settfpa142,title='Direction ',uvalue='FPAFRINGDIR',value='0.',/ret,xsiz=4,uname='fpafringdir',event_func='save_info')
;      fpafringwidth=cw_field(settfpa142,title='Width (pi)',uvalue='FPAFRINGWIDTH',value='1.',/ret,xsiz=4,uname='fpafringwidth',event_func='save_info')
      fpafrinn=cw_field(settfpa142,title='Refrac. index ',uvalue='FPAFRINN',value='1.',/ret,xsiz=5,uname='fpafrinn',event_func='save_info')
      fpafrinref=cw_field(settfpa142,title='Reflectivity ',uvalue='FPAFRINREF',value='0.3',/ret,xsiz=4,uname='fpafrinref',event_func='save_info')
      settfpa14b=widget_base(darkflatgroup,/row,xsiz=10,ysiz=40)
      settfpa143=widget_base(settfpa14b,/row,sensitive=0,map=1,uname='fringesettgr',xoffs=4)
      fpafrinsizd=cw_field(settfpa143,title='            Thickness ',uvalue='FPAFRINSIZD',value='200.',/ret,xsiz=4,uname='fpafrinsizd',event_func='save_info')
      fpafrintheta=cw_field(settfpa143,title='Incidence angle ',uvalue='FPAFRINTHETA',value='0.',/ret,xsiz=4,uname='fpafrintheta',event_func='save_info')
;      settfpa15=widget_base(darkflatgroup,/row,xsiz=470,ysiz=40,sensitive=0,map=1,uname='fringesettgr2')
;;      fpafringtime=cw_bgroup(settfpa13,['Variable fringes'],/nonexclusive,row=1,uvalue='FPAFRINGTIME',/return_index,event_func='save_info')
;      fpafringtime=cw_bgroup(settfpa15,['Time','Lambda','Polariz.'],col=3,/nonexc,/return_index,uvalue='FPAFRINGTIME',space=10,label_left='Variable Fringes: ',event_func='save_info')

;-----------------------------------------------------------------
;accumulation parameters
      accugroup=widget_base(sett,/col,uname='groupaccu',map=0)
      settaccu1=widget_base(accugroup,/row,xsiz=450,ysiz=40)
      nacc=cw_field(settaccu1,title='Number of Accumulations: ',uvalue='NACC',value=12,/ret,xsiz=5,event_func='save_info')
;  wText = WIDGET_TEXT(settaccu1, VALUE='?',/ALL_EVENTS,/CONTEXT_EVENTS,uvalue='contexito',xsize=1,ysize=4,/frame)
;  contextBase1 = WIDGET_BASE(wText,/CONTEXT_MENU,UNAME="contextMenu",uvalue='contex')
;  cb11 = WIDGET_BUTTON(contextBase1, VALUE = 'Text selection 1',EVENT_PRO = 'help')


;-----------------------------------------------------------------
;demodulation parameters
      demogroup=widget_base(sett,/col,uname='groupdemo',map=0)
      settdemo1=widget_base(demogroup,/row,xsiz=450,ysiz=40)
      adhoc=cw_bgroup(settdemo1,['Software correction crosstalk from V to Q, U'],/nonexclusive,row=1,uvalue='ADHOC',/return_index,event_func='save_info')


;-----------------------------------------------------------------
;inversion parameters
      invergroup=widget_base(sett,/col,uname='groupinver',map=0)
      settinver1=widget_base(invergroup,/row,xsiz=450,ysiz=40)
      inversigma=cw_field(settinver1,title='Noise threshold factor: ',uvalue='INVERSIGMA',value=3,/ret,xsiz=5,event_func='save_info')
      iter=cw_field(settinver1,title='Maximum iterations: ',uvalue='ITER',value=15,/ret,xsiz=5,event_func='save_info')
      settinver2=widget_base(invergroup,/row,xsiz=450,ysiz=40)
      invercont=cw_field(settinver2,title='Index position of continuum: ',uvalue='INVERCONT',value=0,/ret,xsiz=5,event_func='save_info')
      settinver3=widget_base(invergroup,/row,xsiz=450,ysiz=40)
      invercent=cw_field(settinver3,title='Central wavelength of inversion: ',uvalue='INVERCENT',value='6173.34',/ret,xsiz=7,event_func='save_info')
      settinver4=widget_base(invergroup,/row,xsiz=450,ysiz=40)
      inverquanlow=cw_field(settinver4,title='Quantic numbers low level (S,L,J): ',uvalue='INVERQUANLOW',value='2,1,1',/ret,xsiz=8,event_func='save_info')
      settinver5=widget_base(invergroup,/row,xsiz=450,ysiz=40)
      inverquanup=cw_field(settinver5,title='Quantic numbers up level (S,L,J): ',uvalue='INVERQUANUP',value='2,2,0',/ret,xsiz=8,event_func='save_info')
      settinver6=widget_base(invergroup,/row,xsiz=450,ysiz=40)
      inversmall=cw_bgroup(settinver6,['Limit FOV for inversion: '],/nonexclusive,row=1,uvalue='INVERSMALL',uname='inversmall');,/return_index,event_func='save_info')
      settinver7=widget_base(invergroup,/row,xsiz=450,ysiz=40,sensitive=0,map=1,uname='fovinvergr')
      inver0=cw_field(settinver7,title='x0,y0 for FOV: ',uvalue='INVER0',value='100,100',/ret,xsiz=8,event_func='save_info')
      inverfov=cw_field(settinver7,title='FOV size in x,y for inversion: ',uvalue='INVERFOV',value='500,500',/ret,xsiz=8,event_func='save_info')


;-----------------------------------------------------------------
;compresssion parameters
      comprgroup=widget_base(sett,/col,uname='groupcompr',map=0)
      settcompr1=widget_base(comprgroup,/row,xsiz=450,ysiz=40)
      comprhw=cw_bgroup(settcompr1,['Simulate Hardware compression'],/nonexclusive,row=1,uvalue='COMPRHW',/return_index,event_func='save_info')
      settcompr2=widget_base(comprgroup,/row,xsiz=450,ysiz=40)
      comprratio=cw_field(settcompr2,title='Compression ratio: ',uvalue='COMPRRATIO',value='2.',/ret,xsiz=5,event_func='save_info')
      settcompr3=widget_base(comprgroup,/row,xsiz=450,ysiz=40)
      comprbpp=cw_field(settcompr3,title='Bits per pixel: ',uvalue='COMPRBPP',value=14,/ret,xsiz=5,event_func='save_info')
      comprbps=cw_field(settcompr3,title='Blocks per segment: ',uvalue='COMPRBPS',value=256,/ret,xsiz=5,event_func='save_info')
      settcompr4=widget_base(comprgroup,/row,xsiz=450,ysiz=40)
      comprskip=cw_bgroup(settcompr4,['Skip Bytes'],/nonexclusive,row=1,uvalue='COMPRSKIP',/return_index,event_func='save_info')
      comprmsb=cw_bgroup(settcompr4,['Most Significant Byte'],/nonexclusive,row=1,uvalue='COMPRMSB',/return_index,event_func='save_info')


;-----------------------------------------------------------------
;report parameters
      repogroup=widget_base(sett,/col,uname='grouprepo',map=0)
      settrepo1=widget_base(repogroup,/row,xsiz=450,ysiz=300)
      repomodules=cw_bgroup(settrepo1,['Input','Jittering','Polarization','Etalon','Optics','Pupil Apod.','FPA','Accumulation','Demodulation','Inversion','Compression'],/nonexclusive,row=11,uvalue='REPOMODULES',/return_index,space=.3,event_func='save_info')


;-----------------------------------------------------------------
      idmapped=globgroup
;create the structure for storing all the settings (except a few added
;when starting the simulation, in a routine above this). Stablish
;default values for the settings
      info={routines:intarr(12),progma:-1,startmed:0,nstart:0,talk:0,compress:1,saves:'../data/tests/',files:['sscene','jitscene','polscene','specscene','otfscene','pupapscene','fpascene','accuscene','demodscene','inverscene','comprscene'],$ 
            freport:'../reports/sophism_report.pdf',fsett:'../settings/settings_sophism.sav',settloadon:0,fsettl:'../settings/settings_sophism.sav',settloadasc:0,fsettasc:'../settings/sophism_ascii.pro',$
            scvel:0.,sclat:0.,scdis:0.28,temp:20.,acqsche:'Fast Polarization',ntim:1,telap:0.14,dhrew:322.,flen:4.125,discar:0,dualbeam:0,obsset:1,$ ;pdset:0,pddefoc:10., $
            datafilo:'../data/examples2/inverted_profs_139000.fits.gz',dataser:0,numfil:0,sscene_mag:1.,sscene_apod:0.,sz0:288,sz:288,sssamp0:20.8,sssamp:20.8,tsampor:1.,tsamp:0.020,wavepointor:992,wavepoint:0,wavesamp:14.127,wl:6173.341d0,geff:2.5,subset:'Subset',$ 
            jitoffset:[0.,0.],jitrms:0.5,jittype:'white noise',$
            filtype:'Hinode',filrms:1.,filcuth:0.1,filcut1:0.,filcut2:10.,filcutoff:10.,fillevelpre:'2.',filrangepre:'0.,10.',filsmooth:0,$
isson:0,$
            issatten:'kis',issunder:100.,issattf:'../data/',$
            transsys:'0.075',$
            otftype:'diffraction',defocopt:0,zerco:fltarr(10),zercoload:0,zercofile:'../data/hrew_zernikes.txt',frontnorm:0,frontrms:0.1, $ 
            pfwl:6173.34,pfhwb:3.,ncav1:2,fixed_pf:'fixed',pfload:0,pfcurve:'../data/prefilter.sav',pffov:0,pffovmax:-1600.,$
            etal_temp:20.,wavevolt:1,voltarrstr:'-500.,500.',wavearrstr:'-40,40',voltjit:0,voltjitrms:1.5,spectsizd:257.,spectsizdrms:0.,fabfin:200.,refl:0.925,reflrms:0.,frat:56.5,tilt:0.,discaretal:0,tvolt:1500.,$;sizdopt:0,sizdfil:'../data/RD_FE1-2final_ThicknessMap1.txt'
            intrange:3200.,phases:0,$
            pol_temp:20.,modtype:'vector. ideal',modnbuff:4,modspar:'V',modclip:0.0,mueller:'Identity',muellerrors:[0,0,0],muellerramp1:0.05,muellerrpha1:5.,muellerrang1:10.,muellerramp2:0.1,muellerrpha2:10.,muellerrang2:20.,retang1:0.,retang2:45.,retdeg1:[225.,225.,315.,315.],retdeg2:[235.,125.,55.,305.],reterrors:[0,0,0],reterramp1:0.05,reterrpha1:5.,reterrang1:10.,reterramp2:0.1,reterrpha2:10.,reterrang2:20.,tdeath1:5,tdeath2:94,tdeath3:94,tdeath4:80,cycrep:1,stokreorder:1,tobs:0.,$ ;,moddmod_type:'basic'
            birefsel:'None',birefloadname:'../data/biref_mat.sav',birefptv:0.05,birefscale:32.,birefvarmodel:0,birefloadmodel:0,birefloadmodelname:'../data/biref_model.sav',$
            fpa_temp:20.,fpaplate:0.5,etim:20.,fpashutt:'rolling',routtim:100.,fpaquan:0.5,fpawell:50.,fpaphotnoi:'Generate',fpaphotnoiname:'../data/fpascene_phot.sav',fparesamp:1,noitype:'normal',noisigma:70.,fpadeadpix:0,fpahotpix:0,fpasatval:65535.,fpafluxconv:0,fpafluxfactor:3e14,cosmic:0,cosmhit:16.8,fpadarksel:'Generate',fpadarkname:'../data/dark.save',fpagainsel:'Generate',fpagainname:'../data/gain.save',fpagain:0.01,fpaflatsel:'None',fpaflatmore:[0,0],fpaflatname:'../data/flat.save',fpaflrang:0.3,fpafring:0,fpafrinn:1.,fpafrinref:0.3,fpafrinsizd:200.,fpafrintheta:0.,$ ;,fpafringamp:1.,fpafringdir:0.,fpafringwidth:1.,fpafringtime:[0,0,0],$ ;,fpapixsiz:10.,0,$ 
            nacc:12,$
            adhoc:0,$
            inversigma:3,iter:15,invercont:0,invercent:6173.34,inverquanlow:[2,1,1],inverquanup:[2,2,0],inversmall:0,inver0:[100,100],inverfov:[500,500],$
            comprhw:0,comprratio:2.,comprbpp:14,comprbps:256,comprskip:0,comprmsb:0, $
            repomodules:intarr(11)$
}
;-----------------------------------------------------------------
;gui realization, settings structure store and event manager
      widget_control,widget_info(main,find_by_uname='voltarrstr'),sensitive=0
      for ee=1,2 do begin
         widget_control,widget_info(main,find_by_uname='reterramp'+strtrim(ee,2)),sensitive=0
         widget_control,widget_info(main,find_by_uname='reterrpha'+strtrim(ee,2)),sensitive=0
         widget_control,widget_info(main,find_by_uname='reterrang'+strtrim(ee,2)),sensitive=0
         widget_control,widget_info(main,find_by_uname='muellerramp'+strtrim(ee,2)),sensitive=0
         widget_control,widget_info(main,find_by_uname='muellerrpha'+strtrim(ee,2)),sensitive=0
         widget_control,widget_info(main,find_by_uname='muellerrang'+strtrim(ee,2)),sensitive=0
      endfor
      widget_control,main,/realize
      widget_control,main,set_uvalue=info
      xmanager,'sophism',main,/no_block

   end ;case0

;***********************************************************
;***********************************************************

   1: begin

;start the simulator directly without widgets, loading some settings file

;case of loading the settings through an ASCII file.
      if keyword_set(asciifile) then begin
;check if settings file provided or not
         if (n_elements(settfile) eq 0) then begin
            print,'No settings file provided'
            retall
         endif 

         if (settfile ne '../settings/sophism_ascii.pro') then spawn,'cp '+settfile+' ../settings/sophism_ascii.pro'
         sophism_loadascii

         restore,'../settings/settings_sophism.sav'

;start the simulation with the selected modules. 
         progctrl=[info.progctrl]
         for tt=0,n_elements(progctrl)-1 do call_procedure,progctrl(tt)
      endif
;---------------
;when loading previous settings file, consider new selection of
;modules and 'starting data'
      if keyword_set(idlfile) then begin
;check if settings file provided or not
         if (n_elements(settfile) eq 0) then begin
            print,'No settings file provided'
            retall
         endif 
         restore,settfile
;in case of producing the report without selecting any modules in its
;option tab, use the same modules selected for simulation
         if (info.routines(11) eq 1 AND max(info.repomodules) eq 0) then info.repomodules=info.routines(0:10) 

;save settings in default file (keep user's one)
         save,filena='../settings/settings_sophism.sav',info

         progctrl=['sophism_input','sophism_jitter','sophism_polmeas','sophism_linbo','sophism_otf','sophism_papo','sophism_fpa','sophism_accu','sophism_demod','sophism_inversion','sophism_compression','sophism_report']
         indprog=where(info.routines eq 1)

         numprogmas=-1
         if (max(info.progma) gt -1) then begin
            testprogma=file_search(info.saves+info.files(info.progma)+'_*.fits*',count=numprogmas)
;in case pupil apodization and loading linbo or otf, it won't
;find any fits but there should be a series, so force it. And read
;header from polarization files
            if (numprogmas eq 0) then begin
               if (info.routines(5) eq 1 AND (info.progma eq 3 OR info.progma eq 4)) then begin
                  info.dataser=1 
                  testprogma=file_search(info.saves+info.files(2)+'_*.fits*',count=numprogmas)
               endif else begin
                  print,'No files found'
                  return
               endelse  
            endif 
         endif 

;check if no input files are found
         if (info.numfil eq 0 AND max(info.progma) eq -1) then begin
            print,'No files found'
            return
         endif
         if (info.numfil gt 1 OR info.dataser eq 1 OR numprogmas gt 1) then begin
            info.ntim=round(info.tobs/info.tsamp)
            save,filena='../settings/settings_sophism.sav',info
         endif
;check if no modules are selected for the simulation
         if (indprog[0] eq -1) then begin
            print, 'No Modules selected'
            return
         endif else begin
;check if saves variable ends with / and add if not
            if (strmid(info.saves,strlen(info.saves)-1,1) ne '/') then info.saves=info.saves+'/'
;check for saves folder. Create it if not existing
            testsaves=file_search(info.saves,count=numsaves,/test_directory)
            if (numsaves eq 0) then spawn,'mkdir '+info.saves
;start the simulation with the selected modules.
            progctrl=progctrl(indprog)
            for tt=0,n_elements(progctrl)-1 do call_procedure,progctrl(tt)
         endelse
      
      endif 

;if class of settings file to be loaded was not specified,
;don't run and ask for it
      if NOT(keyword_set(asciifile)) AND NOT(keyword_set(idlfile)) then print,'No load file class specified. Run again specifying /asciifile or /idlfile and the settings filename to be loaded'


   end ;case1

endcase

END
