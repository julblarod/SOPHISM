PRO sophism_linbo_lambint,wavearr,index,wg,waves,wavesind

; ==================================================================
; DESCRIPTION
;    Determine if input wavelengths are sampled in input data or
;    neighbouring weighting must be used
;
; CALLED FROM 
;    sophism_linbo
;
; SUBROUTINES
;
; MAIN INPUTS
;    wavearr       array of selected wavelengths (AA)
;
; OUTPUT
;    index         integers to identify belonging of sampled
;                     wavelengths to selected wavelengths for future
;                     calculations 
;    wg            weights for measuring contribution of each
;                     wavelength in a position 
;    waves         wavelengths sampled in input data to work with
;    wavesind      indexes of wavelengths sampled to be used
;
; VERSION HISTORY
;    J. Blanco. Nov 2012. v1_0. Split from main sophism_linbo routine.
;
; ==================================================================

progind=3
;loading settings file
restore,'../settings/settings_sophism.sav'

lll=info.lll
nwave=info.etalpos
;initialize variables
index=0.
wg=0.
waves=0.
wavesind=0.
for ww=0,nwave-1 do begin
;check if there is a sampled wavelength coinciding with the selected ones
   diflam=abs(lll-wavearr[ww])
;if there is not, take the two closest ones and store their
;wavelengths and index positions.
   if (min(diflam) gt 1e-3) then begin
;index just to mark wavelengths belonging to selected positions
      index=[index,ww,ww]
      indy=sort(diflam)
      indy=indy[0:1] 
      indy=indy[sort(indy)]
;calculate the weight to give to each sampled position when adding
;them afterwards to produced the selected one
      wgpre=1-abs(wavearr[ww]-lll(indy[0]))/(info.wavesamp*1e-3)
      wg=[wg,wgpre,1.-wgpre]
      waves=[waves,lll(indy[0:1])]
      wavesind=[wavesind,indy]
   endif else begin
;case the selected position is sampled. There is no need to do the
;above two-wavelengths process.
      index=[index,ww]
      indy=sort(diflam)
      indy=indy[0]
      wg=[wg,1.]
      waves=[waves,wavearr[ww]]
      wavesind=[wavesind,indy]
   endelse 
endfor    ;ww1
;remove the first element (the 0 of variable declaration) from arrays
if (n_elements(index) gt 1) then begin
   index=index(1:*)
   wg=wg(1:*)
   waves=waves(1:*)
   wavesind=wavesind(1:*)
endif 




end
