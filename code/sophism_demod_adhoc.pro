pro sophism_demod_adhoc,ima,apod,contpos

; ==================================================================
; DESCRIPTION
;   Corrects crosstalk from Stokes V to Q and U by a linear fitting.
;
; CALLED FROM 
;   sophism_demod
;
; SUBROUTINES
;   
;
; MAIN INPUTS
;    ima       demodulated data with wavelengths and Stokes paramenters
;    apod      percentage of image size to apodize
;
; OUTPUT
;    ima       xtalk corrected data
;
; VERSION HISTORY
;   J. Blanco. Nov 2012. v1_0. Derived after IMaX reduction programs
;   J. Blanco. Jan 2015. v1_1. Included the invercont variable for the
;   continuum position selected by user.
;   J. Blanco. Feb 2016. v1_2. Corrected bug with continuum position
;
; ==================================================================

;program index in SOPHISM code
progind=8

sz=size(ima)

;apodisation border, in case it was selected
apodx=0
apody=0
if (apod gt 0) then begin
   edge=100./apod
   apodx=fix(round(sz(1)/edge))
   apody=fix(round(sz(2)/edge))
endif

;remove Stokes I xtalk obtained from position at continuum.
;considering the continuum is at first lambda position. Too
;restrictive? Could be generalized?
for ll=0,sz(1)-1 do begin
	for pol=1,sz(2)-1 do ima[ll,pol,*,*]=ima[ll,pol,*,*]-mean(ima[contpos,pol,*,*])/mean(ima[contpos,0,*,*])*ima[ll,0,*,*]
endfor

;linear fits between V and Q, and V and U, for every lambda and correction
for ll=0,sz(1)-1 do begin
	co=poly_fit(ima[ll,3,apodx:sz(3)-apodx-1,apody:sz(4)-apody-1],ima[ll,1,apodx:sz(3)-apodx-1,apody:sz(4)-apody-1],1)
	ima[ll,1,*,*]=ima[ll,1,*,*]-co[1]*ima[ll,3,*,*]
	co=poly_fit(ima[ll,3,apodx:sz(3)-apodx-1,apody:sz(4)-apody-1],ima[ll,2,apodx:sz(3)-apodx-1,apody:sz(4)-apody-1],1)
	ima[ll,2,*,*]=ima[ll,2,*,*]-co[1]*ima[ll,3,*,*]
endfor


end
