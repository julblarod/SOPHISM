PRO sophism_input_redim,sscene

;=================================================================
; DESCRIPTION
; Replicates datacubes in spatial dimensions and resamples to detector
; pixel size.
; Has to be done all together to save computation time and disk space.
; If resampling to detector is enabled, it will change the value of
; spatial sampling in the settings file to the new one (i.e. fpaplate).
;
; CALLED FROM 
;   sophism_sscene
;
; SUBROUTINES
;
;
; MAIN INPUTS
;   sscene_mag     Replication factor, the number of times to
;                     replicate the data in both x and y
;   sssamp         Spatial sampling, from original sampling of the
;                     data and distance S/C-Sun [arcsec]
;   fpaplate       Detector pixel size in [arcsec]
;
; OUTPUT
;   Final input data prepared for simulation run
;
; VERSION HISTORY
;   J. Blanco. Sept 2012. v1_0. Split and re-organized from other
;      subroutines
;   J. Blanco. Dec 2012. v1_1. Added printings for talk case
;   J. Blanco. Jan 2013. v1_2. Correction for reading uncompressed
;      data
;   J. Blanco. Jan 2013. v2_0. Reordering of resample and
;      replication. Modified way of procuder: not reading or writing
;      now, just used as a step in global routine.
;   J. Blanco. Apr 2013. v2_1. Moved verbose line to main routine
;
;=================================================================

restore,'../settings/settings_sophism.sav'

if (info.fparesamp eq 1) then begin
; calculate the resize ratio (rounding the new pixel number, since
; ratio won't give fix result)
   sizy=size(sscene)
   sssamp=(info.sssamp0*1e3)*(1./(info.scdis*1.495e11))*!radeg*3600.
   ccdrat=sssamp/info.fpaplate
   dimnew=fix(round(sizy[3]*ccdrat))
   imor=sscene
   sscene=fltarr(sizy[1],sizy[2],dimnew,dimnew)
   for lam=0,sizy[1]-1 do begin
      for pol=0,sizy[2]-1 do sscene[lam,pol,*,*]=congrid(reform(imor[lam,pol,*,*]),dimnew,dimnew,/int,cub=-0.5,/center)
   endfor
endif ;resamp
   
imor=sscene
sizy=size(sscene)
sscene=fltarr(sizy[1],sizy[2],sizy[3]*info.sscene_mag,sizy[4]*info.sscene_mag)
; replicate spatially the data by concatenating the same cube in x and y
; (considers periodic border condition; MHD data prepared like that)
for ix=0,info.sscene_mag-1 do begin
   for iy=0,info.sscene_mag-1 do begin
      sscene[*,*,ix*sizy(3):(ix+1)*sizy(3)-1,iy*sizy(4):(iy+1)*sizy(4)-1] = imor(*,*,*,*)
   endfor
endfor



END
