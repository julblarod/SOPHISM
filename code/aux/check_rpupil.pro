function check_rpupil,nn,rpupil

; Check whether the size of the pupil array is large enough to include the 
; full aperture of radius "rpupil". If not, the program returns a new size
; for the pupil array. This is the case for sub-critically sampled data (aliasing).

;INPUTS:
;     nn     = size of the pupil array (i.e. image-patch-size/2)
;     rpupil = radius of the aperture [in pixel], as computed from the
;              diffraction cutoff, image-patch-size and sampling interval.

;OUTPUT:
;     mm = new size (even number) for the pupil array so that the entire
; aperture fits into it. If 2.*rpupil </= (nn-2) then the dimension "nn" is
; sufficiently large and the output "mm" = "nn"
;
;              mm = check_rpupil(nn,rpupil)
;____________________________________________________________
; 19-Marzo-2011
; Jose A. Bonet
;____________________________________________________________

mm=nn
dif=2.*rpupil-(nn-2) ;dif must be </= 0 otherwise nn must be increased 

if dif gt 0 then begin
  if dif mod 2 eq 0 then begin
    mm=nn+dif
  endif else begin
    mm=nn+2*(fix(dif/2.)+1)
  endelse
endif

return,mm  

end
