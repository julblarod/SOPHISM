;+
; NAME:
;       RESIZE
; PURPOSE:
;       Resize 1- or 2-dimensional arrays.
;       Much like the built-in function REBIN, but the restriction to 
;       integral multiples or fractions of the original size is dropped.
; CATEGORY:
;       Image processing.
; CALLING SEQUENCE:
;       a=resize(b,nx,ny)     ; for 2-dim arrays
;       a=resize(b,n)         ; for 1-dim arrays
; INPUTS:
;       b = data array to be resized
;       nx,ny = new dimensions
; KEYWORD PARAMETERS:
;       /sample : If set, nearest neighborsampling is used (like with REBIN).
; OUTPUTS:
;       function value = resized data array
; COMMON BLOCKS:
;       none
; SIDE EFFECTS:
;       none
; RESTRICTIONS:
;       Only for one or two dimensional arrays.
; PROCEDURE:
;       The built-in function POLY_2D is used for (bi-)linear interpolation.
;       If the new dimension is smaller than the old one, and /sample is
;       not set, the array is first enlarged to the smallest possible
;       integral multiple of the new dimension, and then averaged using REBIN.
; MODIFICATION HISTORY:
;    Author : A. Welz,  Uni. Wuerzburg, Germany , Feb. 1991
;-
function resize,dati,nxnew,nynew,sample=sample
;
on_error,2
;
if n_elements(nxnew) eq 0 then nxnew=1
if n_elements(nynew) eq 0 then nynew=1
nxnew=fix(nxnew)>1 & nynew=fix(nynew)>1
;
data=reform(dati)
s=size(data)
if s(0) eq 1 then begin
    data=reform(data,s(1),1)
    s=size(data)
endif
if s(0) ne 2 then goto,errex
;
type=1. < abs(data(0,0)) > 1.
nxold=s(1)
nyold=s(2)
;
nx=nxnew
ny=nynew
if not keyword_set(sample) then begin
   if nxnew lt nxold then nx=nxnew*long(nxold/nxnew+1)
   if nynew lt nyold then ny=nynew*long(nyold/nynew+1)
endif
;
   d=fltarr(nxold+1,nyold+1)*type
   d(0:nxold-1,0:nyold-1)=data
   data=d
   data(nxold,0:nyold-1)=data(nxold-1,0:nyold-1)
   data(*,nyold)=data(*,nyold-1)
;
xfac=float(nxold)/float(nx)
yfac=float(nyold)/float(ny)
p=fltarr(2,2) & q=p
p(0,1)=xfac
q(1,0)=yfac
;
if keyword_set(sample) then begin
   return,reform(poly_2d(data,p,q,0,nx,ny))
endif else begin
   return,reform(rebin(poly_2d(data,p,q,1,nx,ny),nxnew,nynew))
endelse
;
errex: print,' RESIZE: dimension must be 1 or 2'
return,dati
end
