;+
; NAME:
;	FFT_SHIFT
; PURPOSE:
;	Shift 1D or 2D data by FFT
; CATEGORY:
;	Image Processing
; CALLING SEQUENCE:
;	y=fft_shift(x,shift)
; INPUTS:
;	x     : Input array, may be 1D or 2D
;       shift : a shift scalar, if x is 1D or a two element vector,
;		if x is 2D
; OPTIONAL INPUT PARAMETERS:
;	-
; KEYWORDS:
;	-
; OUTPUTS:
;	y     : Data shifted by 'shift'
; COMMON BLOCKS:
;	-
; SIDE EFFECTS:
;	Non periodic data will not be shifted correctly. 
; RESTRICTIONS:
;	to 1D or 2D arrays. Data should be periodic.
;       number of elements in each dimension must be even
; NOTES:
; MODIFICATION HISTORY:
;	written Feb.92 by Reinhold Kroll
;	Last Update 24.Feb.92 Reinhold Kroll
;-
;
function sign,x 
; defines the sign of a vector
return,(-1*(x lt 0)) + (x ge 0)	; assign -1 if lt, +1 if ge  0 
end
;
;
function phaser,n
; return the phase vector of an	fft of a real function. 
; looks like /\, with the second half complex conjugate 
; n is the number of elements 
i=complex(0,1)			; define i
m=n/2				; one half of the number of elements
in=m-findgen(n)			; thats m,m-1,....-m+2,-m+1
in=(m - abs(in)) * sign(in)	; thats 0,1,...m,-(m-1),...-1
return,-2.*!pi*i*in/n		; and thats the complex shift
end
;
;  THE REAL THING STARTS HERE
;
function fft_shift, x,shft
on_error,2		; return to caller on error
sx=size(x)
nx=sx(0)		; dimension of x
ns=n_elements(shft)	; elements in shift should be = dimension of x
;
; check sanity of data
;
if (nx gt 2) or (nx lt 1)  then begin
	print,'**** Can only process 1D or 2D data ! '
	print,format='("**** Your Data are ",I2,"D !")',nx
	return,-2
	endif
if nx ne ns then begin	; dimensions should be the same
	print,'**** Input dimension mismatch ! '
	print, $
         format='("**** Data are ",I2,"D, shift vector is ",I2,"D !")', $
         nx,ns
	return,-1
	endif
if (nx eq 1) and (sx(1) mod 2) then begin
	print,'**** Data vector must have even elements ! '
	print,'**** Your Data has ',sx(1),' elements'
	return,-3
	endif
if (nx eq 2) and ( (sx(1) mod 2) or (sx(2) mod 2) ) then begin
	print,'**** Data must have even elements in each dimension! '
	print,'**** Your Data has ',sx(1),' x ',sx(2),' elements'
	return,-4
	endif
;
; lets do the work!
;
;
; 1D case
;
if nx eq 1 then begin
	nel=n_elements(x)		; # of elements in x
	sh=shft*phaser(nel)		; the shift vector!
	xf=fft(x,-1)			; transform
	xf1=xf*exp(sh)			; shift in Fourier space
	y=float(fft(xf1,1))		; transform back
	return,y			; return result
	endif 
;
; 2D case
;
if nx eq 2 then begin
	nelx=sx(1)			; size in x
	nely=sx(2)			; size in y
	shx=exp(shft(0)*phaser(nelx))	; shift vector in x
	shy=exp(shft(1)*phaser(nely))	; shift vector in y
	xf=fft(x,-1)			; transform
	for i=0,nely-1 do xf(*,i)=xf(*,i) * shx	; shift in x
	for i=0,nelx-1 do xf(i,*)=xf(i,*) * shy	; shift in y
	y=float(fft(xf,1))		; transform back
	return,y
	endif


end
