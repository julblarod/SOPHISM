FUNCTION Frac_shift, data, dx, dy
;+
; NAME:
;       FRAC_SHIFT
; PURPOSE:
;       Shift a dataset by a fractional amount using fourier techniques
; CALLING SEQUENCE:
;       Result = Frac_Shift(Data, dx [, dy])
; INPUTS:
;       Data  : 1D or 2D dataset
;       dx,dy : shift in x and y
; OUTPUTS:
;       Result: array of same dimension as Data, shifted by given amount
; RESTRICTIONS:
;       
; PROCEDURE:
;       Pad to 2^n dimensions if needed.  FFT.  compute and add linear
;       phase term.  Transform back.  Cut to original dimensions.  
; MODIFICATION HISTORY:
;       06-Dec-2001  P.Suetterlin, SIU
;       13-Dec-2001  Don't pad when size is already 2^n
;       05-Feb-2002  1D case: Add apodisation in the extended area, a
;                    sharp drop at the rim introduces errors.
;                    Change size computation of the extention.
;-

IF n_params() LT 2 THEN $
  message, 'Use:  RESULT = FRAC_SHIFT( DATA, DX [,DY] )'
IF n_params() LT 3 THEN dy = 0.

result = data
sd = size(data)

IF sd(0) EQ 1 THEN BEGIN
      ;;; extend the array by 20% for apodisation
    n = fix(sd(1)*1.2)
      ;;; extend to a nice number
    n_1 = 2^fix(alog10(n)/0.30103)
    IF n_1 LT n THEN $
      n_2=2^(fix(alog10(n-n_1)/0.30103+1) > 2) $
    ELSE $
      n_2 = 0
    n = n_1+n_2
      ;;; generate an apodisation function 0->1 in the extended area
    n_dif = n-sd(1)
    ap = .5-cos((findgen(n_dif)+1)/(n_dif+1)*!pi)/2
    
    pd = replicate(data(0), n)
    pd(0) = data
      ;;; insert the apodisation
    pd(sd(1)) = data(sd(1)-1) + (data(0)-data(sd(1)-1)) * ap
      ;;; FFT, compute and add phase term
    fpd = fft(pd, -1)
    n2 = n/2
    k = shift(-(findgen(n)-n2), n2)
    d = 2*!pi*K/n*dx
    d1 = complex(cos(d), sin(d))
    pd(*) = fft(fpd*d1, 1)
    result(*) = pd(0:sd(1)-1)
    return, result
ENDIF ELSE BEGIN
    nx = 2^fix(alog10(sd(1))/0.30103)
    IF nx LT sd(1) THEN nx = 2*nx
    ny = 2^fix(alog10(sd(2))/0.30103)
    IF ny LT sd(2) THEN ny = 2*ny
    nx2 = nx/2
    ny2 = ny/2
    
    pd = replicate(data(0), nx, ny)
    pd(*) = avg(data)
    pd(0:sd(1)-1, 0:sd(2)-1) = data
    fpd = fft(pd, -1)
    kx = shift(-(findgen(nx)-nx2), nx2)
    kx = rebin(kx, nx, ny)
    ky = shift(-(findgen(ny)-ny2), ny2)
    ky = rebin(transpose(ky), nx, ny)
    
    d = 2*!pi*Kx/nx*dx + 2*!pi*Ky/ny*dy
    d1 = complex(cos(d), sin(d))
    pd(*) = fft(fpd*d1, 1)
    result(*) = pd(0:sd(1)-1, 0:sd(2)-1)
    return, result
ENDELSE

END
