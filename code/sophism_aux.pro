;+
; ==============================================================================
;
; DESCRIPTION
;   Auxiliary routines 
;
; USAGE INSTRUCTIONS
;
; DEPENDENCIES
;
; VERSION HISTORY 
;   Alex Feller   v0.1    2010-09-21
;   J. Blanco. 2010. v1_0. Fit into sophism
;   B. Loptien. Mar 2013. v1_1. Changed indgen to lindgen in
;      sophism_make_freq for because of possible long arrays
;
; ==============================================================================
;-

; ==============================================================================
; Generate frequency array compatible with IDL FFT convention
; ==============================================================================

function sophism_make_freq, nsamp, fsamp

if nsamp mod 2 eq 0 then begin
   nfreq = nsamp/2
   pfreq = lindgen(nfreq+1) * fsamp
   nfreq = -reverse(lindgen(nfreq-1) + 1) * fsamp
endif else begin
   nfreq = (nsamp-1)/2
   pfreq = lindgen(nfreq+1) * fsamp
   nfreq = -reverse(lindgen(nfreq)) * fsamp
endelse

return, [pfreq, nfreq]

end

; ==============================================================================
; Resample image array using FFT resampling and rebin 
; This method emulates sampling by a detector array with adjacent square 
; pixels
; ==============================================================================

function sophism_resample, pim, psamp1, psamp2

; Copy parameters
im = pim
samp1 = psamp1
samp2 = psamp2 

; Round pixel sampling values to 3 digits
a = 10.^(floor(alog10(samp1))-2)
samp1 = round(samp1/a) * a
a = 10.^(floor(alog10(samp2))-2)
samp2 = round(samp2/a) * a

sfac = samp2 / samp1

if sfac gt 1 then begin

   sz = (size(im))[1]
   rfac = ceil(sfac)

; Enlarge array to integer factor of final array dimensions

; TODO: replace bilinear interpolation by FFT resampling. With the current
; implementation, after enlarging, the difference in stddev(im) is in the 3rd
; significant digit, the difference in mean(im) is in the 6th significant
; digit.  

   efac = rfac / sfac
   esz = round(sz * efac)
   im = congrid(im, esz, esz, /interp)

; Rebin array to samp2
   fsz = floor(sz / sfac) * rfac
   if fsz lt esz then im = im[0:fsz-1, 0:fsz-1]
   im = rebin(im, fsz/rfac, fsz/rfac)

endif

if sfac lt 1 then begin
   message, 'Case sfac < 1 not implemented yet!'
endif

return, im

end

; ==============================================================================
; Apodize image array
; Adadpted from PDjab/apo2d.pro
; ==============================================================================

pro sophism_apo2d, im, perc

;  apodization
s=size(im)
edge=100./perc
sd=stdev(im,av)
im=im-av
xmask=fltarr(s(1))+1.
ymask=fltarr(s(2))+1.
smooth_x=fix(round(s(1)/edge))        ; width of the edge in x-dimension
smooth_y=fix(round(s(2)/edge))        ; width of the edge in y-dimension


;  smoothing with a cosine
for i=0,smooth_x do xmask(i)=(1.-cos(!pi*float(i)/float(smooth_x)))/2.
for i=0,smooth_y do ymask(i)=(1.-cos(!pi*float(i)/float(smooth_y)))/2.
xmask(s(1)-smooth_x:s(1)-1)=reverse(xmask(1:smooth_x))
ymask(s(2)-smooth_y:s(2)-1)=reverse(ymask(1:smooth_y))

; Apply smoothing mask 
for i=0,s(1)-1 do im(i,*,*)=im(i,*,*)*xmask(i)
for i=0,s(2)-1 do im(*,i,*)=im(*,i,*)*ymask(i)    
im=im+av

end

; ==============================================================================
; Repeat array
; Used to define modulation/demodulation functions
; ==============================================================================

function sophism_repeat, a, n

b = a
for i=1,n-1 do begin
   b = [b, a]
endfor

return, b

end
