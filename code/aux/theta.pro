function theta,r

; Creates a 2r x 2r matrix with the angular coordinates (in rad) of each pixel
; as measured from the origin in pixel (r,r).

; Programa editado el 16 de Febrero de 1999
; Jose A. Bonet

;--------------------------------------------------------------------------


  d=2*r                              ; d x d es la dimension final de la imagen
  x=(findgen(d)-r) # transpose(replicate(1.,d))   ;x-coordinates of the pixels
  y=transpose(x)   				  ;y-coordinates of the pixels

  w=where(abs(x) lt 1.e-35,num) ;+/- 1.e-38 = floating underflow in the division
  if num ne 0 then x(w)=0
  th=atan(y,x)

  return,th

end

