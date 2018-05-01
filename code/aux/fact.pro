function fact,numax

; La salida es una array DOUBLE de dimension numax+1 que contiene las
; factoriales de todos los numeros desde 0 hasta numax.
;		fac = fact(numax)

; Programa editado el 16 de Febrero de 1999
; Jose A. Bonet
;--------------------------------------------------------------------------


  fac=dblarr(numax+1)
  fac(0)=double(1)                 ;factorial de 0 es 1 by definition
  if numax eq 0 then return,fac
  fac(1)=double(1)
  for i=2,numax do begin
    fac(i)=fac(i-1)*double(i)
  endfor

  return,fac

end
