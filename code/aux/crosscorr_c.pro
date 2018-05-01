function crosscorr_c,f,g,doble=dd

; Calcula la cross-correlation de las funciones 2D complejas, f y g, cada una
; de dimensiones n x n. Utiliza el Tma. de cross-correlation de la transformada
; discreta de Fourier.

;    /double  : if set, all the calculations are performed in Double Precission
;	        although the output will be converted to Single Precission.

; Cuando se trabaja en Double Precision sale al final del programa un warning
; diciendo: "Program caused arithmetic error: Floating underflow"

; El resultado:
;		cross=crosscorr_c(f,g,doble=dd)
; es una funcion compleja (Single Precission) de dimensiones 2n x 2n
;
; Programa editado el 15 de Febrero de 1999.
; Jose A. Bonet
;
;-------------------------------------------------------------------------

; Prolonga las funciones con ceros hasta dimension 2n x 2n

  n2=2*(size(f))(1)
  fex=complexarr(n2,n2)
  if keyword_set(dd) then fex=dcomplexarr(n2,n2)
  gex=fex
  fex(0,0)=f & gex(0,0)=g

; Aplica el Tma. de correlacion

  fex=fft(fex,-1) & gex=fft(gex,-1)
  cross=fft(n2*n2*conj(fex)*gex,1)
  cross=complex(cross)         ; returns a Single Precission function

return,cross

end
