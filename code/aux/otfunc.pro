function otfunc,h,norm

; Calcula la Optical Transfer Function.
; Inputs:
;	   h = generalized pupil complex function of dimensions n x n.
; Outputs:
;	   norm = normalization factor for the OTF
;	   otf  = otfunc(h,norm)
; otf = funcion compleja de dimensiones 2n x 2n, centered at pixel (n,n)
;
; Programa editado el 15 de Febrero de 1999.
; Jose A. Bonet
;
;-------------------------------------------------------------------------

  n=(size(h))(1)

; Calculo de otf = cte.*autocorrelation (generalized pupil function)

  otf=crosscorr_c(h,h,/doble)    ; aqui crosscorr_c ha de trabajar en Double
		                 ; Precission. Si no, la otf que resulta no
				 ; es perfectamente hermitica y esto redunda
				 ; en Q y Q^2 --> ruido al restaurar.

; Normaliza al valor en el origen (correspondiente a shift de h = 0)

  norm=float(otf(0,0))
  otf=otf/norm

; Desplaza para que el maximo de la otf este en el punto de coord. (n,n)

  otf=shift(otf,n,n)

  return,otf

end
