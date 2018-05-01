function zernike_mn,jmax

; La salida es un array de dimensiones 2 x jmax que contiene las parejas
; de coef. [m,n] para cada valor de j desde 1 hasta jmax.
;		mn = zernike_mn(jmax)
; !!! Tener en cuenta que los valores de [m,n] correspondientes a un deter-
; minado valor de j son los elementos mn(0,j-1), y mn(1,j-1) del array !!!

; Programa editado el 16 de Febrero de 1999
; Jose A. Bonet (transcripcion a idl del programa zernike_mn_f.ana en ana
;                de M.L.)
;--------------------------------------------------------------------------

  m=0
  n=0
  i=0
  REPEAT BEGIN
    i=i+1
    d = intarr(i+1)
    IF Odd(i) THEN BEGIN
      FOR ii=0,i/2 DO BEGIN
        d(2*ii) = 2*ii+1
        d(2*ii+1) = 2*ii+1
      ENDFOR
    ENDIF ELSE BEGIN
      FOR ii=1,i/2 DO BEGIN
        d(2*ii-1) = 2*ii
        d(2*ii) = 2*ii
      ENDFOR
    ENDELSE
    m = [m,d]
    n = [n,replicate(i,i+1)]
  END UNTIL (size(n))(1) GE jmax

  mn=transpose([[m],[n]])

  Return, mn

end
