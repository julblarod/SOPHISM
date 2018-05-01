function radius_aper2,r,rc,ap

; Generate two masks, ap and rd, of dimensions 2r x 2r :

;           rd = radius_aper2(r,rc,ap)      (where rc < r)

; rc can be either integer or fractional. This represents an improvement with
; respect to the routine radius_aper.pro that only works with integer values.

; ap = mask defining a circular aperture of radius rc, centered at pixel (r,r). ;      The pixel value is 1 within the circle and 0 outside. Each pixel
;      intersected by rc gets a value that corresponds to the fractional
;      area of the inner portion of the pixel.

; rd = mask defining the radial coordinate of the pixels in the aperture,
;      normalized to rc. Each pixel intersected by rc gets a value that
;      corresponds to the normalized distance from the center to the
;      centroid of the inner portion of the pixel.

; (See Apendix H)

; Programa editado el 6 de Mayo de 1999 como una extension de radius_aper.pro
; para el caso de rc = numero fraccionario o entero.
; Jose A. Bonet

;-----------------------------------------------------------------------------

; Defines the first quadrant in both masks (first approach)

  d = 2*r                          ; d x d is the final dimension of the masks.
  rd = dist(d,d)
  rd1 = rd(0:r-1,0:r-1)/rc
  ap1 = rd1 le 1.  &  ap1=float(ap1)
  rd1 = rd1*ap1

; Define indices (i,j) of those pixels intersected by rc (first octant).

;  ii=nint(rc)
  ii=fix(round(rc))
  i = lonarr(2*ii)
  j = i
  jant = 0
  npix = 0

start:

  xx = ii-0.5
  jult = fix(sqrt(rc*rc-xx*xx)+0.5)
  for k=jant,jult do begin
    i(npix) = ii
    j(npix) = k
    if (i(npix) eq j(npix)) then goto,fin
    npix = npix+1
  endfor
  ii = ii-1
  jant = jult
  goto,start

fin:

  i = i(0:npix)
  j = j(0:npix)

;for h=0,npix do print,i(h),j(h)

; Calculate de centroid and area of the inner portion of those pixels
; intersected by rc (first octant).

  for k=0,npix do begin       ; Absolute coord.system at center of pixel (0,0)

    x0 = i(k)-0.5            ; absolute coords.of the lower left corner of pixel
    y0 = j(k)-0.5            ; (i,j) ==> origin of the new coord. system OXY.
    x1 = sqrt(rc*rc-y0*y0)-x0                ; (x1,y1), (x2,y2) coord. referred
    y1 = 0                                   ; to OXY, of the intersections of
    x2 = 0                                   ; the circunference with the axis
    y2 = sqrt(rc*rc-x0*x0)-y0                ; OX and OY.

    case 1 of
      (x1 gt 0.) and (x1 le 1.) and (y2 le 1.): begin
        xm = x1/3.  &  ym = y2/3.    ; Centroid coord. referred to OXY.
        area = x1*y2/2.
        xxm = x0+xm  &  yym = y0+ym  ; Centroid coord. referred to absolute axis
        end
      (x1 gt 0.) and (x1 le 1.) and (y2 gt 1.): begin
        yup = y0+1.
        x2 = sqrt(rc*rc-yup*yup)-x0  &  y2 = 1.
        aa = x2  &  ab = (x1-x2)/2.  &  area = aa+ab
        xa = x2/2.  &  ya = y2/2.
        xb = x2+(x1-x2)/3.  &  yb = 1./3.
        xm = (aa*xa+ab*xb)/area     ; Centroid coord. referred to OXY.
        ym = (aa*ya+ab*yb)/area
        xxm = x0+xm  & yym = y0+ym  ; Centroid coord. referred to absolute axis
        end
      (x1 gt 1.): begin
        xup = x0+1.  &  yup = y0+1.
        x1 = 1.  &  y1 = sqrt(rc*rc-xup*xup)-y0
        x2 = sqrt(rc*rc-yup*yup)-x0  &  y2=1.
        a = 1.  &  ab = (1.-y1)*(1.-x2)/2.  &  area = a-ab
        x = 0.5  &  y=0.5
        xb = x2+2.*(1.-x2)/3.  &  yb = y1+2.*(1.-y1)/3.
        xxm = x0+(a*x-ab*xb)/area   ; Centroid coord. referred to absolute axis
        yym = y0+(a*y-ab*yb)/area
        end
      (x1 le 0.): begin
        area = 0.
        xxm = 0.  &  yym = 0.
        end
    endcase

; Define distance and aperture in those pixels intersect.by rc (first quadrant)

    rd1(i(k),j(k)) = sqrt(xxm*xxm+yym*yym)/rc
    rd1(j(k),i(k)) = rd1(i(k),j(k))
    ap1(i(k),j(k)) = area
    ap1(j(k),i(k)) = area

  endfor

; Centering of the masks in the center of the matrix (pixel (r,r)) and generate
; the other three quadrants

  rd = rd*0. & ap = rd
  rd(r,r) = rd1
  rd(1,r) = reverse(rd(r+1:*,r:*),1)
  rd(1,1) = reverse(rd(1:*,r+1:*),2)
  ap(r,r) = ap1
  ap(1,r) = reverse(ap(r+1:*,r:*),1)
  ap(1,1) = reverse(ap(1:*,r+1:*),2)


;openw,1,'resultados'
;  for m=0,npix do printf,1,i(m),j(m)
;  for m=0,d-1 do begin
;    printf,1,rd(*,m),format='(1x,10F7.4)'
;  endfor
;  printf,1,''
;  for m=0,d-1 do begin
;    printf,1,ap(*,m),format='(1x,10F7.4)'
;  endfor
;close,1


return,rd

end
