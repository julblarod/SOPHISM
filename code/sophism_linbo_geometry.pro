PRO sophism_linbo_geometry,tam,rpup,llo,n_etalon,incangle_internal,mask

; ==================================================================
; DESCRIPTION
;    Calculates the cone of incident angles over the etalon in
;    telecentric mounting by means of the scalar product of the
;    etalon's directional vector and incident rays vectors.
;    For use in the case of pupil apodisation.
;
; CALLED FROM 
;    sophism_linbo
;
; SUBROUTINES
;    radius_aper2
;
; MAIN INPUTS
;    tam             Size of the pupil , extended if necessary due to
;                       aliasing [pix]
;    rpup            Radius of pupil [pix]
;    llo             Wavelength scanning positions (just in case of
;                       including a wavelength-dependant rpup; not yet)
;    n_etalon        Refractive index of the etalon (array)
;    frat            F-ratio of the instrument
;    tilt            Inclination of the etalon with... [deg]
;
; OUTPUT
;    3D Array (2D spatial + 1D wavelength) of internal incident angles
;    and binary mask
;
;
; VERSION HISTORY
;    J. Hirzberger
;    J. Blanco. April 2012. v1_0. Fit into sophism
;    J. Blanco. September 2012. v2_0. Added pupil definition by radius_aper2
;       and some changes to variables.
;
; ==================================================================



progind=3

restore,'../settings/settings_sophism.sav'

;define pupil aperture and binary mask
rd = radius_aper2(tam/2.,rpup,mask)
;piensa si tiene sentido dejarlo asi o rd deberia ser 1 en toda la
;periferia. Por pensar...

;--------------------------------------------------------------------- 
; ray incidences on etalon
;--------------------------------------------------------------------- 

; maximum external incident ray angle
alpha = 1./(2.d0*info.frat)                      

; 2D array of telecentric incidence angles (assuming no tilt) 
d=rd*alpha



;--------------------------------------------------------------------- 
; ray incidences on tilted etalon
;--------------------------------------------------------------------- 

; etalon tilt angle in radians
tilt = info.tilt*!dtor

;define the etalon tilt in space. z is optical axis, x-y is plane of
;the etalon. The etalon can only tilt in one axis: x axis
x = 0.
z = 1.
y = z*tan(tilt)

; etalon surface vector
v_etalon = double([x,y,z])
;vector magnitude
absetalon = sqrt(abs(transpose(v_etalon)#v_etalon))

;--------------------------------------------------------------------- 
; Initialize parameter settings
;--------------------------------------------------------------------- 

; array of incidence angles
incangle = dblarr(tam,tam)

FOR i=0,tam-1 DO BEGIN
   FOR j=0,tam-1 DO BEGIN
;define the geometry of the rays from pupil to etalon, with z=1
        th = d[i,j]
        c = tan(th)
        beta = atan(double(j-tam/2),float(i-tam/2))
        x = c*cos(beta)
        y = c*sin(beta)
; ray vector
        ray_vec = [x,y,z]
        absray = sqrt(abs(transpose(ray_vec)#ray_vec))
; scalar poduct
        scal_prod = total(v_etalon*ray_vec)
; external incidence angle
        incangle[i,j] = acos(scal_prod/(absetalon*absray))

   ENDFOR ;j
ENDFOR ;i

;the refractive index is function of wavelength, so create and array
;for internal angles of 2D spatial + 1D wavelength
nnsiz=n_elements(n_etalon)
incangle_internal=fltarr(tam,tam,nnsiz)
;internal incidence angle (Snellius)
for nn=0,nnsiz-1 do begin
   incangle_internal[*,*,nn] = asin(sin(double(incangle))/n_etalon[nn])
   incangle_internal[*,*,nn]=incangle_internal[*,*,nn]*mask
endfor

save,incangle,incangle_internal,filename=info.saves+info.files(progind)+'_geomet.sav'


END
