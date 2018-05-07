
pro sophism_polmeas_birefringence_model,CM,MMatrix,MMatrix_ideal,RMatrixB,ds_PtV=ds_PtV;,doplot=doplot,steps=steps,scale=scale,label=label,model=model,ds_PtV=ds_PtV

; ==============================================================================
;
; DESCRIPTION
;    Calculate the birefringence matrixes from models of spatial
;    variation of the EW to be used in sophism_polmeas_birefringence
;
; CALLED FROM 
;    sophism_polmeas_birefringence
;
; SUBROUTINES
;
;
; MAIN INPUTS
;    scale        Reduction factor scale of the data size to speed up
;                    calculations 
;    ds_ptv       Peak to valley of the varition model of the HREW
;    model        Model of the HREW retardance
;
;
; OUTPUT
;    Muodulation matrixes with the birefringence, without it,
;    birefringence matrix and birefringence crosstalk
;
; VERSION HISTORY 
;   D. Orozco   v0.1   
;   J. Blanco. Feb 2016. v1_0. Fit into sophism
;   J. Blanco. Jul 2016. v1_1. Included the produced hrew model in the
;      output save file
;
; ==============================================================================

restore,'../settings/settings_sophism.sav'
progind=2

if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if progma eq -1 then progma=0

;in this program what we want to do is to simulate the effects of the HREW onto the FDT
;ds_PtV Peak-To-Vallew defocus (in waves) to simulate the HREW
;available models are (use defauts 0 for SOPHI simulations):
;0: Delta = Defocus
;1: Delta = 1/2 Defocus
;2: Delta = 1/3 Defocus
;3: Delta = Defocus*(a+bcosnTheta) ;where a=1/2, b=0.02 and n = 3
;4: Delta = Defocus*(a+bcosnTheta) ;where a=1  , b=1 and n = 3
;output is
;CM,MMatrix,MMatrix_ideal,RMatrixB
;the program generates a sav file as well
;save,filename='MATRIX_'+label+'.sav',MMatrix,MMatrix_ideal,RMatrixB,CM
;the save files are used in modelo_ret_lookuptable.pro for SOPHISM

;if not(keyword_set(steps)) then steps = 1.   ;Pixel step
;if not(keyword_set(scale)) then scale = 10.  ;10 times smaller that 2048x2048 (farter calculations, same results)
;if not(keyword_set(label)) then label = '_x' ;for the ps files
;if not(keyword_set(model)) then model = 0    ;different models for the HREW birrefringence pattern (see below)
;if not(keyword_set(ds_PtV)) then ds_PtV = 1/16.    ;Peak To Vallew defocus (see below)

;Telescope FDT characteristics

image_size_or = info.sz ;2048. ;px
image_size = image_size_or / info.birefscale ;px resampled

distance_HREW = info.dhrew ;322d0 ;Distance Between the HREW and the Entrance Pupil in milimeters [mm]
telescope_D = info.telap*1e3 ;17.5d0    ;telescope diameter [mm] (Entrance pupil)
lambda = info.wl*1d-7 ;6173d-7      ;working wavelength [mm]
telescope_f = info.flen*1e3 ;579.0d0 ;585.3d0    ;effective focal length [mm]
;pixCCD = info.fpapixsiz*1d-3 * info.birefscale ;10.D-3 * scale ;physical size of the pixel in the detector [mm]
pixCCD = double(telescope_f)*info.fpaplate/3600./!radeg * info.birefscale ;10.D-3 * scale ;physical size of the pixel in the detector [mm]
pixrad = pixCCD/telescope_f ;angular size of the pixel in radians
pix = pixrad*180*3600d0/!dpi ;angular size of the pixel in arcsec
nu_cutoff=telescope_D/lambda   ; cut frequency rad^-1
deltanu=1d0/image_size/pixrad          ; sampling interval rad^-1
rpupil=nu_cutoff/deltanu/2d0    ;pupil radius imposed by the previous paraemters [px]
print,rpupil,format='("Pupil-radius = ",F7.3," pixel.")'
pix_cut = 2*rpupil             ;cutoff freq expresed in px units (not used)
sizp=2*check_rpupil(image_size/2,rpupil)  ;extend dimen.of the Fourier dom. if needed

;see librery routines
rd = radius_aper2(sizp/2,pix_cut,tsupportp)

;-----------------------------
;define the HREW size depending on rpupil and the excursion of the optical axis
image_size_mm = image_size*pixCCD
;maximum displacement of the center of pupil image in the HREW [mm]
stepsize_xy = image_size_mm/2.*distance_HREW/telescope_f
;tamaño que susciende el pixel (pix en arcsec) a una distancia distance_HREW [mm]
tp = pix*distance_HREW/3600d0/180d0*!dpi
;radius of the image of the CCD at that distance [mm]
tp_ccd = sqrt((image_size*pix/2.)^2.+(image_size*pix/2.)^2.)*distance_HREW/3600d0/180d0*!dpi
;pupil size at the HREW position in pixels
;we suppose that a HREW pixel has the same size as pixCCD (it's easier)
factor_scala = tp / pixccd ; = distance_HREW / telescope_f

;rpupil = telescope_D/2./pixCCD*scale ;*****
;????????
;rpupil_px_hrew = rpupil ;* factor_scala   ;Do we need to scale the px in the HREW by the scale factor?
rpupil_px_hrew = rpupil* factor_scala
;????????

if (info.talk eq 1) then print,rpupil*tp*factor_scala,format='("Pupil-radius at HREW= ",F7.3," mm.")'

;The HREW image shall be, at least, the entrance pupil radius and the
;displacement 
corner = image_size_or/2.*factor_scala
tmin = (rpupil_px_hrew + corner)
tmin = 2.*sqrt(2.*tmin^2.)

cond = 0
factor = 1
while cond eq 0 do begin
HREW_size = 128.*factor
if HREW_size lt tmin then cond = 0 else cond = 1
factor = factor + 1.
endwhile

;The clear aperture of the FDT  HREW corresponds to 20 mm diameter
;HREW_size = 20./pixCCD*scale   ;*****

;Center of the HREW

xc = HREW_size/2.
yc = HREW_size/2.

;Generate the window (twice)
HREW = fltarr(HREW_size,HREW_size)
HREW2 = fltarr(HREW_size,HREW_size)

;generation of the retardance matrix image using different functional forms

;We will assume a defocus term as the radial dependence (other posibilities can be used as well)
;Units of the defocus will be in waves: ds_PtV [waves]

dz = (-1.)*8.*ds_PtV*lambda*(telescope_f/telescope_D)^2  ;Telescope defoculs [mm]
c4 = (-1.)*dz*(telescope_D/telescope_f)^2*!pi/8./sqrt(3.)/lambda ;coef.Zern.desenf.[rad]

;theta matrix (for later) (fast axis orientation. Definition is that the fast axis is oriented with the radial distance)
theta = fltarr(HREW_size,HREW_size)
for i=0.,HREW_size-1. do begin
   for j=0.,HREW_size-1 do begin
      theta[i,j] = atan(i-xc,j-yc)
   endfor 
endfor 
theta = theta*180./!pi+180.  ;between 0-360

;loop over the HREW
;this is a radial-variating function sqrt(x^2+y^2.)
for i=0.,HREW_size-1. do begin
   for j=0.,HREW_size-1. do begin
      if sqrt((i-xc)^2.+(j-yc)^2.) lt HREW_size/2. then begin
         case info.birefvarmodel of
            0: HREW(i,j) = ( 2d0*sqrt(3)*c4*((i-xc)^2.+(j-yc)^2.) ) /1e6 * 1. ;defocus
            1: HREW(i,j) = ( 2d0*sqrt(3)*c4*((i-xc)^2.+(j-yc)^2.) ) /1e6 * 2. ;1/2 Defocus
            2: HREW(i,j) = ( 2d0*sqrt(3)*c4*((i-xc)^2.+(j-yc)^2.) ) /1e6 * 3. ;1/3 Defocus
            3: HREW(i,j) =   2d0*sqrt(3)*c4/1e6 * ((i-xc)^2.+(j-yc)^2.)*(1/2.+ 0.02*cos(3.*theta[i,j]*!pi/180.)*180./!pi)
            4: begin
               HREW(i,j) = ( 2d0*sqrt(3)*c4*((i-xc)^2.+(j-yc)^2.) ) /1e6 * 1.
               HREW(i,j) = HREW(i,j) + 2d0*sqrt(3)*c4/1e6 * ((i-xc)^2.+(j-yc)^2.)*cos(3.*theta[i,j]*!pi/180.)*180./!pi
            end 
            else: HREW(i,j) = ( 2d0*sqrt(3)*c4*((i-xc)^2.+(j-yc)^2.) ) /1e6
         endcase
      endif
   endfor
endfor 
;----------------------------
;retardances matrix integration

;define a mask which is the one to be moved
mask = fltarr(HREW_size,HREW_size)
for i=0.,HREW_size-1. do begin
   for j=0.,HREW_size-1 do begin
      if sqrt((i-xc)^2.+(j-yc)^2.) lt rpupil_px_hrew then mask(i,j) = 1.
  endfor 
endfor 

RMatrixB=fltarr(image_size,image_size,4,4)
cosHREW = cos(HREW*!dpi/180d0)
sinHREW = sin(HREW*!dpi/180d0)
cos2THETA = cos(2d0*theta*!dpi/180d0)
sin2THETA = sin(2d0*theta*!dpi/180d0)

xcc = image_size/2d0   ;camera px
ycc = image_size/2d0


;st = 0
for i=0,image_size-1. do begin ;,steps do begin

   for j=0,image_size-1. do begin ;,steps do begin

      xp = (i-xcc)*factor_scala*info.birefscale
      yp = (j-ycc)*factor_scala*info.birefscale
      mask_sh = shift(mask,[xp,yp])
      
;stop      
      RMatrixB[i,j,3,3] = total(cosHREW*mask_sh) / total(mask_sh)
      RMatrixB[i,j,2,2] = total( (sin2THETA^2.+cos2THETA^2.*cosHREW)*mask_sh) / total(mask_sh)
      RMatrixB[i,j,1,1] = total( (cos2THETA^2.+sin2THETA^2.*cosHREW)*mask_sh) / total(mask_sh)
      RMatrixB[i,j,2,3] = total(cos2THETA*sinHREW*mask_sh) / total(mask_sh) *(-1.)
      RMatrixB[i,j,1,3] = total(sin2THETA*sinHREW*mask_sh) / total(mask_sh)
      RMatrixB[i,j,2,1] = total((sin2THETA*cos2THETA*(1.-cosHREW))*mask_sh) / total(mask_sh)
      
      if (info.talk eq 1) then print,format='(41(%"\b"),"Loop in x and y ", I4," ",I4," ",F7.3," ",F7.3,$)',i,j,xp,yp

   endfor

   
endfor
;to leave space after the talk loop
print,''

RMatrixB[*,*,3,2] = RMatrixB[*,*,2,3]*(-1.)
RMatrixB[*,*,3,1] = RMatrixB[*,*,1,3]*(-1.)
RMatrixB[*,*,2,1] = RMatrixB[*,*,1,2]*(-1.)
RMatrixB[*,*,0,0] = 1.

;me construyo ahora un retardador como en soPHI
ret1=double(info.retdeg1)
ret2=double(info.retdeg2)

;linear polarizer
PL=[[1,1,0,0],$
    [1,1,0,0],$
    [0,0,0,0],$
    [0,0,0,0]]

;MMatrix_ideal = pmp(ret1,ret2)

MMatrix_ideal=dblarr(4,4)
for jj=0,3 do begin
   lcvr1=transpose(device(amr=1.,pha=ret1(jj),ang=info.retang1)) 
   lcvr2=transpose(device(amr=1.,pha=ret2(jj),ang=info.retang2))
   mm1=transpose(PL#lcvr2#lcvr1)
   MMatrix_ideal(0:3,jj)=mm1(0:3)
endfor

MMatrix=fltarr(image_size_or/info.birefscale,image_size_or/info.birefscale,4,4)
CM=fltarr(image_size_or/info.birefscale,image_size_or/info.birefscale,4,4)
ME=fltarr(image_size_or/info.birefscale,image_size_or/info.birefscale,4,4)

for i=0./info.birefscale,image_size_or/info.birefscale-1. do begin
   for j=0./info.birefscale,image_size_or/info.birefscale-1. do begin
      for k=0,3 do begin
         lcvr1=device(amr=1.,pha=ret1[k],ang=double(info.retang1))
         MM = LCVR1##reform(RMatrixB[i,j,*,*])
         lcvr2=device(amr=1.,pha=ret2[k],ang=double(info.retang2))
         MM = LCVR2##MM
         MM = PL##MM
         MMatrix[i,j,*,k] = MM[*,0]
      endfor
      CM[i,j,*,*] = MMatrix_ideal##Invert(reform(MMatrix[i,j,*,*])) ;crosstalk
      ME[i,j,*,*] =(invert(reform(MMatrix_ideal[*,*]))-invert(reform(MMatrix[i,j,*,*])))/invert(reform(MMatrix[i,j,*,*]))*100.
   endfor
;   print,i
endfor
;stop
save,filenam=info.saves+info.files(progind)+'_biref_model.sav',hrew,MMatrix,MMatrix_ideal,RMatrixB,CM


;stop



end
