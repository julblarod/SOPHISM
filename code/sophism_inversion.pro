pro sophism_inversion

; ==================================================================
; DESCRIPTION
;   Simulate the inversion data of SOPHI, using C version of MILOS
;
;
; CALLED FROM
;   sophism
;
; SUBROUTINES
;
;
; EXTERNAL CODE
;   C-Milos          C version of MILOS inversion code
;
;
; MAIN INPUTS
;    invercont        Index position in wavelength dimension where
;                     continuum is observed
;    invercent        Central wavelength
;    inverquanlow     Quantum numbers of lower level
;    inverquanup      Quantum numbers of upper level
;    inversigma       Noise threshold to perform inversion
;    iter             Maximum number of iterations
;    x0,y0,fovx,fovy  Define size of area over which to invert
;    
;
; OUTPUT
;    Fits file with inversion products. Specifically (and in this
;    order in the array): continuum intensity, vlos (km/s), field
;    strength (Gauss), field inclination, field azimuth
;
;
; VERSION HISTORY
;
;   J. Blanco. Mar 2017. v1_0. Inversion using C version of MILOS.
;   J. Blanco. Apr 2017. v1_1. Include option to perform inversion
;      only on a given area of the FOV (to speed things up).
;
; ==================================================================


progind=9

print,''
print,'sophism_inversion'
print,''

;load settings file
restore,'../settings/settings_sophism.sav'
;should check if there is demodulated data before going into work?


;define the first part of the input and output filenames because later
;the current folder will change and it may cause problems
filepref=info.saves+info.files(progind-1);+'_av'
fileoutpref=info.saves+info.files(progind);+'_av'


llo_f=info.wl
dop=info.scvel/299792.5d0
ldop=llo_f*dop
lambda0=llo_f+ldop
;lambda0=info.wl
iteration=info.iter

;use classical estimates? 0=no. This is fixed for the moment
clasest=1
;[FWHM DELTA NPOINTS] use convolution with a gaussian? if the tree
;parameteres are defined yes, else no. Units in A. NPOINTS has to be odd.

ima=readfits(info.saves+info.files(8)+'_'+'0.fits*',head,/comp,/sil)
sizim=size(ima)
;for MILOS, spectral fwhm needed
;if the etalon module was run, use the spectral fwhm calculated there
;If not, go for a 'close' one (technically, it HAS to be run)
headmod=sxpar(head,'HISTORY')
etal=where(headmod eq 'Etalon Module')
if (etal ne -1) then begin
   restore,info.saves+info.files(3)+'.sav'
   lambda=llo
;fwhm_real is in AA
   fwhmet=mean(fwhm_real)*1e3
endif else begin
   lambda=info.lll
   fwhmet=85
endelse

;since header will be common, update it here and use it for all writefits
sxaddpar,head,"HISTORY",'Inversion Module'
sxaddpar,head,"HISTORY",'0-Ic, 1-Vlos, 3-B, 4-Gamma, 5-Azim'

;in OSX to intrepret correctly the non-unicode whatever characters  
;export LC_CTYPE=C                        
;export LANG=C                                                        

;change directory to milos
cd,'../cmilos-master/',current=dirold 
;modify settings depending on the number of spectral positions and
;central wavelength
;spawn,"sed -i '41 c\'$'\n''#define NLAMBDA "+strtrim(sizim(1),2)+" //numero de muestras de entrada'$'\n' defines.h"
;spawn,"sed -i'.orig' 's/\(#define NLAMBDA\).*/\1 "+strtrim(sizim(1),2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CENTRAL_WL\).*/\1 "+strtrim(info.invercent,2)+" /' defines.h"
;for now, only one wavelength allowed
;spawn,"sed -i 's/\(#define CUANTIC_NWL\).*/\1 1 /' defines.h";+strtrim(info.invernwl,2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CUANTIC_SLOI\).*/\1 "+strtrim(info.inverquanlow[0],2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CUANTIC_LLOI\).*/\1 "+strtrim(info.inverquanlow[1],2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CUANTIC_JLOI\).*/\1 "+strtrim(info.inverquanlow[2],2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CUANTIC_SUPI\).*/\1 "+strtrim(info.inverquanup[0],2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CUANTIC_LUPI\).*/\1 "+strtrim(info.inverquanup[1],2)+" /' defines.h"
spawn,"sed -i'.orig' 's/\(#define CUANTIC_JUPI\).*/\1 "+strtrim(info.inverquanup[2],2)+" /' defines.h"
;re-compile the code
spawn,'make clear'
spawn,'make'

;apply inversion only to a limited FOV
if (info.inversmall eq 1) then begin
   x0=info.inver0[0]
   y0=info.inver0[1]
   fovx=info.inverfov[0]
   fovy=info.inverfov[1]
endif else begin
   x0=0
   y0=0
   fovx=sizim(3)
   fovy=sizim(4)
endelse

;array holding the wavelengths and Stokes I,Q,U,V for all the spectral positions
profil=fltarr(5,sizim(1))
profil[0,*]=lambda
print,'Take only strong profiles'
for obs=0,info.obsset-1 do begin
   strength=fltarr(fovx,fovy)
   gamma=strength
   azimuth=strength
   vlos=strength

;input fits file from demodulation
;   file=info.saves+info.files(progind-1)+'_'+strtrim(obs,2)+'.fits*'
   file=filepref+'_'+strtrim(obs,2)+'.fits*'
;output filename
;   fileout=info.saves+info.files(progind)+'_'+strtrim(obs,2)
   fileout=fileoutpref+'_'+strtrim(obs,2)
   ima=readfits(file,/comp,/sil)
   cont=mean(ima[info.invercont,0,*,*])
   rmsQc=stddev(ima[info.invercont,1,*,*]/cont)
   rmsUc=stddev(ima[info.invercont,2,*,*]/cont)
   rmsVc=stddev(ima[info.invercont,3,*,*]/cont)

   for xx=x0,x0+fovx-1 do begin
      for yy=y0,y0+fovy-1 do begin
         profil[1:4,*]=transpose(reform(ima[*,*,xx,yy]))/cont
;el no invertir perfiles por debajo del ruido de la version idl
;probar y, si acaso, meterlo en los teo, xtalk, teo_xtalk
         fixa=1
         fixb=1
         fixc=1
         if  (max(abs(profil[2,*])) lt (info.inversigma*rmsQc) and max(abs(profil[3,*])) lt (info.inversigma*rmsUc)) then begin
;por ahora, saltarse toda la inversion si salen perfiles
;cutrongos. Luego ya veremos si hacemos lo del fix
            goto,noinver
            fixa=0
            if (max(abs(profil[4,*])) lt (info.inversigma*rmsVc)) then begin
               fixb=0
               fixc=0
            endif
;there is only one sentence with 'int fix'. This search may be a bit
;too general for the future. Be careful.
            spawn,"sed -i'.orig' 's/\(int fix\).*/\1[]={1.,"+strtrim(fixc,2)+".,1.,1.,1.,"+strtrim(fixb,2)+".,"+strtrim(fixa,2)+".,1.,1.,0.,0.}; /' milos.c"
            spawn,'make clear'
            spawn,'make'
         endif

         openw,luni,info.saves+'perfilado.per',/get_lun
         printf,luni,profil
         close,/all
;milos NLAMBDA MAX_ITER CLASSICAL_ESTIMATES [FWHM DELTA NPOINTS] prof.txt > output.txt
;         spawn,'./milos '+strtrim(sizim(1),2)+' '+strtrim(iteration,2)+' '+strtrim(clasest,2)+' 0.09 0.04 9 '+info.saves+'perfilado.per > '+info.saves+'perfiladores.txt'
         spawn,'./milos '+strtrim(sizim(1),2)+' '+strtrim(iteration,2)+' '+strtrim(clasest,2)+' '+info.saves+'perfilado.per > '+info.saves+'perfiladores.txt'
         res=dblarr(12)
         openr,lunii,info.saves+'perfiladores.txt',/get_lun
         readf,lunii,res
         close,/all

         strength[xx-x0,yy-y0]=res(2)
         gamma[xx-x0,yy-y0]=res(3)
         azimuth[xx-x0,yy-y0]=res(4)
         vlos[xx-x0,yy-y0]=res(8)       

         if (fixa eq 0) then begin
            fixa=1
            fixb=1
            fixc=1
            spawn,"sed -i'.orig' 's/\(int fix\).*/\1[]={1.,"+strtrim(fixc,2)+".,1.,1.,1.,"+strtrim(fixb,2)+".,"+strtrim(fixa,2)+".,1.,1.,0.,0.}; /' milos.c"
            spawn,'make clear'
            spawn,'make'
         endif 
         noinver:
      endfor
      if (info.talk eq 1) then begin
;         print,''
         print,format='(21(%"\b"),"Inverting column ",I4,$)',xx+1
;         print,''
      endif 
   endfor

;create array to hold the outputs. Keeping the original structure of
;arrays
   imout=fltarr(5,1,fovx,fovy)
   imout[0,0,*,*]=reform(ima[info.invercont,0,x0:x0+fovx-1,y0:y0+fovy-1])
   imout[1,0,*,*]=vlos
   imout[2,0,*,*]=strength
   imout[3,0,*,*]=gamma
   imout[4,0,*,*]=azimuth

   if (info.compress eq 1) then writefits,fileout+'.fits',imout,head,/compress else writefits,fileout+'.fits',imout,head

endfor ;obs

;display some results
if (info.talk eq 1) then begin
;   restore,fileout+'_MILOS_inversion.sav'
   !p.multi=[0,2,2]
   window,title='MILOS Inversion Results',/free
   tvframe,vlos,/asp,/bar,btit='v!Dlos!N (km/s)'
   tvframe,strength,/asp,/bar,btit='B (G)'
   tvframe,gamma,/asp,/bar,btit='gamma (deg)'
   tvframe,azimuth,/asp,/bar,btit='azimuth (deg)'
   !p.multi=0
endif

;****************************************
;include inversion for theoretical(wrong) demodulation
if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
   if (info.talk eq 1) then begin
      print,''
      print,''
      print,'Inversion for theoretical demodulation'
      print,''
      print,''
   endif 
   
   for obs=0,info.obsset-1 do begin 
      strength=fltarr(sizim(3),sizim(4))
      gamma=strength
      azimuth=strength
      vlos=strength
      eta0=strength

;input fits file from demodulation
;      file=info.saves+info.files(progind-1)+'_teo_'+strtrim(obs,2)+'.fits*'
      file=filepref+'_teo_'+strtrim(obs,2)+'.fits*'
;output filename
;      fileout=info.saves+info.files(progind)+'_teo_'+strtrim(obs,2)
      fileout=fileoutpref+'_teo_'+strtrim(obs,2)

      ima=readfits(file,/comp,/sil)
      cont=mean(ima[info.invercont,0,*,*])
      for xx=0,sizim(3)-1 do begin
         for yy=0,sizim(4)-1 do begin
            profil[1:4,*]=transpose(reform(ima[*,*,xx,yy]))/cont
            openw,luni,info.saves+'perfilado.per',/get_lun
            printf,luni,profil
            close,/all
            spawn,'./milos '+strtrim(iteration,2)+' '+info.saves+'perfilado.per > '+info.saves+'perfiladores.txt'
            res=dblarr(12)
            openr,lunii,info.saves+'perfiladores.txt',/get_lun
            readf,lunii,res
            close,/all
            
            strength[xx,yy]=res(2)
            gamma[xx,yy]=res(3)
            azimuth[xx,yy]=res(4)
            vlos[xx,yy]=res(8)       
         endfor
      endfor

;write to fits, for use in compression, or others
;restore output of milos and define size
;      restore,fileout+'_MILOS_inversion.sav'
;      sz=size(vlos)

;create array to hold the outputs. Keeping the original structure of
;arrays
;      imout=fltarr(5,1,sz(1),sz(2))
      imout=fltarr(5,1,sizim(3),sizim(4))
      imout[0,0,*,*]=eta0
      imout[1,0,*,*]=vlos
      imout[2,0,*,*]=strength
      imout[3,0,*,*]=gamma
      imout[4,0,*,*]=azimuth

      if (info.compress eq 1) then writefits,fileout+'.fits',imout,head,/compress else writefits,fileout+'.fits',imout,head

   endfor ;obsteo

   if (info.talk eq 1) then begin
;      restore,fileout+'_MILOS_inversion.sav'
      !p.multi=[0,2,2]
      window,title='MILOS Inversion Results (Theoretical Demod.)',/free
      tvframe,vlos,/asp,/bar,btit='v!Dlos!N (km/s)'
      tvframe,strength,/asp,/bar,btit='B (G)'
      tvframe,gamma,/asp,/bar,btit='gamma (deg)'
      tvframe,azimuth,/asp,/bar,btit='azimuth (deg)'
      !p.multi=0
   endif
endif 

;****************************************
;include inversion for adhoc corrected demodulation
if (info.adhoc eq 1) then begin
   if (info.talk eq 1) then begin
      print,''
      print,''
      print,'Inversion for adhoc corrected demodulation'
      print,''
      print,''
   endif  
   
   for obs=0,info.obsset-1 do begin 
      strength=fltarr(sizim(3),sizim(4))
      gamma=strength
      azimuth=strength
      vlos=strength
      eta0=strength

;input fits file from demodulation
;      file=info.saves+info.files(progind-1)+'_xtalk_'+strtrim(obs,2)+'.fits*'
      file=filepref+'_xtalk_'+strtrim(obs,2)+'.fits*'
;output filename
;      fileout=info.saves+info.files(progind)+'_xtalk_'+strtrim(obs,2)
      fileout=fileoutpref+'_xtalk_'+strtrim(obs,2)

      ima=readfits(file,/comp,/sil)
      cont=mean(ima[info.invercont,0,*,*])
      for xx=0,sizim(3)-1 do begin
         for yy=0,sizim(4)-1 do begin
            profil[1:4,*]=transpose(reform(ima[*,*,xx,yy]))/cont
            openw,luni,info.saves+'perfilado.per',/get_lun
            printf,luni,profil
            close,/all
            spawn,'./milos '+strtrim(iteration,2)+' '+info.saves+'perfilado.per > '+info.saves+'perfiladores.txt'
            res=dblarr(12)
            openr,lunii,info.saves+'perfiladores.txt',/get_lun
            readf,lunii,res
            close,/all
            
            strength[xx,yy]=res(2)
            gamma[xx,yy]=res(3)
            azimuth[xx,yy]=res(4)
            vlos[xx,yy]=res(8)       
         endfor
      endfor

;write to fits, for use in compression, or others
;restore output of milos and define size
;      restore,fileout+'_MILOS_inversion.sav'
;      sz=size(vlos)

;create array to hold the outputs. Keeping the original structure of
;arrays
;      imout=fltarr(5,1,sz(1),sz(2))
      imout=fltarr(5,1,sizim(3),sizim(4))
      imout[0,0,*,*]=eta0
      imout[1,0,*,*]=vlos
      imout[2,0,*,*]=strength
      imout[3,0,*,*]=gamma
      imout[4,0,*,*]=azimuth

      if (info.compress eq 1) then writefits,fileout+'.fits',imout,head,/compress else writefits,fileout+'.fits',imout,head

   endfor ;obsteo

   if (info.talk eq 1) then begin
;      restore,fileout+'_MILOS_inversion.sav'
      !p.multi=[0,2,2]
      window,title='MILOS Inversion Results (Adhoc Demod.)',/free
      tvframe,vlos,/asp,/bar,btit='v!Dlos!N (km/s)'
      tvframe,strength,/asp,/bar,btit='B (G)'
      tvframe,gamma,/asp,/bar,btit='gamma (deg)'
      tvframe,azimuth,/asp,/bar,btit='azimuth (deg)'
      !p.multi=0
   endif 

   if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
      if (info.talk eq 1) then begin
         print,''
         print,''
         print,'Inversion for theoretical+adhoc corrected demodulation'
         print,''
         print,''
      endif 
   
      for obs=0,info.obsset-1 do begin 
         strength=fltarr(sizim(3),sizim(4))
         gamma=strength
         azimuth=strength
         vlos=strength
         eta0=strength

;         file=info.saves+info.files(progind-1)+'_teo_xtalk_'+strtrim(obs,2)+'.fits*'
         file=filepref+'_teo_xtalk_'+strtrim(obs,2)+'.fits*'
;         fileout=info.saves+info.files(progind)+'_teo_xtalk_'+strtrim(obs,2)
         fileout=fileoutpref+'_teo_xtalk_'+strtrim(obs,2)
         
         ima=readfits(file,/comp,/sil)
         cont=mean(ima[info.invercont,0,*,*])
         for xx=0,sizim(3)-1 do begin
            for yy=0,sizim(4)-1 do begin
               profil[1:4,*]=transpose(reform(ima[*,*,xx,yy]))/cont
               openw,luni,info.saves+'perfilado.per',/get_lun
               printf,luni,profil
               close,/all
               spawn,'./milos '+strtrim(iteration,2)+' '+info.saves+'perfilado.per > '+info.saves+'perfiladores.txt'
               res=dblarr(12)
               openr,lunii,info.saves+'perfiladores.txt',/get_lun
               readf,lunii,res
               close,/all
               
               strength[xx,yy]=res(2)
               gamma[xx,yy]=res(3)
               azimuth[xx,yy]=res(4)
               vlos[xx,yy]=res(8)       
            endfor
         endfor

;         restore,fileout+'_MILOS_inversion.sav'
;         sz=size(vlos)
;         imout=fltarr(5,1,sz(1),sz(2))
         imout=fltarr(5,1,sizim(3),sizim(4))
         imout[0,0,*,*]=eta0
         imout[1,0,*,*]=vlos
         imout[2,0,*,*]=strength
         imout[3,0,*,*]=gamma
         imout[4,0,*,*]=azimuth
         
         if (info.compress eq 1) then writefits,fileout+'.fits',imout,head,/compress else writefits,fileout+'.fits',imout,head
         
      endfor ;obsteo

      if (info.talk eq 1) then begin
         !p.multi=[0,2,2]
         window,title='MILOS Inversion Results (Theoretical+adhoc Demod.)',/free
         tvframe,vlos,/asp,/bar,btit='v!Dlos!N (km/s)'
         tvframe,strength,/asp,/bar,btit='B (G)'
         tvframe,gamma,/asp,/bar,btit='gamma (deg)'
         tvframe,azimuth,/asp,/bar,btit='azimuth (deg)'
         !p.multi=0
      endif 
      

   endif ;adhoc+errors

endif ;adhoc  

;****************************************
;include inversion for average and 2D FOV demodulation?

cd,dirold
noinv:

end
