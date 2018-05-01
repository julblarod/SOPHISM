PRO sophism_demod

; ==================================================================
; DESCRIPTION
;    Demodulates the accumulated images of the simulation run with the
;    demodulation matrix calculated at the polarization module.
;    Adds the two beams of dual-beam mode.
;
; CALLED FROM 
;    sophism
;
; SUBROUTINES
;    sophism_demod_adhoc
;
; MAIN INPUTS
;    dmodm/dmodm_teo      Demodulation matrixes, real one (with actual
;                            parameters for LCVRs, including random
;                            errors if enabled) and theoretical one
;                            ('perfect' LCVRs)
;    reterrors            Enabled random errors in LCVRs parameters in
;                            polarization module
;    muellerrors          Enabled random errors in Mueller matrix
;                            polarization module
;
; OUTPUT
;    Final demodulated images from demodulation matrixes with and
;    without errors (in case they were applied), and images corrected
;    from crosstalk from V to Q and U if adhoc correction was enabled.
;
; VERSION HISTORY
;    J. Blanco. 2011.
;    J. Blanco. Feb 2012. v1_1. Demodulation for random errors in matrix 
;    J. Blanco. Nov 2012. v1_2. Added the adhoc correction option
;    J. Blanco. Dec 2012. v1_3. Show demodulated images for talk mode
;    J. Blanco. Dec 2012. v1_4. Rearrange 'show' before ad-hoc correction
;    J. Blanco. Jan 2013. v1_5. Add demodulation with theoretical
;       matrix also in case of errors in Mueller matrix
;    J. Blanco. Jan 2014. v1_6. Demodulation matrixes stored in 
;       ..._scheme_demod.sav starting in sophism_polmeas_modscheme
;       3.0. Adapted to 2D demodulation matrixes and averaged.
;    J. Blanco. Jun 2014. v1_7. Demodulation matrixes with lambda
;       dependence (with or without FOV-dep.). All 4 cases (simple
;       modulation, lambda-modulation, FOV-modulation,
;       FOVandlambda-modulation) considered
;    J. Blanco. Jul 2014. v1_8. Adapted for observation sets.
;    J. Blanco. Aug 2014. v1_9. Readapted modulation cycles' demodulation.
;    J. Blanco. Sep 2014. v2_0. Correction for 'real' cases with
;       FOV/lambda dependencies and the demodulation matrixes they generate
;    J. Blanco. Nov 2014. v2_1. Corrected bug when reading data for
;       adhoc. Changed output names for adhoc data.
;    J. Blanco. Jul 2015. v2_2. Corrected worng variable name for
;       dual-beam wrong-demodulation
;    J. Blanco. Feb 2016. v2_3. Corrected bug continuum position for
;       ad-hoc subroutine
;    J. Blanco. Apr 2017. v2_31. Allow for starting with data other
;       than accumulated (mostly thinking in no data series)
;
; ==================================================================


print,''
print,'sophism_demod'
print,''
print,'use only after accumulation (and polarization)'

restore,'../settings/settings_sophism.sav'
progind=8

;careful!! Everything is done considering 4 modulation states. May be
;useful even in longitudinal (I,V) case, provided it is performed
;twice (2x2) in modulation. Give some thought.
;nperiod=info.modnbuff*info.nacc
;nmcycl=info.cycrep
;ntim=nperiod*nmcycl*info.etalpos

;el -6... alguna forma mas general?
;check for compatibility with older versions where no demod file
dd=file_search(info.saves+info.files(progind-6)+'_scheme_demod.sav',count=nn)
if (nn ne 0) then restore,info.saves+info.files(progind-6)+'_scheme_demod.sav' else restore,info.saves+info.files(progind-6)+'_scheme.sav'

;only accumulated data enter the demodulation
;progma=progind-1
;but now maybe not, because of flexibility and to allow no dataser,
;i.e. no accumulation module (no data to accumulate)
if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
;not sure this line below makes sense here but maybe useful, leave it
if (progma eq -1) then progma=0

for obs=0,info.obsset-1 do begin

   ima=readfits(info.saves+info.files(progma)+'_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)
   sizei=size(ima)

   demodima=ima*0.

;in case the modulation cycles are gt 1, add together the different
;cycles before demodulating
;the accumulated output will conserve all the repetitions in case some
;test want to be done, but to continue with inversion, here is needed
;to add them up 
   if (info.cycrep gt 1 AND info.dataser ne 0) then begin
;considered that the cycles are stored in the stokes dimension, so it
;has to be downed
      buff=sizei(2)/info.cycrep
      imaint=fltarr(sizei(1),buff,sizei(3),sizei(4))
      for cc=0,info.cycrep-1 do imaint=imaint+ima[*,cc*buff:cc*buff+buff-1,*,*]
      ima=imaint
      sizei=size(ima)
      demodima=ima*0.
   endif 


;demodulation in loop form
;check whether to do the loop in FOV or not necessary
   case (size(dmodm))(0) of
;one constant modulation for the whole FOV and wavelength
      2: begin
         for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima(*,j,*,*)=demodima(*,j,*,*)+dmodm(i,j)*ima(*,i,*,*)
;if random errors were introduced in the modulation, perform a second
;demodulation with the theoretical (wrong) demodulation matrix
         if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
            demodima_teo=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo(*,j,*,*)=demodima_teo(*,j,*,*)+dmodm_teo(i,j)*ima(*,i,*,*)
         endif
         end
;case of modulation changing with wavelength
      3: begin
         demodimaav=ima*0.
         for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav(*,j,*,*)=demodimaav(*,j,*,*)+dmodmav(i,j)*ima(*,i,*,*)
;   for l=0,info.etalpos-1 do for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
         for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+dmodm(i,j,l)*reform(ima(l,i,*,*))
         if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
            demodimaav_teo=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav_teo(*,j,*,*)=demodimaav_teo(*,j,*,*)+dmodmav_teo(i,j)*ima(*,i,*,*)
            demodima_teo=ima*0.
            for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo(l,j,*,*)=demodima_teo(l,j,*,*)+dmodm_teo(i,j,l)*reform(ima(l,i,*,*))
         endif
         end
;case of modulation with FOV dependence
      4: begin
         demodimaav=ima*0.
         for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav(*,j,*,*)=demodimaav(*,j,*,*)+dmodmav(i,j)*ima(*,i,*,*)
;   for l=0,info.etalpos-1 do for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
         for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
         if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
            demodimaav_teo=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav_teo(*,j,*,*)=demodimaav_teo(*,j,*,*)+dmodmav_teo(i,j)*ima(*,i,*,*)
            demodima_teo=ima*0.
            for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo(l,j,*,*)=demodima_teo(l,j,*,*)+reform(dmodm_teo(i,j,*,*))*reform(ima(l,i,*,*))
         endif
         end 
;case of modulation dependent on FOV and wavelength
      5: begin
         demodimaav=ima*0.
         for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav(*,j,*,*)=demodimaav(*,j,*,*)+dmodmav(i,j)*ima(*,i,*,*)
;   for l=0,info.etalpos-1 do for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
         for l=0,sizei(1)-1 do begin
;in case the dependency on lambda is not as 'long' as the wavelength
;dimension of the data
            ld=l mod (size(dmodm))(3)
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,ld,*,*))*reform(ima(l,i,*,*))
         endfor 
         if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
            demodimaav_teo=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav_teo(*,j,*,*)=demodimaav_teo(*,j,*,*)+dmodmav_teo(i,j)*ima(*,i,*,*)
            demodima_teo=ima*0.
            for l=0,sizei(1)-1 do begin 
               ld=l mod (size(dmodm))(3)
               for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo(l,j,*,*)=demodima_teo(l,j,*,*)+reform(dmodm_teo(i,j,ld,*,*))*reform(ima(l,i,*,*))
            endfor 
         endif
         end
   endcase 

;##################################################################

   if (info.dualbeam eq 1) then begin
      if (info.talk eq 1 AND obs eq 0) then print,'Starting demodulation of dual-beam'

      ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)
      demodima2=ima*0.
;demodulation in loop form
;      for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima2(*,j,*,*)=demodima2(*,j,*,*)+dmodm2(i,j)*ima(*,i,*,*)
;if random errors were introduced in the modulation, perform a second
;demodulation with the theoretical (wrong) demodulation matrix
;      if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
;         demodima_teo2=ima*0.
;         for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima_teo2(*,j,*,*)=demodima_teo2(*,j,*,*)+dmodm_teo2(i,j)*ima(*,i,*,*)
;      endif
      
      if (info.cycrep gt 1) then begin
         sizei=size(ima)
;considered that the cycles are stored in the stokes dimension, so it
;has to be downed
         buff=sizei(2)/info.cycrep
         imaint=fltarr(sizei(1),buff,sizei(3),sizei(4))
         for cc=0,info.cycrep-1 do imaint=imaint+ima[*,cc*buff:cc*buff+buff-1,*,*]
         ima=imaint
         sizei=size(ima)
         demodima2=ima*0.
      endif 

      case (size(dmodm))(0) of
;one constant modulation for the whole FOV and wavelength
         2: begin
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima2(*,j,*,*)=demodima2(*,j,*,*)+dmodm2(i,j)*ima(*,i,*,*)
;if random errors were introduced in the modulation, perform a second
;demodulation with the theoretical (wrong) demodulation matrix
            if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
               demodima_teo2=ima*0.
               for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo2(*,j,*,*)=demodima_teo2(*,j,*,*)+dmodm2_teo(i,j)*ima(*,i,*,*)
            endif
            end
;case of modulation changing with wavelength
         3: begin
            demodimaav2=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav2(*,j,*,*)=demodimaav2(*,j,*,*)+dmodmav2(i,j)*ima(*,i,*,*)
;   for l=0,info.etalpos-1 do for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
            for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima2(l,j,*,*)=demodima2(l,j,*,*)+dmodm2(i,j,l)*reform(ima(l,i,*,*))
            if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
               demodimaav_teo2=ima*0.
               for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav_teo2(*,j,*,*)=demodimaav_teo2(*,j,*,*)+dmodmav_teo2(i,j)*ima(*,i,*,*)
               demodima_teo2=ima*0.
               for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo2(l,j,*,*)=demodima_teo2(l,j,*,*)+dmodm_teo2(i,j,l)*reform(ima(l,i,*,*))
            endif
            end
;case of modulation with FOV dependence
         4: begin
            demodimaav2=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav2(*,j,*,*)=demodimaav2(*,j,*,*)+dmodmav2(i,j)*ima(*,i,*,*)
;   for l=0,info.etalpos-1 do for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
            for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima2(l,j,*,*)=demodima2(l,j,*,*)+reform(dmodm2(i,j,*,*))*reform(ima(l,i,*,*))
            if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
               demodimaav_teo2=ima*0.
               for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav_teo2(*,j,*,*)=demodimaav_teo2(*,j,*,*)+dmodmav_teo2(i,j)*ima(*,i,*,*)
               demodima_teo2=ima*0.
               for l=0,sizei(1)-1 do for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo2(l,j,*,*)=demodima_teo2(l,j,*,*)+reform(dmodm_teo2(i,j,*,*))*reform(ima(l,i,*,*))
            endif
            end 
;case of modulation dependent on FOV and wavelength
         5: begin
            demodimaav2=ima*0.
            for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav2(*,j,*,*)=demodimaav2(*,j,*,*)+dmodmav2(i,j)*ima(*,i,*,*)
;   for l=0,info.etalpos-1 do for j=0,info.modnbuff-1 do for i=0,info.modnbuff-1 do demodima(l,j,*,*)=demodima(l,j,*,*)+reform(dmodm(i,j,*,*))*reform(ima(l,i,*,*))
            for l=0,sizei(1)-1 do begin
               ld=l mod (size(dmodm))(3)
               for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima2(l,j,*,*)=demodima2(l,j,*,*)+reform(dmodm2(i,j,ld,*,*))*reform(ima(l,i,*,*))
            endfor 
            if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
               demodimaav_teo2=ima*0.
               for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodimaav_teo2(*,j,*,*)=demodimaav_teo2(*,j,*,*)+dmodmav_teo2(i,j)*ima(*,i,*,*)
               demodima_teo2=ima*0.
               for l=0,sizei(1)-1 do begin
                  ld=l mod (size(dmodm))(3)
                  for j=0,info.modnbuff/info.cycrep-1 do for i=0,info.modnbuff/info.cycrep-1 do demodima_teo2(l,j,*,*)=demodima_teo2(l,j,*,*)+reform(dmodm_teo2(i,j,ld,*,*))*reform(ima(l,i,*,*))
               endfor 
            endif
            end
      endcase 


;adding of the two beams for dual-beam mode
      demodima=demodima+demodima2
      if ((size(dmodm))(0) gt 2) then demodimaav=demodimaav+demodimaav2
      if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
         demodima_teo=demodima_teo+demodima_teo2
         if ((size(dmodm))(0) gt 2) then demodimaav_teo=demodimaav_teo+demodimaav_teo2
      endif 

   endif ;dualbeam

;###################################################################

   heador=head
   sxaddpar,head,"HISTORY",'Demodulated images'
   if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',demodima,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',demodima,head

   headth=heador
   if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
      sxaddpar,headth,"HISTORY",'Theoretical (wrong) demodulated images'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_teo_'+strtrim(obs,2)+'.fits',demodima_teo,headth,/compress else writefits,info.saves+info.files(progind)+'_teo_'+strtrim(obs,2)+'.fits',demodima_teo,headth
   endif 

   if ((size(dmodm))(0) gt 2) then begin
      headorth=heador
      sxaddpar,heador,"HISTORY",'Demodulated images with averaged dem. matrix'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_av_'+strtrim(obs,2)+'.fits',demodimaav,heador,/compress else writefits,info.saves+info.files(progind)+'_av_'+strtrim(obs,2)+'.fits',demodimaav,heador
      if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
         sxaddpar,headorth,"HISTORY",'Theoretical (wrong) demodulated images with av. dem. mat.'
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_teo_av_'+strtrim(obs,2)+'.fits',demodimaav_teo,headorth,/compress else writefits,info.saves+info.files(progind)+'_teo_av_'+strtrim(obs,2)+'.fits',demodimaav_teo,headorth
      endif 
   endif 

endfor ;obs

;display demodulated images example
if (info.talk eq 1) then begin
   !p.multi=[0,info.modnbuff/info.cycrep,info.etalpos]
   !x.omargin=[0,80]
   window,xs=700,ys=info.etalpos*135,title='Demodulated Images at selected wavelengths',/free
   for hh=0,(info.etalpos/info.obsset)-1 do begin
      for jj=0,info.modnbuff/info.cycrep-1 do tvframe,demodima[hh,jj,*,*],/asp,charsize=0.8
      xyouts,550,info.etalpos*135-60-hh*135,'!4k!3!D'+strtrim(hh,2)+'!N',/dev,charsiz=1.3
   endfor

   if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
      window,xs=700,ys=info.etalpos*135,title='Wrong-Demodulated Images at selected wavelengths',/free
      for hh=0,info.etalpos-1 do begin
         for jj=0,info.modnbuff/info.cycrep-1 do tvframe,demodima_teo[hh,jj,*,*],/asp,charsize=0.8
         xyouts,550,info.etalpos*135-60-hh*135,'!4k!3!D'+strtrim(hh,2)+'!N',/dev,charsiz=1.3
      endfor
   endif 
   !p.multi=0
   !x.omargin=0
endif 


;apply the adhoc crosstalk correction
if (info.adhoc eq 1) then begin
   print,'Performing ad-hoc crosstalk correction'

   for obs=0,info.obsset-1 do begin

      demodima=readfits(info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)

      sophism_demod_adhoc,demodima,info.sscene_apod,info.invercont

      sxaddpar,head,"HISTORY",'Adhoc corrected'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_xtalk_'+strtrim(obs,2)+'.fits',demodima,head,/compress $
      else writefits,info.saves+info.files(progind)+'_xtalk_'+strtrim(obs,2)+'.fits',demodima,head

      if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
         demodima_teo=readfits(info.saves+info.files(progind)+'_teo_'+strtrim(obs,2)+'.fits*',headth,/comp,/sil)
         sophism_demod_adhoc,demodima_teo,info.sscene_apod,info.invercont
         sxaddpar,headth,"HISTORY",'Adhoc corrected'
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_teo_xtalk_'+strtrim(obs,2)+'.fits',demodima_teo,headth,/compress $
         else writefits,info.saves+info.files(progind)+'_teo_xtalk_'+strtrim(obs,2)+'.fits',demodima_teo,headth
      endif 

      if ((size(dmodm))(0) gt 2) then begin
         demodimaav=readfits(info.saves+info.files(progind)+'_av_'+strtrim(obs,2)+'.fits*',heador,/comp,/sil)
         sophism_demod_adhoc,demodimaav,info.sscene_apod,info.invercont
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_av_xtalk_'+strtrim(obs,2)+'.fits',demodimaav,heador,/compress $
         else writefits,info.saves+info.files(progind)+'_av_xtalk_'+strtrim(obs,2)+'.fits',demodimaav,heador
         if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
            demodimaav_teo=readfits(info.saves+info.files(progind)+'_teo_av_'+strtrim(obs,2)+'.fits*',headorth,/comp,/sil)
            sophism_demod_adhoc,demodimaav_teo,info.sscene_apod,info.invercont
            sxaddpar,headth,"HISTORY",'Adhoc corrected'
            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_teo_av_xtalk_'+strtrim(obs,2)+'.fits',demodimaav_teo,headorth,/compress $
            else writefits,info.saves+info.files(progind)+'_teo_av_xtalk_'+strtrim(obs,2)+'.fits',demodimaav_teo,headorth
         endif       
         
      endif 


   endfor ;obs

endif ;adhoc


END
