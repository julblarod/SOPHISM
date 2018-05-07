pro sophism_accu

; ==================================================================
; DESCRIPTION
;    Performs image accumulation on detector by different ways
;    depending on which modules (polarization/filtergraph) have been
;    used in the simulation.
;    The observational scheme followed is that of first modulation
;    through LCVRs and then wavelength change (when all the modules
;    have been run).
;
; CALLED FROM 
;    sophism
;
; SUBROUTINES
;    headfits
;
; MAIN INPUTS
;
;
; OUTPUT
;    Accumulated data.
;
; VERSION HISTORY
;   Alex Feller
;   J. Blanco. Feb 2012. v1_0. Modified to fit into sophism. 
;   J. Blanco. March 2012. v2_0. Readaptation of accumulation loops
;      to allow different simulation options and rearrangement of
;      summations to fit modulation schemes. Dual-beam mode
;   J. Blanco. Dec 2012. v2_1. Printings and displays for talk option
;   J. Blanco. Jul 2014. v2_2. Corrected ntim according to new cycrep
;      scheme, and cycrep in other calculations. Added obs and
;      according things for the different observation sets.
;   J. Blanco. Sep 2014. v2_3. Corrected bug in ntim for case of
;      polarization but no etalon.
;   J. Blanco May 2015. v2_4. Corrected bug in the progma variable
;      when loading input data directly
;
; ==================================================================

print,''
print,'sophism_accu'
print,''

restore,'../settings/settings_sophism.sav'
progind=7

; number of time exposures within one modulation period
nperiod=info.modnbuff*info.nacc
;nmcycl=info.cycrep

; number of time exposures
;ntim=nperiod*nmcycl*info.etalpos
ntim=nperiod*info.etalpos

if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if (progma eq -1) then progma=0

; ------------------------------------------------------------------------------
; Accumulate
; ------------------------------------------------------------------------------

;initialize counters
j=0     
k=0
istate_old = 0
ilambda_old=0
iobs_old=0
obs=0
kd=0

;determine if Polarization and Etalon modules have been used by
;looking on data header.
;This will divert the program through one of the accumulation paths.
head=headfits(info.saves+info.files(progma)+'_0.fits*')
headmod=sxpar(head,'HISTORY')
etal=where(headmod eq 'Etalon Module')
poli=where(headmod eq 'Polarization Module')


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;first case: Neither etalon module nor polarization module run
if (etal eq -1) then begin

   if (poli eq -1) then begin
 
      if (info.talk eq 1) then begin
         print,'No Etalon & No Modulation'
         print,'i.e., simple accumulation of original Stokes'
      endif 
;      ntim=info.nacc*nmcycl
;here, with no etalon or polarization, obsset must be considered with ntim
      ntim=info.nacc*info.obsset;*nmcycl
      for nf=0,ntim-1 do begin
         ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
;initialize the output variable
         if (nf eq 0) then a=ima*0.

;save in different outputs the different observational sets. The last
;one is out of the loop
         iobs=floor(nf/(ntim/info.obsset))
         if (iobs gt iobs_old) then begin
            sxaddpar,head,"HISTORY",'Accumulation Module'
            sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',a,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',a,head
            iobs_old=iobs
            obs=obs+1
            a=a*0.
          endif ;iobs

;simple adding of each file. No modulation or wavelength selection so
;no change of buffers of lambda
         a[*,*,*,*]=a[*,*,*,*]+ima[*,*,*,*]
         if (info.talk eq 1) then print,format='(18(%"\b"),"Accumulating ",I3," %",$)',float(nf+1)/ntim*100.
      endfor 
      print,''

      sxaddpar,head,"HISTORY",'Accumulation Module'
;      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',a,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',a,head
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',a,head,/compress $
else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',a,head

   endif else begin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;second case: Polarization module run but no spectal module
      if (info.talk eq 1) then begin
         print,'No Etalon, but Modulation'
         print,'i.e., accumulation of modulated intensities, in temporal order'
      endif 

;      ntim=nperiod*nmcycl
      ntim=nperiod*info.obsset
      for nf=0,ntim-1 do begin
         ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
;initialize the accumulation variable (a) and the output variable
;(dmodbuff) modifying the stokes dimension (modulated data have only
;one element)
         if (nf eq 0) then begin
            a=ima*0.
            sz = size(a)
            dmodbuff = fltarr(sz[1],info.modnbuff,sz[3],sz[4])
;not used etalon module, so probably not prepared the obsset
;there. Think about.
;            dmodbuff = fltarr(sz[1]/info.obsset,info.modnbuff,sz[3],sz[4])
         endif 

;save in different outputs the different observational sets. The last
;one is out of the loop
         iobs=floor(nf/(ntim/info.obsset))
         if (iobs gt iobs_old) then begin
            sxaddpar,head,"HISTORY",'Accumulation Module'
            sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head
            iobs_old=iobs
            obs=obs+1
            a=a*0.
            dmodbuff=dmodbuff*0.
         endif ;iobs

;accumulate for each modulation state
         a[*,*,*,*]=a[*,*,*,*]+ima[*,*,*,*]
         istate=floor((nf+1.)/info.nacc)
;when reached the total accumulations for one state, 'save' it in
;output variable, adding to previous existing data in case the
;polarization cycle is repeated
         if (istate gt istate_old) then begin
            istate_old = istate
            dmodbuff[*,j,*,*] += a[*,*,*,*]
;re-initialize accumulation variable
            a = ima*0.

;change modulation buffers and, if more than one cycle, repeat from
;the first one
            j=(j+1) mod info.modnbuff

         endif ;istate
         if (info.talk eq 1) then print,format='(31(%"\b"),"Accumulating. State ",I3,"  ",I3," %",$)',istate mod info.modnbuff,float(nf+1)/ntim*100.
      endfor 
      print,''

      sxaddpar,head,"HISTORY",'Accumulation Module'
;      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',dmodbuff,head
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress $
      else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head

    endelse ;no etalon, yes modul
endif else begin 
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;third case: etalon module run but no polarization module
   if (poli eq -1) then begin

      if (info.talk eq 1) then begin
         print,'Etalon but No Modulation'
         print,'i.e., simple accumulation of original Stokes, in temporal order of etalon positions'
      endif 

;      ntim=info.nacc*nmcycl*info.etalpos
;without polarization, it would have no meaning the cycles (unless for
;fast-wavelength mode, if that is developed), no changes needed
      ntim=info.nacc*info.etalpos
      for nf=0,ntim-1 do begin
         ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
;initialize accumulation (a) and output (dmodbuff) variables. Stokes
;dimension obtained directly from read data, in case data with only I
;(or other option) is used
         if (nf eq 0) then begin
            sz = size(ima)
            a=fltarr(sz[2],sz[3],sz[4])
            dmodbuff = fltarr(sz[1]/info.obsset,sz[2],sz[3],sz[4])
         endif 

;save in different outputs the different observational sets. The last
;one is out of the loop
         iobs=floor(nf/(ntim/info.obsset))
         if (iobs gt iobs_old) then begin
            sxaddpar,head,"HISTORY",'Accumulation Module'
            sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head
            iobs_old=iobs
            obs=obs+1
            a=a*0.
            dmodbuff=dmodbuff*0.
         endif ;iobs

;accumulate for wavelength position
         a=a+reform(ima[j,*,*,*])
         ilambda=floor((nf+1.)/(info.nacc));*nmcycl))
;when finished with the given position, 'save' in output variable
        if (ilambda gt ilambda_old) then begin
            ilambda_old = ilambda
;            dmodbuff[j,*,*,*] = dmodbuff[j,*,*,*]+ a
            dmodbuff[kd,*,*,*] = dmodbuff[kd,*,*,*]+ a
;reinitialize accumulation variable
            a = a*0.
;change wavelength position element in output
            j=(j+1)
            kd=j mod (sz[1]/info.obsset)
         endif ;ilambda
         if (info.talk eq 1) then print,format='(35(%"\b"),"Accumulating. Wavelength ",I2,"  ",I3," %",$)',ilambda,float(nf+1)/ntim*100.
      endfor 
      print,''

      sxaddpar,head,"HISTORY",'Accumulation Module'
;      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',dmodbuff,head
      sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress $
else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head

   endif else begin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
;fourth case: both spectral and polarization enabled 
      if (info.talk eq 1) then begin
         print,'Etalon & Modulation'
         print,'i.e., accumulation of modulated intensities, in temporal order of LCVR and etalon positions'
      endif 

      for nf=0,ntim-1 do begin
         ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
;initialize accumulation and ouput variables.
         if (nf eq 0) then begin
            sz = size(ima)
            a=fltarr(sz[3],sz[4])  
            dmodbuff = fltarr(sz[1]/info.obsset,info.modnbuff,sz[3],sz[4])      
         endif 

;save in different outputs the different observational sets. The last
;one is out of the loop
         iobs=floor(nf/(ntim/info.obsset))
         if (iobs gt iobs_old) then begin
            sxaddpar,head,"HISTORY",'Accumulation Module'
            sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
            if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head
            iobs_old=iobs
            obs=obs+1
            a=a*0.
            dmodbuff=dmodbuff*0.
         endif ;iobs

;given the default observational scheme, first set a wavelength element
;         ilambda=floor(nf/(nperiod*nmcycl))
;la antigua opcion justo encima, con el nmcycl, podria ser util en
;caso de querer sumar los ciclos directamente
         ilambda=floor(nf/(nperiod))
         if (ilambda gt ilambda_old) then begin
            ilambda_old=ilambda
            k=(k+1)
         endif ;ilambda         

; accumulate on detector (Stokes dimension is 0 becuase it has been modulated)
         a=a+reform(ima[k,0,*,*]) 
      
; if the end of one modulation period is reached (when all
; accumulations are done)
         istate=floor((nf+1.)/info.nacc)
         if istate gt istate_old then begin
            istate_old = istate
; in case of obsset gt 1, k must be mod'ed for dmodbuff, so ima
; can be read to next states but dmodbuff kept in its dimensions
            kd=k mod (sz[1]/info.obsset)
     
;oookeeey.Esto tiene sentido en caso de N repeticiones para una misma
;longitud de onda (x ejemplo) del ciclo de modulacion.La diferencia de
;Juanjo entre nA y N.
            dmodbuff[kd,j,*,*] =dmodbuff[kd,j,*,*]+ a
;re-initialize accumulation variable
            ;a = ima*0. 
            a=a*0.
;to change mod buffers and, if more than one cycle, go to the first again
            j=(j+1) mod info.modnbuff

         endif ;istate
         if (info.talk eq 1) then print,format='(52(%"\b"),"Accumulating. Wavelength ",I2," State  ",I3,"  Total ",I3," %",$)',ilambda,istate mod info.modnbuff,float(nf+1)/ntim*100.
      endfor ;nf
      print,''
     
      sxaddpar,head,"HISTORY",'Accumulation Module'
      sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
;      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',dmodbuff,head
;      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+'0.fits',dmodbuff,head,/compress $ ;+strtrim(obs,2)+'.fits',dmodbuff,head,/compress $
;else writefits,info.saves+info.files(progind)+'_'+'0.fits',dmodbuff,head ;+strtrim(obs,2)+'.fits',dmodbuff,head
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(obs,2)+'.fits',dmodbuff,head

   endelse ;etalon&modulation

endelse ;etalon






;#####################################
;Dual-beam
;#####################################

if (info.dualbeam eq 1) then begin

;re-initialize counters
   j=0     
   k=0
   istate_old = 0
   ilambda_old=0
   iobs_old=0
   obs=0
   kd=0

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (etal eq -1) then begin

      if (poli eq -1) then begin

 ;'No Etalon & No Modulation'

         ntim=info.nacc*nmcycl
         for nf=0,ntim-1 do begin
            ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
            if (nf eq 0) then a=ima*0.
            
            a[*,*,*,*]=a[*,*,*,*]+ima[*,*,*,*]
            if (info.talk eq 1) then print,format='(23(%"\b"),"DUAL Accumulating ",I3," %",$)',float(nf+1)/ntim*100.
         endfor 
         print,''

         sxaddpar,head,"HISTORY",'Accumulation Module'
;         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',a,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',a,head
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',a,head,/compress $
else writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',a,head
         
      endif else begin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ;'No Etalon, but Modulation'

;         ntim=nperiod*nmcycl
         ntim=nperiod*info.obsset
         for nf=0,ntim-1 do begin
            ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
            if (nf eq 0) then begin
               a=ima*0.
               sz = size(a)
               dmodbuff = fltarr(sz[1],info.modnbuff,sz[3],sz[4])
            endif 

            iobs=floor(nf/(ntim/info.obsset))
            if (iobs gt iobs_old) then begin
               sxaddpar,head,"HISTORY",'Accumulation Module'
               sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
               if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head
               iobs_old=iobs
               obs=obs+1
               a=a*0.
               dmodbuff=dmodbuff*0.
            endif ;iobs

            a[*,*,*,*]=a[*,*,*,*]+ima[*,*,*,*]
            istate=floor((nf+1.)/info.nacc)
            if (istate gt istate_old) then begin
               istate_old = istate
               dmodbuff[*,j,*,*] += a[*,*,*,*]
               a = ima*0.

               j=(j+1) mod info.modnbuff

            endif ;istate
            if (info.talk eq 1) then print,format='(36(%"\b"),"DUAL Accumulating. State ",I3,"  ",I3," %",$)',istate mod info.modnbuff,float(nf+1)/ntim*100.
         endfor 
         print,''

         sxaddpar,head,"HISTORY",'Accumulation Module'
;         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',dmodbuff,head
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress $
else writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head

      endelse ;no etalon, yes modul
   endif else begin 
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (poli eq -1) then begin

;'Etalon but No Modulation'
         
         ntim=info.nacc*nmcycl*info.etalpos
         for nf=0,ntim-1 do begin
            ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
            if (nf eq 0) then begin
               sz = size(ima)
               a=fltarr(sz[2],sz[3],sz[4])
               dmodbuff = fltarr(sz[1],sz[2],sz[3],sz[4])
            endif 

            a=a+reform(ima[j,*,*,*])
            ilambda=floor((nf+1.)/(info.nacc*nmcycl))
 
            if (ilambda gt ilambda_old) then begin
               ilambda_old = ilambda

               dmodbuff[j,*,*,*] =dmodbuff[j,*,*,*]+ a
               a = a*0.

               j=(j+1)
            endif ;ilambda
            if (info.talk eq 1) then print,format='(40(%"\b"),"DUAL Accumulating. Wavelength ",I2,"  ",I3," %",$)',ilambda,float(nf+1)/ntim*100.
         endfor 
         print,''

         sxaddpar,head,"HISTORY",'Accumulation Module'
;         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',dmodbuff,head
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress $
else writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head

      endif else begin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
 ;Etalon & Modulation

         for nf=0,ntim-1 do begin
            ima=readfits(info.saves+info.files(progma)+'_2_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
            if (nf eq 0) then begin
               sz = size(ima)
               a=fltarr(sz[3],sz[4])
               dmodbuff = fltarr(sz[1]/info.obsset,info.modnbuff,sz[3],sz[4]) 
            endif 

            iobs=floor(nf/(ntim/info.obsset))
            if (iobs gt iobs_old) then begin
               sxaddpar,head,"HISTORY",'Accumulation Module'
               sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
               if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head
               iobs_old=iobs
               obs=obs+1
               a=a*0.
               dmodbuff=dmodbuff*0.
            endif ;iobs

            ilambda=floor(nf/(nperiod));*nmcycl))
            if (ilambda gt ilambda_old) then begin
               ilambda_old=ilambda
               k=(k+1)
            endif ;ilambda         

            a=a+reform(ima[k,0,*,*])
      
            istate=floor((nf+1.)/info.nacc)
            if istate gt istate_old then begin
               istate_old = istate
               kd=k mod (sz[1]/info.obsset)

               dmodbuff[kd,j,*,*] =dmodbuff[kd,j,*,*]+ a
;               a = ima*0. 
               a=a*0.

               j=(j+1) mod info.modnbuff

            endif ;istate
         if (info.talk eq 1) then print,format='(57(%"\b"),"DUAL Accumulating. Wavelength ",I2," State  ",I3,"  Total ",I3," %",$)',ilambda,istate mod info.modnbuff,float(nf+1)/ntim*100.
         endfor ;nf
         print,''

         sxaddpar,head,"HISTORY",'Accumulation Module'
         sxaddpar,head,'NACC',info.nacc,'Number of accumulations'
;         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',dmodbuff,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf,2)+'.fits',dmodbuff,head
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head,/compress $
else writefits,info.saves+info.files(progind)+'_2_'+strtrim(obs,2)+'.fits',dmodbuff,head

      endelse ;etalon&modulation

   endelse ;etalon

endif ;dualbeam


end
