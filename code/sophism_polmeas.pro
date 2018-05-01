pro sophism_polmeas

; ==============================================================================
;
; DESCRIPTION
;    Simulates polarisation measurement, including the generation of
;       modulation and demodulation matrixes according to a given
;       scheme of parameters for the LCVRs and Mueller matrix of the
;       system, the effect of temporal change of LCVRs' states
;       and the actual modulation of the data.
;
; CALLED FROM 
;    sophism
;
; SUBROUTINES
;    sophism_polmeas_modscheme
;    sophism_polmeas_modscheme_retardance
;
; MAIN INPUTS
;    nstate       Number of time samples in a modulation
;                    state. Function of number of accumulations,
;                    exposure time and sample time
;    nperiod      Number of samples in a complete modulation period
;    cycrep       Repetitions of complete modulation cycle. Instead of
;                    doing all the accumulation in each modulation
;                    state, a lower number may be used and then repeat
;                    the cycle
;    tdeaths      Change times between the states of the LCVRs. Data
;                    'falling' into that time window will be discarded
;
; OUTPUT
;    Modulated images and demodulation matrixes
;
; VERSION HISTORY 
;    Alex Feller   v0.1    2010-09-21
;    J. Blanco. 2011. v1_0. Rearrange data and variable definitions to
;       fit into sophism. Add 'real' modulation case
;    J. Blanco. 2011. v2_0. Include change times of LCVRs and
;       corresponding modulation. Rearrange sample variables 
;    J. Blanco. Dec 2012. v2_1. Add plots and prints for talk case
;    J. Blanco. Dec 2012. v2_2. 'Cosmetic' correction of xy.omargin
;    J. Blanco. Jan 2014. v3_0. Change scheme idea for
;       modulation. Done with modulation matrixes, not 'functions' as
;       before. Take into account possible differences of tdeaths due
;       to tdeathetal. Flexibility for 2D modulation implementation.
;    J. Blanco. Jun 2014. v3_1. Added lambda dependence to modulation
;       matrix
;    J. Blanco. Aug 2014. v3_2. Corrected bug in dual-beam after last
;       modifications 
;    J. Blanco. Sep 2014. v3_3. Change scheme for cycle repetitions
;    J. Blanco. Jun 2016. v3_4. Included start at some given sample
;       option (nstart)
;    J. Blanco. Apr 2017. v3_41. Corrected calculations for case of no
;       time series data
;
; ==============================================================================


print,''
print, 'sophism_polmeas'
print,''
restore,'../settings/settings_sophism.sav'
progind=2

if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if progma eq -1 then progma=0

;?????
if (info.startmed eq 1 AND min(where(info.routines eq 1)) eq progind) then begin
   nstart=0
   head=headfits(info.saves+info.files(progma)+'_0.fits*')
   nlam=fix(sxpar(head,'NAXIS1'))
stop
endif else begin
;create the modulation scheme
   sophism_polmeas_modscheme
   nstart=0
endelse
;??????

; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------
  
restore,info.saves+info.files(progind)+'_scheme.sav'
sizm=size(modmat)

; number of time samples in one state
nstate=info.nstate 
;times of LCVRs states change and number of samples 'lost' in the changes
;tdeaths=[info.tdeath1,info.tdeath2,info.tdeath3,info.tdeath4]
;ndeaths=round(tdeaths*1e-3/info.tsamp)

; number of samples in one modulation period
;nperiod=info.nperiod
; number of cycles
;nmcycl=info.cycrep
; number of samples in one observation run
if (info.dataser eq 0) then ntim=1 else ntim=info.ntim

; ------------------------------------------------------------------------------
; Modulate
; ------------------------------------------------------------------------------

print, 'Modulation type: '+info.modtype

;initialize counter to take into account elimination of data during
;LCVRs changes
nout=0

;------
;set state, modulation and wavelength counters to 0
nst=0
modind=0
modind2=0
lamind=0
;no estoy del todo seguro si el total(ndeaths) es correcto
tim=fltarr(ntim-total(info.ndeaths))
modf=fltarr(info.modnbuff/info.cycrep,ntim-total(info.ndeaths))
;------

;????
if (info.startmed eq 1 AND min(where(info.routines eq 1)) eq progind) then begin
   for nf=0,info.nstart-1 do begin
      if (nf eq nst+nstate+info.ndeaths(modind,lamind)) then begin
         nst=nst+nstate+info.ndeaths(modind,lamind)
         modind=modind+1
         modind2=(modind2+1) mod (info.modnbuff/info.cycrep)
         if (modind eq info.modnbuff) then begin
            modind=0
            lamind=lamind+1
         endif 
      endif 
      if (nf lt nst+nstate) then begin
         tim[nf-nout]=info.tsamp*nf
         for ll=0,nlam-1 do begin 
            for mm=0,info.modnbuff/info.cycrep-1 do begin
               case sizm(0) of
                  2: modi=(modmat[mm,modind2])(0)
                  3: begin
                     llm=ll mod sizm(3)
                     modi=(modmat[mm,modind2,llm])(0)
                  end 
                  4: modi=modmat[mm,modind2,*,*]
                  5: begin
                     llm=ll mod sizm(3)
                     modi=modmat[mm,modind2,llm,*,*]
                  end
               endcase 
               modf[mm,nf-nout]=mean(modi)
            endfor 
         endfor 
      endif else begin
         nout+=1 
      endelse
   endfor 
stop
   nstart0=info.nstart
   nout0=nout
   nstart=info.nstart+nout
   while (nstart gt nst+nstate) do begin
      for nf=nstart0,nstart-1 do begin
         if (nf eq nst+nstate+info.ndeaths(modind,lamind)) then begin
            nst=nst+nstate+info.ndeaths(modind,lamind)
            modind=modind+1
            modind2=(modind2+1) mod (info.modnbuff/info.cycrep)
            if (modind eq info.modnbuff) then begin
               modind=0
               lamind=lamind+1
            endif 
         endif 
         if (nf lt nst+nstate) then begin
            tim[nf-nout]=info.tsamp*nf
            for ll=0,nlam-1 do begin 
               for mm=0,info.modnbuff/info.cycrep-1 do begin
                  case sizm(0) of
                     2: modi=(modmat[mm,modind2])(0)
                     3: begin
                        llm=ll mod sizm(3)
                        modi=(modmat[mm,modind2,llm])(0)
                     end 
                     4: modi=modmat[mm,modind2,*,*]
                     5: begin
                        llm=ll mod sizm(3)
                        modi=modmat[mm,modind2,llm,*,*]
                     end
                  endcase 
                  modf[mm,nf-nout]=mean(modi)
               endfor 
            endfor 
         endif else begin
            nout+=1 
         endelse
      endfor
      nstart0=nstart
      nstart=nstart+(nout-nout0)
      nout0=nout
;stop
   endwhile  

endif 
;stop

for nf=nstart,ntim-1 do begin
;?????
   ;----------------
   if (nf eq nst+nstate+info.ndeaths(modind,lamind)) then begin
      nst=nst+nstate+info.ndeaths(modind,lamind)
      modind=modind+1
      modind2=(modind2+1) mod (info.modnbuff/info.cycrep)
      if (modind eq info.modnbuff) then begin
         modind=0
         lamind=lamind+1
;this for cycles?
;         if (lamind eq info.etalpos) then lamind=0           
      endif 
   endif 
   ;----------------
;only do the modulation to the 'non-discarded' data
   if (nf lt nst+nstate) then begin
      ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/compr,/sil)
;time array
      tim[nf-nout]=info.tsamp*nf
;   if (nf eq 0) then begin
      siz=size(ima)
;only 1 in Stokes dimension because each data is a different
;modulation state
      if (info.dataser eq 0) then imamod=fltarr(siz(1),info.modnbuff/info.cycrep,siz(3),siz(4)) else imamod=fltarr(siz(1),1,siz(3),siz(4))
;   endif
;modulate from modulation matrix and
;rearranging already the order of Stokes parameters from input data
;(i.e. index 1 is not Q but V; indexes 2:3 --> Q:U; change to I,Q,U,V)
;   imamod[*,0,*,*]=modf[nf].i*reform(ima(*,0,*,*))+modf[nf].q*reform(ima(*,2,*,*))+modf[nf].u*reform(ima(*,3,*,*))+modf[nf].v*reform(ima(*,1,*,*))
             
;reorder the input data from the I,V,Q,U of MHDsim to the I,Q,U,V
;usual data, if selected. If not, consider that input data are I,Q,U,V
;already
      if (info.stokreorder eq 1) then stokorder=[0,2,3,1] else stokorder=[0,1,2,3]

;in case the mod. matrix is FOV-dependant, must loop through
;wavelength dimension to make the product of the FOV matrix and data.
;Some better way?
;also, must check FOV-dependence. If not, modmat ends up as 1-element
;array which creates problems in the product with the image 
;Now, adding lambda dependence for mod matrix
      for ll=0,siz(1)-1 do begin
;         for mm=0,3 do imamod[ll,0,*,*]=imamod[ll,0,*,*]+modmat[mm,modind,*,*]*reform(ima(ll,stokorder[mm],*,*))
         for mm=0,info.modnbuff/info.cycrep-1 do begin
;            if (sizm(0) eq 4) then modi=modmat[mm,modind,*,*] else modi=(modmat[mm,modind])(0)
            case sizm(0) of
               2: modi=(modmat[mm,modind2])(0)
               3: begin
;in case the dependency on lambda is not as 'long' as the wavelength
;dimension of the data
                  llm=ll mod sizm(3)
                  modi=(modmat[mm,modind2,llm])(0)
                  end 
               4: modi=modmat[mm,modind2,*,*]
               5: begin
                  llm=ll mod sizm(3)
                  modi=modmat[mm,modind2,llm,*,*]
                  end
            endcase 
;stop
            imamod[ll,0,*,*]=imamod[ll,0,*,*]+modi*reform(ima(ll,stokorder[mm],*,*))
;save an average of mod matrix coeff for each time step to plot later
;in case of wavelength dependance, it won't be an average, just
;the last one
            modf[mm,nf-nout]=mean(modi)

         endfor
      endfor

;apply the modulation when no data series. Must include the other
;modulations in the same image array (instead of normal operation)
      if (info.dataser eq 0) then begin
         for modind2=1,info.modnbuff/info.cycrep-1 do begin
            for ll=0,siz(1)-1 do begin
               for mm=0,info.modnbuff/info.cycrep-1 do begin
                  case sizm(0) of
                     2: modi=(modmat[mm,modind2])(0)
                     3: begin
                        llm=ll mod sizm(3)
                        modi=(modmat[mm,modind2,llm])(0)
                     end 
                     4: modi=modmat[mm,modind2,*,*]
                     5: begin
                        llm=ll mod sizm(3)
                        modi=modmat[mm,modind2,llm,*,*]
                     end
                  endcase 
                  imamod[ll,modind2,*,*]=imamod[ll,modind2,*,*]+modi*reform(ima(ll,stokorder[mm],*,*))
               endfor 
            endfor 
         endfor
         modind2=0
      endif  

;only the 'non-death' frames are taken
;   if (max(imamod) ne 0) then begin
      
      sxaddpar,head,"HISTORY",'Polarization Module'
      if (info.dualbeam eq 1) then sxaddpar,head,"HISTORY",'Dual-beam Mode ON'
      sxaddpar,head,"MUELLER",info.mueller,'Mueller Matrix: Identity or Custom'
      for nlc=0,n_elements(info.ndeathlcvr)-1 do sxaddpar,head,"TDEATLC"+strtrim(nlc,2),info.ndeathlcvr(nlc)*info.tsamp,'Time for LCVR change'
      for net=0,n_elements(info.ndeathetal)-1 do sxaddpar,head,"TDEATET"+strtrim(net,2),info.ndeathetal(net)*info.tsamp,'Time for Etalon change'
      sxaddpar,head,"CYCLES",info.cycrep,'Number of Cycle repetitions'
      if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf-nout,2)+'.fits',imamod,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf-nout,2)+'.fits',imamod,head
        

;#################################
;Dual-beam
;#################################
           
      if (info.dualbeam eq 1) then begin
         imamod2=imamod*0.
;         imamod2[*,0,*,*]=modf2[nf].i*reform(ima(*,0,*,*))+modf2[nf].q*reform(ima(*,2,*,*))+modf2[nf].u*reform(ima(*,3,*,*))+modf2[nf].v*reform(ima(*,1,*,*))

        
      for ll=0,siz(1)-1 do begin
;         for mm=0,3 do imamod2[ll,0,*,*]=imamod2[ll,0,*,*]+modmat2[mm,modind,*,*]*reform(ima(ll,stokorder[mm],*,*))
         for mm=0,info.modnbuff/info.cycrep-1 do begin
;            if (sizm(0) eq 4) then modi=modmat2[mm,modind,*,*] else modi=(modmat2[mm,modind])(0)
            case sizm(0) of
               2: modi=(modmat2[mm,modind2])(0)
               3: begin
                  llm=ll mod sizm(3)
                  modi=(modmat2[mm,modind2,llm])(0)
                  end
               4: modi=modmat2[mm,modind2,*,*]
               5: begin
                  llm=ll mod sizm(3)
                  modi=modmat2[mm,modind2,llm,*,*]
                  end
            endcase 
            imamod2[ll,0,*,*]=imamod2[ll,0,*,*]+modi*reform(ima(ll,stokorder[mm],*,*))
         endfor 
      endfor 

      if (info.dataser eq 0) then begin
         for modind2=1,info.modnbuff/info.cycrep-1 do begin
            for ll=0,siz(1)-1 do begin
               for mm=0,info.modnbuff/info.cycrep-1 do begin
                  case sizm(0) of
                     2: modi=(modmat2[mm,modind2])(0)
                     3: begin
                        llm=ll mod sizm(3)
                        modi=(modmat2[mm,modind2,llm])(0)
                     end 
                     4: modi=modmat2[mm,modind2,*,*]
                     5: begin
                        llm=ll mod sizm(3)
                        modi=modmat2[mm,modind2,llm,*,*]
                     end
                  endcase 
                  imamod2[ll,modind2,*,*]=imamod2[ll,modind2,*,*]+modi*reform(ima(ll,stokorder[mm],*,*))
               endfor 
            endfor 
         endfor 
      endif  

;header is the same as above, no need to do anything about it
         if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf-nout,2)+'.fits',imamod2,head,/compress else writefits,info.saves+info.files(progind)+'_2_'+strtrim(nf-nout,2)+'.fits',imamod2,head
      endif ;dualbeam


   endif else begin
;counter to take into account the 'death' images for the file names
      nout+=1 
   endelse  ;imamod
   if (info.talk eq 1) then print,format='(17(%"\b"),"Modulating: ",I3," %",$)',float(nf+1)/float(ntim)*100.
endfor ;nf
print,''

save,filenam=info.saves+info.files(progind)+'_time.sav',tim,modf

;plot modulation scheme
if (info.talk eq 1) then begin
   !p.multi = [0, 1, info.modnbuff/info.cycrep]
   !x.omargin=[0,10]
   !y.omargin=[0,5]
   window,xs=550,ys=750,title='Modulation Scheme',/free
   xtitle = 'time [s]'
   yr = [-1.1,1.1]
;   laten=where(modf.i ne 0)
   
;   delim=[info.nstate,info.nstate+info.ndeaths(0,0)]
;   for dd=2,info.modnbuff do delim=[delim,info.nstate*dd+total(info.ndeaths(0:dd-2,0),1),info.nstate*dd+total(info.ndeaths(0:dd-1,0),1)]
;   delim2=delim
;   for dd=1,info.etalpos-1 do delim2=[delim2,delim+max(delim2)]
;   delim2=delim2*info.tsamp-info.tsamp/2.

prev=0
   delim=[0];[info.nstate,info.nstate+total(info.ndeaths(0,ddet),1)]
   for ddet=0,info.etalpos-1 do begin
      for ddlc=1,info.modnbuff do begin
         delim=[delim,prev+info.nstate,prev+info.nstate+info.ndeaths(ddlc-1,ddet)]
         prev=prev+info.nstate+info.ndeaths(ddlc-1,ddet)
;         delim=[delim,info.nstate*ddlc+total(info.ndeaths(0:ddlc-2,ddet),1),info.nstate*ddlc+total(info.ndeaths(0:ddlc-1,ddet),1)]
      endfor 
   endfor 
   delim=delim[1:*]
   delim=delim*info.tsamp-info.tsamp/2.

   postitle=fltarr(4)
;   postitle(0)=delim2[0]/2.
;   postitle(1)=delim2[1]+(delim2[2]-delim2[1])/2.
;   postitle(2)=delim2[3]+(delim2[4]-delim2[3])/2.
;   postitle(3)=delim2[5]+(delim2[6]-delim2[5])/2.
   postitle(0)=delim[0]/2.
   postitle(1)=delim[1]+(delim[2]-delim[1])/2.
   postitle(2)=delim[3]+(delim[4]-delim[3])/2.
   postitle(3)=delim[5]+(delim[6]-delim[5])/2.
   uptitle=['I!D1!N', 'I!D2!N', 'I!D3!N', 'I!D4!N']
   ytitle = ['I', 'Q', 'U', 'V']

;   for i = 1, info.modnbuff do begin
;      plot, modf(laten).t, modf(laten).(i), yr=yr, /ys, xtitle=xtitle, ytitle=ytitle[i-1],psy=10
;      if (i eq 1) then xyouts,postitle,replicate(1.25,info.modnbuff),uptitle,/data,charsiz=0.8
;      for j=0,n_elements(delim2)-1 do plots,[delim2(j),delim2(j)],yr,lines=1
;   endfor

   for i=0,info.modnbuff/info.cycrep-1 do begin
      plot,tim,modf[i,*],yr=yr,/ys,xtitle=xtitle,ytitle=ytitle[i],psy=10
      if (i eq 0) then xyouts,postitle,replicate(1.25,info.modnbuff/info.cycrep),uptitle,/data,charsiz=0.8
      for j=0,n_elements(delim)-1 do plots,[delim(j),delim(j)],yr,lines=1
   endfor 

   !p.multi=0
   !x.omargin=0
   !y.omargin=0 
endif ;talk


end

