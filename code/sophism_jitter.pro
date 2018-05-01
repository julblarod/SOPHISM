pro sophism_jitter

; ==============================================================================
;
; DESCRIPTION
;    Generates jittering effect (different options) and applies to the
;       data. Image Stabilization System is also possible.
; 
; CALLED FROM 
;    sophism
;
; SUBROUTINES
;    frac_shift
;    sophism_jitter_tpfilt
;    sophism_jitter_atten
;
; MAIN INPUTS
;    jittype        Jittering generation type: ramdon or white noise
;    jitrms         Jittering rms [arcsec]
;    filrms         Flag to renormalize jittering curve after applying
;                      the filtering
;    isson          Flag to enable the image stabilization system
;    jitoffset      A fixed offset to add to the jittering [arcsec]
;
; OUTPUT
;    Data shifted by the effect of jittering.
;    Curves of jittering, jittering filter and attenuation (if enabled)
;
; VERSION HISTORY 
;    Alex Feller          2010-09-28
;    J. Blanco. 2011. v1_0. Adapted into sophism.
;    J. Blanco Oct.2012. v1_1. Change from interpolate to frac_shift when
;       shifting the data.
;    J. Blanco. Dec 2012. v1_2. Plottings and printings for talk case
;    J. Blanco. Dec 2012. v2_0. Change frac_shitf to fft_shift for
;       some cases. Correct mistake in renormalization
;    A. Feller. Feb 2013. v3_0. Modifications to normalization in
;       Hinode filtering.
;    J. Blanco. Jun 2014. v3_1. Added ISSRMS parameter to header. 
;    J. Blanco. Jun 2016. v3_2. Included start at some given sample
;       option (nstart)
; 
; ==============================================================================

print,''
print,'sophism_jitter'
print,''

restore,'../settings/settings_sophism.sav'
progind=1

; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------

; number of time samples
ntim = info.ntim

;????????
if (info.startmed eq 1 AND min(where(info.routines eq 1)) eq progind) then goto,startmit
;????????

jitt = {t: 0.0, x: 0.0, y: 0.0}
jitt = replicate(jitt, ntim)
jitt.t = indgen(ntim) * info.tsamp

; RNG seed
seed = systime(/s)
seed = long((seed-long(seed))*1.e6)

;  FFT of jitter 
fjitt = {nu: 0.0, fx: complex(0, 0), fy: complex(0,0)}
fjitt = replicate(fjitt, ntim)
fjitt.nu = sophism_make_freq(ntim, 1./info.tobs)

if (info.talk eq 1) then print,'Jitter generation type: '+info.jittype

; ------------------------------------------------------------------------------
; randomn
; ------------------------------------------------------------------------------

if info.jittype eq 'randomn' then begin

; Draw random numbers
   rarr = randomn(seed, ntim, 2)

   for i=1,2 do begin
      jitt.(i) = reform(rarr[*,i-1])
 
; normalize to rms jitter and set mean exactly to 0
      jitt.(i) = (jitt.(i) / stddev(jitt.(i))) * info.jitrms 
      jitt.(i) = jitt.(i) - mean(jitt.(i))

      fjitt.(i) = fft(jitt.(i), 1)
   endfor
   
endif

; ------------------------------------------------------------------------------
; White noise
; ------------------------------------------------------------------------------

if info.jittype eq 'white noise' then begin

; Draw random phase
   if ntim mod 2 eq 0 then nfreq = ntim/2-1 else nfreq = (ntim-1)/2
   rph = randomu(seed, nfreq, 2) * 2 * !pi

; Generate white noise
   for i=1,2 do begin
      posf = reform(complex(cos(rph[*,i-1]), sin(rph[*,i-1])))
      negf = reverse(posf)
      negf = complex(real_part(negf), -imaginary(negf))
      if ntim mod 2 eq 0 then $
         fjitt.(i) = [complex(0,0), posf, complex(1.0, 0), negf] $
      else $
         fjitt.(i) = [complex(0,0), posf, negf]
; normalize to rms jitter [arcsec]
      fjitt.(i) = fjitt.(i) * sqrt(ntim) * info.jitrms
   endfor

endif

; ------------------------------------------------------------------------------
; Remember mean power density for later normalization after filtering
; ------------------------------------------------------------------------------

powdens=fltarr(2)
for i=1,2 do powdens[i-1]=mean(abs(fjitt[1:ntim/2].(i)))^2.

; ------------------------------------------------------------------------------
; Prepare and apply the frequency filter
; ------------------------------------------------------------------------------

if (info.talk eq 1) then print,'Frequency filtering with type '+info.filtype

;generate the filter
sophism_jitter_tpfilt
restore, info.saves+info.files(progind)+'_ftpfilt.sav'

;print,'reading Johanns Astrium amplitude spectrum'
;print,'In tsamp=20 ms, 5 Hz empieza en indice 156'
;lee=read_ascii('../data/sc_pdr_jitter_spectrum_interpolated.txt',data_st=8)
;maxi=where(lee.field1(0,*) eq round(max(fjitt.nu)))
;mimi=min(abs(fjitt.nu-lee.field1(0,0)),mini)
;powop=fltarr(ntim/2+1)
;powop(mini:*)=interpol(lee.field1(1,0:maxi),lee.field1(0,0:maxi),fjitt(mini:ntim/2).nu)
;set to 0 everything over 100 Hz
;powop(min(where(fjitt.nu ge 100)):*)=0.
;powo=[powop,rotate(powop(1:ntim/2-1),2)]
;info.filtype='Otro'
;stop
;for i=1,2 do fjitt.(i) = fjitt.(i) * powo

for i=1,2 do fjitt.(i) = fjitt.(i) * sqrt(tpfilt.pow)

; ------------------------------------------------------------------------------
; Normalize to required rms
; ------------------------------------------------------------------------------

; Normalize rms, if required
if info.filrms gt 0 then begin
   if info.filtype eq 'Hinode' then begin
; The normalization is such that the integral of the filtered jitter
; power spectrum between cut and infinity is equal to rms^2. Below
; cut, the power density is set to 0.
      if (info.talk eq 1) then print,'Renormalizing int(power, '+string(format='(F4.2)',info.filcuth)+', infinity) to '+ string(format='(F4.2)', info.jitrms)+ '^2'
      a = ((-info.jitrms^2*(p[1]+1)) / (info.filcuth^(p[1]+1))) * (1./exp(p[0]))
      for i=1,2 do fjitt.(i) *= sqrt(a * (1./powdens[i-1]) * (ntim/(2*info.tsamp)))
   endif else begin
      if (info.talk eq 1) then print,'Renormalizing rms to ' + string(format='(F4.1)', info.jitrms)+ ' arcsec'
      for i=1,2 do fjitt.(i) = fjitt.(i)*info.jitrms*ntim / sqrt(real_part(total(conj(fjitt.(i))*fjitt.(i))))
   endelse
endif

; ------------------------------------------------------------------------------
; Apply frequency attenuation curve
; ------------------------------------------------------------------------------

if (info.isson eq 1) then begin
      if (info.talk eq 1) then print,'Applying attenuation curve: ' + info.issatten
      sophism_jitter_atten
      restore, info.saves+info.files(1)+'_iss.sav'
      w = where(attenarr.f gt 0)
      for i=1,2 do fjitt.(i) = fjitt[w].(i) / attenarr[w].f
endif

; ------------------------------------------------------------------------------
; Transform back to time domain and normalize rms
; ------------------------------------------------------------------------------

for i=1,2 do jitt.(i) = real_part(fft(fjitt.(i),-1))

if (info.talk eq 1 AND info.isson eq 1) then print,'Jitter rms after ISS: ' + string(format='(F6.4)',stddev(jitt.x) )+ ' arcsec' 

; ------------------------------------------------------------------------------
; Add offset
; ------------------------------------------------------------------------------

if (info.talk eq 1 AND max(abs(info.jitoffset)) gt 0) then print,format='("Adding offset: ",F5.3,", ",F5.3," arcsec")',info.jitoffset
for i=1,2 do jitt.(i) = jitt.(i) + info.jitoffset[i-1] 

; ------------------------------------------------------------------------------
; Finalize
; ------------------------------------------------------------------------------

save,jitt,filename=info.saves+info.files(progind)+'.sav'

;?????????
;start at some given sample with startmed
startmit:
if (info.startmed eq 1 AND min(where(info.routines eq 1)) eq progind) then begin
   restore,info.saves+info.files(progind)+'.sav'
   nstart=info.nstart
endif else nstart=0
;?????????

if (max(info.progma) gt 0) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if (progma eq -1) then progma=0

prg0=0.
prg=prg0

for nf=nstart,ntim-1 do begin
;for nf=0,ntim-1 do begin
   ima=readfits(info.saves+info.files(progma)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
   sizy=size(ima)
   imajitt=ima*0.
;apply jittering by shitfing the images, scaling with the spatial
;sampling of the data
;two cases for the actual shifting of the data. If spatial dimensions
;are even (they should be periodic also), use fft_shift which gives
;better results. If not, use frac_shift.
   if (odd(sizy(3)) eq 1 OR odd(sizy(4)) eq 1) then begin
      for sto=0,sizy(2)-1 do begin
         for ww=0,sizy(1)-1 do begin
            imajitt[ww,sto,*,*]=frac_shift(reform(ima[ww,sto,*,*]),jitt[nf].x/info.sssamp,jitt[nf].y/info.sssamp)
            prg =round(100.*float(ww+1)*float(sto+1)*(nf+1)/(sizy(1)*sizy(2)*info.ntim))
            if (info.talk eq 1 AND prg gt prg0) then begin
               print,format='(24(%"\b"),"Applying jittering ",I3," %",$)', prg
               prg0=prg
            endif 
         endfor ; ww
      endfor ; sto
   endif else begin
      for sto=0,sizy(2)-1 do begin
         for ww=0,sizy(1)-1 do begin
            imajitt[ww,sto,*,*]=fft_shift(reform(ima[ww,sto,*,*]),[jitt[nf].x/info.sssamp,jitt[nf].y/info.sssamp])
            prg =round(100.*float(ww+1)*float(sto+1)*(nf+1)/(sizy(1)*sizy(2)*info.ntim))
            if (info.talk eq 1 AND prg gt prg0) then begin
               print,format='(24(%"\b"),"Applying jittering ",I3," %",$)', prg
               prg0=prg
            endif               
         endfor ; ww
      endfor ; sto
   endelse 

   sxaddpar,head,"HISTORY",'Jittering Module'
   sxaddpar,head,"JITTYPE",info.jittype,'Jittering Generation'
   sxaddpar,head,"FILTYPE",info.filtype,'Jittering Filter Type'
   sxaddpar,head,"JITRMS",format='(F6.3)',info.jitrms,'Jittering rms [arcsec]'
   if (info.isson eq 1) then begin
      sxaddpar,head,"HISTORY",'ISS Mode ON'
      sxaddpar,head,"ISSRMS",format='(F6.4)',stddev(jitt.x),'Jitter rms after ISS [arcsec]'
   endif else sxaddpar,head,"HISTORY",'ISS Mode OFF'
   if (info.compress eq 1) then writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',imajitt,head,/compress else writefits,info.saves+info.files(progind)+'_'+strtrim(nf,2)+'.fits',imajitt,head
endfor ;nf
print,''

;plot jittering
if (info.talk eq 1) then begin
   !p.multi = [0, 1, 2]
   window,progind,title='Jittering'
   xtitle = 'time [s]'
   plot, jitt.t, jitt.x, xtitle=xtitle, ytitle='x jitter [arcsec]',psy=-2
   plot, jitt.t, jitt.y, xtitle=xtitle, ytitle='y jitter [arcsec]',psy=-2
   !p.multi=0
endif


end

