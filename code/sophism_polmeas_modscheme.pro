pro sophism_polmeas_modscheme

; ==============================================================================
;
; DESCRIPTION
;    Define modulation scheme, consisting of modulation and demodulation
;       functions and the demodulation matrix
;
; CALLED FROM 
;    sophism
;
; SUBROUTINES
;    sophism_repeat
;    sophism_polmeas_modscheme_retardance
;
; MAIN INPUTS
;    nstate       Number of time samples for a modulation state,
;                    depending on exposure time and number of accumulations.
;    tdeaths      Times for the changings between the different states
;                    of the LCVRs. During these times, no images are stored.
;    nperiod      Number of time samples in a whole modulation period,
;                    including the later discarded ones from tdeaths.
;    cycrep       Repetitions of the whole modulation cycle.
;    tobs         Total observation time of the simulation run
;    tsamp        Time sampling of the input data
;    mueller      Preparation type for Mueller matrix: Identity, No
;                    errors (real matrix without adding unknown
;                    errors), or Errors (adding unknowns)
;    modtype      Type of modulation to be performed: Longitudinal
;                    (just I and V), Ideal (ideal modulation), Real
;                    (actual modulation with the given parameters for
;                    the LCVRs)
;
; OUTPUT
;    Mueller matrix, modulation functions (scheme for modulating the
;       data based on modulation matrix) for single and dual-beam and
;       demodulation matrixes for single and dual-beam and for correct
;       and LCVR random-error-free
;
; VERSION HISTORY 
;   Alex Feller   v0.1    2010-09-21
;   J. Blanco. 2011. v1_0. Fit into sophism. Add 'real' modulation
;      case and associated variables
;   J. Blanco. 2011. v2_0. Add Move complete Mueller matrix calculation 
;      here to be used also for ideal. Still preliminary calculation though
;   J. Blanco. Jan 2013. v2_1. Added random errors to the Mueller
;      matrix calculation. Correction of demodulation matrix for ideal case
;   J. Blanco. Jan 2014. v3_0. Change scheme idea for modulation. Done
;      with modulation matrixes, not 'functions' as before. Include
;      FOV dependence possibility for matrixes. Correct order of
;      Mueller and modmat in products. Included 'dualmat' to make
;      easier the production of dual-beam matrixes.
;   J. Blanco. Jun 2014. v3_1. Added lambda dependence to modulation
;      matrix. In Ideal case at least.
;   J. Blanco. Sep 2014. v4_0. Option trees included in 'ideal' and
;      'real' cases to discriminate mode depending on Mueller
;      and LCVR's phases dependences with FOV and
;      lambda. Changes in inputs of new modscheme_retardance
;   J. Blanco. Jun 2015. v4_1. Change names of modes for better
;      description. Preparation of longit. ideal.
;   J. Blanco. Feb 2016. v4_2. Introducing birefringence through
;      Mueller matrix
;   J. Blanco. Apr 2016. v4_3. Added theoretical averaged demodulation
;      matrixes for the corresponding cases.
;   J. Blanco. Nov 2016. v4_4. Corrected bug? when retardance has FOV
;      or lambda-dep. Not necessary probably since it should be done
;      in setting file but all the same. Added on-screen progress
;      display for verbose case
;   J. Blanco. Mar 2017. v4_5. Completed production of matrices in the
;      general case (theoretical matrices were missing)
;
; ==============================================================================

restore,'../settings/settings_sophism.sav'
progind=2

; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------
if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if progma eq -1 then progma=0

;all this first part, renaming variables (consequence of former
;versions). Mostly useless. Should be removed

;time samples in one modulation state
nstate=info.nstate

;time samples 'thrown' because of delay in LCVRs changing from one
;state to the other
;tdeaths=[info.tdeath1,info.tdeath2,info.tdeath3,info.tdeath4]
;ndeaths=round(tdeaths*1e-3/info.tsamp)
ndeaths=info.ndeaths

; number of time samples within one modulation period
nperiod=info.nperiod

; number of modulation cycles
nmcycl=info.cycrep

; number of time samples
ntim=info.ntim

; demodulation matrix
;dmodm = fltarr(4, info.modnbuff)
;dmodm2=dmodm

;matrix for changing sign of 'magnetic' stokes for dual-beam
dualmat=fltarr(4,4)+1.
for i=1,3 do dualmat[i,*]=-dualmat[i,*]

; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Mueller
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; mirrors & lenses
;mostly a 'template' this section. To be updated when information is up

mmm1=transpose(mirror(45,ri=0.86,ac=6.39))
mmm2=transpose(mirror(45,ri=0.86,ac=6.39))
mmm3=transpose(mirror(45,ri=0.86,ac=6.39))

;random errors in mueller matrix of .95--1 amplitude ratio, -5--5 deg
;phase difference, -90--90 orientation
mmm_ran1=transpose(device(amr=randomu(seed,1)*0.05+0.95, $
pha=randomu(seed,1)*10.-5.,ang=randomu(seed,1)*180.-90.))
mmm_ran2=transpose(device(amr=randomu(seed,1)*0.05+0.95, $
pha=randomu(seed,1)*10.-5.,ang=randomu(seed,1)*180.-90.))
mmm_ran3=transpose(device(amr=randomu(seed,1)*0.05+0.95, $
pha=randomu(seed,1)*10.-5.,ang=randomu(seed,1)*180.-90.))

;selecting mueller matrix case. Ideal--identity, no errors or random
;errors included
case info.mueller of
   'Identity': mmm=mirror(180,ri=1.0,ac=100.)
   'Real': mmm=mmm3#mmm2#mmm1
;   'Errors': mmm=mmm_ran3#mmm3#mmm_ran2#mmm2#mmm1_ran1#mmm1
   else: begin
           print,'Only Identity or Real matrixes'
           break
         end
endcase

;introducing errors in Mueller matrix
;temporal change. When 'Real' implemented, this should probably change again
;keep Mueller matrix without errors
muellmat_teo=mmm/mmm(0)
if (max(info.muellerrors) eq 1) then begin
   ampdif=abs(info.muellerramp2-info.muellerramp1)
   ampoff=(info.muellerramp1<info.muellerramp2)/float(ampdif)
   signpha=randomu(seed,1)
   if (signpha lt 0.5) then signpha=-1 else signpha=1
   phadif=abs(info.muellerrpha2-info.muellerrpha1)
   phaoff=(info.muellerrpha1<info.muellerrpha2)/float(phadif)
   signang=randomu(seed,1)
   if (signang lt 0.5) then signang=-1 else signang=1
   angdif=abs(info.muellerrang2-info.muellerrang1)
   angoff=(info.muellerrang1<info.muellerrang2)/float(angdif)
   mmm1amperr=(randomu(seed,1)+ampoff)*ampdif*info.muellerrors(0)
   mmm1phaerr=(randomu(seed,1)+phaoff)*phadif*signpha*info.muellerrors(1)
   mmm1angerr=(randomu(seed,1)+angoff)*angdif*signang*info.muellerrors(2)

   mmm_err=transpose(device(amr=1.-mmm1amperr,pha=mmm1phaerr,ang=mmm1angerr))
;this is right? Or is it ##?
   mmm=mmm_err#mmm

   save,filenam=info.saves+info.files(progind)+'_mueller_errors.sav',mmm1amperr,mmm1phaerr,mmm1angerr
endif
muellmat=mmm/mmm(0)


if (info.birefsel ne 'None') then begin
   if (info.birefsel eq 'Load') then begin
      restore,info.saves+info.birefloadname
   endif
   if (info.birefsel eq 'Generate') then begin
      sophism_polmeas_birefringence
      restore,info.saves+info.files(progind)+'_biref_mat.sav'
   endif

;Combine the Mueller of the system with the birefringence
;Have to do it in loop through the spatial dimension of the
;birefringence
;Problem if muellmat has wavelength dependence. Should start cases then...
   muellmatpre=muellmat
   muellmatpreteo=muellmat_teo
   muellmat=matrix_r*0.
   muellmat_teo=matrix_r*0.
;reordering the array dimensions. Supposes only 4 dimensions in
;matrix_r, which should be right
   muellmat=transpose(muellmat,[2,3,0,1])
   muellmat_teo=transpose(muellmat_teo,[2,3,0,1])
   for xx=0,info.sz-1 do begin
      for yy=0,info.sz-1 do begin
;         muellmat=muellmat#mmmbi
         muellmat[*,*,xx,yy]=muellmatpre#reform(matrix_r[xx,yy,*,*])
         muellmat_teo[*,*,xx,yy]=muellmatpreteo#reform(matrix_r[xx,yy,*,*])
      endfor
   endfor
endif
;stop

;***************************
;print,''
;print,'Test para crear matriz mueller 2D'
;print,''
;muellmat=fltarr(4,4,info.sz,info.sz)
;for zz=0,info.sz-1 do begin 
;   for zzp=0,info.sz-1 do begin muellmat[*,*,zzp,zz]=device(amr=1,pha=zzp*25./info.sz,ang=20)
;   endfor 
;endfor 
;print,''
;print,'Test para crear matriz mueller lambda-dep'
;print,''
;heado=headfits(info.saves+info.files(0)+'_0.fits*')
;lambdas=info.etalpos;fix(sxpar(heado,'NAXIS1'))
;solo dependencia con lambda 
;muellmat=fltarr(4,4,lambdas)
;phatot=fltarr(lambdas)
;for zzl=0,lambdas-1 do begin
;   signpha=randomu(seed,1)
;   if (signpha lt 0.5) then signpha=-1 else signpha=1
;   phaa=randomu(seed,1)*100.*signpha
;   muellmat[*,*,zzl]=device(amr=1,pha=phaa,ang=20)
;   phatot(zzl)=phaa
;endfor 
;para tener con dependencia de FOV y de lambda tambien
;muellmat=fltarr(4,4,info.etalpos,info.sz,info.sz)
;phatot=fltarr(info.etalpos,info.sz)
;for zzl=0,info.etalpos-1 do begin
;   phaa=randomu(seed,1)*50.*findgen(info.sz)/info.sz-randomu(seed,1)*10.
;   for zz=0,info.sz-1 do begin 
;       muellmat[*,*,zzl,*,zz]=device(amr=1,pha=phaa,ang=20)
;;      for zzp=0,info.sz-1 do begin 
;;         muellmat[*,*,zzl,zzp,zz]=device(amr=1,pha=zzp*50./info.sz,ang=20)
;;      endfor 
;   endfor 
;   phatot(zzl,*)=phaa
;endfor 
;save,phatot,filename=info.saves+info.files(progind)+'_scheme_phasesmuell.sav'
;***************************

;to check later whether we are in 2D case or not
sizmuel=size(muellmat)

; ------------------------------------------------------------------------------
; Longitudinal I-Q, I-U or I-V modulation/demodulation
; ------------------------------------------------------------------------------

;would be the longitudinal case.
;careful!!! Hasn't been checked in a while. Probably wrong now.

if info.modtype eq 'longit. ideal' then begin

;   print,'Not correctly implemented yet'
;   stop

; demodulation function
;   dmodf = {t: fltarr(ntim), f: fltarr(ntim, info.modnbuff)}
;   dmodf.t = modf.t

; modulation function
;   modf.i = replicate(0.5, ntim)
;   a = fltarr(nperiod)
;   a[0:nperiod/2-1] = 0.5
;   a[nperiod/2:nperiod-1] = -0.5
;   b = 0.0
;   for i=0,nmcycl-1 do b = [b, a]
;   b = b[1:*]
;   modf.(where(tag_names(modf) eq info.modspar)) = b

; demodulation function
;   a = fltarr(nperiod)
;   a[0:nperiod/2-1] = 1
;   a[nperiod/2:nperiod-1] = 0
;   b = 0.0
;   for i=0,nmcycl-1 do b = [b, a]
;   b = b[1:*]
;;dmodf.f[*,0] = b
;;dmodf.f[*,1] = shift(b, nperiod/2)

; demodulation matrix
;   dmodm[0,*] = [1,1]
;   case info.modspar of
;      'Q': dmodm[1,*] = [1,-1]
;      'U': dmodm[2,*] = [1,-1]
;      'V': dmodm[3,*] = [1,-1]
;   endcase

;case according to parameter to modulate
   case info.modspar of 
      'Q': aa=[1,0,0]
      'U': aa=[0,1,0]
      'V': aa=[0,0,1]
   endcase 
   
;write a 4x4 modulation matrix for convenience, although only 2 states
   modmat=[[1,aa],$
           [1,-aa],$
           [0,0,0,0],$
           [0,0,0,0]]

;option-tree for Mueller matrix dependencies (none, FOV, lambda, all)
   case sizmuel(0) of
      2: begin
         dmodm_teo=invert(modmat##muellmat_teo)
         dmodm2_teo=dmodm_teo*transpose(dualmat)
         modmat=modmat##muellmat
         modmat2=modmat*dualmat
         dmodm=invert(modmat)
         dmodm2=invert(modmat2)
         end 

      3: begin
         modmator=modmat
         modmat=muellmat*0.
         modmat2=modmat*0.
;this, as long as we want demodulation matrixes as f(lambda)
         dmodm=muellmat*0.
         dmodm2=dmodm*0.
         dmodm_teo=dmodm*0.
         dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
         for ll=0,sizmuel(3)-1 do begin
               modmat[*,*,ll]=modmator##muellmat[*,*,ll]
               modmat2[*,*,ll]=modmat[*,*,ll]*dualmat
               dmodm[*,*,ll]=invert(modmat[*,*,ll])
               dmodm2[*,*,ll]=invert(modmat2[*,*,ll])
         endfor
;prepare other demodulation matrixes
;Average of the whole FOV
         dmodmav=total(dmodm,3)/float(sizmuel(3))
         dmodmav2=total(dmodm2,3)/float(sizmuel(3))
         end

      4: begin
         modmator=modmat
         modmat=muellmat*0.
         modmat2=modmat*0.
;this, as long as we want 2D demodulation matrixes
         dmodm=muellmat*0.
         dmodm2=dmodm*0.
         dmodm_teo=dmodm*0.
         dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
         for xx=0,sizmuel(3)-1 do begin
            for yy=0,sizmuel(4)-1 do begin
               modmat[*,*,xx,yy]=modmator##muellmat[*,*,xx,yy]
               modmat2[*,*,xx,yy]=modmat[*,*,xx,yy]*dualmat
               dmodm[*,*,xx,yy]=invert(modmat[*,*,xx,yy])
               dmodm2[*,*,xx,yy]=invert(modmat2[*,*,xx,yy])
            endfor 
         endfor
;prepare other demodulation matrixes
;Average of the whole FOV
         dmodmav=total(total(dmodm,3),3)/float(sizmuel(3))/float(sizmuel(4))
         dmodmav2=total(total(dmodm2,3),3)/float(sizmuel(3))/float(sizmuel(4))
;Average to a given number of demod. mats. or in a given size of pixs.
;First, fix the number of demodulation matrixes so everything is
;covered and enters the FOV.
;For now, consider square image/mueller matrix
;     submatsiz=info.sizmuel(3)/float(info.demmats)
;     submat=info.sizmuel(3)/round(submatsiz)
;     dmodmsub=fltarr(4,4,submat,submat)
         end

      5: begin
         modmator=modmat
         modmat=muellmat*0.
         modmat2=modmat*0.
;this, as long as we want 2D demodulation matrixes and lambda dependence
         dmodm=muellmat*0.
         dmodm2=dmodm*0.
         dmodm_teo=dmodm*0.
         dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
         for ll=0,sizmuel(3)-1 do begin
            for xx=0,sizmuel(4)-1 do begin
               for yy=0,sizmuel(5)-1 do begin
                  modmat[*,*,ll,xx,yy]=modmator##muellmat[*,*,ll,xx,yy]
                  modmat2[*,*,ll,xx,yy]=modmat[*,*,ll,xx,yy]*dualmat
                  dmodm[*,*,ll,xx,yy]=invert(modmat[*,*,ll,xx,yy])
                  dmodm2[*,*,ll,xx,yy]=invert(modmat2[*,*,ll,xx,yy])
               endfor 
            endfor
         endfor 
;prepare other demodulation matrixes
;Average of the whole FOV and lambdas
         dmodmav=total(total(total(dmodm,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         dmodmav2=total(total(total(dmodm2,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         end
   endcase

endif ;longit. ideal

; ------------------------------------------------------------------------------
; Real longitudinal modulation/demodulation, single beam & dual-beam
; ------------------------------------------------------------------------------

if info.modtype eq 'longit. real' then begin
   print,'paso a paso, llegara pronto'
   stop



endif ;longit. real

; ------------------------------------------------------------------------------
; Ideal modulation/demodulation, single beam & dual-beam
; ------------------------------------------------------------------------------

if info.modtype eq 'vector. ideal' then begin

;generate the ideal modulation matrix and multiply with Mueller matrix
   a=1./sqrt(3.)
;different modulation orders, produce different output rms in time evolution
;   modmat=[[1,a,a,a],$
;           [1,a,-a,-a],$
;           [1,-a,-a,a],$
;           [1,-a,a,-a]]
   modmat=[[1,a,a,a],$
           [1,a,-a,-a],$
           [1,-a,a,-a],$
           [1,-a,-a,a]]
;   modmat=muellmat##modmat
;   modmat=transpose(modmat)

;considered that ideal modulation has no spatial variation. But
;Mueller matrix may have. Looping.
;   if (sizmuel(0) eq 4) then begin
;      modmator=modmat
;      modmat=muellmat*0.
;      modmat2=modmat*0.
;;      dmodmav=modmator##(total(total(muellmat,3),3)/float(sizmuel(3))/float(sizmuel(4)))
;;      dmodmav2=dmodmav*dualmat
;this, as long as we want 2D demodulation matrixes
;      dmodm=muellmat*0.
;      dmodm2=dmodm*0.
;      dmodm_teo=dmodm*0.
;      dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
;      for xx=0,sizmuel(3)-1 do begin
;         for yy=0,sizmuel(4)-1 do begin
;            modmat[*,*,xx,yy]=modmator##muellmat[*,*,xx,yy]
;            modmat2[*,*,xx,yy]=modmat[*,*,xx,yy]*dualmat
;            dmodm[*,*,xx,yy]=invert(modmat[*,*,xx,yy])
;            dmodm2[*,*,xx,yy]=invert(modmat2[*,*,xx,yy])
;         endfor 
;      endfor
;prepare other demodulation matrixes
;Average of the whole FOV
;      dmodmav=total(total(dmodm,3),3)/float(sizmuel(3))/float(sizmuel(4))
;      dmodmav2=total(total(dmodm2,3),3)/float(sizmuel(3))/float(sizmuel(4))
;Average to a given number of demod. mats. or in a given size of pixs.
;First, fix the number of demodulation matrixes so everything is
;covered and enters the FOV.
;For now, consider square image/mueller matrix
;;     submatsiz=info.sizmuel(3)/float(info.demmats)
;;     submat=info.sizmuel(3)/round(submatsiz)
;;     dmodmsub=fltarr(4,4,submat,submat)
;   endif else begin
;      dmodm_teo=invert(modmat##muellmat_teo)
;      dmodm2_teo=dmodm_teo*transpose(dualmat)
;      modmat=modmat##muellmat
;      modmat2=modmat*dualmat
;      dmodm=invert(modmat)
;      dmodm2=invert(modmat2)
;   endelse 




   case sizmuel(0) of
      2: begin
         dmodm_teo=invert(modmat##muellmat_teo)
         dmodm2_teo=dmodm_teo*transpose(dualmat)
         modmat=modmat##muellmat
         modmat2=modmat*dualmat
         dmodm=invert(modmat)
         dmodm2=invert(modmat2)
         end

      3: begin
         modmator=modmat
         modmat=muellmat*0.
         modmat2=modmat*0.
;this, as long as we want demodulation matrixes as f(lambda)
         dmodm=muellmat*0.
         dmodm2=dmodm*0.
         dmodm_teo=dmodm*0.
         dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
         for ll=0,sizmuel(3)-1 do begin
               modmat[*,*,ll]=modmator##muellmat[*,*,ll]
               modmat2[*,*,ll]=modmat[*,*,ll]*dualmat
               dmodm[*,*,ll]=invert(modmat[*,*,ll])
               dmodm2[*,*,ll]=invert(modmat2[*,*,ll])
         endfor
;prepare other demodulation matrixes
;Average of the whole FOV
         dmodmav=total(dmodm,3)/float(sizmuel(3))
         dmodmav2=total(dmodm2,3)/float(sizmuel(3))
         end

      4: begin
         modmator=modmat
         modmat=muellmat*0.
         modmat2=modmat*0.
;this, as long as we want 2D demodulation matrixes
         dmodm=muellmat*0.
         dmodm2=dmodm*0.
         dmodm_teo=dmodm*0.
         dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
         for xx=0,sizmuel(3)-1 do begin
            for yy=0,sizmuel(4)-1 do begin
               modmat[*,*,xx,yy]=modmator##muellmat[*,*,xx,yy]
               modmat2[*,*,xx,yy]=modmat[*,*,xx,yy]*dualmat
               dmodm[*,*,xx,yy]=invert(modmat[*,*,xx,yy])
               dmodm2[*,*,xx,yy]=invert(modmat2[*,*,xx,yy])
            endfor 
         endfor
;prepare other demodulation matrixes
;Average of the whole FOV
         dmodmav=total(total(dmodm,3),3)/float(sizmuel(3))/float(sizmuel(4))
         dmodmav2=total(total(dmodm2,3),3)/float(sizmuel(3))/float(sizmuel(4))
;Average to a given number of demod. mats. or in a given size of pixs.
;First, fix the number of demodulation matrixes so everything is
;covered and enters the FOV.
;For now, consider square image/mueller matrix
;     submatsiz=info.sizmuel(3)/float(info.demmats)
;     submat=info.sizmuel(3)/round(submatsiz)
;     dmodmsub=fltarr(4,4,submat,submat)
         end

      5: begin
         modmator=modmat
         modmat=muellmat*0.
         modmat2=modmat*0.
;this, as long as we want 2D demodulation matrixes and lambda dependence
         dmodm=muellmat*0.
         dmodm2=dmodm*0.
         dmodm_teo=dmodm*0.
         dmodm2_teo=dmodm*0.
;for now, without the 'teo' consideration, just 'normal' modulation
         for ll=0,sizmuel(3)-1 do begin
            for xx=0,sizmuel(4)-1 do begin
               for yy=0,sizmuel(5)-1 do begin
                  modmat[*,*,ll,xx,yy]=modmator##muellmat[*,*,ll,xx,yy]
                  modmat2[*,*,ll,xx,yy]=modmat[*,*,ll,xx,yy]*dualmat
                  dmodm[*,*,ll,xx,yy]=invert(modmat[*,*,ll,xx,yy])
                  dmodm2[*,*,ll,xx,yy]=invert(modmat2[*,*,ll,xx,yy])
               endfor 
            endfor
         endfor 
;prepare other demodulation matrixes
;Average of the whole FOV and lambdas
         dmodmav=total(total(total(dmodm,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         dmodmav2=total(total(total(dmodm2,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         end
   endcase
   
   
; ideal demodulation matrix (this just 4x4)
;   a = 0.25
;   b = 0.25*sqrt(3.)
;   dmodm[*,0] = a * [1,1,1,1]
;   dmodm[*,1] = b * [1,1,-1,-1]
;   dmodm[*,2] = b * [1,-1,-1,1]
;   dmodm[*,3] = b * [1,-1,1,-1]



endif ;vector. ideal


; ------------------------------------------------------------------------------
; REAL modulation/demodulation, single beam & dual beam
; ------------------------------------------------------------------------------

if info.modtype eq 'vector. real' then begin
;generate the modulation matrixes with the input settings
;   sophism_polmeas_modscheme_retardance,muellmat_teo,muellmat,tt1_teo,tt2_teo,tt1,tt2
;with the single-beam scheme, the polarizer is set just after the
;PMP. Thus, no mueller addition can be induced between both, as
;prepared here originally. Then, muellmat_teo,muellmat must be
;identity not to interfere
;   sophism_polmeas_modscheme_retardance,muellmat_teo,muellmat,modmat_teo,modmat2_teo,modmat,modmat2

;creating FOV-dependence of the LCVR?
;could take the input retardance and run a routine to create a 2D
;array with some variation
;retardance_2d,info.retdeg1,info.retdeg2,pha1,pha2
;this below, just provisional
;I'm considering phases for both LCVRs have dependencies on the
;same dimensions (both on lambda, or both on FOV, or everything)
;for j=0,119 do begin ret[*,j]=sin((findgen(120)+j*(-20)/45.)/(!pi*2.5)+!pi*0)*1/100.
   pha1=info.retdeg1
   pha2=info.retdeg2

   sizpha=size(pha1)
;I'm supposing that, in the LCVR, only the retardance may
;change with FOV. The angle?

;to introduce the FOV-dependence, create a three possibility 'if'
;tree, based on the dimension number of muellmat and retardance?
;the retardance would have 1 dim for the 4 modulation states and
;another 2 for the x/y variation: 3 in total
;   if (sizmuel(0) ne 4 AND sizpha(0) ne 3) then begin
;Mueller matrix could be 3 dimensions or 5, when lambda-dependent
   if (sizmuel(0) eq 2 AND sizpha(0) eq 1) then begin
;no FOV-dep or lambda-dep
;      sophism_polmeas_modscheme_retardance,identity(4),identity(4),modmat_teo,modmat
      sophism_polmeas_modscheme_retardance,identity(4),identity(4),info.reterrors,lcvrerror,pha1,pha2,info.retang1,info.retang2,modmat_teo,modmat
      modmat=modmat##muellmat
      modmat2=modmat*dualmat
      dmodm=invert(modmat)
;      dmodm2=invert(modmat2)
      dmodm2=dmodm*transpose(dualmat)
      dmodm_teo=invert(modmat_teo##muellmat_teo)
;      dmodm2_teo=invert((modmat_teo##muellmat_teo)*dualmat)
      dmodm2_teo=dmodm_teo*transpose(dualmat)
   endif else begin
;FOV-dep
      if (sizpha(0) eq 1) then begin
;Mueller FOV-dep or lambda-dep, modulation matrix not
;         sophism_polmeas_modscheme_retardance,identity(4),identity(4),modmat_teo,modmator
         sophism_polmeas_modscheme_retardance,identity(4),identity(4),info.reterrors,lcvrerror,pha1,pha2,info.retang1,info.retang2,modmat_teo,modmator
         modmat=muellmat*0.
         modmat2=modmat
         dmodm=modmat
         dmodm2=dmodm
         dmodm_teo=dmodm
         dmodm2_teo=dmodm
;maybe create a fixed 5-dimensions matrix, even with 1-element
;dimensions and then works without splitting in cases? 
;         for ll=0,sizmuel(3)-1 do begin
;            for xx=0,sizmuel(4)-1 do begin
;               for yy=0,sizmuel(5)-1 do begin
;                  modmat[*,*,ll,xx,yy]=modmator##muellmat[*,*,ll,xx,yy]
;                  modmat2[*,*,ll,xx,yy]=modmat[*,*,ll,xx,yy]*dualmat
;                  dmodm[*,*,ll,xx,yy]=invert(modmat[*,*,ll,xx,yy])
;                  dmodm2[*,*,ll,xx,yy]=dmodm[*,*,ll,x,y]*dualmat
;                  dmodm_teo[*,*,ll,xx,yy]=invert(modmator##muellmat_teo[*,*,ll,xx,yy])
;                  dmodm2_teo[*,*,ll,xx,yy]=dmodm_teo[*,*,ll,xx,yy]*dualmat
;               endfor
;            endfor
;         endfor

;split according to Mueller dep: lambda, FOV, or both
         case sizmuel(0) of
         
            3: begin
               for ll=0,sizmuel(3)-1 do begin
                  modmat[*,*,ll]=modmator##muellmat[*,*,ll]
                  modmat2[*,*,ll]=modmat[*,*,ll]*dualmat
                  dmodm[*,*,ll]=invert(modmat[*,*,ll])
                  dmodm2[*,*,ll]=invert(modmat2[*,*,ll])
                  dmodm_teo[*,*,ll]=invert(modmator##muellmat_teo[*,*,ll])
                  dmodm2_teo[*,*,ll]=dmodm_teo[*,*,ll]*transpose(dualmat)
               endfor
               dmodmav=total(dmodm,3)/float(sizmuel(3))
               dmodmav2=total(dmodm2,3)/float(sizmuel(3))
               dmodmav_teo=total(dmodm_teo,3)/float(sizmuel(3))
               dmodmav_teo2=total(dmodm2_teo,3)/float(sizmuel(3))
               end 

            4: begin
               for xx=0,sizmuel(3)-1 do begin
                  for yy=0,sizmuel(4)-1 do begin
                     modmat[*,*,xx,yy]=modmator##muellmat[*,*,xx,yy]
                     modmat2[*,*,xx,yy]=modmat[*,*,xx,yy]*dualmat
                     dmodm[*,*,xx,yy]=invert(modmat[*,*,xx,yy])
                     dmodm2[*,*,xx,yy]=invert(modmat2[*,*,xx,yy])
                     dmodm_teo[*,*,xx,yy]=invert(modmator##muellmat_teo[*,*,xx,yy])
                     dmodm2_teo[*,*,xx,yy]=dmodm_teo[*,*,xx,yy]*transpose(dualmat)
                  endfor
               endfor
               dmodmav=total(total(dmodm,3),3)/float(sizmuel(3))/float(sizmuel(4))
               dmodmav2=total(total(dmodm2,3),3)/float(sizmuel(3))/float(sizmuel(4))
               dmodmav_teo=total(total(dmodm_teo,3),3)/float(sizmuel(3))/float(sizmuel(4))
               dmodmav_teo2=total(total(dmodm2_teo,3),3)/float(sizmuel(3))/float(sizmuel(4))
               end

            5: begin
               for ll=0,sizmuel(3)-1 do begin
                  for xx=0,sizmuel(4)-1 do begin
                     for yy=0,sizmuel(5)-1 do begin
                        modmat[*,*,ll,xx,yy]=modmator##muellmat[*,*,ll,xx,yy]
                        modmat2[*,*,ll,xx,yy]=modmat[*,*,ll,xx,yy]*dualmat
                        dmodm[*,*,ll,xx,yy]=invert(modmat[*,*,ll,xx,yy])
                        dmodm2[*,*,ll,xx,yy]=invert(modmat2[*,*,ll,xx,yy])
                        dmodm_teo[*,*,ll,xx,yy]=invert(modmator##muellmat_teo[*,*,ll,xx,yy])
                        dmodm2_teo[*,*,ll,xx,yy]=dmodm_teo[*,*,ll,xx,yy]*transpose(dualmat)
                     endfor
                  endfor
               endfor
               dmodmav=total(total(total(dmodm,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
               dmodmav2=total(total(total(dmodm2,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
               dmodmav_teo=total(total(total(dmodm_teo,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
               dmodmav_teo2=total(total(total(dmodm2_teo,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
               end 

         endcase 

      endif else begin
;FOV-dep or lambda-dep of everything
;Select which case Mueller and phase are in and fit accordingly
         case (5-sizmuel(0)) of
            1: begin
               muellmat=reform(muellmat,sizmuel(1),sizmuel(2),1,sizmuel(3),sizmuel(4))
               muellmat_teo=reform(muellmat_teo,sizmuel(1),sizmuel(2),1,sizmuel(3),sizmuel(4))
               end 
            2: begin
               muellmat=reform(muellmat,sizmuel(1),sizmuel(2),sizmuel(3),1,1)
               muellmat_teo=reform(muellmat_teo,sizmuel(1),sizmuel(2),sizmuel(3),1,1)
               end 
            3: begin
               muellmat=reform(muellmat,sizmuel(1),sizmuel(2),1,1,1)
               muellmat_teo=reform(muellmat_teo,sizmuel(1),sizmuel(2),1,1,1)
               end
         endcase
         sizmuel=size(muellmat)

         case (4-sizpha(0)) of 
            1: begin
               pha1=reform(pha1,sizpha(1),1,sizpha(2),sizpha(3))
               pha2=reform(pha2,sizpha(1),1,sizpha(2),sizpha(3))
               end
            2: begin
               pha1=reform(pha1,sizpha(1),sizpha(2),1,1)
               pha2=reform(pha2,sizpha(1),sizpha(2),1,1)
               end
         endcase 
         sizpha=size(pha1)
;now, rebin muellmat and phases so that both have the same size in all
;dimensions
;take only the appropiate dimensions of lambda, x and y
         modsiz=sizpha[2:4]>sizmuel[3:5]
         muellmat=rebin(muellmat,[sizmuel(1),sizmuel(2),modsiz])
         muellmat_teo=rebin(muellmat_teo,[sizmuel(1),sizmuel(2),modsiz])
         pha1=rebin(pha1,[sizpha(1),modsiz])
         pha2=rebin(pha2,[sizpha(1),modsiz])
         sizmuel=size(muellmat)
         sizpha=size(pha1)
;stop
         modmat=muellmat*0.
         modmat2=modmat
         dmodm=modmat
         dmodm2=dmodm
         dmodm_teo=dmodm
         dmodm2_teo=dmodm
         prg0=0 ;this just for info.talk display
         for ll=0,sizmuel(3)-1 do begin
            for xx=0,sizmuel(4)-1 do begin
               for yy=0,sizmuel(5)-1 do begin
;This way of preparing would create random errors each time, i.e. for
;every pixel. Too much?
                  sophism_polmeas_modscheme_retardance,identity(4),identity(4),info.reterrors,lcvrerror,pha1[*,ll,xx,yy],pha2[*,ll,xx,yy],info.retang1,info.retang2,modmat_teo,modmator
;stop
                  modmat[*,*,ll,xx,yy]=modmator##muellmat[*,*,ll,xx,yy]
                  modmat2[*,*,ll,xx,yy]=modmat[*,*,ll,xx,yy]*dualmat
                  dmodm[*,*,ll,xx,yy]=invert(modmat[*,*,ll,xx,yy])
                  dmodm2[*,*,ll,xx,yy]=invert(modmat2[*,*,ll,xx,yy])
                  dmodm_teo[*,*,ll,xx,yy]=invert(modmator##muellmat_teo[*,*,ll,xx,yy])
                  dmodm2_teo[*,*,ll,xx,yy]=dmodm_teo[*,*,ll,xx,yy]*transpose(dualmat)
                  prg=round(100.*float(ll+1)*float(xx+1)*float(yy+1)/(sizmuel(3)*sizmuel(4)*sizmuel(5)))
                  if (info.talk eq 1 and prg gt prg0) then begin
                     print,format='(37(%"\b"),"Generating Modulation matrices ",I3," %",$)',prg
                     prg0=prg
                  endif 
               endfor 
            endfor  
         endfor 
         dmodmav=total(total(total(dmodm,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         dmodmav2=total(total(total(dmodm2,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         dmodmav_teo=total(total(total(dmodm_teo,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
         dmodmav_teo2=total(total(total(dmodm2_teo,3),3),3)/float(sizmuel(3))/float(sizmuel(4))/float(sizmuel(5))
      endelse
   endelse 

;   sophism_polmeas_modscheme_retardance,identity(4),identity(4),modmat_teo,modmat2_teo,modmat,modmat2

;obtain the demodulation matrixes
;;   dmodm=transpose(invert(transpose(tt1),stat))
;;   dmodm2=transpose(invert(transpose(tt2),stat))
;   dmodm=invert(modmat1,stat)
;   dmodm2=invert(modmat2,stat)
; this below, it's done always although only useful for the case of random errors in LCVRs
;;   dmodm_teo=transpose(invert(transpose(tt1_teo),stat))
;;   dmodm_teo2=transpose(invert(transpose(tt2_teo),stat))
;   dmodm_teo=invert(modmat1_teo,stat)
;   dmodm_teo2=invert(modmat2_teo,stat)

;scheme for the modulation function. Repeat the appropriate modulation
;from the modulation matrix for all the samples of a modulation
;state. Include zeros for the change of LCVRs state.
;   modstate=fltarr(nperiod,4)
;   modmat=transpose(tt1)



   endif ;vector. real
stop

; ------------------------------------------------------------------------------
; Finalize
; ------------------------------------------------------------------------------

;restore,'../../SOPHISM2/data/flat_polari/polscene_flaterr_scheme.sav'
print,'Mueller matrixes (theoretical and real - one pix only in case of FOV-dep):'
print,' '
print,muellmat_teo[*,*,0,0,0]
print,' '
print,muellmat[*,*,0,0,0]
print,' '

;split the save file in two, one for modulation one for demodulation,
;due to the possible large size of matrixes with FOV-dependence
save,muellmat_teo,muellmat,modmat,modmat2,filename=info.saves+info.files(progind)+'_scheme.sav'
save,dmodm,dmodm2,dmodm_teo,dmodm2_teo,dmodmav,dmodmav2,dmodmav_teo,dmodmav_teo2,filename=info.saves+info.files(progind)+'_scheme_demod.sav'


end

