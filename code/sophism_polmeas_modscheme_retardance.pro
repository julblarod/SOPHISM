
pro sophism_polmeas_modscheme_retardance,mmm_teo,mmm,reterr,lcvrerr,retar1,retar2,ang1,ang2,tt1_teo,tt1

; ==============================================================================
;
; DESCRIPTION
;    Defines modulation scheme, consisting of modulation and demodulation
;    functions and the demodulation matrix, from the information for
;    LCVRs and the Mueller matrix of the system.
;    Calculates the same including random errors, in case they are enabled.
;    Performs also the orthogonal states for dual-beam.
;
; CALLED FROM 
;    sophism_polmeas_modscheme
;
; SUBROUTINES
;    mirror
;    rota
;    device
;
; MAIN INPUTS
;    retang1/2      Orientations of axis of LCVR1 and LCVR2 (deg)
;    retdeg1/2      Retardances for all modulation states of LCVR1 and
;                      LCVR2 (deg)
;    mmm_teo        Mueller matrix of the system without random errors
;    mmm            Mueller matrix of the system
;    reterrors      Array enabling random errors in the parameters of
;                      LCVRs (amplitude ratio, orientation, retardances)
;
; OUTPUT
;    Modulation matrixes for both single and dual-beam in the correct
;    case (including random errors) and theoretical case (withour them)
;
; VERSION HISTORY 
;    V. Martinez Pillet
;    J. Blanco. 2011. v1_0. Fit into sophism
;    J. Blanco. Feb 2012. v2_0. Modified random errors and added dual-beam
;       option
;    J. Blanco. Dec 2012. v2_1. Added varibles from input with values
;       for LCVR errors
;    J. Blanco. Jan 2013. v2_2. Added theoretical (error-free) Mueller
;       matrix. Save error values for LCVRs
;    J. Blanco. Jan 2014. v2_3. Removed dual-beam part, because it is
;       done in modscheme.
;    J. Blanco. Sep 2014. v2_4. Take retardandces, errors and angles
;       as input from command line, better for new cases of FOV and
;       lambda dependence.
;    J. Blanco. Feb 2016. v2_5. Comment the random sign generation for
;       errors, now fixed to positive. Generate different random
;       errors for each modulation state.
;    J. Blanco. Nov 2016. v2_55. Corrected the on-screen info about
;       error signs to appear only when errors enabled
;
; ==============================================================================

restore,'../settings/settings_sophism.sav'

progind=2

plmas=transpose(device(amr=0.,pha=0.,ang=0.))

volt1=[2.54,2.54,3.12,3.12]
volt2=[2.46,9,4.3,2.9]

; temperature dependence.
; ojo!!! esto por ahora pero...!
cht=0.
ret1=(-1.16+0.305*volt1-0.02*volt1^2)*cht
ret2=(-1.16+0.305*volt2-0.02*volt2^2)*cht
ret2[1]=0. ; this is by hand; because for V>8 there is no change with T
; this is a cl of delta ret/delta T (raquel's paper)
; times cht deg celsius change

;random errors in amplitude ratio, phases and orientations of the LCVRs.
;up to now, only in limited fixed ranges. Change in future versions.
ampdif=abs(info.reterramp2-info.reterramp1)
ampoff=(info.reterramp1<info.reterramp2)/float(ampdif)
if (max(info.reterrors) gt 0) then begin
   print,'Error signs fixed to positive, not random'
   print,'Errors generated randomly for each of 4 mod. states'
endif
signpha=1;randomu(seed,1)
if (signpha lt 0.5) then signpha=-1 else signpha=1
phadif=abs(info.reterrpha2-info.reterrpha1)
phaoff=(info.reterrpha1<info.reterrpha2)/float(phadif)
signang=1;randomu(seed,1)
if (signang lt 0.5) then signang=-1 else signang=1
angdif=abs(info.reterrang2-info.reterrang1)
angoff=(info.reterrang1<info.reterrang2)/float(angdif)
lcvr1amperr=(randomu(seed,4)+ampoff)*ampdif*info.reterrors(0)
lcvr1phaerr=(randomu(seed,4)+phaoff)*phadif*signpha*info.reterrors(1)
lcvr1angerr=(randomu(seed,4)+angoff)*angdif*signang*info.reterrors(2)

lcvr2amperr=(randomu(seed,4)+ampoff)*ampdif*info.reterrors(0)
signpha=1;randomu(seed,1)
if (signpha lt 0.5) then signpha=-1 else signpha=1
signang=1;randomu(seed,1)
if (signang lt 0.5) then signang=-1 else signang=1
lcvr2phaerr=(randomu(seed,4)+phaoff)*phadif*signpha*info.reterrors(1)
lcvr2angerr=(randomu(seed,4)+angoff)*angdif*signang*info.reterrors(2)

;parameters defining LCVRs behaviour: orientation and phases
ang1=info.retang1
ang2=info.retang2
rett1=retar1;info.retdeg1
rett2=retar2;info.retdeg2
;stop
		; mod. state 1 
		lcvr1=transpose(device(amr=1.,pha=rett1(0),ang=ang1)) 
		lcvr2=transpose(device(amr=1.,pha=rett2(0),ang=ang2))

		mm1=transpose(plmas#mmm_teo#lcvr2#lcvr1)

		tt1=fltarr(4,4)
		tt1(0:3,0)=mm1(0:3)

		; mod. state 2 
		lcvr1=transpose(device(amr=1.,pha=rett1(1),ang=ang1)) 
		lcvr2=transpose(device(amr=1.,pha=rett2(1),ang=ang2))

		mm2=transpose(plmas#mmm_teo#lcvr2#lcvr1)

		tt1(0:3,1)=mm2(0:3)

		; mod. state 3		
		lcvr1=transpose(device(amr=1.,pha=rett1(2),ang=ang1)) 
		lcvr2=transpose(device(amr=1.,pha=rett2(2),ang=ang2))

		mm3=transpose(plmas#mmm_teo#lcvr2#lcvr1)

		tt1(0:3,2)=mm3(0:3)

		; mod. state 4 
		lcvr1=transpose(device(amr=1.,pha=rett1(3),ang=ang1)) 
		lcvr2=transpose(device(amr=1.,pha=rett2(3),ang=ang2))

		mm4=transpose(plmas#mmm_teo#lcvr2#lcvr1)
					
		tt1(0:3,3)=mm4(0:3)

                tt1_teo=tt1

;this probably for nothing. Coming from old times
reterr=0.
reterr1=ret1[0]*reterr
reterr2=ret1[1]*reterr
reterr3=ret2[0]*reterr
reterr4=ret2[1]*reterr
reterr5=ret2[2]*reterr
reterr6=ret2[3]*reterr

rett1=[retar1(0)+reterr1,retar1(1)+reterr1,retar1(2)+reterr2,retar1(3)+reterr2]
rett2=[retar2(0)+reterr3,retar2(1)+reterr4,retar2(2)+reterr5,retar2(3)+reterr6]

		; mod. state 1 
		lcvr1=transpose(device(amr=1.-lcvr1amperr(0),pha=rett1(0)+lcvr1phaerr(0),ang=ang1+lcvr1angerr(0)))
		lcvr2=transpose(device(amr=1.-lcvr2amperr(0),pha=rett2(0)+lcvr2phaerr(0),ang=ang2+lcvr2angerr(0)))

		mm1=transpose(plmas#mmm#lcvr2#lcvr1)

		tt1=fltarr(4,4)
		tt1(0:3,0)=mm1(0:3)

		; mod. state 2 
		lcvr1=transpose(device(amr=1.-lcvr1amperr(1),pha=rett1(1)+lcvr1phaerr(1),ang=ang1+lcvr1angerr(1)))
		lcvr2=transpose(device(amr=1.-lcvr2amperr(1),pha=rett2(1)+lcvr2phaerr(1),ang=ang2+lcvr2angerr(1)))

		mm2=transpose(plmas#mmm#lcvr2#lcvr1)

		tt1(0:3,1)=mm2(0:3)

		; estado 3	
		lcvr1=transpose(device(amr=1.-lcvr1amperr(2),pha=rett1(2)+lcvr1phaerr(2),ang=ang1+lcvr1angerr(2)))
		lcvr2=transpose(device(amr=1.-lcvr2amperr(2),pha=rett2(2)+lcvr2phaerr(2),ang=ang2+lcvr2angerr(2)))

		mm3=transpose(plmas#mmm#lcvr2#lcvr1)

		tt1(0:3,2)=mm3(0:3)

		; estado 4 
		lcvr1=transpose(device(amr=1.-lcvr1amperr(3),pha=rett1(3)+lcvr1phaerr(3),ang=ang1+lcvr1angerr(3)))
		lcvr2=transpose(device(amr=1.-lcvr2amperr(3),pha=rett2(3)+lcvr2phaerr(3),ang=ang2+lcvr2angerr(3)))

		mm4=transpose(plmas#mmm#lcvr2#lcvr1)
					
		tt1(0:3,3)=mm4(0:3)
	

;save random errors introduced in LCVRs
if (max(info.reterrors) eq 1) then save,filenam=info.saves+info.files(progind)+'_lcvr_errors.sav',lcvr1amperr,lcvr2amperr,lcvr1phaerr,lcvr2phaerr,lcvr1angerr,lcvr2angerr

                
end
