pro sophism_jitter_tpfilt

; ==============================================================================
;
; DESCRIPTION
;   Defines a jittering frequency filter from a selection of possibilities.
;
; CALLED FROM
;    sophism_jitter
;
; SUBROUTINES
;    sophism_make_freq
; 
; MAIN INPUTS
;    filtype      Determines which type of filtering will be applied
;
; OUTPUT
;    Filtering curve for jittering power spectrum
;
; VERSION HISTORY 
;    Alex Feller   v0.1    2010-09-21
;    J. Blanco. 2011. v1_0. Minor variable changes to fit into SOPHISM structure
;    A. Feller. Feb 2013. v1_1. Hinode filter modifications
;
; ==============================================================================

restore,'../settings/settings_sophism.sav'
progind=1
;dummy, to save later the Hinode filter parameters
p=0

; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------

ntim = info.ntim 
tpfilt = {nu: 0.0, pow: 0.0}
tpfilt = replicate(tpfilt, ntim)
tpfilt.nu = sophism_make_freq(ntim, 1./info.tobs)
;stop
; ------------------------------------------------------------------------------
; Cutoff filter
; ------------------------------------------------------------------------------

if info.filtype eq 'cutoff' then begin
   w = where(abs(tpfilt.nu) le info.filcutoff)
   tpfilt[w].pow = 1.0
endif

; ------------------------------------------------------------------------------
; Cutout filter
; ------------------------------------------------------------------------------

if info.filtype eq 'cutout' then begin
   a = abs(tpfilt.nu)
   w = where(a ge info.filcut1 and a le info.filcut2)
   tpfilt[w].pow = 1.0
endif

; ------------------------------------------------------------------------------
; Hinode/SOT filter
;
; Power-law approximation of the power spectrum of the Hinode/SOT pointing
; error. 
; 
; Ref.: Y. Katsukawa et al. 2010, "Pointing Stability of Hinode and Requirements
; for the Next Solar Mission Solar-C", in Intern. Conference on Space Optics,
; Rhodes, Greece, 4-8/10/2010
;
; Normalization to different RMS: see sophism_jitter.pro
; ------------------------------------------------------------------------------

if info.filtype eq 'Hinode' then begin

; Parameters describing the approximation of the Hinode power spectrum
;
; p[2]: threshold frequency in [Hz]
; In the low-frequency region (frequencies < p[2]):
; power density [arcsec^2/Hz] = exp(p[0]) * frequency [Hz] ^ p[1]
; In the high-frequency region (frequencies >= p[2]):
; power density [arcsec^2/Hz] = p[3]

;NOT in the new version. Now no limit in high-frequency:
;       For frequency [Hz] >= filcuth:
;         power density [arcsec^2/Hz] = exp(p[0]) * frequency [Hz] ^ p[1]
;       For frequency [Hz] < filcuth:
;         power density [arcsec^2/Hz] = 0

   p = [-10.1774, -1.97677, 45.0000, 2.10000e-08]

   w = where(abs(tpfilt.nu) ge info.filcuth, complement=wc)
   tpfilt[w].pow = exp(p[0]) * abs(tpfilt[w].nu)^p[1]
;stop
;this below to use Johanns spectra of jittering
;restore,'../data/ESA_jitter_spectrum_Johann.sav'
;spectrum2=spectrum2>0
;tpfilt.pow=interpol(spectrum2,f,abs(tpfilt.nu))
   tpfilt[wc].pow = 0.

;   w = where( (tpfilt.nu ne 0) and (abs(tpfilt.nu) ge p[2]) )
;careful! 'cuidado! aqui cambio porque las freq son pequenhas
;muchas veces para casos de pocos datos'
;   if (max(w) ne -1) then tpfilt[w].pow = p[3]

; Set value at nu=0 to next non-zero nu value
;   w = where(tpfilt.nu eq 0)
;   tpfilt[w].pow = tpfilt[max(w)+1].pow

; Clip frequency range
;   a = abs(tpfilt.nu)
;   w = where( (a lt info.filcut1) or (a gt info.filcut2) )
;   if (max(w) ne -1) then tpfilt[w].pow = 0.0

endif

; ------------------------------------------------------------------------------
; Levels filter
; ------------------------------------------------------------------------------

if info.filtype eq 'levels' then begin
   a = abs(tpfilt.nu)
   nrange = n_elements(info.filrange)/2
   for i=0,nrange-1 do begin
      level = info.fillevel[i]^2 / (info.filrange[2*i+1]-info.filrange[2*i])
      w = where(a ge info.filrange[2*i] and a lt info.filrange[2*i+1])
      tpfilt[w].pow = level
   endfor

; smoothing
   if info.filsmooth gt 0 then begin
      dnu = abs(tpfilt[1].nu-tpfilt[0].nu)
      tpfilt.pow = smooth(tpfilt.pow, round(info.filsmooth/dnu))
   endif

endif

; ------------------------------------------------------------------------------
; Finalize
; ------------------------------------------------------------------------------

save, tpfilt, p, filename=info.saves+info.files(progind)+'_ftpfilt.sav'


end

