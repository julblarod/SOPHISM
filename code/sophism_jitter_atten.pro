pro sophism_jitter_atten

; ==============================================================================
;
; DESCRIPTION
;    Generates attenuation curve from an image stabilisation system
;       for jittering effects
; 
; CALLED FROM
;    sophism_jitter
;
; SUBROUTINES
;    sophism_make_freq
; 
; MAIN INPUTS
;    
;
; OUTPUT
;    Attenuation curve
;
; VERSION HISTORY 
;    Alex Feller          2010-12-10
;    J. Blanco. 2011. v1_0. Minor variable changes to fit into SOPHISM
;       structure 
;    J. Blanco. Feb 2013. v1_1. Added parameter to underperform ISS
;       (implemented in sophism_2.5 but not included in header)
;
; ==============================================================================

restore,'../settings/settings_sophism.sav'
progind=1

; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------

; number of time samples
ntim = info.ntim 
;define variable for frequency and curve
attenarr = {nu: 0.0, f: complex(0, 0)}
attenarr = replicate(attenarr, ntim)
attenarr.nu = sophism_make_freq(ntim, 1./info.tobs)

; ------------------------------------------------------------------------------
; Main
; ------------------------------------------------------------------------------

;selection of attenuation type. Now, only KIS curve
case info.issatten of
   'kis' : begin
 ; support points
      x = [2.5,5,7.5,10,30]
      y = [100.,27,13,8,1]*info.issunder/100.

; fit power law
      p = linfit(alog(x), alog(y))
      w = where(abs(attenarr.nu) gt 0)
      attenarr[w].f = exp(p[0]) * abs(attenarr[w].nu)^p[1]
      w = where(abs(attenarr.f) lt 1)
      attenarr[w].f = 1
   end
   else: message, 'Attenuation type ' + atten + ' not implemented!'
 
endcase

; ------------------------------------------------------------------------------
; Finalize
; ------------------------------------------------------------------------------

save, attenarr, filename=info.saves+info.files(1)+'_iss.sav'

end

