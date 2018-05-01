;=================================================================
; DESCRIPTION
;   Compilation of auxiliary routines
;
; CALLED FROM 
;   sophism
;
; SUBROUTINES
;
;
; MAIN INPUTS
;
; OUTPUT
;
; VERSION HISTORY
;   J. Blanco. Dec 2012. v1_0.
;   J. Blanco. Dec 2012. v1_1. Reorder of functions
;   J. Blanco. Dec 2012. v1_2. Added new functions to compile (tvframe,avg,...)
;   J. Blanco. Feb 2013. v1_3. New functions related with check_sum of
;      fits files
;   J. Blanco. Apr 2017. v1_4. Included daycnv.pro
;
;=================================================================

@aux/gettok.pro
@aux/valid_num.pro
@aux/sxaddpar.pro
@aux/sxdelpar.pro
@aux/sxpar.pro
@aux/fxpar.pro
@aux/fxparpos.pro
@aux/fxaddpar.pro
@aux/fxmove.pro
@aux/fxposit.pro
@aux/get_date.pro
@aux/mkhdr.pro
@aux/mrd_hread.pro
@aux/mrd_skip.pro
@aux/check_fits
@aux/headfits.pro
@aux/readfits.pro
@aux/is_ieee_big.pro
@aux/host_to_ieee.pro
@aux/fits_ascii_encode.pro
@aux/n_bytes.pro
@aux/checksum32.pro
@aux/fits_add_checksum.pro
@aux/writefits.pro

@aux/avg
@aux/frac_shift
@aux/fft_shift

@aux/crosscorr_c
@aux/otfunc.pro
@aux/radius_aper2.pro
@aux/theta.pro
@aux/check_rpupil
@aux/fact
@aux/odd
@aux/zernike_mn.pro
@aux/zernike2.pro

@aux/sign.pro
@aux/mirror.pro
@aux/rota.pro
@aux/device

@aux/minmax.pro
@aux/resize.pro
@aux/tvframe.pro
@aux/daycnv.pro
