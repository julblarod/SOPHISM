;function modelo_ret_lookuptable,PtV
pro sophism_polmeas_birefringence

; ==============================================================================
;
; DESCRIPTION
;    Obtain the birefringence Mueller matrix from the EW variation
;
; CALLED FROM 
;    sophism_polmeas_birefringence
;
; SUBROUTINES
;    sophism_polmeas_birefringence_model
;
; MAIN INPUTS
;    scale        Reduction factor scale of the data size to speed up
;                    calculations 
;
;
; OUTPUT
;    Mueller matrix of the birefringence
;
; VERSION HISTORY 
;   D. Orozco   v0.1   
;   J. Blanco. Feb 2016. v1_0. Fit into sophism.
;
; ==============================================================================

restore,'../settings/settings_sophism.sav'
progind=2

ptv=info.birefptv

if (info.birefloadmodel eq 1) then begin
   restore,info.saves+info.birefloadmodelnam ;'HREW_mueller_matrices.sav'
endif else begin
;generate the look up tables used below
;.r HREW_model_lt.pro before runing gem_mod_ret?
   ds_PtV = [0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20]
   dsnum=n_elements(ds_ptv)

   scale = info.birefscale ;32.
   x_size=info.sz ;2048.
   y_size=info.sz ;2048.

;   RMatrixB_T = fltarr(x_size/scale,y_size/scale,4,4,10)
   RMatrixB_T = fltarr(x_size/scale,y_size/scale,4,4,dsnum)
   
   for i=0,dsnum-1 do begin     ;n_elements(ds_PtV)-1 do begin    
;      HREW_model_lt,CM,MMatrix,MMatrix_ideal,RMatrixB,ds_ptv=ds_ptv[i];,scale=scale,ds_PtV=ds_PtV[i];,/doplot ;,label='_test'
      sophism_polmeas_birefringence_model,CM,MMatrix,MMatrix_ideal,RMatrixB,ds_ptv=ds_ptv[i]
      RMatrixB_T(*,*,*,*,i) = RMatrixB
   endfor

   M13 = Reform(RMatrixB_T[*,*,1,3,*])
   M23 = Reform(RMatrixB_T[*,*,2,3,*])
   M31 = Reform(RMatrixB_T[*,*,3,1,*])
   M32 = Reform(RMatrixB_T[*,*,3,2,*])
   M11 = Reform(RMatrixB_T[*,*,1,1,*])
   M22 = Reform(RMatrixB_T[*,*,2,2,*])
   M33 = Reform(RMatrixB_T[*,*,3,3,*])
;+ II
   save,filename=info.saves+info.files(progind)+'_HREW_mueller_matrices.sav',M13,M23,M31,M32,M11,M22,M33,ds_PtV,scale,x_size,y_size

;STOP
endelse

;find closest ds_ptv in sav file and returns the Retardance Matrix of
;the HREW

;Shall be between 0.02 and 0.2

Matrix_R = fltarr(x_size,y_size,4,4)
for i=0,3 do Matrix_R[*,*,i,i] = 1d0

;Scale px
;howmany = n_elements(ds_ptv)
if (PtV lt min(ds_ptv) or ptv gt max(ds_ptv)) then begin ;ds_ptv(howmany-1) then begin
   print,''
   print,'Out of ps_ptv range ',ds_ptv
   print,'Assuming identity matrix'
   print,''
;   return,Matrix_R
;   stop
   goto,elfin
endif

;busco el mas cercano
;below = min(where(ds_ptv le ds_ptv))
;above = max(where(ds_ptv ge ds_ptv))
below = max(where(ds_ptv le ptv))
above = min(where(ds_ptv ge ptv))

if below eq above then begin
   Matrix_R[*,*,1,3] = Rebin(M13[*,*,below],x_size,y_size)
   Matrix_R[*,*,2,3] = Rebin(M23[*,*,below],x_size,y_size)
   Matrix_R[*,*,3,1] = Rebin(M31[*,*,below],x_size,y_size)
   Matrix_R[*,*,3,2] = Rebin(M32[*,*,below],x_size,y_size)
   Matrix_R[*,*,1,1] = Rebin(M11[*,*,below],x_size,y_size)
   Matrix_R[*,*,2,2] = Rebin(M22[*,*,below],x_size,y_size)
   Matrix_R[*,*,3,3] = Rebin(M33[*,*,below],x_size,y_size)
endif else begin
;interpolation (linear)
   for i=0,x_size/scale-1 do begin
      for j=0,x_size/scale-1 do begin
         M13[i,j,below] = interpol( M13[i,j,*], ds_ptv,ptv)
         M23[i,j,below] = interpol( M23[i,j,*], ds_ptv,ptv)
         M31[i,j,below] = interpol( M31[i,j,*], ds_ptv,ptv)
         M32[i,j,below] = interpol( M32[i,j,*], ds_ptv,ptv)
         M11[i,j,below] = interpol( M11[i,j,*], ds_ptv,ptv)
         M22[i,j,below] = interpol( M22[i,j,*], ds_ptv,ptv)
         M33[i,j,below] = interpol( M33[i,j,*], ds_ptv,ptv)
      endfor
   endfor
   Matrix_R[*,*,1,3] = Rebin(M13[*,*,below],x_size,y_size)
   Matrix_R[*,*,2,3] = Rebin(M23[*,*,below],x_size,y_size)
   Matrix_R[*,*,3,1] = Rebin(M31[*,*,below],x_size,y_size)
   Matrix_R[*,*,3,2] = Rebin(M32[*,*,below],x_size,y_size)
   Matrix_R[*,*,1,1] = Rebin(M11[*,*,below],x_size,y_size)
   Matrix_R[*,*,2,2] = Rebin(M22[*,*,below],x_size,y_size)
   Matrix_R[*,*,3,3] = Rebin(M33[*,*,below],x_size,y_size)
   
endelse

;diagonal elements equal 1. Rest of elements are zero

;return,Matrix_R
elfin:
save,filenam=info.saves+info.files(progind)+'_biref_mat.sav',matrix_r

end





