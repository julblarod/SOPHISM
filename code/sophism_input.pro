pro sophism_input

; ==============================================================================
;
; DESCRIPTION
;   Prepares a set of synthetic Stokes images from MHD simulations.
;   Writes the filenames in appropriate given format.
;   Selects a subset of 1/2 or 1/4 of the wavelength dimension from
;   each data if selected.
;   Generates a time series if needed, by replication (when only one
;   input data is used) or from a series of given input data (with or
;   without temporal interpolation in the series). 
;   Replicates data and resamples to detector pixel if selected.
;
; CALLED FROM 
;   sophism
;
; SUBROUTINES
;   sophism_input_redim
;
; MAIN INPUTS
;   dataf              Filenames of original input data
;   ntim               Number of time samples needed for the
;                         simulation run; calculated in base to other
;                         settings (exposure time, etalon positions,...)
;   subset             Selection of a wavelength subset from the
;                         original data. Can be chosen from Complete
;                         set, a Subset for half the original
;                         wavelength points or Minimal for a
;                         fourth. All centered as the original one.
;   numfil             Number of files selected as input data
;   dataser            Option to replicate a single input data to
;                         create artificial time series
;   tsampor/tsamp      Original and new time samplings [seconds]. If
;                         different, interpolation is performed
;
; OUTPUT
;   Final input data [fits format] prepared for simulation run
;
; LIMITATION (solved?)
;   Due to reading/writing scheme based on openw/writeu (to save
;   memory), only so many input files can be managed.
;   May divise an alternative solution with read/writefits using nslice/append.
;
; VERSION HISTORY 
;   Alex Feller   v0.1    2010-09-21
;   J. Blanco. 2010. v1_0. Reorganization of whole process. Change to openw
;      routines to have all data opened and interpolate them.
;   J. Blanco. Sept 2012. v1_1. Split subroutine sophism_input_redim to perform
;      spatial replication and resize to detector pixel. Few other
;      changes for further generalization of read arrays.
;   J. Blanco. Oct 2012. v2_0. Reformed reading/writing for time
;      interpolation option and numfil gt 1. Should overcome idl-units
;      limitation problem.
;   J. Blanco. Dec 2012. v2_1. Added printing and displaying options
;      for talk case
;   J. Blanco. Jan 2013. v3_0. Reordering  of resampling, replication,
;      time series-interpolation
;   J. Blanco. Apr 2013. v3_1. Modified verbose lines for better display
;
; ==============================================================================


restore,'../settings/settings_sophism.sav'
print,'sophism_input'

dataf=info.dataf
;number of time samples
ntim=info.ntim 

;         rotrate=0.1            ;degrees/sec
;         rotdata=rotrate*info.tsamp
;         rotdeg=findgen(ntim)*rotdata

;first, version for only one input data, i.e. no time evolution
if (info.numfil eq 1) then begin
   datacub=readfits(info.dataf[0],head,/sil)

;extract subsets of wavelengths (due to huge original dimension of 992 points)
   if (info.subset eq 'Subset') then datacub=datacub[info.wavepointor/4:info.wavepointor-1-(info.wavepointor/4),*,*,*] 
   if (info.subset eq 'Minimal') then datacub=datacub[3*info.wavepointor/8:info.wavepointor-1-(3*info.wavepointor/8),*,*,*] 

   if (info.talk eq 1) then print,'Replicating and/or resampling data'
   sophism_input_redim,datacub

; re-adjust the variables of size and sampling into the settings file
   sizy=size(datacub)
   info.sz=sizy[3]
   if (info.talk eq 1) then print,'Final data size: ',sizy[1:4]
   resamp=info.sssamp/info.fpaplate
   if (info.fparesamp eq 1) then info.sssamp=info.fpaplate
   save,filenam='../settings/settings_sophism.sav',info
   save,filena=info.fsett,info

;write data. First, case of no time series
   if (info.dataser eq 0) then begin
      sxaddpar,head,"HISTORY",'Input Module'
      sxaddpar,head,"HISTORY",'No time series'
      sxaddpar,head,"MAGNIF",info.sscene_mag,'Replication factor'
      if (info.fparesamp eq 1) then sxaddpar,head,"RESAMP",resamp,'Resampling factor'
      sxaddpar,head,"SPATSAMP",info.sssamp,'Spatial sampling'
      if (info.compress eq 1) then writefits,info.saves+info.files(0)+'_0.fits',datacub,head,/compress else writefits,info.saves+info.files(0)+'_0.fits',datacub,head
      if (info.talk eq 1) then print,'No time series generated'
   endif else begin

;create an artificial time series just by repeating the only input data
      if (info.talk eq 1) then print,'Time series with no solar evolution'

      sxaddpar,head,"HISTORY",'Input Module'
      sxaddpar,head,"HISTORY",'Time series, no solar evolution'
      sxaddpar,head,"MAGNIF",info.sscene_mag,'Replication factor'
      if (info.fparesamp eq 1) then sxaddpar,head,"RESAMP",resamp,'Resampling factor'
      sxaddpar,head,"SPATSAMP",info.sssamp,'Spatial sampling'
      for nf=0,info.ntim-1 do begin

;         sophism_input_rot,datacub,rotdeg[nf],cut=info.rotcut

         if (info.compress eq 1) then writefits,info.saves+info.files(0)+'_'+strtrim(nf,2)+'.fits',datacub,head,/compress else writefits,info.saves+info.files(0)+'_'+strtrim(nf,2)+'.fits',datacub,head

         if (info.talk eq 1) then print,format='(25(%"\b"),"Generating series: ",I3," %",$)',float(nf)/(info.ntim-1)*100.
      endfor

   endelse
   print,''
endif

;case of a list of input files
if (info.numfil gt 1) then begin
   if (info.talk eq 1) then print,'Time series, solar evolution'
   for nf=0,info.numfil-1 do begin
      datacub=readfits(info.dataf[nf],head,/sil)

;extract subsets of wavelengths (due to huge original dimension of 992 points)
      if (info.subset eq 'Subset') then datacub=datacub[info.wavepointor/4:info.wavepointor-1-(info.wavepointor/4),*,*,*] 
      if (info.subset eq 'Minimal') then datacub=datacub[3*info.wavepointor/8:info.wavepointor-1-(3*info.wavepointor/8),*,*,*] 

      if (info.talk eq 1) then print,format='(42(%"\b"),"Replicating and/or resampling data: ",I3," %",$)',float(nf)/(info.numfil-1)*100.
      sophism_input_redim,datacub

;      if (info.tsamp eq info.tsampor AND...) then sophism_input_rot,datacub,rotdeg[nf],cut=info.rotcut


;write fits data. This will be output for no time interpolation case
      sxaddpar,head,"HISTORY",'Input Module'
      sxaddpar,head,"HISTORY",'Time series solar evolution'
      sxaddpar,head,"TSAMP",info.tsamp,'Time sampling (s)'
      sxaddpar,head,"MAGNIF",info.sscene_mag,'Replication factor'
      if (info.fparesamp eq 1) then begin
         sxaddpar,head,"RESAMP",info.sssamp/info.fpaplate,'Resampling factor'
         sxaddpar,head,"SPATSAMP",info.fpaplate,'Spatial sampling'
      endif else sxaddpar,head,"SPATSAMP",info.sssamp,'Spatial sampling'
      sxaddpar,head,"TSAMP",format='(F8.3)',info.tsamp,'Time sampling (s)'
      if (info.compress eq 1) then writefits,info.saves+info.files(0)+'_'+strtrim(nf,2)+'.fits',datacub,head,/compress else writefits,info.saves+info.files(0)+'_'+strtrim(nf,2)+'.fits',datacub,head
   endfor ;nf

   sizy=size(datacub)
   info.sz=sizy[3]
   if (info.talk eq 1) then print,'Final data size: ',sizy[1:4]
   if (info.fparesamp eq 1) then info.sssamp=info.fpaplate
   save,filenam='../settings/settings_sophism.sav',info
   save,filena=info.fsett,info

;case of time interpolation
   if (info.tsamp ne info.tsampor) then begin
      if (info.talk eq 1) then print,'Time interpolation'
;prepare array of interpolation positions.
      ninterp=indgen(ntim)*info.tsamp/info.tsampor

;work with list of files row by row: read one for each file - 'store'
;them in a new array - interpolate that row-temporal dimension to new
;temporal sampling - save each new row in new files
      for rr=0,sizy[3]-1 do begin
         datarow=fltarr(sizy[1],sizy[2],sizy[3],info.numfil)
         for nf=0,info.numfil-1 do begin
            strnam=info.saves+info.files(0)+'_'+strtrim(nf,2)
            datamed=readfits(strnam+'.fits*',head,nsli=rr,/comp,/silent)
            datarow[*,*,*,nf]=datamed
         endfor ;nf
         datarow=interpolate(datarow,ninterp)
;writing of new interpolated files done with openw because of its
;append option (writefits should work but...)
         if (rr eq 0) then begin
            for ni=0,n_elements(ninterp)-1 do begin
               strnamint=info.saves+info.files(0)+strtrim(ni,2)
               openw,lun,strnamint,/get_lun
               writeu,lun,datarow[*,*,*,ni]
               close,/all
            endfor
         endif else begin

            for ni=0,n_elements(ninterp)-1 do begin
               strnamint=info.saves+info.files(0)+strtrim(ni,2)
               openw,lun,strnamint,/get_lun,/append
               writeu,lun,datarow[*,*,*,ni]
               close,/all
            endfor 
         endelse ;rr eq 0

         if (info.talk eq 1) then print,format='(20(%"\b"),"Interpolating: ",I3," %",$)',float(rr)/(info.sz-1)*100.
      endfor ;rr
      print,''
      close,/all

;convert final data from idl-write to writefits with header
      sxaddpar,head,"TSAMP",format='(F8.3)',info.tsamp,'Time sampling (s)'
;      sxaddpar,head,"MAGNIF",format='(F6.2)',info.sscene_mag,'Replication factor'
;      if (info.fparesamp eq 1) then sxaddpar,head,"RESAMP",info.sssamp/info.fpaplate,'Resampling factor'

      for ni=0,n_elements(ninterp)-1 do begin
         datafin=fltarr(info.wavepoint,4,info.sz,info.sz)
         strnamint=info.saves+info.files(0)+strtrim(ni,2)
         openr,5,strnamint
         readu,5,datafin

;         sophism_input_rot,datafin,rotdeg[ni],cut=info.rotcut

         if (info.compress eq 1) then writefits,info.saves+info.files(0)+'_'+strtrim(ni,2)+'.fits',datafin,head,/compress else writefits,info.saves+info.files(0)+'_'+strtrim(ni,2)+'.fits',datafin,head
         close,5
;remove idl-write files
         spawn,'rm '+info.saves+info.files(0)+strtrim(ni,2)
         if (info.talk eq 1) then print,format='(25(%"\b"),"Converting to fits: ",I3," %",$)',float(ni)/(n_elements(ninterp)-1)*100.
      endfor  ;ni
      print,''
      close,/all

   endif ;tsamp ne tsampor


endif ;numfil and tsamp

;show some images if selected, at -400 mAA
if (info.talk eq 1) then begin
   cc=info.wavepoint/2.-round(400./info.wavesamp)
   datafin=readfits(info.saves+info.files(0)+'_'+strtrim(0,2)+'.fits*',head,/comp,/sil)
   sz=size(datafin)
   if (cc lt 0) then cc=0
   cont=reform(datafin(cc,0,*,*))
   si = reform(datafin(cc,0,*,*))/mean(cont)
   if (sz(2) lt 4) then begin
      window,0,title='Input data at -'+strtrim(fix((round(400./info.wavesamp) < info.wavepoint/2.)*info.wavesamp),2)+' mAA'
      tvframe,si,/bar,/aspect,btitle='I!DC!N'
   endif else begin
      sv = reform(datafin(cc,1,*,*))/cont
      mm=max(abs(sv))
      sv[0,0]=mm
      sv[0,1]=-mm
      sq = reform(datafin(cc,2,*,*))/cont
      mm=max(abs(sq))
      sq[0,0]=mm
      sq[0,1]=-mm
      su = reform(datafin(cc,3,*,*))/cont
      mm=max(abs(su))
      su[0,0]=mm
      su[0,1]=-mm
      !p.multi = [0,2,2]
      window,0,title='Input data at -'+strtrim(fix((round(400./info.wavesamp) < info.wavepoint/2.)*info.wavesamp),2)+' mAA'
      tvframe,si,/bar,/aspect,btitle='I!DC!N'
      tvframe,sq,/bar,/aspect,btitle='Q/I!DC!N'
      tvframe,su,/bar,/aspect,btitle='U/I!DC!N'
      tvframe,sv,/bar,/aspect,btitle='V/I!DC!N'
      !p.multi=0
   endelse 

endif


end

