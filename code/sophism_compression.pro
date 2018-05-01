pro sophism_compression

; ==============================================================================
;
; DESCRIPTION
;   Simulate the data compression.
;   Right now (June 2016), the compression code is only available in
;   windows executable. Thus, it has to be run through wine command. 
;   In principle, any wine linux distribution should work.
;
; CALLED FROM
;   sophism
;
; SUBROUTINES
;    
;
; MAIN INPUTS
;   comprhw        Simulate hardware compression (1) or software only (0)
;   comprratio     Target compression ratio (1=lossless)
;   comprbpp       Bits per pixel. Up to 16
;   comprbps       Blocks per segment. Default: 256 (software), 32 (hw)
;   comprskip      Skip possible header bytes in input file before
;                  compression (0)
;   comprmsb       Most Significant Byte first for all input and
;                  output files (1). 0: LSB is used (e.g. windows systems)
;
; OUTPUT
;    Compressed and decompressed data
;
; VERSION HISTORY 
;   A. Lagg   May 2014. Original routines
;   J. Blanco May 2014. v1_0. Fit into sophism
;   J. Blanco Jul 2014. v1_1. Adapted for observation sets. Included
;      cases for theoretical (random errors in modulation), average
;      (FOV and/or lambda dependent modulation), and dual (if
;      it's not demodulated data which are going to be used)
;   J. Blanco Nov 2014. v1_2. Corrected the output decompressed files
;      size, depending whether it has been replicated or not
;   J. Blanco May 2015. v1_3. Corrected bug in the progma variable
;      when loading input data directly
;   J. Blanco Jun 2016. v1_4. Modified replication for cases where
;      input data size < 128 pix. Now replicates the necessary
;      factor. Re-arranged histogram display. Remove intermediate files
;
; ==============================================================================


restore,'../settings/settings_sophism.sav'
progind=10

; number of time exposures
;info.modnbuff=1
;nperiod=info.modnbuff*info.nacc
;nmcycl=info.cycrep
;ntim=nperiod*nmcycl*info.etalpos


if (max(info.progma) gt 0 AND min(where(info.routines eq 1)) eq progind) then progma=info.progma else progma=max(where(info.routines(0:progind-1) eq 1))
if (progma eq -1) then progma=0

for obs=0,info.obsset-1 do begin

   fileout=info.saves+info.files(progind)+'_'+strtrim(obs,2)

;only one output now, from the demodulation module
;(should have to be generalized? And the inversion output?)

   data=readfits(info.saves+info.files(progma)+'_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)
   sizd=size(data)

;if (info.talk eq 1) then begin
;   print,'Compression parameters: '
;   print,'Ratio: ',ratio
;   print,'Bits per pixel: ',bpp
;   print,'Blocks per segment: ',bps
;endif

;variable to store the results' compression ratios
;rat=fltarr(info.etalpos*info.modnbuff)
   rat=fltarr(sizd(1)*sizd(2))  ;info.etalpos*info.modnbuff)
;storing the ranges of the images to 'recover' after decompressing
;something similar to the input
   rangetot=fltarr(2,sizd(1)*sizd(2)) ;info.etalpos*info.modnbuff)

;do a loop to compress each lambda/stokes?
;Creo que si. Pero para inversion? Como agruparlo? En un solo fits
;especificando en la cabecera a que corresponde cada dimension? Puede.
   for ll=0,sizd(1)-1 do begin  ;info.etalpos-1 do begin 
      for ss=0,sizd(2)-1 do begin ;info.modnbuff-1 do begin
         img=reform(data[ll,ss,*,*])
;)))))
;crazy tests
;img=rotate(img,3)
;img[0:511,*]=0
;img[*,0:511]=0
;img[2047-512:2047,*]=0
;img[*,2047-512:2047]=0
;img=img[512:1535,512:1535]
;img=img*(img gt 9e3)
;((((
;***********************
;these lines below are needed for too small input data. Less than 128
;pixels seems not possible for compression 
;creates problems when flats/fringes present because of the intesity
;gradients, it is no longer periodic
         siz=size(img,/dim)
         flagsiz=0
         if (siz(0) lt 128 OR siz(1) lt 128) then begin
            if (ll eq 0 AND ss eq 0) then print,'Data smaller than 128 pixels. Replicating...'
            flagsiz=1
            img0=img
            facrep=ceil(128./min(siz))
;            img=fltarr(siz[0]*2,siz[1]*2)
;            img[0:siz[0]-1,0:siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,0:siz[1]-1]=img0
;            img[0:siz[0]-1,siz[1]:2*siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,siz[1]:2*siz[1]-1]=img0
            img=fltarr(siz[0]*facrep,siz[1]*facrep)
            for jjy=0,facrep-1 do begin
               for jjx=0,facrep-1 do begin
                  img[siz[0]*jjx:siz[0]*(jjx+1)-1,siz[1]*jjy:siz[1]*(jjy+1)-1]=img0
               endfor 
            endfor 
;***********************
            siz=size(img,/dim)
         endif 
;      if (info.modnbuff*ll+ss eq 3 OR info.modnbuff*ll+ss eq 15) then img=reform(data[ll+1,ss,*,*]) ;img/10.
;stop
;      img=img*(img gt 6e3)
;stop

;define image range
;depending if there are negative values or not, use 0 for min or
;proper minimum
;      range=[min(img),max(img)]
         if (min(img) ge 0) then range=[0,max(img)] else range=[min(img),max(img)] 
         rangetot[*,sizd(2)*ll+ss]=range

;image must be cut to 2^n dimensions. Cut image to closest 2^n size
;cortar desde la esquina? o alrededor del centro?
         siz=2l^fix(alog(siz)/alog(2))
         img=img[0:siz[0]-1,0:siz[1]-1]
;stop
         siz=size(img,/dim)
         rows=siz[0]
         columns=siz[1]

;define the output variable with the new size
;no, keep the original size and cut the extra replication
;CAREFUL!! no, use the new size. In case of replication, it will go further than
;original but if not, it will be smaller
;put another flagsiz and no problem
;         if (ll eq 0 AND ss eq 0) then decompima=fltarr(sizd(1),sizd(2),rows,columns) ;info.etalpos,info.modnbuff,rows,columns)
         if (ll eq 0 AND ss eq 0) then begin
            if (flagsiz eq 1) then decompima=fltarr(sizd(1),sizd(2),sizd(3),sizd(4)) else decompima=fltarr(sizd(1),sizd(2),rows,columns)
         endif 

;convert image to 16 bit
;         scimg=(img-range[0])/(range[1]-range[0]) ;0...1
         scimg=(img-range[0])/(range[1]-range[0]+(range[1]-range[0])*1e-7) ;0...1
         img16=fix(scimg*2l^16,type=12)
         img=uint(img16/2l^(16-info.comprbpp))
         if (ll eq 0 AND ss eq 0) then imgref=img
;stop
;write image compression-ready
         openw,unit,fileout+'_16bit',/get_lun
         writeu,unit,img
         free_lun,unit

;run the compression code through wine
         cmd=strcompress('wine DataComp.exe '+fileout+'_16bit '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw)+' '+string(rows)+' '+string(columns)+' '+string(info.comprratio)+' '+string(info.comprbpp)+' '+string(info.comprbps))
;stop
         spawn,cmd,result,error

;store compression ratio from results
         ratind=strpos(result,'ratio')
         if (max(ratind) ne -1) then rat[sizd(2)*ll+ss]=float(strmid(result[where(ratind ne -1)],ratind[where(ratind ne -1)]+7)) else rat[sizd(2)*ll+ss]=0
;stop
         if (info.talk ge 1) then print,format='(35(%"\b"),"Compressing...",I3," %   ratio ",F5.3,$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.),rat[sizd(2)*ll+ss] ;cmd
;      if (info.talk ge 1) then print, ;transpose(result),transpose(error)
         if max(strpos(strlowcase(result),'success')) eq -1 then message,'COMPRESS-Error: '+strjoin([cmd,error],' / ')
      endfor ;ss
   endfor ;ll

   print,''

   for ll=0,sizd(1)-1 do begin  ;info.etalpos-1 do begin 
      for ss=0,sizd(2)-1 do begin ;info.modnbuff-1 do begin

;run the decompression
;      if (info.modnbuff*ll+ss eq 3 OR info.modnbuff*ll+ss eq 15) then goto,nodec
         cmd=strcompress('wine DataComp.exe '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+fileout+'_decomp '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw))
         spawn,cmd,result,error
         if (info.talk ge 1) then print,format='(35(%"\b"),"Decompressing...",I3," %",$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.) ;cmd
;      if info.talk ge 1 then print,transpose(result),transpose(error)
;habria que hacer que imprimiese el 'ERROR' de la variable result, que
;es mas descriptivo
         if max(strpos(strlowcase(result),'success')) eq -1 then message,'DECOMPRESS-Error: '+strjoin([cmd,error],' / ')

;read decompressed image
         uimg=uintarr(rows,columns)
         openr,unit,fileout+'_decomp',/get_lun
         readu,unit,uimg
         free_lun,unit

;undo the 'image preparation'
         uimgor=float(uimg)*2l^(16-info.comprbpp)/2l^16*(rangetot[1,sizd(2)*ll+ss]-rangetot[0,sizd(2)*ll+ss])+rangetot[0,sizd(2)*ll+ss]

;recover original dimensions in case of replication because size<128 pix
         if (flagsiz eq 1) then uimgor=uimgor[0:sizd(3)-1,0:sizd(4)-1]

         decompima[ll,ss,*,*]=uimgor
;nodec:
         if (info.talk eq 1 AND ll eq 0 AND ss eq 0) then begin
            window,/free
            zer=where(imgref ne 0)
            diff=(float(uimg[zer])-float(imgref[zer]))/float(imgref[zer])
            h=histogram(diff,nbins=100,location=x)
            plot,psym=10,x,h,/xst,/yst,tit='Histogram of (decompressed image - 16bit orig image) / 16bit orig image'
         endif

      endfor ;ss
   endfor ;ll
   print,''

;save decompressed data in fits format
   sxaddpar,head,"HISTORY",'Compressed-decompressed images'
   if (info.comprhw eq 1) then sxaddpar,head,'COMPHW  ','Hardware','Hardware/software compression simulated' else sxaddpar,head,'COMPHW  ','Software','Hardware/software compression simulated'
   sxaddpar,head,'BITSPIX ',format='(F8.3)',info.comprbpp,'Bits per pixel'
   sxaddpar,head,'COMPRAT ',format='(F8.3)',mean(rat),'Compression ratio'
   writefits,fileout+'.fits',decompima,head,/compress

;------------------------------
   if (max(info.reterrors) eq 1 OR max(info.muellerrors) eq 1) then begin
      fileout=info.saves+info.files(progind)+'_teo_'+strtrim(obs,2)

      data=readfits(info.saves+info.files(progma)+'_teo_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)
      sizd=size(data)

      rat=fltarr(sizd(1)*sizd(2)) 
      rangetot=fltarr(2,sizd(1)*sizd(2)) 

      for ll=0,sizd(1)-1 do begin 
         for ss=0,sizd(2)-1 do begin 
            img=reform(data[ll,ss,*,*])
;***********************
;these lines below are needed for too small input data. Less than 128
;pixels seems not possible for compression 
;creates problems when flats/fringes present because of the intesity
;gradients, it is no longer periodic
            siz=size(img,/dim)
            if (siz(0) lt 128 OR siz(1) lt 128) then begin
               if (ll eq 0 AND ss eq 0) then print,'Data smaller than 128 pixels. Replicating...'
               flagsiz=1
               img0=img
               facrep=ceil(128./min(siz))
;            img=fltarr(siz[0]*2,siz[1]*2)
;            img[0:siz[0]-1,0:siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,0:siz[1]-1]=img0
;            img[0:siz[0]-1,siz[1]:2*siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,siz[1]:2*siz[1]-1]=img0
               img=fltarr(siz[0]*facrep,siz[1]*facrep)
               for jjy=0,facrep-1 do begin
                  for jjx=0,facrep-1 do begin
                     img[siz[0]*jjx:siz[0]*(jjx+1)-1,siz[1]*jjy:siz[1]*(jjy+1)-1]=img0
                  endfor 
               endfor 
;***********************
               siz=size(img,/dim)
            endif 

            if (min(img) ge 0) then range=[0,max(img)] else range=[min(img),max(img)] 
            rangetot[*,sizd(2)*ll+ss]=range

            siz=2l^fix(alog(siz)/alog(2))
            img=img[0:siz[0]-1,0:siz[1]-1]
            siz=size(img,/dim)
            rows=siz[0]
            columns=siz[1]

;            if (ll eq 0 AND ss eq 0) then decompima=fltarr(sizd(1),sizd(2),rows,columns) 
            if (ll eq 0 AND ss eq 0) then begin
               if (flagsiz eq 1) then decompima=fltarr(sizd(1),sizd(2),sizd(3),sizd(4)) else decompima=fltarr(sizd(1),sizd(2),rows,columns)
            endif 

            scimg=(img-range[0])/(range[1]-range[0]) ;0...1
            img16=fix(scimg*2l^16,type=12)
            img=uint(img16/2l^(16-info.comprbpp))
            openw,unit,fileout+'_16bit',/get_lun
            writeu,unit,img
            free_lun,unit

            cmd=strcompress('wine DataComp.exe '+fileout+'_16bit '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw)+' '+string(rows)+' '+string(columns)+' '+string(info.comprratio)+' '+string(info.comprbpp)+' '+string(info.comprbps))
            spawn,cmd,result,error

            ratind=strpos(result,'ratio')
            if (max(ratind) ne -1) then rat[sizd(2)*ll+ss]=float(strmid(result[where(ratind ne -1)],ratind[where(ratind ne -1)]+7)) else rat[sizd(2)*ll+ss]=0
            if (info.talk ge 1) then print,format='(35(%"\b"),"Compressing teo...",I3," %   ratio ",F5.3,$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.),rat[sizd(2)*ll+ss] ;cmd
            if max(strpos(strlowcase(result),'success')) eq -1 then message,'COMPRESS-Error: '+strjoin([cmd,error],' / ')
         endfor ;ss
      endfor ;ll
      print,''

      for ll=0,sizd(1)-1 do begin 
         for ss=0,sizd(2)-1 do begin 
            cmd=strcompress('wine DataComp.exe '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+fileout+'_decomp '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw))
            spawn,cmd,result,error
            if (info.talk ge 1) then print,format='(35(%"\b"),"Decompressing teo...",I3," %",$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.) ;cmd
            if max(strpos(strlowcase(result),'success')) eq -1 then message,'DECOMPRESS-Error: '+strjoin([cmd,error],' / ')

            uimg=uintarr(rows,columns)
            openr,unit,fileout+'_decomp',/get_lun
            readu,unit,uimg
            free_lun,unit
            uimgor=float(uimg)*2l^(16-info.comprbpp)/2l^16*(rangetot[1,sizd(2)*ll+ss]-rangetot[0,sizd(2)*ll+ss])+rangetot[0,sizd(2)*ll+ss]
;recover original dimensions in case of replication because size<128 pix
            if (flagsiz eq 1) then uimgor=uimgor[0:sizd(3)-1,0:sizd(4)-1]
            decompima[ll,ss,*,*]=uimgor
         endfor ;ss
      endfor ;ll
      print,''

      sxaddpar,head,"HISTORY",'Compressed-decompressed images'
      if (info.comprhw eq 1) then sxaddpar,head,'COMPHW  ','Hardware','Hardware/software compression simulated' else sxaddpar,head,'COMPHW  ','Software','Hardware/software compression simulated'
      sxaddpar,head,'BITSPIX ',format='(F8.3)',info.comprbpp,'Bits per pixel'
      sxaddpar,head,'COMPRAT ',format='(F8.3)',mean(rat),'Compression ratio'
      writefits,fileout+'.fits',decompima,head,/compress  

   endif ;teo

;------------------------------
   if ((size(dmodm))(0) gt 2) then begin
      fileout=info.saves+info.files(progind)+'_av_'+strtrim(obs,2)

      data=readfits(info.saves+info.files(progma)+'_av_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)
      sizd=size(data)

      rat=fltarr(sizd(1)*sizd(2)) 
      rangetot=fltarr(2,sizd(1)*sizd(2)) 

      for ll=0,sizd(1)-1 do begin 
         for ss=0,sizd(2)-1 do begin 
            img=reform(data[ll,ss,*,*])
;***********************
;these lines below are needed for too small input data. Less than 128
;pixels seems not possible for compression 
;creates problems when flats/fringes present because of the intesity
;gradients, it is no longer periodic
            siz=size(img,/dim)
            if (siz(0) lt 128 OR siz(1) lt 128) then begin
               if (ll eq 0 AND ss eq 0) then print,'Data smaller than 128 pixels. Replicating...'
               flagsiz=1
               img0=img
               facrep=ceil(128./min(siz))
;            img=fltarr(siz[0]*2,siz[1]*2)
;            img[0:siz[0]-1,0:siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,0:siz[1]-1]=img0
;            img[0:siz[0]-1,siz[1]:2*siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,siz[1]:2*siz[1]-1]=img0
               img=fltarr(siz[0]*facrep,siz[1]*facrep)
               for jjy=0,facrep-1 do begin
                  for jjx=0,facrep-1 do begin
                     img[siz[0]*jjx:siz[0]*(jjx+1)-1,siz[1]*jjy:siz[1]*(jjy+1)-1]=img0
                  endfor 
               endfor 
;***********************
               siz=size(img,/dim)
            endif 

            if (min(img) ge 0) then range=[0,max(img)] else range=[min(img),max(img)] 
            rangetot[*,sizd(2)*ll+ss]=range

            siz=2l^fix(alog(siz)/alog(2))
            img=img[0:siz[0]-1,0:siz[1]-1]
            siz=size(img,/dim)
            rows=siz[0]
            columns=siz[1]

;            if (ll eq 0 AND ss eq 0) then decompima=fltarr(sizd(1),sizd(2),rows,columns) 
            if (ll eq 0 AND ss eq 0) then begin
               if (flagsiz eq 1) then decompima=fltarr(sizd(1),sizd(2),sizd(3),sizd(4)) else decompima=fltarr(sizd(1),sizd(2),rows,columns)
            endif 

            scimg=(img-range[0])/(range[1]-range[0]) ;0...1
            img16=fix(scimg*2l^16,type=12)
            img=uint(img16/2l^(16-info.comprbpp))
            openw,unit,fileout+'_16bit',/get_lun
            writeu,unit,img
            free_lun,unit

            cmd=strcompress('wine DataComp.exe '+fileout+'_16bit '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw)+' '+string(rows)+' '+string(columns)+' '+string(info.comprratio)+' '+string(info.comprbpp)+' '+string(info.comprbps))
            spawn,cmd,result,error

            ratind=strpos(result,'ratio')
            if (max(ratind) ne -1) then rat[sizd(2)*ll+ss]=float(strmid(result[where(ratind ne -1)],ratind[where(ratind ne -1)]+7)) else rat[sizd(2)*ll+ss]=0
            if (info.talk ge 1) then print,format='(35(%"\b"),"Compressing av...",I3," %   ratio ",F5.3,$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.),rat[sizd(2)*ll+ss] ;cmd
            if max(strpos(strlowcase(result),'success')) eq -1 then message,'COMPRESS-Error: '+strjoin([cmd,error],' / ')
         endfor ;ss
      endfor ;ll
      print,''

      for ll=0,sizd(1)-1 do begin  
         for ss=0,sizd(2)-1 do begin 
            cmd=strcompress('wine DataComp.exe '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+fileout+'_decomp '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw))
            spawn,cmd,result,error
            if (info.talk ge 1) then print,format='(35(%"\b"),"Decompressing av...",I3," %",$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.) ;cmd
            if max(strpos(strlowcase(result),'success')) eq -1 then message,'DECOMPRESS-Error: '+strjoin([cmd,error],' / ')

            uimg=uintarr(rows,columns)
            openr,unit,fileout+'_decomp',/get_lun
            readu,unit,uimg
            free_lun,unit
            uimgor=float(uimg)*2l^(16-info.comprbpp)/2l^16*(rangetot[1,sizd(2)*ll+ss]-rangetot[0,sizd(2)*ll+ss])+rangetot[0,sizd(2)*ll+ss]
;recover original dimensions in case of replication because size<128 pix
            if (flagsiz eq 1) then uimgor=uimgor[0:sizd(3)-1,0:sizd(4)-1]
            decompima[ll,ss,*,*]=uimgor
         endfor ;ss
      endfor ;ll
      print,''

      sxaddpar,head,"HISTORY",'Compressed-decompressed images'
      if (info.comprhw eq 1) then sxaddpar,head,'COMPHW  ','Hardware','Hardware/software compression simulated' else sxaddpar,head,'COMPHW  ','Software','Hardware/software compression simulated'
      sxaddpar,head,'BITSPIX ',format='(F8.3)',info.comprbpp,'Bits per pixel'
      sxaddpar,head,'COMPRAT ',format='(F8.3)',mean(rat),'Compression ratio'
      writefits,fileout+'.fits',decompima,head,/compress  

   endif ;av

;to check if input data has been demodulated or not
headmod=sxpar(head,'HISTORY')
demhead=where(headmod eq 'Demodulated images')
;------------------------------
   if (info.dualbeam eq 1 AND demhead eq -1) then begin
;if there are separated data for dual-beam, there has been no
;demodulation, so no need for theoretical and average cases      
      fileout=info.saves+info.files(progind)+'_2_'+strtrim(obs,2)

      data=readfits(info.saves+info.files(progma)+'_2_'+strtrim(obs,2)+'.fits*',head,/comp,/sil)
      sizd=size(data)

      rat=fltarr(sizd(1)*sizd(2)) 
      rangetot=fltarr(2,sizd(1)*sizd(2)) 

      for ll=0,sizd(1)-1 do begin ;info.etalpos-1 do begin 
         for ss=0,sizd(2)-1 do begin ;info.modnbuff-1 do begin
            img=reform(data[ll,ss,*,*])
;***********************
;these lines below are needed for too small input data. Less than 128
;pixels seems not possible for compression 
;creates problems when flats/fringes present because of the intesity
;gradients, it is no longer periodic
            siz=size(img,/dim)
            if (siz(0) lt 128 OR siz(1) lt 128) then begin
               if (ll eq 0 AND ss eq 0) then print,'Data smaller than 128 pixels. Replicating...'
               flagsiz=1
               img0=img
               facrep=ceil(128./min(siz))
;            img=fltarr(siz[0]*2,siz[1]*2)
;            img[0:siz[0]-1,0:siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,0:siz[1]-1]=img0
;            img[0:siz[0]-1,siz[1]:2*siz[1]-1]=img0
;            img[siz[0]:2*siz[0]-1,siz[1]:2*siz[1]-1]=img0
               img=fltarr(siz[0]*facrep,siz[1]*facrep)
               for jjy=0,facrep-1 do begin
                  for jjx=0,facrep-1 do begin
                     img[siz[0]*jjx:siz[0]*(jjx+1)-1,siz[1]*jjy:siz[1]*(jjy+1)-1]=img0
                  endfor 
               endfor 
;***********************
               siz=size(img,/dim)
            endif 

            if (min(img) ge 0) then range=[0,max(img)] else range=[min(img),max(img)] 
            rangetot[*,sizd(2)*ll+ss]=range

            siz=2l^fix(alog(siz)/alog(2))
            img=img[0:siz[0]-1,0:siz[1]-1]
            siz=size(img,/dim)
            rows=siz[0]
            columns=siz[1]

;            if (ll eq 0 AND ss eq 0) then decompima=fltarr(sizd(1),sizd(2),rows,columns) 
            if (ll eq 0 AND ss eq 0) then begin
               if (flagsiz eq 1) then decompima=fltarr(sizd(1),sizd(2),sizd(3),sizd(4)) else decompima=fltarr(sizd(1),sizd(2),rows,columns)
            endif 

            scimg=(img-range[0])/(range[1]-range[0]) ;0...1
            img16=fix(scimg*2l^16,type=12)
            img=uint(img16/2l^(16-info.comprbpp))
            openw,unit,fileout+'_16bit',/get_lun
            writeu,unit,img
            free_lun,unit

            cmd=strcompress('wine DataComp.exe '+fileout+'_16bit '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw)+' '+string(rows)+' '+string(columns)+' '+string(info.comprratio)+' '+string(info.comprbpp)+' '+string(info.comprbps))
            spawn,cmd,result,error

            ratind=strpos(result,'ratio')
            if (max(ratind) ne -1) then rat[sizd(2)*ll+ss]=float(strmid(result[where(ratind ne -1)],ratind[where(ratind ne -1)]+7)) else rat[sizd(2)*ll+ss]=0
            if (info.talk ge 1) then print,format='(35(%"\b"),"Compressing dual...",I3," %   ratio ",F5.3,$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.),rat[sizd(2)*ll+ss] ;cmd
            if max(strpos(strlowcase(result),'success')) eq -1 then message,'COMPRESS-Error: '+strjoin([cmd,error],' / ')
         endfor ;ss
      endfor ;ll
      print,''

      for ll=0,sizd(1)-1 do begin 
         for ss=0,sizd(2)-1 do begin 
            cmd=strcompress('wine DataComp.exe '+fileout+'_comp_'+strtrim(sizd(2)*ll+ss,2)+' '+fileout+'_decomp '+string(info.comprskip)+' '+string(info.comprmsb)+' '+string(info.comprhw))
            spawn,cmd,result,error
            if (info.talk ge 1) then print,format='(35(%"\b"),"Decompressing dual...",I3," %",$)',fix(float((sizd(2)*ll+ss+1))/(sizd(1)*sizd(2))*100.) ;cmd
            if max(strpos(strlowcase(result),'success')) eq -1 then message,'DECOMPRESS-Error: '+strjoin([cmd,error],' / ')

            uimg=uintarr(rows,columns)
            openr,unit,fileout+'_decomp',/get_lun
            readu,unit,uimg
            free_lun,unit
            uimgor=float(uimg)*2l^(16-info.comprbpp)/2l^16*(rangetot[1,sizd(2)*ll+ss]-rangetot[0,sizd(2)*ll+ss])+rangetot[0,sizd(2)*ll+ss]
;recover original dimensions in case of replication because size<128 pix
            if (flagsiz eq 1) then uimgor=uimgor[0:sizd(3)-1,0:sizd(4)-1]
            decompima[ll,ss,*,*]=uimgor

         endfor ;ss
      endfor ;ll
      print,''

      sxaddpar,head,"HISTORY",'Compressed-decompressed images'
      if (info.comprhw eq 1) then sxaddpar,head,'COMPHW  ','Hardware','Hardware/software compression simulated' else sxaddpar,head,'COMPHW  ','Software','Hardware/software compression simulated'
      sxaddpar,head,'BITSPIX ',format='(F8.3)',info.comprbpp,'Bits per pixel'
      sxaddpar,head,'COMPRAT ',format='(F8.3)',mean(rat),'Compression ratio'
      writefits,fileout+'.fits',decompima,head,/compress   

   endif ;dual

;remove all the intermediate files
spawn,'rm '+fileout+'_16bit*'
spawn,'rm '+fileout+'_comp*'
spawn,'rm '+fileout+'_decomp*'

endfor ;obs


;stop
;if (info.talk eq 1) then begin
;   window,/free
;   zer=where(img ne 0)
;   diff=(float(uimg[zer])-float(img[zer]))/float(img[zer])
;;   imgn=data[0,0,*,*]
;;   zer=where(imgn ne 0)
;;   diff=()
;   h=histogram(diff,nbins=100,location=x)
;   plot,psym=10,x,h,/xst,/yst,tit='Histogram of (decompressed image - 16bit orig image) / 16bit orig image'
;endif

end
