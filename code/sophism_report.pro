pro sophism_report

; ==============================================================================
;
; DESCRIPTION
;    Generates a simulation report, including parameters, graphs,
;    images,...
;
; CALLED FROM
;    sophism
;
; SUBROUTINTES
;
; MAIN INPUTS
;    Settings file, and save files and fits from each module run.
;
; OUTPUT
;    Report in pdf file.
; 
; VERSION HISTORY 
;    Alex Feller   v0.1    2010-09-21
;    J. Blanco. 2012. v1_0. Adapted for sophism. Changed to printf mode
;       adapting from J. Hirzberger 
;    J. Blanco. Nov 2012. v2_0. Major changes in most parts. Images from
;       inversion part included now
;    J. Blanco. Dec 2012. v2_1. Modified to use the 'report modules' settings
;    J. Blanco. Sep 2014. v3_0. Updating whole routine to reflect
;       latest changes in simulator
;    J. Blanco. Feb 2016. v3_1. Corrected for case of various
;       observation sets in spectral positions list and demodulation 
;
; ==============================================================================


; ------------------------------------------------------------------------------
; Initialize
; ------------------------------------------------------------------------------

restore,'../settings/settings_sophism.sav'
print,''
print,'sophism_report'
print,''
;change the image plotting to ps
set_plot,'ps'

;start the text file with the name given on settings
barr=strpos(info.freport,'/',/reverse_search)+1
dott=strpos(info.freport,'.',/reverse_search)
dir_rep=strmid(info.freport,0,barr);+'/'
fil_rep=strmid(info.freport,barr,dott-barr)
get_lun,unit
openw,unit,dir_rep+fil_rep+'.tex'

;start with the header of the tex document. Title page as well.
printf,unit,'\documentclass[a4paper,10pt]{article}'
printf,unit,'\usepackage{graphicx}'

printf,unit,'\addtolength{\textwidth}{3cm}'
printf,unit,'\addtolength{\textheight}{2cm}'
printf,unit,'\addtolength{\topmargin}{-1.5cm}'
printf,unit,'\addtolength{\oddsidemargin}{-2cm}'

printf,unit,'\title{SOPHISM REPORT}'
;get the author directly from the session. Not the easiest way probably
spawn,'who -m > auth'
openr,unit2,'auth',/get_lun
auth=strarr(1)
readf,unit2,auth
close,unit2
spawn,'rm auth'
auth=strmid(auth,0,strpos(auth,'pts'))
auth=strmid(auth,0,strpos(auth,' '))
printf,unit,'\author{'+auth+'}'

printf,unit,'\begin{document}'
printf,unit,'\maketitle'

; ------------------------------------------------------------------------------
; Simulation settings
; ------------------------------------------------------------------------------
printf,unit,'\section*{Global Settings}'

printf,unit,'Distance S/C - Sun: ' + string(format='(F5.3)', info.scdis) + ' AU \\'
printf,unit,'S/C Velocity: ' + string(format='(F6.3)', info.scvel) + ' km/s \\'
;printf,unit,'Latitude: ',string(format='(F6.2)',info.sclat),' deg \\'
;printf,unit,' Temperature: '+string(format='(F5.2)',info.temp)+' $^{\circ}$C \\'

;printf,unit, 'Telescope and detector: \\'
printf,unit,'Telescope aperture: ' + string(format='(F4.2)', info.telap) + ' m \\'

printf,unit,'Acquisition Scheme: '+info.acqsche+' \\'

printf,unit,'Temporal settings: \\'
printf,unit,'\indent Exposure time: ',string(format='(F6.2)',info.etim),' ms \\'
printf,unit,'\indent Detector readout time: ',string(format='(F6.2)',info.routtim),' ms \\'
if (info.discar eq 0) then discopt='Samples' else discopt='Frames'
printf,unit,'\indent Discarding scheme: '+discopt+' \\'
printf,unit,'\indent Observation Sets: ',string(format='(F6.2)',info.obsset),' \\'
printf,unit,'Full observation time: ' + string(format='(F6.2)', info.tobs) + ' s \\'

if (info.dualbeam eq 1) then db='ON' else db='OFF'
printf,unit,'Dual-beam mode: '+db+' \\'

if (info.routines(7) eq 1) then printf,unit,'Number of accumulations: ',strtrim(info.nacc,2),' \\'

; ------------------------------------------------------------------------------
; Original solar scene
; ------------------------------------------------------------------------------

printf,unit,'\newpage'

printf,unit,'\section*{Original Solar Scene}'
printf,unit,'Original MHD solar scene: \\' 

;printf,unit,'Quiet sun, B\_mean=30G, blue wing of Fe I 630.2 nm (N. Vitas) \\'
datos=strsplit(info.datafilo,'_',/extract)
if (n_elements(datos) gt 1) then begin
   name=datos(0)
   for dd=1,n_elements(datos)-1 do name=name+'\_'+datos(dd)
endif 
printf,unit,'Data file(s): '+name+' \\'

;maybe not a good idea. Not necessarily run the etalon module
;restore,info.saves+info.files(3)+'.sav'
printf,unit,'Wavelength range (\AA): '+string(format='(F9.3)',min(info.lll))+' - '+string(format='(F9.3)',max(info.lll))+' (Doppler shift from S/C velocity not included) \\'
printf,unit,'Wavelength sampling (m\AA): ',info.wavesamp,' \\'

printf,unit,'Spatial dimensions: \\'
printf,unit,'\indent Replication: '+strtrim(fix(info.sscene_mag),2)+' \\'
printf,unit,'\indent ',string(format='("Image size: ", I3, " x ", I3, " pixels")', info.sz, info.sz),' \\'
printf,unit,'Spatial sampling: \\'
printf,unit,'\indent On the Sun: ' + string(format='(F5.1)', info.sssamp0) + ' km/pixel \\'
printf,unit,'\indent On the sky: ' + string(format='(F5.3)', (info.sssamp0*1e3)*(1./(info.scdis*1.495e11))*!radeg*3600.) + ' arcsec \\'
if (info.fparesamp eq 1) then printf,unit,'\indent On the detector: '+string(format='(F5.3)', info.sssamp)+' arcsec \\'

printf,unit,'Time Series: \\'
if (info.ntim gt 1) then begin
   if (info.numfil gt 1) then  begin
      printf,unit,'\indent Time Evolution \\'
      printf,unit,'\indent Time Sampling: \\'
      printf,unit,'\indent \indent ',string(format='("  Original: ", F5.3, " s")', info.tsampor),' \\'
      if (info.tsamp ne info.tsampor) then printf,unit,'\indent \indent ',string(format='("  Interpolated: ", F5.3, " s")', info.tsamp),' \\'
   endif else begin
      printf,unit,'\indent No Time Evolution'
      printf,unit,'\indent Artificial Time Sampling: '+string(format='(F4.2)',info.tsamp)+' ms \\'
   endelse 
endif else printf,unit,'\indent Single time sample'


;esto solo para el primero de la serie, haya mas o no.
nf=0
ima=readfits(info.saves+info.files(0)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
sz=size(ima)
;muy lejos estamos en continuo puro y no hay stokes, solo I
cc=info.wavepoint/2.-round(200./info.wavesamp)
if (cc lt 0) then cc=0
cont=reform(ima(cc,0,*,*))

loadct, 0, /s
device,filename=dir_rep+'input.eps',bit=64,/color,/encapsulated
si = reform(ima(cc,0,*,*))/mean(cont)
if (sz(2) lt 4) then tvframe,si,/bar,/aspect,btitle='I!DC!N' else begin
   sv = reform(ima(cc,1,*,*))/cont
   mm=max(abs(sv))
   sv[0,0]=mm
   sv[0,1]=-mm
   sq = reform(ima(cc,2,*,*))/cont
   mm=max(abs(sq))
   sq[0,0]=mm
   sq[0,1]=-mm
   su = reform(ima(cc,3,*,*))/cont
   mm=max(abs(su))
   su[0,0]=mm
   su[0,1]=-mm
   !p.multi = [0,2,2]
   tvframe,si,/bar,/aspect,btitle='I!DC!N'
   tvframe,sq,/bar,/aspect,btitle='Q/I!DC!N'
   tvframe,su,/bar,/aspect,btitle='U/I!DC!N'
   tvframe,sv,/bar,/aspect,btitle='V/I!DC!N'
   device,/close
endelse

printf,unit,'\begin{center}'
printf,unit,'\begin{figure}[!h]'
printf,unit,'\includegraphics[width=13cm]{'+dir_rep+'input.eps}'
printf,unit,'\caption{Stokes parameters from original input data at $-$'+strtrim(fix((round(400./info.wavesamp) < info.wavepoint/2.)*info.wavesamp),2)+' m\AA\ from center.}'
printf,unit,'\end{figure}'
printf,unit,'\end{center}'
!p.multi=0


; ------------------------------------------------------------------------------
; Jitter
; ------------------------------------------------------------------------------

;if (info.routines(1) eq 1) then begin
if (info.repomodules(1) eq 1) then begin
   printf,unit,'\newpage'
   printf,unit,'\section*{Jitter}'

   restore,info.saves+info.files(1)+'.sav'
   restore,info.saves+info.files(1)+'_ftpfilt.sav'

   printf,unit,'Jittering generation: '+info.jittype+' \\'

   printf,unit,'Frequency filter: \\'
   printf,unit,'\indent Type: ' + info.filtype+' \\'

   case info.filtype of
      'cutoff': begin
         printf,unit, '\indent Max. frequency: ' + string(format='(F4.1)', info.filcutoff) + ' Hz \\'
      end ;cutoff

      'cutout': begin
         printf,unit,'\indent Min. frequency: ' + string(format='(F4.1)', info.filcut1) + ' Hz \\'
         printf,unit,'\indent Max. frequency: ' + string(format='(F4.1)', info.filcut2) + ' Hz \\'
      end ;cutout

      'levels': begin
         nrange = n_elements(info.filrange)/2
         for nn=0,nrange-1 do begin
            printf,unit,'\indent Range from '+info.filrange[2*nn]+' to '+info.filrange[2*nn+1]+' Hz -- Power: '+info.fillevel[nn]^2 / (info.filrange[2*nn+1]-info.filrange[2*nn])
         endfor 
      end ;levels

      'Hinode': begin
            printf,unit,'\indent Min. frequency: ' + string(format='(F4.1)', info.filcut1) + ' Hz \\'
            printf,unit,'\indent Max. frequency: ' + string(format='(F4.1)', info.filcut2) + ' Hz \\'
      end ;hinode

   endcase 
   
   if info.filrms then a = 'AFTER filtering' else a = 'BEFORE filtering'
   printf,unit,'Jittering rms '+a+': '+string(format='(F5.2)',info.jitrms)+' arcsec \\'

   if (info.isson eq 1) then begin
      printf,unit,'ISS active \\'
      printf,unit,'\indent Attenuation: ' + info.issatten+' \\'
   endif else printf,unit,'ISS INactive \\'

   loadct, 0, /s
   !p.multi = [0, 1, 3]
; Plot jitter in x direction
   xtitle = 'time [s]'
   ytitle = 'x [arcsec]'

   device,filename=dir_rep+'jitter.eps',/encapsulated,bit=64,/color
   
   plot, jitt.t, jitt.x, xtitle=xtitle, ytitle=ytitle,psy=-2

; Plot jitter in y direction
   xtitle = 'time [s]'
   ytitle = 'y [arcsec]'

   plot, jitt.t, jitt.y, xtitle=xtitle, ytitle=ytitle,psy=-2

; Plot power spectrum filter
   ptitle = 'Filter'
   if (info.isson eq 1) then ptitle = 'Filter (solid) and attenuation (dashed)'
   xtitle = 'frequency [Hz]'
   ytitle = 'power'

   s = sort(tpfilt.nu)

   plot, tpfilt[s].nu, tpfilt[s].pow, subtitle=ptitle, xtitle=xtitle, ytitle=ytitle, /xs, /ylog,yr=[1e-8,5], /ys

   if (info.isson eq 1) then begin
      restore,info.saves+info.files(1)+'_iss.sav'
      s = sort(attenarr.nu)
      attenarr.nu = attenarr[s].nu
      attenarr.f = attenarr[s].f
      w = where(abs(attenarr.nu) gt 0)
      y = 1./attenarr[w].f^2.
      oplot, attenarr[w].nu, y, linestyle=2
   endif
   device,/close
   printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'jitter.eps} \\'

   printf,unit,'RMS values in arcsec: \\'
   printf,unit, '\indent x = ',stddev(jitt.x),' \\'
   printf,unit, '\indent y = ',stddev(jitt.y),' \\'

   !p.multi=0

; ------------------------------------------------------------------------------
; Jitter power spectra
; ------------------------------------------------------------------------------

   printf,unit,'\subsection*{Jitter Power Spectra}'
; Compute jitter power spectra
;goto,powerfin
   jittpow = {nu: 0., x: 0., y: 0.}
   ntim = info.ntim
   if ntim mod 2 eq 0 then nfreq = ntim/2 else nfreq = (ntim-1)/2
   jittpow = replicate(jittpow, nfreq)
   jittpow.nu = indgen(nfreq) * (1./info.tobs)
   for i=1,2 do begin
      a = jitt.(i)
      a = (fft(a,-1))[0:nfreq-1]
      a = abs(conj(a)*a)
; normalization to int(a) = stddev(jitter)^2
      b = int_tabulated(jittpow.nu, a)
      c = stddev(jitt.(i))^2.
      a = a * c / b
      jittpow.(i) = a
   endfor

; Plot jitter power spectra
   !p.multi = [0, 1, 2]
   !p.charsize = 0
   xtitle = 'Frequency [Hz]'
   ytitle = ['x power [arcsec^2/Hz]', 'y power [arcsec^2/Hz]']
   device,filename=dir_rep+'jitterpower.eps',/encapsulated
   for i=1,2 do plot, jittpow.nu, jittpow.(i),xtitle=xtitle, ytitle=ytitle[i-1], /ylog, yr=[1e-6,1]
   device,/close

   printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'jitterpower.eps}'
powerfin:

   !p.multi=0
endif ;jitter--info.routines(1)


; ------------------------------------------------------------------------------
; Polarization modulation
; ------------------------------------------------------------------------------

;if (info.routines(2) eq 1) then begin
if (info.repomodules(2) eq 1) then begin
   restore,info.saves+info.files(2)+'_scheme.sav'
   restore,info.saves+info.files(2)+'_scheme_demod.sav'
   restore,info.saves+info.files(2)+'_time.sav'

   printf,unit,'\newpage'
   printf,unit,'\section*{Polarization Modulation}'
;Mueller matrix
   printf,unit,'Number of cycle repetitions: ',info.cycrep,' \\'

   printf,unit,'Mueller matrix type: ',info.mueller,' \\'
   if (info.mueller ne 'Identity') then begin
      printf,unit,'\[ \left( \begin{array}{cccc}'
      for mm=0,2 do begin 
         printf,unit,string(format='(3(F5.3,"&"),F5.3,"\\")',muellmat[*,mm])
      endfor
      printf,unit,string(format='(3(F5.3,"&"),F5.3,"\end{array} \right)\]")',muellmat[*,3])
   endif

   printf,unit,'Modulation Type: ' + info.modtype+' \\'
   spar = ['I', 'Q', 'U', 'V']

   printf,unit,'  Modulation period: ' + string(format='(F5.2)', info.nperiod(0)*info.tsamp) + ' s \\'

;   printf,unit,'  Demodulation type: ' + info.moddmod_type+' \\'
   printf,unit,'Demodulation matrix: \\'

;longitudinal case
   if (info.modtype eq 'longitudinal') then begin
      a = ['I', 'V']
      printf,unit,'\[ \left( \begin{array}{cc}'
      printf,unit,string(format='(F6.3,"&",F6.3,"\\")',dmodm[*,0])
      printf,unit,string(format='(F6.3,"&",F6.3,"\end{array} \right)\]")',dmodm[*,1])
   endif else begin
;ideal and real cases
      a = ['I', 'Q', 'U', 'V']
      printf,unit,'\[ \left( \begin{array}{cccc}'
      for mm=0,2 do printf,unit,string(format='(3(F6.3,"&"),F6.3,"\\")',dmodm[*,mm])
 printf,unit,string(format='(3(F6.3,"&"),F6.3,"\end{array} \right)\]")',dmodm[*,3])
   endelse ;modnbuff 4

;in the real case, when also random errors are enabled, we have the
;theoretical demodulation matrix
   if (info.modtype eq 'real' AND max(info.reterrors) eq 1) then begin
      printf,unit,'Theoretical demodulation matrix (without considering random errors): \\'
      printf,unit,'\[ \left( \begin{array}{cccc}'
      for mm=0,2 do printf,unit,string(format='(3(F6.3,"&"),F6.3,"\\")',dmodm_teo[*,mm])
 printf,unit,string(format='(3(F6.3,"&"),F6.3,"\end{array} \right)\]")',dmodm_teo[*,3])
   endif 

   loadct, 0, /s
   !p.multi = [0, 1, info.modnbuff/info.cycrep]
   !x.omargin=[0,10]
   !y.omargin=[0,5]
   xtitle = 'time [s]'
   yr = [-1.1,1.1]
   prev=0
   delim=[0]
   for ddlc=1,info.modnbuff do begin
      delim=[delim,prev+info.nstate,prev+info.nstate+info.ndeaths(ddlc-1,0)]
      prev=prev+info.nstate+info.ndeaths(ddlc-1,0)
   endfor 
   delim=delim[1:*]
   delim=delim*info.tsamp-info.tsamp/2.

   postitle=fltarr(4)
   postitle(0)=delim[0]/2.
   postitle(1)=delim[1]+(delim[2]-delim[1])/2.
   postitle(2)=delim[3]+(delim[4]-delim[3])/2.
   postitle(3)=delim[5]+(delim[6]-delim[5])/2.
   uptitle=['I!D1!N', 'I!D2!N', 'I!D3!N', 'I!D4!N']
   ytitle = ['I', 'Q', 'U', 'V']
   device,filename=dir_rep+'modscheme.eps',/encapsulated
   for i=0,info.modnbuff/info.cycrep-1 do begin
      plot,tim,modf[i,*],yr=yr,/ys,xtitle=xtitle,ytitle=ytitle[i],psy=10
;      plot,tim,tim,yr=yr,/ys,xtitle=xtitle,ytitle=ytitle[i],psy=10
      if (i eq 0) then xyouts,postitle,replicate(1.25,info.modnbuff/info.cycrep),uptitle,/data,charsiz=0.8
      for j=0,n_elements(delim)-1 do plots,[delim(j),delim(j)],yr,lines=1
   endfor 
   device,/close
   printf,unit,'\includegraphics[width=12cm]{'+dir_rep+'modscheme.eps} \\'
   !p.multi=0
   !x.omargin=0
   !y.omargin=0 




;   !p.multi = [0, 1, info.modnbuff]
;   xtitle = 'time [s]'
;   yr = [-1.1,1.1]
;   laten=where(modf.i ne 0)
;   tdeaths=[info.tdeath1,info.tdeath2,info.tdeath3,info.tdeath4]
;   ndeaths=fix(round(tdeaths*1e-3/info.tsamp))
;   delim=[info.nstate,info.nstate+ndeaths(0)]
;   for dd=2,info.modnbuff do delim=[delim,info.nstate*dd+total(ndeaths(0:dd-2)),info.nstate*dd+total(ndeaths(0:dd-1))]
;   delim2=delim
;   for dd=1,info.etalpos-1 do delim2=[delim2,delim+max(delim2)]
;   delim2=delim2*info.tsamp-info.tsamp/2.
;   if (info.modnbuff eq 4) then begin
;      postitle=fltarr(4)
;      postitle(0)=delim2[0]/2.
;      postitle(1)=delim2[1]+(delim2[2]-delim2[1])/2.
;      postitle(2)=delim2[3]+(delim2[4]-delim2[3])/2.
;      postitle(3)=delim2[5]+(delim2[6]-delim2[5])/2.
;      uptitle=['I!D1!N', 'I!D2!N', 'I!D3!N', 'I!D4!N']
;      ytitle = ['I', 'Q', 'U', 'V']
;   endif else begin
;      postitle=fltarr(2)
;      postitle(0)=delim2[0]/2.
;      postitle(1)=delim2[1]+(delim2[2]-delim2[1])/2.
;      uptitle=['I!D1!N', 'I!D2!N']
;      ytitle = ['I', 'V']
;   endelse 
;   printf,unit,'Modulation scheme: \\'
;   printf,unit,''
;   device,filename=dir_rep+'modscheme.eps',/encapsulated
;   for i = 1, info.modnbuff do begin
;      plot, modf(laten).t, modf(laten).(i), yr=yr, /ys, xtitle=xtitle, ytitle=ytitle[i-1],psy=10
;      if (i eq 1) then xyouts,postitle,replicate(1.25,info.modnbuff),uptitle,/data,charsiz=0.8
;      for j=0,n_elements(delim2)-1 do plots,[delim2(j),delim2(j)],yr,lines=1
;   endfor
;   device,/close
;;podria hacerse lo de unir, cuando haya ndeath=/=0, con linea no
;;histograma. Haciendo histograma antes, histograma despues (con
;;noerase), y linea normal entre medias. Un poco pesado para hacerlo bonito
;   printf,unit,'\includegraphics[width=12cm]{'+dir_rep+'modscheme.eps} \\'
;   !p.multi=0

endif ;info.routines(2)




; ------------------------------------------------------------------------------
; Spectral curves & parameters
; ------------------------------------------------------------------------------

;if (info.routines(3) eq 1) then begin
if (info.repomodules(3) eq 1) then begin
   restore,info.saves+info.files(3)+'.sav'
   szff=size(fff)
   printf,unit,'\newpage'
   printf,unit,'\section*{Spectral Profile}'
   printf,unit,'  Prefilter FWHM (A): '+string(format='(F5.3)',info.pfhwb)+' \\'
   printf,unit,'  Prefilter Cavities: '+strtrim(info.ncav1,2)+' \\'
   printf,unit,'  Etalon Thickness (micron): '+strtrim(info.spectsizd,2)+' \\'

;if pupil apodization was enabled, the transmission curve has four
;dimensions
;   if (info.routines(5) eq 1) then begin
   if (szff(0) eq 4) then begin
;just to show different pupil intensities. In order not to
;overcomplicate it, just taking the first selected wavelength (without
;even interpolating in case of neighbouring wavelengths) and
;positions at +-fwhm, +-fwhm/2.
      fwhmpos=round(fwhm_real(0)/(info.wavesamp*1d-3))
      siz=size(fff)
      pupint=fltarr(siz(1)*5,siz(2))
      pupint[2*siz(2):3*siz(2)-1,*]=fff[*,*,wavesind(0),0]
      pupint[0:siz(2)-1,*]=fff[*,*,wavesind(0)-fwhmpos,0]
      pupint[siz(2):2*siz(2)-1,*]=fff[*,*,wavesind(0)-fwhmpos/2.,0]
      pupint[3*siz(2):4*siz(2)-1,*]=fff[*,*,wavesind(0)+fwhmpos/2.,0]
      pupint[4*siz(2):5*siz(2)-1,*]=fff[*,*,wavesind(0)+fwhmpos,0]
;for all plotting purposes, use the pupil-averaged etalon transmission
      fff=fffav
      loadct,39
      device,filename=dir_rep+'pupilintens.eps',/color,bit=64,/encapsulated,xsiz=20,ysiz=5
      tvframe,pupint/max(pupint),/asp,/bar
      device,/close
      
      printf,unit,'Pupil Intensities: \\'
      printf,unit,'\begin{figure}[!h]'
      printf,unit,'\begin{center}'
      printf,unit,'\includegraphics[width=15cm]{'+dir_rep+'pupilintens.eps}\\'
      printf,unit,'\vspace{-0.7cm}'
      printf,unit,'\caption{Pupil intensities at -FWHM, -FWHM/2, center, FWHM/2 and FWHM of etalon transmission.}'
      printf,unit,'\end{center}'
      printf,unit,'\end{figure}'
      printf,unit,''
      loadct,0
   endif 


;   llopre=llo
;   llo=fltarr(info.etalpos)
;   fwhmpre=fwhm_real
;   fwhm_real=fltarr(info.etalpos)
   fsrpre=fsr0
   fsr0=fltarr(info.etalpos)
   ffpre=ff0
   ff0=fltarr(info.etalpos)
   voltpre=volt
   volt=fltarr(info.etalpos)
   fffpre=fff
   fff=fltarr(info.wavepoint,info.etalpos)
   mind=max(index)
   for ii=0,mind do begin
      indy=where(index eq ii,nindy)
      for bb=0,nindy-1 do begin
;         llo[ii]=llo[ii]+double(wg[indy(bb)]*llopre[indy(bb)])
;         fwhm_real[ii]=fwhm_real[ii]+double(wg[indy(bb)]*fwhmpre[indy(bb)])
         fsr0[ii]=fsr0[ii]+double(wg[indy(bb)]*fsrpre[indy(bb)])
         ff0[ii]=ff0[ii]+double(wg[indy(bb)]*ffpre[indy(bb)])
         volt[ii]=volt[ii]+double(wg[indy(bb)]*voltpre[indy(bb)])
         fff[*,ii]=fff[*,ii]+double(wg[indy(bb)]*fffpre[*,indy(bb)])
      endfor
   endfor

   printf,unit,'\noindent Wavelength positions: '+strtrim(info.etalpos,2)+' \\'
;   for ww=0,n_elements(info.wavearr)-1 do begin
   for ww=0,n_elements(info.wavearr/info.obsset)-1 do begin
      printf,unit,'\textbf{Position '+strtrim(ww+1,2)+'} \\'
      printf,unit,'\indent Applied Voltage: '+strtrim(volt[ww],2)+' \\' 
      if (info.voltjit eq 1) then printf,unit,'\indent \indent including tuning jittering (rms): '+string(format='(F5.3)',info.voltjitrms)+' V \\'
      printf,unit,'\indent Effective transmitted lambda (A)='+strtrim(llo[ww],2)+' \\'
      printf,unit,'\indent Total Finesse: '+strtrim(ff0[ww],2),'\\'
      printf,unit,'\indent Free Spectral Range (A): '+strtrim(fsr0[ww],2)+' \\'
      printf,unit,'\indent Real FWHM (A)='+strtrim(fwhm_real[ww],2)+' \\'
      ddb=where(lll ge llo[ww]-fwhm_real[ww]*2 and lll le llo[ww]+fwhm_real[ww]*2)
      ddl=where(lll lt llo[ww]-fwhm_real[ww]*2)         
      ddg=where(lll gt llo[ww]+fwhm_real[ww]*2)        
      printf,unit,'\indent Spectral Stray-light (\%)=',100*total([fff(ddl,ww),fff(ddg,ww)])/total(fff(ddb,ww)),' \\'
   endfor

   llo_f=info.wl
   dop=info.scvel/299792.5d0
   ldop=llo_f*dop
   llo_f=llo_f+ldop
   xtitle='!7k!3 ('+STRING("305B)+')'
   ytitle=''

   restore,'../data/fts6173.save'
   iii=iii*1.d-4
   iii=fft_shift(iii,ldop/mean(deriv(xxx)))
   iip=interpol(iii,xxx,lll)

   device,filename=dir_rep+'specplot.eps',/encapsulated,xsiz=17.78,ysiz=12.7
   plot,lll,fff[*,0],xrange=[fix(llo_f-4.),fix(llo_f+4.)],yrange=[0,1],xtitle=xtitle, ytitle=ytitle
   oplot,lll,iip
   oplot,lll,filter,lines=1
   plots,[llo(0),llo(0)],[0,1],lines=2
   plots,[llo(0)-fsr0(0),llo(0)-fsr0(0)],[0,1],lines=2
   plots,[llo(0)+fsr0(0),llo(0)+fsr0(0)],[0,1],lines=2
   capt='etalon transmission curve'
   if (info.etalpos gt 1) then begin
      oplot,lll,fff[*,info.etalpos-1],lines=3
      plots,[llo(info.etalpos-1),llo(info.etalpos-1)],[0,1],lines=3
      plots,[llo(info.etalpos-1)-fsr0(info.etalpos-1),llo(info.etalpos-1)-fsr0(info.etalpos-1)],[0,1],lines=3
      plots,[llo(info.etalpos-1)+fsr0(info.etalpos-1),llo(info.etalpos-1)+fsr0(info.etalpos-1)],[0,1],lines=3
      capt='etalon transmission curves from two tuned wavelengths'
   endif
   device,/close
   printf,unit,'\begin{figure}[!h]'
   printf,unit,'\begin{center}'
   printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'specplot.eps}'
   printf,unit,'\caption{Plot of FTS-atlas spectrum of Fe {\sc i} 6173 \AA\ line (Doppler-shifted), prefilter curve (dotted line) and  '+capt+'. Vertical lines mark transmission maximum and secondary peaks.}'
   printf,unit,'\end{center}'
   printf,unit,'\end{figure}'

endif ;info.routines(3)


; ------------------------------------------------------------------------------
; OTF-convolved images
; ------------------------------------------------------------------------------

;if (info.routines(4) eq 1) then begin
if (info.repomodules(4) eq 1) then begin
   printf,unit,'\newpage'
   printf,unit,'\section*{OTF}'
   printf,unit,'Image apodisation (\%): ',string(format='(F5.2)',info.sscene_apod),' \\'

   llo_f=info.wl
   dop=info.scvel/299792.5d0
   ldop=llo_f*dop
   llo_f=llo_f+ldop

   printf,unit,'Type: ' + info.otftype+' \\'
   cutfreq=(info.telap/(llo_f*1.e-10))*(1./(!radeg*3600.))
   deltanu=1./(info.sz*info.sssamp)
   rpup=cutfreq/(2.*deltanu)
   printf,unit,'Cutoff frequency: ' +string(format='(F5.3)',cutfreq)+ ' arcsec$^{-1}$ \\'
   if (2.*rpup gt (info.sz/2.-2)) then printf,unit,'Aliasing \\'

   if (info.zercoload eq 1) then begin
      datos=strsplit(info.zercofile,'_',/extract)
      if (n_elements(datos) gt 1) then begin
         name=datos(0)
         for dd=1,n_elements(datos)-1 do name=name+'\_'+datos(dd)
      endif 
      printf,unit,'Zernikes loaded from: '+name+' \\'
   endif 

   if (max(abs(info.zerco)) gt 0) then printf,unit,'Zernike coefficients input (rad): ',string(format='(F6.2)',info.zerco),' \\'
   if (info.frontnorm eq 1) then printf,unit,'Wavefront normalization rms: '+string(format='(F5.3)',info.frontrms)+' lambda \\'

   xpsf=(findgen(info.sz)-info.sz/2)*info.sssamp

   restore,info.saves+info.files(4)+'.sav'
   siz=size(otfs)

;if pupil apodization was enabled, the transmission curve has four
;dimensions
;case no-pupil apod
;   if (info.routines(5) ne 1) then begin 
   if (siz(0) ne 4) then begin
      nf=0
      ima=readfits(info.saves+info.files(4)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)

      !p.multi=[0,2,1]
      device,filename=dir_rep+'psf_ima.eps',/encapsulated,bit=64,/color
      psf=float(shift(fft(otfs,1),info.sz/2,info.sz/2))
      tvframe,psf/max(psf)<0.3,/bar,/aspect,btit='PSF norm'
      if (cc gt (size(ima))(1)) then cc=0
      tvframe,ima(cc,0,*,*),/bar,/aspect,btit='I!DC!N'
      device,/close
      printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'psf_ima.eps}\\'     
      !p.multi=0
print,'mirar interpolacion para psf, que se vea bien'

      device,filename=dir_rep+'psf_plot.eps',/encapsulated
      plot,xpsf,psf[*,info.sz/2]/max(psf),psy=-2,xr=[-2,2],yr=[0,1],ytit='PSF'
      device,/close
      printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'psf_plot.eps}'


;case pupil apod

   endif else begin

      if (info.phases eq 1) then printf,unit,'Phases from etalon included\\'

      restore,info.saves+info.files(3)+'.sav'
     
;monochromatic PSFs
print,'mirar interpolacion para psf, que se vea bien'

      device,filename=dir_rep+'psf_plot_mono.eps',/encapsulated
      for oo=0,n_elements(wavesind)-1 do begin
         mpsf=float(shift(fft(otfs[*,*,wavesind(oo),oo],1),info.sz/2,info.sz/2))
         plot,xpsf,mpsf(*,info.sz/2)/max(mpsf),psy=-2,xr=[-2,2],yr=[0,1],xtit='arcsec',ytit='PSF!Dnorm!N',tit='Monochromatic PSFs',lines=oo,/noerase
      endfor
      device,/close
      printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'psf_plot_mono.eps}\\'


;polychromatic PSFs
      ppsf = fltarr(info.sz,info.sz,n_elements(wavesind))
      prof=readfits(info.saves+info.files(0)+'_'+strtrim(0,2)+'.fits*',/comp,/sil)
      prof=total(total(prof[*,0,*,*],3),3)/float(info.sz)/float(info.sz)
      int_range=fix(round(info.intrange/info.wavesamp))
      wsteps=wavesind
      int_range=wsteps < int_range
      wrange=where(wsteps+int_range gt info.wavepoint-1)
      if (min(wrange) ne -1) then int_range[wrange]=info.wavepoint-1-wsteps[wrange]
      device,filename=dir_rep+'psf_plot_poly.eps',/encapsulated
      for oo=0,n_elements(wavesind)-1 do begin 
         otf = fltarr(info.sz,info.sz)
         for i=-int_range[oo],int_range[oo] do begin 
            otf = otf+reform(otfs[*,*,wavesind[oo]+i,oo])*prof[wavesind[oo]+i]
         endfor 
         ppsf[*,*,oo]=float(shift(fft(otf,1),info.sz/2,info.sz/2))
         plot,xpsf,ppsf(*,info.sz/2,oo)/max(ppsf[*,*,oo]),psy=-2,xr=[-2,2],yr=[0,1],xtit='arcsec',ytit='PSF!Dnorm!N',tit='Polychromatic PSFs',lines=oo,/noerase
      endfor
      device,/close
      printf,unit,'\includegraphics[width=10cm]{'+dir_rep+'psf_plot_poly.eps}'



   endelse 


endif ;info.routines(4)


; ------------------------------------------------------------------------------
; Pupil Apodization
; ------------------------------------------------------------------------------

;if (info.routines(5) eq 1) then begin
;   printf,unit,'\section*{SOPHISM Pupil Apodization Report}'

 
;endif ;info.routines(5)


; ------------------------------------------------------------------------------
; Detectors
; ------------------------------------------------------------------------------

;if (info.routines(6) eq 1) then begin
if (info.repomodules(6) eq 1) then begin
   printf,unit,'\newpage'
   printf,unit,'\section*{Focus Plane Assembly}'
   
   printf,unit,'FPA Shutter type: ',info.fpashutt,' \\'
   
   printf,unit,'Quantum efficiency: ',string(format='(F5.2)',info.fpaquan),' \\'
   printf,unit,'Full well: ',string(format='(F6.2)',info.fpawell*100.),' \% \\'
   
   printf,unit,'Photon Noise: ',info.fpaphotnoi,' \\'
   if (info.fpaphotnoi eq 'Load') then printf,unit,'\indent ',info.fpaphotnoiname,' \\'

;   printf,unit,'\subsection*{Readout Noise}'
   printf,unit,'Readout Noise: '
;   printf,unit,'Noise Type: ',info.noitype,' \\'
   printf,unit,'\indent Noise sigma: '+string(format='(F7.3)',info.noisigma)+' \\'

   printf,unit,'Darks: ',info.fpadarksel,' \\'
   printf,unit,'Gain Table: ',info.fpagainsel,' \\'
   if (info.fpagainsel eq 'Load') then printf,unit,'\indent ',info.fpagainname,' \\'
   if (info.fpagainsel eq 'Generate') then printf,unit,'\indent Gain Table rms: ',string(format='(F6.2)',info.fpagain),' \\'
   printf,unit,'Flats: ',info.fpaflatsel,' \\'
   if (info.fpaflatsel eq 'Load') then printf,unit,'\indent ',info.fpaflatname,' \\'
   if (info.fpaflatsel eq 'Generate') then begin
      if (max(info.fpaflatmore) gt 0) then begin
         printf,unit,'\indent for: '
         if (info.fpaflatmore(0) eq 1) then printf,unit,'each wavelength '
         if (info.fpaflatmore(1) eq 1) then printf,unit,'each polarization '
         printf,unit,' \\'
      endif 
      printf,unit,'\indent Max-min intensity variation: ',info.fpaflrang,' \\'
   endif 
   if (info.fpafring eq 1) then begin
      printf,unit,'Fringes: ',info.fpafrinsel,' \\'
      if (max(info.fpafringtime) gt 0) then begin
         printf,unit,'\indent variable in: '
         if (info.fpafringtime(0) eq 1) then printf,unit,'time '
         if (info.fpafringtime(1) eq 1) then pritnf,unit,'lambda '
         if (info.fpafringtime(2) eq 1) then printf,unit,'polarization '
         printf,unit,' \\'
      endif 
      printf,unit,'\indent Amplitude (\%): ',string(format='(F5.2)',info.fpafringamp),' \\'
      printf,unit,'\indent Width: ',string(format='(F5.2)',info.fpafringwidth),' \\'
      printf,unit,'\indent Direction: ',string(format='(F5.2)',info.fpafringdir),' \\'
   endif 

;dark and gain table
   if (info.dualbeam eq 0) then begin
      restore,info.saves+info.files(6)+'.sav'
      restore,info.saves+info.files(6)+'_flat.sav'
      !p.multi=[0,3,1]
      device,filename=dir_rep+'darkgain.eps',bit=64,/color,/encapsulated
      if (max(dc_im) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Dark' else tvframe,dc_im[0,*,*],/bar,/aspect,btit='Dark'
      if (max(gain_table) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Gain Table' else tvframe,gain_table,/bar,/aspect,btit='Gain Table'
      if (max(flat[0,0,*,*]*fring[0,0,0,*,*]) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Flat and fringes' else tvframe,reform(flat[0,0,*,*]*fring[0,0,0,*,*]),/bar,/aspect,btit='Flat and fringes'
      device,/close
      capt='detector'
   endif else begin
      restore,info.saves+info.files(6)+'.sav'
      restore,info.saves+info.files(6)+'_flat.sav'
      !p.multi=[0,3,2]
;stop
      device,filename=dir_rep+'darkgain.eps',bit=64,/color,/encapsulated
      if (max(dc_im) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Dark' else tvframe,dc_im[0,*,*],/bar,/aspect,btit='Dark'
      if (max(gain_table) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Gain Table' else tvframe,gain_table,/bar,/aspect,btit='Gain Table'
      if (max(flat[0,0,*,*]*fring[0,0,*,*]) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Flat and fringes' else tvframe,reform(flat[0,0,*,*]*fring[0,0,0,*,*]),/bar,/aspect,btit='Flat and fringes'
      restore,info.saves+info.files(6)+'_2.sav'
      restore,info.saves+info.files(6)+'_flat_2.sav'
      if (max(dc_im) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Dark' else tvframe,dc_im[0,*,*],/bar,/aspect,btit='Dark'
      if (max(gain_table) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Gain Table' else tvframe,gain_table,/bar,/aspect,btit='Gain Table'
      if (max(flat[0,0,*,*]*fring[0,0,*,*]) eq 0) then tvframe,fltarr(info.sz,info.sz),/bar,/aspect,btit='Flat and fringes' else tvframe,reform(flat[0,0,*,*]*fring[0,0,0,*,*]),/bar,/aspect,btit='Flat and fringes'
      device,/close
      capt='detector 1 (up) and detector 2 (down)'
   endelse 
   !p.multi=0

   ima=readfits(info.saves+info.files(6)+'_0.fits*',head,/comp,/sil)
   device,filename=dir_rep+'detect_ima.eps',bit=64,/color,/encapsulated
   tvframe,reform(ima[0,0,*,*])/mean(ima[0,0,*,*]),/asp,/bar,btit='I!DC!N'
   device,/close

   printf,unit,'\begin{figure}[!h]'
   printf,unit,'\begin{center}'
   printf,unit,'\includegraphics[width=12cm]{'+dir_rep+'darkgain.eps}\\'
   printf,unit,'\caption{Dark and gain table from '+capt+'.}'
   printf,unit,'\end{center}'
   printf,unit,'\end{figure}'
   printf,unit,'\begin{figure}[!h]'
   printf,unit,'\begin{center}'
   printf,unit,'\includegraphics[width=8cm]{'+dir_rep+'detect_ima.eps}'
   printf,unit,'\caption{Image example at detector.}'
   printf,unit,'\end{center}'
   printf,unit,'\end{figure}'

endif ;info.routines(6)


; ------------------------------------------------------------------------------
; Demodulation
; ------------------------------------------------------------------------------

;if (info.routines(8) eq 1) then begin
if (info.repomodules(8) eq 1) then begin

   printf,unit,'\newpage'
   printf,unit,'\section*{Demodulated images}'

;   nperiod=info.modnbuff*info.nacc
;   nmcycl=info.cycrep
;   ntim=nperiod*nmcycl*info.etalpos
;   nf=fix(ntim)
;   ima=readfits(info.saves+info.files(8)+'_'+strtrim(nf,2)+'.fits*',head,/comp,/sil)
   ima=readfits(info.saves+info.files(8)+'_0.fits*',head,/comp,/sil)

   !p.multi=[0,info.modnbuff/info.cycrep,info.etalpos]
   !x.omargin=[0,80]
   device,filename=dir_rep+'demod_ima.eps',bit=64,/color,/encapsulated
;no normalization. Use contind variable or make some supposition?
   for hh=0,info.etalpos/info.obsset-1 do begin
;   for hh=0,info.etalpos-1 do begin
      for jj=0,info.modnbuff/info.cycrep-1 do tvframe,ima[hh,jj,*,*],/asp,/bar,charsize=0.8
   endfor
   device,/close
   printf,unit,'\begin{center}'
   printf,unit,'\includegraphics[width=22cm]{'+dir_rep+'demod_ima.eps} \\'
   printf,unit,'\end{center}'

   avfil=findfile(info.saves+info.files(8)+'_av_0.fits*',count=nav)
   if (nav gt 0) then begin
      printf,unit,'\newpage'
      device,filename=dir_rep+'demod_ima_av.eps',bit=64,/color,/encapsulated
      ima=readfits(info.saves+info.files(8)+'_av_0.fits*',head,/comp,/sil)
      for hh=0,info.etalpos-1 do begin
         for jj=0,info.modnbuff/info.cycrep-1 do tvframe,ima[hh,jj,*,*],/asp,/bar,charsize=0.8
      endfor 
      device,/close
      printf,unit,'Data demodulated with averaged demodulation matrix: \\'
      printf,unit,'\begin{center}'
      printf,unit,'\includegraphics[width=22cm]{'+dir_rep+'demod_ima_av.eps}'
      printf,unit,'\end{center}'

   endif 

   if (max(info.reterrors) eq 1) then begin
      printf,unit,'\newpage'
      device,filename=dir_rep+'demod_ima_teo.eps',bit=64,/color,/encapsulated
      ima=readfits(info.saves+info.files(8)+'_teo_0.fits*',head,/comp,/sil)
      for hh=0,info.etalpos-1 do begin
         for jj=0,info.modnbuff/info.cycrep-1 do tvframe,ima[hh,jj,*,*],/asp,/bar,charsize=0.8
      endfor
      device,/close
      printf,unit,'Data demodulated with theoretical (wrong) demodulation matrix: \\'
      printf,unit,'\begin{center}'
      printf,unit,'\includegraphics[width=22cm]{'+dir_rep+'demod_ima_teo.eps}'
      printf,unit,'\end{center}'
   endif 

   !p.multi=0
   !x.omargin=0

endif ;demod--info.routines(8)

; ------------------------------------------------------------------------------
; Inversion
; ------------------------------------------------------------------------------

;if (info.routines(9) eq 1) then begin
if (info.repomodules(9) eq 1) then begin

   printf,unit,'\newpage'
   printf,unit,'\section*{Inversion results}'

   printf,unit,'Maximum number of iterations: '+strtrim(info.iter,2)+' \\'
   printf,unit,'Noise threshold factor: '+strtrim(info.inversigma,2)+' \\'

;   nperiod=info.modnbuff*info.nacc
;   nmcycl=info.cycrep
;   ntim=nperiod*nmcycl*info.etalpos
;   nf=fix(ntim)
;   restore,info.saves+info.files(9)+'_'+strtrim(nf,2)+'MILOS_inversion.sav'
;we choose the first of the observation sets (or the only one)
   ima=readfits(info.saves+info.files(9)+'_0.fits*',head,/comp,/sil)

   !p.multi=[0,2,2]

   device,filename=dir_rep+'invers_ima.eps',bit=64,/color,/encapsulated
   tvframe,ima(1,0,*,*),/asp,/bar,btit='v!Dlos!N (km/s)'
   tvframe,ima(2,0,*,*),/asp,/bar,btit='B (G)'
   tvframe,ima(3,0,*,*),/asp,/bar,btit='gamma (deg)'
   tvframe,ima(4,0,*,*),/asp,/bar,btit='azimuth (deg)'
   device,/close

   printf,unit,'\begin{center}'
   printf,unit,'\includegraphics[width=15cm]{'+dir_rep+'invers_ima.eps}'
   printf,unit,'\end{center}'

   avfil=findfile(info.saves+info.files(9)+'_av_0.fits*',count=nav)
   if (nav gt 0) then begin
      printf,unit,'\newpage'
;      restore,info.saves+info.files(9)+'_teo_'+strtrim(nf,2)+'MILOS_inversion.sav'
      printf,unit,'Inversion results from averaged demodulation: '

      ima=readfits(info.saves+info.files(9)+'_av_0.fits*',head,/comp,/sil)

      device,filename=dir_rep+'invers_ima_av.eps',bit=64,/color,/encapsulated
      tvframe,ima(1,0,*,*),/asp,/bar,btit='v!Dlos!N (km/s)'
      tvframe,ima(2,0,*,*),/asp,/bar,btit='B (G)'
      tvframe,ima(3,0,*,*),/asp,/bar,btit='gamma (deg)'
      tvframe,ima(4,0,*,*),/asp,/bar,btit='azimuth (deg)'
      device,/close
      printf,unit,'\begin{center}'
      printf,unit,'\includegraphics[width=15cm]{'+dir_rep+'invers_ima_av.eps}'
      printf,unit,'\end{center}'
   endif 

   if (max(info.reterrors) eq 1) then begin
      printf,unit,'\newpage'
;      restore,info.saves+info.files(9)+'_teo_'+strtrim(nf,2)+'MILOS_inversion.sav'
      printf,unit,'Inversion results from theoretical random-error-free modulation: '

      ima=readfits(info.saves+info.files(9)+'_teo_0.fits*',head,/comp,/sil)

      device,filename=dir_rep+'invers_ima_teo.eps',bit=64,/color,/encapsulated
      tvframe,ima(1,0,*,*),/asp,/bar,btit='v!Dlos!N (km/s)'
      tvframe,ima(2,0,*,*),/asp,/bar,btit='B (G)'
      tvframe,ima(3,0,*,*),/asp,/bar,btit='gamma (deg)'
      tvframe,ima(4,0,*,*),/asp,/bar,btit='azimuth (deg)'
      device,/close
      printf,unit,'\begin{center}'
      printf,unit,'\includegraphics[width=15cm]{'+dir_rep+'invers_ima_teo.eps}'
      printf,unit,'\end{center}'
   endif 

endif ;inversion--info.routines(9)

; ------------------------------------------------------------------------------
; Compression
; ------------------------------------------------------------------------------

if (info.repomodules(10) eq 1) then begin

   printf,unit,'\newpage'
   printf,unit,'\section*{Data Compression}'

   printf,unit,'Compression Ratio: '
   printf,unit,'Bits per pixel: '

   ima=readfits(info.saves+info.files(10)+'_0.fits*',head,/comp,/sil)
   headmod=sxpar(head,'HISTORY')
   invi=where(headmod eq 'Inversion Module')

   !p.multi=[0,2,2]

   device,filename=dir_rep+'compress_ima.eps',bit=64,/color,/encapsulated
   if (max(invi) eq -1) then begin
;demodulated data, not normalizing by continuum. Show the first
;spectral position, whatever that is
      tvframe,ima[0,0,*,*],/asp,/bar,btit='I'
      tvframe,ima[0,1,*,*],/asp,/bar,btit='Q'
      tvframe,ima[0,2,*,*],/asp,/bar,btit='U'
      tvframe,ima[0,3,*,*],/asp,/bar,btit='V'
   endif else begin
      tvframe,ima[1,0,*,*],/asp,/bar,btit='v!Dlos!N (km/s)'
      tvframe,ima[2,0,*,*],/asp,/bar,btit='B (G)'
      tvframe,ima[3,0,*,*],/asp,/bar,btit='gamma (deg)'
      tvframe,ima[4,0,*,*],/asp,/bar,btit='azimuth (deg)'
   endelse 
   device,/close

   printf,unit,'\begin{center}'
   printf,unit,'\includegraphics[width=15cm]{'+dir_rep+'compress_ima.eps}'
   printf,unit,'\end{center}'


endif ;compression


; ------------------------------------------------------------------------------
; Finalize
; ------------------------------------------------------------------------------

final:

set_plot, 'X'

printf,unit,'\end{document}'

free_lun,unit

spawn,'latex '+dir_rep+fil_rep+'.tex'
spawn,'dvipdf '+fil_rep+'.dvi '+fil_rep+'.pdf'
;spawn,'rm '+strmid(info.freport,barr+1,dott-barr-1)+'.dvi  '+strmid(info.freport,barr+1,dott-barr-1)+'.aux  '+strmid(info.freport,barr+1,dott-barr-1)+'.log  '
spawn,'rm '+fil_rep+'.dvi '+fil_rep+'.log '+fil_rep+'.aux'
spawn,'mv '+fil_rep+'.pdf '+dir_rep

; ------------------------------------------------------------------------------

end

