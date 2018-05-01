
;
	function DEVICE,AMR=amr,PHA=pha,ANG=ang
;
; INPUT:
;	amr	amplitude ratio, (P_{ort}/P_{par} where P_{ort} is the
;		transmission in the extinction direction and P_{par} is
;		the transmission in the linear polarization direction.
;
;	phase	phase difference introduced by the device (phase difference
;		in the perpendicular directrion minus phase difference in the 
;		parallel direction)
;	
;	ang   orientation of the device with respect to x-axis
;
; PURPOSE:
;	It computes the Mueller matrix of a general device with an amplitude 
;	ratio, a phase difference and an arbitrary orientation
;
;
;
; CALLING SEQUENCE: 
;	example: IDL>t=device(amr=0,pha=0,ang=0) lin. polarizer at 0 degrees
;	example: IDL>t=device(amr=0,pha=0,ang=90)lin. polarizer at 90 degrees
;	example: IDL>t=device(amr=0,pha=0,ang=45)lin. polarizer at 45 degrees
;	example: IDL>t=device(amr=0,pha=0,ang=135)lin. polarizer at 135 degrees
;	example: IDL>t=device(amr=1,pha=90,ang=0) lambda/4 at 0 degrees 
;	example: IDL>t=device(amr=1,pha=180,ang=0) lambda/2 at 0 degrees 
;
; CATEGORY:
;	Instrumental Polarization
; 
; MODIFICATION HISTORY:
;	November, 19, 1992	Valentin Martinez Pillet


;	Previous computations

	if (amr lt 0 or amr gt 1) then stop

	ang=double(ang)*double(!pi)/180.d0
	pha=double(pha)*double(!pi)/180.d0
	cpha=cos(pha)
	spha=sin(pha)
	cang=cos(2.d0*ang)
	sang=sin(2.d0*ang)
	pmas=0.5d0*(1+amr*amr)
	pmin=0.5d0*(1-amr*amr)
;	
;	Computing the Mueller matrix (IDL matrix are transpose of algebraic
;	matrix)

	t=dblarr(4,4)
	t(0,0)=pmas
	t(1,0)=cang*pmin
	t(2,0)=sang*pmin
	t(0,1)=cang*pmin
	t(1,1)=cang*cang*pmas+sang*sang*amr*cpha
	t(2,1)=cang*sang*pmas-cang*sang*amr*cpha
	t(3,1)=-sang*amr*spha
	t(0,2)=sang*pmin
	t(1,2)=cang*sang*pmas-cang*sang*amr*cpha
	t(2,2)=sang*sang*pmas+cang*cang*amr*cpha
	t(3,2)=cang*amr*spha
	t(1,3)=sang*amr*spha
	t(2,3)=-cang*amr*spha
	t(3,3)=amr*cpha
	t=t/t(0,0)
	z=where(abs(t) le 1.d-6)
	if(z(0) ne -1) then t(z)=0.d0
;	print,t
;	End
	ang=double(ang)*180.d0/double(!pi)
	pha=double(pha)*180.d0/double(!pi)

	return,t
	end
















