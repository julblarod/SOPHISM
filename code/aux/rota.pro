
;
	function ROTA,angle
;
; INPUT:
;	angle   Rotation angle of the reference system
;
; PURPOSE:
;	It computes the Mueller matrix of a rotation of the reference
;	system. The sign convention is that of Capitani et al SP 120, 173
;
;
; CALLING SEQUENCE: 
;	example: IDL>t=rota(45)
;
; CATEGORY:
;	Instrumental Polarization
; 
; MODIFICATION HISTORY:
;	October, 23, 1992	Valentin Martinez Pillet


	angle=double(angle)
	ang=angle*double(!pi)/180.d0

;	Computing the Mueller matrix (IDL matrix are transpose of algebraic
;	matrix)

	t=dblarr(4,4)
	t(0,0)=1.d0
	t(1,1)=cos(2.d0*ang)
	t(2,2)=cos(2.d0*ang)
	t(1,2)=-sin(2.d0*ang)
	t(2,1)=sin(2.d0*ang)
	t(3,3)=1.d0
	

;	End

	return,t
	end
















