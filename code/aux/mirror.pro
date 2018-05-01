;
;
	function MIRROR,angle,WAVE=wave,RI=ri,AC=ac
;
; INPUT:
;	angle	Incidence angle with respect to the mirror normal (degrees)
;
;       WAVE	wavelength
;
;	RI	refractive index
;
;	AC	absorption coefficient
;
; PURPOSE:
;	It computes the Mueller matrix of a mirror reflection.
;	The equations are given in Capitani et al SP 120, 173
;
;	If WAVE is given, then RI and AC are calculated whithin the 
;	function and RI and AC must not have assigned value
;
; CALLING SEQUENCE: 
;	example: IDL>t=mirror(45,ri=1.13,ac=6.39)
;
; CATEGORY:
;	Instrumental Polarization
; 
; MODIFICATION HISTORY:
;	October, 23, 1992	Valentin Martinez Pillet



	if(n_elements(wave) eq 0) then begin
		
		if(n_elements(ri) eq 0 and n_elements(ac) eq 0) then stop
		angle=double(angle)
		ang=angle*double(!pi)/180.d0
		ri=double(ri)
		ac=double(ac)

	endif else begin

		stop

;	this choice is still not implemented. It will use the wavelength
;	dependence of the refrective index of aluminium.

	endelse

;	Computing f and g

	v1=ri^2-ac^2-(sin(ang)^2)	
	v2=ac^2-ri^2+(sin(ang)^2)	
	v3=4.d0*ri^2*ac^2
	
	f=sqrt(0.5d0*(v1+sqrt(v1^2+v3)))
	g=sqrt(0.5d0*(v2+sqrt(v1^2+v3)))

;	Computing x and tau

	v4=sin(ang)*tan(ang)
	v5=f^2+g^2-2.d0*f*v4+v4^2
	v6=f^2+g^2+2.d0*f*v4+v4^2
	v7=v4^2-f^2-g^2
	x=sqrt(v5/v6)
	if (v4 eq 0.d0) then begin 
		tau=double(!pi)
	endif else begin
		tau=atan((2.d0*g*v4)/v7)
		if(tau lt 0.d0) then tau=tau+double(!pi)
		if(angle gt 89.9999) then tau=0.d0
	endelse
;	print,x,tau*180.d0/(!pi)
;	Computing the Mueller matrix (IDL matrix are transposed of algebraic
;	matrix!!!!!)
	t=dblarr(4,4)
	t(0,0)=x^2+1
	t(1,1)=x^2+1
	t(1,0)=x^2-1
	t(0,1)=x^2-1
	t(2,2)=2.d0*x*cos(tau)
	t(3,3)=2.d0*x*cos(tau)
	t(2,3)=-2.d0*x*sin(tau)
	t(3,2)=2.d0*x*sin(tau)
	t=t/t(0,0)


;	End

	return,t
	end
