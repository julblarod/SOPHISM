FUNCTION ODD,A
;+
; NAME:
;	ODD
;
; PURPOSE:
;	Search for odd numbers in A.
;
; CALLING SEQUENCE:
;	Result = ODD(A)
;
; INPUT:
;	A = array of any size; type must be byte, integer or longword.
;
; OUTPUT:
;	Result = a byte array of same size as A, filled with 1 and 0.
;		 0 is even and 1 means odd.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	Input should be of byte, integer or longword integer type.
;
; PROCEDURE:
;	Straightforward. Output = 1 is odd; 0 means even.
;
; MODIFICATION HISTORY:
;	Written by Roberto Luis Molowny Horas, 1991.
;-
ON_ERROR,2

	s = SIZE(a)
	array_type = s(s(0)+1)
	IF array_type NE 1 AND array_type NE 2 AND array_type NE 3 THEN $
		MESSAGE,'Input must be byte, integer or longword integer'

	dumb = a / 2.
	RETURN,LONG(dumb) NE dumb	;Easy algorithm

END