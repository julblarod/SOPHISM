function sign, var
;+
;
;	function:  sign
;
;	purpose:  return the sign of a number
;
;	author:  rob@ncar, 11/93
;
;==============================================================================
;
;	Check number of parameters.
;
if n_params() ne 1 then begin
	print
	print, "usage:  ret = sign(var)"
	print
	print, "	Return the sign of a number"
	print, "	(-1 if var lt 0; +1 if var ge 0).
	print
	print, "	Arguments"
	print, "		var	- number to return the sign of"
	print
	print
	print, "   ex:  v2 = v1 * sign(v0)
	print
	return, 0
endif
;-
;
;	Return I.
;
return, 2*(var ge 0) - 1
end
