FUNCTION Avg, data
;+
; NAME:
;       AVG
; PURPOSE:
;       Compute the average of the data
; CALLING SEQUENCE:
;       RESULT = AVG ( DATA )
; INPUTS:
;       DATA : Array. No limitations
; OUTPUTS:
;       RESULT: Average of data
; PROCEDURE:
; 	Return total/n_elements
; MODIFICATION HISTORY:
;	10-Jan-1990  P.Suetterlin, KIS
;-

return,total(data)/n_elements(data)

end
