function NOLABEL, axis, index, value
; a hack to get no axis labels when using the axis command, use
; axis, XTICKFORMAT='NOLABEL'
;
; function nolabel returns an empty string for no tick annotation
; 
; sometimes using this routine in this manner gives an underflow warning
return, ''
end
