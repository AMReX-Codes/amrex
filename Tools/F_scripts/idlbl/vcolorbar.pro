; NAME:
;       vcolorbar
;
; PURPOSE:
;       This procedure produces a colorbar legend for a colormap plot
;
; CATEGORY:
;       Plotting, Contour
;
; CALLING SEQUENCE:
;       colorbar, CTR_LEVELS, COLORS, POS
;
; INPUTS:
;       CTR_LEVELS: The contour levels to make the legend with 
;
;       COLORS:     The colors associated with the contour levels
;
;       POS:        The position of the colorbar, specified in normal
;                   coordinates, POS = [x1, y1, x2, y2], where (x1,y1)
;                   are the coordinates of the lower left corner, and
;                   (x2,y2) are the coordinates of the upper right
;                   corner. 
;
; EXAMPLE:
;         Generate some data and make a contour plot and legend
;
;         x = findgen(20) # replicate(1,30)
;         y = replicate(1,20) # findgen(30)
;         z = sin(x^2 + y^2)
;
;         ctr_levels = [-1, -.5, 0, .5, 1]
;         
;         tek_color
;         colors = [2, 3, 4, 5, 6]
;
;         contour, z, x, y, /fill, levels = ctr_levels, c_color = colors, $   
;            pos = [.15, .30, .95, .95]     
;         
;         colorbar, ctr_levels, colors, [.15, .16, .95, .21]
;
; MODIFICATION HISTORY:
;       MZ -- 9-4-98
;-

PRO vcolorbar, min, max, pos, COLORMIN = colormin, COLORMAX = colormax, $
              LABEL_STEP = label_step, CHARSIZE = charsize

forward_function sci_notat, number


if n_elements(charsize) EQ 0 then charsize = 1
if n_elements(colormin) EQ 0 then colormin = 0
if n_elements(colormax) EQ 0 then colormax = 255

; produce a decent sized array containing the contour levels
ctr_levels = indgen(256)

bar = ctr_levels # replicate(1b,10)
bar = rotate(bar,1)
npts = 100
slope = (max - min)/(npts)

ctr_bar = indgen(npts)*slope # replicate(1b,10)
ctr_bar = rotate(ctr_bar,1)

; we need to tell IDL not to erase the page when we make the colorbar
!p.multi = [!p.multi[0] + 1]

; use axis style 4 to hide the axis and 1 to force the range
contour, ctr_bar,  position = pos,  $
  yrange = [min,max], xrange = [0,9], $
  xstyle = 5, ystyle = 5

; find the extent of this legend so we can label it
lower = convert_coord([0,min], /data, /to_normal) 
upper = convert_coord([9,max], /data, /to_normal)

ixsize = upper[0] - lower[0] + 1
iysize = upper[1] - lower[1] + 1

if colormax GT colormin then begin
    bar = bytscl(bar, top = colormax - colormin) + colormin
endif else begin
    bar = colormin - bytscl(bar, top = colormin - colormax)
endelse

tvimage, bar, pos = [lower[0], lower[1], upper[0], upper[1]], /overplot

; restore the axes
axis, lower[0],lower[1], yaxis = 0, /normal, ytickformat = 'nolabel', $
  ystyle = 1
axis, upper[0],lower[0], yaxis = 1, /normal, $
  ystyle = 1

end







