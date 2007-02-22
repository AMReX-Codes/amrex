function flamelength, data, x, y, threshold

contour, data, x, y, level=[threshold], path_xy=xy, path_info=info, closed=0, /path_double
contour, data, x, y, level=[threshold]

; the coordinate information in xy is stored in normalized units [0:1]
; transform it to physical space
coords = convert_coord(xy, /normal, /data)

xcoords = reform(coords[0,*])
ycoords = reform(coords[1,*])

; info.type tells us if it is closed (1) or open (0)
; info.level tells us which contour level it is
; info.n tells us the number of xy pairs that make up the contour
; offset tells us where in the xy[] array the pairs begin

dist = 0.0

; loop over all the contours
for i = 0, (n_elements(info) - 1) do begin

; compute an array giving the relative indices of all the points that
; make up the current contour.  If we are closed, we need to add a 0
; to the end so that we wrap around
    if (info[i].type EQ 1) then begin
        s = [indgen(info[i].n), 0]
    endif else begin
        s = [indgen(info[i].n)]
    endelse

; the points that make up the current contour begin are just
; info[i].offset + s in the xy array.  Loop over all the points and
; compute the length.
    ldist = 0.0

    ind = info[i].offset + s

    for n = 1, info[i].n - 1 do begin
        ldist = ldist + sqrt((xcoords[ind[n]] - xcoords[ind[n-1]])^2 + $
                             (ycoords[ind[n]] - ycoords[ind[n-1]])^2)
    endfor

    dist = dist + ldist


    for n = 0, info[i].n-1 do begin
        plots, [xcoords[ind[n]], ycoords[ind[n]]], /continue
    endfor
        
endfor

return, dist

end
 
