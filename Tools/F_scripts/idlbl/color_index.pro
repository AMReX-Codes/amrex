; store the names of the colortables defined for in
; lbl_colors.tbl.  Given the string name of the color, return the
; index into the table

function color_index, name, GET_NAMES=colors, $
                      MIN_VALUE=colorMin, MAX_VALUE=colorMax

colors = ['Red Temp', $
          'Blue -> Red', $
          'Rainbow', $
          'Fire', $
          'Menacing', $
          'RT', $
          'Grayscale', $
          'vort']

index = (where(name EQ colors))[0]

; the first few colors in the table are reserved as standard
; colors.   Set the range of usable colors in the table.

case index of 
    0: begin
        colorMin = 255
        colorMax = 120
    end
    1: begin
        colorMin = 255
        colorMax = 20
    end
    6: begin
        colorMin = 255
        colorMax = 60
    end
    else: begin
        colorMin = 255
        colorMax = 12
    end
endcase

return, index

end

