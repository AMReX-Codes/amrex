basec12 = "narrow_yc12_"
basetemp = "narrow_temp_"
basedyc12dt = "narrow_dyc12dt_"

extension = ".raw"

startfile = 0
endfile = 15400
skip = 50

width = 96.0
height = 3072.0

!p.multi = [0,1,3]

!p.charsize = 2.0

window, XSIZE=700, YSIZE=800

for i = startfile, endfile, skip do begin
    inum = string(i, format='(i5.5)')

    filename = basec12 + inum + extension
    yc12 = rawread(filename, width, height, XCOORDS=x, YCOORDS=y, TIME=time)

    filename = basetemp + inum + extension
    temp = rawread(filename, width, height, XCOORDS=x, YCOORDS=y, TIME=time)

    filename = basedyc12dt + inum + extension
    dyc12dt = rawread(filename, width, height, XCOORDS=x, YCOORDS=y, TIME=time)

    avg_yc12 = total(yc12,1)/(size(yc12))[1]
    avg_temp = total(temp,1)/(size(temp))[1]
    avg_dyc12dt = total(dyc12dt,1)/(size(dyc12dt))[1]
    
    plot, y, avg_yc12, xrange=[800,2400], xstyle=1, yrange=[0,.6], title = "Y(C12)", xtitle = "height [cm]", thick=2
    oplot, y, max(yc12, dimension=1), linestyle=1
    oplot, y, min(yc12, dimension=1), linestyle=1

    plot, y, avg_temp, xrange=[800,2400], xstyle=1, yrange=[0,2.5e9], title="temperature", xtitle = "height [cm]",thick=2
    oplot, y, max(temp, dimension=1), linestyle=1
    oplot, y, min(temp, dimension=1), linestyle=1

    plot, y, avg_dyc12dt, xrange=[800,2400], xstyle=1, yrange=[-60.0,0], title = "dY(C12)/dt", xtitle = "height [cm]", thick=2
    oplot, y, min(dyc12dt, dimension=1), linestyle=1
    oplot, y, max(dyc12dt, dimension=1), linestyle=1

    print, min(avg_dyc12dt)

    xyouts, 0.65, 0.95, "time = " + string(time), /normal, charsize=1.5

    color_bitmap, "frame_" + inum + '.png'

endfor

end
