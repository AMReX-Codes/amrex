
index = 100


file = "rt_1.5e7_3d_" + string(index, format="(i4.4)") + ".Y(C12).raw"
print, file

data = rawread3d(file,/single)

data = temporary(data[*,*,256:767])
help, data

;s = size(data)
;scale3, xrange=[0, s[1]], yrange=[0, s[2]], zrange=[0, s[3]]
;shade_volume, data, 0.1, v, p
;help, v
;help, p
;tv, polyshade(v,p,/t3d)

isosurface, data, 0.1, v, conn
help, v
help, conn

; we want the connectivity 1 based
conn = conn + 1

; print out the vertices
n = (size(v))[2]
m = (size(conn))[1]

openw,lun,'triangles.dat',/get_lun

printf, lun, '# of vertices = ', n
printf, lun, '# of triangles = ', m/4
for i = 0l, n-1 do begin
    printf, lun, format='(f,f,f)', v[*,i]
endfor

; now the connectivity
for i = 0l, m-1, 4 do begin
    printf, lun, format='(f,f,f)',conn[i+1], conn[i+2], conn[i+3]
endfor

close, lun


end
