base = "narrow"
extension = ".C21.raw"

startfile = 0
endfile = 6900
skip = 50

width = 96.0
height = 3072.0

for i = startfile, endfile, skip do begin
    inum = string(i, format='(i4.4)')
    filename = base + '_' + inum + extension

    data = rawread(filename, width, height, XCOORDS=x, YCOORDS=y)

    length = flamelength(data,x,y,0.25)

    print, i, length
endfor

end

