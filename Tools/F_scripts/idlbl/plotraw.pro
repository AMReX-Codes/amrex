pro plotraw, filename, xrange=xrange, yrange=yrange, $
             postscript = postscript, log=log, colorbar=colorbar, $
             xvelfile=xvelfile, yvelfile=yvelfile, plotvel = plotvel, $
             lamspeed=lamspeed, varrange=varrange, yaxis=yaxis, xaxis=xaxis, yshift=yshift, $
             forceportrait=forceportrait

set_plot, 'X'

print, filename

if n_elements(postscript) EQ 0 then postscript = 0
if n_elements(log) EQ 0 then log = 0
if n_elements(colorbar) EQ 0 then colorbar = 1
if n_elements(plotvel) EQ 0 then plotvel = 0
if n_elements(lamspeed) EQ 0 then lamspeed = 0
if n_elements(varrange) EQ 0 then begin
    have_varrange = 0
endif else begin
    have_varrange = 1
endelse
if n_elements(yaxis) EQ 0 then yaxis = 1
if n_elements(xaxis) EQ 0 then xaxis = 1
if n_elements(yshift) EQ 0 then yshift = 0.0
if n_elements(forceportrait) EQ 0 then forceportrait = 0

hsize = 500
vsize = 600

if NOT postscript then window, XSIZE=hsize, YSIZE=vsize


; load the colortable -- we use special colortables, the first 12
; colors are reserved and constant across the different tables
idlbl_dir = get_idlbl_path()

colormap = "Menacing"
;colormap = "Blue -> Red"

; get the colortable
iclrmap = color_index(colormap, MIN_VALUE=colorMin, MAX_VALUE=colorMax)
loadct, iclrmap, FILE = idlbl_dir + 'flash_colors.tbl', /SILENT

;temp = colorMin
;colorMin = colorMax
;colorMax = temp

!p.background = color('white')
!p.color      = color('black')
!p.charsize = 1.0
!p.charthick = 2.0

data = rawread(filename, XCOORDS=x, YCOORDS=y, $
               XMIN=xmin, YMIN=ymin, XMAX=xmax, YMAX=ymax)

y = y + yshift
ymin = ymin + yshift
ymax = ymax + yshift

if n_elements(xrange) EQ 0 then begin
    xrange = [xmin, xmax]
    imin = 0
    imax = (size(x))[1]-1

endif else begin
    imin = (where(xrange[0] LT x))[0]
    imax = (where(xrange[1] LT x))[0]
    xrange[0] = x[imin]
    xrange[1] = x[imax]
endelse
   
if n_elements(yrange) EQ 0 then begin
    yrange = [ymin, ymax]
    jmin = 0
    jmax = (size(y))[1]-1

endif else begin
    print, yrange
    print, where(yrange[1] GT y)
    jmin = (where(yrange[0] LT y))[0]
    jmax = (where(yrange[1] LT y))[0]
    print, jmin, jmax
    yrange[0] = y[jmin]
    yrange[1] = y[jmax]
endelse

print, imin, imax
print, jmin, jmax
minvar = min(data[imin:imax,jmin:jmax])
maxvar = max(data[imin:imax,jmin:jmax])
print, 'varrange : ', minvar, maxvar
print, have_varrange


dy = yrange[1] - yrange[0]
dx = xrange[1] - xrange[0]

dataAspectRatio = dy/dx

; compute the fraction of the total domain that is being plotted
xfrac = dx/(xmax - xmin)
yfrac = dy/(ymax - ymin)


; ---- initialize the device --------------------------------------------------

; choose orientation based on dataAspectRatio
if postscript then begin
    current_device = !d.name
    set_plot, 'PS'

    if dataAspectRatio GE 1. then begin
        iorient = 0
    endif else begin
        iorient = 1
    endelse

    print, 'forceportrait', forceportrait
    if forceportrait then iorient = 0

    case iorient of

; portrait orientation
        0: begin

            xsize = 5.5
            ysize = 7.0

            outfile = filename + '_br.eps'
            device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
              XOFF = 0.5, YOFF = 0.5, /INCH, /COLOR, $
              BITS_PER_PIXEL = 8, LANDSCAPE=0, /encapsul, language_level=2

        end

; landscape orientation
        1: begin

            xsize = 10.
            ysize = 7.5

            outfile = filename + '_br.eps'
            device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
              XOFF = 0.5, YOFF = 10.5, /INCH, /COLOR, $
              BITS_PER_PIXEL = 8, /LANDSCAPE, /encapsul, language_level=2
        end
    endcase

    deviceAspect = (ysize)/(xsize)


endif else begin
    deviceAspect = float(vsize)/float(hsize)
endelse

; compute the aspect ratio of the domain and window (paper) together
aspect = dataAspectRatio/deviceAspect

; ---- determine the bounds for the plot --------------------------------------

; set the normal coordinates of the portion of the display/page you
; wish to use -- leave a little margin so the plots don't run to the
; edge of the page
page_nx1 = 0.05
page_nx2 = 0.95

page_ny1 = 0.05
page_ny2 = 0.90

; compute the bounding normal coordinates for this subplot
; the current subplot occupies the rectangle with opposite
; corners (nx1, ny1) to (nx2, ny2).  (note, 0,0 is in the
; lower left corner in normal coords).

dpagex = page_nx2 - page_nx1
dpagey = page_ny2 - page_ny1

nx1 = (dpagex) * (0.15) + page_nx1
nx2 = (dpagex) * (1.0 - 0.05)  + page_nx1

ny1 = page_ny2 - (dpagey) * (1.0 - 0.15)
ny2 = page_ny2 - (dpagey) * (0.05)

; aspect compares the aspect ratio of the data to that of the device,
; aspect = (dy_data/dy_device)/(dx_data/dx_device).  If aspect > 1,
; it means that the y-coordinate is the one that is going to set the
; overall scaling.  If aspect < 1, then the x-coordinate will set the
; scaling.  Consider each case separately below.
cb_factor = 0.90

if (aspect GE 1.) then begin

; the y size of the data sets the scaling of the plot

; set the initial values of the y min and max normal coordinates.  We
; leave some room for the axis labels
    py1 = ny1
    py2 = ny2

; compute the x size, using the aspect ratio
    dpy = py2 - py1
    dpx = dpy/aspect < cb_factor 


; recompute dpy, in case dpx was adjusted
    dpy = aspect*dpx

; compute the plot coordinates
    px1 = nx1
    px2 = px1 + dpx

    py_center = 0.5*(ny1 + ny2)
    py1 = py_center - .5*dpy
    py2 = py_center + .5*dpy

; set the plot and legend bounding box -- the legend will be vertical
    plot_pos   = [px1, py1, px2, py2]

; the legend goes in the remaining space in the x direction -- figure
; out how much room is left there
    dx_legend = dpagex  +page_nx1 - px2
    lwidth = dx_legend/4 < 0.25*(px2 - px1)
    lcenter = 0.5* (dpagex + page_nx1 + px2)
    legend_pos = [lcenter-0.5*lwidth, py1, lcenter+0.5*lwidth, py2]

endif else begin

; the x size of the data sets the scaling of the plot


; set the initial x min and max normal coordiantes
    px1 = nx1
    px2 = nx2

    dpx = px2 - px1
    dpy = aspect*dpx < cb_factor 

; recompute dpx, in case dpy was adjusted
    dpx = dpy/aspect

; recompute the plot coordinates
    px_center = 0.5*(nx1 + nx2)
    px1 = px_center - .5*dpx
    px2 = px_center + .5*dpx

    py2 = ny2
    py1 = py2 - dpy

; set the plot and legend bounding box -- the legend will be horizontal
    plot_pos   = [px1, py1, px2, py2]

; the legend goes in the remaining space in the lower y direction --
; figure out how much room is left there
    dy_legend = py1 - (page_ny2 - dpagey)
    lheight = dy_legend/4 < 0.25*(py2 - py1)
    lcenter = 0.5*(py1 + (page_ny2 - dpagey))
    legend_pos = [px1, lcenter-0.5*lheight, px2, lcenter+0.5*lheight]
endelse

if have_varrange then begin
    minvar = varrange[0]
    maxvar = varrange[1]
endif

case log of 
    0: sdata = scale_color(data, VARMIN=minvar, VARMAX=maxvar, $
                           COLORMAP_MIN = colorMin, COLORMAP_MAX = colorMax)
    1: sdata = scale_color(alog10(data), VARMIN=alog10(minvar), VARMAX=alog10(maxvar), $
                           COLORMAP_MIN = colorMin, COLORMAP_MAX = colorMax)
endcase

plot, [xrange[0], xrange[1]], [yrange[0], yrange[1]], POS = plot_pos, $
  XSTYLE=5, YSTYLE=5

lower = convert_coord([xrange[0], yrange[0]], /DATA, /TO_NORMAL)
upper = convert_coord([xrange[1], yrange[1]], /DATA, /TO_NORMAL)

tvimage, sdata[imin:imax,jmin:jmax], position = [lower[0], lower[1], upper[0], upper[1]], $
  /OVERPLOT, /HALF_HALF

if plotvel then begin

    xpts = (size(x))[1]
    ypts = (size(y))[1]

    xt = x # replicate(1,ypts)
    yt = replicate(1,xpts) # y

    xvel = rawread(xvelfile)
    yvel = rawread(yvelfile) - lamspeed

    partvelvec, xvel, yvel, xt, yt, /over, color=color('verygray'), $
      typvel = 5.e4, minmag=0, maxmag=1.e10, xskip=16, yskip=16
    
endif

if (xaxis) then begin
    axis, xrange[0], yrange[0], xaxis = 0, xstyle = 1, xtitle = 'x (cm)', xthick=2.0
endif else begin
    axis, xrange[0], yrange[0], xaxis = 0, xstyle = 1, xtickformat = 'nolabel', xthick=2.0
endelse

if (yaxis) then begin
    axis, xrange[0], yrange[0], yaxis = 0, ystyle = 1, ytitle = 'r (cm)', ythick=2.0
endif else begin
    axis, xrange[0], yrange[0], yaxis = 0, ystyle = 1, ytickformat = 'nolabel',ythick=2.0
endelse

axis, xrange[0], yrange[1], xaxis = 1, xstyle = 1, xtickformat = 'nolabel', xthick=2.0
axis, xrange[1], yrange[0], yaxis = 1, ystyle = 1, ytickformat = 'nolabel', ythick=2.0

if colorbar then begin
    if aspect LT 1 then begin

        case log of
            0: colorbar2, float(minvar), float(maxvar), $
              legend_pos, $
              COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5
            
            1: colorbar2, alog10(minvar), alog10(maxvar), $
              legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, $
              CHARSIZE = 1.5
        endcase
        
    endif else begin
        
        case log of
            
            0: vcolorbar, float(minvar), float(maxvar), $
              legend_pos, $
              COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5
            
            1: vcolorbar, alog10(minvar), alog10(maxvar), $
              legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, $
              CHARSIZE = 1.5
            
        endcase
    endelse
endif

if postscript then device, /close

end
