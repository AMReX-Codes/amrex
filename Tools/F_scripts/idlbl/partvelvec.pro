; NAME:
;       partvelvec
;
; PURPOSE:
;       This procedure plots the velocity vectors of particles (at the
;       positions of the particles).
;
; CATEGORY:
;       Plotting, Two-dimensional.
;
; CALLING SEQUENCE:
;       PARTVELVEC, VELX, VELY, POSX, POSY [, X, Y]
;
; INPUTS:
;       VELX:  An array of any dimension, containing the x-components
;              of the particle velocities.
;       VELY:  An array of the same dimension as velx, containing the
;              y-components of the particle velocities.
;       POSX:  An array of the same dimension as velx, containing the
;              x-components of the particle positions.
;       POSY:  An array of the same dimension as velx, containing the
;              y-components of the particle positions.
;
; OPTIONAL INPUTS:
;       X:   Optional abcissae values. X must be a vector.
;       Y:   Optional ordinate values. Y must be a vector. If only X
;            is specified, then Y is taken equal to be equal to X.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       FRACTION:   The fraction of the vectors to plot. They are
;                   taken at random from the complete sample.    
;                   Default is FRACTION = 1.0, use all vectors
;
;       LENGTH:     The maximum vectorlength relative to the plot data
;                   window.   
;                   Default = 0.08
;
;       COLOR:      The color for the vectors, axes and titles.
;                   Default=!P.COLOR
;
;       OVER:       Plot over the previous plot
;
;       MINMAG:     The minimum magnitude vector to plot
;                   Default is plot all
;
;       MAXMAG:     The maximum magnitude vector to plot
;                   Default is 1.e30
;
;       TYPVEL:     The magnitude to scale all velocities to.  This is
;                   useful if you want the arrow length to mean the
;                   same thing across plots.  
;                   The default is the maximum velocity
;
;       XSKIP &     Plot every xskip or yskip vector in x or y respectively.
;       YSKIP:      This allows you to thin out the velocity field in
;                   the x and y directions independently. 
;      
;       LEGEND:     The x & y NORMAL coordinates to plot a velocity legend
;
;       OUTLINE:    Draw a white outline around the arrows --
;                   increases contrast
;              
;       Plot        All other keywords available to PLOT are also used
;       Keywords:   by this procedure.
;
; OUTPUTS:
;       This procedure plots the velocity vectors (VELX,VELY) at the
;       positions of the particles, (POSX,POSY). If X and Y are not
;       specified, then the size of the plot is such that all vectors
;       just fit within in the plot data window.
;
; SIDE EFFECTS:
;       Plotting on the current device is performed.
;
; EXAMPLE:
;       Generate some particle positions and velocities.
;
;         POSX=RANDOMU(seed,200)
;         POSY=RANDOMU(seed,200)
;         VELX=RANDOMU(seed,200)-0.5
;         VELY=RANDOMU(seed,200)-0.5
;
;       Plot the particle velocities.
;
;         PARTVELVEC, VELX, VELY, POSX, POSY
;
; MODIFICATION HISTORY:
;       Written by:  Joop Schaye (jschaye@astro.rug.nl), Sep 1996.
;
;       Modified:    Theo Brauers (th.brauers@fz-juelich.de) Oct. 1997
;                    use with maps, incl. debug
;             
;                    Michael Zingale (zingale@oddjob.uchicago.edu)
;                    Aug. 1998, added minimum and maximum cutoff for 
;                    velocity, legend plotting, clipping of velocity
;                    vectors outside the plot window, option of
;                    skipping in x and y directions, and scaling to a
;                    typical velocity instead of the max if desired.
;-

PRO partvelvec, orig_velx, orig_vely, orig_posx, orig_posy, x, y, $
                OVER=over, FRACTION=fraction, LENGTH=length, COLOR=color, $
                _EXTRA=extra, MINMAG = orig_minmag, MAXMAG = orig_maxmag, $ 
                XSKIP = xskip, YSKIP = yskip, TYPVEL = orig_typvel, $
                LEGEND = legend, LEGCLR = legclr, PARTICLES=particles,$
                TAG=tag, SYM_SIZE=sym_size, SHOWVECTORS=showvectors, $
                SHOWTAGS=showtags, OUTLINE=outline

;               debug, '1.10 T.B. 1997-OCT-20'

;---------------------------------------------
; Various settings, modify these to customize
;---------------------------------------------


forward_function sci_notat, number

c={customize, $
   length: 0.08, $     ; Maximum vector length relative to plot region. (*)
   lengtharrow: 0.3, $ ; Length of arrowhead legs relative to vectorlength.
   angle: 22.5 }       ; 1/2 times the angle between the arrowhead legs.

; (*) Not used if keyword LENGTH is present

; ---------------------------------------
; create working copies of the arguments
;----------------------------------------
posx   = orig_posx
posy   = orig_posy
velx   = orig_velx
vely   = orig_vely
typvel = orig_typvel

if (n_elements(orig_minmag) EQ 0) then orig_minmag = 0
if (n_elements(orig_maxmag) EQ 0) then orig_maxmag = 1.e30

if (n_elements(xskip) EQ 0) then xskip = 0
if (n_elements(yskip) EQ 0) then yskip = 0

if (n_elements(particles) EQ 0) then particles = 0
if (n_elements(showvectors) EQ 0) then showvectors = 0
if (n_elements(showtags) EQ 0) then showtags = 0

if (n_elements(sym_size) EQ 0) then sym_size=0

if (n_elements(outline) EQ 0) then outline = 0

minmag = orig_minmag
max_max = orig_maxmag

;---------------------
; Some error handling
;---------------------

on_error,2  ; Return to caller if an error occurs.

nparams=n_params()
IF nparams NE 4 THEN BEGIN
    IF (nparams NE 5 AND nparams NE 6) THEN BEGIN
        message,'Wrong number of parameters!',/continue
        message,'Syntax: PARTVELVEC, VELX, VELY, POSX, POSY [, X, Y]', $
          /noname,/noprefix
    ENDIF
    IF nparams EQ 5 THEN y=x
    sizex=size(x)
    sizey=size(y)
    IF (sizex(0) NE 1 OR sizey(0) NE 1) THEN $
      message,'X and Y must be vectors!'
ENDIF

sizevelx=size(velx)
sizevely=size(vely)
sizeposx=size(posx)
sizeposy=size(posy)

IF (total(sizevelx(0:sizevelx(0))-sizevely(0:sizevelx(0))) NE 0 $
    OR total(sizevelx(0:sizevelx(0))-sizeposx(0:sizevelx(0))) NE 0 $
    OR total(sizevelx(0:sizevelx(0))-sizeposy(0:sizevelx(0))) NE 0) THEN $
  message,'All arguments must have the same dimension and size!'

IF n_elements(fraction) GT 0 THEN $
  IF (fraction LT 0.0 OR fraction GT 1.0) THEN $
  message,'Fraction has to be between 0.0 and 1.0.'

; modified -- MZ
;---------------------------------------------------
; thin out the velocities if xskip or yskip are set
;---------------------------------------------------

if xskip LT 1 then xskip = 1
if yskip LT 1 then yskip = 1

print, 'skipping = ', xskip, yskip

if xskip GT 1 then begin
    xskip = fix(xskip)

    help, velx

    xindices = lindgen(n_elements(posx[*,0])/xskip)*xskip

    posx = temporary(posx[xindices,*])
    posy = temporary(posy[xindices,*])

    velx = temporary(velx[xindices,*])
    vely = temporary(vely[xindices,*])
    help, velx
endif

if yskip GT 1 then begin
    yskip = fix(yskip)
    
    yindices = lindgen(n_elements(posx[0,*])/yskip)*yskip

    posx = temporary(posx[*,yindices])
    posy = temporary(posy[*,yindices])

    velx = temporary(velx[*,yindices])
    vely = temporary(vely[*,yindices])
endif

;--------------
; Prepare plot
;--------------

nvecs=n_elements(velx)  ; Number of particles.
vel=sqrt(velx^2+vely^2)  ; Total velocity.

maxvel=max(vel)  ; Maximum velocity.
print, '<<< maxvel = ', maxvel, ' >>>'
if n_elements(typvel) EQ 0 then typvel = maxvel

; Compute maximum length of vectors.
IF n_elements(length) LE 0 THEN length=c.length
minposx=min(posx)
maxposx=max(posx)
minposy=min(posy)
maxposy=max(posy)
length=length*((maxposx-minposx) > (maxposy-minposy))


; Convert velocities.
vx=length*temporary(velx)/typvel
vy=length*temporary(vely)/typvel
vel=length*temporary(vel)/typvel

; modified -- MZ
; check to see if a minimum or maximum vector length was specified
IF n_elements(minmag) EQ 0 then begin
    minmag = -1
endif else begin
    minmag = length*minmag/typvel
endelse

IF n_elements(maxmag) EQ 0 then begin
    maxmag = 1.e30
endif else begin
    maxmag = length*maxmag/typvel
endelse


; Make sure no vectors extend beyond the plot data window.
x1=posx+vx  ; End of vector.
y1=posy+vy
IF (nparams EQ 4 and n_elements(over) EQ 0) THEN BEGIN
    minposx=min(x1)<minposx
    maxposx=max(x1)>maxposx
    minposy=min(y1)<minposy
    maxposy=max(y1)>maxposy
ENDIF


angle=c.angle*!dtor  ; Convert from degrees to radians.
sinangle=sin(angle)  ; Need these.
cosangle=cos(angle)


;-----------
; Plot axes
;-----------

IF n_elements(color) EQ 0 THEN color=!p.color

IF n_elements(over) EQ 0 THEN BEGIN
  IF nparams EQ 4 THEN $
    plot,[minposx,maxposx],[minposy,maxposy], $
     /nodata,/xstyle,/ystyle,COLOR=color,_EXTRA=extra $
  ELSE plot,x,y,/nodata,/xstyle,/ystyle,COLOR=color,_EXTRA=extra
ENDIF

;--------------
; Plot vectors
;--------------

IF n_elements(fraction) GT 0 THEN BEGIN
    IF fraction EQ 1.0 THEN GOTO,plotall
    nrgood=long(fraction*nvecs)  ; # of vectors to plot.
    IF nrgood EQ 0 THEN return
    ; Compute indices of vectors to plot. I use two lines to get more
    ; random "random numbers".
    good=long(randomu(seed,nrgood+1)*(nvecs-1.0))
    good=good(1:*)
    vx=temporary(vx(good))
    vy=temporary(vy(good))
    px=posx(good)  ; Can't use temporary if we wan't to keep the data.
    py=posy(good)
    x1=temporary(x1(good))
    y1=temporary(y1(good))
    nvecs=nrgood
ENDIF ELSE BEGIN
plotall:
    px=posx
    py=posy
ENDELSE

print, '**** in partvelvec, nvecs = ', nvecs
print, 'max x = ', max(posx)

FOR i=0l,nvecs-1l DO BEGIN  ; Loop over particles.

; Note that we cannot put the next three lines outside the loop,
; because we want the arrow size to be relative to the vector length.
    r=c.lengtharrow*vel(i)  ; Length of arrow head.
    rsin=r*sinangle
    rcos=r*cosangle

; Draw basis, arrow leg, same arrow leg, other arrow leg.
; One arrow leg is drawn twice, because we need to return to the end
; of the vector to draw the other leg.

; modified -- MZ

; define a user symbol that is a circle
    NPTS = 24
    tsym = findgen(NPTS)*2.*!pi/NPTS
    xsym = cos(tsym)
    ysym = sin(tsym)

    usersym, xsym, ysym, /fill
    if ((vel(i) GE minmag) and (vel(i) LE maxmag)) then begin

        if (particles EQ 0 OR (particles EQ 1 AND showvectors EQ 1)) then begin

            if (outline EQ 1) then begin
                plots,[px(i),x1(i),x1(i)-(vx(i)*rcos+vy(i)*rsin)/vel(i), $
                       x1(i),x1(i)-(vx(i)*rcos-vy(i)*rsin)/vel(i)], $
                  [py(i),y1(i),y1(i)-(vy(i)*rcos-vx(i)*rsin)/vel(i), $
                   y1(i),y1(i)-(vy(i)*rcos+vx(i)*rsin)/vel(i)], $
                  COLOR=255, thick=2, noclip = 0
            endif

            plots,[px(i),x1(i),x1(i)-(vx(i)*rcos+vy(i)*rsin)/vel(i), $
                   x1(i),x1(i)-(vx(i)*rcos-vy(i)*rsin)/vel(i)], $
              [py(i),y1(i),y1(i)-(vy(i)*rcos-vx(i)*rsin)/vel(i), $
               y1(i),y1(i)-(vy(i)*rcos+vx(i)*rsin)/vel(i)], $
              COLOR=color, noclip = 0

        endif

        if (particles EQ 1) then begin
            oplot, [px(i)], [py(i)], psym=8, symsize=sym_size
            if (showtags EQ 1) then $
              xyouts, px(i), py(i), ' '+ strtrim(string(tag(i)),2)
        endif
    endif

ENDFOR

if n_elements(legend) EQ 2 then begin

    if n_elements(legclr) EQ 0 then legclr = 0

    ptemp = convert_coord(legend[0], legend[1], /normal, /to_data)
    px = ptemp[0]
    py = ptemp[1]

    vx = length
    vy = 0
    
    vel = vx

    r = c.lengtharrow*vel
    rcos = r*cosangle
    rsin = r*sinangle

    x1 = px + vx
    y1 = py + vy

    plots,[px,x1,x1-(vx*rcos+vy*rsin)/vel,x1,x1-(vx*rcos-vy*rsin)/vel], $
      [py,y1,y1-(vy*rcos-vx*rsin)/vel, y1,y1-(vy*rcos+vx*rsin)/vel], $
      COLOR=legclr


    vx = 0
    vy = length
    
    vel = vy

    x1 = px + vx
    y1 = py + vy

    plots,[px,x1,x1-(vx*rcos+vy*rsin)/vel,x1,x1-(vx*rcos-vy*rsin)/vel], $
      [py,y1,y1-(vy*rcos-vx*rsin)/vel, y1,y1-(vy*rcos+vx*rsin)/vel], $
      COLOR = legclr


    xyouts, x1 + length, py + vy/2, sci_notat(typvel) + ' cm/s'
endif

END  ; End of procedure PARTVELVEC.







