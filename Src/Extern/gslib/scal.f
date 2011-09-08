        subroutine scal(xmin,xmax,ymin,ymax,xaxmin,xaxmax,yaxmin,
     +                  yaxmax,ilog,i45)
c-----------------------------------------------------------------------
c
c Draws a reasonable graph axes for a PostScript plot.  The appropriate
c labelling and tic mark interval are established.
c
c INPUT VARIABLES:
c       xmin   - the minimum of the x axis (labeled on axis)
c       xmax   - the maximum of the x axis (labeled on axis)
c       ymin   - the minimum of the y axis (labeled on axis)
c       ymax   - the maximum of the y axis (labeled on axis)
c       xaxmin - the minimum of the x axis (on PostScript window)
c       xaxmax - the maximum of the x axis (on PostScript window)
c       yaxmin - the minimum of the y axis (on PostScript window)
c       yaxmax - the maximum of the y axis (on PostScript window)
c       ilog   - scale option: 0 - both cartesian
c                              1 - semi log with x being log scale
c                              2 - semi log with y being log scale
c                              3 - log log with both being log scale
c       i45    - 45 degree line option: 0 - no, 1 - if axes are the same
c              
c
c
c-----------------------------------------------------------------------
      parameter (EPS=0.001)
      real       xloc(5),yloc(5)
      character  label*8,lfmt*8
c
c Common Block for Postscript Output Unit and Scaling:
c
      common /psdata/ psscl,pxmin,pxmax,pymin,pymax,wxmin,
     +                wxmax,wymin,wymax,lpsout
c
c Check to make sure that the scale can be plotted:
c
      if((xmax-xmin).le.0.0001.or.(ymax-ymin).le.0.0001) return
c
c Set up some of the parameters:
c
      tlng =       0.013 * ((xaxmax-xaxmin) + (yaxmax-yaxmin))
      tsht =       0.007 * ((xaxmax-xaxmin) + (yaxmax-yaxmin))
      psz  = 4.0 + 0.060 *  (xaxmax-xaxmin)
      pl1  = 0.6 + 0.005 *  (xaxmax-xaxmin)
      pl2  = 0.3 + 0.003 *  (xaxmax-xaxmin)
c
c Draw the axis:
c
      xloc(1) = xaxmin
      yloc(1) = yaxmax
      xloc(2) = xaxmin
      yloc(2) = yaxmin
      xloc(3) = xaxmax
      yloc(3) = yaxmin
      call psline(3,xloc,yloc,pl1,0)
c
c Show a 45 degree line?
c
      if(i45.eq.1) then
         if(abs(xmin-ymin).le.0.0001.and.abs(xmax-ymax).le.0.0001) then
            xloc(1) = xaxmin
            yloc(1) = yaxmin
            xloc(2) = xaxmax
            yloc(2) = yaxmax
            call psline(2,xloc,yloc,pl1,0)
         end if
      end if
c
c CONSTRUCT THE X AXIS:
c
c
c Log scale?
c
      if(ilog.eq.1.or.ilog.eq.3) then
c
c      The start, end, number of log(10) cycles, the tic mark start,
c      and the number of points in defining a tic mark:
c
            tminx   = alog10(xmin)
            tmaxx   = alog10(xmax)
            ncyc    = tmaxx - tminx
            cbas    = xmin/10
            yloc(1) = yaxmin
            num     = 2
c
c      Loop along the axis drawing the tic marks and labels:
c
            do icyc=1,ncyc+1
                  cbas = cbas * 10
                  do i=1,9
                  t1   = alog10(cbas*real(i))
                  xloc(1) = resc(tminx,tmaxx,xaxmin,xaxmax,t1)
                  xloc(2) = xloc(1)
                  if(i.eq.1) then
c
c            First point - long tic mark:
c
                        yloc(2) = yloc(1) - tlng
                        call psline(num,xloc,yloc,pl2,0)
                        yloc(2) = yloc(1) - 2.5*tlng
                        if(abs(t1+9.).le.EPS) label = '1.0e-9  '
                        if(abs(t1+8.).le.EPS) label = '1.0e-8  '
                        if(abs(t1+7.).le.EPS) label = '1.0e-7  '
                        if(abs(t1+6.).le.EPS) label = '1.0e-6  '
                        if(abs(t1+5.).le.EPS) label = '0.00001 '
                        if(abs(t1+4.).le.EPS) label = '0.0001  '
                        if(abs(t1+3.).le.EPS) label = '0.001   '
                        if(abs(t1+2.).le.EPS) label = '0.01    '
                        if(abs(t1+1.).le.EPS) label = '0.1     '
                        if(abs(t1)   .le.EPS) label = '1       '
                        if(abs(t1-1.).le.EPS) label = '10      '
                        if(abs(t1-2.).le.EPS) label = '100     '
                        if(abs(t1-3.).le.EPS) label = '1000    '
                        if(abs(t1-4.).le.EPS) label = '10000   '
                        if(abs(t1-5.).le.EPS) label = '100000  '
                        if(abs(t1-6.).le.EPS) label = '1.0e+6  '
                        if(abs(t1-7.).le.EPS) label = '1.0e+7  '
                        if(abs(t1-8.).le.EPS) label = '1.0e+8  '
                        if(abs(t1-9.).le.EPS) label = '1.0e+9  '
                        call pstext(xloc(1),yloc(2),8,label,
     +                                      psz,1,0.d0,1)
                  else
c
c            Not first point - short tic mark:
c
                        if(icyc.le.ncyc) then
                              yloc(2) = yloc(1) - tsht
                              call psline(num,xloc,yloc,pl2,0)
                        endif
                  endif
                    end do
              end do
      else
c
c Arithmetic Scale:
c
            do i=1,20
                  test = (xmax-xmin)/(10.0**(6-i))
                  if(test.gt.0.9) go to 1
            end do
 1          if(test.gt.3.0)                 zval = 1.0
            if(test.le.3.0.and.test.gt.2.0) zval = 0.5
            if(test.le.2.0.and.test.gt.1.2) zval = 0.4
            if(test.le.1.2)                 zval = 0.2
            nval = 5
            if(zval.eq.0.4.or.zval.eq.0.2)  nval = 4
            zval = zval * 10.0**(6-i)
            tval = zval / real(nval)
            if(i.ge.12) lfmt = '(f8.8)'
            if(i.eq.11) lfmt = '(f8.7)'
            if(i.eq.10) lfmt = '(f8.6)'
            if(i.eq.9)  lfmt = '(f8.5)'
            if(i.eq.8)  lfmt = '(f8.4)'
            if(i.eq.7)  lfmt = '(f8.3)'
            if(i.eq.6)  lfmt = '(f8.2)'
            if(i.eq.5)  lfmt = '(f8.1)'
            if(i.le.4)  lfmt = '(f8.0)'
c
c      Loop along the axis drawing the tic marks and labels:
c
            yloc(1) = yaxmin
            pos     = xmin
            num     = 2
            do i=1,100
                  yloc(2) = yaxmin - tlng
                  xloc(1) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                  xloc(2) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                  call psline(num,xloc,yloc,pl2,0)
                  yloc(2) = yloc(1) - 2.5*tlng
                  write(label,lfmt) pos                        
                  call pstext(xloc(1),yloc(2),8,label,psz,1,0.d0,1)
                  yloc(2) = yaxmin - tsht
                  do j=1,nval-1
                       pos     = pos + tval
                       if(pos.gt.xmax) go to 2
                       xloc(1) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                       xloc(2) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                       call psline(num,xloc,yloc,pl2,0)
                    end do
                  pos = pos + tval
                  if(pos.gt.xmax) go to 2
            end do
 2          continue
      endif
c
c CONSTRUCT THE Y AXIS:
c
c
c Log scale?
c
      if(ilog.eq.2.or.ilog.eq.3) then
c
c      The start, end, number of log(10) cycles, the tic mark start,
c      and the number of points in defining a tic mark:
c
            tminy   = alog10(ymin)
            tmaxy   = alog10(ymax)
            ncyc    = tmaxy - tminy
            cbas    = ymin/10
            xloc(1) = xaxmin
            num     = 2
c
c      Loop along the axis drawing the tic marks and labels:
c
            do icyc=1,ncyc+1
                  cbas = cbas * 10
                  do i=1,9
                  t1   = alog10(cbas*real(i))
                  yloc(1) = resc(tminy,tmaxy,yaxmin,yaxmax,t1)
                  yloc(2) = yloc(1)
                  if(i.eq.1) then
c
c            First point - long tic mark:
c
                        xloc(2) = xloc(1) - tlng
                        call psline(num,xloc,yloc,pl2,0)
                        xloc(2) = xloc(2) - 0.1*tlng
                        if(abs(t1+9.).le.EPS) label = '1.0e-9  '
                        if(abs(t1+8.).le.EPS) label = '1.0e-8  '
                        if(abs(t1+7.).le.EPS) label = '1.0e-7  '
                        if(abs(t1+6.).le.EPS) label = '1.0e-6  '
                        if(abs(t1+5.).le.EPS) label = '0.00001 '
                        if(abs(t1+4.).le.EPS) label = '0.0001  '
                        if(abs(t1+3.).le.EPS) label = '0.001   '
                        if(abs(t1+2.).le.EPS) label = '0.01    '
                        if(abs(t1+1.).le.EPS) label = '0.1     '
                        if(abs(t1)   .le.EPS) label = '1       '
                        if(abs(t1-1.).le.EPS) label = '10      '
                        if(abs(t1-2.).le.EPS) label = '100     '
                        if(abs(t1-3.).le.EPS) label = '1000    '
                        if(abs(t1-4.).le.EPS) label = '10000   '
                        if(abs(t1-5.).le.EPS) label = '100000  '
                        if(abs(t1-6.).le.EPS) label = '1.0e+6  '
                        if(abs(t1-7.).le.EPS) label = '1.0e+7  '
                        if(abs(t1-8.).le.EPS) label = '1.0e+8  '
                        if(abs(t1-9.).le.EPS) label = '1.0e+9  '
                        call pstext(xloc(2),yloc(2),8,label,
     +                                      psz,1,0.0,2)
                  else
c
c            Not first point - short tic mark:
c
                        if(icyc.le.ncyc) then
                              xloc(2) = xloc(1) - tsht
                              call psline(num,xloc,yloc,pl2,0)
                        endif
                  endif
                  end do
            end do
      else
c
c      Determine a labelling and tic mark increment:
c
            do i=1,20
                  test = (ymax-ymin)/(10.0**(6-i))
                  if(test.gt.0.9) go to 11
            end do
 11         if(test.ge.3.0)                 zval = 1.0
            if(test.le.3.0.and.test.gt.2.0) zval = 0.5
            if(test.le.2.0.and.test.gt.1.2) zval = 0.4
            if(test.le.1.2)                 zval = 0.2
            nval = 5
            if(zval.eq.0.4.or.zval.eq.0.2)  nval = 4
            zval = zval * 10.0**(6-i)
            tval = zval / real(nval)
            if(i.ge.12) lfmt = '(f8.8)'
            if(i.eq.11) lfmt = '(f8.7)'
            if(i.eq.10) lfmt = '(f8.6)'
            if(i.eq.9)  lfmt = '(f8.5)'
            if(i.eq.8)  lfmt = '(f8.4)'
            if(i.eq.7)  lfmt = '(f8.3)'
            if(i.eq.6)  lfmt = '(f8.2)'
            if(i.eq.5)  lfmt = '(f8.1)'
            if(i.le.4)  lfmt = '(f8.0)'
c
c      Loop along the axis drawing the tic marks and labels:
c
            xloc(1) = xaxmin
            pos     = ymin
            num     = 2
            do i=1,100
                  xloc(2) = xaxmin - tlng
                  yloc(1) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                  yloc(2) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                  call psline(num,xloc,yloc,pl2,0)
                  xloc(2) = xloc(2) - 0.2*tlng
                  write(label,lfmt) pos                        
                  call pstext(xloc(2),yloc(2),8,label,psz,1,0.0,2)
                  xloc(2) = xaxmin - tsht
                  do j=1,nval-1
                       pos     = pos + tval
                       if(pos.gt.ymax) go to 12
                       yloc(1) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                       yloc(2) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                       call psline(num,xloc,yloc,pl2,0)
                     end do
                  pos = pos + tval
                  if(pos.gt.ymax) go to 12
            end do
 12         continue
      endif
c
c Return to calling program:
c
      return
      end
