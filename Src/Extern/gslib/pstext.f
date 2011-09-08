      subroutine pstext(xs,ys,lostr,str,tsiz,ifont,rot,iadj)
c-----------------------------------------------------------------------
c
c              Write Postscript Text commands to a file
c              ****************************************
c
c
c CALLING ARGUMENTS:
c
c  xs            starting value of x in the range xmin to xmax
c  ys            starting value of y in the range ymin to ymax
c  lostr      number of characters in str to print
c  str            the character string
c  tsiz            Text size in 1/72 of an inch
c  ifont      Font Number: See font number below
c  rot             Rotation Angle to post the text (default to 0.0)
c  iadj            Adjustment: 0=left adjusted, 1=centre, 2=right
c
c
c-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      character str*80,fnnt(10)*32,line*132,size*4,part1*1,part2*7
c
c Common Block for Postscript Output Unit and Scaling:
c
      common /psdata/ pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                xmax,ymin,ymax,lpsout
      save    fnnt,ifold,tsold,izip
c
c Preset 10 different fonts:
c
      data fnnt/'/Helvetica             findfont ',
     +          '/Helvetica-Bold        findfont ',
     +          '/Helvetica-BoldOblique findfont ',
     +          '/Times-Roman           findfont ',
     +          '/Times-Bold            findfont ',
     +          '/Times-Italic          findfont ',
     +          '/Times-BoldItalic      findfont ',
     +          '/Courier               findfont ',
     +          '/Courier-Bold          findfont ',
     +          '/Courier-BoldOblique   findfont '/
      data ifold/0/,tsold/0.0/,izip/0/
      part1 = '('
      part2 = ')  text'
c
c Remove leading and trailing blanks:
c
      lost = lostr
      do i=1,lostr
            if(str(1:1).eq.' ') then
                  lost = lost - 1
                  do j=1,lost
                        k = j + 1
                        str(j:j) = str(k:k)
                  end do
            else
                  go to 1
            endif
      end do
 1    k = lost
      do i=1,k
            ix = k - i + 1
            if(str(ix:ix).ne.' ') go to 2
            lost = lost - 1
      end do
 2    if(lost.le.0) return
c
c Create line to set the text size and type:
c
      if(ifont.ne.ifold.or.tsiz.ne.tsold) then
            isiz=int(tsiz/pscl)
            write(size,'(i4)') isiz
            line=fnnt(ifont)//size//' scalefont setfont'      
            write(lpsout,'(a)')line(1:54)
            ifold = ifont
            tsold = tsiz
      endif
c
c Set the correct adjustment:
c
      part2(3:3) = 'l'
      if(iadj.eq.1) part2(3:3) = 'c'
      if(iadj.eq.2) part2(3:3) = 'r'
c
c Write the lines and position to the Postscript file:
c                  
      ix = int((resc(xmin,xmax,pxmin,pxmax,xs))/pscl)
      iy = int((resc(ymin,ymax,pymin,pymax,ys))/pscl)
c
c Rotate if Necessary:
c
      line = part1//str(1:lost)//part2
      if(rot.ne.0.0) then
            irot = int(rot)
            write(lpsout,102)   ix,iy
            write(lpsout,103)   irot 
            write(lpsout,100)   izip,izip
            write(lpsout,'(a)') line(1:lost+8)
            ix   = -1.0 * ix
            iy   = -1.0 * iy
            irot = -1.0 * irot
            write(lpsout,103)   irot 
            write(lpsout,102)   ix,iy
      else
c
c Just write out the text if no rotation:
c
            write(lpsout,100)   ix,iy
            write(lpsout,'(a)') line(1:lost+8)
      endif
 100  format(i5,1x,i5,1x,'m')
 102  format(i5,1x,i5,1x,'translate')
 103  format(i5,1x,'rotate')
c
c Finished - Return to calling program:
c
      return
      end
