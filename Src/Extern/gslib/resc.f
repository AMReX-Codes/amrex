      real function resc(xmin1,xmax1,xmin2,xmax2,x111)

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      real*8 rsc
c
c Simple linear rescaling (get a value in coordinate system "2" given
c a value in "1"):
c
      rsc  = dble((xmax2-xmin2)/(xmax1-xmin1))
      resc = xmin2 + real( dble(x111 - xmin1) * rsc )
      return
      end
