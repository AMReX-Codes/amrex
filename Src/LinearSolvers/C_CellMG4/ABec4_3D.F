#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>

#include "ABec4_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "LO_BCTYPES.H"

c-----------------------------------------------------------------------
c
c     Fill in a matrix x vector operator here
c
      subroutine FORT_ADOTX(
     $     y,DIMS(y),
     $     x,DIMS(x),
     $     alpha, beta,
     $     a,DIMS(a),
     $     b,DIMS(b),
     $     lo,hi,nc,
     $     h
     $     )

      implicit none

      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(y)
      integer DIMDEC(x)
      integer DIMDEC(a)
      integer DIMDEC(b)
      REAL_T  y(DIMV(y),nc)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T  b(DIMV(b))
      REAL_T h(BL_SPACEDIM)
      integer i,j,n
      REAL_T b1L, b1R, b1B, b1T, g1L, g1R, g1B, g1T, FL, FR, FB, FT
      REAL_T b2L, b2R, b2B, b2T, g2L, g2R, g2B, g2T
      REAL_T b3L, b3R, b3B, b3T, g3L, g3R, g3B, g3T
      REAL_T btm2, btm1, btp1, btp2, btzL, btzR, btzB, btzT
      REAL_T gtm2, gtm1, gtp1, gtp2, gtzL, gtzR, gtzB, gtzT
      REAL_T i144, i80, i12, bh1i, bh2i
      REAL_T h11, h12, h13, h14
      REAL_T h21, h22, h23, h24
      REAL_T h1i1_12, h1i1_48, h1i2_16, h1i3_12, h12_12
      REAL_T h2i1_12, h2i1_48, h2i2_16, h2i3_12, h22_12

      h11 = h(1)
      h12 = h11*h11
      h13 = h12*h11
      h14 = h13*h11
      h1i1_12 = 1.d0 / (12.d0*h11)
      h1i1_48 = 1.d0 / (48.d0*h11)
      h1i2_16 = 1.d0 / (16.d0*h12)
      h1i3_12 = 1.d0 / (12.d0*h13)
      h12_12  = h12 / 12.d0

      h21 = h(2)
      h22 = h21*h21
      h23 = h22*h21
      h24 = h23*h21
      h2i1_12 = 1.d0 / (12.d0*h21)
      h2i1_48 = 1.d0 / (48.d0*h21)
      h2i2_16 = 1.d0 / (16.d0*h22)
      h2i3_12 = 1.d0 / (12.d0*h23)
      h22_12  = h22 / 12.d0

      i12  = 1.d0/12.d0
      i144 = 1.d0/144.d0
      i80  = 1.d0/80.d0
      bh1i = beta / h11
      bh2i = beta / h21

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            btm2 = btildei(i-1,j-2)
            btm1 = btildei(i-1,j-1)
            btzL = btildei(i-1,j  )
            btp1 = btildei(i-1,j+1)
            btp2 = btildei(i-1,j+2)

            b1L = ( 5*(btm2 - btp2) + 34*(btp1 - btm1)          ) * h1i1_48
            b2L = (-1*(btm2 + btp2) + 12*(btp1 + btm1) - 22*btzL) * h1i2_16
            b3L = (  - btm2 + btp2  -  2*(btp1 - btm1)          ) * h1i3_12

            btm2 = btildei(i  ,j-2)
            btm1 = btildei(i  ,j-1)
            btzR = btildei(i  ,j  )
            btp1 = btildei(i  ,j+1)
            btp2 = btildei(i  ,j+2)

            b1R = ( 5*(btm2 - btp2) + 34*(btp1 - btm1)          ) * h1i1_48
            b2R = (-1*(btm2 + btp2) + 12*(btp1 + btm1) - 22*btzR) * h1i2_16
            b3R = (  - btm2 + btp2  -  2*(btp1 - btm1)          ) * h1i3_12

            btm2 = btildej(i-2,j-1)
            btm1 = btildej(i-1,j-1)
            btzB = btildej(i  ,j-1)
            btp1 = btildej(i+1,j-1)
            btp2 = btildej(i+2,j-1)

            b1B = ( 5*(btm2 - btp2) + 34*(btp1 - btm1)          ) * h2i1_48
            b2B = (-1*(btm2 + btp2) + 12*(btp1 + btm1) - 22*btzB) * h2i2_16
            b3B = (  - btm2 + btp2  -  2*(btp1 - btm1)          ) * h2i3_12

            btm2 = btildej(i-2,j  )
            btm1 = btildej(i-1,j  )
            btzT = btildej(i  ,j  )
            btp1 = btildej(i+1,j  )
            btp2 = btildej(i+2,j  )

            b1T = ( 5*(btm2 - btp2) + 34*(btp1 - btm1)          ) * h2i1_48
            b2T = (-1*(btm2 + btp2) + 12*(btp1 + btm1) - 22*btzT) * h2i2_16
            b3T = (  - btm2 + btp2  -  2*(btp1 - btm1)          ) * h2i3_12

            do n = 1, nc

               gtm2 = gxtildei(i-1,j-2,n)
               gtm1 = gxtildei(i-1,j-1,n)
               gtzL = gxtildei(i-1,j  ,n)
               gtp1 = gxtildei(i-1,j+1,n)
               gtp2 = gxtildei(i-1,j+2,n)

               g1L = ( 5*(gtm2 - gtp2) + 34*(gtp1 - gtm1)          ) * h1i1_48
               g2L = (-1*(gtm2 + gtp2) + 12*(gtp1 + gtm1) - 22*gtzL) * h1i2_16
               g3L = (  - gtm2 + gtp2  -  2*(gtp1 - gtm1)          ) * h1i3_12

               gtm2 = gxtildei(i  ,j-2,n)
               gtm1 = gxtildei(i  ,j-1,n)
               gtzR = gxtildei(i  ,j  ,n)
               gtp1 = gxtildei(i  ,j+1,n)
               gtp2 = gxtildei(i  ,j+2,n)

               g1R = ( 5*(gtm2 - gtp2) + 34*(gtp1 - gtm1)          ) * h1i1_48
               g2R = (-1*(gtm2 + gtp2) + 12*(gtp1 + gtm1) - 22*gtzR) * h1i2_16
               g3R = (  - gtm2 + gtp2  -  2*(gtp1 - gtm1)          ) * h1i3_12

               gtm2 = gytildej(i-2,j-1,n)
               gtm1 = gytildej(i-1,j-1,n)
               gtzB = gytildej(i  ,j-1,n)
               gtp1 = gytildej(i+1,j-1,n)
               gtp2 = gytildej(i+2,j-1,n)

               g1B = ( 5*(gtm2 - gtp2) + 34*(gtp1 - gtm1)          ) * h2i1_48
               g2B = (-1*(gtm2 + gtp2) + 12*(gtp1 + gtm1) - 22*gtzB) * h2i2_16
               g3B = (  - gtm2 + gtp2  -  2*(gtp1 - gtm1)          ) * h2i3_12

               gtm2 = gytildej(i-2,j  ,n)
               gtm1 = gytildej(i-1,j  ,n)
               gtzT = gytildej(i  ,j  ,n)
               gtp1 = gytildej(i+1,j  ,n)
               gtp2 = gytildej(i+2,j  ,n)

               g1T = ( 5*(gtm2 - gtp2) + 34*(gtp1 - gtm1)          ) * h2i1_48
               g2T = (-1*(gtm2 + gtp2) + 12*(gtp1 + gtm1) - 22*gtzT) * h2i2_16
               g3T = (  - gtm2 + gtp2  -  2*(gtp1 - gtm1)          ) * h2i3_12

               FL = gtzL*btzL + g1L*b1L*h12_12 + h14*( - g2L*b2L*i144 + ( g1L*b3L + g2L*b2L + g3L*b1L )*i80 )
               FR = gtzR*btzR + g1R*b1R*h12_12 + h14*( - g2R*b2R*i144 + ( g1R*b3R + g2R*b2R + g3R*b1R )*i80 )
               FB = gtzB*btzB + g1B*b1B*h22_12 + h24*( - g2B*b2B*i144 + ( g1B*b3B + g2B*b2B + g3B*b1B )*i80 )
               FT = gtzT*btzT + g1T*b1T*h22_12 + h24*( - g2T*b2T*i144 + ( g1T*b3T + g2T*b2T + g3T*b1T )*i80 )

               y(i,j,n) = alpha*a(i,j)*x(i,j,n) - ( bh1i*(FR - FL) + bh2i*(FT - FB) )

            end do
         end do
      end do

      contains

c     4th order approx to b at (i+1/2,j), given cell-averaged values at ([i-1,i,i+1,i+2],j)
      REAL_T function btildei(i,j)
      integer, intent(in) :: i, j
      btildei = (-b(i-1,j)+7*(b(i,j)+b(i+1,j))-b(i+2,j))*i12
      end function

      REAL_T function btildej(i,j)
      integer, intent(in) :: i,j
      btildej = (-b(i,j-1)+7*(b(i,j)+b(i,j+1))-b(i,j+2))*i12
      end function

c     Average over edge (i+1/2,j) of derivative normal to the edge, based on cell-averaged values
      REAL_T function gxtildei(i,j,n)
      integer, intent(in) :: i,j,n
      gxtildei = (x(i-1,j,n)-x(i+2,j,n)+15.d0*(x(i+1,j,n)-x(i,j,n)))*h1i1_12
      end function

      REAL_T function gytildej(i,j,n)
      integer, intent(in) :: i,j,n
      gytildej = (x(i,j-1,n)-x(i,j+2,n)+15.d0*(x(i,j+1,n)-x(i,j,n)))*h2i1_12
      end function

      end

c-----------------------------------------------------------------------
c
c     Fill in fluxes
c
      subroutine FORT_FLUX(
     $     x,DIMS(x),
     $     alpha, beta,
     $     a, DIMS(a),
     $     b, DIMS(b),
     $     lo,hi,nc,
     $     h,
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux)
     $     )
      implicit none
      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(x)
      integer DIMDEC(a)
      integer DIMDEC(b)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T  b(DIMV(b))
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T h(BL_SPACEDIM)

      integer i,j,n
      REAL_T btm2, btm1, btz, btp1, btp2, b1, b2, b3
      REAL_T gtm2, gtm1, gtz, gtp1, gtp2, g1, g2, g3

      REAL_T i144, i80, i12, bh1i, bh2i
      REAL_T h11, h12, h13, h14
      REAL_T h21, h22, h23, h24
      REAL_T h1i1_12, h1i1_48, h1i2_16, h1i3_12, h12_12
      REAL_T h2i1_12, h2i1_48, h2i2_16, h2i3_12, h22_12

      h11 = h(1)
      h12 = h11*h11
      h13 = h12*h11
      h14 = h13*h11
      h1i1_12 = 1.d0 / (12.d0*h11)
      h1i1_48 = 1.d0 / (48.d0*h11)
      h1i2_16 = 1.d0 / (16.d0*h12)
      h1i3_12 = 1.d0 / (12.d0*h13)
      h12_12  = h12 / 12.d0

      h21 = h(2)
      h22 = h21*h21
      h23 = h22*h21
      h24 = h23*h21
      h2i1_12 = 1.d0 / (12.d0*h21)
      h2i1_48 = 1.d0 / (48.d0*h21)
      h2i2_16 = 1.d0 / (16.d0*h22)
      h2i3_12 = 1.d0 / (12.d0*h23)
      h22_12  = h22 / 12.d0

      i12  = 1.d0/12.d0
      i144 = 1.d0/144.d0
      i80  = 1.d0/80.d0
      bh1i = beta / h11
      bh2i = beta / h21

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            
            btm2 = btildei(i-1,j-2)
            btm1 = btildei(i-1,j-1)
            btz  = btildei(i-1,j  )
            btp1 = btildei(i-1,j+1)
            btp2 = btildei(i-1,j+2)

            b1 = ( 5*(btm2 - btp2) + 34*(btp1 - btm1)         ) * h1i1_48
            b2 = (-1*(btm2 + btp2) + 12*(btp1 + btm1) - 22*btz) * h1i2_16
            b3 = (  - btm2 + btp2  -  2*(btp1 - btm1)         ) * h1i3_12

            do n = 1, nc

               gtm2 = gxtildei(i-1,j-2,n)
               gtm1 = gxtildei(i-1,j-1,n)
               gtz  = gxtildei(i-1,j  ,n)
               gtp1 = gxtildei(i-1,j+1,n)
               gtp2 = gxtildei(i-1,j+2,n)

               g1 = ( 5*(gtm2 - gtp2) + 34*(gtp1 - gtm1)         ) * h1i1_48
               g2 = (-1*(gtm2 + gtp2) + 12*(gtp1 + gtm1) - 22*gtz) * h1i2_16
               g3 = (  - gtm2 + gtp2  -  2*(gtp1 - gtm1)         ) * h1i3_12

               xflux(i,j,n) = gtz*btz + g1*b1*h12_12 + h14*( - g2*b2*i144 + ( g1*b3 + g2*b2 + g3*b1 )*i80 ) 

            end do
         end do
      end do

      do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            
            btm2 = btildej(i-2,j-1)
            btm1 = btildej(i-1,j-1)
            btz  = btildej(i  ,j-1)
            btp1 = btildej(i+1,j-1)
            btp2 = btildej(i+2,j-1)

            b1 = ( 5*(btm2 - btp2) + 34*(btp1 - btm1)         ) * h2i1_48
            b2 = (-1*(btm2 + btp2) + 12*(btp1 + btm1) - 22*btz) * h2i2_16
            b3 = (  - btm2 + btp2  -  2*(btp1 - btm1)         ) * h2i3_12

            do n = 1, nc

               gtm2 = gytildej(i-2,j-1,n)
               gtm1 = gytildej(i-1,j-1,n)
               gtz  = gytildej(i  ,j-1,n)
               gtp1 = gytildej(i+1,j-1,n)
               gtp2 = gytildej(i+2,j-1,n)

               g1 = ( 5*(gtm2 - gtp2) + 34*(gtp1 - gtm1)         ) * h2i1_48
               g2 = (-1*(gtm2 + gtp2) + 12*(gtp1 + gtm1) - 22*gtz) * h2i2_16
               g3 = (  - gtm2 + gtp2  -  2*(gtp1 - gtm1)         ) * h2i3_12

               yflux(i,j,n) = gtz*btz + g1*b1*h22_12 + h24*( - g2*b2*i144 + ( g1*b3 + g2*b2 + g3*b1 )*i80 ) 

            end do
         end do
      end do

      contains

c     4th order approx to b at (i+1/2,j), given cell-averaged values at ([i-1,i,i+1,i+2],j)
      REAL_T function btildei(i,j)
      integer, intent(in) :: i, j
      btildei = (-b(i-1,j)+7*(b(i,j)+b(i+1,j))-b(i+2,j))*i12
      end function

      REAL_T function btildej(i,j)
      integer, intent(in) :: i,j
      btildej = (-b(i,j-1)+7*(b(i,j)+b(i,j+1))-b(i,j+2))*i12
      end function

c     Average over edge (i+1/2,j) of derivative normal to the edge, based on cell-averaged values
      REAL_T function gxtildei(i,j,n)
      integer, intent(in) :: i,j,n
      gxtildei = (x(i-1,j,n)-x(i+2,j,n)+15.d0*(x(i+1,j,n)-x(i,j,n)))*h1i1_12
      end function

      REAL_T function gytildej(i,j,n)
      integer, intent(in) :: i,j,n
      gytildej = (x(i,j-1,n)-x(i,j+2,n)+15.d0*(x(i,j+1,n)-x(i,j,n)))*h2i1_12
      end function

      end

c-----------------------------------------------------------------------
      subroutine FORT_APPLYBC4 (
     $     flagden, flagbc, maxorder,
     $     phi,   DIMS(phi),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     den,   DIMS(den),
     $     lo, hi, nc,
     $     h
     $     )

      implicit none
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c      
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.
c      
      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den))
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer bct
      REAL_T bcl, tmp, b, d
      parameter (b = 1.d0/24.d0)
      parameter (d = 11.d0/12.d0)
      REAL_T h(BL_SPACEDIM)
c
      integer i
      integer j
      integer n
      logical is_dirichlet
      logical is_neumann
c
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
c
c
c     The Left face of the grid
c
      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  if (mask(lo(1)-1,j) .gt. 0) then
c                     phi(lo(1)-1,j,n) =  (9*phi(lo(1),j,n) +  3*phi(lo(1)+1,j,n) - phi(lo(1)+2,j,n))/11.d0
                     phi(lo(1)-1,j,n) =  (5*phi(lo(1),j,n) +  9*phi(lo(1)+1,j,n) - 5*phi(lo(1)+2,j,n) + phi(lo(1)+3,j,n))/10.d0
                  endif
                  if (mask(lo(1)-2,j) .gt. 0) then
c                     phi(lo(1)-2,j,n) = (- 30*phi(lo(1),j,n) + 56*phi(lo(1)+1,j,n) - 15*phi(lo(1)+2,j,n))/11.d0
                     phi(lo(1)-2,j,n) = (- 75*phi(lo(1),j,n) + 145*phi(lo(1)+1,j,n) - 75*phi(lo(1)+2,j,n) + 15*phi(lo(1)+3,j,n))/10.d0
                  endif
               end do
            end do
            if ( flagden .eq. 1) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = 1.0D0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     if (mask(lo(1)-1,j) .gt. 0) then
                        if (mask(lo(1)-1,j-1).gt.0 .and. mask(lo(1)-1,j+1).gt.0) then
                           tmp = b*(bcval(lo(1)-1,j+1,n)+bcval(lo(1)-1,j-1,n)) + d*bcval(lo(1)-1,j,n)
                        else
                           tmp = bcval(lo(1)-1,j,n)
                        endif
c                        phi(lo(1)-1,j,n) = (12*tmp - 13*phi(lo(1),j,n) + 5*phi(lo(1)+1,j,n) - phi(lo(1)+2,j,n))/3.d0
                        phi(lo(1)-1,j,n) = (60*tmp - 77*phi(lo(1),j,n) + 43*phi(lo(1)+1,j,n) - 17*phi(lo(1)+2,j,n) + 3*phi(lo(1)+3,j,n))/12.d0
                        if (mask(lo(1)-2,j) .gt. 0) then
c                           phi(lo(1)-2,j,n) = (48*tmp - 70*phi(lo(1),j,n) + 32*phi(lo(1)+1,j,n) - 7*phi(lo(1)+2,j,n))/3.d0
                           phi(lo(1)-2,j,n) = (300*tmp - 505*phi(lo(1),j,n) + 335*phi(lo(1)+1,j,n) - 145*phi(lo(1)+2,j,n) + 27*phi(lo(1)+3,j,n))/12.d0
                        endif
                     endif
                  end do
               else
                  do j = lo(2), hi(2)
                     if (mask(lo(1)-1,j) .gt. 0) then
                        phi(lo(1)-1, j, n) = 0.d0
                     endif
                     if (mask(lo(1)-2,j) .gt. 0) then
                        phi(lo(1)-2, j, n) = 0.d0
                     endif
                  end do
               end if
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = merge(66.D0, 0.0D0,
     $                 mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(lo(1)-1, j, n) = merge(
     $                -phi(lo(1),j,n),
     $                 phi(lo(1)-1, j, n),
     $                 mask(lo(1)-1,j) .gt. 0)
                  phi(lo(1)-2, j, n) = merge(
     $                -phi(lo(1)+1,j,n),
     $                 phi(lo(1)-2, j, n),
     $                 mask(lo(1)-2,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = merge(-1.0D0, 0.0D0,
     $                 mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else 
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c     
c     The Right face of the grid
c
      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  if (mask(hi(1)+1,j) .gt. 0) then
c                     phi(hi(1)+1,j,n) =  (9*phi(hi(1),j,n) +  3*phi(hi(1)-1,j,n) - phi(hi(1)-2,j,n))/11.d0
                     phi(hi(1)+1,j,n) =  (5*phi(hi(1)  ,j,n) +  9*phi(hi(1)-1,j,n) - 5*phi(hi(1)-2,j,n) + phi(hi(1)-3,j,n))/10.d0
                  endif
                  if (mask(hi(1)+2,j) .gt. 0) then
c                     phi(hi(1)+2,j,n) = (- 30*phi(hi(1),j,n) + 56*phi(hi(1)-1,j,n) - 15*phi(hi(1)-2,j,n))/11.d0
                     phi(hi(1)+2,j,n) = (- 75*phi(hi(1),j,n) + 145*phi(hi(1)-1,j,n) - 75*phi(hi(1)-2,j,n) + 15*phi(hi(1)-3,j,n))/10.d0
                  endif
               end do
            end do
	    if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j) = 1.0D0
               end do
	    end if
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     if (mask(hi(1)+1,j) .gt. 0) then
                        if (mask(hi(1)+1,j-1).gt.0 .and. mask(hi(1)+1,j+1).gt.0) then
                           tmp = b*(bcval(hi(1)+1,j+1,n)+bcval(hi(1)+1,j-1,n)) + d*bcval(hi(1)+1,j,n)
                        else
                           tmp = bcval(hi(1)+1,j,n)
                        endif
c                        phi(hi(1)+1,j,n) = (12*tmp - 13*phi(hi(1),j,n) + 5*phi(hi(1)-1,j,n) - phi(hi(1)-2,j,n))/3.d0
                        phi(hi(1)+1,j,n) = (60*tmp - 77*phi(hi(1),j,n) + 43*phi(hi(1)-1,j,n) - 17*phi(hi(1)-2,j,n) + 3*phi(hi(1)-3,j,n))/12.d0
                        if (mask(hi(1)+2,j) .gt. 0) then
c                            phi(hi(1)+2,j,n) = (48*tmp - 70*phi(hi(1),j,n) + 32*phi(hi(1)-1,j,n) - 7*phi(hi(1)-2,j,n))/3.d0
                            phi(hi(1)+2,j,n) = (300*tmp - 505*phi(hi(1),j,n) + 335*phi(hi(1)-1,j,n) - 145*phi(hi(1)-2,j,n) + 27*phi(hi(1)-3,j,n))/12.d0
                        endif
                     endif
                  end do
               else
                  do j = lo(2), hi(2)
                     if (mask(hi(1)+1,j) .gt. 0) then
                        phi(hi(1)+1, j, n) = 0.d0
                     endif
                     if (mask(hi(1)+2,j) .gt. 0) then
                        phi(hi(1)+2, j, n) = 0.d0
                     endif
                  end do
               end if
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j)   = merge(66.D0, 0.0D0,
     $                 mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(hi(1)+1, j, n) = merge(
     $                -phi(hi(1),j,n),
     $                 phi(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
                  phi(hi(1)+2, j, n) = merge(
     $                -phi(hi(1)-1,j,n),
     $                 phi(hi(1)+2, j, n),
     $                 mask(hi(1)+2,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j) = merge(-1.0D0, 0.0D0,
     $                 mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
c     The Bottom of the Grid
c
      if(cdir .eq. 1) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1),hi(1)
                  if (mask(i,lo(2)-1) .gt. 0) then
c                     phi(i,lo(2)-1,n) =  (9*phi(i,lo(2),n) +  3*phi(i,lo(2)+1,n) - phi(i,lo(2)+2,n))/11.d0
                     phi(i,lo(2)-1,n) =  (5*phi(i,lo(2)  ,n) +  9*phi(i,lo(2)+1,n) - 5*phi(i,lo(2)+2,n) + phi(i,lo(2)+3,n))/10.d0
                  endif
                  if (mask(i,lo(2)-2) .gt. 0) then
c                     phi(i,lo(2)-2,n) = (- 30*phi(i,lo(2),n) + 56*phi(i,lo(2)+1,n) - 15*phi(i,lo(2)+2,n))/11.d0
                     phi(i,lo(2)-2,n) = (- 75*phi(i,lo(2),n) + 145*phi(i,lo(2)+1,n) - 75*phi(i,lo(2)+2,n) + 15*phi(i,lo(2)+3,n))/10.d0
                  endif
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1),hi(1)
                  den(i,lo(2))   = 1.0D0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     if (mask(i,lo(2)-1) .gt. 0) then
                        if (mask(i-1,lo(2)-1).gt.0 .and. mask(i+1,lo(2)-1).gt.0) then
                           tmp = b*(bcval(i+1,lo(2)-1,n)+bcval(i-1,lo(2)-1,n)) + d*bcval(i,lo(2)-1,n)
                        else
                           tmp = bcval(i,lo(2)-1,n)
                        endif
c                        phi(i,lo(2)-1,n) = (12*tmp - 13*phi(i,lo(2),n) + 5*phi(i,lo(2)+1,n) - phi(i,lo(2)+2,n))/3.d0
                        phi(i,lo(2)-1,n) = (60*tmp - 77*phi(i,lo(2),n) + 43*phi(i,lo(2)+1,n) - 17*phi(i,lo(2)+2,n) + 3*phi(i,lo(2)+3,n))/12.d0
                        if (mask(i,lo(2)-2) .gt. 0) then
c                           phi(i,lo(2)-2,n) =   (48*tmp - 70*phi(i,lo(2),n) + 32*phi(i,lo(2)+1,n) - 7*phi(i,lo(2)+2,n))/3.d0
                           phi(i,lo(2)-2,n) = (300*tmp - 505*phi(i,lo(2),n) + 335*phi(i,lo(2)+1,n) - 145*phi(i,lo(2)+2,n) + 27*phi(i,lo(2)+3,n))/12.d0
                        endif
                     endif
                  end do
               else
                  do i = lo(1), hi(1)
                     if (mask(i,lo(2)-1) .gt. 0) then
                        phi(i,lo(2)-1, n) = 0.d0
                     endif
                     if (mask(i,lo(2)-2) .gt. 0) then
                        phi(i,lo(2)-2, n) = 0.d0
                     endif
                  end do
               end if
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i, lo(2))   = merge(66.D0, 0.0D0,
     $                 mask(i, lo(2)-1) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,lo(2)-1,n) = merge(
     $                -phi(i,lo(2),n),
     $                 phi(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
                  phi(i,lo(2)-2,n) = merge(
     $                -phi(i,lo(2)+1,n),
     $                 phi(i,lo(2)-2,n),
     $                 mask(i,lo(2)-2) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,lo(2)) = merge(-1.0D0, 0.0D0,
     $                 mask(i,lo(2)-1) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
c     The top of the grid
c
      if (cdir .eq. 3) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1), hi(1)
                  if (mask(i,hi(2)+1) .gt. 0) then
c                     phi(i,hi(2)+1,n) =  (9*phi(i,hi(2),n) +  3*phi(i,hi(2)-1,n) - phi(i,hi(2)-2,n))/11.d0
                     phi(i,hi(2)+1,n) =  (5*phi(i,hi(2),n) +  9*phi(i,hi(2)-1,n) - 5*phi(i,hi(2)-2,n) + phi(i,hi(2)-3,n))/10.d0
                  endif
                  if (mask(i,hi(2)+2) .gt. 0) then
c                     phi(i,hi(2)+2,n) = (- 30*phi(i,hi(2),n) + 56*phi(i,hi(2)-1,n) - 15*phi(i,hi(2)-2,n))/11.d0
                     phi(i,hi(2)+2,n) = (- 75*phi(i,hi(2),n) + 145*phi(i,hi(2)-1,n) - 75*phi(i,hi(2)-2,n) + 15*phi(i,hi(2)-3,n))/10.d0
                  endif
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2))   = 1.0D0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     if (mask(i,hi(2)+1) .gt. 0) then
                        if (mask(i-1,hi(2)+1).gt.0 .and. mask(i+1,hi(2)+1).gt.0) then
                           tmp = b*(bcval(i+1,hi(2)+1,n)+bcval(i-1,hi(2)+1,n)) + d*bcval(i,hi(2)+1,n)
                        else
                           tmp = bcval(i,hi(2)+1,n)
                        endif
c                        phi(i,hi(2)+1,n) = (12*tmp - 13*phi(i,hi(2),n) + 5*phi(i,hi(2)-1,n) - phi(i,hi(2)-2,n))/3.d0
                        phi(i,hi(2)+1,n) = (60*tmp - 77*phi(i,hi(2),n) + 43*phi(i,hi(2)-1,n) - 17*phi(i,hi(2)-2,n) + 3*phi(i,hi(2)-3,n))/12.d0
                        if (mask(i,hi(2)+2) .gt. 0) then
c                           phi(i,hi(2)+2,n) =   (48*tmp -  70*phi(i,hi(2),n) + 32*phi(i,hi(2)-1,n) - 7*phi(i,hi(2)-2,n))/3.d0
                           phi(i,hi(2)+2,n) = (300*tmp - 505*phi(i,hi(2),n) + 335*phi(i,hi(2)-1,n) - 145*phi(i,hi(2)-2,n) + 27*phi(i,hi(2)-3,n))/12.d0
                        endif
                     endif
                  end do
               else
                  do i = lo(1), hi(1)
                     if (mask(i,hi(2)+1) .gt. 0) then
                        phi(i,hi(2)+1, n) = 0.d0
                     endif
                     if (mask(i,hi(2)+2) .gt. 0) then
                        phi(i,hi(2)+2, n) = 0.d0
                     endif
                  end do
               end if
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2))   = merge(66.D0, 0.0D0,
     $                 mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,hi(2)+1,n) = merge(
     $                -phi(i,hi(2),n),
     $                 phi(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
                  phi(i,hi(2)+2,n) = merge(
     $                -phi(i,hi(2)-2,n),
     $                 phi(i,hi(2)+2,n),
     $                 mask(i,hi(2)+2) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2)) = merge(-1.0D0, 0.0D0,
     $                 mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
      end

c-----------------------------------------------------------------------

      subroutine FORT_APPLYBC4_TOUCHUP (
     $     phi,   DIMS(phi),
     $     cdir, lo, hi, nc)

      implicit none
      integer nc, cdir
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
c
      integer i, j, n

      if(cdir .eq. 0) then
         do n = 1, nc
            do i = lo(1)-1, ARG_L1(phi), -1
               do j = lo(2)-1, ARG_L2(phi), -1
                  phi(i,j,n) = 0.5d0 *
     &                 ( (4*phi(i,j+1,n)-6*phi(i,j+2,n)+4*phi(i,j+3,n)-phi(i,j+4,n))
     &                 + (4*phi(i+1,j,n)-6*phi(i+2,j,n)+4*phi(i+3,j,n)-phi(i+4,j,n)) )
               enddo
               do j = hi(2)+1, ARG_H2(phi)
                  phi(i,j,n) = 0.5d0 *
     &                 ( (4*phi(i,j-1,n)-6*phi(i,j-2,n)+4*phi(i,j-3,n)-phi(i,j-4,n))
     &                 + (4*phi(i+1,j,n)-6*phi(i+2,j,n)+4*phi(i+3,j,n)-phi(i+4,j,n)) )
               enddo
            enddo

            do i = hi(1)+1, ARG_H1(phi)
               do j = lo(2)-1, ARG_L2(phi), -1
                  phi(i,j,n) = 0.5d0 *
     &                 ( (4*phi(i,j+1,n)-6*phi(i,j+2,n)+4*phi(i,j+3,n)-phi(i,j+4,n))
     &                 + (4*phi(i-1,j,n)-6*phi(i-2,j,n)+4*phi(i-3,j,n)-phi(i-4,j,n)) )
               enddo
               do j = hi(2)+1, ARG_H2(phi)
                  phi(i,j,n) = 0.5d0 *
     &                 ( (4*phi(i,j-1,n)-6*phi(i,j-2,n)+4*phi(i,j-3,n)-phi(i,j-4,n))
     &                 + (4*phi(i-1,j,n)-6*phi(i-2,j,n)+4*phi(i-3,j,n)-phi(i-4,j,n)) )
               enddo
            enddo
         enddo
      else
         print *,'In ABec4_2D::FORT_APPLYBC4_TOUCHUP: BAD cdir = ',cdir
         call bl_pd_abort()
      endif
      end

c-----------------------------------------------------------------------

      subroutine FORT_CA2CC(lo, hi, ca, DIMS(ca), cc, DIMS(cc), nc)

      implicit none
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(ca)
      integer DIMDEC(cc)
      integer nc
      REAL_T ca(DIMV(ca),nc)
      REAL_T cc(DIMV(cc),nc)

      integer i,j,n
      REAL_T one24th, seven6th
      parameter (one24th = 1.d0 / 24.d0)
      parameter (seven6th = 7.d0 / 6.d0)

      do n=1,nc
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,n) = -one24th*(
     &              ca(i,j-1,n)+ca(i-1,j,n)+ca(i+1,j,n)+ca(i,j+1,n))
     &              + seven6th*ca(i,j,n)
            enddo
         enddo
      enddo

      end

c-----------------------------------------------------------------------

      subroutine FORT_CC2CA(lo, hi, cc, DIMS(cc), ca, DIMS(ca), nc)

      implicit none
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(cc)
      integer DIMDEC(ca)
      integer nc
      REAL_T cc(DIMV(cc),nc)
      REAL_T ca(DIMV(ca),nc)

      integer i,j,n
      REAL_T one24th, five6th
      parameter (one24th = 1.d0 / 24.d0)
      parameter (five6th = 5.d0 / 6.d0)

      do n=1,nc
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ca(i,j,n) = one24th*(
     &              cc(i,j-1,n)+cc(i-1,j,n)+cc(i+1,j,n)+cc(i,j+1,n))
     &              + five6th*cc(i,j,n)
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------

      subroutine FORT_LO_CC2EC(lo, hi, 
     &     cfab, DIMS(cfab),
     &     efab, DIMS(efab), nc, dir,
     &     isharm)
      implicit none
      integer lo(2), hi(2), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.0.d0)then
                        efab(i,j,n)
     &                       = 2*(cfab(i,j,n) * cfab(i-1,j,n))/
     &                       (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=0.d0
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.0.d0)then
                        efab(i,j,n)
     &                       = 2*(cfab(i,j,n) * cfab(i,j-1,n))/
     &                       (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=0.d0
                     endif
                  end do
               end do
            end do
         end if
      end if
      end
