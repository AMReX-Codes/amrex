
module amrex_abec4_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  include 'AMReX_lo_bctypes.fi'

contains

#undef BC_ORDER_4
#define BC_ORDER_5

!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
      subroutine amrex_ab4_adotx( &
          y,y_l1,y_l2,y_l3,y_h1,y_h2,y_h3, &
          x,x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, &
          alpha, beta, &
          a,a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
          b,b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, &
          lo,hi,nc, &
          h &
          ) bind(c,name='amrex_ab4_adotx')
      implicit none

      real(amrex_real) alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer y_l1,y_l2,y_l3,y_h1,y_h2,y_h3
      integer x_l1,x_l2,x_l3,x_h1,x_h2,x_h3
      integer a_l1,a_l2,a_l3,a_h1,a_h2,a_h3
      integer b_l1,b_l2,b_l3,b_h1,b_h2,b_h3
      real(amrex_real)  y(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3,nc)
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,x_l3:x_h3,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2,a_l3:a_h3)
      real(amrex_real)  b(b_l1:b_h1,b_l2:b_h2,b_l3:b_h3)
      real(amrex_real) h(BL_SPACEDIM)

      integer i,j,k,n
      real(amrex_real) bh1i, bh2i, bh3i
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer zlo(BL_SPACEDIM), zhi(BL_SPACEDIM)

      real(amrex_real), allocatable :: xflux(:,:,:,:), yflux(:,:,:,:), zflux(:,:,:,:)

      xlo(:) = lo(:)
      xhi(:) = hi(:)
      ylo(:) = lo(:)
      yhi(:) = hi(:)
      zlo(:) = lo(:)
      zhi(:) = hi(:)
      xhi(1) = hi(1) + 1
      yhi(2) = hi(2) + 1
      zhi(3) = hi(3) + 1

      allocate(xflux(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),nc))
      allocate(yflux(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),nc))
      allocate(zflux(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3),nc))

      call flux_dir(x, x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, alpha, beta, &
           b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, nc, &
           h(1), xlo, xhi, xflux, xlo(1), xlo(2), xlo(3), xhi(1), xhi(2), xhi(3), 1)
      call flux_dir(x, x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, alpha, beta, &
           b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, nc, &
           h(2), ylo, yhi, yflux, ylo(1), ylo(2), ylo(3), yhi(1), yhi(2), yhi(3), 2)
      call flux_dir(x, x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, alpha, beta, &
           b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, nc, &
           h(3), zlo, zhi, zflux, zlo(1), zlo(2), zlo(3), zhi(1), zhi(2), zhi(3), 3)

      bh1i = beta / h(1)
      bh2i = beta / h(2)
      bh3i = beta / h(3)

      do n=1,nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) &
     &                 - bh1i*(xflux(i+1,j,k,n) - xflux(i,j,k,n)) &
     &                 - bh2i*(yflux(i,j+1,k,n) - yflux(i,j,k,n)) &
     &                 - bh3i*(zflux(i,j,k+1,n) - zflux(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      deallocate(xflux)
      deallocate(yflux)
      deallocate(zflux)

    end subroutine amrex_ab4_adotx

!-----------------------------------------------------------------------
!>
!>     Fill in fluxes
!>
      subroutine amrex_ab4_flux( &
          x,x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, &
          alpha, beta, &
          a, a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
          b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, &
          nc, &
          h, &
          xlo,xhi,xflux,xflux_l1,xflux_l2,xflux_l3,xflux_h1,xflux_h2,xflux_h3, &
          ylo,yhi,yflux,yflux_l1,yflux_l2,yflux_l3,yflux_h1,yflux_h2,yflux_h3, &
          zlo,zhi,zflux,zflux_l1,zflux_l2,zflux_l3,zflux_h1,zflux_h2,zflux_h3 &
          ) bind(c,name='amrex_ab4_flux')
      implicit none
      real(amrex_real) alpha, beta
      integer nc
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer zlo(BL_SPACEDIM), zhi(BL_SPACEDIM)
      integer x_l1,x_l2,x_l3,x_h1,x_h2,x_h3
      integer a_l1,a_l2,a_l3,a_h1,a_h2,a_h3
      integer b_l1,b_l2,b_l3,b_h1,b_h2,b_h3
      integer xflux_l1,xflux_l2,xflux_l3,xflux_h1,xflux_h2,xflux_h3
      integer yflux_l1,yflux_l2,yflux_l3,yflux_h1,yflux_h2,yflux_h3
      integer zflux_l1,zflux_l2,zflux_l3,zflux_h1,zflux_h2,zflux_h3
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,x_l3:x_h3,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2,a_l3:a_h3)
      real(amrex_real)  b(b_l1:b_h1,b_l2:b_h2,b_l3:b_h3)
      real(amrex_real) xflux(xflux_l1:xflux_h1,xflux_l2:xflux_h2,xflux_l3:xflux_h3,nc)
      real(amrex_real) yflux(yflux_l1:yflux_h1,yflux_l2:yflux_h2,yflux_l3:yflux_h3,nc)
      real(amrex_real) zflux(zflux_l1:zflux_h1,zflux_l2:zflux_h2,zflux_l3:zflux_h3,nc)
      real(amrex_real) h(BL_SPACEDIM)

      call flux_dir(x, x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, alpha, beta, &
           b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, nc, &
           h(1), xlo, xhi, xflux, xflux_l1,xflux_l2,xflux_l3,xflux_h1,xflux_h2,xflux_h3, 1)
      call flux_dir(x, x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, alpha, beta, &
           b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, nc, &
           h(2), ylo, yhi, yflux, yflux_l1,yflux_l2,yflux_l3,yflux_h1,yflux_h2,yflux_h3, 2)
      call flux_dir(x, x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, alpha, beta, &
           b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, nc, &
           h(3), zlo, zhi, zflux, zflux_l1,zflux_l2,zflux_l3,zflux_h1,zflux_h2,zflux_h3, 3)

    end subroutine amrex_ab4_flux

!-----------------------------------------------------------------------

      subroutine flux_dir( &
          x,x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, &
          alpha, beta, &
          b, b_l1,b_l2,b_l3,b_h1,b_h2,b_h3, &
          nc, &
          h, &
          lo,hi,flux,flux_l1,flux_l2,flux_l3,flux_h1,flux_h2,flux_h3, &
          dir)
      implicit none
      real(amrex_real) alpha, beta
      integer nc, dir
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer x_l1,x_l2,x_l3,x_h1,x_h2,x_h3
      integer b_l1,b_l2,b_l3,b_h1,b_h2,b_h3
      integer flux_l1,flux_l2,flux_l3,flux_h1,flux_h2,flux_h3
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,x_l3:x_h3,nc)
      real(amrex_real)  b(b_l1:b_h1,b_l2:b_h2,b_l3:b_h3)
      real(amrex_real) flux(flux_l1:flux_h1,flux_l2:flux_h2,flux_l3:flux_h3,nc)
      real(amrex_real) h

      integer i,j,k,n
      real(amrex_real) i12, i48
      real(amrex_real) hi1_12

      real(amrex_real), allocatable :: gt(:,:,:)
      real(amrex_real), allocatable :: bt(:,:,:)
      real(amrex_real) g1, g2, b1, b2

      i12 = 1.d0/12.d0
      i48 = 1.d0/48.d0
      hi1_12 = 1.d0 / (12.d0*h)

      if (dir .eq. 1) then
         allocate(gt(lo(1):hi(1),lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2))
         allocate(bt(lo(1):hi(1),lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2))

         do k = lo(3)-2, hi(3)+2
            do j = lo(2)-2, hi(2)+2
               do i = lo(1), hi(1)
                  bt(i,j,k) = (-b(i-2,j,k)+7*(b(i-1,j,k)+b(i,j,k))-b(i+1,j,k))*i12
               enddo
            enddo
         enddo
         do n=1,nc
            do k = lo(3)-2, hi(3)+2
               do j = lo(2)-2, hi(2)+2
                  do i = lo(1), hi(1)
                     gt(i,j,k) = (x(i-2,j,k,n)-x(i+1,j,k,n)+15.d0*(x(i,j,k,n)-x(i-1,j,k,n)))*hi1_12
                  enddo
               enddo
            enddo
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     g1 = (34.d0*(gt(i,j+1,k)-gt(i,j-1,k))+5.d0*(gt(i,j-2,k)-gt(i,j+2,k)))*i48
                     b1 = (34.d0*(bt(i,j+1,k)-bt(i,j-1,k))+5.d0*(bt(i,j-2,k)-bt(i,j+2,k)))*i48
                     g2 = (34.d0*(gt(i,j,k+1)-gt(i,j,k-1))+5.d0*(gt(i,j,k-2)-gt(i,j,k+2)))*i48
                     b2 = (34.d0*(bt(i,j,k+1)-bt(i,j,k-1))+5.d0*(bt(i,j,k-2)-bt(i,j,k+2)))*i48
                     flux(i,j,k,n) = bt(i,j,k)*gt(i,j,k)  +  i12*(b1*g1 + b2*g2)
                  enddo
               enddo
            enddo
         enddo
      else if (dir .eq. 2) then
         allocate(gt(lo(1)-2:hi(1)+2,lo(2):hi(2),lo(3)-2:hi(3)+2))
         allocate(bt(lo(1)-2:hi(1)+2,lo(2):hi(2),lo(3)-2:hi(3)+2))

         do k = lo(3)-2, hi(3)+2
            do j = lo(2), hi(2)
               do i = lo(1)-2, hi(1)+2
                  bt(i,j,k) = (-b(i,j-2,k)+7*(b(i,j-1,k)+b(i,j,k))-b(i,j+1,k))*i12
               enddo
            enddo
         enddo
         do n=1,nc
            do k = lo(3)-2, hi(3)+2
               do j = lo(2), hi(2)
                  do i = lo(1)-2, hi(1)+2
                     gt(i,j,k) = (x(i,j-2,k,n)-x(i,j+1,k,n)+15.d0*(x(i,j,k,n)-x(i,j-1,k,n)))*hi1_12
                  enddo
               enddo
            enddo
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     g1 = (34.d0*(gt(i+1,j,k)-gt(i-1,j,k))+5.d0*(gt(i-2,j,k)-gt(i+2,j,k)))*i48
                     b1 = (34.d0*(bt(i+1,j,k)-bt(i-1,j,k))+5.d0*(bt(i-2,j,k)-bt(i+2,j,k)))*i48
                     g2 = (34.d0*(gt(i,j,k+1)-gt(i,j,k-1))+5.d0*(gt(i,j,k-2)-gt(i,j,k+2)))*i48
                     b2 = (34.d0*(bt(i,j,k+1)-bt(i,j,k-1))+5.d0*(bt(i,j,k-2)-bt(i,j,k+2)))*i48
                     flux(i,j,k,n) = bt(i,j,k)*gt(i,j,k)  +  i12*(b1*g1 + b2*g2)
                  enddo
               enddo
            enddo
         enddo
      else
         allocate(gt(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3):hi(3)))
         allocate(bt(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3):hi(3)))

         do k = lo(3), hi(3)
            do j = lo(2)-2, hi(2)+2
               do i = lo(1)-2, hi(1)+2
                  bt(i,j,k) = (-b(i,j,k-2)+7*(b(i,j,k-1)+b(i,j,k))-b(i,j,k+1))*i12
               enddo
            enddo
         enddo
         do n=1,nc
            do k = lo(3), hi(3)
               do j = lo(2)-2, hi(2)+2
                  do i = lo(1)-2, hi(1)+2
                     gt(i,j,k) = (x(i,j,k-2,n)-x(i,j,k+1,n)+15.d0*(x(i,j,k,n)-x(i,j,k-1,n)))*hi1_12
                  enddo
               enddo
            enddo
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     g1 = (34.d0*(gt(i+1,j,k)-gt(i-1,j,k))+5.d0*(gt(i-2,j,k)-gt(i+2,j,k)))*i48
                     b1 = (34.d0*(bt(i+1,j,k)-bt(i-1,j,k))+5.d0*(bt(i-2,j,k)-bt(i+2,j,k)))*i48
                     g2 = (34.d0*(gt(i,j+1,k)-gt(i,j-1,k))+5.d0*(gt(i,j-2,k)-gt(i,j+2,k)))*i48
                     b2 = (34.d0*(bt(i,j+1,k)-bt(i,j-1,k))+5.d0*(bt(i,j-2,k)-bt(i,j+2,k)))*i48
                     flux(i,j,k,n) = bt(i,j,k)*gt(i,j,k)  +  i12*(b1*g1 + b2*g2)
                  enddo
               enddo
            enddo
         enddo
      endif

      deallocate(gt)
      deallocate(bt)

    end subroutine flux_dir

!-----------------------------------------------------------------------
      subroutine amrex_ab4_applybc4 ( &
          flagden, flagbc, maxorder, &
          phi,   phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
          cdir, bct, bcl, &
          bcval, bcval_l1,bcval_l2,bcval_l3,bcval_h1,bcval_h2,bcval_h3, &
          mask,  mask_l1,mask_l2,mask_l3,mask_h1,mask_h2,mask_h3, &
          den,   den_l1,den_l2,den_l3,den_h1,den_h2,den_h3, &
          lo, hi, nc, &
          h &
          ) bind(c,name='amrex_ab4_applybc4')

      implicit none
!
!     If the boundary is of Neumann type, set the ghost cell value to
!     that of the outermost point in the valid data (2nd order accurate)
!     and then fill the "den" array with the value "1"
!
!
!     If flagbc==1:
!
!     If the boundary is of Dirichlet type, construct a polynomial
!     interpolation through the boundary location and internal points
!     (at locations x(-1:len-2) that generates the ghost cell value (at
!     location xInt).  Then fill the ghost cell with the interpolated value.
!     If flagden==1, load the "den" array with the interpolation
!     coefficient corresponding to outermost point in the valid region
!     ( the coef(0) corresponding to the location x(0) )
!
!     Note:
!     The bc type = LO_REFLECT_ODD is a special type of boundary condition.
!
      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer den_l1,den_l2,den_l3,den_h1,den_h2,den_h3
      integer bcval_l1,bcval_l2,bcval_l3,bcval_h1,bcval_h2,bcval_h3
      integer mask_l1,mask_l2,mask_l3,mask_h1,mask_h2,mask_h3
      real(amrex_real)  phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3,nc)
      real(amrex_real)  den(den_l1:den_h1,den_l2:den_h2,den_l3:den_h3)
      real(amrex_real)  bcval(bcval_l1:bcval_h1,bcval_l2:bcval_h2,bcval_l3:bcval_h3,nc)
      integer mask(mask_l1:mask_h1,mask_l2:mask_h2,mask_l3:mask_h3)
      integer bct
      real(amrex_real) bcl, tmp, b, d
      parameter (b = 1.d0/24.d0)
      parameter (d = 5.d0/6.d0)
      real(amrex_real)  h(BL_SPACEDIM)
      integer i, j, k, n
      logical is_dirichlet, is_neumann

      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )

!
!
!     The Left face of the grid
!
      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     if (mask(lo(1)-1,j,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(lo(1)-1,j,k,n) = (9*phi(lo(1),j,k,n)+3*phi(lo(1)+1,j,k,n) &
     &                       -phi(lo(1)+2,j,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(lo(1)-1,j,k,n) = (5*phi(lo(1),j,k,n)+9*phi(lo(1)+1,j,k,n) &
     &                       -5*phi(lo(1)+2,j,k,n)+phi(lo(1)+3,j,k,n))/10.d0
#endif
                     endif
                     if (mask(lo(1)-2,j,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(lo(1)-2,j,k,n) = (-30*phi(lo(1),j,k,n)+56*phi(lo(1)+1,j,k,n) &
     &                       -15*phi(lo(1)+2,j,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(lo(1)-2,j,k,n) = (-75*phi(lo(1),j,k,n)+145*phi(lo(1)+1,j,k,n) &
     &                       -75*phi(lo(1)+2,j,k,n)+15*phi(lo(1)+3,j,k,n))/10.d0
#endif
                     endif
                  enddo
               enddo
            enddo
            if (flagden .eq. 1) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(lo(1),j,k) = 1.0D0
                  enddo
               enddo
            endif
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        if (mask(lo(1)-1,j,k) .gt. 0) then
!     Convert ec BC data to ea BC data -> tmp
                           if (mask(lo(1)-1,j-1,k).gt.0 .and. mask(lo(1)-1,j+1,k).gt.0 &
                               .and. mask(lo(1)-1,j,k-1).gt.0 .and. mask(lo(1)-1,j,k+1).gt.0) then
                              tmp = b*(bcval(lo(1)-1,j+1,k,n)+bcval(lo(1)-1,j-1,k,n) &
     &                             +bcval(lo(1)-1,j,k+1,n)+bcval(lo(1)-1,j,k-1,n))+d*bcval(lo(1)-1,j,k,n)
                           else
                              tmp = bcval(lo(1)-1,j,k,n)
                           endif
#if defined(BC_ORDER_4)
                           phi(lo(1)-1,j,k,n) = (12*tmp-13*phi(lo(1),j,k,n) &
     &                          +5*phi(lo(1)+1,j,k,n)-phi(lo(1)+2,j,k,n))/3.d0
#elif defined(BC_ORDER_5)
                           phi(lo(1)-1,j,k,n) = (60*tmp-77*phi(lo(1),j,k,n)+43*phi(lo(1)+1,j,k,n) &
     &                          -17*phi(lo(1)+2,j,k,n)+3*phi(lo(1)+3,j,k,n))/12.d0
#endif
                           if (mask(lo(1)-2,j,k) .gt. 0) then
#if defined(BC_ORDER_4)
                             phi(lo(1)-2,j,k,n) = (48*tmp - 70*phi(lo(1),j,k,n) + &
     &                             32*phi(lo(1)+1,j,k,n) - 7*phi(lo(1)+2,j,k,n))/3.d0
#elif defined(BC_ORDER_5)
                              phi(lo(1)-2,j,k,n) = (300*tmp-505*phi(lo(1),j,k,n)+ &
     &                            335*phi(lo(1)+1,j,k,n)-145*phi(lo(1)+2,j,k,n)+27*phi(lo(1)+3,j,k,n))/12.d0
#endif
                           endif
                        endif
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        if (mask(lo(1)-1,j,k) .gt. 0) then
                           phi(lo(1)-1,j,k,n) = 0.d0
                        endif
                        if (mask(lo(1)-2,j,k) .gt. 0) then
                           phi(lo(1)-2,j,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(lo(1),j,k) = merge(66.D0, 0.0D0, &
                         mask(lo(1)-1,j,k) .gt. 0)
                  enddo
               enddo
            endif

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     phi(lo(1)-1,j,k,n) = merge( &
                         -phi(lo(1),j,k,n), &
                         phi(lo(1)-1,j,k,n), &
                         mask(lo(1)-1,j,k) .gt. 0)
                     phi(lo(1)-2,j,k,n) = merge( &
                         -phi(lo(1)+1,j,k,n), &
                         phi(lo(1)-2,j,k,n), &
                         mask(lo(1)-2,j,k) .gt. 0)
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(lo(1),j,k) = merge(-1.0D0, 0.0D0, &
                         mask(lo(1)-1,j,k) .gt. 0)
                  enddo
               enddo
            endif
         else
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The Right face of the grid
!
      if(cdir .eq. 3) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     if (mask(hi(1)+1,j,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(hi(1)+1,j,k,n) =  (9*phi(hi(1),j,k,n) &
     &                       +3*phi(hi(1)-1,j,k,n)-phi(hi(1)-2,j,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(hi(1)+1,j,k,n) =  (5*phi(hi(1)  ,j,k,n) +  9*phi(hi(1)-1,j,k,n) &
     &                       - 5*phi(hi(1)-2,j,k,n) + phi(hi(1)-3,j,k,n))/10.d0
#endif
                     endif
                     if (mask(hi(1)+2,j,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(hi(1)+2,j,k,n) = (- 30*phi(hi(1),j,k,n) + 56*phi(hi(1)-1,j,k,n) &
     &                       - 15*phi(hi(1)-2,j,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(hi(1)+2,j,k,n) = (- 75*phi(hi(1),j,k,n) + 145*phi(hi(1)-1,j,k,n) &
     &                       - 75*phi(hi(1)-2,j,k,n) + 15*phi(hi(1)-3,j,k,n))/10.d0
#endif
                     endif
                  enddo
               enddo
            enddo
	    if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(hi(1),j,k) = 1.0D0
                  enddo
               enddo
	    endif
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        if (mask(hi(1)+1,j,k) .gt. 0) then
!     Convert ec BC data to ea BC data -> tmp
                           if (mask(hi(1)+1,j-1,k).gt.0 .and. mask(hi(1)+1,j+1,k).gt.0 &
                               .and. mask(hi(1)+1,j,k-1).gt.0 .and. mask(hi(1)+1,j,k+1).gt.0) then
                              tmp = b*(bcval(hi(1)+1,j+1,k,n)+bcval(hi(1)+1,j-1,k,n)+bcval(hi(1)+1,j,k+1,n) &
     &                             +bcval(hi(1)+1,j,k-1,n)) + d*bcval(hi(1)+1,j,k,n)
                           else
                              tmp = bcval(hi(1)+1,j,k,n)
                           endif
#if defined(BC_ORDER_4)
                           phi(hi(1)+1,j,k,n) = (12*tmp - 13*phi(hi(1),j,k,n) &
     &                          + 5*phi(hi(1)-1,j,k,n) - phi(hi(1)-2,j,k,n))/3.d0
#elif defined(BC_ORDER_5)
                           phi(hi(1)+1,j,k,n) = (60*tmp - 77*phi(hi(1),j,k,n) + 43*phi(hi(1)-1,j,k,n) &
     &                          - 17*phi(hi(1)-2,j,k,n) + 3*phi(hi(1)-3,j,k,n))/12.d0
#endif
                           if (mask(hi(1)+2,j,k) .gt. 0) then
#if defined(BC_ORDER_4)
                              phi(hi(1)+2,j,k,n) = (48*tmp - 70*phi(hi(1),j,k,n) &
     &                             + 32*phi(hi(1)-1,j,k,n) - 7*phi(hi(1)-2,j,k,n))/3.d0
#elif defined(BC_ORDER_5)
                              phi(hi(1)+2,j,k,n) = (300*tmp - 505*phi(hi(1),j,k,n) + 335*phi(hi(1)-1,j,k,n) &
     &                             - 145*phi(hi(1)-2,j,k,n) + 27*phi(hi(1)-3,j,k,n))/12.d0
#endif
                           endif
                        endif
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        if (mask(hi(1)+1,j,k) .gt. 0) then
                           phi(hi(1)+1,j,k,n) = 0.d0
                        endif
                        if (mask(hi(1)+2,j,k) .gt. 0) then
                           phi(hi(1)+2,j,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(hi(1),j,k)   = merge(66.D0, 0.0D0, &
                         mask(hi(1)+1,j,k) .gt. 0)
                  enddo
               enddo
            endif

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,k,n) = merge( &
                         -phi(hi(1),j,k,n), &
                         phi(hi(1)+1,j,k,n), &
                         mask(hi(1)+1,j,k) .gt. 0)
                     phi(hi(1)+2,j,k,n) = merge( &
                         -phi(hi(1)-1,j,k,n), &
                         phi(hi(1)+2,j,k,n), &
                         mask(hi(1)+2,j,k) .gt. 0)
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(hi(1),j,k) = merge(-1.0D0, 0.0D0, &
                         mask(hi(1)+1,j,k) .gt. 0)
                  enddo
               enddo
            endif
         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         endif
      endif
!
!     The Bottom of the Grid
!
      if(cdir .eq. 1) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do k = lo(3),hi(3)
                  do i = lo(1),hi(1)
                     if (mask(i,lo(2)-1,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,lo(2)-1,k,n) =  (9*phi(i,lo(2),k,n) +  3*phi(i,lo(2)+1,k,n) &
     &                       - phi(i,lo(2)+2,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,lo(2)-1,k,n) =  (5*phi(i,lo(2)  ,k,n) +  9*phi(i,lo(2)+1,k,n) &
     &                       - 5*phi(i,lo(2)+2,k,n) + phi(i,lo(2)+3,k,n))/10.d0
#endif
                     endif
                     if (mask(i,lo(2)-2,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,lo(2)-2,k,n) = (- 30*phi(i,lo(2),k,n) + 56*phi(i,lo(2)+1,k,n) &
     &                       - 15*phi(i,lo(2)+2,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,lo(2)-2,k,n) = (- 75*phi(i,lo(2),k,n) + 145*phi(i,lo(2)+1,k,n) &
     &                       - 75*phi(i,lo(2)+2,k,n) + 15*phi(i,lo(2)+3,k,n))/10.d0
#endif
                     endif
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3),hi(3)
                  do i = lo(1),hi(1)
                     den(i,lo(2),k)   = 1.0D0
                  enddo
               enddo
            endif
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        if (mask(i,lo(2)-1,k) .gt. 0) then
!     Convert ec BC data to ea BC data -> tmp
                           if (mask(i-1,lo(2)-1,k).gt.0 .and. mask(i+1,lo(2)-1,k).gt.0 &
                               .and. mask(i,lo(2)-1,k-1).gt.0 .and. mask(i,lo(2)-1,k+1).gt.0) then
                              tmp = b*(bcval(i+1,lo(2)-1,k,n)+bcval(i-1,lo(2)-1,k,n)+bcval(i,lo(2)-1,k+1,n) &
     &                             +bcval(i,lo(2)-1,k-1,n)) + d*bcval(i,lo(2)-1,k,n)
                           else
                              tmp = bcval(i,lo(2)-1,k,n)
                           endif
#if defined(BC_ORDER_4)
                           phi(i,lo(2)-1,k,n) = (12*tmp - 13*phi(i,lo(2),k,n) + 5*phi(i,lo(2)+1,k,n) &
     &                          - phi(i,lo(2)+2,k,n))/3.d0
#elif defined(BC_ORDER_5)
                           phi(i,lo(2)-1,k,n) = (60*tmp - 77*phi(i,lo(2),k,n) + 43*phi(i,lo(2)+1,k,n) &
     &                          - 17*phi(i,lo(2)+2,k,n) + 3*phi(i,lo(2)+3,k,n))/12.d0
#endif
                           if (mask(i,lo(2)-2,k) .gt. 0) then
#if defined(BC_ORDER_4)
                              phi(i,lo(2)-2,k,n) =   (48*tmp - 70*phi(i,lo(2),k,n) + 32*phi(i,lo(2)+1,k,n) &
     &                             - 7*phi(i,lo(2)+2,k,n))/3.d0
#elif defined(BC_ORDER_5)
                              phi(i,lo(2)-2,k,n) = (300*tmp - 505*phi(i,lo(2),k,n) + 335*phi(i,lo(2)+1,k,n) &
     &                             - 145*phi(i,lo(2)+2,k,n) + 27*phi(i,lo(2)+3,k,n))/12.d0
#endif
                           endif
                        endif
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        if (mask(i,lo(2)-1,k) .gt. 0) then
                           phi(i,lo(2)-1,k,n) = 0.d0
                        endif
                        if (mask(i,lo(2)-2,k) .gt. 0) then
                           phi(i,lo(2)-2,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,lo(2),k)   = merge(66.D0, 0.0D0, &
                         mask(i,lo(2)-1,k) .gt. 0)
                  enddo
               enddo
            endif

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     phi(i,lo(2)-1,k,n) = merge( &
                         -phi(i,lo(2),k,n), &
                         phi(i,lo(2)-1,k,n), &
                         mask(i,lo(2)-1,k) .gt. 0)
                     phi(i,lo(2)-2,k,n) = merge( &
                         -phi(i,lo(2)+1,k,n), &
                         phi(i,lo(2)-2,k,n), &
                         mask(i,lo(2)-2,k) .gt. 0)
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,lo(2),k) = merge(-1.0D0, 0.0D0, &
                         mask(i,lo(2)-1,k) .gt. 0)
                  enddo
               enddo
            endif

         else
            print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The top of the grid
!
      if (cdir .eq. 4) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     if (mask(i,hi(2)+1,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,hi(2)+1,k,n) =  (9*phi(i,hi(2),k,n) +  3*phi(i,hi(2)-1,k,n) &
     &                       - phi(i,hi(2)-2,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,hi(2)+1,k,n) =  (5*phi(i,hi(2),k,n) +  9*phi(i,hi(2)-1,k,n) &
     &                       - 5*phi(i,hi(2)-2,k,n) + phi(i,hi(2)-3,k,n))/10.d0
#endif
                     endif
                     if (mask(i,hi(2)+2,k) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,hi(2)+2,k,n) = (- 30*phi(i,hi(2),k,n) + 56*phi(i,hi(2)-1,k,n) &
     &                       - 15*phi(i,hi(2)-2,k,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,hi(2)+2,k,n) = (- 75*phi(i,hi(2),k,n) + 145*phi(i,hi(2)-1,k,n) &
     &                       - 75*phi(i,hi(2)-2,k,n) + 15*phi(i,hi(2)-3,k,n))/10.d0
#endif
                     endif
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,hi(2),k)   = 1.0D0
                  enddo
               enddo
            endif
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        if (mask(i,hi(2)+1,k) .gt. 0) then
!     Convert ec BC data to ea BC data -> tmp
                           if (mask(i-1,hi(2)+1,k).gt.0 .and. mask(i+1,hi(2)+1,k).gt.0 &
                               .and. mask(i,hi(2)+1,k-1).gt.0 .and. mask(i,hi(2)+1,k+1).gt.0) then
                              tmp = b*(bcval(i+1,hi(2)+1,k,n)+bcval(i-1,hi(2)+1,k,n)+bcval(i,hi(2)+1,k+1,n) &
     &                             +bcval(i,hi(2)+1,k-1,n)) + d*bcval(i,hi(2)+1,k,n)
                           else
                              tmp = bcval(i,hi(2)+1,k,n)
                           endif
#if defined(BC_ORDER_4)
                           phi(i,hi(2)+1,k,n) = (12*tmp - 13*phi(i,hi(2),k,n) + 5*phi(i,hi(2)-1,k,n) &
     &                          - phi(i,hi(2)-2,k,n))/3.d0
#elif defined(BC_ORDER_5)
                           phi(i,hi(2)+1,k,n) = (60*tmp - 77*phi(i,hi(2),k,n) + 43*phi(i,hi(2)-1,k,n) &
     &                          - 17*phi(i,hi(2)-2,k,n) + 3*phi(i,hi(2)-3,k,n))/12.d0
#endif
                           if (mask(i,hi(2)+2,k) .gt. 0) then
#if defined(BC_ORDER_4)
                              phi(i,hi(2)+2,k,n) =   (48*tmp -  70*phi(i,hi(2),k,n) + 32*phi(i,hi(2)-1,k,n) &
     &                             - 7*phi(i,hi(2)-2,k,n))/3.d0
#elif defined(BC_ORDER_5)
                              phi(i,hi(2)+2,k,n) = (300*tmp - 505*phi(i,hi(2),k,n) + 335*phi(i,hi(2)-1,k,n) &
     &                             - 145*phi(i,hi(2)-2,k,n) + 27*phi(i,hi(2)-3,k,n))/12.d0
#endif
                           endif
                        endif
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        if (mask(i,hi(2)+1,k) .gt. 0) then
                           phi(i,hi(2)+1,k,n) = 0.d0
                        endif
                        if (mask(i,hi(2)+2,k) .gt. 0) then
                           phi(i,hi(2)+2,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,hi(2),k)   = merge(66.D0, 0.0D0, &
                         mask(i,hi(2)+1,k) .gt. 0)
                  enddo
               enddo
            endif

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     phi(i,hi(2)+1,k,n) = merge( &
                         -phi(i,hi(2),k,n), &
                         phi(i,hi(2)+1,k,n), &
                         mask(i,hi(2)+1,k) .gt. 0)
                     phi(i,hi(2)+2,k,n) = merge( &
                         -phi(i,hi(2)-2,k,n), &
                         phi(i,hi(2)+2,k,n), &
                         mask(i,hi(2)+2,k) .gt. 0)
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,hi(2),k) = merge(-1.0D0, 0.0D0, &
                         mask(i,hi(2)+1,k) .gt. 0)
                  enddo
               enddo
            endif

         else
            print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
            call bl_error("stop")
         end if
      endif
!
!
!     The Front of the Grid
!
      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     if (mask(i,j,lo(3)-1) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,j,lo(3)-1,n) =  (9*phi(i,j,lo(3),n) +  3*phi(i,j,lo(3)+1,n) &
     &                       - phi(i,j,lo(3)+2,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,j,lo(3)-1,n) =  (5*phi(i,j,lo(3),n) +  9*phi(i,j,lo(3)+1,n) &
     &                       - 5*phi(i,j,lo(3)+2,n) + phi(i,j,lo(3)+3,n))/10.d0
#endif
                     endif
                     if (mask(i,k,lo(3)-2) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,j,lo(3)-2,n) = (- 30*phi(i,j,lo(3),n) + 56*phi(i,j,lo(3)+1,n) &
     &                       - 15*phi(i,j,lo(3)+2,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,j,lo(3)-2,n) = (- 75*phi(i,j,lo(3),n) + 145*phi(i,j,lo(3)+1,n) &
     &                       - 75*phi(i,j,lo(3)+2,n) + 15*phi(i,j,lo(3)+3,n))/10.d0
#endif
                     endif
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     den(i,j,lo(3))   = 1.0D0
                  enddo
               enddo
            endif
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if (mask(i,j,lo(3)-1) .gt. 0) then
!     Convert ec BC data to ea BC data -> tmp
                           if (mask(i-1,j,lo(3)-1).gt.0 .and. mask(i+1,j,lo(3)-1).gt.0 &
                               .and. mask(i,j-1,lo(3)-1).gt.0 .and. mask(i,j+1,lo(3)-1).gt.0) then
                              tmp = b*(bcval(i+1,j,lo(3)-1,n)+bcval(i-1,j,lo(3)-1,n)+bcval(i,j+1,lo(3)-1,n) &
     &                             +bcval(i,j-1,lo(3)-1,n)) + d*bcval(i,j,lo(3)-1,n)
                           else
                              tmp = bcval(i,j,lo(3)-1,n)
                           endif
#if defined(BC_ORDER_4)
                           phi(i,j,lo(3)-1,n) = (12*tmp - 13*phi(i,j,lo(3),n) + 5*phi(i,j,lo(3)+1,n) &
     &                          - phi(i,j,lo(3)+2,n))/3.d0
#elif defined(BC_ORDER_5)
                           phi(i,j,lo(3)-1,n) = (60*tmp - 77*phi(i,j,lo(3),n) + 43*phi(i,j,lo(3)+1,n) &
     &                          - 17*phi(i,j,lo(3)+2,n) + 3*phi(i,j,lo(3)+3,n))/12.d0
#endif
                           if (mask(i,j,lo(3)-2) .gt. 0) then
#if defined(BC_ORDER_4)
                              phi(i,j,lo(3)-2,n) =   (48*tmp - 70*phi(i,j,lo(3),n) + 32*phi(i,j,lo(3)+1,n) &
     &                             - 7*phi(i,j,lo(3)+2,n))/3.d0
#elif defined(BC_ORDER_5)
                              phi(i,j,lo(3)-2,n) = (300*tmp - 505*phi(i,j,lo(3),n) + 335*phi(i,j,lo(3)+1,n) &
     &                             - 145*phi(i,j,lo(3)+2,n) + 27*phi(i,j,lo(3)+3,n))/12.d0
#endif
                           endif
                        endif
                     enddo
                  enddo
               else
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if (mask(i,j,lo(3)-1) .gt. 0) then
                           phi(i,j,lo(3)-1,n) = 0.d0
                        endif
                        if (mask(i,j,lo(3)-2) .gt. 0) then
                           phi(i,j,lo(3)-2,n) = 0.d0
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,lo(3))   = merge(66.D0, 0.0D0, &
                         mask(i,j,lo(3)-1) .gt. 0)
                  enddo
               enddo
            endif

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     phi(i,j,lo(3)-1,n) = merge( &
                         -phi(i,j,lo(3),n), &
                         phi(i,j,lo(3)-1,n), &
                         mask(i,j,lo(3)-1) .gt. 0)
                     phi(i,j,lo(3)-2,n) = merge( &
                         -phi(i,j,lo(3)+1,n), &
                         phi(i,j,lo(3)-2,n), &
                         mask(i,j,lo(3)-2) .gt. 0)
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,lo(3)) = merge(-1.0D0, 0.0D0, &
                         mask(i,j,lo(3)-1) .gt. 0)
                  enddo
               enddo
            endif

         else
            print *,'UNKNOWN BC ON FRONT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The back of the grid
!
      if (cdir .eq. 5) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if (mask(i,j,hi(2)+1) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,j,hi(3)+1,n) =  (9*phi(i,j,hi(3),n) +  3*phi(i,j,hi(3)-1,n) &
     &                       - phi(i,j,hi(3)-2,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,j,hi(3)+1,n) =  (5*phi(i,j,hi(3),n) +  9*phi(i,j,hi(3)-1,n) &
     &                       - 5*phi(i,j,hi(3)-2,n) + phi(i,j,hi(3)-3,n))/10.d0
#endif
                     endif
                     if (mask(i,j,hi(3)+2) .gt. 0) then
#if defined(BC_ORDER_4)
                        phi(i,j,hi(3)+2,n) = (- 30*phi(i,j,hi(3),n) + 56*phi(i,j,hi(3)-1,n) &
     &                       - 15*phi(i,j,hi(3)-2,n))/11.d0
#elif defined(BC_ORDER_5)
                        phi(i,j,hi(3)+2,n) = (- 75*phi(i,j,hi(3),n) + 145*phi(i,j,hi(3)-1,n) &
     &                       - 75*phi(i,j,hi(3)-2,n) + 15*phi(i,j,hi(3)-3,n))/10.d0
#endif
                     endif
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,hi(3))   = 1.0D0
                  enddo
               enddo
            endif
         else if (is_dirichlet(bct)) then
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if (mask(i,j,hi(3)+1) .gt. 0) then
!     Convert ec BC data to ea BC data -> tmp
                           if (mask(i-1,j,hi(3)+1).gt.0 .and. mask(i+1,j,hi(3)+1).gt.0 &
                               .and. mask(i,j-1,hi(3)+1).gt.0 .and. mask(i,j+1,hi(3)+1).gt.0) then
                              tmp = b*(bcval(i+1,j,hi(3)+1,n)+bcval(i-1,j,hi(3)+1,n)+bcval(i,j+1,hi(3)+1,n) &
     &                             +bcval(i,j-1,hi(3)+1,n)) + d*bcval(i,j,hi(3)+1,n)
                           else
                              tmp = bcval(i,j,hi(3)+1,n)
                           endif
#if defined(BC_ORDER_4)
                           phi(i,j,hi(3)+1,n) = (12*tmp - 13*phi(i,j,hi(3),n) + 5*phi(i,j,hi(3)-1,n) &
     &                          - phi(i,j,hi(3)-2,n))/3.d0
#elif defined(BC_ORDER_5)
                           phi(i,j,hi(3)+1,n) = (60*tmp - 77*phi(i,j,hi(3),n) + 43*phi(i,j,hi(3)-1,n) &
     &                          - 17*phi(i,j,hi(3)-2,n) + 3*phi(i,j,hi(3)-3,n))/12.d0
#endif
                           if (mask(i,j,hi(3)+2) .gt. 0) then
#if defined(BC_ORDER_4)
                              phi(i,j,hi(3)+2,n) =   (48*tmp -  70*phi(i,j,hi(3),n) + 32*phi(i,j,hi(3)-1,n) &
     &                             - 7*phi(i,j,hi(3)-2,n))/3.d0
#elif defined(BC_ORDER_5)
                              phi(i,j,hi(3)+2,n) = (300*tmp - 505*phi(i,j,hi(3),n) + 335*phi(i,j,hi(3)-1,n) &
     &                             - 145*phi(i,j,hi(3)-2,n) + 27*phi(i,j,hi(3)-3,n))/12.d0
#endif
                           endif
                        endif
                     enddo
                  enddo
               else
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if (mask(i,j,hi(3)+1) .gt. 0) then
                           phi(i,j,hi(3)+1,n) = 0.d0
                        endif
                        if (mask(i,j,hi(3)+2) .gt. 0) then
                           phi(i,j,hi(3)+2,n) = 0.d0
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,hi(3))   = merge(66.D0, 0.0D0, &
                         mask(i,j,hi(3)+1) .gt. 0)
                  enddo
               enddo
            endif

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     phi(i,j,hi(3)+1,n) = merge( &
                         -phi(i,j,hi(3),n), &
                         phi(i,j,hi(3)+1,n), &
                         mask(i,j,hi(3)+1) .gt. 0)
                     phi(i,j,hi(3)+2,n) = merge( &
                         -phi(i,j,hi(3)-2,n), &
                         phi(i,j,hi(3)+2,n), &
                         mask(i,j,hi(3)+2) .gt. 0)
                  enddo
               enddo
            enddo
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,hi(3)) = merge(-1.0D0, 0.0D0, &
                         mask(i,j,hi(3)+1) .gt. 0)
                  enddo
               enddo
            endif

         else
            print *,'UNKNOWN BC ON BACK FACE IN APPLYBC'
            call bl_error("stop")
         end if
      endif
!
      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_applybc4_touchup ( &
          phi,   phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
          lo, hi, nc) bind(c,name='amrex_ab4_applybc4_touchup')

      implicit none
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      real(amrex_real) phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3,nc)
!
      integer i, j, k, n

      do n = 1, nc
         do k = lo(3),hi(3)
            do i = lo(1)-1, phi_l1, -1
               do j = lo(2)-1, phi_l2, -1
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) )
               enddo
               do j = hi(2)+1, phi_h2
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) )
               enddo
            enddo

            do i = hi(1)+1, phi_h1
               do j = lo(2)-1, phi_l2, -1
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) )
               enddo
               do j = hi(2)+1, phi_h2
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) )
               enddo
            enddo

         enddo

         do j = lo(2),hi(2)
            do i = lo(1)-1, phi_l1, -1
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo

            do i = hi(1)+1, phi_h1
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo
         enddo


         do i = lo(1),hi(1)
            do j = lo(2)-1, phi_l2, -1
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo

            do j = hi(2)+1, phi_h2
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = half * &
     &                 ( (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo
         enddo


         do i = lo(1)-1, phi_l1, -1
            do j = lo(2)-1, phi_l2, -1
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo
            do j = hi(2)+1, phi_h2
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i+1,j,k,n)-6*phi(i+2,j,k,n)+4*phi(i+3,j,k,n)-phi(i+4,j,k,n)) &
     &                 + (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo
         enddo

         do i = hi(1)+1, phi_h1
            do j = lo(2)-1, phi_l2, -1
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j+1,k,n)-6*phi(i,j+2,k,n)+4*phi(i,j+3,k,n)-phi(i,j+4,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo
            do j = hi(2)+1, phi_h2
               do k = lo(3)-1, phi_l3, -1
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) &
     &                 + (4*phi(i,j,k+1,n)-6*phi(i,j,k+2,n)+4*phi(i,j,k+3,n)-phi(i,j,k+4,n)) )
               enddo
               do k = hi(3)+1, phi_h3
                  phi(i,j,k,n) = third * &
     &                 ( (4*phi(i-1,j,k,n)-6*phi(i-2,j,k,n)+4*phi(i-3,j,k,n)-phi(i-4,j,k,n)) &
     &                 + (4*phi(i,j-1,k,n)-6*phi(i,j-2,k,n)+4*phi(i,j-3,k,n)-phi(i,j-4,k,n)) &
     &                 + (4*phi(i,j,k-1,n)-6*phi(i,j,k-2,n)+4*phi(i,j,k-3,n)-phi(i,j,k-4,n)) )
               enddo
            enddo
         enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_ca2cc(lo, hi, ca, ca_l1,ca_l2,ca_l3,ca_h1,ca_h2,ca_h3, &
           cc, cc_l1,cc_l2,cc_l3,cc_h1,cc_h2,cc_h3, nc) bind(c,name='amrex_ab4_ca2cc')

      implicit none
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer ca_l1,ca_l2,ca_l3,ca_h1,ca_h2,ca_h3
      integer cc_l1,cc_l2,cc_l3,cc_h1,cc_h2,cc_h3
      integer nc
      real(amrex_real) ca(ca_l1:ca_h1,ca_l2:ca_h2,ca_l3:ca_h3,nc)
      real(amrex_real) cc(cc_l1:cc_h1,cc_l2:cc_h2,cc_l3:cc_h3,nc)

      integer i,j,k,n
      real(amrex_real) b, d
      parameter (b = -1.d0/24.d0)
      parameter (d = 1.25d0)

      do n=1,nc
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  cc(i,j,k,n) = b*( &
     &                  ca(i-1,j,k,n)+ca(i,j-1,k,n)+ca(i,j,k-1,n) &
     &                 +ca(i+1,j,k,n)+ca(i,j+1,k,n)+ca(i,j,k+1,n)) &
     &                 + d*ca(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_cc2ca(lo, hi, cc, cc_l1,cc_l2,cc_l3,cc_h1,cc_h2,cc_h3, &
           ca, ca_l1,ca_l2,ca_l3,ca_h1,ca_h2,ca_h3, nc) bind(c,name='amrex_ab4_cc2ca')

      implicit none
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cc_l1,cc_l2,cc_l3,cc_h1,cc_h2,cc_h3
      integer ca_l1,ca_l2,ca_l3,ca_h1,ca_h2,ca_h3
      integer nc
      real(amrex_real) cc(cc_l1:cc_h1,cc_l2:cc_h2,cc_l3:cc_h3,nc)
      real(amrex_real) ca(ca_l1:ca_h1,ca_l2:ca_h2,ca_l3:ca_h3,nc)

      integer i,j,k,n
      real(amrex_real) b, d
      parameter (b = 1.d0/24.d0)
      parameter (d = 0.75d0)

      do n=1,nc
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  ca(i,j,k,n) = b*( &
     &                  cc(i-1,j,k,n)+cc(i,j-1,k,n)+cc(i,j,k-1,n) &
     &                 +cc(i+1,j,k,n)+cc(i,j+1,k,n)+cc(i,j,k+1,n)) &
     &                 + d*cc(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_lo_cc2ec(lo, hi, &
     &     cfab, cfab_l1,cfab_l2,cfab_l3,cfab_h1,cfab_h2,cfab_h3, &
     &     efab, efab_l1,efab_l2,efab_l3,efab_h1,efab_h2,efab_h3, &
     &     nc, dir, &
     &     isharm) bind(c,name='amrex_ab4_lo_cc2ec')
      implicit none
      integer lo(3), hi(3), nc, dir, isharm
      integer cfab_l1,cfab_l2,cfab_l3,cfab_h1,cfab_h2,cfab_h3
      integer efab_l1,efab_l2,efab_l3,efab_h1,efab_h2,efab_h3
      real(amrex_real)  cfab(cfab_l1:cfab_h1,cfab_l2:cfab_h2,cfab_l3:cfab_h3, nc)
      real(amrex_real)  efab(efab_l1:efab_h1,efab_l2:efab_h2,efab_l3:efab_h3, nc)

      integer i,j,k,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     enddo
                  enddo
               enddo
            enddo
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     enddo
                  enddo
               enddo
            enddo
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                     enddo
                  enddo
               enddo
            enddo
         endif
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)).ne.0.d0) then
                           efab(i,j,k,n) &
     &                          = 2*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/ &
     &                          (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               enddo
            enddo
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).ne.0.d0) then
                           efab(i,j,k,n) &
     &                          = 2*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/ &
     &                          (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               enddo
            enddo
         else
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).ne.0.d0) then
                           efab(i,j,k,n) &
     &                          = 2*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/ &
     &                          (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n) = 0.d0
                        endif
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
      end

end module amrex_abec4_module
