
module amrex_abec4_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  include 'AMReX_lo_bctypes.fi'

contains

!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
      subroutine amrex_ab4_adotx( &
          y,y_l1,y_l2,y_h1,y_h2, &
          x,x_l1,x_l2,x_h1,x_h2, &
          alpha, beta, &
          a,a_l1,a_l2,a_h1,a_h2, &
          b,b_l1,b_l2,b_h1,b_h2, &
          lo,hi,nc, &
          h &
          ) bind(c,name='amrex_ab4_adotx')

      implicit none

      real(amrex_real) alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer y_l1,y_l2,y_h1,y_h2
      integer x_l1,x_l2,x_h1,x_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer b_l1,b_l2,b_h1,b_h2
      real(amrex_real)  y(y_l1:y_h1,y_l2:y_h2,nc)
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real)  b(b_l1:b_h1,b_l2:b_h2)
      real(amrex_real) h(BL_SPACEDIM)

      integer i,j,n
      real(amrex_real) bh1i, bh2i
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)

      real(amrex_real), allocatable :: xflux(:,:,:), yflux(:,:,:)

      xlo(:) = lo(:)
      xhi(:) = hi(:)
      ylo(:) = lo(:)
      yhi(:) = hi(:)
      xhi(1) = hi(1) + 1
      yhi(2) = hi(2) + 1

      allocate(xflux(xlo(1):xhi(1),xlo(2):xhi(2),nc))
      allocate(yflux(ylo(1):yhi(1),ylo(2):yhi(2),nc))

      bh1i = beta / h(1)
      bh2i = beta / h(2)

      call flux_dir(x, x_l1,x_l2,x_h1,x_h2, alpha, beta, &
           a, a_l1,a_l2,a_h1,a_h2, &
           b, b_l1,b_l2,b_h1,b_h2, &
           nc, h(1), xlo, xhi, xflux, xlo(1), xlo(2), xhi(1), xhi(2), 1)
      call flux_dir(x, x_l1,x_l2,x_h1,x_h2, alpha, beta, &
           a, a_l1,a_l2,a_h1,a_h2, &
           b, b_l1,b_l2,b_h1,b_h2, &
           nc, h(2), ylo, yhi, yflux, ylo(1), ylo(2), yhi(1), yhi(2), 2)

      do n=1,nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               y(i,j,n) = alpha*a(i,j)*x(i,j,n) &
     &          - bh1i*(xflux(i+1,j,n) - xflux(i,j,n)) &
     &          - bh2i*(yflux(i,j+1,n) - yflux(i,j,n))
            enddo
         enddo
      enddo

      deallocate(xflux)
      deallocate(yflux)

      end

!-----------------------------------------------------------------------
!>
!>     Fill in fluxes
!>
      subroutine amrex_ab4_flux( &
          x,x_l1,x_l2,x_h1,x_h2, &
          alpha, beta, &
          a, a_l1,a_l2,a_h1,a_h2, &
          b, b_l1,b_l2,b_h1,b_h2, &
          nc, &
          h, &
          xlo,xhi,xflux,xflux_l1,xflux_l2,xflux_h1,xflux_h2, &
          ylo,yhi,yflux,yflux_l1,yflux_l2,yflux_h1,yflux_h2 &
          ) bind(c,name='amrex_ab4_flux')
      implicit none
      real(amrex_real) alpha, beta
      integer nc
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer x_l1,x_l2,x_h1,x_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer b_l1,b_l2,b_h1,b_h2
      integer xflux_l1,xflux_l2,xflux_h1,xflux_h2
      integer yflux_l1,yflux_l2,yflux_h1,yflux_h2
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real)  b(b_l1:b_h1,b_l2:b_h2)
      real(amrex_real) xflux(xflux_l1:xflux_h1,xflux_l2:xflux_h2,nc)
      real(amrex_real) yflux(yflux_l1:yflux_h1,yflux_l2:yflux_h2,nc)
      real(amrex_real) h(BL_SPACEDIM)

      call flux_dir(x, x_l1,x_l2,x_h1,x_h2, alpha, beta, &
           a, a_l1,a_l2,a_h1,a_h2, &
           b, b_l1,b_l2,b_h1,b_h2, &
           nc, h(1), xlo, xhi, xflux, xflux_l1,xflux_l2,xflux_h1,xflux_h2, 1)
      call flux_dir(x, x_l1,x_l2,x_h1,x_h2, alpha, beta, &
           a, a_l1,a_l2,a_h1,a_h2, &
           b, b_l1,b_l2,b_h1,b_h2, &
           nc,h(2), ylo, yhi, yflux, yflux_l1,yflux_l2,yflux_h1,yflux_h2, 2)

      end

!-----------------------------------------------------------------------
      subroutine flux_dir( &
          x,x_l1,x_l2,x_h1,x_h2, &
          alpha, beta, &
          a, a_l1,a_l2,a_h1,a_h2, &
          b, b_l1,b_l2,b_h1,b_h2, &
          nc, &
          h, &
          lo,hi,flux,flux_l1,flux_l2,flux_h1,flux_h2, &
          dir)
      implicit none
      real(amrex_real) alpha, beta
      integer nc, dir
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer x_l1,x_l2,x_h1,x_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer b_l1,b_l2,b_h1,b_h2
      integer flux_l1,flux_l2,flux_h1,flux_h2
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real)  b(b_l1:b_h1,b_l2:b_h2)
      real(amrex_real) flux(flux_l1:flux_h1,flux_l2:flux_h2,nc)
      real(amrex_real) h

      integer i,j,n
      real(amrex_real) i12, i48
      real(amrex_real) hi1_12

      real(amrex_real), allocatable :: gt(:,:)
      real(amrex_real), allocatable :: bt(:,:)
      real(amrex_real) g1, b1

      i12 = 1.d0/12.d0
      i48 = 1.d0/48.d0
      hi1_12 = 1.d0 / (12.d0*h)

      if (dir .eq. 1) then
         allocate(gt(lo(1):hi(1),lo(2)-2:hi(2)+2))
         allocate(bt(lo(1):hi(1),lo(2)-2:hi(2)+2))

         do j = lo(2)-2, hi(2)+2
            do i = lo(1), hi(1)
               bt(i,j) = (-b(i-2,j)+7*(b(i-1,j)+b(i,j))-b(i+1,j))*i12
            enddo
         enddo
         do n=1,nc
            do j = lo(2)-2, hi(2)+2
               do i = lo(1), hi(1)
                  gt(i,j) = (x(i-2,j,n)-x(i+1,j,n)+15.d0*(x(i,j,n)-x(i-1,j,n)))*hi1_12
               enddo
            enddo
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  g1 = (34.d0*(gt(i,j+1)-gt(i,j-1))+5.d0*(gt(i,j-2)-gt(i,j+2)))*i48
                  b1 = (34.d0*(bt(i,j+1)-bt(i,j-1))+5.d0*(bt(i,j-2)-bt(i,j+2)))*i48
                  flux(i,j,n) = bt(i,j)*gt(i,j)  +  i12*b1*g1
               enddo
            enddo
         enddo
      else
         allocate(gt(lo(1)-2:hi(1)+2,lo(2):hi(2)))
         allocate(bt(lo(1)-2:hi(1)+2,lo(2):hi(2)))

         do j = lo(2), hi(2)
            do i = lo(1)-2, hi(1)+2
               bt(i,j) = (-b(i,j-2)+7*(b(i,j-1)+b(i,j))-b(i,j+1))*i12
            enddo
         enddo
         do n=1,nc
            do j = lo(2), hi(2)
               do i = lo(1)-2, hi(1)+2
                  gt(i,j) = (x(i,j-2,n)-x(i,j+1,n)+15.d0*(x(i,j,n)-x(i,j-1,n)))*hi1_12
               enddo
            enddo
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  g1 = (34.d0*(gt(i+1,j)-gt(i-1,j))+5.d0*(gt(i-2,j)-gt(i+2,j)))*i48
                  b1 = (34.d0*(bt(i+1,j)-bt(i-1,j))+5.d0*(bt(i-2,j)-bt(i+2,j)))*i48
                  flux(i,j,n) = bt(i,j)*gt(i,j)  +  i12*b1*g1
               enddo
            enddo
         enddo
      endif

      deallocate(gt)
      deallocate(bt)

      end

!-----------------------------------------------------------------------
      subroutine amrex_ab4_applybc4 ( &
          flagden, flagbc, maxorder, &
          phi,   phi_l1,phi_l2,phi_h1,phi_h2, &
          cdir, bct, bcl, &
          bcval, bcval_l1,bcval_l2,bcval_h1,bcval_h2, &
          mask,  mask_l1,mask_l2,mask_h1,mask_h2, &
          den,   den_l1,den_l2,den_h1,den_h2, &
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
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer phi_l1,phi_l2,phi_h1,phi_h2
      real(amrex_real) phi(phi_l1:phi_h1,phi_l2:phi_h2,nc)
      integer den_l1,den_l2,den_h1,den_h2
      real(amrex_real) den(den_l1:den_h1,den_l2:den_h2)
      integer bcval_l1,bcval_l2,bcval_h1,bcval_h2
      real(amrex_real) bcval(bcval_l1:bcval_h1,bcval_l2:bcval_h2,nc)
      integer mask_l1,mask_l2,mask_h1,mask_h2
      integer mask(mask_l1:mask_h1,mask_l2:mask_h2)
      integer bct
      real(amrex_real) bcl, tmp, b, d
      parameter (b = 1.d0/24.d0)
      parameter (d = 11.d0/12.d0)
      real(amrex_real) h(BL_SPACEDIM)
!
      integer i, j, n
      logical is_dirichlet, is_neumann
!
!
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
!
!
!     The Left face of the grid
!
      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  if (mask(lo(1)-1,j) .gt. 0) then
!                     phi(lo(1)-1,j,n) =  (9*phi(lo(1),j,n) +  3*phi(lo(1)+1,j,n) - phi(lo(1)+2,j,n))/11.d0
                     phi(lo(1)-1,j,n) =  (5*phi(lo(1),j,n) +  9*phi(lo(1)+1,j,n) - 5*phi(lo(1)+2,j,n) + phi(lo(1)+3,j,n))/10.d0
                  endif
                  if (mask(lo(1)-2,j) .gt. 0) then
!                     phi(lo(1)-2,j,n) = (- 30*phi(lo(1),j,n) + 56*phi(lo(1)+1,j,n) - 15*phi(lo(1)+2,j,n))/11.d0
                     phi(lo(1)-2,j,n) = (- 75*phi(lo(1),j,n) + 145*phi(lo(1)+1,j,n) - 75*phi(lo(1)+2,j,n) &
                          + 15*phi(lo(1)+3,j,n))/10.d0
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
!                        phi(lo(1)-1,j,n) = (12*tmp - 13*phi(lo(1),j,n) + 5*phi(lo(1)+1,j,n) - phi(lo(1)+2,j,n))/3.d0
                        phi(lo(1)-1,j,n) = (60*tmp - 77*phi(lo(1),j,n) + 43*phi(lo(1)+1,j,n) - 17*phi(lo(1)+2,j,n) &
                             + 3*phi(lo(1)+3,j,n))/12.d0
                        if (mask(lo(1)-2,j) .gt. 0) then
!                           phi(lo(1)-2,j,n) = (48*tmp - 70*phi(lo(1),j,n) + 32*phi(lo(1)+1,j,n) - 7*phi(lo(1)+2,j,n))/3.d0
                           phi(lo(1)-2,j,n) = (300*tmp - 505*phi(lo(1),j,n) + 335*phi(lo(1)+1,j,n) - 145*phi(lo(1)+2,j,n) &
                                + 27*phi(lo(1)+3,j,n))/12.d0
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
                  den(lo(1),j) = merge(66.D0, 0.0D0, &
                      mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(lo(1)-1, j, n) = merge( &
                     -phi(lo(1),j,n), &
                      phi(lo(1)-1, j, n), &
                      mask(lo(1)-1,j) .gt. 0)
                  phi(lo(1)-2, j, n) = merge( &
                     -phi(lo(1)+1,j,n), &
                      phi(lo(1)-2, j, n), &
                      mask(lo(1)-2,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = merge(-1.0D0, 0.0D0, &
                      mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The Right face of the grid
!
      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  if (mask(hi(1)+1,j) .gt. 0) then
!                     phi(hi(1)+1,j,n) =  (9*phi(hi(1),j,n) +  3*phi(hi(1)-1,j,n) - phi(hi(1)-2,j,n))/11.d0
                     phi(hi(1)+1,j,n) =  (5*phi(hi(1)  ,j,n) +  9*phi(hi(1)-1,j,n) - 5*phi(hi(1)-2,j,n) &
                          + phi(hi(1)-3,j,n))/10.d0
                  endif
                  if (mask(hi(1)+2,j) .gt. 0) then
!                     phi(hi(1)+2,j,n) = (- 30*phi(hi(1),j,n) + 56*phi(hi(1)-1,j,n) - 15*phi(hi(1)-2,j,n))/11.d0
                     phi(hi(1)+2,j,n) = (- 75*phi(hi(1),j,n) + 145*phi(hi(1)-1,j,n) - 75*phi(hi(1)-2,j,n) &
                          + 15*phi(hi(1)-3,j,n))/10.d0
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
!                        phi(hi(1)+1,j,n) = (12*tmp - 13*phi(hi(1),j,n) + 5*phi(hi(1)-1,j,n) - phi(hi(1)-2,j,n))/3.d0
                        phi(hi(1)+1,j,n) = (60*tmp - 77*phi(hi(1),j,n) + 43*phi(hi(1)-1,j,n) - 17*phi(hi(1)-2,j,n) &
                             + 3*phi(hi(1)-3,j,n))/12.d0
                        if (mask(hi(1)+2,j) .gt. 0) then
!                            phi(hi(1)+2,j,n) = (48*tmp - 70*phi(hi(1),j,n) + 32*phi(hi(1)-1,j,n) - 7*phi(hi(1)-2,j,n))/3.d0
                            phi(hi(1)+2,j,n) = (300*tmp - 505*phi(hi(1),j,n) + 335*phi(hi(1)-1,j,n) - 145*phi(hi(1)-2,j,n) &
                                 + 27*phi(hi(1)-3,j,n))/12.d0
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
                  den(hi(1),j)   = merge(66.D0, 0.0D0, &
                      mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(hi(1)+1, j, n) = merge( &
                     -phi(hi(1),j,n), &
                      phi(hi(1)+1, j, n), &
                      mask(hi(1)+1,j) .gt. 0)
                  phi(hi(1)+2, j, n) = merge( &
                     -phi(hi(1)-1,j,n), &
                      phi(hi(1)+2, j, n), &
                      mask(hi(1)+2,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j) = merge(-1.0D0, 0.0D0, &
                      mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The Bottom of the Grid
!
      if(cdir .eq. 1) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1),hi(1)
                  if (mask(i,lo(2)-1) .gt. 0) then
!                     phi(i,lo(2)-1,n) =  (9*phi(i,lo(2),n) +  3*phi(i,lo(2)+1,n) - phi(i,lo(2)+2,n))/11.d0
                     phi(i,lo(2)-1,n) =  (5*phi(i,lo(2)  ,n) +  9*phi(i,lo(2)+1,n) - 5*phi(i,lo(2)+2,n) &
                          + phi(i,lo(2)+3,n))/10.d0
                  endif
                  if (mask(i,lo(2)-2) .gt. 0) then
!                     phi(i,lo(2)-2,n) = (- 30*phi(i,lo(2),n) + 56*phi(i,lo(2)+1,n) - 15*phi(i,lo(2)+2,n))/11.d0
                     phi(i,lo(2)-2,n) = (- 75*phi(i,lo(2),n) + 145*phi(i,lo(2)+1,n) - 75*phi(i,lo(2)+2,n) &
                          + 15*phi(i,lo(2)+3,n))/10.d0
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
!                        phi(i,lo(2)-1,n) = (12*tmp - 13*phi(i,lo(2),n) + 5*phi(i,lo(2)+1,n) - phi(i,lo(2)+2,n))/3.d0
                        phi(i,lo(2)-1,n) = (60*tmp - 77*phi(i,lo(2),n) + 43*phi(i,lo(2)+1,n) - 17*phi(i,lo(2)+2,n) &
                             + 3*phi(i,lo(2)+3,n))/12.d0
                        if (mask(i,lo(2)-2) .gt. 0) then
!                           phi(i,lo(2)-2,n) =   (48*tmp - 70*phi(i,lo(2),n) + 32*phi(i,lo(2)+1,n) - 7*phi(i,lo(2)+2,n))/3.d0
                           phi(i,lo(2)-2,n) = (300*tmp - 505*phi(i,lo(2),n) + 335*phi(i,lo(2)+1,n) - 145*phi(i,lo(2)+2,n) &
                                + 27*phi(i,lo(2)+3,n))/12.d0
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
                  den(i, lo(2))   = merge(66.D0, 0.0D0, &
                      mask(i, lo(2)-1) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,lo(2)-1,n) = merge( &
                     -phi(i,lo(2),n), &
                      phi(i,lo(2)-1,n), &
                      mask(i,lo(2)-1) .gt. 0)
                  phi(i,lo(2)-2,n) = merge( &
                     -phi(i,lo(2)+1,n), &
                      phi(i,lo(2)-2,n), &
                      mask(i,lo(2)-2) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,lo(2)) = merge(-1.0D0, 0.0D0, &
                      mask(i,lo(2)-1) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The top of the grid
!
      if (cdir .eq. 3) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1), hi(1)
                  if (mask(i,hi(2)+1) .gt. 0) then
!                     phi(i,hi(2)+1,n) =  (9*phi(i,hi(2),n) +  3*phi(i,hi(2)-1,n) - phi(i,hi(2)-2,n))/11.d0
                     phi(i,hi(2)+1,n) =  (5*phi(i,hi(2),n) +  9*phi(i,hi(2)-1,n) - 5*phi(i,hi(2)-2,n) &
                          + phi(i,hi(2)-3,n))/10.d0
                  endif
                  if (mask(i,hi(2)+2) .gt. 0) then
!                     phi(i,hi(2)+2,n) = (- 30*phi(i,hi(2),n) + 56*phi(i,hi(2)-1,n) - 15*phi(i,hi(2)-2,n))/11.d0
                     phi(i,hi(2)+2,n) = (- 75*phi(i,hi(2),n) + 145*phi(i,hi(2)-1,n) - 75*phi(i,hi(2)-2,n) &
                          + 15*phi(i,hi(2)-3,n))/10.d0
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
!                        phi(i,hi(2)+1,n) = (12*tmp - 13*phi(i,hi(2),n) + 5*phi(i,hi(2)-1,n) - phi(i,hi(2)-2,n))/3.d0
                        phi(i,hi(2)+1,n) = (60*tmp - 77*phi(i,hi(2),n) + 43*phi(i,hi(2)-1,n) - 17*phi(i,hi(2)-2,n) &
                             + 3*phi(i,hi(2)-3,n))/12.d0
                        if (mask(i,hi(2)+2) .gt. 0) then
!                           phi(i,hi(2)+2,n) =   (48*tmp -  70*phi(i,hi(2),n) + 32*phi(i,hi(2)-1,n) - 7*phi(i,hi(2)-2,n))/3.d0
                           phi(i,hi(2)+2,n) = (300*tmp - 505*phi(i,hi(2),n) + 335*phi(i,hi(2)-1,n) - 145*phi(i,hi(2)-2,n) &
                                + 27*phi(i,hi(2)-3,n))/12.d0
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
                  den(i,hi(2))   = merge(66.D0, 0.0D0, &
                      mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,hi(2)+1,n) = merge( &
                     -phi(i,hi(2),n), &
                      phi(i,hi(2)+1,n), &
                      mask(i,hi(2)+1) .gt. 0)
                  phi(i,hi(2)+2,n) = merge( &
                     -phi(i,hi(2)-2,n), &
                      phi(i,hi(2)+2,n), &
                      mask(i,hi(2)+2) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2)) = merge(-1.0D0, 0.0D0, &
                      mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_applybc4_touchup ( &
          phi,   phi_l1,phi_l2,phi_h1,phi_h2, &
          lo, hi, nc) bind(c,name='amrex_ab4_applybc4_touchup')

      implicit none
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer phi_l1,phi_l2,phi_h1,phi_h2
      real(amrex_real) phi(phi_l1:phi_h1,phi_l2:phi_h2,nc)
!
      integer i, j, n

      do n = 1, nc
         do i = lo(1)-1, phi_l1, -1
            do j = lo(2)-1, phi_l2, -1
               phi(i,j,n) = 0.5d0 * &
     &              ( (4*phi(i,j+1,n)-6*phi(i,j+2,n)+4*phi(i,j+3,n)-phi(i,j+4,n)) &
     &              + (4*phi(i+1,j,n)-6*phi(i+2,j,n)+4*phi(i+3,j,n)-phi(i+4,j,n)) )
            enddo
            do j = hi(2)+1, phi_h2
               phi(i,j,n) = 0.5d0 * &
     &              ( (4*phi(i,j-1,n)-6*phi(i,j-2,n)+4*phi(i,j-3,n)-phi(i,j-4,n)) &
     &              + (4*phi(i+1,j,n)-6*phi(i+2,j,n)+4*phi(i+3,j,n)-phi(i+4,j,n)) )
            enddo
         enddo

         do i = hi(1)+1, phi_h1
            do j = lo(2)-1, phi_l2, -1
               phi(i,j,n) = 0.5d0 * &
     &              ( (4*phi(i,j+1,n)-6*phi(i,j+2,n)+4*phi(i,j+3,n)-phi(i,j+4,n)) &
     &              + (4*phi(i-1,j,n)-6*phi(i-2,j,n)+4*phi(i-3,j,n)-phi(i-4,j,n)) )
            enddo
            do j = hi(2)+1, phi_h2
               phi(i,j,n) = 0.5d0 * &
     &              ( (4*phi(i,j-1,n)-6*phi(i,j-2,n)+4*phi(i,j-3,n)-phi(i,j-4,n)) &
     &              + (4*phi(i-1,j,n)-6*phi(i-2,j,n)+4*phi(i-3,j,n)-phi(i-4,j,n)) )
            enddo
         enddo
      enddo
      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_ca2cc(lo, hi, ca, ca_l1,ca_l2,ca_h1,ca_h2, cc, cc_l1,cc_l2,cc_h1,cc_h2, nc) &
           bind(c,name='amrex_ab4_ca2cc')

      implicit none
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer ca_l1,ca_l2,ca_h1,ca_h2
      integer cc_l1,cc_l2,cc_h1,cc_h2
      integer nc
      real(amrex_real) ca(ca_l1:ca_h1,ca_l2:ca_h2,nc)
      real(amrex_real) cc(cc_l1:cc_h1,cc_l2:cc_h2,nc)

      integer i,j,n
      real(amrex_real) one24th, seven6th
      parameter (one24th = 1.d0 / 24.d0)
      parameter (seven6th = 7.d0 / 6.d0)

      do n=1,nc
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,n) = -one24th*( &
     &              ca(i,j-1,n)+ca(i-1,j,n)+ca(i+1,j,n)+ca(i,j+1,n)) &
     &              + seven6th*ca(i,j,n)
            enddo
         enddo
      enddo

      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_cc2ca(lo, hi, cc, cc_l1,cc_l2,cc_h1,cc_h2, &
           ca, ca_l1,ca_l2,ca_h1,ca_h2, nc) bind(c,name='amrex_ab4_cc2ca')

      implicit none
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cc_l1,cc_l2,cc_h1,cc_h2
      integer ca_l1,ca_l2,ca_h1,ca_h2
      integer nc
      real(amrex_real) cc(cc_l1:cc_h1,cc_l2:cc_h2,nc)
      real(amrex_real) ca(ca_l1:ca_h1,ca_l2:ca_h2,nc)

      integer i,j,n
      real(amrex_real) one24th, five6th
      parameter (one24th = 1.d0 / 24.d0)
      parameter (five6th = 5.d0 / 6.d0)

      if (cc_h1.eq.31) then
         print *,'hello'
      endif

      do n=1,nc
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ca(i,j,n) = one24th*( &
     &              cc(i,j-1,n)+cc(i-1,j,n)+cc(i+1,j,n)+cc(i,j+1,n)) &
     &              + five6th*cc(i,j,n)
            enddo
         enddo
      enddo
      if (cc_h1.eq.31) then
         print *,'goodbye'
      endif
      end

!-----------------------------------------------------------------------

      subroutine amrex_ab4_lo_cc2ec(lo, hi, &
     &     cfab, cfab_l1,cfab_l2,cfab_h1,cfab_h2, &
     &     efab, efab_l1,efab_l2,efab_h1,efab_h2, nc, dir, &
     &     isharm) bind(c,name='amrex_ab4_lo_cc2ec')
      implicit none
      integer lo(2), hi(2), nc, dir, isharm
      integer cfab_l1,cfab_l2,cfab_h1,cfab_h2
      integer efab_l1,efab_l2,efab_h1,efab_h2
      real(amrex_real)  cfab(cfab_l1:cfab_h1,cfab_l2:cfab_h2, nc)
      real(amrex_real)  efab(efab_l1:efab_h1,efab_l2:efab_h2, nc)

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
                     if((cfab(i,j,n) * cfab(i-1,j,n)).ne.0.d0)then
                        efab(i,j,n) &
     &                       = 2*(cfab(i,j,n) * cfab(i-1,j,n))/ &
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
                     if((cfab(i,j,n) * cfab(i,j-1,n)).ne.0.d0)then
                        efab(i,j,n) &
     &                       = 2*(cfab(i,j,n) * cfab(i,j-1,n))/ &
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

end module amrex_abec4_module
