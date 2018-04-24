
module amrex_lo_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_reflect_odd

  implicit none

contains

!-----------------------------------------------------------------------
    subroutine amrex_lo_harmonic_averageec ( &
           c, c_l1,c_l2,c_h1,c_h2, &
           f, f_l1,f_l2,f_h1,f_h2, &
           lo, hi, nc, &
           cdir &
           ) bind(c,name='amrex_lo_harmonic_averageec')

      implicit none

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer f_l1,f_l2,f_h1,f_h2
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,nc)
      integer c_l1,c_l2,c_h1,c_h2
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,nc)

      real(amrex_real) factor, den
      parameter(factor=2.00D0)
      integer n
      integer i
      integer j

      if ( cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  den = f(2*i,2*j,n) + f(2*i,2*j+1,n)
                  if (den .ne. 0.0D0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i,2*j+1,n)/den
                  else
                    c(i,j,n) =  0.0D0
                  end if
               end do
            end do
         end do
      else if (cdir .eq. 1 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  den = f(2*i,2*j,n) + f(2*i+1,2*j,n)
                  if (den .ne. 0.0D0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i+1,2*j,n)/den
                  else
                    c(i,j,n) =  0.0D0
                  end if
               end do
            end do
         end do
      end if

    end subroutine amrex_lo_harmonic_averageec
!-----------------------------------------------------------------------
    subroutine amrex_lo_averageec ( &
           c, c_l1,c_l2,c_h1,c_h2, &
           f, f_l1,f_l2,f_h1,f_h2, &
           lo, hi, nc, &
           cdir &
           ) bind(c,name='amrex_lo_averageec')

      implicit none

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer f_l1,f_l2,f_h1,f_h2
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,nc)
      integer c_l1,c_l2,c_h1,c_h2
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,nc)

      integer n
      integer i
      integer j
      real(amrex_real) denom
      parameter(denom=half)

      if (cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i,2*j+1,n))*denom
               end do
            end do
         end do
      else if (cdir .eq. 1) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i+1,2*j,n))*denom
               end do
            end do
         end do
      end if

    end subroutine amrex_lo_averageec
!-----------------------------------------------------------------------
    subroutine amrex_lo_averagecc ( &
           c, c_l1,c_l2,c_h1,c_h2, &
           f, f_l1,f_l2,f_h1,f_h2, &
           lo, hi, nc &
           ) bind(c,name='amrex_lo_averagecc')

      implicit none

      integer nc
      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,nc)
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,nc)

      integer i
      integer j
      integer n
      real(amrex_real) denom
      parameter(denom=fourth)

      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) =  ( &
                    f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n) &
                    + f(2*i+1,2*j  ,n) + f(2*i  ,2*j  ,n))*denom
            end do
         end do
      end do

    end subroutine amrex_lo_averagecc
!-----------------------------------------------------------------------
    subroutine amrex_lo_applybc ( &
           flagden, flagbc, maxorder, &
           phi,   phi_l1,phi_l2,phi_h1,phi_h2, &
           cdir, bct, bcl, &
           bcval, bcval_l1,bcval_l2,bcval_h1,bcval_h2, &
           mask,  mask_l1,mask_l2,mask_h1,mask_h2, &
           den,   den_l1,den_l2,den_h1,den_h2, &
           lo, hi, nc, &
           h &
           ) bind(c,name='amrex_lo_applybc')

      use amrex_lo_util_module, only : polyInterpCoeff
      implicit none

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
!     The bc type = amrex_lo_reflect_odd is a special type of boundary condition.
 
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
      real(amrex_real) bcl
      real(amrex_real) h(BL_SPACEDIM)

      integer i
      integer j
      integer n
      logical is_dirichlet
      logical is_neumann
!      real(amrex_real) xb

      integer lenx
      integer leny
      integer m

      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      real(amrex_real) x(-1:maxmaxorder-2)
      real(amrex_real) coef(-1:maxmaxorder-2)
      real(amrex_real) xInt

      is_dirichlet(i) = ( i .eq. amrex_lo_dirichlet )
      is_neumann(i)   = ( i .eq. amrex_lo_neumann )

      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)

!     TODO:
!     In order for this to work with growing multigrid, must
!     sort xa[] because it is possible for the xb value to lay
!     within this range.
!     
!     The Left face of the grid

      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(lo(1)-1,j,n) = merge( &
                       phi(lo(1),j,n), &
                       phi(lo(1)-1,j,n), &
                       mask(lo(1)-1,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = 1.0D0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     phi(lo(1)-1, j, n) = merge( &
                          bcval(lo(1)-1,j,n)*coef(-1), &
                          phi(lo(1)-1, j, n), &
                          mask(lo(1)-1,j) .gt. 0)
                  end do
               else
                  do j = lo(2), hi(2)
                     phi(lo(1)-1, j, n) = merge( &
                          0.0D0, &
                          phi(lo(1)-1, j, n), &
                          mask(lo(1)-1,j) .gt. 0)
                  end do
               end if
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     phi(lo(1)-1,j,n) = merge( &
                          phi(lo(1)-1,j,n) &
                          + phi(lo(1)+m, j, n)*coef(m), &
                          phi(lo(1)-1,j,n), &
                          mask(lo(1)-1,j) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = merge(coef(0), 0.0D0, &
                       mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. amrex_lo_reflect_odd ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(lo(1)-1, j, n) = merge( &
                      -phi(lo(1),j,n), &
                       phi(lo(1)-1, j, n), &
                       mask(lo(1)-1,j) .gt. 0)
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

!     The Right face of the grid

      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(hi(1)+1,j,n) = merge( &
                       phi(hi(1), j, n), &
                       phi(hi(1)+1, j, n), &
                       mask(hi(1)+1,j) .gt. 0)
               end do
            end do
	    if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j) = 1.0D0
               end do
	    end if
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,n) = merge( &
                          bcval(hi(1)+1,j,n)*coef(-1), &
                          phi(hi(1)+1,j,n), &
                          mask(hi(1)+1,j) .gt. 0)
                  end do
               else
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,n) = merge( &
                          0.0D0, &
                          phi(hi(1)+1,j,n), &
                          mask(hi(1)+1,j) .gt. 0)
                  end do
               end if
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,n) = merge( &
                          phi(hi(1)+1,j,n) &
                          + phi(hi(1)-m,j,n)*coef(m), &
                          phi(hi(1)+1,j,n), &
                          mask(hi(1)+1,j) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j)   = merge(coef(0), 0.0D0, &
                       mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. amrex_lo_reflect_odd ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(hi(1)+1, j, n) = merge( &
                      -phi(hi(1),j,n), &
                       phi(hi(1)+1, j, n), &
                       mask(hi(1)+1,j) .gt. 0)
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

!     The Bottom of the Grid

      if(cdir .eq. 1) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1),hi(1)
                  phi(i,lo(2)-1,n) = merge( &
                       phi(i,lo(2),n), &
                       phi(i,lo(2)-1,n), &
                       mask(i,lo(2)-1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1),hi(1)
                  den(i,lo(2))   = 1.0D0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do m=0,leny
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(2)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, leny+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     phi(i,lo(2)-1,n) = merge( &
                          bcval(i,lo(2)-1,n)*coef(-1), &
                          phi(i,lo(2)-1,n), &
                          mask(i,lo(2)-1) .gt. 0)
                  end do
               else
                  do i = lo(1), hi(1)
                     phi(i,lo(2)-1,n) = merge( &
                          0.0D0, &
                          phi(i,lo(2)-1,n), &
                          mask(i,lo(2)-1) .gt. 0)
                  end do
               end if
               do m = 0, leny
                  do i = lo(1), hi(1)
                     phi(i, lo(2)-1, n) = merge( &
                          phi(i, lo(2)-1,n) &
                          + phi(i, lo(2)+m,n)*coef(m), &
                          phi(i, lo(2)-1, n), &
                          mask(i, lo(2)-1) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i, lo(2))   = merge(coef(0), 0.0D0, &
                       mask(i, lo(2)-1) .gt. 0)
               end do
            end if

         else if ( bct .eq. amrex_lo_reflect_odd ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,lo(2)-1,n) = merge( &
                      -phi(i,lo(2),n), &
                       phi(i,lo(2)-1,n), &
                       mask(i,lo(2)-1) .gt. 0)
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

!     The top of the grid

      if (cdir .eq. 3) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,hi(2)+1,n) = merge( &
                       phi(i,hi(2),n), &
                       phi(i,hi(2)+1,n), &
                       mask(i,hi(2)+1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2))   = 1.0D0
               end do
            end if
         else if (is_dirichlet(bct)) then
            if ( bct .eq. amrex_lo_reflect_odd ) leny = 0
            do m=0,leny
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(2)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, leny+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     phi(i,hi(2)+1,n) = merge( &
                          bcval(i,hi(2)+1,n)*coef(-1), &
                          phi(i,hi(2)+1,n), &
                          mask(i,hi(2)+1) .gt. 0)
                  end do
               else
                  do i = lo(1), hi(1)
                     phi(i,hi(2)+1,n) = merge( &
                          0.0D0, &
                          phi(i,hi(2)+1,n), &
                          mask(i,hi(2)+1) .gt. 0)
                  end do
               end if
               do m = 0, leny
                  do i = lo(1), hi(1)
                     phi(i, hi(2)+1,n) = merge( &
                          phi(i,hi(2)+1,n) &
                          + phi(i, hi(2)-m,n)*coef(m), &
                          phi(i,hi(2)+1,n), &
                          mask(i,hi(2)+1) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2))   = merge(coef(0), 0.0D0, &
                       mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else if ( bct .eq. amrex_lo_reflect_odd ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,hi(2)+1,n) = merge( &
                      -phi(i,hi(2),n), &
                       phi(i,hi(2)+1,n), &
                       mask(i,hi(2)+1) .gt. 0)
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

    end subroutine amrex_lo_applybc

end module amrex_lo_module
