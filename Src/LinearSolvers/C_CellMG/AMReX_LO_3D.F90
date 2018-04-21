
module amrex_lo_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_reflect_odd

  implicit none

contains

!-----------------------------------------------------------------------
    subroutine amrex_lo_harmonic_averageec ( &
           c, c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
           f, f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
           lo, hi, nc, &
           cdir &
           ) bind(c,name='amrex_lo_harmonic_averageec')

      implicit none

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nc)
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nc)

      integer n, i, j, k

      select case(cdir)

      case (0)
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = four/( &
                          + 1.0D0/f(2*i,2*j  ,2*k  ,n) &
                          + 1.0D0/f(2*i,2*j+1,2*k  ,n) &
                          + 1.0D0/f(2*i,2*j  ,2*k+1,n) &
                          + 1.0D0/f(2*i,2*j+1,2*k+1,n) )

                  end do
               end do
            end do
         end do

      case (1)
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = four/( &
                          + 1.0D0/f(2*i  ,2*j,2*k  ,n) &
                          + 1.0D0/f(2*i+1,2*j,2*k  ,n) &
                          + 1.0D0/f(2*i  ,2*j,2*k+1,n) &
                          + 1.0D0/f(2*i+1,2*j,2*k+1,n) )

                  end do
               end do
            end do
         end do
      case (2)
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = four/( &
                           + 1.0D0/f(2*i  ,2*j  ,2*k,n) &
                          + 1.0D0/f(2*i+1,2*j  ,2*k,n) &
                          + 1.0D0/f(2*i  ,2*j+1,2*k,n) &
                          + 1.0D0/f(2*i+1,2*j+1,2*k,n) )

                  end do
               end do
            end do
         end do

      end select

    end subroutine amrex_lo_harmonic_averageec
!-----------------------------------------------------------------------
    subroutine amrex_lo_averageec ( &
           c, c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
           f, f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
           lo, hi, nc, &
           cdir &
           ) bind(c,name='amrex_lo_averageec')

      implicit none

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nc)
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nc)

      integer n, i, j, k

      select case(cdir)

      case (0)
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = fourth*( &
                          + f(2*i,2*j  ,2*k  ,n) &
                          + f(2*i,2*j+1,2*k  ,n) &
                          + f(2*i,2*j  ,2*k+1,n) &
                          + f(2*i,2*j+1,2*k+1,n) )

                  end do
               end do
            end do
         end do

      case (1)
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = fourth*( &
                          + f(2*i  ,2*j,2*k  ,n) &
                          + f(2*i+1,2*j,2*k  ,n) &
                          + f(2*i  ,2*j,2*k+1,n) &
                          + f(2*i+1,2*j,2*k+1,n) )

                  end do
               end do
            end do
         end do

      case (2)
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = fourth*( &
                          + f(2*i  ,2*j  ,2*k,n) &
                          + f(2*i+1,2*j  ,2*k,n) &
                          + f(2*i  ,2*j+1,2*k,n) &
                          + f(2*i+1,2*j+1,2*k,n) )

                  end do
               end do
            end do
         end do

      end select

    end subroutine amrex_lo_averageec
!-----------------------------------------------------------------------
    subroutine amrex_lo_averagecc ( &
           c, c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
           f, f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
           lo, hi, nc &
           ) bind(c,name='amrex_lo_averagecc')

      implicit none

      integer nc
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nc)
      real(amrex_real) c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nc)

      integer i, j, k, n

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  c(i,j,k,n) =  eighth*( &
                       + f(2*i+1,2*j+1,2*k  ,n) &
                       + f(2*i  ,2*j+1,2*k  ,n) &
                       + f(2*i+1,2*j  ,2*k  ,n) &
                       + f(2*i  ,2*j  ,2*k  ,n) &
                       + f(2*i+1,2*j+1,2*k+1,n) &
                       + f(2*i  ,2*j+1,2*k+1,n) &
                       + f(2*i+1,2*j  ,2*k+1,n) &
                       + f(2*i  ,2*j  ,2*k+1,n) )

               end do
            end do
         end do
      end do

    end subroutine amrex_lo_averagecc
!-----------------------------------------------------------------------
!
! Don't thread this.  We instead thread LinOp::applyBC() across faces.
!
    subroutine amrex_lo_applybc ( &
           flagden, flagbc, maxorder, &
           phi, phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
           cdir, bct, bcl, &
           bcval, bcval_l1,bcval_l2,bcval_l3,bcval_h1,bcval_h2,bcval_h3, &
           mask, mask_l1,mask_l2,mask_l3,mask_h1,mask_h2,mask_h3, &
           den, den_l1,den_l2,den_l3,den_h1,den_h2,den_h3, &
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
!     The bc type = amrex_lo_reflect_odd is a special type of dirichlet condition,
!     in that we want a "zeroth" order interpolant to fill the ghost cell.
!     If this were treated in the normal way, then ALL boundaries would be
!     low order.
      
      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      real(amrex_real) phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3,nc)
      integer den_l1,den_l2,den_l3,den_h1,den_h2,den_h3
      real(amrex_real) den(den_l1:den_h1,den_l2:den_h2,den_l3:den_h3)
      integer bcval_l1,bcval_l2,bcval_l3,bcval_h1,bcval_h2,bcval_h3
      real(amrex_real) bcval(bcval_l1:bcval_h1,bcval_l2:bcval_h2,bcval_l3:bcval_h3,nc)
      integer mask_l1,mask_l2,mask_l3,mask_h1,mask_h2,mask_h3
      integer mask(mask_l1:mask_h1,mask_l2:mask_h2,mask_l3:mask_h3)
      integer bct
      real(amrex_real) bcl
      real(amrex_real) h(BL_SPACEDIM)

      integer i, j, k, n
      logical is_dirichlet, is_neumann

      integer lenx, leny, lenz, m

      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      real(amrex_real) x(-1:maxmaxorder-2)
      real(amrex_real) coef(-1:maxmaxorder-2)
      real(amrex_real) xInt
      parameter(xInt = -0.5D0)

      is_dirichlet(i) = ( i .eq. amrex_lo_dirichlet )
      is_neumann(i) = (i .eq. amrex_lo_neumann)

      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
      lenz = MIN(hi(3)-lo(3), Lmaxorder-2)

      do m=0,maxmaxorder-2
         x(m) = m + 0.5D0
      end do

!     TODO:
!     In order for this to work with growing multigrid, must
!     sort xa[] because it is possible for the xb value to lay
!     within this range.

      select case (cdir)

      case (0)
         !     
         ! The Left face of the grid
         !
         if (is_neumann(bct)) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     phi(lo(1)-1,j,k,n) = merge( &
                          phi(lo(1),j,k,n), &
                          phi(lo(1)-1,j,k,n), &
                          mask(lo(1)-1,j,k) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(lo(1),j,k) = 1.0D0
                  end do
               end do
            end if
         else if (is_dirichlet(bct)) then
            x(-1) = - bcl/h(1)
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        phi(lo(1)-1, j, k, n) = merge( &
                             bcval(lo(1)-1,j,k,n)*coef(-1), &
                             phi(lo(1)-1, j,k, n), &
                             mask(lo(1)-1,j,k) .gt. 0)
                     end do
                  end do
               else
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        phi(lo(1)-1, j, k, n) = merge( &
                             0.0D0, &
                             phi(lo(1)-1, j, k, n), &
                             mask(lo(1)-1,j, k) .gt. 0)
                     end do
                  end do
               end if
               do m = 0, lenx
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        phi(lo(1)-1,j,k,n) = merge( &
                             phi(lo(1)-1,j,k,n) &
                             + phi(lo(1)+m, j, k, n)*coef(m), &
                             phi(lo(1)-1,j,k,n), &
                             mask(lo(1)-1,j,k) .gt. 0)
                     end do
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(lo(1),j,k) = merge(coef(0), 0.0D0, &
                          mask(lo(1)-1,j,k) .gt. 0)
                  end do
               end do
            end if
         else if ( bct .eq. amrex_lo_reflect_odd ) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     phi(lo(1)-1, j, k, n) = merge( &
                         -phi(lo(1),j,k,n), &
                          phi(lo(1)-1,j,k,n), &
                          mask(lo(1)-1,j,k) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(lo(1),j,k) = merge(-1.0D0, 0.0D0, &
                          mask(lo(1)-1,j,k) .gt. 0)
                  end do
               end do
            end if
         else
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if

      case (3)
         !
         ! The Right face of the grid
         ! 
         if(is_neumann(bct)) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,k,n) = merge( &
                          phi(hi(1), j, k, n), &
                          phi(hi(1)+1, j, k, n), &
                          mask(hi(1)+1,j,k) .gt. 0)
                  end do
               end do
            end do
	    if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(hi(1),j,k) = 1.0D0
                  end do
               end do
	    end if
         else if (is_dirichlet(bct)) then
            x(-1) = - bcl/h(1)
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        phi(hi(1)+1,j,k,n) = merge( &
                             bcval(hi(1)+1,j,k,n)*coef(-1), &
                             phi(hi(1)+1,j,k,n), &
                             mask(hi(1)+1,j,k) .gt. 0)
                     end do
                  end do
               else
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        phi(hi(1)+1,j,k,n) = merge( &
                             0.0D0, &
                             phi(hi(1)+1,j,k,n), &
                             mask(hi(1)+1,j,k) .gt. 0)
                     end do
                  end do
               end if
               do m = 0, lenx
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        phi(hi(1)+1,j,k,n) = merge( &
                             phi(hi(1)+1,j,k,n) &
                             + phi(hi(1)-m,j,k,n)*coef(m), &
                             phi(hi(1)+1,j,k,n), &
                             mask(hi(1)+1,j,k) .gt. 0)
                     end do
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(hi(1),j,k)   = merge(coef(0), 0.0D0, &
                          mask(hi(1)+1,j,k) .gt. 0)
                  end do
               end do
            end if
         else if ( bct .eq. amrex_lo_reflect_odd ) then
            do n = 1, nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     phi(hi(1)+1, j, k, n) = merge( &
                         -phi(hi(1),j,k,n), &
                          phi(hi(1)+1,j,k,n), &
                          mask(hi(1)+1,j,k) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     den(hi(1),j,k) = merge(-1.0D0, 0.0D0, &
                          mask(hi(1)+1,j,k) .gt. 0)
                  end do
               end do
            end if
         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if

         case (1)
            !
            ! The Bottom of the Grid
            !
            if(is_neumann(bct)) then
               do n = 1, nc
                  do k = lo(3), hi(3)
                     do i = lo(1),hi(1)
                        phi(i,lo(2)-1,k,n) = merge( &
                             phi(i,lo(2),k,n), &
                             phi(i,lo(2)-1,k,n), &
                             mask(i,lo(2)-1,k) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1),hi(1)
                        den(i,lo(2),k)   = 1.0D0
                     end do
                  end do
               end if
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(2)
               call polyInterpCoeff(xInt, x, leny+2, coef)
               do n = 1, nc
                  if ( flagbc .eq. 1 ) then
                     do k = lo(3), hi(3)
                        do i = lo(1), hi(1)
                           phi(i,lo(2)-1,k,n) = merge( &
                                bcval(i,lo(2)-1,k,n)*coef(-1), &
                                phi(i,lo(2)-1,k,n), &
                                mask(i,lo(2)-1,k) .gt. 0)
                        end do
                     end do
                  else
                     do k = lo(3), hi(3)
                        do i = lo(1), hi(1)
                           phi(i,lo(2)-1,k,n) = merge( &
                               0.0D0, &
                               phi(i,lo(2)-1,k,n), &
                               mask(i,lo(2)-1,k) .gt. 0)
                        end do
                     end do
                  end if
                  do m = 0, leny
                     do k = lo(3), hi(3)
                        do i = lo(1), hi(1)
                           phi(i, lo(2)-1, k, n) = merge( &
                               phi(i, lo(2)-1,k,n) &
                               + phi(i, lo(2)+m, k,n)*coef(m), &
                               phi(i, lo(2)-1, k, n), &
                               mask(i, lo(2)-1, k) .gt. 0)
                        end do
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i, lo(2),k)   = merge(coef(0), 0.0D0, &
                             mask(i, lo(2)-1,k) .gt. 0)
                     end do
                  end do
               end if
            else if ( bct .eq. amrex_lo_reflect_odd ) then
               do n = 1, nc
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        phi(i, lo(2)-1, k, n) = merge( &
                             -phi(i,lo(2),k,n), &
                             phi(i,lo(2)-1,k,n), &
                             mask(i,lo(2)-1,k) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,lo(2),k) = merge(-1.0D0, 0.0D0, &
                             mask(i,lo(2)-1,k) .gt. 0)
                     end do
                  end do
               end if
            else
               print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
               call bl_error("stop")
            end if

         case (4)
            !
            ! The top of the grid
            !
            if(is_neumann(bct)) then
               do n = 1, nc
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        phi(i,hi(2)+1,k,n) = merge( &
                             phi(i,hi(2),k,n), &
                             phi(i,hi(2)+1,k,n), &
                             mask(i,hi(2)+1,k) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,hi(2),k)   = 1.0D0
                     end do
                  end do
               end if
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(2)
               call polyInterpCoeff(xInt, x, leny+2, coef)
               do n = 1, nc
                  if ( flagbc .eq. 1 ) then
                     do k = lo(3), hi(3)
                        do i = lo(1), hi(1)
                           phi(i,hi(2)+1,k,n) = merge( &
                                bcval(i,hi(2)+1,k,n)*coef(-1), &
                                phi(i,hi(2)+1,k,n), &
                                mask(i,hi(2)+1,k) .gt. 0)
                        end do
                     end do
                  else
                     do k = lo(3), hi(3)
                        do i = lo(1), hi(1)
                           phi(i,hi(2)+1,k,n) = merge( &
                                0.0D0, &
                                phi(i,hi(2)+1,k,n), &
                                mask(i,hi(2)+1,k) .gt. 0)
                        end do
                     end do
                  end if
                  do m = 0, leny
                     do k = lo(3), hi(3)
                        do i = lo(1), hi(1)
                           phi(i, hi(2)+1,k,n) = merge( &
                                phi(i,hi(2)+1,k,n) &
                                + phi(i, hi(2)-m,k,n)*coef(m), &
                                phi(i,hi(2)+1,k,n), &
                                mask(i,hi(2)+1,k) .gt. 0)
                        end do
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,hi(2),k)   = merge(coef(0), 0.0D0, &
                             mask(i,hi(2)+1,k) .gt. 0)
                     end do
                  end do
               end if
            else if ( bct .eq. amrex_lo_reflect_odd ) then
               do n = 1, nc
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        phi(i, hi(2)+1, k, n) = merge( &
                            -phi(i,hi(2),k,n), &
                             phi(i,hi(2)+1,k,n), &
                             mask(i,hi(2)+1,k) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,hi(2),k) = merge(-1.0D0, 0.0D0, &
                             mask(i,hi(2)+1,k) .gt. 0)
                     end do
                  end do
               end if
            else
               print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
               call bl_error("stop")
            end if

         case (2)
            !
            ! The Front of the Grid
            !
            if(is_neumann(bct)) then
               do n = 1, nc
                  do j = lo(2), hi(2)
                     do i = lo(1),hi(1)
                        phi(i,j,lo(3)-1,n) = merge( &
                             phi(i,j,lo(3),n), &
                             phi(i,j,lo(3)-1,n), &
                             mask(i,j,lo(3)-1) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1),hi(1)
                        den(i,j,lo(3))   = 1.0D0
                     end do
                  end do
               end if
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(3)
               call polyInterpCoeff(xInt, x, lenz+2, coef)
               do n = 1, nc
                  if ( flagbc .eq. 1 ) then
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           phi(i,j,lo(3)-1,n) = merge( &
                                bcval(i,j,lo(3)-1,n)*coef(-1), &
                                phi(i,j,lo(3)-1,n), &
                                mask(i,j,lo(3)-1) .gt. 0)
                        end do
                     end do
                  else
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           phi(i,j,lo(3)-1,n) = merge( &
                                0.0D0, &
                                phi(i,j,lo(3)-1,n), &
                                mask(i,j,lo(3)-1) .gt. 0)
                        end do
                     end do
                  end if
                  do m = 0, lenz
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           phi(i, j, lo(3)-1, n) = merge( &
                                phi(i, j, lo(3)-1,n) &
                                + phi(i, j, lo(3)+m, n)*coef(m), &
                                phi(i, j, lo(3)-1,n), &
                                mask(i, j, lo(3)-1) .gt. 0)
                        end do
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i, j, lo(3))   = merge(coef(0), 0.0D0, &
                             mask(i, j, lo(3)-1) .gt. 0)
                     end do
                  end do
               end if
            else if ( bct .eq. amrex_lo_reflect_odd ) then
               do n = 1, nc
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        phi(i, j, lo(3)-1, n) = merge( &
                             -phi(i,j,lo(3),n), &
                             phi(i,j,lo(3)-1,n), &
                             mask(i,j,lo(3)-1) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j,lo(3)) = merge(-1.0D0, 0.0D0, &
                             mask(i,j,lo(3)-1) .gt. 0)
                     end do
                  end do
               end if
            else
               print *,'UNKNOWN BC ON FRONT FACE IN APPLYBC'
               call bl_error("stop")
            end if

         case (5)
            !
            ! The back of the grid
            !
            if(is_neumann(bct)) then
               do n = 1, nc
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        phi(i,j, hi(3)+1,n) = merge( &
                             phi(i,j, hi(3),n), &
                             phi(i,j, hi(3)+1,n), &
                             mask(i,j, hi(3)+1) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j, hi(3))   = 1.0D0
                     end do
                  end do
               end if
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(3)
               call polyInterpCoeff(xInt, x, lenz+2, coef)
               do n = 1, nc
                  if ( flagbc .eq. 1 ) then
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           phi(i,j, hi(3)+1,n) = merge( &
                                bcval(i,j, hi(3)+1,n)*coef(-1), &
                                phi(i,j, hi(3)+1,n), &
                                mask(i,j, hi(3)+1) .gt. 0)
                        end do
                     end do
                  else
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           phi(i,j, hi(3)+1,n) = merge( &
                                0.0D0, &
                                phi(i,j, hi(3)+1,n), &
                                mask(i,j, hi(3)+1) .gt. 0)
                        end do
                     end do
                  end if
                  do m = 0, lenz
                     do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                           phi(i, j, hi(3)+1,n) = merge( &
                                phi(i,j, hi(3)+1,n) &
                                + phi(i, j, hi(3)-m,n)*coef(m), &
                                phi(i,j, hi(3)+1,n), &
                                mask(i,j, hi(3)+1) .gt. 0)
                        end do
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j, hi(3))   = merge(coef(0), 0.0D0, &
                             mask(i,j, hi(3)+1) .gt. 0)
                     end do
                  end do
               end if
            else if ( bct .eq. amrex_lo_reflect_odd ) then
               do n = 1, nc
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        phi(i, j, hi(3)+1, n) = merge( &
                             -phi(i,j,hi(3),n), &
                             phi(i,j,hi(3)+1,n), &
                             mask(i,j,hi(3)+1) .gt. 0)
                     end do
                  end do
               end do
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j,hi(3)) = merge(-1.0D0, 0.0D0, &
                            mask(i,j,hi(3)+1) .gt. 0)
                     end do
                  end do
               end if
            else
               print *,'UNKNOWN BC ON BACK FACE IN APPLYBC'
               call bl_error("stop")
            end if
         end select

       end subroutine amrex_lo_applybc

end module amrex_lo_module
