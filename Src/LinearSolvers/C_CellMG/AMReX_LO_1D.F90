
module amrex_lo_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  include 'AMReX_lo_bctypes.fi'

contains

!-----------------------------------------------------------------------
    subroutine amrex_lo_harmonic_averageec ( &
           c, c_l1,c_h1, &
           f, f_l1,f_h1, &
           lo, hi, nc, &
           cdir &
           ) bind(c,name='amrex_lo_harmonic_averageec')

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer f_l1,f_h1
      real(amrex_real) f(f_l1:f_h1,nc)
      integer c_l1,c_h1
      real(amrex_real) c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) =  f(2*i,n)
         end do
      end do

    end subroutine amrex_lo_harmonic_averageec
!-----------------------------------------------------------------------
    subroutine amrex_lo_averageec ( &
           c, c_l1,c_h1, &
           f, f_l1,f_h1, &
           lo, hi, nc, &
           cdir &
           ) bind(c,name='amrex_lo_averageec')

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer f_l1,f_h1
      real(amrex_real) f(f_l1:f_h1,nc)
      integer c_l1,c_h1
      real(amrex_real) c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) = f(2*i,n)
         end do
      end do

    end subroutine amrex_lo_averageec
!-----------------------------------------------------------------------
    subroutine amrex_lo_averagecc ( &
           c, c_l1,c_h1, &
           f, f_l1,f_h1, &
           lo, hi, nc &
           ) bind(c,name='amrex_lo_averagecc')

      integer nc
      integer f_l1,f_h1
      integer c_l1,c_h1
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      real(amrex_real) f(f_l1:f_h1,nc)
      real(amrex_real) c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) =  (f(2*i+1,n) + f(2*i,n))*half
         end do
      end do

    end subroutine amrex_lo_averagecc
!-----------------------------------------------------------------------
    subroutine amrex_lo_applybc ( &
           flagden, flagbc, maxorder, &
           phi,   phi_l1,phi_h1, &
           cdir, bct, bcl, &
           bcval, bcval_l1,bcval_h1, &
           mask,  mask_l1,mask_h1, &
           den,   den_l1,den_h1, &
           lo, hi, nc, &
           h &
           ) bind(c,name='amrex_lo_applybc')

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

      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer phi_l1,phi_h1
      real(amrex_real) phi(phi_l1:phi_h1,nc)
      integer den_l1,den_h1
      real(amrex_real) den(den_l1:den_h1)
      integer bcval_l1,bcval_h1
      real(amrex_real) bcval(bcval_l1:bcval_h1,nc)
      integer mask_l1,mask_h1
      integer mask(mask_l1:mask_h1)
      integer bct
      real(amrex_real) bcl
      real(amrex_real) h(BL_SPACEDIM)

      integer i,n
      logical is_dirichlet
      logical is_neumann

      integer lenx
      integer m

      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      real(amrex_real) x(-1:maxmaxorder-2)
      real(amrex_real) coef(-1:maxmaxorder-2)
      real(amrex_real) xInt

      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )

      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)

!     TODO:
!     In order for this to work with growing multigrid, must
!     sort xa[] because it is possible for the xb value to lay
!     within this range.
!     
!     The Left face of the grid

      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               phi(lo(1)-1,n) = merge( &
                    phi(lo(1),n), &
                    phi(lo(1)-1,n), &
                    mask(lo(1)-1) .gt. 0)
            end do
            if ( flagden .eq. 1) then
               den(lo(1)) = 1.0D0
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
                  phi(lo(1)-1,n) = merge( &
                       bcval(lo(1)-1,n)*coef(-1), &
                       phi(lo(1)-1,n), &
                       mask(lo(1)-1) .gt. 0)
               else
                  phi(lo(1)-1,n) = merge( &
                       0.0D0, &
                       phi(lo(1)-1,n), &
                       mask(lo(1)-1) .gt. 0)
               end if
               do m = 0, lenx
                  phi(lo(1)-1,n) = merge(&
                       phi(lo(1)-1,n) &
                       + phi(lo(1)+m,n)*coef(m), &
                       phi(lo(1)-1,n), &
                       mask(lo(1)-1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               den(lo(1)) = merge(coef(0), 0.0D0, &
                    mask(lo(1)-1) .gt. 0)
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               phi(lo(1)-1,n) = merge( &
                   -phi(lo(1),n), &
                    phi(lo(1)-1,n), &
                    mask(lo(1)-1) .gt. 0)
            end do
            if ( flagden .eq. 1 ) then
               den(lo(1)) = merge(-1.0D0, 0.0D0,&
                    mask(lo(1)-1) .gt. 0)
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
               phi(hi(1)+1,n) = merge( &
                    phi(hi(1),n), &
                    phi(hi(1)+1,n), &
                    mask(hi(1)+1) .gt. 0)
            end do
	    if ( flagden .eq. 1 ) then
                den(hi(1)) = 1.0D0
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
                  phi(hi(1)+1,n) = merge( &
                       bcval(hi(1)+1,n)*coef(-1), &
                       phi(hi(1)+1,n), &
                       mask(hi(1)+1) .gt. 0)
               else
                  phi(hi(1)+1,n) = merge( &
                       0.0D0, &
                       phi(hi(1)+1,n), &
                       mask(hi(1)+1) .gt. 0)
               end if
               do m = 0, lenx
                  phi(hi(1)+1,n) = merge( &
                       phi(hi(1)+1,n) &
                       + phi(hi(1)-m,n)*coef(m), &
                       phi(hi(1)+1,n), &
                       mask(hi(1)+1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               den(hi(1))   = merge(coef(0), 0.0D0, &
                    mask(hi(1)+1) .gt. 0)
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               phi(hi(1)+1,n) = merge( &
                   -phi(hi(1),n), &
                    phi(hi(1)+1,n), &
                    mask(hi(1)+1) .gt. 0)
            end do
            if ( flagden .eq. 1 ) then
               den(hi(1)) = merge(-1.0D0, 0.0D0, &
                    mask(hi(1)+1) .gt. 0)
            end if

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if

    end subroutine amrex_lo_applybc

end module amrex_lo_module
