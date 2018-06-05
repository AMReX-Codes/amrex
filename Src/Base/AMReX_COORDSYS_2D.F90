
module amrex_coordsys_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

! :: ----------------------------------------------------------
! :: SETVOL
! ::             Compute the volume of each cell
! ::
! :: INPUTS / OUTPUTS:
! ::  vol         <=  volume array
! ::  vlo,vhi      => index limits of vol array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ, 2 = RTHETA)
! :: ----------------------------------------------------------

  subroutine AMREX_SETVOL(reg_l1,reg_l2,reg_h1,reg_h2,vol,vol_l1,vol_l2,vol_h1,vol_h2,offset,dx,coord) &
       bind(c,name='amrex_setvol')

    implicit none

    integer    reg_l1,reg_l2,reg_h1,reg_h2
    integer    vol_l1,vol_l2,vol_h1,vol_h2
    integer    coord
    real(amrex_real)     dx(2), offset(2)
    real(amrex_real)     vol(vol_l1:vol_h1,vol_l2:vol_h2)

    integer    i, j
    real(amrex_real)     ri, ro, pi, po, v
    real(amrex_real)     RZFACTOR
    parameter (RZFACTOR = two*M_PI)
       
    if (coord .eq. 0) then

       ! cartesian

       v = dx(1)*dx(2)
       do j = reg_l2, reg_h2
          do i = reg_l1, reg_h1
             vol(i,j) = v
          end do
       end do

    elseif(coord .eq. 1) then

       ! R-Z

       do i = reg_l1, reg_h1
          ri = offset(1) + dx(1)*i
          ro = ri + dx(1)
          v = (half*RZFACTOR)*dx(2)*dx(1)*(ro + ri)
          do j = reg_l2, reg_h2
             vol(i,j) = abs(v)
          end do
       end do

    elseif(coord .eq. 2) then

       ! R-THETA

       do i = reg_l1, reg_h1
          ri = offset(1) + dx(1)*i
          ro = ri + dx(1)
          do j = reg_l2, reg_h2
             pi = offset(2) + dx(2)*j
             po = pi + dx(2)
             v = RZFACTOR*(ro - ri)*(ro**2 + ro*ri + ri**2)*(cos(pi)-cos(po))/three
             vol(i,j) = abs(v)
          enddo
       enddo

    end if
       
  end subroutine AMREX_SETVOL

!========================================================

! :: ----------------------------------------------------------
! :: SETDLOGA
! ::             Compute  d(log(A))/dr in each cell
! ::
! :: INPUTS / OUTPUTS:
! ::  dloga        <=  dloga array
! ::  dlo,dhi      => index limits of dloga array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
! :: ----------------------------------------------------------

  subroutine AMREX_SETDLOGA(dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2,offset,dx,dir,coord) &
       bind(c,name='amrex_setdloga')

    implicit none

    integer    dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer    coord
    real(amrex_real)     dx(2), offset(2)
    real(amrex_real)     dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    integer dir
       
    integer    i, j
    real(amrex_real)     rc, dlga, po, pi
       
    if (coord .eq. 0) then

       do j = dloga_l2, dloga_h2
          do i = dloga_l1, dloga_h1
             dloga(i,j) = zero
          end do
       end do

    else if( coord .eq. 1 ) then

       if (dir .eq. 0) then

          do i = dloga_l1, dloga_h1
             rc = offset(1) + dx(1)*(dble(i)+half)
             dlga = 1.d0/rc
             do j = dloga_l2, dloga_h2
                dloga(i,j) = dlga
             end do
          end do

       else if (dir .eq. 1) then

          do i = dloga_l1, dloga_h1
             do j = dloga_l2, dloga_h2
                dloga(i,j) = zero
             end do
          end do

       else

          call bl_abort('setdloga: illegal direction')

       end if

    else if( coord .eq. 2) then

       if (dir .eq. 0) then
          do i = dloga_l1, dloga_h1
             rc = offset(1) + dx(1)*(dble(i)+half)
             dlga = 2.d0/rc
             do j = dloga_l2, dloga_h2
                dloga(i,j) = dlga
             enddo
          enddo

       else if (dir .eq. 1) then

          do i = dloga_l1, dloga_h1
             rc = offset(1) + dx(1)*(dble(i)+half)
             dlga = 1.d0/rc
             do j = dloga_l2, dloga_h2
                pi = offset(2) + dx(2)*j
                po = pi + dx(2)
                dloga(i,j) = dlga/tan(half*(pi+po))
             enddo
          enddo

       else

          call bl_abort('setdloga: illegal coordinate system')

       endif

    end if
   
  end subroutine AMREX_SETDLOGA

! :: ----------------------------------------------------------
! :: SETAREA
! ::             Compute the area of given cell face
! ::
! :: INPUTS / OUTPUTS:
! ::  area        <=  area array
! ::  alo,ahi      => index limits of area array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 =cartesian, 1 = RZ)
! :: ----------------------------------------------------------

  subroutine AMREX_SETAREA(reg_l1,reg_l2,reg_h1,reg_h2,area,area_l1,area_l2,area_h1,area_h2,offset,dx,dir,coord) &
       bind(c,name='amrex_setarea')

    implicit none
    integer    reg_l1,reg_l2,reg_h1,reg_h2
    integer    area_l1,area_l2,area_h1,area_h2
    integer    coord, dir
    real(amrex_real)     dx(2), offset(2)
    real(amrex_real)     area(area_l1:area_h1,area_l2:area_h2)

    integer    i, j
    real(amrex_real)     rc, ri, ro, a, pi, po
    real(amrex_real)     RZFACTOR
    parameter (RZFACTOR = two*M_PI)
       
    if (coord .eq. 0) then

       ! cartesian

       if (dir .eq. 0) then

          do j = reg_l2, reg_h2
             do i = reg_l1, reg_h1
                area(i,j) = dx(2)
             end do
          end do

       else

          do j = reg_l2, reg_h2
             do i = reg_l1, reg_h1
                area(i,j) = dx(1)
             end do
          end do

       end if

    else if (coord .eq. 1) then

       ! R-Z

       if (dir .eq. 0) then

          do i = reg_l1, reg_h1
             ri = offset(1) + dx(1)*i
             a = abs(RZFACTOR*ri*dx(2))
             do j = reg_l2, reg_h2
                area(i,j) = a
             end do
          end do

       else

          do i = reg_l1, reg_h1
             rc = offset(1) + dx(1)*(dble(i)+half)
             a = abs(dx(1)*RZFACTOR*rc)
             do j = reg_l2, reg_h2
                area(i,j) = a
             end do
          end do

       end if

    elseif(coord .eq. 2) then

       if (dir .eq. 0) then

          do i = reg_l1, reg_h1
             ri = offset(1) + dx(1)*i
             do j = reg_l2, reg_h2
                pi = offset(2) + dx(2)*j
                po = pi + dx(2)
                a = RZFACTOR*ri*ri*(cos(pi)-cos(po))
                area(i,j) = abs(a)
             enddo
          enddo

       elseif(dir .eq. 1) then

          do i = reg_l1, reg_h1
             ri = offset(1) + dx(1)*i
             ro = ri + dx(1)
             do j = reg_l2, reg_h2
                pi = offset(2) + dx(2)*j
                a = RZFACTOR*sin(pi)*(ro - ri)*(ro + ri)/two
                area(i,j) = abs(a)
             enddo
          enddo

       else

          write(6,*)' bogus dir ', dir
          call bl_abort(" ")

       endif

    end if
       
  end subroutine AMREX_SETAREA

end module amrex_coordsys_module

