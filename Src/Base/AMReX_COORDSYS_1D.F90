
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_COORDSYS_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 1

! :: ----------------------------------------------------------
! :: SETVOL
! ::             Compute the volume of each cell
! ::
! :: INPUTS / OUTPUTS:
! ::  vol         <=  volume array
! ::  vlo,vhi      => index limits of vol array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ, 2 = spherical)
! :: ----------------------------------------------------------

  subroutine FORT_SETVOL(reg_l1,reg_h1,vol,vol_l1,vol_h1,offset,dx,coord)

    implicit none

    integer    reg_l1,reg_h1
    integer    vol_l1,vol_h1
    integer    coord
    REAL_T     dx(SDIM), offset(SDIM)
    REAL_T     vol(vol_l1:vol_h1)

    integer    i
    REAL_T     ri, ro, v
    REAL_T     RZFACTOR
    parameter (RZFACTOR = two*Pi)

    if (coord .eq. 0) then

       v = dx(1)
       do i = ARG_L1(reg), ARG_H1(reg)
          vol(i) = v
       enddo

    else if (coord .eq. 1) then

       do i = ARG_L1(reg), ARG_H1(reg)
          ri = offset(1) + dx(1)*i
          ro = ri + dx(1)
          v = half*(RZFACTOR)*(ro - ri)*(ro + ri)
          vol(i) = abs(v)
       enddo

    else if (coord .eq. 2) then

       do i = ARG_L1(reg), ARG_H1(reg)
          ri = offset(1) + dx(1)*i
          ro = ri + dx(1)
          v = (two3rd*RZFACTOR)*(ro - ri)*(ro**2 + ro*ri + ri**2)
          vol(i) = abs(v)
       enddo

    else

       call bl_abort('bogus value of coord... bndrylib::SETVOL')

    endif

  end subroutine FORT_SETVOL



  subroutine FORT_SETVOLPT(vol, volloi1, volhii1, &
       ro, roloi1, rohii1, ri, riloi1, rihii1, dx, coord)

    integer volloi1, volhii1
    integer roloi1, rohii1, riloi1, rihii1
    integer coord
    REAL_T dx(SDIM)
    REAL_T vol(volloi1:volhii1)
    REAL_T ro(roloi1:rohii1)
    REAL_T ri(riloi1:rihii1)

    integer i
    REAL_T     RZFACTOR
    parameter (RZFACTOR = two*Pi)

    !  note that dx is usually unity.  dx not unity is used by the nfluid
    !  slic reconstruction

    if (coord .eq. 0) then

       do i = roloi1, rohii1
          vol(i) = (ro(i)-ri(i))
       enddo

    else if (coord .eq. 1) then

       do i = roloi1, rohii1
          vol(i) = half*RZFACTOR*(ro(i) - ri(i))*(ro(i) + ri(i))
          vol(i) = abs(vol(i))
       enddo

    else if (coord .eq. 2) then

       do i = roloi1, rohii1
          vol(i) = two3rd*RZFACTOR*(ro(i) - ri(i))*(ro(i)**2 + ro(i)*ri(i) + ri(i)**2)
          vol(i) = abs(vol(i))
       enddo

    else

       call bl_abort('bogus value of coord ... bndrylib::SETVOLPT')
    endif

  end subroutine FORT_SETVOLPT

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

  subroutine FORT_SETDLOGA(dloga,dloga_l1,dloga_h1,offset,dx,dir,coord)

    integer    dloga_l1,dloga_h1
    integer    coord
    REAL_T     dx(SDIM), offset(SDIM)
    REAL_T     dloga(dloga_l1:dloga_h1)
    integer dir

    integer    i
    REAL_T     rc

    if (coord .eq. 0) then

       do i = ARG_L1(dloga), ARG_H1(dloga)
          dloga(i) = zero
       enddo

    else if (coord .eq. 1) then

       do i = ARG_L1(dloga), ARG_H1(dloga)
          rc = offset(1) + dx(1)*(dble(i) + half)
          dloga(i) = one / rc
       enddo

    else if (coord .eq. 2) then

       do i = ARG_L1(dloga), ARG_H1(dloga)
          rc = offset(1) + dx(1)*(dfloat(i) + half)
          dloga(i) = two/rc
       enddo

    else

       call bl_abort('SETDLOGA: illegal coordinate system')

    endif

  end subroutine FORT_SETDLOGA

! :: ----------------------------------------------------------
! :: SETAREA
! ::             Compute the area of given cell face
! ::
! :: INPUTS / OUTPUTS:
! ::  area        <=  area array
! ::  alo,ahi      => index limits of area array
! ::  offset       => shift to origin of computational domain
! ::  dx           => cell size
! ::  coord        => coordinate flag (0 =cartesian, 1 = RZ, 2 = spherical)
! :: ----------------------------------------------------------

  subroutine FORT_SETAREA(reg_l1,reg_h1,area,area_l1,area_h1,offset,dx,dir,coord)

    integer    reg_l1,reg_h1
    integer    area_l1,area_h1
    integer    coord, dir
    REAL_T     dx(SDIM), offset(SDIM)
    REAL_T     area(area_l1:area_h1)

    integer    i
    REAL_T     ri, a
    REAL_T     RZFACTOR
    parameter (RZFACTOR = two*Pi)

    if (coord .eq. 0) then

       do i = ARG_L1(reg), ARG_H1(reg)
          area(i) = one
       enddo

    else if (coord .eq. 1) then

       do i = ARG_L1(reg), ARG_H1(reg)
          ri = offset(1) + dx(1)*dble(i)
          a = RZFACTOR*ri
          area(i) = abs(a)
       enddo

    else if( coord .eq. 2) then

       do i = ARG_L1(reg), ARG_H1(reg)
          ri = offset(1) + dx(1)*dble(i)
          a = two*RZFACTOR*ri*ri
          area(i) = abs(a)
       enddo

    else

       call bl_abort('bogus value for coord... SETAREA')

    endif

  end subroutine FORT_SETAREA
