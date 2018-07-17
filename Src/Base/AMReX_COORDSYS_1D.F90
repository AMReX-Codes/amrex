
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
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ, 2 = spherical)
! :: ----------------------------------------------------------

  subroutine AMREX_SETVOL(reg_l1,reg_h1,vol,vol_l1,vol_h1,offset,dx,coord) &
       bind(c,name='amrex_setvol')

    implicit none

    integer    reg_l1,reg_h1
    integer    vol_l1,vol_h1
    integer    coord
    real(amrex_real)     dx(1), offset(1)
    real(amrex_real)     vol(vol_l1:vol_h1)

    integer    i
    real(amrex_real)     ri, ro, v
    real(amrex_real)     RZFACTOR
    parameter (RZFACTOR = two*M_PI)

    if (coord .eq. 0) then

       v = dx(1)
       do i = reg_l1, reg_h1
          vol(i) = v
       enddo

    else if (coord .eq. 1) then

       do i = reg_l1, reg_h1
          ri = offset(1) + dx(1)*i
          ro = ri + dx(1)
          v = half*(RZFACTOR)*(ro - ri)*(ro + ri)
          vol(i) = abs(v)
       enddo

    else if (coord .eq. 2) then

       do i = reg_l1, reg_h1
          ri = offset(1) + dx(1)*i
          ro = ri + dx(1)
          v = (two3rd*RZFACTOR)*(ro - ri)*(ro**2 + ro*ri + ri**2)
          vol(i) = abs(v)
       enddo

    else

       call bl_abort('bogus value of coord... bndrylib::SETVOL')

    endif

  end subroutine AMREX_SETVOL



  subroutine AMREX_SETVOLPT(vol, volloi1, volhii1, &
       ro, roloi1, rohii1, ri, riloi1, rihii1, dx, coord) bind(c,name='amrex_setvolpt')

    integer volloi1, volhii1
    integer roloi1, rohii1, riloi1, rihii1
    integer coord
    real(amrex_real) dx(1)
    real(amrex_real) vol(volloi1:volhii1)
    real(amrex_real) ro(roloi1:rohii1)
    real(amrex_real) ri(riloi1:rihii1)

    integer i
    real(amrex_real)     RZFACTOR
    parameter (RZFACTOR = two*M_PI)

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

  end subroutine AMREX_SETVOLPT

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

  subroutine AMREX_SETDLOGA(dloga,dloga_l1,dloga_h1,offset,dx,dir,coord) bind(c,name='amrex_setdloga')

    integer    dloga_l1,dloga_h1
    integer    coord
    real(amrex_real)     dx(1), offset(1)
    real(amrex_real)     dloga(dloga_l1:dloga_h1)
    integer dir

    integer    i
    real(amrex_real)     rc

    if (coord .eq. 0) then

       do i = dloga_l1, dloga_h1
          dloga(i) = zero
       enddo

    else if (coord .eq. 1) then

       do i = dloga_l1, dloga_h1
          rc = offset(1) + dx(1)*(dble(i) + half)
          dloga(i) = one / rc
       enddo

    else if (coord .eq. 2) then

       do i = dloga_l1, dloga_h1
          rc = offset(1) + dx(1)*(dfloat(i) + half)
          dloga(i) = two/rc
       enddo

    else

       call bl_abort('SETDLOGA: illegal coordinate system')

    endif

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
! ::  coord        => coordinate flag (0 =cartesian, 1 = RZ, 2 = spherical)
! :: ----------------------------------------------------------

  subroutine AMREX_SETAREA(reg_l1,reg_h1,area,area_l1,area_h1,offset,dx,dir,coord) bind(c,name='amrex_setarea')

    integer    reg_l1,reg_h1
    integer    area_l1,area_h1
    integer    coord, dir
    real(amrex_real)     dx(1), offset(1)
    real(amrex_real)     area(area_l1:area_h1)

    integer    i
    real(amrex_real)     ri, a
    real(amrex_real)     RZFACTOR
    parameter (RZFACTOR = two*M_PI)

    if (coord .eq. 0) then

       do i = reg_l1, reg_h1
          area(i) = one
       enddo

    else if (coord .eq. 1) then

       do i = reg_l1, reg_h1
          ri = offset(1) + dx(1)*dble(i)
          a = RZFACTOR*ri
          area(i) = abs(a)
       enddo

    else if( coord .eq. 2) then

       do i = reg_l1, reg_h1
          ri = offset(1) + dx(1)*dble(i)
          a = two*RZFACTOR*ri*ri
          area(i) = abs(a)
       enddo

    else

       call bl_abort('bogus value for coord... SETAREA')

    endif

  end subroutine AMREX_SETAREA

end module amrex_coordsys_module
