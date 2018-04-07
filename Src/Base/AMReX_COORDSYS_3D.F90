
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
! ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
! :: ----------------------------------------------------------

  subroutine FORT_SETVOL(reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
       vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3,offset,dx,coord) &
       bind(c,name='amrex_setvol')

    implicit none

    integer    reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
    integer    vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
    integer    coord
    real(amrex_real)     dx(3), offset(3)
    real(amrex_real)     vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
       
    integer    i, j, k
    real(amrex_real)     v
       
    if (coord .eq. 0) then

       ! cartesian

       v = dx(1)*dx(2)*dx(3)
       do k = reg_l3, reg_h3
          do j = reg_l2, reg_h2
             do i = reg_l1, reg_h1
                vol(i,j,k) = v
             end do
          end do
       end do

    else

       write(6,*) "FORT_SETVOLUME not define for coord = ",coord

       call bl_abort(" ")

    end if
       
  end subroutine FORT_SETVOL

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

  subroutine FORT_SETAREA(reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
       area,area_l1,area_l2,area_l3,area_h1,area_h2,area_h3,offset,dx,dir,coord) &
       bind(c,name='amrex_setarea')

    implicit none

    integer    reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
    integer    area_l1,area_l2,area_l3,area_h1,area_h2,area_h3
    integer    coord, dir
    real(amrex_real)     dx(3), offset(3)
    real(amrex_real)     area(area_l1:area_h1,area_l2:area_h2,area_l3:area_h3)
       
    integer    i, j, k
    real(amrex_real)     fa

    fa = 1.d0
       
    if (coord .eq. 0) then

       if (dir .eq. 0) then
          fa = dx(2)*dx(3)
       else if (dir .eq. 1) then
          fa = dx(1)*dx(3)
       else if (dir .eq. 2) then
          fa = dx(1)*dx(2)
       else
          write(6,*) "FORT_SETAREA: invalid dir = ",dir
          call bl_abort(" ")
       end if

       do k = reg_l3, reg_h3
          do j = reg_l2, reg_h2
             do i = reg_l1, reg_h1
                area(i,j,k) = fa
             end do
          end do
       end do

    else

       write(6,*) "FORT_SETAREA not define for coord = ",coord
       call bl_abort(" ")

    end if
       
  end subroutine FORT_SETAREA

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


  subroutine FORT_SETDLOGA(dloga,dloga_l1,dloga_l2,dloga_l3,dloga_h1,dloga_h2,dloga_h3,offset,dx,dir,coord) &
       bind(c,name='amrex_setdloga')

    implicit none

    integer    dloga_l1,dloga_l2,dloga_l3,dloga_h1,dloga_h2,dloga_h3
    integer    coord
    real(amrex_real)     dx(3), offset(3)
    real(amrex_real)     dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2,dloga_l3:dloga_h3)
    integer dir
       
    integer    i, j, k
       
    if (coord .eq. 0) then

       ! cartesian

       do k = dloga_l3, dloga_h3
          do j = dloga_l2, dloga_h2
             do i = dloga_l1, dloga_h1
                dloga(i,j,k) = zero
             end do
          end do
       enddo

    else 

       write(6,*)' non-cartesian not allowed in 3D yet'
       call bl_abort(" ")

    endif

  end subroutine FORT_SETDLOGA

end module amrex_coordsys_module
