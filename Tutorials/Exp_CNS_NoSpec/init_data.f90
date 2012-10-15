module init_data_module

  use bl_error_module
  use multifab_module
  use layout_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(data,dx,plo,phi)

    type(multifab),   intent(inout) :: data
    double precision, intent(in   ) :: dx(data%dim)
    double precision, intent(in   ) :: plo(data%dim), phi(data%dim)

    integer                   :: lo(data%dim), hi(data%dim), ng, i
    double precision, pointer :: dp(:,:,:,:)

    ng = data%ng

    do i=1,nfabs(data)

       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))

       select case(data%dim)
       case (2)
          call bl_error('We only support 3-D')
       case (3)
          call init_data_3d(lo,hi,ng,dx,dp,plo,phi)
       end select
    end do

  end subroutine init_data

  subroutine init_data_3d(lo,hi,ng,dx,cons,plo,phi)

    use advance_module ! For: irho imx imy imz iene.

    integer,          intent(in   ) :: lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3),plo(3),phi(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)

    integer          :: i,j,k
    double precision :: xloc,yloc,zloc,rholoc,eloc,uvel,vvel,wvel,scale(3)

    double precision, parameter :: twopi = 2.0d0 * 3.141592653589793238462643383279502884197d0

    do i=1,3
       scale(i) = (phi(i)-plo(i))/twopi
    end do

    !$OMP PARALLEL DO PRIVATE(i,j,k,zloc,yloc,xloc,uvel,vvel,wvel,rholoc,eloc)
    do k=lo(3),hi(3)
       zloc = dfloat(k)*dx(3)/scale(3)
       do j=lo(2),hi(2)
          yloc = dfloat(j)*dx(2)/scale(2)
          do i=lo(1),hi(1)
             xloc = dfloat(i)*dx(1)/scale(1)

             uvel   = 1.1d4*sin(1*xloc)*sin(2*yloc)*sin(3*zloc)
             vvel   = 1.0d4*sin(2*xloc)*sin(4*yloc)*sin(1*zloc)
             wvel   = 1.2d4*sin(3*xloc)*cos(2*yloc)*sin(2*zloc)
             rholoc = 1.0d-3 + 1.0d-5*sin(1*xloc)*cos(2*yloc)*cos(3*zloc)
             eloc   = 2.5d9  + 1.0d-3*sin(2*xloc)*cos(2*yloc)*sin(2*zloc)

             cons(i,j,k,irho) = rholoc
             cons(i,j,k,imx)  = rholoc*uvel
             cons(i,j,k,imy)  = rholoc*vvel
             cons(i,j,k,imz)  = rholoc*wvel
             cons(i,j,k,iene) = rholoc*(eloc + (uvel**2+vvel**2+wvel**2)/2)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_data_3d

end module init_data_module
