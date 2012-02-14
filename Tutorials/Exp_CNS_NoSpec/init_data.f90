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

    do i=1,nboxes(data)
       if ( multifab_remote(data,i) ) cycle

       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))

       select case(data%dim)
       case (2)
          call bl_error('We only support 3-D')
       case (3)
          call init_data_3d(lo,hi,ng,dx,dp)
       end select
    end do

  end subroutine init_data

  subroutine init_data_3d(lo,hi,ng,dx,cons)

    use advance_module ! For: irho imx imy imz iene.

    integer,          intent(in   ) :: lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)

    double precision :: xloc,yloc,zloc,rholoc,eloc,uvel,vvel,wvel,scale

    integer :: i,j,k

    scale = 0.02d0

    do k=lo(3)-ng,hi(3)+ng
       zloc = dfloat(k)*dx(3)
       do j=lo(2)-ng,hi(2)+ng
          yloc = dfloat(j)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             xloc = dfloat(i)*dx(1)

             uvel   = sin(xloc/scale)*sin(2*yloc/scale)*sin(3.d0*zloc/scale)
             vvel   = sin(2*xloc/scale)*sin(4.d0*yloc/scale)*sin(1.d0*zloc/scale)
             wvel   = sin(3.d0*xloc/scale)*cos(2*yloc/scale)*sin(2*zloc/scale)
             rholoc = 1.0d-3+1.0d-5*sin(xloc/scale)*cos(2*yloc/scale)*cos(3*zloc/scale)
             eloc   = 2.5d9+1.0d-3*sin(2*xloc/scale)*cos(2*yloc/scale)*sin(2*zloc/scale)

             cons(i,j,k,irho) = rholoc
             cons(i,j,k,imx)  = rholoc*uvel
             cons(i,j,k,imy)  = rholoc*vvel
             cons(i,j,k,imz)  = rholoc*wvel
             cons(i,j,k,iene) = rholoc*(eloc+(uvel**2+vvel**2+wvel**2)/2)

          enddo
       enddo
    enddo

  end subroutine init_data_3d

end module init_data_module
