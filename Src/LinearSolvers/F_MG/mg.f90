module mg_defect_module

  use bl_types
  use stencil_module
  use stencil_nodal_module

  implicit none

  interface grid_laplace
     module procedure grid_laplace_1d
     module procedure grid_laplace_2d
     module procedure grid_laplace_3d
  end interface

  real(dp_t), private, parameter :: ZERO   = 0.0_dp_t
  real(dp_t), private, parameter :: EIGHTH = 0.125_dp_t
  real(dp_t), private, parameter :: FOURTH = 0.25_dp_t
  real(dp_t), private, parameter :: HALF   = 0.5_dp_t
  real(dp_t), private, parameter :: ONE    = 1.0_dp_t
  real(dp_t), private, parameter :: TWO    = 2.0_dp_t
  real(dp_t), private, parameter :: FOUR   = 4.0_dp_t

contains

  subroutine grid_laplace_1d(ss, dd, ff, uu, mm, ng, nodal_ng)
    integer, intent(in) :: ng
    integer,            intent(in)  :: nodal_ng
    real (kind = dp_t), intent(in)  :: ff(1-nodal_ng:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    real (kind = dp_t), intent(in)  :: ss(:,0:)
    real (kind = dp_t), intent(out) :: dd(1-nodal_ng:)
    integer,            intent(in)  :: mm(:)

    integer :: i,nx
    nx = size(ss,dim=1)-1

    i = 1
    dd(i) = HALF*ff(i) - (ss(i,0)*uu(i) + ss(i,1)*(uu(i+1)-uu(i)))

    do i = 2,nx
      dd(i) = ff(i) - &
              (ss(i,0)*uu(i) + ss(i,1) * uu(i+1) + ss(i,2) * uu(i-1))
    end do

    i = nx+1
    dd(i) = HALF*ff(i) - (ss(i,2)*(uu(i-1)-uu(i)))

  end subroutine grid_laplace_1d

  subroutine grid_laplace_2d(ss, dd, ff, uu, mm, ng, nodal_ng, face_type)
    integer, intent(in) :: ng
    integer,            intent(in   ) :: nodal_ng
    real (kind = dp_t), intent(in   ) :: ff(0:,0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(  out) :: dd(0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,0:)
    integer,            intent(in   ) :: mm(:,:)
    integer,            intent(in ) :: face_type(:,:)
    integer :: i,j,nx,ny,lo(2)
    integer :: istart,iend,jstart,jend

    nx = size(ss,dim=1)-1
    ny = size(ss,dim=2)-1

    lo = 1
    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    if (face_type(1,1) .eq. BC_NEU) then
      istart = 1
    else
      istart = 2
    end if
    if (face_type(1,2) .eq. BC_NEU) then
      iend = nx+1
    else
      iend = nx
    end if
    if (face_type(2,1) .eq. BC_NEU) then
      jstart = 1
    else
      jstart = 2
    end if
    if (face_type(2,2) .eq. BC_NEU) then
      jend = ny+1
    else
      jend = ny
    end if

!   Corners
    i = 1
    j = 1
    dd(i,j) = ss(i,j,8)*(uu(i+1,j+1)+HALF*uu(i,j+1)+HALF*uu(i+1,j)-TWO*uu(i,j))
    dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

    i = 1
    j = ny+1
    dd(i,j) = ss(i,j,3)*(uu(i+1,j-1)+HALF*uu(i,j-1)+HALF*uu(i+1,j)-TWO*uu(i,j))
    dd(i,j) = FOURTH*ff(i,j) - dd(i,j)
 
    i = nx+1
    j = 1
    dd(i,j) = ss(i,j,6)*(uu(i-1,j+1)+HALF*uu(i,j+1)+HALF*uu(i-1,j)-TWO*uu(i,j))
    dd(i,j) = FOURTH*ff(i,j) - dd(i,j)

    i = nx+1
    j = ny+1
    dd(i,j) = ss(i,j,1)*(uu(i-1,j-1)+HALF*uu(i,j-1)+HALF*uu(i-1,j)-TWO*uu(i,j))
    dd(i,j) = FOURTH*ff(i,j) - dd(i,j)
 
!   Lo-x edge
    i = 1
    do j = jstart,jend
       dd(i,j) = ss(i,j,3)*(uu(i+1,j-1)+HALF*uu(i,j-1)+HALF*uu(i+1,j)-TWO*uu(i,j)) &
                +ss(i,j,8)*(uu(i+1,j+1)+HALF*uu(i,j+1)+HALF*uu(i+1,j)-TWO*uu(i,j))
       dd(i,j) = HALF*ff(i,j) - dd(i,j)
    end do

!   Hi-x edge
    i = nx+1
    do j = jstart,jend
       dd(i,j) = ss(i,j,1)*(uu(i-1,j-1)+HALF*uu(i,j-1)+HALF*uu(i-1,j)-TWO*uu(i,j)) &
                +ss(i,j,6)*(uu(i-1,j+1)+HALF*uu(i,j+1)+HALF*uu(i-1,j)-TWO*uu(i,j))
       dd(i,j) = HALF*ff(i,j) - dd(i,j)
    end do
 
!   Lo-y edge
    j = 1
    do i = istart,iend
       dd(i,j) = ss(i,j,6)*(uu(i-1,j+1)+HALF*uu(i,j+1)+HALF*uu(i-1,j)-TWO*uu(i,j)) &
                +ss(i,j,8)*(uu(i+1,j+1)+HALF*uu(i,j+1)+HALF*uu(i+1,j)-TWO*uu(i,j))
       dd(i,j) = HALF*ff(i,j) - dd(i,j)
    end do

!   Hi-y edge
    j = ny+1
    do i = istart,iend
       dd(i,j) = ss(i,j,1)*(uu(i-1,j-1)+HALF*uu(i,j-1)+HALF*uu(i-1,j)-TWO*uu(i,j)) &
                +ss(i,j,3)*(uu(i+1,j-1)+HALF*uu(i,j-1)+HALF*uu(i+1,j)-TWO*uu(i,j))
       dd(i,j) = HALF*ff(i,j) - dd(i,j)
    end do
 
!   Interior
    do j = jstart,jend
    do i = istart,iend
       dd(i,j) = ss(i,j,0)*uu(i,j) + ss(i,j,1) * uu(i-1,j-1) &
                                   + ss(i,j,2) * uu(i  ,j-1) &
                                   + ss(i,j,3) * uu(i+1,j-1) &
                                   + ss(i,j,4) * uu(i-1,j  ) &
                                   + ss(i,j,5) * uu(i+1,j  ) &
                                   + ss(i,j,6) * uu(i-1,j+1) &
                                   + ss(i,j,7) * uu(i  ,j+1) &
                                   + ss(i,j,8) * uu(i+1,j+1)
       dd(i,j) = ff(i,j) - dd(i,j)
    end do
    end do


  end subroutine grid_laplace_2d

  subroutine grid_laplace_3d(ss, dd, ff, uu, mm, ng, nodal_ng, face_type)
    integer, intent(in) :: ng
    integer,            intent(in   ) :: nodal_ng
    real (kind = dp_t), intent(in   ) :: ff(0:,0:,0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(  out) :: dd(0:,0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(:,:,:,0:)
    integer,            intent(in   ) :: mm(:,:,:)
    integer,            intent(in   ) :: face_type(:,:)

    integer :: i,j,k,lo(3)
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: nx,ny,nz

    nx = size(ss,dim=1)-1
    ny = size(ss,dim=2)-1
    nz = size(ss,dim=3)-1

    lo = 1
    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    if (face_type(1,1) .eq. BC_NEU) then
      istart = 1
    else
      istart = 2
    end if
    if (face_type(1,2) .eq. BC_NEU) then
      iend = nx+1
    else
      iend = nx
    end if
    if (face_type(2,1) .eq. BC_NEU) then
      jstart = 1
    else
      jstart = 2
    end if
    if (face_type(2,2) .eq. BC_NEU) then
      jend = ny+1
    else
      jend = ny
    end if
    if (face_type(3,1) .eq. BC_NEU) then
      kstart = 1
    else
      kstart = 2
    end if
    if (face_type(3,2) .eq. BC_NEU) then
      kend = nz+1
    else
      kend = nz
    end if

!   Corners
    i = 1
    j = 1
    k = 1
    dd(i,j,k) = ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                             +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

    i = 1
    j = ny+1
    k = 1
    dd(i,j,k) = ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                             +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
 
    i = nx+1
    j = 1
    k = 1
    dd(i,j,k) = ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                             +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

    i = nx+1
    j = ny+1
    k = 1
    dd(i,j,k) = ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i-1,j-1,k  ) &
                             +uu(i-1,j  ,k+1) + uu(i  ,j-1,k+1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

    i = 1
    j = 1
    k = nz+1
    dd(i,j,k) = ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                             +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

    i = 1
    j = ny+1
    k = nz+1
    dd(i,j,k) = ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                             +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
 
    i = nx+1
    j = 1
    k = nz+1
    dd(i,j,k) = ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                             +uu(i-1,j  ,k-1) + uu(i  ,j+1,k-1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)

    i = nx+1
    j = ny+1
    k = nz+1
    dd(i,j,k) = ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                             +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                       - FOUR*uu(i  ,j  ,k) )
    dd(i,j,k) = EIGHTH*ff(i,j,k) - dd(i,j,k)
 
 
!   Lo-x / Lo-y edge
    i = 1
    j = 1
    do k = kstart,kend
       dd(i,j,k) = ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Hi-x / Lo-y edge
    i = nx+1
    j = 1
    do k = kstart,kend
       dd(i,j,k) = ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                 +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do
 
!   Lo-x / Hi-y edge
    i = 1
    j = ny+1
    do k = kstart,kend
       dd(i,j,k) = ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                                +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Hi-x / Hi-y edge
    i = nx+1
    j = ny+1
    do k = kstart,kend
       dd(i,j,k) = ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k+1) + uu(i  ,j-1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do
 
!   Lo-x / Lo-z edge
    i = 1
    k = 1
    do j = jstart,jend
       dd(i,j,k) = ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Hi-x / Lo-z edge
    i = nx+1
    k = 1
    do j = jstart,jend
       dd(i,j,k) = ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i-1,j  ,k+1) &
                                +uu(i-1,j-1,k  ) + uu(i  ,j-1,k+1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do
 
!   Lo-x / Hi-z edge
    i = 1
    k = nz+1
    do j = jstart,jend
       dd(i,j,k) = ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j  ,k-1) &
                                +uu(i+1,j+1,k  ) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                                +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Hi-x / Hi-z edge
    i = nx+1
    k = nz+1
    do j = jstart,jend
       dd(i,j,k) = ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j  ,k-1) &
                                +uu(i-1,j+1,k  ) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do
 
!   Lo-y / Lo-z edge
    j = 1
    k = 1
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                +uu(i  ,j+1,k+1) + uu(i-1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Hi-y / Lo-z edge
    j = ny+1
    k = 1
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                +uu(i  ,j-1,k+1) + uu(i+1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i  ,j-1,k+1) &
                                +uu(i-1,j-1,k  ) + uu(i-1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do
 
!   Lo-y / Hi-z edge
    j = 1
    k = nz+1
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j  ,k-1) &
                                +uu(i+1,j+1,k  ) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                +uu(i  ,j+1,k-1) + uu(i-1,j  ,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Hi-y / Hi-z edge
    j = ny+1
    k = nz+1
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i  ,j-1,k-1) &
                                +uu(i+1,j-1,k  ) + uu(i+1,j  ,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = FOURTH*ff(i,j,k) - dd(i,j,k)
    end do

!   Lo-x face
    i = 1
    do k = kstart,kend
    do j = jstart,jend
       dd(i,j,k) = ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j-1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i+1,j-1,k  ) &
                                +uu(i+1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
    end do
    end do
 
!   Hi-x face
    i = nx+1
    do k = kstart,kend
    do j = jstart,jend
       dd(i,j,k) = ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i-1,j  ,k+1) &
                                +uu(i-1,j-1,k  ) + uu(i  ,j-1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j  ,k-1) &
                                +uu(i-1,j+1,k  ) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
    end do
    end do

!   Lo-y face
    j = 1
    do k = kstart,kend
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k-1) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                +uu(i-1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
    end do
    end do

!   Hi-y face
    j = ny+1
    do k = kstart,kend
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i  ,j-1,k-1) &
                                +uu(i+1,j-1,k  ) + uu(i+1,j  ,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                +uu(i  ,j-1,k+1) + uu(i+1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i  ,j-1,k+1) &
                                +uu(i-1,j-1,k  ) + uu(i-1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
    end do
    end do

!   Lo-z face
    k = 1
    do j = jstart,jend
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k,15)*(uu(i+1,j-1,k+1) + uu(i+1,j-1,k  ) &
                                +uu(i  ,j-1,k+1) + uu(i+1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,13)*(uu(i-1,j-1,k+1) + uu(i  ,j-1,k+1) &
                                +uu(i-1,j-1,k  ) + uu(i-1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,20)*(uu(i+1,j+1,k+1) + uu(i+1,j+1,k  ) &
                                +uu(i+1,j  ,k+1) + uu(i  ,j+1,k+1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k,18)*(uu(i-1,j+1,k+1) + uu(i-1,j+1,k  ) &
                                +uu(i  ,j+1,k+1) + uu(i-1,j  ,k+1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
    end do
    end do

!   Hi-z face
    k = nz+1
    do j = jstart,jend
    do i = istart,iend
       dd(i,j,k) = ss(i,j,k, 3)*(uu(i+1,j-1,k-1) + uu(i  ,j-1,k-1) &
                                +uu(i+1,j-1,k  ) + uu(i+1,j  ,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 1)*(uu(i-1,j-1,k-1) + uu(i-1,j-1,k  ) &
                                +uu(i-1,j  ,k-1) + uu(i  ,j-1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 8)*(uu(i+1,j+1,k-1) + uu(i+1,j  ,k-1) &
                                +uu(i+1,j+1,k  ) + uu(i  ,j+1,k-1) &
                          - FOUR*uu(i  ,j  ,k) ) &
                  +ss(i,j,k, 6)*(uu(i-1,j+1,k-1) + uu(i-1,j+1,k  ) &
                                +uu(i  ,j+1,k-1) + uu(i-1,j  ,k-1) &
                          - FOUR*uu(i  ,j  ,k) )
       dd(i,j,k) = HALF*ff(i,j,k) - dd(i,j,k)
    end do
    end do

!   Interior
    do k = kstart,kend
    do j = jstart,jend
    do i = istart,iend

        dd(i,j,k) = ss(i,j,k,0)*uu(i,j,k) &
          + ss(i,j,k, 1) * uu(i-1,j-1,k-1) + ss(i,j,k, 2) * uu(i  ,j-1,k-1) &
          + ss(i,j,k, 3) * uu(i+1,j-1,k-1) + ss(i,j,k, 4) * uu(i-1,j  ,k-1) &
          + ss(i,j,k, 5) * uu(i+1,j  ,k-1) + ss(i,j,k, 6) * uu(i-1,j+1,k-1) &
          + ss(i,j,k, 7) * uu(i  ,j+1,k-1) + ss(i,j,k, 8) * uu(i+1,j+1,k-1) &
          + ss(i,j,k, 9) * uu(i-1,j-1,k  ) + ss(i,j,k,10) * uu(i+1,j-1,k  ) &
          + ss(i,j,k,11) * uu(i-1,j+1,k  ) + ss(i,j,k,12) * uu(i+1,j+1,k  ) &
          + ss(i,j,k,13) * uu(i-1,j-1,k+1) + ss(i,j,k,14) * uu(i  ,j-1,k+1) &
          + ss(i,j,k,15) * uu(i+1,j-1,k+1) + ss(i,j,k,16) * uu(i-1,j  ,k+1) &
          + ss(i,j,k,17) * uu(i+1,j  ,k+1) + ss(i,j,k,18) * uu(i-1,j+1,k+1) &
          + ss(i,j,k,19) * uu(i  ,j+1,k+1) + ss(i,j,k,20) * uu(i+1,j+1,k+1)

!       Add faces (only non-zero for non-uniform dx)
        dd(i,j,k) = dd(i,j,k) + &
            ss(i,j,k,21) * uu(i-1,j  ,k  ) + ss(i,j,k,22) * uu(i+1,j  ,k  ) &
          + ss(i,j,k,23) * uu(i  ,j-1,k  ) + ss(i,j,k,24) * uu(i  ,j+1,k  ) &
          + ss(i,j,k,25) * uu(i  ,j  ,k-1) + ss(i,j,k,26) * uu(i  ,j  ,k+1)

        dd(i,j,k) = ff(i,j,k) - dd(i,j,k)
    end do
    end do
    end do
  end subroutine grid_laplace_3d

end module mg_defect_module

module mg_prolongation_module

  use bl_types

  implicit none

  interface pc_c_prolongation
     module procedure pc_c_prolongation_1d
     module procedure pc_c_prolongation_2d
     module procedure pc_c_prolongation_3d
  end interface

  interface nodal_prolongation
     module procedure nodal_prolongation_1d
     module procedure nodal_prolongation_2d
     module procedure nodal_prolongation_3d
  end interface

contains

  subroutine pc_c_prolongation_1d(ff, cc, ir)
    real (kind = dp_t), intent(inout) :: ff(0:)
    real (kind = dp_t), intent(in) :: cc(0:)
    integer, intent(in) :: ir(:)
    integer nx, i, l
    nx = size(cc,dim=1)
    do l = 0, ir(1)-1
       do i = 0, nx - 1
          ff(ir(1)*i+l) = ff(ir(1)*i+l) + cc(i)
       end do
    end do
  end subroutine pc_c_prolongation_1d

  subroutine pc_c_prolongation_2d(ff, cc, ir)
    real (kind = dp_t), intent(inout) :: ff(0:,0:)
    real (kind = dp_t), intent(in) :: cc(0:,0:)
    integer, intent(in) :: ir(:)
    integer nx, ny, i, j, l, m
    nx = size(cc,dim=1)
    ny = size(cc,dim=2)
    do m = 0, ir(2)-1
       do l = 0, ir(1)-1
          !$OMP PARALLEL DO PRIVATE(j,i)
          do j = 0, ny - 1
             do i = 0, nx - 1
                ff(ir(1)*i+l, ir(2)*j+m) = ff(ir(1)*i+l, ir(2)*j+m) + cc(i,j)
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do
  end subroutine pc_c_prolongation_2d

  subroutine pc_c_prolongation_3d(ff, cc, ir)
    real (kind = dp_t), intent(inout) :: ff(0:,0:,0:)
    real (kind = dp_t), intent(in) :: cc(0:,0:,0:)
    integer, intent(in) :: ir(:)
    integer nx, ny, nz, i, j, k, l, m, n
    nx = size(cc,dim=1)
    ny = size(cc,dim=2)
    nz = size(cc,dim=3)
    do n = 0, ir(3)-1
       do m = 0, ir(2)-1
          do l = 0, ir(1)-1
             !$OMP PARALLEL DO PRIVATE(j,i,k)
             do k = 0, nz - 1
                do j = 0, ny - 1
                   do i = 0, nx - 1
                      ff(ir(1)*i+l, ir(2)*j+m, ir(3)*k+n) = ff(ir(1)*i+l, ir(2)*j+m, ir(3)*k+n) + cc(i,j,k)
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end do
       end do
    end do
  end subroutine pc_c_prolongation_3d

  subroutine nodal_prolongation_1d(ff, cc, ir)
    real (kind = dp_t), intent(inout) :: ff(0:)
    real (kind = dp_t), intent(in) :: cc(0:)
    integer, intent(in) :: ir(:)
    integer nx, i, l
    real (kind = dp_t) :: fac_left, fac_rght

    nx = size(cc,dim=1)-1

    do i = 0, nx
      ff(ir(1)*i) = ff(ir(1)*i) + cc(i)
    end do

    do l = 1, ir(1)-1
      fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
      fac_rght = 1.0_dp_t - fac_left
      do i = 0, nx - 1
        ff(ir(1)*i+l) = ff(ir(1)*i+l) + fac_left*cc(i) + fac_rght*cc(i+1)
      end do
    end do

  end subroutine nodal_prolongation_1d

  subroutine nodal_prolongation_2d(ff, cc, ir)
    real (kind = dp_t), intent(inout) :: ff(0:,0:)
    real (kind = dp_t), intent(inout) :: cc(0:,0:)
    integer, intent(in) :: ir(:)
    integer nx, ny, i, j, l, m

    real (kind = dp_t) :: fac_left, fac_rght
    real (kind = dp_t) :: temp(0:size(ff,dim=1)-1,0:size(ff,dim=2)-1)

    nx = size(cc,dim=1)-1
    ny = size(cc,dim=2)-1

!   Interpolate at fine nodes on top of coarse nodes
    do j = 0,ny
      do i = 0,nx
        temp(ir(1)*i, ir(2)*j) = cc(i,j)
      end do
    end do

!   Interpolate at fine nodes between coarse nodes in the i-direction only
    do j = 0,ny
      do l = 1, ir(1)-1
        fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
        fac_rght = 1.0_dp_t - fac_left
        do i = 0,nx-1
          temp(ir(1)*i+l, ir(2)*j) = fac_left*cc(i,j) + fac_rght*cc(i+1,j)
        end do
      end do
    end do

!   Interpolate in the j-direction using previously interpolated "temp"
    do m = 1, ir(2)-1
      fac_left = real(m,kind=dp_t) / real(ir(2),kind=dp_t)
      fac_rght = 1.0_dp_t - fac_left
      do j = 0,ny-1
      do i = 0,ir(1)*nx
          temp(i, ir(2)*j+m) = fac_left*temp(i,ir(2)*(j  )) &
                              +fac_rght*temp(i,ir(2)*(j+1))
      end do
      end do
    end do

    do j = 0,ir(2)*ny
    do i = 0,ir(1)*nx
      ff(i,j) = ff(i,j) + temp(i,j)
    end do
    end do

  end subroutine nodal_prolongation_2d

  subroutine nodal_prolongation_3d(ff, cc, ir)
    real (kind = dp_t), intent(inout) :: ff(0:,0:,0:)
    real (kind = dp_t), intent(in) :: cc(0:,0:,0:)
    integer, intent(in) :: ir(:)
    integer nx, ny, nz, i, j, k, l, m, n

    real (kind = dp_t) :: fac_left, fac_rght
    real (kind = dp_t) :: temp(0:size(ff,dim=1)-1,0:size(ff,dim=2)-1,0:size(ff,dim=3)-1)

    nx = size(cc,dim=1)-1
    ny = size(cc,dim=2)-1
    nz = size(cc,dim=3)-1

!   Interpolate at fine nodes on top of coarse nodes
    do k = 0,nz
    do j = 0,ny
    do i = 0,nx
      temp(ir(1)*i,ir(2)*j,ir(3)*k) = cc(i,j,k)
    end do
    end do
    end do

!   Interpolate at fine nodes between coarse nodes in the i-direction only
    do k = 0,nz
    do j = 0,ny
      do l = 1, ir(1)-1
        fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
        fac_rght = 1.0_dp_t - fac_left
        do i = 0,nx-1
          temp(ir(1)*i+l,ir(2)*j,ir(3)*k) = fac_left*cc(i  ,j, k) &
                                           +fac_rght*cc(i+1,j, k)
        end do
      end do
    end do
    end do

!   Interpolate in the j-direction using previously interpolated "temp"
    do m = 1, ir(2)-1
      fac_left = real(m,kind=dp_t) / real(ir(2),kind=dp_t)
      fac_rght = 1.0_dp_t - fac_left
      do k = 0,nz
      do j = 0,ny-1
      do i = 0,ir(1)*nx
          temp(i,ir(2)*j+m,ir(3)*k) = fac_left*temp(i,ir(2)*(j  ),ir(3)*k) &
                                     +fac_rght*temp(i,ir(2)*(j+1),ir(3)*k)
      end do
      end do
      end do
    end do

!   Interpolate in the k-direction using previously interpolated "temp"
    do n = 1, ir(3)-1
      fac_left = real(n,kind=dp_t) / real(ir(3),kind=dp_t)
      fac_rght = 1.0_dp_t - fac_left
      do k = 0,nz-1
      do j = 0,ir(2)*ny
      do i = 0,ir(1)*nx
          temp(i, j, ir(3)*k+n) = fac_left*temp(i,j,ir(3)*(k  )) &
                                 +fac_rght*temp(i,j,ir(3)*(k+1))
      end do
      end do
      end do
    end do

    do k = 0,ir(3)*nz
    do j = 0,ir(2)*ny
    do i = 0,ir(1)*nx
      ff(i,j,k) = ff(i,j,k) + temp(i,j,k)
    end do
    end do
    end do

  end subroutine nodal_prolongation_3d

end module mg_prolongation_module

module mg_restriction_module

  use bl_types
  use stencil_module
  use stencil_nodal_module

  implicit none

  interface cc_restriction
     module procedure cc_restriction_1d
     module procedure cc_restriction_2d
     module procedure cc_restriction_3d
  end interface

  interface nodal_restriction
     module procedure nodal_restriction_1d
     module procedure nodal_restriction_2d
     module procedure nodal_restriction_3d
  end interface

  real(kind=dp_t), private, parameter :: ZERO = 0.0_dp_t
  real(kind=dp_t), private, parameter :: ONE  = 1.0_dp_t
  real(kind=dp_t), private, parameter :: TWO  = 2.0_dp_t
  real(kind=dp_t), private, parameter :: HALF = 0.5_dp_t

contains
  subroutine cc_restriction_1d(cc, loc, ff, lof, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(out) :: cc(loc(1):)
    real (kind = dp_t), intent(in) :: ff(lof(1):)
    integer, intent(in) :: ir(:)
    real (kind = dp_t) :: fac
    integer :: i, l

    fac = one/real(product(ir),kind=dp_t)

    do i = lo(1),hi(1)
       cc(i) = zero
       do l = 0, ir(1)-1
          cc(i) = cc(i) + ff(ir(1)*i+l)
       end do
       cc(i) = cc(i)*fac
    end do

  end subroutine cc_restriction_1d

  subroutine cc_restriction_2d(cc, loc, ff, lof, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(out) :: cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in) :: ff(lof(1):,lof(2):)
    integer, intent(in) :: ir(:)
    real (kind = dp_t) :: fac
    integer :: i, j, l, m

    fac = one/real(product(ir),kind=dp_t)

    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
      cc(i,j) = zero
      do m = 0, ir(2)-1
         do l = 0, ir(1)-1
             cc(i,j) = cc(i,j) + ff(ir(1)*i+l,ir(2)*j+m)
          end do
      end do
      cc(i,j) = cc(i,j)*fac
    end do
    end do

  end subroutine cc_restriction_2d

  subroutine cc_restriction_3d(cc, loc, ff, lof, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:),hi(:)
    real (kind = dp_t), intent(out) :: cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in) :: ff(lof(1):,lof(2):,lof(3):)
    integer, intent(in) :: ir(:)
    real (kind = dp_t) :: fac
    integer :: i, j, k, l, m, n

    fac = one/real(product(ir),kind=dp_t)

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
      cc(i,j,k) = zero
      do n = 0, ir(3)-1
        do m = 0, ir(2)-1
          do l = 0, ir(1)-1
            cc(i,j,k) = cc(i,j,k) + ff(ir(1)*i+l,ir(2)*j+m,ir(3)*k+n)
          end do
        end do
      end do
      cc(i,j,k) = cc(i,j,k)*fac
    end do
    end do
    end do

  end subroutine cc_restriction_3d

  subroutine nodal_restriction_1d(cc, loc, ff, lof, &
                                  mm_fine, lom_fine, mm_crse, lom_crse, &
                                  lo, hi, ir, inject, mg_restriction_mode)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lom_fine(:)
    integer, intent(in) :: lom_crse(:)
    integer, intent(in) :: lo(:), hi(:)
    real(kind=dp_t), intent(out) :: cc(loc(1):)
    real(kind=dp_t), intent(in)  :: ff(lof(1):)
    integer        , intent(in)  :: mm_fine(lom_fine(1):)
    integer        , intent(in)  :: mm_crse(lom_crse(1):)
    integer :: ir(:)
    integer :: hif
    logical,intent(in) :: inject
    integer,intent(in) :: mg_restriction_mode

    integer :: i,ifine,m

    real(kind=dp_t) :: fac,fac0

    hif = lof(1)+size(ff,dim=1)-1

    if ( inject ) then

       do i = lo(1),hi(1)
          cc(i) = ff(ir(1)*i)
       end do

    else if (mg_restriction_mode == 1) then

       fac0 = 1.d0 
       do m = 0, ir(1)-1
         fac = (ir(1)-m) * fac0
         if (m .eq. 0) fac = HALF * fac
         do i = lo(1),hi(1)
           ifine = i*ir(1)
           if (.not.bc_dirichlet(mm_fine(ifine),1,0)) &
             cc(i) = cc(i) + fac * (ff(ifine-m) + ff(ifine+m))
         end do
       end do

    else 

       fac0 = 1.d0 
       do m = 0, ir(1)-1
         fac = (ir(1)-m) * fac0
         if (m .eq. 0) fac = HALF * fac
         do i = lo(1),hi(1)
           ifine = i*ir(1)
           if (.not.bc_dirichlet(mm_fine(ifine),1,0)) then
              if (ifine == lof(1)+1) then
                 cc(i) = cc(i) + fac * ff(ifine+m)
              else if (ifine == hif-1) then
                 cc(i) = cc(i) + fac * ff(ifine-m)
              else
                 cc(i) = cc(i) + fac * (ff(ifine-m) + ff(ifine+m))
              end if
           end if
         end do
       end do

    end if

    i = lo(1)
    if (bc_dirichlet(mm_crse(i),1,0)) then
      cc(i) = ZERO
    else if (bc_neumann(mm_crse(i),1,-1)) then
      cc(i) = TWO*cc(i)
    end if

    i = hi(1)
    if (bc_dirichlet(mm_crse(i),1,0)) then
      cc(i) = ZERO
    else if (bc_neumann(mm_crse(i),1,+1)) then
      cc(i) = TWO*cc(i)
    end if

  end subroutine nodal_restriction_1d

  subroutine nodal_restriction_2d(cc, loc, ff, lof, &
                                  mm_fine, lom_fine, mm_crse, lom_crse, &
                                  lo, hi, ir, inject, mg_restriction_mode)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lom_fine(:)
    integer, intent(in) :: lom_crse(:)
    integer, intent(in) :: lo(:), hi(:)
    real(kind=dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real(kind=dp_t), intent(  out) :: cc(loc(1):,loc(2):)
    integer        , intent(in)  :: mm_fine(lom_fine(1):,lom_fine(2):)
    integer        , intent(in)  :: mm_crse(lom_crse(1):,lom_crse(2):)
    integer :: ir(:)
    logical, intent(in) :: inject
    integer,intent(in) :: mg_restriction_mode

    integer :: i, j, ifine, jfine, m, n, ng
    integer :: hif(2)
    real(kind=dp_t) :: fac,fac0,fac1
    logical :: add_lo_x, add_lo_y, add_hi_x, add_hi_y

    hif(1) = lof(1)+size(ff,dim=1)-1
    hif(2) = lof(2)+size(ff,dim=2)-1
    
    ng = lom_fine(1) - lof(1)
    call impose_neumann_bcs_2d(ff,mm_fine,lom_fine,ng)

    if ( inject ) then

       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          cc(i,j) = ff(ir(1)*i,ir(2)*j)
       end do
       end do

    else if (mg_restriction_mode == 1) then

       ng = lom_fine(1) - lof(1)
       call impose_neumann_bcs_2d(ff,mm_fine,lom_fine,ng)

       fac0 = 1.0_dp_t / (ir(1)*ir(2))
       do n = 0, ir(2)-1
         fac1 = (ir(2)-n) * fac0
         if (n .eq. 0) fac1 = HALF * fac1
         do m = 0, ir(1)-1
            fac = (ir(1)-m) * fac1
            if (m .eq. 0) fac = HALF * fac
            do j = lo(2),hi(2)
               jfine = j*ir(2)
               do i = lo(1),hi(1)
                  ifine = i*ir(1)
                  if (.not. bc_dirichlet(mm_fine(ifine,jfine),1,0)) then
                    cc(i,j) = cc(i,j) &
                       + fac*ff(ifine-m,jfine-n) &
                       + fac*ff(ifine+m,jfine-n) &
                       + fac*ff(ifine-m,jfine+n) &
                       + fac*ff(ifine+m,jfine+n)
                  end if
               end do
            end do
         end do
       end do

    else

       ng = lom_fine(1) - lof(1)
       call impose_neumann_bcs_2d(ff,mm_fine,lom_fine,ng)
     
       fac0 = 1.0_dp_t / (ir(1)*ir(2))
       do n = 0, ir(2)-1
         fac1 = (ir(2)-n) * fac0
         if (n .eq. 0) fac1 = HALF * fac1
         do m = 0, ir(1)-1
            fac = (ir(1)-m) * fac1
            if (m .eq. 0) fac = HALF * fac
            do j = lo(2),hi(2)
              jfine = j*ir(2)
              do i = lo(1),hi(1)
               add_lo_x = .true.
               add_lo_y = .true.
               add_hi_x = .true.
               add_hi_y = .true.
               ifine = i*ir(1)
               if (.not. bc_dirichlet(mm_fine(ifine,jfine),1,0)) then
                 
                 if (ifine == lof(1)+1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine),1,-1)) add_lo_x = .false. 
                 if (jfine == lof(2)+1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine),2,-1)) add_lo_y = .false. 
                 if (ifine == hif(1)-1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine),1,+1)) add_hi_x = .false.
                 if (jfine == hif(2)-1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine),2,+1)) add_hi_y = .false.

                 if (add_lo_x .and. add_lo_y) cc(i,j) = cc(i,j) + fac*ff(ifine-m,jfine-n)
                 if (add_hi_x .and. add_lo_y) cc(i,j) = cc(i,j) + fac*ff(ifine+m,jfine-n)
                 if (add_lo_x .and. add_hi_y) cc(i,j) = cc(i,j) + fac*ff(ifine-m,jfine+n)
                 if (add_hi_x .and. add_hi_y) cc(i,j) = cc(i,j) + fac*ff(ifine+m,jfine+n)

               end if
            end do
           end do
         end do
       end do

    end if

    do j = lo(2),hi(2)
      if (bc_dirichlet(mm_crse(lo(1),j),1,0)) cc(lo(1),j) = ZERO
      if (bc_dirichlet(mm_crse(hi(1),j),1,0)) cc(hi(1),j) = ZERO
    end do

    do i = lo(1),hi(1)
      if (bc_dirichlet(mm_crse(i,lo(2)),1,0)) cc(i,lo(2)) = ZERO
      if (bc_dirichlet(mm_crse(i,hi(2)),1,0)) cc(i,hi(2)) = ZERO
    end do

  end subroutine nodal_restriction_2d

  subroutine nodal_restriction_3d(cc, loc, ff, lof, &
                                  mm_fine, lom_fine, mm_crse, lom_crse, &
                                  lo, hi, ir, inject, mg_restriction_mode)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lom_fine(:)
    integer, intent(in) :: lom_crse(:)
    integer, intent(in) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real(kind=dp_t), intent(  out) :: cc(loc(1):,loc(2):,loc(3):)
    integer        , intent(in   ) :: mm_fine(lom_fine(1):,lom_fine(2):,lom_fine(3):)
    integer        , intent(in   ) :: mm_crse(lom_crse(1):,lom_crse(2):,lom_crse(3):)
    integer :: ir(:)
    logical, intent(in) :: inject
    integer,intent(in) :: mg_restriction_mode

    integer :: i, j, k, l, m, n, ng
    integer :: ifine,jfine,kfine
    integer :: hif(3)
    real(kind=dp_t) :: fac,fac0,fac1,fac2
    logical :: add_lo_x, add_lo_y, add_lo_z, add_hi_x, add_hi_y, add_hi_z

    hif(1) = lof(1)+size(ff,dim=1)-1
    hif(2) = lof(2)+size(ff,dim=2)-1
    hif(3) = lof(3)+size(ff,dim=3)-1

    if ( inject ) then

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          cc(i,j,k) = ff(ir(1)*i,ir(2)*j,ir(3)*k)
       end do
       end do
       end do

    else if (mg_restriction_mode == 1) then

       ng = lom_fine(1) - lof(1)
       call impose_neumann_bcs_3d(ff,mm_fine,lom_fine,ng)

       fac0 = 1.0_dp_t / (ir(1)*ir(2)*ir(3))
       do l = 0, ir(3)-1
         fac2 = (ir(3)-l) * fac0
         if (l .eq. 0) fac2 = HALF * fac2
         do n = 0, ir(2)-1
           fac1 = (ir(2)-n) * fac2
           if (n .eq. 0) fac1 = HALF * fac1
           do m = 0, ir(1)-1
             fac = (ir(1)-m) * fac1
             if (m .eq. 0) fac = HALF * fac
             do k = lo(3),hi(3)
               kfine = k*ir(3)
               do j = lo(2),hi(2)
                 jfine = j*ir(2)
                 do i = lo(1),hi(1)
                   ifine = i*ir(1)
                   if (.not. bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0)) then
                     cc(i,j,k) = cc(i,j,k) &
                        + fac*ff(ifine-m,jfine-n,kfine-l) &
                        + fac*ff(ifine+m,jfine-n,kfine-l) &
                        + fac*ff(ifine-m,jfine+n,kfine-l) &
                        + fac*ff(ifine+m,jfine+n,kfine-l) &
                        + fac*ff(ifine-m,jfine-n,kfine+l) &
                        + fac*ff(ifine+m,jfine-n,kfine+l) &
                        + fac*ff(ifine-m,jfine+n,kfine+l) &
                        + fac*ff(ifine+m,jfine+n,kfine+l)
                   end if
                 end do
               end do
             end do
           end do
         end do
       end do

    else

       ng = lom_fine(1) - lof(1)
       call impose_neumann_bcs_3d(ff,mm_fine,lom_fine,ng)

       fac0 = 1.0_dp_t / (ir(1)*ir(2)*ir(3))
       do l = 0, ir(3)-1
         fac2 = (ir(3)-l) * fac0
         if (l .eq. 0) fac2 = HALF * fac2
         do n = 0, ir(2)-1
           fac1 = (ir(2)-n) * fac2
           if (n .eq. 0) fac1 = HALF * fac1
           do m = 0, ir(1)-1
             fac = (ir(1)-m) * fac1
             if (m .eq. 0) fac = HALF * fac
             do k = lo(3),hi(3)
               kfine = k*ir(3)
               do j = lo(2),hi(2)
                 jfine = j*ir(2)
                 do i = lo(1),hi(1)
                   ifine = i*ir(1)
                   add_lo_x = .true.
                   add_lo_y = .true.
                   add_lo_z = .true.
                   add_hi_x = .true.
                   add_hi_y = .true.
                   add_hi_z = .true.
                   if (.not. bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0)) then
                 
                    if (ifine == lof(1)+1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine,kfine),1,-1)) add_lo_x = .false. 
                    if (jfine == lof(2)+1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine,kfine),2,-1)) add_lo_y = .false. 
                    if (kfine == lof(3)+1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine,kfine),3,-1)) add_lo_z = .false. 
                    if (ifine == hif(1)-1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine,kfine),1,+1)) add_hi_x = .false.
                    if (jfine == hif(2)-1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine,kfine),2,+1)) add_hi_y = .false.
                    if (kfine == hif(3)-1 .and. &
                     .not. bc_neumann(mm_fine(ifine,jfine,kfine),3,+1)) add_hi_z = .false.

                    if (add_lo_z) then
                      if (add_lo_x .and. add_lo_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine-m,jfine-n,kfine-l)
                      if (add_hi_x .and. add_lo_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine+m,jfine-n,kfine-l)
                      if (add_lo_x .and. add_hi_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine-m,jfine+n,kfine-l)
                      if (add_hi_x .and. add_hi_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine+m,jfine+n,kfine-l)
                    end if
                    if (add_hi_z) then
                      if (add_lo_x .and. add_lo_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine-m,jfine-n,kfine+l)
                      if (add_hi_x .and. add_lo_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine+m,jfine-n,kfine+l)
                      if (add_lo_x .and. add_hi_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine-m,jfine+n,kfine+l)
                      if (add_hi_x .and. add_hi_y) &
                          cc(i,j,k) = cc(i,j,k) + fac*ff(ifine+m,jfine+n,kfine+l)
                    end if

                   end if
                 end do
               end do
             end do
           end do
         end do
       end do

    end if

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      if (bc_dirichlet(mm_crse(lo(1),j,k),1,0)) cc(lo(1),j,k) = ZERO
      if (bc_dirichlet(mm_crse(hi(1),j,k),1,0)) cc(hi(1),j,k) = ZERO
    end do
    end do

    do k = lo(3),hi(3)
    do i = lo(1),hi(1)
      if (bc_dirichlet(mm_crse(i,lo(2),k),1,0)) cc(i,lo(2),k) = ZERO
      if (bc_dirichlet(mm_crse(i,hi(2),k),1,0)) cc(i,hi(2),k) = ZERO
    end do
    end do

    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
      if (bc_dirichlet(mm_crse(i,j,lo(3)),1,0)) cc(i,j,lo(3)) = ZERO
      if (bc_dirichlet(mm_crse(i,j,hi(3)),1,0)) cc(i,j,hi(3)) = ZERO
    end do
    end do

  end subroutine nodal_restriction_3d

end module mg_restriction_module

module mg_module

  use multifab_module
  use stencil_module
  use stencil_nodal_module
  use itsol_module
  use sparse_solve_module
  use bl_timer_module

  implicit none

  integer, parameter :: MG_SMOOTHER_GS_RB  = 1
  integer, parameter :: MG_SMOOTHER_JACOBI = 2
  integer, parameter :: MG_SMOOTHER_GS_LEX = 3

  integer, parameter :: MG_FCycle = 1
  integer, parameter :: MG_WCycle = 2
  integer, parameter :: MG_VCycle = 3

  type mg_tower

     integer :: dim = 0

     ! defaults appropriate for abec, laplacian
     integer :: nc = 1
     integer :: ng = 1
     integer :: ns = 1

     ! defaults
     integer :: solver  = 1

     ! gsrb
     integer :: smoother = MG_SMOOTHER_GS_RB
     integer :: nu1 = 2
     integer :: nu2 = 2
     integer :: gamma = 1
     integer :: cycle = MG_Vcycle
     real(kind=dp_t) :: omega = 1.0_dp_t

     ! bottom solver defaults good for bicg
     integer :: bottom_solver = 1
     integer :: bottom_max_iter = 100
     real(kind=dp_t) :: bottom_solver_eps = 1.0e-4_dp_t

     integer :: nboxes =  0
     integer :: nlevels =  0

     ! let MG pick the maximum number of levels
     integer :: max_nlevel = Huge(1)
     integer :: min_width  = 1

     ! good for many problems
     integer :: max_iter = 20
     real(kind=dp_t) :: eps = 1.0e-10_dp_t

     type(box), pointer :: pd(:) => Null()
     real(kind=dp_t), pointer :: dh(:,:) => Null()

     type(multifab), pointer :: cc(:) => Null()
     type(multifab), pointer :: ff(:) => Null()
     type(multifab), pointer :: dd(:) => Null()
     type(multifab), pointer :: uu(:) => Null()
     type(multifab), pointer :: ss(:) => Null()
     type(imultifab), pointer :: mm(:) => Null()

     type(stencil), pointer :: st(:) => Null()

     integer, pointer :: face_type(:,:,:)

     type(timer), pointer :: tm(:) => Null()

     type(sparse) sparse_object
     type(multifab) :: rh1
     type(multifab) :: uu1
     type(multifab) :: ss1

     integer :: verbose = 0
  end type mg_tower

  real(kind=dp_t), parameter, private :: zero = 0.0_dp_t

contains

  subroutine mg_tower_build(mgt, la, pd, domain_bc, &
       solver, nu1, nu2, gamma, cycle, &
       smoother, omega, &
       dh, &
       ns, &
       nc, ng, &
       max_nlevel, min_width, &
       max_iter, eps, &
       bottom_solver, bottom_max_iter, bottom_solver_eps, &
       st_type, &
       verbose, nodal)

    type(mg_tower), intent(inout) :: mgt
    type(layout), intent(inout) :: la
    type(box), intent(in) :: pd
    integer, intent(in) :: domain_bc(:,:)

    integer, intent(in), optional :: ns
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: ng
    integer, intent(in), optional :: solver, nu1, nu2, gamma, cycle
    integer, intent(in), optional :: smoother
    logical, intent(in), optional :: nodal(:)
    real(dp_t), intent(in), optional :: omega
    real(kind=dp_t), intent(in), optional :: dh(:)
    real(kind=dp_t), intent(in), optional :: eps
    real(kind=dp_t), intent(in), optional :: bottom_solver_eps
    integer, intent(in), optional :: max_nlevel
    integer, intent(in), optional :: min_width
    integer, intent(in), optional :: max_iter
    integer, intent(in), optional :: bottom_solver
    integer, intent(in), optional :: bottom_max_iter
    integer, intent(in), optional :: verbose
    integer, intent(in), optional :: st_type

    integer :: lo_grid,hi_grid,lo_dom,hi_dom
    integer :: ng_for_res
    integer :: n, i, id
    type(layout) :: la1, la2
    logical :: nodal_flag

    ! Paste in optional arguments
    if ( present(ng)                ) mgt%ng                = ng
    if ( present(nc)                ) mgt%nc                = nc
    if ( present(max_nlevel)        ) mgt%max_nlevel        = max_nlevel
    if ( present(max_iter)          ) mgt%max_iter          = max_iter
    if ( present(eps)               ) mgt%eps               = eps
    if ( present(solver)            ) mgt%solver            = solver
    if ( present(smoother)          ) mgt%smoother          = smoother
    if ( present(nu1)               ) mgt%nu1               = nu1
    if ( present(nu2)               ) mgt%nu2               = nu2
    if ( present(gamma)             ) mgt%gamma             = gamma
    if ( present(omega)             ) mgt%omega             = omega
    if ( present(cycle)             ) mgt%cycle             = cycle
    if ( present(bottom_solver)     ) mgt%bottom_solver     = bottom_solver
    if ( present(bottom_solver_eps) ) mgt%bottom_solver_eps = bottom_solver_eps
    if ( present(bottom_max_iter)   ) mgt%bottom_max_iter   = bottom_max_iter
    if ( present(min_width)         ) mgt%min_width         = min_width
    if ( present(verbose)           ) mgt%verbose           = verbose

    nodal_flag = .false.
    if ( present(nodal) ) then
       if ( all(nodal) ) then
          nodal_flag = .true.
       else if ( any(nodal) ) then
          call bl_error("mixed nodal/not nodal forbidden")
       end if
    end if

    if ( present(ns) ) then
       mgt%ns = ns
    else
       if ( nodal_flag ) then
          mgt%ns = 3**mgt%dim
       else
          mgt%ns = 1 + 3*mgt%dim
       end if
    end if

    if ( .not. nodal_flag ) then
       if ( .not. present(omega) ) then
          select case ( mgt%dim )
          case (2)
             select case (smoother)
             case ( MG_SMOOTHER_JACOBI )
                mgt%omega = 4.0_dp_t/5.0_dp_t
             end select
          case (3)
             select case (smoother)
             case ( MG_SMOOTHER_JACOBI )
                mgt%omega = 6.0_dp_t/7.0_dp_t
             case ( MG_SMOOTHER_GS_RB )
                mgt%omega = 1.15_dp_t
             end select
          end select
       end if
    end if
    ng_for_res = 0; if ( nodal_flag ) ng_for_res = 1

    n = max_mg_levels(layout_boxarray(la), mgt%min_width)
    mgt%nlevels = min(n, mgt%max_nlevel)
    n = mgt%nlevels
    mgt%nboxes = layout_nboxes(la)
    mgt%dim    = layout_dim(la)

    allocate(mgt%face_type(mgt%nboxes,mgt%dim,2))

    allocate(mgt%cc(n), mgt%ff(n), mgt%dd(n), mgt%uu(n-1), mgt%ss(n), mgt%mm(n))
    allocate(mgt%pd(n),mgt%dh(mgt%dim,n))
    allocate(mgt%tm(n))
    if ( n == 1 ) then
       mgt%tm(n)%name = "SINGLE LEVEL"
    else
       mgt%tm(n)%name = "FINEST LEVEL"
       mgt%tm(1)%name = "COARSEST LEVEL"
    end if

    la1 = la
    do i = mgt%nlevels, 1, -1
       call  multifab_build(mgt%cc(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%ff(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%dd(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%ss(i), la1, ns, 0, nodal)

       call setval(mgt%cc(i), zero, all = .TRUE.)
       call setval(mgt%ff(i), zero, all = .TRUE.)
       call setval(mgt%dd(i), zero, all = .TRUE.)
       call setval(mgt%ss(i), zero, all = .TRUE.)

       call imultifab_build(mgt%mm(i), la1, 1, 0, nodal)
       if ( i /= mgt%nlevels ) then
          call multifab_build(mgt%uu(i), la1, mgt%nc, mgt%ng, nodal)
       end if
       mgt%pd(i) = coarsen(pd, 2**(mgt%nlevels-i))
       if ( i > 1 ) call layout_build_coarse(la2, la1, (/(2,i=1,mgt%dim)/))
       la1 = la2
    end do

    if ( mgt%bottom_solver == 3 .and. parallel_nprocs() > 0 ) then
       ! FIXME
       la2 = mgt%cc(1)%la
       call layout_build_derived(la1, la2)
       call multifab_build(mgt%rh1, la1, mgt%nc, 0)
       call multifab_build(mgt%uu1, la1, mgt%nc, 0)
       call multifab_build(mgt%ss1, la1, ns,     0)
    end if

    if ( present(dh) ) then
       mgt%dh(:,mgt%nlevels) = dh(:)
    else
       mgt%dh(:,mgt%nlevels) = 1.0_dp_t
    end if

    do i = mgt%nlevels-1, 1, -1
       mgt%dh(:,i) = mgt%dh(:,i+1)*2.0_dp_t
    end do

    if ( present(st_type) ) then
       allocate(mgt%st(mgt%nlevels))
       do i = mgt%nlevels, 1, -1
          call stencil_build(mgt%st(i), la1, mgt%dh(:,i), &
               type = st_type, nc = ns, nodal = nodal)
       end do
    end if

    !   Set the face_type array to be BC_DIR or BC_NEU depending on domain_bc
    mgt%face_type = BC_INT
    do id = 1,mgt%dim
       lo_dom = lwb(pd,id)
       hi_dom = upb(pd,id)
       do i = 1,mgt%nboxes
          lo_grid =  lwb(get_box(mgt%ss(mgt%nlevels), i),id)
          if (lo_grid .eq. lo_dom) mgt%face_type(i,id,1) = domain_bc(id,1)

          hi_grid = upb(get_box(mgt%ss(mgt%nlevels), i),id)
          if (hi_grid .eq. hi_dom) mgt%face_type(i,id,2) = domain_bc(id,2)
       end do
    end do

    if ( mgt%bottom_solver == 0 .and. .not. present(bottom_max_iter) ) mgt%bottom_max_iter = 20

    if ( mgt%cycle == MG_WCycle ) mgt%gamma = 2

    ! if only the bottom solver is 'solving' make sure that its eps is in effect
    if ( mgt%nlevels == 1 ) mgt%bottom_solver_eps = mgt%eps

    if ( .false. ) then
       if ( mgt%verbose > 0 .AND. parallel_IOProcessor() ) then
          write(unit=*, fmt=*) 'MG Solver'
          write(unit=*, fmt=*) 'nlevels           = ', mgt%nlevels
          write(unit=*, fmt=*) 'ng                = ', mgt%ng
          write(unit=*, fmt=*) 'nc                = ', mgt%nc
          write(unit=*, fmt=*) 'min_width         = ', mgt%min_width
          write(unit=*, fmt=*) 'max_nlevel        = ', mgt%max_nlevel
          write(unit=*, fmt=*) 'max_iter          = ', mgt%max_iter
          write(unit=*, fmt=*) 'eps               = ', mgt%eps
          write(unit=*, fmt=*) 'solver            = ', mgt%solver
          write(unit=*, fmt=*) 'smoother          = ', mgt%smoother
          write(unit=*, fmt=*) 'nu1               = ', mgt%nu1
          write(unit=*, fmt=*) 'nu2               = ', mgt%nu2
          write(unit=*, fmt=*) 'gamma             = ', mgt%gamma
          write(unit=*, fmt=*) 'omega             = ', mgt%omega
          write(unit=*, fmt=*) 'cycle             = ', mgt%cycle
          write(unit=*, fmt=*) 'bottom_solver     = ', mgt%bottom_solver
          write(unit=*, fmt=*) 'bottom_solver_eps = ', mgt%bottom_solver_eps
          write(unit=*, fmt=*) 'bottom_max_iter   = ', mgt%bottom_max_iter
       end if
    end if
  end subroutine mg_tower_build

  subroutine mg_tower_destroy(mgt)
    type(mg_tower), intent(inout) :: mgt
    integer :: i
    do i = 1, mgt%nlevels
       call destroy(mgt%cc(i))
       call destroy(mgt%ff(i))
       call destroy(mgt%dd(i))
       call destroy(mgt%ss(i))
       call destroy(mgt%mm(i))
       if ( i /= mgt%nlevels ) then
          call destroy(mgt%uu(i))
       end if
    end do
    if ( associated(mgt%st) ) then
       do i = 1, mgt%nlevels
          call stencil_destroy(mgt%st(i))
       end do
       deallocate(mgt%st)
    end if
    deallocate(mgt%cc, mgt%ff, mgt%dd, mgt%uu, mgt%mm, mgt%ss)
    deallocate(mgt%dh, mgt%pd)
    deallocate(mgt%tm)
    deallocate(mgt%face_type)
    if ( built_q(mgt%sparse_object) ) call destroy(mgt%sparse_object)
    if ( built_q(mgt%rh1)           ) call destroy(mgt%rh1)
    if ( built_q(mgt%uu1)           ) call destroy(mgt%uu1)
    if ( built_q(mgt%ss1)           ) call destroy(mgt%ss1)
  end subroutine mg_tower_destroy

  function max_mg_levels(ba, min_size) result(r)
    type(boxarray), intent(in) :: ba
    integer, intent(in), optional :: min_size
    integer :: r
    integer, parameter :: rrr = 2
    type(box) :: bx, bx1
    integer :: i, rr, lmn, dm
    lmn = 1; if ( present(min_size) ) lmn = min_size
    r = 1
    rr = rrr
    dm = ba%dim
    do
       do i = 1, size(ba%bxs)
          bx = ba%bxs(i)
          bx1 = coarsen(bx, rr)
          if ( any(extent(bx1) < lmn) ) return
          if ( bx /= refine(bx1, rr) ) then
             return
          end if
       end do
       rr = rr*rrr
       r  = r + 1
    end do
  end function max_mg_levels

  subroutine mg_tower_v_cycle(mgt, uu, rh)
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in)    :: rh

    call mg_tower_cycle(mgt, MG_VCYCLE, mgt%nlevels, mgt%ss(mgt%nlevels), &
         uu, rh, &
         mgt%mm(mgt%nlevels), mgt%nu1, mgt%nu2, mgt%gamma)
  end subroutine mg_tower_v_cycle

  subroutine mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)
    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in) :: rh
    type( multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer stat
    integer i
    call bl_assert(lev==1, "bottom_solve only at level 1")

    stat = 0
    select case ( mgt%bottom_solver )
    case (0)
       do i = 1, mgt%bottom_max_iter
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (1)
       call bicgstab_solve(ss, uu, rh, mm, &
                           mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%verbose, stat = stat)
    case (2)
       call cg_solve(ss, uu, rh, mm, &
                     mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%verbose, stat = stat)
       do i = 1, 2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (3)

       if ( parallel_nprocs() > 0 ) then
          ! call bl_error("can't do this yet")
          ! Copy into a local mf the rhs, copy out the solution to the distributed uu
          call copy(mgt%rh1, rh)
          call sparse_solve(mgt%sparse_object, mgt%uu1, mgt%rh1, &
               mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%verbose, stat)
          call copy(uu, mgt%uu1)
       else
          call sparse_solve(mgt%sparse_object, uu, rh, &
               mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%verbose, stat)
       end if

    case default
       call bl_error("MG_TOWER_BOTTOM_SOLVE: no such solver: ", mgt%bottom_solver)
    end select
    if ( stat /= 0 ) then
       call bl_warn("BREAKDOWN in bottom_solver: trying smoother")
       do i = 1, 2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    end if

  end subroutine mg_tower_bottom_solve

  subroutine mg_defect(ss, dd, ff, uu, mm)

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm

    call stencil_apply(ss, dd, uu, mm)

    call saxpy(dd, ff, -1.0_dp_t, dd)

  end subroutine mg_defect

  subroutine grid_res(mgt, lev, ss, dd, ff, uu, mm, face_type)
    use mg_defect_module
    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    type(mg_tower), intent(inout) :: mgt
    integer, intent(in) :: lev
    integer, intent(in) :: face_type(:,:,:)
    integer :: i, n
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: nodal_ng

    nodal_ng = 0; if ( nodal_q(uu) ) nodal_ng = 1

    call multifab_fill_boundary(uu)

    do i = 1, mgt%nboxes
       if ( multifab_remote(dd, i) ) cycle
       dp => dataptr(dd, i)
       fp => dataptr(ff, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, mgt%nc
          select case(mgt%dim)
          case (1)
             call grid_laplace(sp(:,1,1,:), dp(:,1,1,n), fp(:,1,1,n), up(:,1,1,n), &
                               mp(:,1,1,1), mgt%ng, nodal_ng)
          case (2)
             call grid_laplace(sp(:,:,1,:), dp(:,:,1,n), fp(:,:,1,n), up(:,:,1,n), &
                               mp(:,:,1,1), mgt%ng, nodal_ng, face_type(i,:,:))
          case (3)
             call grid_laplace(sp(:,:,:,:), dp(:,:,:,n), fp(:,:,:,n), up(:,:,:,n), &
                               mp(:,:,:,1), mgt%ng, nodal_ng, face_type(i,:,:))
          end select
       end do
    end do
  end subroutine grid_res

  subroutine mg_tower_restriction(mgt, lev, crse, fine, mm_fine, mm_crse)
    use mg_restriction_module
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    type(imultifab), intent(in)   :: mm_fine
    type(imultifab), intent(in)   :: mm_crse
    type(mg_tower), intent(inout) :: mgt
    integer, intent(in) :: lev
    integer :: loc(mgt%dim), lof(mgt%dim)
    integer :: lom_fine(mgt%dim)
    integer :: lom_crse(mgt%dim)
    integer :: lo(mgt%dim), hi(mgt%dim)
    integer :: i, n, ir(mgt%dim)
    logical :: nodal_flag
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp_fine(:,:,:,:)
    integer        , pointer :: mp_crse(:,:,:,:)

    integer :: mg_restriction_mode

    ir = 2

    nodal_flag = nodal_q(crse)

    if ( nodal_flag ) then
       call multifab_fill_boundary(fine)
       call setval(crse, 0.0_dp_t, all=.true.)
       mg_restriction_mode = 1
    end if 

    do i = 1, mgt%nboxes
       if ( multifab_remote(crse, i) ) cycle
       cp => dataptr(crse, i)
       fp => dataptr(fine, i)
       mp_fine => dataptr(mm_fine, i)
       mp_crse => dataptr(mm_crse, i)

       loc(:) = lwb(get_ibox(crse,i)) - crse%ng
       lof(:) = lwb(get_ibox(fine,i)) - fine%ng
       lom_fine(:) = lwb(get_ibox(mm_fine,i)) - mm_fine%ng
       lom_crse(:) = lwb(get_ibox(mm_crse,i)) - mm_crse%ng
       lo(:)  = lwb(get_ibox(crse,i))
       hi(:)  = upb(get_ibox(crse,i))

       do n = 1, mgt%nc
          select case ( mgt%dim)
          case (1)
             if ( .not.nodal_flag ) then
               call cc_restriction(cp(:,1,1,n), loc, fp(:,1,1,n), lof, lo, hi, ir)
             else
               call nodal_restriction(cp(:,1,1,n), loc, fp(:,1,1,n), lof, &
                    mp_fine(:,1,1,1), lom_fine, &
                    mp_crse(:,1,1,1), lom_crse, lo, hi, ir, .false., mg_restriction_mode)
             end if  
          case (2)
             if ( .not.nodal_flag ) then
               call cc_restriction(cp(:,:,1,n), loc, fp(:,:,1,n), lof, lo, hi, ir)
             else
               call nodal_restriction(cp(:,:,1,n), loc, fp(:,:,1,n), lof, &
                    mp_fine(:,:,1,1), lom_fine, &
                    mp_crse(:,:,1,1), lom_crse, lo, hi, ir, .false., mg_restriction_mode)
             end if
          case (3)
             if ( .not.nodal_flag ) then
               call cc_restriction(cp(:,:,:,n), loc, fp(:,:,:,n), lof, lo, hi, ir)
             else
               call nodal_restriction(cp(:,:,:,n), loc, fp(:,:,:,n), lof, &
                    mp_fine(:,:,:,1), lom_fine, &
                    mp_crse(:,:,:,1), lom_crse, lo, hi, ir, .false., mg_restriction_mode)
             end if
          end select
       end do
    end do
  end subroutine mg_tower_restriction

  subroutine mg_tower_smoother(mgt, lev, ss, uu, ff, mm)
    use mg_smoother_module
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ff
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, n, nn
    integer :: lo(mgt%dim)
    logical :: skewed

    if ( cell_centered_q(uu) ) then
       select case ( mgt%smoother )
       case ( MG_SMOOTHER_GS_RB )
          do nn = 0, 1
             call multifab_fill_boundary(uu)
             do i = 1, mgt%nboxes
                if ( multifab_remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                ! skewed = skewed_q(mp)
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (1)
                      call gs_rb_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                           mp(:,1,1,1), lo, mgt%ng, nn)
                   case (2)
                      call gs_rb_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                           mp(:,:,1,1), lo, mgt%ng, nn)
                   case (3)
                      call gs_rb_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                           mp(:,:,:,1), lo, mgt%ng, nn)
                   end select
                end do
             end do
          end do
       case ( MG_SMOOTHER_JACOBI )
          call multifab_fill_boundary(uu)
          do i = 1, mgt%nboxes
             if ( multifab_remote(ff, i) ) cycle
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case ( mgt%dim)
                case (1)
                   call jac_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                        mp(:,1,1,1), mgt%ng)
                case (2)
                   call jac_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                        mp(:,:,1,1), mgt%ng)
                case (3)
                   call jac_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                        mp(:,:,:,1), mgt%ng)
                end select
             end do
          end do
       case ( MG_SMOOTHER_GS_LEX )
          call multifab_fill_boundary(uu)
          do i = 1, mgt%nboxes
             if ( multifab_remote(ff, i) ) cycle
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case ( mgt%dim)
                case (1)
                   call gs_lex_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                        mp(:,1,1,1), mgt%ng)
                case (2)
                   call gs_lex_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                        mp(:,:,1,1), mgt%ng)
                case (3)
                   call gs_lex_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                        mp(:,:,:,1), mgt%ng)
                end select
             end do
          end do
       case default
          call bl_error("MG_TOWER_SMOOTHER: no such smoother")
       end select
    else 
       call multifab_fill_boundary(uu)
       do i = 1, mgt%nboxes
          if ( multifab_remote(ff, i) ) cycle
          up => dataptr(uu, i)
          fp => dataptr(ff, i)
          sp => dataptr(ss, i)
          mp => dataptr(mm, i)
          lo =  lwb(get_box(ss, i))
          do n = 1, mgt%nc
             select case ( mgt%dim)
             case (1)
                call nodal_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                     mp(:,1,1,1), lo, mgt%ng)
             case (2)
                call nodal_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                     mp(:,:,1,1), lo, mgt%ng)
             case (3)
                call nodal_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                     mp(:,:,:,1), lo, mgt%ng)
             end select
          end do
       end do

       call multifab_internal_sync(uu)
    end if

  end subroutine mg_tower_smoother

  subroutine mg_tower_prolongation(mgt, lev, uu, uu1)
    use mg_prolongation_module
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu, uu1
    integer, intent(in) :: lev
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    type(box) :: nbox, nbox1
    integer :: i, n, ir(mgt%dim)
    logical :: nodal_flag
    ir = 2

    nodal_flag = nodal_q(uu)

    if ( .not.nodal_flag ) then
      do i = 1, mgt%nboxes
         if ( multifab_remote(mgt%ff(lev), i) ) cycle
         fp => dataptr(uu,  i, get_box(uu,i))
         cp => dataptr(uu1, i, get_box(uu1,i))
         do n = 1, mgt%nc
            select case ( mgt%dim)
            case (1)
               call pc_c_prolongation(fp(:,1,1,n), cp(:,1,1,n), ir)
            case (2)
               call pc_c_prolongation(fp(:,:,1,n), cp(:,:,1,n), ir)
            case (3)
               call pc_c_prolongation(fp(:,:,:,n), cp(:,:,:,n), ir)
            end select
         end do
      end do
    else
      do i = 1, mgt%nboxes
         if ( multifab_remote(mgt%ff(lev), i) ) cycle
         nbox  = box_grow_n_f(get_box(uu,i),1,1)
         nbox1 = box_grow_n_f(get_box(uu1,i),1,1)
         fp => dataptr(uu,  i, nbox )
         cp => dataptr(uu1, i, nbox1)
         do n = 1, mgt%nc
            select case ( mgt%dim)
            case (1)
               call nodal_prolongation(fp(:,1,1,n), cp(:,1,1,n), ir)
            case (2)
               call nodal_prolongation(fp(:,:,1,n), cp(:,:,1,n), ir)
            case (3)
               call nodal_prolongation(fp(:,:,:,n), cp(:,:,:,n), ir)
            end select
         end do
      end do
    endif

  end subroutine mg_tower_prolongation

  function mg_tower_converged(mgt, lev, dd, uu, Anorm, Ynorm) result(r)
    logical :: r
    type(mg_tower), intent(inout) :: mgt
    integer, intent(in) :: lev
    real(dp_t), intent(in) :: Anorm, Ynorm
    type(multifab), intent(in) :: dd
    type(multifab), intent(in) :: uu
    r = itsol_converged(dd, uu, Anorm, Ynorm, mgt%eps)
  end function mg_tower_converged

  recursive subroutine mg_tower_cycle(mgt, cyc, lev, ss, uu, rh, mm, nu1, nu2, gamma)
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(in) :: rh
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer, intent(in) :: nu1, nu2
    integer, intent(inout) :: gamma
    integer, intent(in) :: cyc
    integer i
    logical do_diag
    real(dp_t) :: nrm

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    call timer_start(mgt%tm(lev))
    if ( lev == 1 ) then
       if (do_diag) then
         call mg_defect(ss, mgt%cc(lev), rh, uu, mm)
         nrm = norm_inf(mgt%cc(lev))
         if ( parallel_IOProcessor() ) then
            print *,'DN: NORM BEFORE BOTTOM ',lev, nrm
         end if
       end if
       call mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)
       if ( cyc == MG_FCycle ) gamma = 1
       if (do_diag) then
         call mg_defect(ss, mgt%cc(lev), rh, uu, mm)
         nrm = norm_inf(mgt%cc(lev))
         if ( parallel_IOProcessor() ) then
            print *,'DN: NORM AFTER BOTTOM ',lev, nrm
         end if
       end if
    else 

       if (do_diag) then
         nrm = norm_inf(rh)
         if ( parallel_IOProcessor() ) then
            print *,'DN: NORM BEFORE RELAX ',lev, nrm
         end if
       end if
       do i = 1, nu1
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
       call mg_defect(ss, mgt%cc(lev), rh, uu, mm)

       if (do_diag) then
         nrm = norm_inf(mgt%cc(lev))
         if ( parallel_IOProcessor() ) then
            print *,'DN: NORM AFTER RELAX ',lev, nrm
         end if
       end if

       call mg_tower_restriction(mgt, lev, mgt%dd(lev-1), mgt%cc(lev), &
                                 mgt%mm(lev),mgt%mm(lev-1))
       call setval(mgt%uu(lev-1), zero, all = .TRUE.)
       do i = gamma, 1, -1
          call mg_tower_cycle(mgt, cyc, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), &
                              mgt%dd(lev-1), mgt%mm(lev-1), nu1, nu2, gamma)
       end do
       ! uu  += cc, done, by convention, using the prolongation routine.
       call mg_tower_prolongation(mgt, lev, uu, mgt%uu(lev-1))

       if (do_diag) then
         call mg_defect(ss, mgt%cc(lev), rh, uu, mm)
         nrm = norm_inf(mgt%cc(lev))
         if ( parallel_IOProcessor() ) then
            print *,'UP: NORM AFTER INTERP ',lev, nrm
         end if
       end if

       do i = 1, nu2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

       if (do_diag) then
         call mg_defect(ss, mgt%cc(lev), rh, uu, mm)
         nrm = norm_inf(mgt%cc(lev))
         if ( parallel_IOProcessor() ) then
            print *,'UP: NORM AFTER RELAX ',lev, nrm
         end if
       end if

       ! if at top of tower and doing an FCycle reset gamma
       if ( lev == mgt%nlevels .AND. cyc == MG_FCycle ) gamma = 2
    end if

    call timer_stop(mgt%tm(lev))

  end subroutine mg_tower_cycle

  subroutine mini_cycle(mgt, cyc, lev, ss, uu, rh, mm, nu1, nu2, gamma)
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(in) :: rh
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer, intent(in) :: nu1, nu2
    integer, intent(inout) :: gamma
    integer, intent(in) :: cyc
    integer i

    call timer_start(mgt%tm(lev))

    if ( lev == 1 ) then

       call mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)

    else 

       call mg_defect(ss, mgt%cc(lev), rh, uu, mm)

       call mg_tower_restriction(mgt, lev, mgt%dd(lev-1), mgt%cc(lev), &
                                 mgt%mm(lev),mgt%mm(lev-1))

       call setval(mgt%uu(lev-1), zero, all = .TRUE.)

       call mg_tower_bottom_solve(mgt, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), &
                                  mgt%dd(lev-1), mgt%mm(lev-1))

       call mg_tower_prolongation(mgt, lev, uu, mgt%uu(lev-1))

       do i = 1, nu2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

    end if

    call timer_stop(mgt%tm(lev))

  end subroutine mini_cycle

  subroutine mg_tower_solve(mgt, uu, rh, qq, num_iter, defect_history, defect_dirname, stat)
    use fabio_module
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu, rh
    integer, intent(out), optional :: stat
    real(dp_t), intent(out), optional :: qq(0:)
    integer, intent(out), optional :: num_iter
    character(len=*), intent(in), optional :: defect_dirname
    logical, intent(in), optional :: defect_history
    integer :: gamma
    integer :: it, cyc
    real(dp_t) :: ynorm, Anorm, nrm(3)
    integer :: n_qq, i_qq
    logical :: ldef
    character(len=128) :: defbase

    ldef = .false.; if ( present(defect_history) ) ldef = defect_history
    if ( ldef ) then
       if ( .not. present(defect_dirname) ) then
          call bl_error("MG_TOWER_SOLVE: defect_history but no defect_dirname")
       end if
    end if
       
    if ( present(stat) ) stat = 0
    if ( mgt%cycle == MG_FCycle ) then
       gamma = 2
    else
       gamma = mgt%gamma
    end if

    n_qq = 0
    i_qq = 0
    if ( present(qq) ) then
       n_qq = size(qq)
    end if
    cyc   = mgt%cycle
    Anorm = stencil_norm(mgt%ss(mgt%nlevels))
    ynorm = norm_inf(rh)
    call mg_defect(mgt%ss(mgt%nlevels), &
                   mgt%dd(mgt%nlevels), rh, uu, mgt%mm(mgt%nlevels))
    if ( i_qq < n_qq ) then
       qq(i_qq) = norm_l2(mgt%dd(mgt%nlevels))
       i_qq = i_qq + 1
    end if
    if ( ldef ) then
       write(unit = defbase, fmt='("def",I3.3)') 0
       call fabio_write(mgt%dd(mgt%nlevels), defect_dirname, defbase)
    end if
       
    if ( mgt%verbose > 0 ) then
       nrm(3) = norm_inf(mgt%dd(mgt%nlevels))
       nrm(1) = norm_inf(uu)
       nrm(2) = norm_inf(rh)
       if ( parallel_IOProcessor() ) then
          write(unit=*, &
               fmt='(i3,": Unorm= ",g15.8,", Rnorm= ",g15.8,", Ninf(defect) = ",g15.8, ", Anorm=",g15.8)') &
               0, nrm, Anorm
       end if
    end if
    if ( mg_tower_converged(mgt, mgt%nlevels, mgt%dd(mgt%nlevels), uu, Anorm, Ynorm) ) then
       if ( present(stat) ) stat = 0
       if ( mgt%verbose > 0 .AND. parallel_IOProcessor() ) then
          write(unit=*, fmt='("MG finished at on input")') 
       end  if
       return
    end if
    do it = 1, mgt%max_iter
       call mg_tower_cycle(mgt, cyc, mgt%nlevels, mgt%ss(mgt%nlevels), &
                           uu, rh, mgt%mm(mgt%nlevels), mgt%nu1, mgt%nu2, &
                           gamma)
       call mg_defect(mgt%ss(mgt%nlevels), &
                            mgt%dd(mgt%nlevels), rh, uu, mgt%mm(mgt%nlevels))
       if ( mgt%verbose > 0 ) then
          nrm(1) = norm_inf(mgt%dd(mgt%nlevels))
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Ninf(defect) = ",g15.8)') it, nrm(1)
          end if
       end if
       if ( i_qq < n_qq ) then
          qq(i_qq) = norm_l2(mgt%dd(mgt%nlevels))
          i_qq = i_qq + 1
       end if
       if ( ldef ) then
          write(unit = defbase,fmt='("def",I3.3)') it
          call fabio_write(mgt%dd(mgt%nlevels), defect_dirname, defbase)
       end if
       if ( mg_tower_converged(mgt, mgt%nlevels, mgt%dd(mgt%nlevels), uu, Anorm, Ynorm) ) exit
    end do
    if ( mgt%verbose > 0 .AND. parallel_IOProcessor() ) then
       write(unit=*, fmt='("MG finished at ", i3, " iterations")') it
    end  if

    if ( it > mgt%max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("MG_TOWER_SOLVE: failed to converge in max_iter iterations")
       end if
    end if

    if ( present(num_iter) ) num_iter = it

  end subroutine mg_tower_solve

end module mg_module

module new_interp_module

  use bl_types

  implicit none

  real(kind=dp_t), parameter :: ONE = 1.0_dp_t
  real(kind=dp_t), parameter :: HALF = 0.5_dp_t
  real(kind=dp_t), parameter :: TWO  = 2.0_dp_t
  real(kind=dp_t), parameter :: ZERO = 0.0_dp_t

contains

  subroutine bilinear_interp_expensive_1d(fot, fin, bc)
    real(kind=dp_t),intent(in)  :: fin(0:)
    real(kind=dp_t),intent(out) :: fot(0:)
    integer, intent(in) :: bc(2,1)

    real(kind=dp_t) dhi(1), dho(1)
    integer i
    real(kind=dp_t) wil, wih
    real(kind=dp_t) gamma, dii, doi
    integer nxi, nxo
    integer ii
    nxi = size(fin,dim=1)-2
    nxo = size(fot,dim=1)-2
    dhi(1) = ONE/nxi
    dho(1) = ONE/nxo
    do i = 1, nxo
       doi  = (i-HALF)*dho(1)
       ii   = int(doi/dhi(1) + HALF)
       dii  = (ii-HALF)*dhi(1)
       gamma = (doi-dii)/dhi(1)
       wil = ONE-gamma
       wih = gamma
       if ( ii < 1 ) then
          if ( bc(1,1) == 1 ) then
             wil = ZERO
             wih = doi/(dhi(1)/TWO)
          else if ( bc(1,1) == 2 ) then
          else if ( bc(1,1) == 3 ) then
             stop 'no neumann yet:x:1'
          else
             stop 'bad bc:x:2'
          end if
       end if
       if ( ii > nxi - 1 ) then
          if ( bc(2,1) == 1 ) then
             wil = ONE-(doi-dii)/(dhi(1)/TWO)
             wih = ZERO
          else if ( bc(2,1) == 2 ) then
          else if ( bc(2,1) == 3 ) then
             stop 'no neumann yet:x:2'
          else
             stop 'bad bc:x:2'
          end if
       end if
       fot(i) = &
            + wil*fin(ii  ) &
            + wih*fin(ii+1) 
    end do
  end subroutine bilinear_interp_expensive_1d

  subroutine bilinear_interp_expensive_2d(fot, fin, bc)
    real(kind=dp_t),intent(in)  :: fin(0:,0:)
    real(kind=dp_t),intent(out) :: fot(0:,0:)
    integer, intent(in) :: bc(2,2)

    real(kind=dp_t) dhi(2), dho(2)
    integer i, j
    real(kind=dp_t) wil, wih, wjl, wjh
    real(kind=dp_t) gamma, dii, doi, dij, doj
    integer nxi, nyi, nxo, nyo
    integer ii, jj
    nxi = size(fin,dim=1)-2
    nyi = size(fin,dim=2)-2
    nxo = size(fot,dim=1)-2
    nyo = size(fot,dim=2)-2
    dhi(1) = ONE/nxi
    dhi(2) = ONE/nyi
    dho(1) = ONE/nxo
    dho(2) = ONE/nyo
    do j = 1, nyo
       doj = (j-HALF)*dho(2)
       jj  = int(doj/dhi(2) + HALF)
       dij = (jj-HALF)*dhi(2)
       gamma = (doj-dij)/dhi(2)
       wjl = ONE-gamma
       wjh = gamma
       if ( jj < 1 ) then
          if ( bc(1,2) == 1 ) then
             wjl = ZERO
             wjh = doj/(dhi(2)/TWO)
          else if ( bc(1,2) == 2 ) then
             !     Keep things the same
          else if ( bc(1,2) == 3 ) then
             stop 'no neumann yet:y:1'
          else
             stop 'bad bc:y:1'
          end if
       endif
       if ( jj > nyi - 1) then
          if ( bc(1,2) == 1 ) then
             wjl = ONE-(doj-dij)/(dhi(2)/TWO)
             wjh = ZERO
          else if ( bc(2,2) == 2 ) then
          else if ( bc(2,2) == 3 ) then
             stop 'no neumann yet:y:2'
          else
             stop 'bad bc:y:2'
          endif
       end if
       do i = 1, nxo
          doi  = (i-HALF)*dho(1)
          ii   = int(doi/dhi(1) + HALF)
          dii  = (ii-HALF)*dhi(1)
          gamma = (doi-dii)/dhi(1)
          wil = ONE-gamma
          wih = gamma
          if ( ii < 1 ) then
             if ( bc(1,1) == 1 ) then
                wil = ZERO
                wih = doi/(dhi(1)/TWO)
             else if ( bc(1,1) == 2 ) then
             else if ( bc(1,1) == 3 ) then
                stop 'no neumann yet:x:1'
             else
                stop 'bad bc:x:2'
             end if
          end if
          if ( ii > nxi - 1 ) then
             if ( bc(2,1) == 1 ) then
                wil = ONE-(doi-dii)/(dhi(1)/TWO)
                wih = ZERO
             else if ( bc(2,1) == 2 ) then
             else if ( bc(2,1) == 3 ) then
                stop 'no neumann yet:x:2'
             else
                stop 'bad bc:x:2'
             end if
          end if
          fot(i,j) = &
               + wil*wjl*fin(ii  ,  jj) &
               + wil*wjh*fin(ii  ,jj+1) &
               + wih*wjl*fin(ii+1,  jj) &
               + wih*wjh*fin(ii+1,jj+1) 
       end do
    end do
  end subroutine bilinear_interp_expensive_2d

  subroutine trilinear_interp_expensive_3d(fot, fin, bc)
    integer, intent(in) ::  bc(2,3)
    real(kind=dp_t) dhi(3), dho(3)
    real(kind=dp_t), intent(in) ::  fin(0:,0:,0:)
    real(kind=dp_t), intent(out) ::  fot(0:,0:,0:)
    integer i, j, k
    real(kind=dp_t) wil, wih, wjl, wjh, wkl, wkh
    real(kind=dp_t) gamma, dii, doi, dij, doj, dok, dik
    integer ii, jj, kk
    integer nxi, nyi, nxo, nyo, nzo, nzi

    nxi = size(fin,dim=1)-2
    nyi = size(fin,dim=2)-2
    nzi = size(fin,dim=3)-2
    nxo = size(fot,dim=1)-2
    nyo = size(fot,dim=2)-2
    nzo = size(fot,dim=3)-2
    dhi(1) = ONE/nxi
    dhi(2) = ONE/nyi
    dhi(3) = ONE/nzi
    dho(1) = ONE/nxo
    dho(2) = ONE/nyo
    dho(3) = ONE/nzo
    do k = 1, nzo
       !     The Z direction
       dok = (k-HALF)*dho(3)
       kk  = int(dok/dhi(3) + HALF)
       dik = (kk-HALF)*dhi(3)
       gamma = (dok-dik)/dhi(3)
       wkl = ONE-gamma
       wkh = gamma
       if ( kk < 1 ) then
          wkl = ZERO
          wkh = dok/(dhi(3)/TWO)
       end if
       if ( kk + 1 > nzi ) then
          wkl = ONE-(dok-dik)/(dhi(3)/TWO)
          wkh = ZERO
       end if
       do j = 1, nyo
          ! The Y direction
          doj = (j-HALF)*dho(2)
          jj  = int(doj/dhi(2) + HALF)
          dij = (jj-HALF)*dhi(2)
          gamma = (doj-dij)/dhi(2)
          wjl = ONE-gamma
          wjh = gamma
          if ( jj < 1 ) then
             wjl = ZERO
             wjh = doj/(dhi(2)/TWO)
          end if
          if ( jj + 1 > nyi ) then
             wjl = ONE-(doj-dij)/(dhi(2)/TWO)
             wjh = ZERO
          end if
          do i = 1, nxo
             ! The X direction
             doi  = (i-HALF)*dho(1)
             ii   = int(doi/dhi(1) + HALF)
             dii  = (ii-HALF)*dhi(1)
             gamma = (doi-dii)/dhi(1)
             wil = ONE-gamma
             wih = gamma
             if ( ii < 1 ) then
                wil = ZERO
                wih = doi/(dhi(1)/TWO)
             end if
             if ( ii + 1 > nxi ) then
                wil = ONE-(doi-dii)/(dhi(1)/TWO)
                wih = ZERO
             end if
             fot(i,j,k) = &
                  + wil*wjl*wkl*fin(ii  ,  jj,  kk) &
                  + wil*wjl*wkh*fin(ii  ,  jj,kk+1) &
                  + wil*wjh*wkl*fin(ii  ,jj+1,  kk) &
                  + wil*wjh*wkh*fin(ii  ,jj+1,kk+1) &
                  + wih*wjl*wkl*fin(ii+1,  jj,  kk) &
                  + wih*wjl*wkh*fin(ii+1,  jj,kk+1) &
                  + wih*wjh*wkl*fin(ii+1,jj+1,  kk) &
                  + wih*wjh*wkh*fin(ii+1,jj+1,kk+1)
          end do
       end do
    end do
  end subroutine trilinear_interp_expensive_3d

end module new_interp_module
