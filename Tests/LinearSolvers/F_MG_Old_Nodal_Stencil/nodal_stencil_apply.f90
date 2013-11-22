module nodal_stencil_apply_module

  use bl_types
  use multifab_module
  use bc_functions_module
  use impose_neumann_bcs_module
  use stencil_types_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t
  real (kind = dp_t), private, parameter :: THREE = 3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR  = 4.0_dp_t
  real (kind = dp_t), private, parameter :: FIVE  = 5.0_dp_t
  real (kind = dp_t), private, parameter :: SIX   = 6.0_dp_t
  real (kind = dp_t), private, parameter :: SEVEN = 7.0_dp_t
  real (kind = dp_t), private, parameter :: EIGHT = 8.0_dp_t
  real (kind = dp_t), private, parameter :: TEN   = 10.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t
  real (kind = dp_t), private, parameter :: FOURTH= 0.25_dp_t
  real (kind = dp_t), private, parameter :: THIRD = 1.0_dp_t/3.0_dp_t
  real (kind = dp_t), private, parameter :: SIXTH = 1.0_dp_t/6.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR_THIRD = 4.0_dp_t/3.0_dp_t

contains

  subroutine stencil_apply_1d_nodal(ss, dd, uu, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in   ) :: ss(0:,:)
    real (kind = dp_t), intent(inout) :: dd(0:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:)
    integer           , intent(in   ) :: mm(:)

    integer :: i,lo(1)
 
    dd = ZERO

    lo = 1
    call impose_neumann_bcs_1d(uu,mm,lo,ng)
   
    i = 1
    if (.not. bc_dirichlet(mm(i),1,0)) then
      dd(i) = ss(0,i)*uu(i) + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1)
    end if

    do i = 2,size(ss,dim=2)-1
      dd(i) = ss(0,i)*uu(i) + ss(1,i) * uu(i+1) + ss(2,i) * uu(i-1)
    end do

    i = size(ss,dim=2)
    if (.not. bc_dirichlet(mm(i),1,0)) then
      dd(i) = ss(0,i)*uu(i) + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1)
    end if

  end subroutine stencil_apply_1d_nodal

  subroutine stencil_apply_2d_nodal(ss, dd, uu, mm, ng, stencil_type)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(0:,:,:)
    integer           , intent(in   ) :: mm(:,:)
    integer           , intent(in   ) :: stencil_type

    integer :: i,j,lo(2),nx,ny
    logical :: zeroit,iface,jface

    lo = 1
    call impose_neumann_bcs_2d(uu,mm,lo,ng)
 
    nx = size(ss,dim=2)
    ny = size(ss,dim=3)

    if (stencil_type .eq. ND_CROSS_STENCIL) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                dd(i,j) = ss(0,i,j)*uu(i,j) + ss(1,i,j) * uu(i+1,j  ) &
                     + ss(2,i,j) * uu(i-1,j  ) &
                     + ss(3,i,j) * uu(i  ,j+1) &
                     + ss(4,i,j) * uu(i  ,j-1) 
             end if
          end do
       end do

    else if (stencil_type .eq. ND_DENSE_STENCIL) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                dd(i,j) = ss(0,i,j)*uu(i,j) + ss(1,i,j) * uu(i-1,j-1) &
                     + ss(2,i,j) * uu(i  ,j-1) &
                     + ss(3,i,j) * uu(i+1,j-1) &
                     + ss(4,i,j) * uu(i-1,j  ) &
                     + ss(5,i,j) * uu(i+1,j  ) &
                     + ss(6,i,j) * uu(i-1,j+1) &
                     + ss(7,i,j) * uu(i  ,j+1) &
                     + ss(8,i,j) * uu(i+1,j+1)
             end if
          end do
       end do

    else
       call bl_error("stencil_apply_2d_nodal: dont know this stencil_type")
    end if

  end subroutine stencil_apply_2d_nodal

  subroutine stencil_apply_3d_nodal(ss, dd, uu, mm, ng, stencil_type, uniform_dh, bottom_solver)
    integer           , intent(in   ) :: ng
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(inout) :: dd(0:,0:,0:)
    real (kind = dp_t), intent(in   ) :: ss(0:,:,:,:)
    integer           , intent(in   ) :: mm(:,:,:)
    integer           , intent(in   ) :: stencil_type
    logical           , intent(in   ) :: uniform_dh, bottom_solver

    integer :: i,j,k,lo(3),nx,ny,nz
    logical :: jface,kface

    lo = 1

    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    nz = size(ss,dim=4)
    ny = size(ss,dim=3)
    nx = size(ss,dim=2)

    if (stencil_type .eq. ND_CROSS_STENCIL) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,jface,kface) IF(.not.bottom_solver)
       do k = 1,nz
          kface = .false. ; if ( (k.eq.1) .or. (k.eq.nz) ) kface = .true.

          do j = 1,ny
             jface = .false. ; if ( (j.eq.1) .or. (j.eq.ny) ) jface = .true.

             do i = 1,nx

                if ( (jface .or. kface .or. (i.eq.1) .or. (i.eq.nx)) .and. bc_dirichlet(mm(i,j,k),1,0) ) then
                   dd(i,j,k) = ZERO
                else
                   dd(i,j,k) = &
                        ss(0,i,j,k) * uu(i,j,k)       + &
                        ss(1,i,j,k) * uu(i+1,j  ,k  ) + &
                        ss(2,i,j,k) * uu(i-1,j  ,k  ) + &
                        ss(3,i,j,k) * uu(i  ,j+1,k  ) + &
                        ss(4,i,j,k) * uu(i  ,j-1,k  ) + &
                        ss(5,i,j,k) * uu(i  ,j  ,k+1) + &
                        ss(6,i,j,k) * uu(i  ,j  ,k-1)
                end if

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if (stencil_type .eq. ND_DENSE_STENCIL) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,jface,kface) IF(.not.bottom_solver)
       do k = 1,nz
          kface = .false. ; if ( (k.eq.1) .or. (k.eq.nz) ) kface = .true.

          do j = 1,ny
             jface = .false. ; if ( (j.eq.1) .or. (j.eq.ny) ) jface = .true.

             do i = 1,nx

                if ( (jface .or. kface .or. (i.eq.1) .or. (i.eq.nx)) .and. bc_dirichlet(mm(i,j,k),1,0) ) then
                   dd(i,j,k) = ZERO
                else
                   dd(i,j,k) = ss(0,i,j,k)*uu(i,j,k) &
                        + ss( 1,i,j,k) * uu(i-1,j-1,k-1) + ss( 2,i,j,k) * uu(i  ,j-1,k-1) &
                        + ss( 3,i,j,k) * uu(i+1,j-1,k-1) + ss( 4,i,j,k) * uu(i-1,j  ,k-1) &
                        + ss( 5,i,j,k) * uu(i+1,j  ,k-1) + ss( 6,i,j,k) * uu(i-1,j+1,k-1) &
                        + ss( 7,i,j,k) * uu(i  ,j+1,k-1) + ss( 8,i,j,k) * uu(i+1,j+1,k-1) &
                        + ss( 9,i,j,k) * uu(i-1,j-1,k  ) + ss(10,i,j,k) * uu(i+1,j-1,k  ) &
                        + ss(11,i,j,k) * uu(i-1,j+1,k  ) + ss(12,i,j,k) * uu(i+1,j+1,k  ) &
                        + ss(13,i,j,k) * uu(i-1,j-1,k+1) + ss(14,i,j,k) * uu(i  ,j-1,k+1) &
                        + ss(15,i,j,k) * uu(i+1,j-1,k+1) + ss(16,i,j,k) * uu(i-1,j  ,k+1) &
                        + ss(17,i,j,k) * uu(i+1,j  ,k+1) + ss(18,i,j,k) * uu(i-1,j+1,k+1) &
                        + ss(19,i,j,k) * uu(i  ,j+1,k+1) + ss(20,i,j,k) * uu(i+1,j+1,k+1)

                   if ((size(ss,dim=1) .eq. 27) .and. (.not. uniform_dh)) then
                      !
                      ! Add faces (only non-zero for non-uniform dx)
                      !
                      dd(i,j,k) = dd(i,j,k) + &
                           ss(21,i,j,k) * uu(i-1,j  ,k  ) + ss(22,i,j,k) * uu(i+1,j  ,k  ) &
                           + ss(23,i,j,k) * uu(i  ,j-1,k  ) + ss(24,i,j,k) * uu(i  ,j+1,k  ) &
                           + ss(25,i,j,k) * uu(i  ,j  ,k-1) + ss(26,i,j,k) * uu(i  ,j  ,k+1)
                   end if

                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       call bl_error("stencil_apply_3d_nodal: dont know this stencil_type")
    end if

  end subroutine stencil_apply_3d_nodal

end module nodal_stencil_apply_module
