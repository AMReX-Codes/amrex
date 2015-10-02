module stencil_util_module

  use bl_error_module
  use bl_types
  use bl_constants_module
  use multifab_module

  implicit none

  private
  public :: is_ibc_stencil, make_ibc_stencil_fab, simple_ib_const,  stencil_multifab_copy

contains

  ! Is i'th fab of ss a special stencil for interior box with const coefficients?
  pure function is_ibc_stencil(ss,i) result(r)
    logical :: r
    integer, intent(in) :: i
    type(multifab), intent(in) :: ss
    r = fab_is_0d(ss%fbs(i))
  end function is_ibc_stencil

  subroutine make_ibc_stencil_fab(ss, i, dm)
    integer, intent(in) :: i, dm
    type(multifab), intent(inout) :: ss
    call fab_destroy(ss%fbs(i))
    call fab_build_0d(ss%fbs(i), dm+1)
  end subroutine make_ibc_stencil_fab

  ! interior box, const coefficient 
  ! simple cross stencil
  subroutine simple_ib_const(ss, alpha_const, beta_const, dh, dm)
    real (kind = dp_t), intent(  out) ::   ss(0:)
    real (kind = dp_t), intent(in   ) :: alpha_const, beta_const
    real (kind = dp_t), intent(in   ) :: dh(:)
    integer           , intent(in   ) :: dm
    real (kind = dp_t) :: f1(dm)
    f1 = ONE/dh**2
    ss(0) = TWO*beta_const*sum(f1) + alpha_const
    ss(1:dm) = -beta_const*f1
  end subroutine simple_ib_const

  subroutine stencil_multifab_copy(dst, src)
    type(multifab), intent(in) :: src
    type(multifab), intent(inout) :: dst
    integer :: i, dm
    double precision, pointer :: pdst(:,:,:,:), psrc(:,:,:,:)
    dm = get_dim(dst)
    do i = 1, nfabs(dst)
       if (is_ibc_stencil(src,i)) then
          call make_ibc_stencil_fab(dst, i, dm)
       end if
    end do
    !$omp parallel do private(i,pdst,psrc)
    do i = 1, nfabs(dst)
       pdst => dataptr(dst, i)
       psrc => dataptr(src, i)
       call cpy_d(pdst, psrc)
    end do
    !$omp end parallel do
  end subroutine stencil_multifab_copy

  ! subroutine copy_from_ibc(dst, src)
  !   type(multifab), intent(in) :: src
  !   type(multifab), intent(inout) :: dst
  !   integer :: i, dm
  !   double precision, pointer :: pdst(:,:,:,:), psrc(:,:,:,:)
  !   dm = get_dim(dst)
  !   !$omp parallel do private(i,pdst,psrc)
  !   do i = 1, nfabs(dst)
  !      pdst => dataptr(dst, i)
  !      psrc => dataptr(src, i)
  !      if (is_ibc_stencil(src,i)) then
  !         select case (dm)
  !         case (1)
  !            call bl_error("stencil_util::copy_from_ibc: 1d is not supported")
  !         case (2)
  !            call copy_from_ibc_to_simple_2d_const(pdst(:,:,:,1), psrc(:,1,1,1))
  !         case (3)
  !            call copy_from_ibc_to_simple_3d_const(pdst(:,:,:,:), psrc(:,1,1,1))
  !         end select
  !      else
  !         call cpy_d(pdst, psrc)
  !      end if
  !   end do
  !   !$omp end parallel do
  ! end subroutine copy_from_ibc

  ! subroutine copy_from_ibc_to_simple_2d_const(sdst, ssrc)
  !   real (kind = dp_t), intent(in ) :: ssrc(0:)
  !   real (kind = dp_t), intent(out) :: sdst(0:,1:,1:)
  !   integer :: i, j, nx, ny
  !   nx = size(sdst,dim=2)
  !   ny = size(sdst,dim=3)
  !   ! see cc_stencil.f90 for simple_2d_const
  !   do j = 1, ny
  !      do i = 1, nx
  !         sdst(0,i,j) = ssrc(0)
  !         sdst(1,i,j) = ssrc(1)
  !         sdst(2,i,j) = ssrc(1)
  !         sdst(3,i,j) = ssrc(2)
  !         sdst(4,i,j) = ssrc(2)
  !         sdst(5:,i,j) = zero
  !      end do
  !   end do
  ! end subroutine copy_from_ibc_to_simple_2d_const

  ! subroutine copy_from_ibc_to_simple_3d_const(sdst, ssrc)
  !   real (kind = dp_t), intent(in ) :: ssrc(0:)
  !   real (kind = dp_t), intent(out) :: sdst(0:,1:,1:,1:)
  !   integer :: i, j, k, nx, ny, nz
  !   nx = size(sdst,dim=2)
  !   ny = size(sdst,dim=3)
  !   nz = size(sdst,dim=4)
  !   ! see cc_stencil.f90 for simple_3d_const
  !   do k = 1, nz
  !      do j = 1, ny
  !         do i = 1, nx
  !            sdst(0,i,j,k) = ssrc(0)
  !            sdst(1,i,j,k) = ssrc(1)
  !            sdst(2,i,j,k) = ssrc(1)
  !            sdst(3,i,j,k) = ssrc(2)
  !            sdst(4,i,j,k) = ssrc(2)
  !            sdst(5,i,j,k) = ssrc(3)
  !            sdst(6,i,j,k) = ssrc(3)
  !            sdst(7:,i,j,k) = zero
  !         end do
  !      end do
  !   end do
  ! end subroutine copy_from_ibc_to_simple_3d_const

end module stencil_util_module
