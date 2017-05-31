module multifab_physbc_module

  use bl_constants_module
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: multifab_physbc

contains

  subroutine multifab_physbc(s,start_scomp,start_bccomp,ncomp,the_bc_level, &
                             time_in,dx_in,prob_lo_in,prob_hi_in)

    integer        , intent(in   )           :: start_scomp,start_bccomp,ncomp
    type(multifab) , intent(inout)           :: s
    type(bc_level) , intent(in   )           :: the_bc_level
    real(kind=dp_t), intent(in   ), optional :: time_in,dx_in(:),prob_lo_in(:),prob_hi_in(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng,dm,scomp,bccomp
    real(kind=dp_t)          :: time,dx(get_dim(s))
    real(kind=dp_t)          :: prob_lo(get_dim(s)),prob_hi(get_dim(s))
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (nghost(s) == 0) return
    
    ! set optional arguments
    time    = ZERO
    dx      = ZERO
    prob_lo = ZERO
    prob_hi = ZERO
    if (present(time_in))       time = time_in
    if (present(dx_in))           dx = dx_in
    if (present(prob_lo_in)) prob_lo = prob_lo_in
    if (present(prob_hi_in)) prob_hi = prob_hi_in

    ng = nghost(s)
    dm = get_dim(s)

    !$OMP PARALLEL DO PRIVATE(i,sp,lo,hi,scomp,bccomp)
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (1)
          do scomp=start_scomp,start_scomp+ncomp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_1d(sp(:,1,1,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                            time, dx, prob_lo, prob_hi)
          end do          
       case (2)
          do scomp=start_scomp,start_scomp+ncomp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                            time, dx, prob_lo, prob_hi)
          end do
       case (3)
          do scomp=start_scomp,start_scomp+ncomp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                            time, dx, prob_lo, prob_hi)
          end do
       end select
    end do
    !$OMP END PARALLEL DO
 
  end subroutine multifab_physbc


  !============================================================================
  ! physbc_1d
  !============================================================================
  subroutine physbc_1d(s,lo,hi,ng,bc,time,dx,prob_lo,prob_hi)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: time,dx(:),prob_lo(:),prob_hi(:)

    ! Local variables
    integer :: i

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) == EXT_DIR) then
       ! set all ghost cell values to a prescribed dirichlet
       ! value; in this example, we have chosen 1
       s(lo(1)-ng:lo(1)-1) = ONE

    else if (bc(1,1) == FOEXTRAP) then
       ! set all ghost cell values to first interior value
       s(lo(1)-ng:lo(1)-1) = s(lo(1))

    else if (bc(1,1) == HOEXTRAP) then
       ! set all ghost cell values equal to quadratic interpolate
       ! to the physical location of the domain boundary, NOT the
       ! physical location of the first ghost cell-center
       s(lo(1)-ng:lo(1)-1) = EIGHTH* &
            (FIFTEEN*s(lo(1)) - TEN*s(lo(1)+1) + THREE*s(lo(1)+2))

    else if (bc(1,1) == REFLECT_EVEN) then
       ! mirror the interior
       do i = 1, ng
          s(lo(1)-i) = s(lo(1)+i-1)
       end do

    else if (bc(1,1) == REFLECT_ODD) then
       ! mirror the interior with opposite sign
       do i = 1, ng
          s(lo(1)-i) = -s(lo(1)+i-1)
       end do

    else if (bc(1,1) == INTERIOR) then
       ! do nothing, these ghost cell values should be filled by something else
       ! multifab_fill_boundary fills from periodic or neighboring cells at
       ! the same level of refinement, multifab_fill_ghost_cells fills fine
       ! ghost cells by interpolating coarse values in space

    else 
       print *,'bc(1,1) = ',bc(1,1)
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------

    if (bc(1,2) == EXT_DIR) then
          s(hi(1)+1:hi(1)+ng) = ONE

    else if (bc(1,2) == FOEXTRAP) then
       s(hi(1)+1:hi(1)+ng) = s(hi(1))

    else if (bc(1,2) == HOEXTRAP) then
       s(hi(1)+1:hi(1)+ng) = EIGHTH* &
            (FIFTEEN*s(hi(1)) - TEN*s(hi(1)-1) + THREE*s(hi(1)-2))

    else if (bc(1,2) == REFLECT_EVEN) then
       do i = 1, ng
          s(hi(1)+i) = s(hi(1)-i+1)
       end do

    else if (bc(1,2) == REFLECT_ODD) then
       do i = 1, ng
          s(hi(1)+i) = -s(hi(1)-i+1)
       end do

    else if (bc(1,2) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if
    
  end subroutine physbc_1d


  !============================================================================
  ! physbc_2d
  !============================================================================
  subroutine physbc_2d(s,lo,hi,ng,bc,time,dx,prob_lo,prob_hi)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: time,dx(:),prob_lo(:),prob_hi(:)

    ! Local variables
    integer :: i,j


    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) == EXT_DIR) then
       ! set all ghost cell values to a prescribed dirichlet
       ! value; in this example, we have chosen 1
       do j = lo(2)-ng, hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = ONE
       end do

    else if (bc(1,1) == FOEXTRAP) then
       ! set all ghost cell values to first interior value
       do j = lo(2)-ng, hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do

    else if (bc(1,1) == HOEXTRAP) then
       ! set all ghost cell values equal to quadratic interpolate
       ! to the physical location of the domain boundary, NOT the
       ! physical location of the first ghost cell-center
       do j = lo(2)-ng, hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = EIGHTH* &
               (FIFTEEN*s(lo(1),j) - TEN*s(lo(1)+1,j) + THREE*s(lo(1)+2,j))
       end do

    else if (bc(1,1) == REFLECT_EVEN) then
       ! mirror the interior
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(lo(1)-i,j) = s(lo(1)+i-1,j)
          end do
       end do

    else if (bc(1,1) == REFLECT_ODD) then
       ! mirror the interior with opposite sign
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(lo(1)-i,j) = -s(lo(1)+i-1,j)
          end do
       end do

    else if (bc(1,1) == INTERIOR) then
       ! do nothing, these ghost cell values should be filled by something else
       ! multifab_fill_boundary fills from periodic or neighboring cells at
       ! the same level of refinement, multifab_fill_ghost_cells fills fine
       ! ghost cells by interpolating coarse values in space

    else 
       print *,'bc(1,1) = ',bc(1,1)
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) == EXT_DIR) then
       do j = lo(2)-ng, hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = ONE
       end do

    else if (bc(1,2) == FOEXTRAP) then
       do j = lo(2)-ng, hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do

    else if (bc(1,2) == HOEXTRAP) then
       do j = lo(2)-ng, hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = EIGHTH* &
               (FIFTEEN*s(hi(1),j) - TEN*s(hi(1)-1,j) + THREE*s(hi(1)-2,j))
       end do

    else if (bc(1,2) == REFLECT_EVEN) then
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(hi(1)+i,j) = s(hi(1)-i+1,j)
          end do
       end do

    else if (bc(1,2) == REFLECT_ODD) then
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(hi(1)+i,j) = -s(hi(1)-i+1,j)
          end do
       end do

    else if (bc(1,2) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------
    if (bc(2,1) == EXT_DIR) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = ONE
       end do

    else if (bc(2,1) == FOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do

    else if (bc(2,1) == HOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = EIGHTH* &
               (FIFTEEN*s(i,lo(2)) - TEN*s(i,lo(2)+1) + THREE*s(i,lo(2)+2))
       end do

    else if (bc(2,1) == REFLECT_EVEN) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,lo(2)-j) = s(i,lo(2)+j-1)
          end do
       end do

    else if (bc(2,1) == REFLECT_ODD) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,lo(2)-j) = -s(i,lo(2)+j-1)
          end do
       end do

    else if (bc(2,1) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(2,1) = ',bc(2,1)
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------
    if (bc(2,2) == EXT_DIR) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = ONE
       end do

    else if (bc(2,2) == FOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do

    else if (bc(2,2) == HOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = EIGHTH* &
               (FIFTEEN*s(i,hi(2)) - TEN*s(i,hi(2)-1) + THREE*s(i,hi(2)-2))
       end do

    else if (bc(2,2) == REFLECT_EVEN) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,hi(2)+j) = s(i,hi(2)-j+1)
          end do
       end do

    else if (bc(2,2) == REFLECT_ODD) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,hi(2)+j) = -s(i,hi(2)-j+1)
          end do
       end do

    else if (bc(2,2) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(2,2) = ',bc(2,2)
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

  end subroutine physbc_2d


  !============================================================================
  ! physbc_3d
  !============================================================================
  subroutine physbc_3d(s,lo,hi,ng,bc,time,dx,prob_lo,prob_hi)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: time,dx(:),prob_lo(:),prob_hi(:)

    ! Local variables
    integer :: i,j,k

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) == EXT_DIR) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = ONE
          end do
       end do

    else if (bc(1,1) == FOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = s(lo(1),j,k)
          end do
       end do

    else if (bc(1,1) == HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = EIGHTH* &
                  (FIFTEEN*s(lo(1),j,k) - TEN*s(lo(1)+1,j,k) + THREE*s(lo(1)+2,j,k))
          end do
       end do

    else if (bc(1,1) == REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(lo(1)-i,j,k) = s(lo(1)+i-1,j,k)
             end do
          end do
       end do

    else if (bc(1,1) == REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(lo(1)-i,j,k) = -s(lo(1)+i-1,j,k)
             end do
          end do
       end do

    else if (bc(1,1) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(1,1) = ',bc(1,1)
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) == EXT_DIR) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = ONE
          end do
       end do

    else if (bc(1,2) == FOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k)
          end do
       end do

    else if (bc(1,2) == HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = EIGHTH* &
                  (FIFTEEN*s(hi(1),j,k) - TEN*s(hi(1)-1,j,k) + THREE*s(hi(1)-2,j,k))
          end do
       end do

    else if (bc(1,2) == REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(hi(1)+i,j,k) = s(hi(1)-i+1,j,k)
             end do
          end do
       end do

    else if (bc(1,2) == REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(hi(1)+i,j,k) = -s(hi(1)-i+1,j,k)
             end do
          end do
       end do

    else if (bc(1,2) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------
    if (bc(2,1) == EXT_DIR) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = ONE
          end do
       end do

    else if (bc(2,1) == FOEXTRAP .or. bc(2,1) == REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k)
          end do
       end do

    else if (bc(2,1) == HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = EIGHTH* &
                  (FIFTEEN*s(i,lo(2),k) - TEN*s(i,lo(2)+1,k) + THREE*s(i,lo(2)+2,k))
          end do
       end do

    else if (bc(2,1) == REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,lo(2)-j,k) = s(i,lo(2)+j-1,k)
             end do
          end do
       end do

    else if (bc(2,1) == REFLECT_ODD) then
       do k = lo(3)-1,hi(3)+1
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,lo(2)-j,k) = -s(i,lo(2)+j-1,k)
             end do
          end do
       end do

    else if (bc(2,1) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(2,1) = ',bc(2,1)
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------
    if (bc(2,2) == EXT_DIR) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = ONE
          end do
       end do

    else if (bc(2,2) == FOEXTRAP .or. bc(2,2) == REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k)
          end do
       end do

    else if (bc(2,2) == HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = EIGHTH* &
                  (FIFTEEN*s(i,hi(2),k) - TEN*s(i,hi(2)-1,k) + THREE*s(i,hi(2)-2,k))
          end do
       end do

    else if (bc(2,2) == REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,hi(2)+j,k) = s(i,hi(2)-j+1,k)
             end do
          end do
       end do

    else if (bc(2,2) == REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,hi(2)+j,k) = -s(i,hi(2)-j+1,k)
             end do
          end do
       end do

    else if (bc(2,2) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(2,2) = ',bc(2,2)
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! lower Z
    !--------------------------------------------------------------------------
    if (bc(3,1) == EXT_DIR) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = ONE
          end do
       end do

    else if (bc(3,1) == FOEXTRAP .or. bc(3,1) == REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do

    else if (bc(3,1) == HOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = EIGHTH* &
                  (FIFTEEN*s(i,j,lo(3)) - TEN*s(i,j,lo(3)+1) + THREE*s(i,j,lo(3)+2))
          end do
       end do

    else if (bc(3,1) == REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,lo(3)-k) = s(i,j,lo(3)+k-1)
             end do
          end do
       end do

    else if (bc(3,1) == REFLECT_ODD) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,lo(3)-k) = -s(i,j,lo(3)+k-1)
             end do
          end do
       end do

    else if (bc(3,1) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(3,1) = ',bc(3,1)
       call bl_error('BC(3,1) = NOT YET SUPPORTED')
    end if


    !--------------------------------------------------------------------------
    ! upper Z
    !--------------------------------------------------------------------------
    if (bc(3,2) == EXT_DIR) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = ONE
          end do
       end do

    else if (bc(3,2) == FOEXTRAP .or. bc(3,2) == REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do

    else if (bc(3,2) == HOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = EIGHTH* &
                  (FIFTEEN*s(i,j,hi(3)) - TEN*s(i,j,hi(3)-1) + THREE*s(i,j,hi(3)-2))
          end do
       end do

    else if (bc(3,2) == REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,hi(3)+k) = s(i,j,hi(3)-k+1)
             end do
          end do
       end do

    else if (bc(3,2) == REFLECT_ODD) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,hi(3)+k) = -s(i,j,hi(3)-k+1)
             end do
          end do
       end do

    else if (bc(3,2) == INTERIOR) then
       ! do nothing - see comment above

    else 
       print *,'bc(3,2) = ',bc(3,2)
       call bl_error('BC(3,2) = NOT YET SUPPORTED')
    end if

  end subroutine physbc_3d

end module multifab_physbc_module
