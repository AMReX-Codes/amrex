program unittests

  use BoxLib
  use bl_types
  use multifab_module
  use ml_layout_module
  use ml_boxarray_module
  use box_util_module

  implicit none

  real(kind=dp_t), pointer :: dx(:,:)

  type(ml_layout) :: mla
  type(ml_boxarray) :: mba

  type(multifab) :: a

  character (len=256) :: test_set

  logical :: pmask(3)

  real(dp_t) :: prob_lo(3), prob_hi(3)

  integer :: n

  !--------------------------------------------------------------------------
  ! initialize the grid
  !--------------------------------------------------------------------------
  call boxlib_initialize()

  test_set = "gr0_3d.1level"

  pmask(:) = [.true., .true., .true.]
  prob_lo(:) = [0.0_dp_t, 0.0_dp_t, 0.0_dp_t]
  prob_hi(:) = [10.0_dp_t, 10.0_dp_t, 10.0_dp_t]

  call read_a_hgproj_grid(mba, test_set)
  call ml_layout_build(mla, mba, pmask)

  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('ERROR: fixed_grids not properly nested')
  end if

  call initialize_dx(dx, mba, mla%nlevel, prob_lo, prob_hi)


  !--------------------------------------------------------------------------
  ! start of unit tests
  !--------------------------------------------------------------------------
  do n = 1, mla%nlevel
     call build(a, mla%la(n), 1, 0)
  enddo

  call setval(a, 1.0_dp_t)



  call boxlib_finalize()


contains

  subroutine initialize_dx(dx, mba, num_levs, prob_lo, prob_hi)

    use ml_boxarray_module
    use bl_types

    real(kind=dp_t),   intent(inout), pointer :: dx(:,:)

    type(ml_boxarray), intent(in ) :: mba
    integer          , intent(in ) :: num_levs
    real(dp_t)       , intent(in ) :: prob_lo(3), prob_hi(3)

    integer :: n,d,dm

    dm = mba%dim

    allocate(dx(num_levs,dm))

    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mba%pd(1),d),kind=dp_t)
    end do
    do n=2,num_levs
       dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
    end do
    
  end subroutine initialize_dx

end program unittests
