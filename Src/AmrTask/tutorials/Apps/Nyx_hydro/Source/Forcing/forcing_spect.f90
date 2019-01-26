! This module initializes and sets the Fourier modes for the computation of the physical force field

module forcing_spect_module

  use amrex_fort_module, only : bl_spacedim, rt => amrex_real

  implicit none

  integer, parameter :: LENVEC = 4  ! assumed minimal vector length

  integer :: num_modes
  integer :: num_modes_ext
  integer, allocatable :: wavevectors(:,:)
  real(rt), allocatable :: modes_even(:,:), modes_odd(:,:)

  public fort_alloc_spect, fort_set_wavevector, fort_set_modes

contains

  ! allocate arrays
  subroutine fort_alloc_spect(length) &
             bind(C, name="fort_alloc_spect")

    integer,  intent(in) :: length

    integer alloc, modvec

    if (length > 1) then
       num_modes = length

       ! determine extend array length as integer multiple of 4 
       ! to support vectorization of loops over all modes
       modvec = mod(length, LENVEC)
       if (modvec > 0) then
          num_modes_ext = num_modes + LENVEC - modvec
       else
          num_modes_ext = num_modes
       end if

       allocate (modes_even(num_modes_ext,bl_spacedim), modes_odd(num_modes_ext,bl_spacedim), &
                 wavevectors(bl_spacedim,num_modes_ext), STAT=alloc)
       if (alloc > 0) call bl_abort('failed to allocate arrays for forcing modes')
    else
       call bl_abort('number of forcing modes must be positive')
    end if

    modes_even(:,:) = 0.d0
    modes_even(:,:) = 0.d0
    wavevectors(:,:) = 0.d0

  end subroutine fort_alloc_spect

  ! set m-th nonzero wavevector of spectrum
  subroutine fort_set_wavevector(kvect, m) &
             bind(C, name="fort_set_wavevector")

    integer, intent(in) :: m
    integer, intent(in) :: kvect(bl_spacedim)

    if ((m.lt.0).or.(m.ge.num_modes)) then
	call bl_abort('invalid index of wavevector')
    end if

    wavevectors(:,m+1) = kvect(:) ! index conversion from C to Fortran

  end subroutine fort_set_wavevector

  ! set modes
  subroutine fort_set_modes(even, odd, length, comp) &
             bind(C, name="fort_set_modes")

    integer,  intent(in) :: length, comp
    real(rt), intent(in) :: even(length), odd(length)

    integer m

    if ((length.ne.num_modes).or.(comp.ge.bl_spacedim)) then
	call bl_abort('dimensions of input arrays do not match')
    end if

    modes_even(1:num_modes,comp+1) = even(:)
    modes_odd(1:num_modes,comp+1) = odd(:)

  end subroutine fort_set_modes

end module forcing_spect_module

