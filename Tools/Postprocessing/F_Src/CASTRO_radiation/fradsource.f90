! Analysis routine for the radiation source test.  
!
! This problem is a thermal relaxiation problem.  The domain is 
! completely uniform, so we just need to look at the state variables
! in a single zone.
!
! Take a list of files and print out (rho e) and the total radiation 
! energy density in the first zone as a function of time.

program fradsphere

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) pf
  integer :: i, j

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)

  character(len=256) :: pltfile
  integer :: unit

  integer :: narg, farg, f

  integer :: irad_begin, irhoe_begin

  pltfile  = ''

  ! all the arguments are assumed to be plotfiles
  narg = command_argument_count()
  farg = 1

  ! loop over all the plotfiles
  do f = farg, narg

     call get_command_argument(f, value = pltfile)

     unit = unit_new()
     call build(pf, pltfile, unit)

     ! find the index of the first radiation energy
     irad_begin = plotfile_var_index(pf, "rad")

     ! find the index of the gas internal energy density
     irhoe_begin = plotfile_var_index(pf, "rho_e")

     if (pf%dim /= 1) call bl_error("ERROR: fradsphere only works for dim = 1")

     ! only work with the finest level's data
     i = pf%flevel

     ! we only care about a single zone, so take the first box
     j = 1

     ! read in the data patch
     call fab_bind(pf, i, j)

     lo = lwb(get_box(pf, i, j))
     hi = upb(get_box(pf, i, j))

     p => dataptr(pf, i, j)

     if (f == farg) then
        print *, "# time, (rho e), rad"
     endif

     print *, pf%tm, p(lo(1),1,1,irhoe_begin), p(lo(1),1,1,irad_begin)

     call fab_unbind(pf, i, j)

     call destroy(pf)

  enddo

end program fradsphere
