! Analysis routine for RHD_shocktube

program frhdshocktube

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) pf
  integer :: i, j, ig

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  character(len=256) :: pltfile
  character(len=256) :: groupfile
  integer :: unit, ung, uno

  integer :: narg, farg, f

  integer :: rad0_comp, pres_comp, dens_comp, velx_comp, ngroups

  character(len=256) :: header_line
  integer :: ipos

  integer :: nx
  real(kind=dp_t), allocatable :: rad(:), pres(:), dens(:), velx(:), x(:)

  pltfile  = ''
  groupfile = "group_structure.dat"

  ! all the arguments are assumed to be plotfiles
  narg = command_argument_count()
  farg = 1

  call get_command_argument(1, value = pltfile)
  
  unit = unit_new()
  call build(pf, pltfile, unit)
  
  ! find the index of the first radiation energy
  rad0_comp = plotfile_var_index(pf, "rad0")
  dens_comp = plotfile_var_index(pf, "density")
  pres_comp = plotfile_var_index(pf, "pressure")
  velx_comp = plotfile_var_index(pf, "x_velocity")
  
  ! open up the group file and read in the group information
  open(unit=ung,file=groupfile)

  ! read in the number of groups
  read(ung,fmt='(a256)') header_line
  ipos = index(header_line, "=") + 1
  read (header_line(ipos:),*) ngroups
  
  close(ung)
  
  if (pf%dim /= 1) call bl_error("ERROR: frhdshocktube only works for dim = 1")
  
  ! Figure out how many levels in this plotfile
  if (plotfile_nlevels(pf) /= 1) call bl_error("ERROR: frhdshocktube only works for single level") 
  
  ! get the index bounds and dx for the coarse level.  Note, lo and hi are
  ! ZERO based indicies
  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))
  dx = plotfile_get_dx(pf, 1)
  
  nx = hi(1) - lo(1)
  allocate(rad(lo(1):hi(1)))
  allocate(pres(lo(1):hi(1)))
  allocate(dens(lo(1):hi(1)))
  allocate(velx(lo(1):hi(1)))
  allocate(x(lo(1):hi(1)))

  do i=lo(1), hi(1)
     x(i) = (i+0.5d0) * dx(1)
  end do
  
  rad = 0.d0
  pres = 0.d0
  dens = 0.d0
  velx = 0.d0
  
  i = pf%flevel
  call fab_bind_level(pf, i)
  
  do j = 1, nboxes(pf,i)
     
     lo = lwb(get_box(pf, i, j))
     hi = upb(get_box(pf, i, j))
     
     ! get a pointer to the current patch
     p => dataptr(pf, i, j)
     
     do ig = 0, ngroups-1
        rad(lo(1):hi(1)) = rad(lo(1):hi(1)) + p(lo(1):hi(1), 1, 1, rad0_comp+ig)
     end do
     pres(lo(1):hi(1)) = pres(lo(1):hi(1)) +  p(lo(1):hi(1), 1, 1, pres_comp)
     dens(lo(1):hi(1)) = dens(lo(1):hi(1)) +  p(lo(1):hi(1), 1, 1, dens_comp)
     velx(lo(1):hi(1)) = velx(lo(1):hi(1)) +  p(lo(1):hi(1), 1, 1, velx_comp)
     
  end do
  
  call destroy(pf)

1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  uno =  unit_new()
  open(unit=uno, file=trim(pltfile)//'.dat', status = 'replace')

  ! write the header
  write(uno,1000) "x", "density", "velocity", "pressure", "rad"

  do i=0, nx-1
     write(uno,1001) x(i), dens(i), velx(i), pres(i), rad(i)
  end do

  close(unit=uno)
  
end program frhdshocktube
