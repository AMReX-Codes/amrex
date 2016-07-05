!
! This program takes a 2-d plotfile and calculates the adiabatic excess
! along the vertical direction.  See MAESTRO/docs/thermo_notes for details.
!

program fad_excess

  use bl_space, only: MAX_SPACEDIM
  use bl_error_module
  use bl_constants_module, only: ZERO
  use bl_IO_module
  use plotfile_module
  use multifab_module
  use network
  use eos_module

  implicit none

  ! argument variables
  character(len=256) :: pltfile, outputfile
  real(kind=dp_t) :: low_cutoff

  ! f2kcli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local variables
  integer :: chk_int, ipos
  type(plotfile) :: pf
  type(fab) :: fb
  type(layout) :: la
  type(boxarray) :: ba
  type(list_box) :: bl
  type(box) :: bx,pd
  integer, allocatable :: ref_ratio(:)
  real(kind=dp_t),dimension(MAX_SPACEDIM) :: prob_lo, prob_hi
  real(kind=dp_t) :: time
  type(multifab), allocatable :: ad_excess(:)

  integer :: uin, uout
  integer :: dim, nlevs
  integer :: i, j, ii, jj
  real(kind=dp_t), dimension(MAX_SPACEDIM) :: dx
  integer, dimension(MAX_SPACEDIM) :: flo, fhi, lo, hi
  real(kind=dp_t), pointer :: p(:,:,:,:), ap(:,:,:,:)
  real(kind=dp_t), allocatable :: pres(:,:), nabla_ad(:,:)
  real(kind=dp_t) :: chi_rho, chi_t, dp, dt, nabla

  integer :: dens_comp, spec_comp, temp_comp

  character(len=20) :: plot_name(1)

  logical :: use_eos_coulomb = .true., do_diag = .true.
  
  real(kind=dp_t), parameter :: small = 1.e-14


  uin = unit_new()

  ! defaults
  pltfile = ''
  outputfile = 'ad_excess'
  low_cutoff = 1.e4

  plot_name(1) = "adiabatic excess"

  ! parse the arguements
  narg = command_argument_count()

  farg = 1
  do while (farg <= narg) 
     call get_command_argument(farg, value = fname)

     select case(fname)
        
     case ('-i', '--input')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-o', '--output')
        farg = farg + 1
        call get_command_argument(farg, value = outputfile)

     case ('-l', '--low_cutoff')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) low_cutoff

     case default
        exit

     end select
     farg = farg + 1
  enddo

  ! sanity check
  if (pltfile == '') then
     call print_usage()
     stop
  endif

  print *, 'working on pltfile:', pltfile
  print *, 'dumping data to file: ', outputfile
  print *, 'using low density cutoff:', low_cutoff

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

  ! build the input plotfile
  call build(pf, pltfile, uin)

  nlevs = pf%flevel
  dim = pf%dim

  dens_comp = plotfile_var_index(pf, "density")
  spec_comp = plotfile_var_index(pf, "X(He4)")
  temp_comp = plotfile_var_index(pf, "tfromp")

  if (dens_comp < 0 .or. spec_comp < 0 .or. temp_comp < 0) &
       call bl_error("Variables not found")

  allocate(ref_ratio(nlevs-1),ad_excess(nlevs))

  ! get the index bounds for the finest level
  flo(1:dim) = lwb(plotfile_get_pd_box(pf, nlevs))
  fhi(1:dim) = upb(plotfile_get_pd_box(pf, nlevs))

  allocate(pres(flo(1):fhi(1),flo(2):fhi(2)))
  allocate(nabla_ad(flo(1):fhi(1),flo(2):fhi(2)))

  ! define the ref ratios
  do i = 1, nlevs-1
     ref_ratio(i) = pf%refrat(i,1)
  enddo

  do i = 1, pf%flevel
     do j = 1, nboxes(pf,i)
        call push_back(bl,get_box(pf,i,j))
     enddo

     call build(ba,bl)

     call layout_build_ba(la,ba)
     
     call destroy(bl)
     call destroy(ba)

     ! build the mutifab with 0 ghost cells and 1 component
     call multifab_build(ad_excess(i),la,1,0)

  enddo
     
  ! loop over the plotfile data starting at the first
  do i = pf%flevel, 1, -1

     ! loop over each box at this level
     do j = 1, nboxes(pf, i)

        ! read in the data 1 patch at a time
        call fab_bind(pf, i, j)

        ! get the integer bounds of the current box, in terms of this 
        ! level's index space
        lo(1:dim) = lwb(get_box(pf,i,j))
        hi(1:dim) = upb(get_box(pf,i,j))

        ! pointers to the plotfile data and the ad_excess data
        p => dataptr(pf, i, j)
        ap => dataptr(ad_excess(i), j)!, get_box(ad_excess(i),j))

        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)

              den_eos(1) = p(ii,jj,1,dens_comp)
              temp_eos(1) = p(ii,jj,1,temp_comp)
              xn_eos(1,:) = p(ii,jj,1,spec_comp:spec_comp+nspec-1)

              call eos(eos_input_rt, den_eos, temp_eos, &
                          npts, &
                          xn_eos, &
                          p_eos, h_eos, e_eos, &
                          cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                          dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                          dpdX_eos, dhdX_eos, &
                          gam1_eos, cs_eos, s_eos, &
                          dsdt_eos, dsdr_eos, do_diag)

              pres(ii,jj) = p_eos(1)

              chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
              chi_t   = temp_eos(1) * dpdt_eos(1) / p_eos(1)
              nabla_ad(ii,jj) = (gam1_eos(1) - chi_rho) / (gam1_eos(1)*chi_t)

           enddo
        enddo

        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)

              ! forward difference
              if (jj==lo(2)) then
                 dt = p(ii,jj+1,1,temp_comp) - p(ii,jj,1,temp_comp)
                 dp = pres(ii,jj+1) - pres(ii,jj)
              ! backward difference
              else if (jj == hi(2)) then
                 dt = p(ii,jj,1,temp_comp) - p(ii,jj-1,1,temp_comp)
                 dp = pres(ii,jj) - pres(ii,jj-1)
              ! centered difference
              else
                 dt = p(ii,jj+1,1,temp_comp) - p(ii,jj-1,1,temp_comp)
                 dp = pres(ii,jj+1) - pres(ii,jj-1)
              endif

              if (p(ii,jj,1,dens_comp) <= low_cutoff .or. abs(dp) <= small) then
                 nabla = ZERO
              else
                 nabla = pres(ii,jj) * dt / (dp*p(ii,jj,1,temp_comp))
              endif
              
              ap(ii,jj,1,1) = nabla - nabla_ad(ii,jj)

           enddo
        enddo

        call fab_unbind(pf, i,j)

     enddo
  enddo

  pd = plotfile_get_pd_box(pf,1)
  dx(1:dim) = plotfile_get_dx(pf,1)
  prob_lo(1:dim) = pf%plo
  prob_hi(1:dim) = pf%phi
  time = plotfile_time(pf)

  call fabio_ml_multifab_write_d(ad_excess, ref_ratio, trim(outputfile), &
                                 plot_name, pd, &
                                 prob_lo(1:dim), prob_hi(1:dim), time, &
                                 dx(1:dim))

  call destroy(pf)

contains
  subroutine print_usage()
    implicit none

    print *, 'This program takes a 2-d plotfile and calculates the adiabatic'
    print *, 'excess which is dumped as a 2-d plotfile.  See '
    print *, 'MAESTRO/docs/thermo_notes for details.'
    print *, ''
    print *, ' calling sequence: '
    print *, ' fad_excess -i <inputfile> [-o <outputfile>, -l <cutoff>]'
    print *, ''
    print *, ' arguments:'
    print *, '-i | --input:'
    print *, '     specify the input plotfile used to calculate the adiabatic'
    print *, '     excess.  required.'
    print *, '-o | --output:'
    print *, '     specify the output filename for the adiabatic excess.'
    print *, '     defaults to "ad_excess"'
    print *, '-l | --low_cutoff:'
    print *, '     specifies the low density cutoff.  for densities below the'
    print *, '     low density cutoff, the actual thermal gradient is set to'
    print *, '     ZERO.  defaults to 1.e4.'
    print *, ''

  end subroutine print_usage

end program fad_excess
