! Take 2 plotfiles as input.  We assume that the first plotfile is 
! the earlier time.
!
! Make sure the grids and variables in each plotfile agree with
! one-another.
!
! Loop over all the patches in the first plotfile.  Find the zone
! with the maximum temperature and find the zones corresponding to
! the center -- store the grid location information for these.
! 
! Now get T, rho, and X for the hottest and center zones for each
! plotfile, compute c_p, and report the change in T and c_p over
! the time interval between the 2 plotfiles.

program fthermo

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module
  use network
  use eos_module

  implicit none

  type grid_info
     integer :: level = -1
     integer :: box = -1
     integer :: coords(3)
  end type grid_info

  type(plotfile) :: pf_a, pf_b
  character (len=256) :: plotfile_a, plotfile_b
  integer :: unit_a, unit_b

  real(kind=dp_t), pointer :: p_a(:,:,:,:), p_b(:,:,:,:)

  integer :: lo_a(MAX_SPACEDIM), hi_a(MAX_SPACEDIM)
  integer :: lo_b(MAX_SPACEDIM), hi_b(MAX_SPACEDIM)

  integer :: flo_a(MAX_SPACEDIM), fhi_a(MAX_SPACEDIM)
  integer :: flo_b(MAX_SPACEDIM), fhi_b(MAX_SPACEDIM)

  real(kind=dp_t) :: dx_a(MAX_SPACEDIM), dx_b(MAX_SPACEDIM)

  integer :: nboxes_a, nboxes_b

  integer :: n_a, n_b
  integer, allocatable :: ivar_b(:)

  integer :: narg, farg
  character (len=256) :: fname

  integer :: i, j
  integer :: ii, jj, kk

  real(kind=dp_t) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(kind=dp_t) :: x, y, z, xc, yc, zc

  type(grid_info) :: T_max_loc, T_center_loc(8)

  integer :: dens_comp, temp_comp, XC12_comp, XO16_comp, Xash_comp
  integer :: ic12, io16, iash

  real(kind=dp_t) :: dens_cutoff, T_ignite

  real(kind=dp_t) :: T_hot_a, T_hot_b, T_center_a, T_center_b
  real(kind=dp_t) :: rho_hot_a, rho_hot_b, rho_center_a, rho_center_b
  real(kind=dp_t) :: X_hot_a(nspec), X_hot_b(nspec)
  real(kind=dp_t) :: X_center_a(nspec), X_center_b(nspec)

  real(kind=dp_t) :: T_max

  real(kind=dp_t) :: cp_hot_a, cp_hot_b, cp_center_a, cp_center_b

  real(kind=dp_t) :: dist

  integer :: ncenter
  integer :: n


  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  plotfile_a = ""
  plotfile_b = ""

  dens_cutoff = 3.e6_dp_t
  T_ignite = 8.e8_dp_t

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)
     
     select case (fname)

     case ('--infile1')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_a)

     case ('--infile2')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_b)

     case default
        exit

     end select
     farg = farg + 1
  enddo

  if (len_trim(plotfile_a) == 0 .OR. len_trim(plotfile_b) == 0) then
     print *, " "
     print *, "Compare two plotfiles, zone by zone, to machine precision"
     print *, "and report the maximum absolute and relative errors for each"
     print *, "variable."
     print *, " "
     print *, "usage:"
     print *, "   fcompare --infile1 file1 --infile2 file2"
     print *, " "
     stop
  endif

  !---------------------------------------------------------------------------
  ! build the plotfiles and do initial comparisons
  
  unit_a = unit_new()
  call build(pf_a, plotfile_a, unit_a)

  unit_b = unit_new()
  call build(pf_b, plotfile_b, unit_b)

  
  ! check if they are the same dimensionality
  if (pf_a%dim /= pf_b%dim) then
     call bl_error("ERROR: plotfiles have different numbers of spatial dimensions")
  endif


  ! 3-d only
  if (pf_a%dim /= 3) then
     call bl_error("ERROR: this routine works on 3-d files only")
  endif


  ! check if they have the same number of levels
  if (pf_a%flevel /= pf_b%flevel) then
     call bl_error("ERROR: number of levels do not match")
  endif


  ! check if the finest domains are the same size
  flo_a = lwb(plotfile_get_pd_box(pf_a, pf_a%flevel))
  fhi_a = upb(plotfile_get_pd_box(pf_a, pf_a%flevel))

  flo_b = lwb(plotfile_get_pd_box(pf_b, pf_b%flevel))
  fhi_b = upb(plotfile_get_pd_box(pf_b, pf_b%flevel))

  if ( (flo_a(1) /= flo_b(1) .OR. fhi_a(1) /= fhi_b(1)) .OR. &
       (flo_a(2) /= flo_b(2) .OR. fhi_a(2) /= fhi_b(2)) .OR. &
       (flo_a(3) /= flo_b(3) .OR. fhi_a(3) /= fhi_b(3)) ) then
     call bl_error("ERROR: grids do not match")
  endif

  
  ! check if they have the same number of variables
  if (pf_a%nvars /= pf_b%nvars) then
     print *, " "
     print *, "WARNING: number of variables do not match"
  endif

  
  ! check if the variables are in the same order
  do n_a = 1, pf_a%nvars
     if (pf_a%names(n_a) /= pf_b%names(n_a)) then
        call bl_error("ERROR: variable ordering does not agree")
     endif
  enddo


  ! get the indices for the variables we care about 

  ! density
  dens_comp = plotfile_var_index(pf_a, "density")

  ! temperature (here we used T from p0)
  temp_comp = plotfile_var_index(pf_a, "tfromp")

  ! species
  XC12_comp = plotfile_var_index(pf_a, "X(C12)")
  XO16_comp = plotfile_var_index(pf_a, "X(O16)")
  Xash_comp = plotfile_var_index(pf_a, "X(ash)")


  if ( dens_comp < 0 .or. temp_comp < 0 .or. &
       XC12_comp < 0 .or. XO16_comp < 0 .or. Xash_comp < 0) then
     call bl_error("ERROR: varaible(s) not defined")
  endif



  ! initialize the network and eos
  call network_init()
  call eos_init(use_eos_coulomb=.true.)

  
  ! we are assuming a specific composition here -- make sure the right
  ! network is in place
  if (nspec /= 3) then
     call bl_error("ERROR: nspec incorrect")
  endif


  ! get the indices that the network uses to address the species
  ic12  = network_species_index("carbon-12")
  io16  = network_species_index("oxygen-16")
  iash = network_species_index("ash")



  ! get the domain extrema and compute the center
  xmin = pf_a%plo(1)
  xmax = pf_a%phi(1)
  xc = HALF*(xmin + xmax)

  ymin = pf_a%plo(2)
  ymax = pf_a%phi(2)
  yc = HALF*(ymin + ymax)

  zmin = pf_a%plo(3)
  zmax = pf_a%phi(3)
  zc = HALF*(zmin + zmax)
  


  !---------------------------------------------------------------------------
  ! go level-by-level and patch-by-patch and find the location of the max
  ! T and the zones that correspond to the center of the star.

  T_max = 0.0
  ncenter = 0

  do i = 1, pf_a%flevel


     ! make sure the dx's agree
     dx_a = plotfile_get_dx(pf_a, i)
     dx_b = plotfile_get_dx(pf_b, i)  

     if ((dx_a(1) /= dx_b(1)) .OR. &
         (dx_a(2) /= dx_b(2)) .OR. &
         (dx_a(3) /= dx_b(3))) then
        call bl_error("ERROR: grid dx does not match")
     endif


     ! make sure the number of boxes agree
     nboxes_a = nboxes(pf_a, i)
     nboxes_b = nboxes(pf_b, i)
     
     if (nboxes_a /= nboxes_b) then
        call bl_error("ERROR: number of boxes do not match")
     endif

     do j = 1, nboxes_a


        ! make sure that the grids match
        lo_a = lwb(get_box(pf_a, i, j))
        hi_a = upb(get_box(pf_a, i, j))

        lo_b = lwb(get_box(pf_b, i, j))
        hi_b = upb(get_box(pf_b, i, j))

        if ( (lo_a(1) /= lo_b(1) .OR. hi_a(1) /= hi_b(1)) .OR. &
             (lo_a(2) /= lo_b(2) .OR. hi_a(2) /= hi_b(2)) .OR. &
             (lo_a(3) /= lo_b(3) .OR. hi_a(3) /= hi_b(3)) ) then
           call bl_error("ERROR: grids do not match")
        endif
        

        ! get the data for the current box
        call fab_bind(pf_a, i, j)
        p_a => dataptr(pf_a, i, j)


        ! loop over all the cells in the current box
        do kk = lo_a(3), hi_a(3)
           z = zmin + dble(kk + HALF)*dx_a(3)

           do jj = lo_a(2), hi_a(2)
              y = ymin + dble(jj + HALF)*dx_a(2)

              do ii = lo_a(1), hi_a(1)
                 x = xmin + dble(ii + HALF)*dx_a(1)


                 dist = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
                 

                 ! if this the maximum temperature?
                 if (p_a(ii,jj,kk,temp_comp) > T_max .and. p_a(ii,jj,kk,dens_comp) > dens_cutoff) then
                    T_max = p_a(ii,jj,kk,temp_comp)
                    
                    T_max_loc%level = i
                    T_max_loc%box = j
                    T_max_loc%coords(1) = ii
                    T_max_loc%coords(2) = jj
                    T_max_loc%coords(3) = kk
                 endif

                 ! are we at the center (on the finest grid)?
                 if (i == pf_a%flevel .and. dist < dx_a(1)) then
                    ncenter = ncenter + 1
                    
                    T_center_loc(ncenter)%level = i
                    T_center_loc(ncenter)%box = j
                    T_center_loc(ncenter)%coords(1) = ii
                    T_center_loc(ncenter)%coords(2) = jj
                    T_center_loc(ncenter)%coords(3) = kk
                 endif                   

              enddo
           enddo
        enddo

        call fab_unbind(pf_a, i, j)

     enddo  ! boxes loop
     
  enddo  ! levels loop


  ! only consider pre-ignition
  if (T_max > T_ignite) then
     call bl_error("ERROR: ignited")
  endif

  ! we should have found 8 center zones
  if (ncenter /= 8) then
     call bl_error("ERROR: ncenter /= 8")
  endif

  !-------------------------------------------------------------------------
  ! compute the difference in T and c_p for the hottest zone from
  ! plotfile a to b.
  !-------------------------------------------------------------------------

  ! a
  call fab_bind(pf_a, T_max_loc%level, T_max_loc%box)
  p_a => dataptr(pf_a, T_max_loc%level, T_max_loc%box)

  T_hot_a    = p_a(T_max_loc%coords(1), &
                   T_max_loc%coords(2), &
                   T_max_loc%coords(3), temp_comp)

  rho_hot_a  = p_a(T_max_loc%coords(1), &
                   T_max_loc%coords(2), &
                   T_max_loc%coords(3), dens_comp)

  X_hot_a(ic12) = p_a(T_max_loc%coords(1), &
                      T_max_loc%coords(2), &
                      T_max_loc%coords(3), XC12_comp)

  X_hot_a(io16) = p_a(T_max_loc%coords(1), &
                      T_max_loc%coords(2), &
                      T_max_loc%coords(3), XO16_comp)

  X_hot_a(iash) = p_a(T_max_loc%coords(1), &
                      T_max_loc%coords(2), &
                      T_max_loc%coords(3), Xash_comp)

  ! call the EOS
  temp_eos(1) = T_hot_a
  den_eos(1)  = rho_hot_a
  xn_eos(1,:) = X_hot_a(:)

  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)

  cp_hot_a = cp_eos(1)

  call fab_unbind(pf_a, T_max_loc%level, T_max_loc%box)


  ! b
  call fab_bind(pf_b, T_max_loc%level, T_max_loc%box)
  p_b => dataptr(pf_b, T_max_loc%level, T_max_loc%box)

  T_hot_b    = p_b(T_max_loc%coords(1), &
                   T_max_loc%coords(2), &
                   T_max_loc%coords(3), temp_comp)

  rho_hot_b  = p_b(T_max_loc%coords(1), &
                   T_max_loc%coords(2), &
                   T_max_loc%coords(3), dens_comp)

  X_hot_b(ic12) = p_b(T_max_loc%coords(1), &
                      T_max_loc%coords(2), &
                      T_max_loc%coords(3), XC12_comp)

  X_hot_b(io16) = p_b(T_max_loc%coords(1), &
                      T_max_loc%coords(2), &
                      T_max_loc%coords(3), XO16_comp)

  X_hot_b(iash) = p_b(T_max_loc%coords(1), &
                      T_max_loc%coords(2), &
                      T_max_loc%coords(3), Xash_comp)

  ! call the EOS
  temp_eos(1) = T_hot_b
  den_eos(1)  = rho_hot_b
  xn_eos(1,:) = X_hot_b(:)

  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)

  cp_hot_b = cp_eos(1)

  call fab_unbind(pf_b, T_max_loc%level, T_max_loc%box)


  !-------------------------------------------------------------------------  
  ! compute the thermodynamic state at the domain center by averaging and
  ! find c_p there
  !-------------------------------------------------------------------------

  ! a
  T_center_a   = ZERO
  rho_center_a = ZERO
  X_center_a   = ZERO

  do n = 1, 8
     call fab_bind(pf_a, T_center_loc(n)%level, T_center_loc(n)%box)
     p_a => dataptr(pf_a, T_center_loc(n)%level, T_center_loc(n)%box)

     T_center_a   = T_center_a   + p_a(T_center_loc(n)%coords(1), &
                                       T_center_loc(n)%coords(2), &
                                       T_center_loc(n)%coords(3), temp_comp)

     rho_center_a = rho_center_a + p_a(T_center_loc(n)%coords(1), &
                                       T_center_loc(n)%coords(2), &
                                       T_center_loc(n)%coords(3), dens_comp)

     X_center_a(ic12) = X_center_a(ic12) + p_a(T_center_loc(n)%coords(1), &
                                               T_center_loc(n)%coords(2), &
                                               T_center_loc(n)%coords(3), XC12_comp)

     X_center_a(io16) = X_center_a(io16) + p_a(T_center_loc(n)%coords(1), &
                                               T_center_loc(n)%coords(2), &
                                               T_center_loc(n)%coords(3), XO16_comp)

     X_center_a(iash) = X_center_a(iash) + p_a(T_center_loc(n)%coords(1), &
                                               T_center_loc(n)%coords(2), &
                                               T_center_loc(n)%coords(3), Xash_comp)

     call fab_unbind(pf_a, T_center_loc(n)%level, T_center_loc(n)%box)
  enddo

  ! compute the average central state
  T_center_a    = T_center_a/8
  rho_center_a  = rho_center_a/8
  X_center_a(:) = X_center_a(:)/8

  temp_eos(1) = T_center_a
  den_eos(1)  = rho_center_a
  xn_eos(1,:) = X_center_a(:)

  ! call the EOS
  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)

  cp_center_a = cp_eos(1)
  
  
  ! b
  T_center_b   = ZERO
  rho_center_b = ZERO
  X_center_b   = ZERO

  do n = 1, 8
     call fab_bind(pf_b, T_center_loc(n)%level, T_center_loc(n)%box)
     p_b => dataptr(pf_b, T_center_loc(n)%level, T_center_loc(n)%box)

     T_center_b   = T_center_b   + p_b(T_center_loc(n)%coords(1), &
                                       T_center_loc(n)%coords(2), &
                                       T_center_loc(n)%coords(3), temp_comp)

     rho_center_b = rho_center_b + p_b(T_center_loc(n)%coords(1), &
                                       T_center_loc(n)%coords(2), &
                                       T_center_loc(n)%coords(3), dens_comp)

     X_center_b(ic12) = X_center_b(ic12) + p_b(T_center_loc(n)%coords(1), &
                                               T_center_loc(n)%coords(2), &
                                               T_center_loc(n)%coords(3), XC12_comp)

     X_center_b(io16) = X_center_b(io16) + p_b(T_center_loc(n)%coords(1), &
                                               T_center_loc(n)%coords(2), &
                                               T_center_loc(n)%coords(3), XO16_comp)

     X_center_b(iash) = X_center_b(iash) + p_b(T_center_loc(n)%coords(1), &
                                               T_center_loc(n)%coords(2), &
                                               T_center_loc(n)%coords(3), Xash_comp)

     call fab_unbind(pf_b, T_center_loc(n)%level, T_center_loc(n)%box)
  enddo

  ! compute the average central state
  T_center_b    = T_center_b/8
  rho_center_b  = rho_center_b/8
  X_center_b(:) = X_center_b(:)/8

  temp_eos(1) = T_center_b
  den_eos(1)  = rho_center_b
  xn_eos(1,:) = X_center_b(:)

  ! call the EOS
  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)

  cp_center_b = cp_eos(1)


  !-------------------------------------------------------------------------
  ! output
  !-------------------------------------------------------------------------
1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  write (*,1000), "time (file1)", "time (file2)", "T center (file1)", "cp center (file1)", "T center (file2)", "cp center (file2)", "peak T (file1)", "cp @ peak T (file1)", "cp @ peak T (file2)"
  write (*,1001) pf_a%tm, pf_b%tm, T_center_a, cp_center_a, T_center_b, cp_center_b, T_hot_a, cp_hot_a, cp_hot_b


end program fthermo
