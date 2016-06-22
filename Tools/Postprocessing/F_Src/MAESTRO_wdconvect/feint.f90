! Analysis routines for the 3-d WD convect MAESTRO problem 
!
! This version is only computes the internal energy.
!
! Note: it will assume that C, O, and Mg are the only species
!
program feint

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use eos_module
  use network
  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj, kk
  real(kind=dp_t) :: xx, yy, zz
  integer :: rr, r1
  integer :: uno

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  real(kind=dp_t), pointer :: p(:,:,:,:)

  real(kind=dp_t) :: internal_energy

  integer :: dens_comp, temp_comp, XC12_comp, XO16_comp, XMg24_comp
  integer :: ic12, io16, img24

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  character(len=256) :: pltfile
  real(kind=dp_t) :: rho_cutoff

  integer :: narg, farg
  character(len=256) :: fname

  real(kind=dp_t) :: time

  unit = unit_new()
  uno =  unit_new()


  ! set the defaults
  pltfile  = ''
  rho_cutoff = 3.d6

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('--rho_cutoff')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) rho_cutoff

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0) then
     print *, "usage: fwdconvect args"
     print *, "args [-p|--pltfile]   plotfile   : plot file directory (required)"
     print *, "     --rho_cutoff cutoff value   : low density cutoff (optional)"
     stop
  end if

  print *, 'pltfile    = "', trim(pltfile), '"'
  print *, 'rho_cutoff = "', rho_cutoff

  call build(pf, pltfile, unit)

  time = pf%tm

!  do i = 1, pf%flevel
!     call fab_bind_level(pf, i)
!  end do


  ! figure out the variable indices from the plotfile
  
  ! density
  dens_comp = plotfile_var_index(pf, "density")
  
  ! temperature (here we used T from p0)
  temp_comp = plotfile_var_index(pf, "tfromp") 

  ! species
  XC12_comp = plotfile_var_index(pf, "X(C12)")
  XO16_comp = plotfile_var_index(pf, "X(O16)")
  XMg24_comp = plotfile_var_index(pf, "X(Mg24)")


  if ( dens_comp < 0 .or. temp_comp < 0 .or. &
       XC12_comp < 0 .or. XO16_comp < 0 .or. XMg24_comp < 0) then
     call bl_error("ERROR: varaible(s) not defined")
  endif



  ! get the indices that the network uses to address the species
  call network_init()
  call eos_init(use_eos_coulomb=.true.)



  ic12  = network_species_index("carbon-12")
  io16  = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  print *, ic12, io16, img24
  

  ! get dx for the coarse level.  
  dx = plotfile_get_dx(pf, 1)


  ! get the index bounds for the finest level.  Note, lo and hi are
  ! ZERO based indicies
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  print *, 'Size of domain (zones): ', fhi(1)-flo(1)+1, fhi(2)-flo(2)+1, fhi(3)-flo(3)+1



  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)))


  !----------------------------------------------------------------------------
  ! compute the total internal energy
  !----------------------------------------------------------------------------

  internal_energy = ZERO

  imask(:,:,:) = .true.


  ! loop over the data, starting at the finest grid, and if we haven't
  ! already stored data in that grid location (according to imask),
  ! store it.  


  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

     print *, 'processing level ', i, ' rr = ', rr

     do j = 1, nboxes(pf, i)
        
        ! read in the data 1 patch at a time, and only read in
        ! those components that we are going to use.
        call fab_bind_comp_vec(pf, i, j, &
                               (/ dens_comp,    &  ! density       = 1
                                  temp_comp,    &  ! temperature   = 2
                                  XC12_comp,    &  ! X(C12)        = 3
                                  XO16_comp,    &  ! X(O16)        = 4
                                  XMg24_comp /))   ! X(Mg24)       = 5

        lo(:) = 1
        hi(:) = 1
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        ! get a pointer to the current patch
        p => dataptr(pf, i, j)
        
        ! loop over all of the zones in the patch.  Here, we convert
        ! the cell-centered indices at the current level into the
        ! corresponding RANGE on the finest level, and test if we've
        ! stored data in any of those locations.  If we haven't then
        ! we store this level's data and mark that range as filled.
        do kk = lo(3), hi(3)
           zz = (kk + HALF)*dx(3)/rr

           do jj = lo(2), hi(2)
              yy = (jj + HALF)*dx(2)/rr

              do ii = lo(1), hi(1)
                 xx = (ii + HALF)*dx(1)/rr

                 if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1, &
                                kk*r1:(kk+1)*r1-1) ) ) then


                    ! weight the zone's data by its size

                    ! note, for p(:,:,:,n), n refers to the order
                    ! used in the call to fab_bind_comp_vec


                    ! call the EOS
                    temp_eos(1) = p(ii,jj,kk,2)
                    den_eos(1)  = p(ii,jj,kk,1)
                    xn_eos(1,ic12) = p(ii,jj,kk,3)
                    xn_eos(1,io16) = p(ii,jj,kk,4)
                    xn_eos(1,img24) = p(ii,jj,kk,5)

                    call eos(eos_input_rt, den_eos, temp_eos, &
                             npts, &
                             xn_eos, &
                             p_eos, h_eos, e_eos, &
                             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                             dpdX_eos, dhdX_eos, &
                             gam1_eos, cs_eos, s_eos, &
                             dsdt_eos, dsdr_eos, &
                             do_diag)


                    ! internal energy is sum { rho * e * dV } 
                    if (p(ii,jj,kk,1) >= rho_cutoff) then
                       internal_energy = internal_energy + &
                            p(ii,jj,kk,1) * e_eos(1) * &
                            (dx(1)/rr)*(dx(2)/rr)*(dx(3)/rr)
                    endif


                    imask(ii*r1:(ii+1)*r1-1, &
                          jj*r1:(jj+1)*r1-1, &
                          kk*r1:(kk+1)*r1-1) = .false.
                 end if

              end do
           enddo
        enddo

        call fab_unbind(pf, i, j)
                
     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do


  print *, "internal energy = ", internal_energy


  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

end program feint
