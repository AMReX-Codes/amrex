! Analysis routines for the 3-d WD convect MAESTRO problem 
! These are based on fsedov3d_sph.f90
!
! Here we compute angle-averaged quantities, <q>, and RMS quantities,
! q' = sqrt { < (q - <q>)**2 > }, where the averaging is done at constant
! radius.
!
! For density, rho_0 = <rho>, and in the plotfiles we store rhopert =
! rho - rho0, so we can compute the RMS density directly from this.
!
! Similarly, for temperature, we store tpert = T - <T>, so we will use
! tpert directly from the plotfiles for computing the RMS value.
!
! we also compute some global quantities, including T_peak and its
! location, enuc_peak and its location, and the total kinetic energy
!
! for the total kinetic energy inside the star, we only consider those
! zones where rho > rho_cutoff
!
! --globals_only will just print T_max, it's location, and kinetic energy
!
program fwdconvect

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj, kk
  real(kind=dp_t) :: xx, yy, zz
  integer :: rr, r1
  integer :: uno

  integer :: nbins
  real(kind=dp_t), allocatable :: r(:)
  real(kind=dp_t) :: maxdist, x_maxdist, y_maxdist, z_maxdist
  real(kind=dp_t) :: xctr, yctr, zctr

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  real(kind=dp_t) :: r_zone
  integer :: index

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer, allocatable :: ncount(:)
  real(kind=dp_t), allocatable :: dens_avg_bin(:), dens_rms_bin(:)
  real(kind=dp_t), allocatable :: temp_avg_bin(:), temp_rms_bin(:)
  real(kind=dp_t), allocatable :: XC12_avg_bin(:), XO16_avg_bin(:), Xash_avg_bin(:)
  real(kind=dp_t), allocatable :: entropy_avg_bin(:), entropy_rms_bin(:)
  real(kind=dp_t) :: kinetic_energy

  integer :: dens_comp, temp_comp, XC12_comp, XO16_comp, Xash_comp
  integer :: rhopert_comp, tpert_comp
  integer :: magvel_comp, enuc_comp, xvel_comp, yvel_comp, zvel_comp
  integer :: s_comp, spert_comp

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  character(len=256) :: slicefile
  character(len=256) :: pltfile
  real(kind=dp_t) :: rho_cutoff

  integer :: narg, farg
  character(len=256) :: fname

  real(kind=dp_t) :: time

  real(kind=dp_t) :: T_peak
  real(kind=dp_t) :: xloc_Tpeak, yloc_Tpeak, zloc_Tpeak, R_Tpeak
  real(kind=dp_t) :: vx_Tpeak, vy_Tpeak, vz_Tpeak, vr_Tpeak

  real(kind=dp_t) :: enuc_peak
  real(kind=dp_t) :: xloc_enucpeak, yloc_enucpeak, zloc_enucpeak, R_enucpeak

  logical :: globals_only

  unit = unit_new()
  uno =  unit_new()


  ! set the defaults
  slicefile = ''
  pltfile  = ''
  rho_cutoff = 3.d6
  globals_only = .false.

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-s', '--slicefile')
        farg = farg + 1
        call get_command_argument(farg, value = slicefile)

     case ('--rho_cutoff')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) rho_cutoff

     case ('--globals_only')
        globals_only = .true.

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 .OR. (len_trim(slicefile) == 0 .and. .not. globals_only)) then
     print *, "usage: fwdconvect args"
     print *, "args [-p|--pltfile]   plotfile   : plot file directory (required)"
     print *, "     [-s|--slicefile] slice file : slice file (required if not globals_only)"
     print *, "     --rho_cutoff cutoff value   : low density cutoff (optional)"
     print *, "     --globals_only              : only output global quantities (optional) "
     stop
  end if

  if (.not. globals_only) then
     print *, 'pltfile    = "', trim(pltfile), '"'
     print *, 'slicefile  = "', trim(slicefile), '"'
     print *, 'rho_cutoff = "', rho_cutoff
  endif

  call build(pf, pltfile, unit)

  time = pf%tm

!  do i = 1, pf%flevel
!     call fab_bind_level(pf, i)
!  end do


  ! figure out the variable indices
  
  ! density
  dens_comp = plotfile_var_index(pf, "density")
  
  ! temperature (here we used T from p0)
  temp_comp = plotfile_var_index(pf, "tfromp") 

  ! X(C12)
  XC12_comp = plotfile_var_index(pf, "X(C12)")

  ! X(O16)
  XO16_comp = plotfile_var_index(pf, "X(O16)")

  ! Xash -- look for ether X(Mg24) or X(ash)
  Xash_comp = plotfile_var_index(pf, "X(Mg24)")
  if (Xash_comp < 1) then
     Xash_comp = plotfile_var_index(pf, "X(ash)")
  endif

  ! rhopert
  rhopert_comp = plotfile_var_index(pf, "rhopert")

  ! tpert
  tpert_comp = plotfile_var_index(pf, "tpert") 

  ! magvel
  magvel_comp = plotfile_var_index(pf, "magvel") 

  ! enuc
  enuc_comp = plotfile_var_index(pf, "enucdot") 

  ! xvel
  xvel_comp = plotfile_var_index(pf, "x_vel")

  ! yvel
  yvel_comp = plotfile_var_index(pf, "y_vel") 

  ! zvel
  zvel_comp = plotfile_var_index(pf, "z_vel") 

  ! entropy
  s_comp = plotfile_var_index(pf, "entropy")
  
  ! entropy perturbation
  spert_comp = plotfile_var_index(pf, "entropypert")



  if ( dens_comp < 0 .or. temp_comp < 0 .or. &
       XC12_comp < 0 .or. XO16_comp < 0 .or. Xash_comp < 0 .or. &
       rhopert_comp < 0 .or. tpert_comp < 0 .or. &
       magvel_comp < 0 .or. enuc_comp < 0 .or. &
       xvel_comp < 0 .or. yvel_comp < 0 .or. zvel_comp < 0 .or. &
       s_comp < 0 .or. spert_comp < 0) then
     call bl_error("ERROR: varaible(s) not defined")
  endif


  ! get dx for the coarse level.  
  dx = plotfile_get_dx(pf, 1)


  ! get the index bounds for the finest level.  Note, lo and hi are
  ! ZERO based indicies
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  if (.not. globals_only) then
     print *, 'Size of domain (zones): ', fhi(1)-flo(1)+1, fhi(2)-flo(2)+1, fhi(3)-flo(3)+1
  endif

  ! the default for the center of the star will be the geometric center
  ! of the domain
  xctr = HALF*(pf%phi(1) + pf%plo(1))
  yctr = HALF*(pf%phi(2) + pf%plo(2))
  zctr = HALF*(pf%phi(3) + pf%plo(3))

  if (.not. globals_only) then
     print *, 'Center of the star: ', xctr, yctr, zctr
  endif

  ! compute the size of the radially-binned array -- we'll do it to
  ! the furtherest corner of the domain
  x_maxdist = max(abs(pf%phi(1) - xctr), abs(pf%plo(1) - xctr))
  y_maxdist = max(abs(pf%phi(2) - yctr), abs(pf%plo(2) - yctr))
  z_maxdist = max(abs(pf%phi(3) - zctr), abs(pf%plo(3) - zctr))
  
  maxdist = sqrt(x_maxdist**2 + y_maxdist**2 + z_maxdist**2)

  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))
  nbins = int(maxdist/dx_fine)

  allocate(r(0:nbins-1))

  do i = 0, nbins-1
     r(i) = (dble(i) + HALF)*dx_fine
  enddo


  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)))

  ! ncount will keep track of how many fine zones were written into
  ! each bin
  allocate(   ncount(0:nbins-1))


  !----------------------------------------------------------------------------
  ! compute the angle averaged quantities
  !----------------------------------------------------------------------------

  ! allocate storage for the data 
  allocate(   dens_avg_bin(0:nbins-1))
  allocate(   temp_avg_bin(0:nbins-1))
  allocate(   XC12_avg_bin(0:nbins-1))
  allocate(   XO16_avg_bin(0:nbins-1))
  allocate(   Xash_avg_bin(0:nbins-1))
  allocate(   dens_rms_bin(0:nbins-1))
  allocate(   temp_rms_bin(0:nbins-1))
  allocate(entropy_avg_bin(0:nbins-1))
  allocate(entropy_rms_bin(0:nbins-1))

  ncount(:) = 0
  dens_avg_bin(:) = ZERO
  temp_avg_bin(:) = ZERO
  XC12_avg_bin(:) = ZERO
  XO16_avg_bin(:) = ZERO
  Xash_avg_bin(:) = ZERO
  dens_rms_bin(:) = ZERO
  temp_rms_bin(:) = ZERO
  entropy_avg_bin(:) = ZERO
  entropy_rms_bin(:) = ZERO

  kinetic_energy = ZERO

  imask(:,:,:) = .true.

  T_peak = 0.0
  xloc_Tpeak = 0.0
  yloc_Tpeak = 0.0
  zloc_Tpeak = 0.0

  vx_Tpeak = 0.0
  vy_Tpeak = 0.0
  vz_Tpeak = 0.0
  vr_Tpeak = 0.0

  enuc_peak = 0.0
  xloc_enucpeak = 0.0
  yloc_enucpeak = 0.0
  zloc_enucpeak = 0.0

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

     if (.not. globals_only) then
        print *, 'processing level ', i, ' rr = ', rr
     endif

     do j = 1, nboxes(pf, i)
        
        ! read in the data 1 patch at a time -- read in all the variables
        call fab_bind(pf, i, j)

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

                    r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)

                    index = r_zone/dx_fine

                    ! weight the zone's data by its size

                    ! note, for p(:,:,:,n), n refers to index of the
                    ! variable as found via plotfile_var_index

                    dens_avg_bin(index) = dens_avg_bin(index) + &
                         p(ii,jj,kk,dens_comp) * r1**3

                    temp_avg_bin(index) = temp_avg_bin(index) + &
                         p(ii,jj,kk,temp_comp) * r1**3

                    entropy_avg_bin(index) = entropy_avg_bin(index) + &
                         p(ii,jj,kk,s_comp) * r1**3

                    ! do the Favre-average here, < rho * X(C12) > / < rho >
                    XC12_avg_bin(index) = XC12_avg_bin(index) + &
                         p(ii,jj,kk,dens_comp)*p(ii,jj,kk,XC12_comp) * r1**3

                    ! do the Favre-average here, < rho * X(O16) > / < rho >
                    XO16_avg_bin(index) = XO16_avg_bin(index) + &
                         p(ii,jj,kk,dens_comp)*p(ii,jj,kk,XO16_comp) * r1**3

                    ! do the Favre-average here, < rho * X(ash) > / < rho >
                    Xash_avg_bin(index) = Xash_avg_bin(index) + &
                         p(ii,jj,kk,dens_comp)*p(ii,jj,kk,Xash_comp) * r1**3

                    ! for the RMS quantities, we use the perturbational quantities
                    ! already stored in the plotfile
                    dens_rms_bin(index) = dens_rms_bin(index) + &
                         p(ii,jj,kk,rhopert_comp)*p(ii,jj,kk,rhopert_comp) * r1**3

                    temp_rms_bin(index) = temp_rms_bin(index) + &
                         p(ii,jj,kk,tpert_comp)*p(ii,jj,kk,tpert_comp) * r1**3

                    entropy_rms_bin(index) = entropy_rms_bin(index) + &
                         p(ii,jj,kk,spert_comp)*p(ii,jj,kk,spert_comp) * r1**3

                    ! kinetic energy is integral of rho U U dV
                    if (p(ii,jj,kk,1) >= rho_cutoff) then
                       kinetic_energy = kinetic_energy + &
                            p(ii,jj,kk,dens_comp) * p(ii,jj,kk,magvel_comp)**2 * &
                            (dx(1)/rr)*(dx(2)/rr)*(dx(3)/rr)
                    endif

                    ncount(index) = ncount(index) + r1**3


                    ! store the location and value of the peak temperature
                    if (p(ii,jj,kk,temp_comp) > T_peak) then
                       T_peak = p(ii,jj,kk,temp_comp)
                       xloc_Tpeak = xx
                       yloc_Tpeak = yy
                       zloc_Tpeak = zz

                       R_Tpeak = sqrt( (xctr - xloc_Tpeak)**2 + &
                                       (yctr - yloc_Tpeak)**2 + &
                                       (zctr - zloc_Tpeak)**2 )

                       vx_Tpeak = p(ii,jj,kk,xvel_comp)
                       vy_Tpeak = p(ii,jj,kk,yvel_comp)
                       vz_Tpeak = p(ii,jj,kk,zvel_comp)

                       vr_Tpeak = vx_Tpeak*(xx - xctr)/R_Tpeak + &
                                  vy_Tpeak*(yy - yctr)/R_Tpeak + &
                                  vz_Tpeak*(zz - zctr)/R_Tpeak
                       

                    endif


                    ! store the location and value of the peak enucdot
                    if (p(ii,jj,kk,enuc_comp) > enuc_peak) then
                       enuc_peak = p(ii,jj,kk,enuc_comp)

                       xloc_enucpeak = xx
                       yloc_enucpeak = yy
                       zloc_enucpeak = zz

                       R_enucpeak = sqrt( (xctr - xloc_enucpeak)**2 + &
                                          (yctr - yloc_enucpeak)**2 + &
                                          (zctr - zloc_enucpeak)**2 )
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


  ! normalize
  do i = 0, nbins-1
     if (ncount(i) /= 0) then

        ! simple averages
        dens_avg_bin(i) = dens_avg_bin(i)/ncount(i)
        temp_avg_bin(i) = temp_avg_bin(i)/ncount(i)
        entropy_avg_bin(i) = entropy_avg_bin(i)/ncount(i)

        ! Favre averaged composition
        XC12_avg_bin(i) = (XC12_avg_bin(i)/ncount(i)) / dens_avg_bin(i)
        XO16_avg_bin(i) = (XO16_avg_bin(i)/ncount(i)) / dens_avg_bin(i)
        Xash_avg_bin(i) = (Xash_avg_bin(i)/ncount(i)) / dens_avg_bin(i)

        ! RMS quantities
        dens_rms_bin(i) = sqrt(dens_rms_bin(i)/ncount(i))        
        temp_rms_bin(i) = sqrt(temp_rms_bin(i)/ncount(i))
        entropy_rms_bin(i) = sqrt(entropy_rms_bin(i)/ncount(i))

        ! we already normalized the kinetic energy by multiplying 
        ! each zone by dV

     endif
  enddo



990  format("# time = ",g24.12)
991  format("# ---------------------------------------------------------------------------")
992  format("# total kinetic energy = ", g24.12)
994  format("# peak temperature = ", g24.12)
995  format("# peak temp loc (x,y,z) = ", 3(g24.12,1x))
996  format("# peak temp radius = ", g24.12)
9961 format("# velocity @ peak T loc (vx, vy, vz) = ", 3(g24.12,1x))
9962 format("# radial velocity @ peak T loc = ", g24.12,1x)
997  format("# peak enucdot = ", g24.12)
998  format("# peak enucdot loc (x,y,z) = ", 3(g24.12,1x))
999  format("# peak enucdot radius = ", g24.12)
1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  ! slicefile
  if (.not. globals_only) then
     open(unit=uno, file=slicefile, status = 'replace')

     ! write the header
     write(uno,990) time
     write(uno,991)

     write(uno,992) kinetic_energy
     write(uno,991)

     write(uno,994) T_peak
     write(uno,995) xloc_Tpeak, yloc_Tpeak, zloc_Tpeak     
     write(uno,996) R_Tpeak
     write(uno,9961) vx_Tpeak, vy_Tpeak, vz_Tpeak
     write(uno,9962) vr_Tpeak

     write(uno,997) enuc_peak
     write(uno,998) xloc_enucpeak, yloc_enucpeak, zloc_enucpeak
     write(uno,999) R_enucpeak
     write(uno,991)

     write(uno,1000) "r", "density", "temperature", "entropy", "X(C12)", "X(O16)", "X(ash)", &
          "RMS density", "RMS temperature", "RMS entropy"

     ! write the data in columns
     do i = 0, nbins-1
        ! Use this to protect against a number being xx.e-100
        !   which will print without the "e"
        if (abs(   dens_avg_bin(i)) .lt. 1.d-99)    dens_avg_bin(i) = 0.d0
        if (abs(   temp_avg_bin(i)) .lt. 1.d-99)    temp_avg_bin(i) = 0.d0
        if (abs(entropy_avg_bin(i)) .lt. 1.d-99) entropy_avg_bin(i) = 0.d0
        if (abs(   XC12_avg_bin(i)) .lt. 1.d-99)    XC12_avg_bin(i) = 0.d0
        if (abs(   XO16_avg_bin(i)) .lt. 1.d-99)    XO16_avg_bin(i) = 0.d0
        if (abs(   Xash_avg_bin(i)) .lt. 1.d-99)    Xash_avg_bin(i) = 0.d0
        if (abs(   dens_rms_bin(i)) .lt. 1.d-99)    dens_rms_bin(i) = 0.d0
        if (abs(   temp_rms_bin(i)) .lt. 1.d-99)    temp_rms_bin(i) = 0.d0
        if (abs(entropy_rms_bin(i)) .lt. 1.d-99) entropy_rms_bin(i) = 0.d0
        
        write(uno,1001) r(i), dens_avg_bin(i), temp_avg_bin(i), &
             entropy_avg_bin(i), &
             XC12_avg_bin(i), XO16_avg_bin(i), Xash_avg_bin(i), &
             dens_rms_bin(i), temp_rms_bin(i), &
             entropy_rms_bin(i)
     end do

     close(unit=uno)
  endif

  if (.not. globals_only) then
     print *, 'Total kinetic energy = ', kinetic_energy
     print *, 'Peak temperature = ', T_peak
     print *, 'Peak temperature location (x,y,z) = ', xloc_Tpeak, yloc_Tpeak, zloc_Tpeak
     print *, 'Peak temperature radius = ', R_Tpeak
     print *, 'Peak enucdot = ', enuc_peak
     print *, 'Peak enucdot location (x,y,z) = ', xloc_enucpeak, yloc_enucpeak, zloc_enucpeak
     print *, 'Peak enucdot radius = ', R_enucpeak
  else
     write (*,1000) "time", "T_peak", "x(T_peak)", "y(T_peak)", "z(T_peak)", "R(T_peak)", &
          "enuc_peak", "x(enuc_peak)", "y(enuc_peak)", "z(enuc_peak)", "R(enuc_peak)", "KE"
     write (*,1001) time, T_peak, xloc_Tpeak, yloc_Tpeak, zloc_Tpeak, R_Tpeak, &
          enuc_peak, xloc_enucpeak, yloc_enucpeak, zloc_enucpeak, R_enucpeak, kinetic_energy
  endif

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

end program fwdconvect
