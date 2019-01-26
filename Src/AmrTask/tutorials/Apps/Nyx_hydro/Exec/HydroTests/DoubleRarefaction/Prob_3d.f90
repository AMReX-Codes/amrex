
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      use meth_params_module, only : gamma_minus_1
      use network   , only : network_init
      implicit none

      integer init, namlen
      integer name(namlen)
      real(rt) problo(3), probhi(3)

      real(rt) vctr
      integer untin,i

      namelist /fortin/ &
           p_l, u_l, rho_l, p_r, u_r, rho_r, frac, idir, &
           denerr,  dengrad,  max_denerr_lev,  max_dengrad_lev, &
           velgrad,  max_velgrad_lev, &
           presserr,pressgrad,max_presserr_lev,max_pressgrad_lev

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      call network_init()

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

! set namelist defaults

      p_l = 1.0               ! left pressure (erg/cc)
      u_l = 0.0               ! left velocity (cm/s)
      rho_l = 1.0             ! left density (g/cc)

      p_r = 0.1               ! right pressure (erg/cc)
      u_r = 0.0               ! right velocity (cm/s)
      rho_r = 0.125           ! right density (g/cc)

      idir = 1                ! direction across which to jump
      frac = 0.5              ! fraction of the domain for the interface

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = -1
      max_dengrad_lev = -1

      presserr = 1.d20
      pressgrad = 1.d20
      max_presserr_lev = -1
      max_pressgrad_lev = -1

      velgrad = 1.d20
      max_velgrad_lev = -1

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      center(1) = frac*(problo(1)+probhi(1))
      center(2) = frac*(problo(2)+probhi(2))
      center(3) = frac*(problo(3)+probhi(3))

! compute the internal energy (erg/cc) for the left and right state
      rhoe_l = p_l/gamma_minus_1
      rhoe_r = p_r/gamma_minus_1

      end subroutine amrex_probinit

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

      subroutine fort_initdata(level,time,lo,hi, &
                               ns, state   ,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                               nd, diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                               delta,xlo,xhi)  &
                               bind(C, name="fort_initdata")

     use amrex_fort_module, only : rt => amrex_real
     use probdata_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
     use atomic_rates_module, only: XHYDROGEN
     implicit none

     integer level, ns, nd
     integer lo(3), hi(3)
     integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
     integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
     real(rt) xlo(3), xhi(3), time, delta(3)
     real(rt)    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
     real(rt) diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)

     real(rt) xcen,ycen,zcen
     integer i,j,k

      do k = lo(3), hi(3)
         zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)

         do j = lo(2), hi(2)
            ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)

            do i = lo(1), hi(1)
               xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

               if (idir == 1) then
                  if (xcen <= center(1)) then
                     state(i,j,k,URHO) = rho_l
                     state(i,j,k,UMX) = rho_l*u_l
                     state(i,j,k,UMY) = 0.d0
                     state(i,j,k,UMZ) = 0.d0
                     state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                     state(i,j,k,UEINT) = rhoe_l
                  else
                     state(i,j,k,URHO) = rho_r
                     state(i,j,k,UMX) = rho_r*u_r
                     state(i,j,k,UMY) = 0.d0
                     state(i,j,k,UMZ) = 0.d0
                     state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                     state(i,j,k,UEINT) = rhoe_r
                  endif

               else if (idir == 2) then
                  if (ycen <= center(2)) then
                     state(i,j,k,URHO) = rho_l
                     state(i,j,k,UMX) = 0.d0
                     state(i,j,k,UMY) = rho_l*u_l
                     state(i,j,k,UMZ) = 0.d0
                     state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                     state(i,j,k,UEINT) = rhoe_l
                  else
                     state(i,j,k,URHO) = rho_r
                     state(i,j,k,UMX) = 0.d0
                     state(i,j,k,UMY) = rho_r*u_r
                     state(i,j,k,UMZ) = 0.d0
                     state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                     state(i,j,k,UEINT) = rhoe_r
                  endif

               else if (idir == 3) then
                  if (zcen <= center(3)) then
                     state(i,j,k,URHO) = rho_l
                     state(i,j,k,UMX) = 0.d0
                     state(i,j,k,UMY) = 0.d0
                     state(i,j,k,UMZ) = rho_l*u_l
                     state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                     state(i,j,k,UEINT) = rhoe_l
                  else
                     state(i,j,k,URHO) = rho_r
                     state(i,j,k,UMX) = 0.d0
                     state(i,j,k,UMY) = 0.d0
                     state(i,j,k,UMZ) = rho_r*u_r
                     state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                     state(i,j,k,UEINT) = rhoe_r
                  endif

               else
                  call bl_abort('invalid idir')
               endif
 
               if (UFS .gt. -1) then
                   state(i,j,k,UFS  ) =         XHYDROGEN  * state(i,j,k,URHO)
                   state(i,j,k,UFS+1) = (1.d0 - XHYDROGEN) * state(i,j,k,URHO)
               end if

            enddo
         enddo
      enddo

      end subroutine fort_initdata
