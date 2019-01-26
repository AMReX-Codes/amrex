
     subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      use turbinit_module

      implicit none
      integer init, namlen
      integer name(namlen)
      real(rt) problo(3), probhi(3)
 
      integer untin,i

      namelist /fortin/ &
           denerr, dengrad,  max_denerr_lev,max_dengrad_lev, &
           alpha, rho0, temp0, &
           turb_scale, force_scale, forcing_type, spectrum_type, &
           mode_start, nmodes, forcing_time_scale_min, forcing_time_scale_max, &
           stop_forcing, do_mode_division
!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

!     Set namelist defaults

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = 10
      max_dengrad_lev = 10

      alpha = 0.d0        ! adiabatic
      rho0 = 1.d0
      temp0 = 10.d0

      stop_forcing = 1.d9 ! never stops
      force_scale = 0.0   ! force is zero

!     Read namelists

      untin = 9 
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      ! These are stored in probdata_module
      prob_lo(:) = problo(:)
      prob_hi(:) = probhi(:)

      call turbforce_init(problo,probhi)

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

      use probdata_module
      use turbforce_module
      use bl_constants_module, only : TWO, ONE, HALF, ZERO, M_PI
      use atomic_rates_module, only: XHYDROGEN
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, TEMP_COMP, NE_COMP, &
	                             small_dens, small_temp, small_pres
      use eos_module, only: nyx_eos_given_RT
      use eos_params_module
 
      implicit none

      integer level, ns, nd
      integer lo(3), hi(3)
      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt) xlo(3), xhi(3), time, delta(3)
      real(rt)    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      real(rt) diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)

      integer          :: i,j,k
      real(rt) :: fact,twicePi
      real(rt) :: eint0,rhoe0,ne0,pres0,a

      a=1.d0
      ne0=1.d0
      twicePi  = TWO*M_PI

      if (turb_scale.gt.ZERO) then
         fact=turb_scale
      else
         fact=ONE
      end if

      ! Call EOS to get the internal energy for constant initial temperature
      call nyx_eos_given_RT(eint0, pres0, rho0, temp0, ne0, a)

      rhoe0 = rho0*eint0

      small_dens = 1.d-9*rhoe0
      small_temp = 1.d-3*temp0
      small_pres = 1.d-12*pres0
 
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               diag_eos(i,j,k,TEMP_COMP) = temp0
               diag_eos(i,j,k,  NE_COMP) = ne0

               state(i,j,k,URHO) = rho0
               state(i,j,k,UMX:UMZ) = 0.0d0
               state(i,j,k,UEINT) = rhoe0

               state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                  0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

               if (UFS .gt. -1) then
                  state(i,j,k,UFS  ) = XHYDROGEN
                  state(i,j,k,UFS+1) = (1.d0 - XHYDROGEN)
               end if
            enddo
         enddo
      enddo

      end subroutine fort_initdata
