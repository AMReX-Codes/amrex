
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      implicit none
      integer init, namlen
      integer name(namlen)
      real(rt) :: problo(3), probhi(3)

      integer untin,i

      namelist /fortin/ probtype, p_ambient, dens_ambient, exp_energy, &
                        r_init, nsub,  &
                        denerr,dengrad,max_denerr_lev,max_dengrad_lev, &
                        presserr,pressgrad,max_presserr_lev,max_pressgrad_lev

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
         
!    set namelist defaults

      p_ambient = 1.d-5        ! ambient pressure (in erg/cc)
      dens_ambient = 1.d0      ! ambient density (in g/cc)
      exp_energy = 1.d0        ! absolute energy of the explosion (in erg)
      r_init = 0.05d0          ! initial radius of the explosion (in cm)
      nsub = 4

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = -1
      max_dengrad_lev = -1

      presserr = 1.d20
      pressgrad = 1.d20
      max_presserr_lev = -1
      max_pressgrad_lev = -1

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

! set local variable defaults
      center(1) = (problo(1)+probhi(1))/2.d0
      center(2) = (problo(2)+probhi(2))/2.d0
      center(3) = (problo(3)+probhi(3))/2.d0

      end

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
      use bl_constants_module, only: M_PI, FOUR3RD
      use atomic_rates_module, only: XHYDROGEN
      use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                     gamma_minus_1

      implicit none
      integer level, ns, nd
      integer lo(3), hi(3)
      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt) :: xlo(3), xhi(3), time, delta(3)
      real(rt) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      real(rt) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)

      real(rt) :: xmin,ymin,zmin
      real(rt) :: xx, yy, zz
      real(rt) :: dist
      real(rt) :: eint, p_zone
      real(rt) :: vctr, p_exp

      integer i,j,k, ii, jj, kk
      integer :: npert, nambient

      if (probtype .eq. 32) then

!       set explosion pressure -- we will convert the point-explosion energy into
!       a corresponding pressure distributed throughout the perturbed volume
        vctr  = M_PI*r_init**2
        p_exp = gamma_minus_1*exp_energy/vctr

        do k = lo(3), hi(3)
         zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

         do j = lo(2), hi(2)
            ymin = xlo(2) + delta(2)*dble(j-lo(2))

            do i = lo(1), hi(1)
               xmin = xlo(1) + delta(1)*dble(i-lo(1))

               npert = 0
               nambient = 0

               do jj = 0, nsub-1
                  yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)

                  do ii = 0, nsub-1
                     xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)

                     dist = (center(1)-xx)**2 + (center(2)-yy)**2

                     if(dist <= r_init**2) then
                        npert = npert + 1
                     else
                        nambient = nambient + 1
                     endif

                  enddo
               enddo

               p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient) / &
                        (dble(npert) + dble(nambient))

               eint = p_zone/gamma_minus_1

               state(i,j,k,URHO) = dens_ambient
               state(i,j,k,UMX) = 0.d0
               state(i,j,k,UMY) = 0.d0
               state(i,j,k,UMZ) = 0.d0

               state(i,j,k,UEDEN) = eint +  &
                    0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                           state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                           state(i,j,k,UMZ)**2/state(i,j,k,URHO))

               state(i,j,k,UEINT) = eint

            enddo
         enddo
        enddo

      else if (probtype .eq. 33) then

!       set explosion pressure -- we will convert the point-explosion energy into
!       a corresponding pressure distributed throughout the perturbed volume
        vctr  = FOUR3RD*M_PI*r_init**3
        p_exp = gamma_minus_1*exp_energy/vctr

        do k = lo(3), hi(3)
         zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

         do j = lo(2), hi(2)
            ymin = xlo(2) + delta(2)*dble(j-lo(2))

            do i = lo(1), hi(1)
               xmin = xlo(1) + delta(1)*dble(i-lo(1))

               npert = 0
               nambient = 0

               do kk = 0, nsub-1
                  zz = zmin + (delta(3)/dble(nsub))*(kk + 0.5d0)

                  do jj = 0, nsub-1
                     yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)

                     do ii = 0, nsub-1
                        xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)

                        dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2

                        if(dist <= r_init**2) then
                           npert = npert + 1
                        else
                           nambient = nambient + 1
                        endif

                     enddo
                  enddo
               enddo

               p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient)/  &
                    dble(nsub*nsub*nsub)

               eint = p_zone/gamma_minus_1

               state(i,j,k,URHO) = dens_ambient
               state(i,j,k,UMX) = 0.d0
               state(i,j,k,UMY) = 0.d0
               state(i,j,k,UMZ) = 0.d0

               state(i,j,k,UEDEN) = eint + &
                    0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                           state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                           state(i,j,k,UMZ)**2/state(i,j,k,URHO))

               state(i,j,k,UEINT) = eint
 
               if (UFS .gt. -1) then
                   state(i,j,k,UFS  ) =         XHYDROGEN  * state(i,j,k,URHO)
                   state(i,j,k,UFS+1) = (1.d0 - XHYDROGEN) * state(i,j,k,URHO)
               end if

            enddo
         enddo
        enddo

      else 

         call bl_error('Dont know this probtype in initdata')

      end if

      end subroutine fort_initdata
