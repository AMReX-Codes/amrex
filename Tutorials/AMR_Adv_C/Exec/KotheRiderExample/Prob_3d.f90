
      subroutine PROBINIT (init,name,namlen,problo,probhi)

      use probdata_module

      implicit none

      integer init, namlen
      integer name(namlen)
      double precision problo(2), probhi(2)

      integer untin,i

      namelist /fortin/ probtype,adverr,specerr,specgrad,max_specerr_lev,max_specgrad_lev,&
                        diff_coeff

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
         
!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

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
! ::: dx        => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
      subroutine initdata(level,time,lo,hi,nscal, &
                          state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                          dx,xlo,xhi)

      use probdata_module
      use meth_params_module , only: NVAR, URHO, UX, UY, UZ, UFS, UFA
      use network            , only: nspec
      implicit none

      integer level, nscal
      integer lo(3), hi(3)
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      double precision xlo(3), xhi(3), time, dx(3)
      double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

      integer          :: i,j,k,idir
      double precision :: x,y,z,r1,r2

      idir = 3

      state(:,:,:,UFA  ) = 0.d0
      state(:,:,:,UFS  ) = 1.d0
      state(:,:,:,UFS+1) = 0.d0
      state(:,:,:,UFS+2:UFS+nspec-1) = 0.d0

      state(:,:,:,URHO)  = 1.d0
      state(:,:,:,UX:UZ) = 1.d0

      do k = lo(3), hi(3)
        z = xlo(3) + (k-lo(2)+0.5d0)*dx(3) 

        do j = lo(2), hi(2)
           y = xlo(2) + (j-lo(2)+0.5d0)*dx(2) 

           do i = lo(1), hi(1)
             x = xlo(1) + (i-lo(1)+0.5d0)*dx(1) 

             ! Define two blobs 
             if (idir.eq.3) then 
                r1 = sqrt( (x-0.25d0)**2+(y-0.25d0)**2 )
             else if (idir.eq.2) then
                r1 = sqrt( (x-0.25d0)**2+(z-0.25d0)**2 )
             else if (idir.eq.1) then
                r1 = sqrt( (y-0.25d0)**2+(z-0.25d0)**2 )
             end if

             if (r1.lt.0.2d0) then
                state(i,j,k,UFS  ) = 1.d0 - 0.5d0 * (1.d0-tanh(30.*(r1-0.1)))
                state(i,j,k,UFS+1) = 1.d0 - state(i,j,k,UFS)
                state(i,j,k,UFA  ) = 1.d0
             end if

             if (idir.eq.3) then 
                r2 = sqrt( (x-0.75d0)**2+(y-0.75d0)**2 )
             else if (idir.eq.2) then
                r2 = sqrt( (x-0.75d0)**2+(z-0.75d0)**2 )
             else if (idir.eq.1) then
                r2 = sqrt( (y-0.75d0)**2+(z-0.75d0)**2 )
             end if

             if (r2.lt.0.2d0) then
                state(i,j,k,UFS  ) = 1.d0 - 0.5d0 * (1.d0-tanh(30.*(r2-0.1))) 
                state(i,j,k,UFS+1) = 1.d0 - state(i,j,k,UFS)
                state(i,j,k,UFA  ) = 2.d0
             end if

             state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO)*state(i,j,k,UFS:UFS+nspec-1)
             state(i,j,k,UFA            ) = state(i,j,k,URHO)*state(i,j,k,UFA              ) 

         enddo
      enddo
      enddo

      end subroutine initdata

