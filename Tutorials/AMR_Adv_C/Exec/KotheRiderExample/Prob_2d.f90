
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
                          state,state_l1,state_l2,state_h1,state_h2, &
                          dx,xlo,xhi)

      use probdata_module
      use meth_params_module , only: NVAR, URHO, UX, UY, UFS, UFA
      use network            , only: nspec
      implicit none

      integer level, nscal
      integer lo(2), hi(2)
      integer state_l1,state_l2,state_h1,state_h2
      double precision xlo(2), xhi(2), time, dx(2)
      double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

      integer          :: i,j
      double precision :: x,y,r1
      double precision :: pi

      state(:,:,UFA  ) = 0.d0
      state(:,:,UFS  ) = 1.d0
      state(:,:,UFS+1:UFS+nspec-1) = 0.d0

      state(:,:,URHO)  = 1.d0
      state(:,:,UX  )  = 1.d0
      state(:,:,UY  )  = 1.d0

      pi = 4.d0 * datan(1.d0)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
             x = xlo(1) + (i-lo(1)+0.5d0)*dx(1) 
             y = xlo(2) + (j-lo(2)+0.5d0)*dx(2) 
             state(i,j,UX) =  -2.d0 * sin(pi * (x+0.5d0))**2 * sin(pi*(y+0.5d0)) * cos(pi*(y+0.5d0)) 
             state(i,j,UY) =   2.d0 * sin(pi * (x+0.5d0)) * cos(pi * (x+0.5d0)) * sin(pi*(y+0.5d0))**2
         enddo
      enddo

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

             x = xlo(1) + (i-lo(1)+0.5d0)*dx(1) 
             y = xlo(2) + (j-lo(2)+0.5d0)*dx(2) 

             ! Define one blob
             r1 = sqrt( (x-0.5d0)**2+(y-0.75d0)**2 )
             if (r1.lt.0.15d0) then
                state(i,j,UFA) = 1.d0 - 0.5d0 * (1.d0-tanh(30.*(r1-0.1)))
             end if

             state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO)*state(i,j,UFS:UFS+nspec-1)
             state(i,j,UFA            ) = state(i,j,URHO)*state(i,j,UFA            ) 

         enddo
      enddo

      end subroutine initdata

