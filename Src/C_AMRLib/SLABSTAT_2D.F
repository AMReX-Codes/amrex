c
c $Id: SLABSTAT_2D.F,v 1.3 2000-05-22 19:54:00 sstanley Exp $
c
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"
#include "SLABSTAT_F.H"

#define SDIM 2

c
c ::: -----------------------------------------------------------
c ::: This is a general routine to calculate the running statistics
c ::: required to calculate the Reynolds averaged first four moments
c ::: of the state components passed in.  This is a general routine 
c ::: and will work for any number of state components.  If the state
c ::: is passed in in the order (a, b, c, d), then the running averages
c ::: which are returned are in the order,
c :::   <a>,   <b>,   <c>,   <d>,
c :::   <a^2>, <b^2>, <c^2>, <d^2>,
c :::   <a^3>, <b^3>, <c^3>, <d^3>,
c :::   <a^4>, <b^4>, <c^4>, <d^4>
c :::
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: runStats      <=  Array to contain the running statistics
c ::: DIMS(runStats) => Dimensions of the running statistics array
c ::: nStats         => Number of components in the statistics array
c ::: state          =) Array containing the state data
c ::: DIMS(state)    => Dimensions of the state data array
c ::: nStateComp     => Number of components in the state array
c ::: dt             => Time step
c ::: delta          => Grid spacing
c ::: -----------------------------------------------------------
c
      subroutine FORT_REYNAVG_FIRST4MOMS(state, DIMS(state), nStateComp,
     $                                   runStats, DIMS(runStats), nStats,
     $                                   dt, delta)

      implicit none

c
c     ::::: passed variables
c
      integer nStats, nStateComp
      integer DIMDEC(state)
      integer DIMDEC(runStats)
      REAL_T dt
      REAL_T delta(SDIM)
      REAL_T runStats(DIMV(runStats),nStats)
      REAL_T state(DIMV(state),nStateComp)

c
c     ::::: local variables
c
      integer i, j, n, nghost
      integer loStats(SDIM), hiStats(SDIM)
      integer loState(SDIM), hiState(SDIM)


c
c     ===================================
c     ::: Set the Values of LO and HI :::
c     ===================================
c
      call SET_LOHI(DIMS(runStats), loStats, hiStats)
      call SET_LOHI(DIMS(state), loState, hiState)

#ifndef NDEBUG
      if (nStats .NE. 4*nStateComp) then
        write(*,1000) nStats, nStateComp
 1000   format('Incorrect relationship between nStats and nStateComp in',
     $        /'FORT_REYNAVG_FIRST4MOMS.  Shoud have nStats = 4*nStateComp',
     $        /'   nStats = ', I5, 5x, 'nStateComp = ', I5)
        stop
      endif

      nghost = 0
      do n = 1, SDIM
        if (loState(n) .GT. loStats(n) - nghost .OR.
     $      hiState(n) .LT. hiStats(n) + nghost) then
          write(*,1010) n, nghost, loStats, hiStats, loState, hiState
 1010     format('Incorrect number of ghost cells in the state date in',
     $          /'FORT_REYNAVG_FIRST4MOMS.',
     $          /'   Direction = ', I2, 5x, 'nghost = ', I2,
     $          /'   loStats = ', SDIM (I2,1x), 5x, 'hiStats = ', SDIM (I2,1x),
     $          /'   loState = ', SDIM (I2,1x), 5x, 'hiState = ', SDIM (I2,1x))
          stop
        endif
      enddo
#endif


c
c     ====================================
c     ::: Calculate Running Statistics :::
c     ====================================
c
      do n = 1, nStateComp
        do j = loStats(2), hiStats(2)
          do i = loStats(1), hiStats(1)
            runStats(i,j,n) = runStats(i,j,n)            + dt * state(i,j,n)
            runStats(i,j,n+nStateComp) = runStats(i,j,n+nStateComp)
     $                                                   + dt * state(i,j,n)**2
            runStats(i,j,n+2*nStateComp) = runStats(i,j,n+2*nStateComp)
     $                                                   + dt * state(i,j,n)**3
            runStats(i,j,n+3*nStateComp) = runStats(i,j,n+3*nStateComp)
     $                                                   + dt * state(i,j,n)**4
          enddo
        enddo
      enddo

c
c
      return
      end




c
c ::: -----------------------------------------------------------
c ::: This is a general routine to calculate the running statistics
c ::: required to calculate the Favre averaged first four moments
c ::: of the state components passed in.  This is a general routine 
c ::: and will work for any number of state components.  The first state 
c ::: component is assumed to be density, so it is used in the Favre averages.
c ::: If the state is passed in in the order (Rho, a, b, c, d),
c ::: then the running averages which are returned are in the order,
c :::   <Rho>,   <Rho a>,   <Rho b>,   <Rho c>,   <Rho d>,
c :::   <Rho^2>, <Rho a^2>, <Rho b^2>, <Rho c^2>, <Rho d^2>,
c :::   <Rho^3>, <Rho a^3>, <Rho b^3>, <Rho c^3>, <Rho d^3>,
c :::   <Rho^4>, <Rho a^4>, <Rho b^4>, <Rho c^4>, <Rho d^4>
c :::
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: runStats      <=  Array to contain the running statistics
c ::: DIMS(runStats) => Dimensions of the running statistics array
c ::: nStats         => Number of components in the statistics array
c ::: state          =) Array containing the state data
c ::: DIMS(state)    => Dimensions of the state data array
c ::: nStateComp     => Number of components in the state array
c ::: dt             => Time step
c ::: delta          => Grid spacing
c ::: -----------------------------------------------------------
c
      subroutine FORT_FAVREAVG_FIRST4MOMS(state, DIMS(state), nStateComp,
     $                                    runStats, DIMS(runStats), nStats,
     $                                    dt, delta)

      implicit none

c
c     ::::: passed variables
c
      integer nStats, nStateComp
      integer DIMDEC(state)
      integer DIMDEC(runStats)
      REAL_T dt
      REAL_T delta(SDIM)
      REAL_T runStats(DIMV(runStats),nStats)
      REAL_T state(DIMV(state),nStateComp)

c
c     ::::: local variables
c
      integer i, j, n, nghost
      integer loStats(SDIM), hiStats(SDIM)
      integer loState(SDIM), hiState(SDIM)


c
c     ===================================
c     ::: Set the Values of LO and HI :::
c     ===================================
c
      call SET_LOHI(DIMS(runStats), loStats, hiStats)
      call SET_LOHI(DIMS(state), loState, hiState)

#ifndef NDEBUG
      if (nStats .NE. 4*nStateComp) then
        write(*,1000) nStats, nStateComp
 1000   format('Incorrect relationship between nStats and nStateComp in',
     $        /'FORT_FAVREAVG_FIRST4MOMS.  Shoud have nStats = 4*nStateComp',
     $        /'   nStats = ', I5, 5x, 'nStateComp = ', I5)
        stop
      endif

      nghost = 0
      do n = 1, SDIM
        if (loState(n) .GT. loStats(n) - nghost .OR.
     $      hiState(n) .LT. hiStats(n) + nghost) then
          write(*,1010) n, nghost, loStats, hiStats, loState, hiState
 1010     format('Incorrect number of ghost cells in the state date in',
     $          /'FORT_FAVREAVG_FIRST4MOMS.',
     $          /'   Direction = ', I2, 5x, 'nghost = ', I2,
     $          /'   loStats = ', SDIM (I2,1x), 5x, 'hiStats = ', SDIM (I2,1x),
     $          /'   loState = ', SDIM (I2,1x), 5x, 'hiState = ', SDIM (I2,1x))
          stop
        endif
      enddo
#endif


c
c     ====================================
c     ::: Calculate Running Statistics :::
c     ====================================
c
      do j = loStats(2), hiStats(2)
        do i = loStats(1), hiStats(1)
          runStats(i,j,1) = runStats(i,j,1)              + dt * state(i,j,1)
          runStats(i,j,1+nStateComp) = runStats(i,j,1+nStateComp)
     $                                                   + dt * state(i,j,1)**2
          runStats(i,j,1+2*nStateComp) = runStats(i,j,1+2*nStateComp)
     $                                                   + dt * state(i,j,1)**3
          runStats(i,j,1+3*nStateComp) = runStats(i,j,1+3*nStateComp)
     $                                                   + dt * state(i,j,1)**4
        enddo
      enddo

      do n = 2, nStateComp
        do j = loStats(2), hiStats(2)
          do i = loStats(1), hiStats(1)
            runStats(i,j,n) = runStats(i,j,n)
     $                                  + dt * (state(i,j,1) * state(i,j,n))
            runStats(i,j,n+nStateComp) = runStats(i,j,n+nStateComp)
     $                                  + dt * (state(i,j,1) * state(i,j,n))**2
            runStats(i,j,n+2*nStateComp) = runStats(i,j,n+2*nStateComp)
     $                                  + dt * (state(i,j,1) * state(i,j,n))**3
            runStats(i,j,n+3*nStateComp) = runStats(i,j,n+3*nStateComp)
     $                                  + dt * (state(i,j,1) * state(i,j,n))**4
          enddo
        enddo
      enddo

c
c
      return
      end

