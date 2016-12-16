
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "MACOPERATOR_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 2

c :: ----------------------------------------------------------
c :: MACCOEF
c ::             Compute the coefficents for MAC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  cx,cy       <=  edge coef arrays
c ::  DIMS(cx)     => index limits for cx
c ::  DIMS(cy)     => index limits for cy
c ::  lo,hi        => index limits for rhs
c ::  ax,ay        => edge based area arrays
c ::  DIMS(ax)     => index limits for ax
c ::  DIMS(ay)     => index limits for ay
c ::  rho          => cell centered density array
c ::  DIMS(rho)    => index limits for rho
c ::  dx           => cell size
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACCOEF (cx,DIMS(cx),cy,DIMS(cy),
     &                          ax,DIMS(ax),ay,DIMS(ay),
     &                          rho,DIMS(rho),lo,hi,dx)
       implicit none
       integer DIMDEC(cx)
       integer DIMDEC(cy)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(rho)
       integer lo(SDIM), hi(SDIM)
       REAL_T  dx(SDIM)
       REAL_T  cx(DIMV(cx))
       REAL_T  cy(DIMV(cy))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  rho(DIMV(rho))

       integer i, j
       REAL_T rhoavg

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            if (rho(i,j) .lt. zero) then
	      print *,' '
	      print *,'TESTING in MACCOEF '
	      print *,'RHO HAS GONE NEGATIVE AT ',i,j,rho(i,j)
              call bl_abort(" ")
            end if
          end do
       end do
c
c      ::::: finish coef in X direction (part 2)
c
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             rhoavg = half * (rho(i,j) + rho(i-1,j))
             cx(i,j) = dx(1) * ax(i,j) / rhoavg
          end do
       end do
c
c      ::::: finish coef in Y direction (part 2)
c
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             rhoavg = half * (rho(i,j) + rho(i,j-1))
             cy(i,j) = dx(2) * ay(i,j) / rhoavg
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: MACRHS
c ::             Compute the RHS for MAC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  ux,uy       <=  edge velocity arrays
c ::  DIMS(ux)     => index limits for ux
c ::  DIMS(uy)     => index limits for uy
c ::  lo,hi        => index limits for rhs
c ::  ax,ay        => edge based area arrays
c ::  DIMS(ax)     => index limits for ax
c ::  DIMS(ay)     => index limits for ay
c ::  vol          => cell centered volume array
c ::  vlo,vhi      => index limits of vol array
c ::  DIMS(vol)    => index limits for vol
c ::  rhs         <=> cell centered rhs array
c ::  DIMS(rhs)    => index limits for rhs
c ::  scale        => scale factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACRHS (ux,DIMS(ux),uy,DIMS(uy),
     &                         ax,DIMS(ax),ay,DIMS(ay),
     &                         vol,DIMS(vol),rhs,DIMS(rhs),
     &                         lo,hi,scale)
       implicit none
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(vol)
       integer DIMDEC(rhs)
       integer lo(SDIM), hi(SDIM)
       REAL_T  scale
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  vol(DIMV(vol))
       REAL_T  rhs(DIMV(rhs))

       integer i, j
       REAL_T  divu
c
c      ::::: rhs holds the divergence condition (possibly zero)
c
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             divu = ax(i+1,j)*ux(i+1,j) - ax(i,j)*ux(i,j)
     &            + ay(i,j+1)*uy(i,j+1) - ay(i,j)*uy(i,j)
             rhs(i,j) = scale*(divu - vol(i,j)*rhs(i,j))
          end do
       end do
 
       end

c :: ----------------------------------------------------------
c :: MACUPDATE
c ::             Compute the update to velocity field to
c ::             make it divergence free
c ::
c :: INPUTS / OUTPUTS:
c ::  ux,uy       <=  edge based velocity arrays
c ::  DIMS(ux)     => index limits for ux
c ::  DIMS(uy)     => index limits for uy
c ::  phi          => soln from MAC project
c ::  DIMS(phi)    => index limits for phi
c ::  rho          => density at time N
c ::  DIMS(rho)    => index limits for rho
c ::  dx           => cell size
c ::  mult         => scalar multiplier
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACUPDATE(init,ux,DIMS(ux),uy,DIMS(uy),
     &                           phi,DIMS(phi),rho,DIMS(rho),
     &                           lo,hi,dx,mult)

       implicit none
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(phi)
       integer DIMDEC(rho)
       integer lo(SDIM), hi(SDIM)
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  phi(DIMV(phi))
       REAL_T  rho(DIMV(rho))
       REAL_T  dx(SDIM), mult
       integer init

       integer i, j
       REAL_T  rhoavg, gp
c
c     set values to 0.0 if initializing
c
       if ( init .eq. 1 ) then
          do j = ARG_L2(ux), ARG_H2(ux)
             do i = ARG_L1(ux), ARG_H1(ux)
                ux(i,j) = zero
             end do
          end do
          do j = ARG_L2(uy), ARG_H2(uy)
             do i = ARG_L1(uy), ARG_H1(uy)
                uy(i,j) = zero
             end do
          end do
       end if
c
c     compute the x mac gradient
c
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             rhoavg = half*(rho(i,j) + rho(i-1,j))
             gp = (phi(i,j)-phi(i-1,j))/dx(1)
             ux(i,j) = ux(i,j) + mult * gp / rhoavg
          end do
       end do
c
c      compute the y mac gradient
c
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             rhoavg = half*(rho(i,j) + rho(i,j-1))
             gp = (phi(i,j)-phi(i,j-1))/dx(2)
             uy(i,j) = uy(i,j) + mult * gp / rhoavg
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: MACSYNCRHS
c ::        Modify the RHS for MAC SYNC solve
c ::
c :: INPUTS / OUTPUTS:
c ::  rhs         <=  right hand side array
c ::  DIMS(rhs)    => index limits for rhs
c ::  vol          => cell centered volume array
c ::  DIMS(vol)    => index limits of vol array
c ::  rhsscale     => const multiplier to rhs
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACSYNCRHS(rhs,DIMS(rhs),lo,hi,
     &                            vol,DIMS(vol),rhsscale)
       implicit none
       integer DIMDEC(rhs)
       integer DIMDEC(vol)
       integer lo(SDIM), hi(SDIM)
       REAL_T  rhsscale
       REAL_T  rhs(DIMV(rhs))
       REAL_T  vol(DIMV(vol))

       integer i, j
c
c      ::::: multiply by volume since reflux step (which computed rhs)
c      ::::: divided by volume.
c
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhs(i,j) = rhsscale*vol(i,j)*rhs(i,j)
          end do
       end do

       end
