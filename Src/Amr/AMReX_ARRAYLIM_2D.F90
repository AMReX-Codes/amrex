
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_ArrayLim.H>

#define SDIM 2

! ::: -----------------------------------------------------------
! ::: This routine sets the values for the lo() and hi() arrays
! ::: from the ARG_L1, ARG_H1, ... macros.  This is done since
! ::: it is more convenient to use the lo() and hi() arrays.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: DIMS(holder)=> index extent of place holder array
! ::: lo(SDIM)    <= lower index limits
! ::: hi(SDIM)    <= upper index limits
! ::: -----------------------------------------------------------

    subroutine SET_LOHI(DIMS(holder), lo, hi)

      implicit none

!     :::: Passed Variables ::::

      integer DIMDEC(holder)
      integer lo(SDIM), hi(SDIM)

!     --------------------------------------
!     :::: Set Values for lo() and hi() ::::
!     --------------------------------------

      lo(1) = ARG_L1(holder)
      hi(1) = ARG_H1(holder)
      lo(2) = ARG_L2(holder)
      hi(2) = ARG_H2(holder)

    end subroutine SET_LOHI


! ::: -----------------------------------------------------------
! ::: This routine sets the values for the ARG_L1, ARG_H1, ... macros
! ::: from the lo() and hi() arrays.  This is done since
! ::: it is more convenient to use the macros to dimension arrays.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: FF_DIMS(holder) <=  index extent of place holder array
! ::: lo(SDIM)         => lower index limits
! ::: hi(SDIM)         => upper index limits
! ::: -----------------------------------------------------------

    subroutine SET_ARGS(DIMS(holder), lo, hi)

      implicit none

!     :::: Passed Variables ::::

      integer DIMDEC(holder)
      integer lo(SDIM), hi(SDIM)

!     --------------------------------------
!     :::: Set Values for lo() and hi() ::::
!     --------------------------------------

      ARG_L1(holder) = lo(1)
      ARG_H1(holder) = hi(1)
      ARG_L2(holder) = lo(2)
      ARG_H2(holder) = hi(2)

    end subroutine SET_ARGS
