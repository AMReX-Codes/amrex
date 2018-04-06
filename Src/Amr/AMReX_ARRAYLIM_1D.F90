
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_ArrayLim.H>

#define SDIM 1

! ::: -----------------------------------------------------------
! ::: This routine sets the values for the lo() and hi() arrays
! ::: from the ARG_L1, ARG_H1, ... macros.  This is done since
! ::: it is more convenient to use the lo() and hi() arrays.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: holder_l1,holder_h1=> index extent of place holder array
! ::: lo(SDIM)    <= lower index limits
! ::: hi(SDIM)    <= upper index limits
! ::: -----------------------------------------------------------

    subroutine SET_LOHI(holder_l1,holder_h1, lo, hi)

      implicit none

!     :::: Passed Variables ::::

      integer holder_l1,holder_h1
      integer lo(SDIM), hi(SDIM)

!     --------------------------------------
!     :::: Set Values for lo() and hi() ::::
!     --------------------------------------

      lo(1) = holder_l1
      hi(1) = holder_h1

    end subroutine SET_LOHI


! ::: -----------------------------------------------------------
! ::: This routine sets the values for the ARG_L1, ARG_H1, ... macros
! ::: from the lo() and hi() arrays.  This is done since
! ::: it is more convenient to use the macros to dimension arrays.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: FF_holder_l1,holder_h1 <=  index extent of place holder array
! ::: lo(SDIM)         => lower index limits
! ::: hi(SDIM)         => upper index limits
! ::: -----------------------------------------------------------

    subroutine SET_ARGS(holder_l1,holder_h1, lo, hi)

      implicit none

!     :::: Passed Variables ::::

      integer holder_l1,holder_h1
      integer lo(SDIM), hi(SDIM)

!     --------------------------------------
!     :::: Set Values for lo() and hi() ::::
!     --------------------------------------

      holder_l1 = lo(1)
      holder_h1 = hi(1)

    end subroutine SET_ARGS
