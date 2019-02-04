! ------------------------------------------------------------------
! $Revision$
! $Date$
! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2014, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! ------------------------------------------------------------------
! This module implements the Fortran 2003 interface to the SUNDIALS
! dense linear solver (dls) matrix structure.
!-------------------------------------------------------------------

module sundials_fdlsmat

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  type, bind(C) :: sundlsmat
     integer(c_int)  :: type  ! dense matrix type
     integer(c_long) :: M     ! number of rows
     integer(c_long) :: N     ! number of columns
     integer(c_long) :: ldim  ! leading dimension
     integer(c_long) :: mu    ! upper bandwidth
     integer(c_long) :: ml    ! lower bandwidth
     integer(c_long) :: s_mu  ! storage upper bandwidth
     type(c_ptr)     :: data  ! pointer to matrix data
     integer(c_long) :: ldata ! length of data array
     type(c_ptr)     :: cols  ! array of pointers to matrix columns
  end type sundlsmat

contains

  subroutine sundlsmat_GetData_Dense(matrix, mdata)
    ! ----------------------------------------------------------------
    ! Description: subroutine to extract data from a dense matrix
    ! ----------------------------------------------------------------

    !======= Declarations =========
    implicit none

    type(sundlsmat)         :: matrix     ! SUNDIALS dense matrix
    real(c_double), pointer :: mdata(:,:) ! dense matrix data

    !======= Internals ============

    ! extract and reshape 1D data array
    call c_f_pointer(matrix%data, mdata, (/matrix%M,matrix%N/))

  end subroutine sundlsmat_GetData_Dense

end module sundials_fdlsmat
