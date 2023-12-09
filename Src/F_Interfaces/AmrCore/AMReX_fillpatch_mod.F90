#include <AMReX_Config.H>

module amrex_fillpatch_module

  use amrex_base_module

  implicit none
  private

  public :: amrex_fillpatch, amrex_fillcoarsepatch, amrex_interp_hook_proc

  interface amrex_fillpatch
     module procedure amrex_fillpatch_single
     module procedure amrex_fillpatch_two
     module procedure amrex_fillpatch_dgconservative_two
     module procedure amrex_fillpatch_dgpointwise_two
     module procedure amrex_fillpatch_two_faces
  end interface amrex_fillpatch

  interface amrex_fillcoarsepatch
     module procedure amrex_fillcoarsepatch_default
     module procedure amrex_fillcoarsepatch_dgconservative
     module procedure amrex_fillcoarsepatch_dgpointwise
     module procedure amrex_fillcoarsepatch_faces
  end interface amrex_fillcoarsepatch

  interface
     subroutine amrex_interp_hook_proc (lo, hi, d, dlo, dhi, nd, icomp, ncomp) bind(c)
       import
       implicit none
       integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
       integer(c_int), intent(in), value :: nd, icomp, ncomp
       real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nd)
     end subroutine amrex_interp_hook_proc

     subroutine amrex_interp_hook_arr_proc (lo, hi, dx, dxlo, dxhi, &
#if (AMREX_SPACEDIM > 1)
          &                                         dy, dylo, dyhi, &
#endif
#if (AMREX_SPACEDIM > 2)
          &                                         dz, dzlo, dzhi, &
#endif
          &                                 nd, icomp, ncomp) bind(c)
       import
       implicit none
       integer(c_int), intent(in) :: lo(3), hi(3), dxlo(3), dxhi(3)
       integer(c_int), intent(in), value :: nd, icomp, ncomp
       real(amrex_real), intent(inout) :: dx(dxlo(1):dxhi(1),dxlo(2):dxhi(2),dxlo(3):dxhi(3),nd)
#if (AMREX_SPACEDIM > 1)
       integer(c_int), intent(in) :: dylo(3), dyhi(3)
       real(amrex_real), intent(inout) :: dy(dylo(1):dyhi(1),dylo(2):dyhi(2),dylo(3):dyhi(3),nd)
#endif
#if (AMREX_SPACEDIM > 2)
       integer(c_int), intent(in) :: dzlo(3), dzhi(3)
       real(amrex_real), intent(inout) :: dz(dzlo(1):dzhi(1),dzlo(2):dzhi(2),dzlo(3):dzhi(3),nd)
#endif
     end subroutine amrex_interp_hook_arr_proc
  end interface

  interface
     subroutine amrex_fi_fillpatch_single(mf, time, smf, stime, ns, scomp, dcomp, ncomp, &
          geom, fill) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom
       type(c_ptr), intent(in) :: smf(*)
       real(amrex_real), value :: time
       real(amrex_real), intent(in) :: stime(*)
       integer(c_int), value :: scomp, dcomp, ncomp, ns
       type(c_funptr), value :: fill
     end subroutine amrex_fi_fillpatch_single

     SUBROUTINE amrex_fi_fillpatch_two &
       ( pMF, Time, &
         pCrseMF, CrseTime, nCrse, pFineMF, FineTime, nFine, &
         sComp, dComp, nComp, &
         pCrseGeom, pFineGeom, fpCrseFillPhysBC, fpFineFillPhysBC, &
         RefRatio, interp, pLoBC, pHiBC, pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE      :: pMF, pCrseGeom, pFineGeom
       TYPE(C_PTR)     , INTENT(in) :: pCrseMF(*), pFineMF(*), pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE      :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                       pre_interp, post_interp
       REAL(amrex_real), VALUE      :: Time
       REAL(amrex_real), INTENT(in) :: CrseTime(*), FineTime(*)
       INTEGER         , VALUE      :: nCrse, nFine, sComp, dComp, nComp, RefRatio, interp
     END SUBROUTINE amrex_fi_fillpatch_two

     SUBROUTINE amrex_fi_fillpatch_dgconservative_two &
       ( pMF, pMF_G, Time, &
         pCrseMF, pCrseMF_G, CrseTime, nCrse, &
         pFineMF, pFineMF_G, FineTime, nFine, &
         sComp, dComp, nComp, &
         pCrseGeom, pFineGeom, fpCrseFillPhysBC, fpFineFillPhysBC, &
         RefRatio, interp, pLoBC, pHiBC, nFineV, nDOFX, &
         vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE      :: pMF, pMF_G, pCrseGeom, pFineGeom, &
                                       vpCoarseToFineProjectionMatrix
       TYPE(C_PTR)     , INTENT(in) :: pCrseMF(*), pCrseMF_G(*), &
                                       pFineMF(*), pFineMF_G(*), &
                                       pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE      :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                       pre_interp, post_interp
       REAL(amrex_real), VALUE      :: Time
       REAL(amrex_real), INTENT(in) :: CrseTime(*), FineTime(*)
       INTEGER         , VALUE      :: nCrse, nFine, sComp, dComp, nComp, &
                                       RefRatio, interp, nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillpatch_dgconservative_two

     SUBROUTINE amrex_fi_fillpatch_dgpointwise_two &
       ( pMF, Time, &
         pCrseMF, CrseTime, nCrse, &
         pFineMF, FineTime, nFine, &
         sComp, dComp, nComp, &
         pCrseGeom, pFineGeom, fpCrseFillPhysBC, fpFineFillPhysBC, &
         RefRatio, interp, pLoBC, pHiBC, nFineV, nDOFX, &
         vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE      :: pMF, pCrseGeom, pFineGeom, &
                                       vpCoarseToFineProjectionMatrix
       TYPE(C_PTR)     , INTENT(in) :: pCrseMF(*), &
                                       pFineMF(*), &
                                       pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE      :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                       pre_interp, post_interp
       REAL(amrex_real), VALUE      :: Time
       REAL(amrex_real), INTENT(in) :: CrseTime(*), FineTime(*)
       INTEGER         , VALUE      :: nCrse, nFine, sComp, dComp, nComp, &
                                       RefRatio, interp, nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillpatch_dgpointwise_two

     subroutine amrex_fi_fillpatch_two_faces(mf, time, &
          cmf, ctime, nc, fmf, ftime, nf, scomp, dcomp, ncomp, &
          cgeom, fgeom, cfill, ffill, rr, interp, lo_bc, hi_bc, pre_interp, post_interp) &
          bind(c)
       import
       implicit none
       type(c_ptr), value :: cgeom, fgeom
       type(c_ptr), intent(in) :: mf(*)
       type(c_ptr), intent(in) :: cmf(*), fmf(*), lo_bc(*), hi_bc(*)
       type(c_funptr), intent(in) :: cfill(*), ffill(*)
       type(c_funptr), value :: pre_interp, post_interp
       real(amrex_real), value :: time
       real(amrex_real), intent(in) :: ctime(*), ftime(*)
       integer, value :: nc, nf, scomp, dcomp, ncomp, rr, interp
     end subroutine amrex_fi_fillpatch_two_faces

     SUBROUTINE amrex_fi_fillcoarsepatch &
       ( pMF, Time, &
         pCrseMF, sComp, dComp, nComp, pCrseGeom, pFineGeom, &
         fpCrseFillPhysBC, fpFineFillPhysBC, RefRatio, interp, pLoBC, pHiBC, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: pMF, pCrseMF, pCrseGeom, pFineGeom
       TYPE(C_PTR), INTENT(in) :: pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                  pre_interp, post_interp
       REAL(amrex_real), VALUE :: Time
       INTEGER         , VALUE :: sComp, dComp, nComp, RefRatio, interp
     END SUBROUTINE amrex_fi_fillcoarsepatch

     SUBROUTINE amrex_fi_fillcoarsepatch_dgconservative &
       ( pMF, pMF_G, Time, &
         pCrseMF, pCrseMF_G, sComp, dComp, nComp, pCrseGeom, pFineGeom, &
         fpCrseFillPhysBC, fpFineFillPhysBC, RefRatio, interp, pLoBC, pHiBC, &
         nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: pMF, pMF_G, pCrseMF, pCrseMF_G, &
                                  pCrseGeom, pFineGeom, &
                                  vpCoarseToFineProjectionMatrix
       TYPE(C_PTR), INTENT(in) :: pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                  pre_interp, post_interp
       REAL(amrex_real), VALUE :: Time
       INTEGER         , VALUE :: sComp, dComp, nComp, RefRatio, interp, &
                                  nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillcoarsepatch_dgconservative

     SUBROUTINE amrex_fi_fillcoarsepatch_dgpointwise &
       ( pMF, Time, &
         pCrseMF, sComp, dComp, nComp, pCrseGeom, pFineGeom, &
         fpCrseFillPhysBC, fpFineFillPhysBC, RefRatio, interp, pLoBC, pHiBC, &
         nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: pMF, pCrseMF, &
                                  pCrseGeom, pFineGeom, &
                                  vpCoarseToFineProjectionMatrix
       TYPE(C_PTR), INTENT(in) :: pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                  pre_interp, post_interp
       REAL(amrex_real), VALUE :: Time
       INTEGER         , VALUE :: sComp, dComp, nComp, RefRatio, interp, &
                                  nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillcoarsepatch_dgpointwise

     subroutine amrex_fi_fillcoarsepatch_faces(mf, time, cmf, scomp, dcomp, ncomp, &
          cgeom, fgeom, cfill, ffill, rr, interp, lo_bc, hi_bc, pre_interp, post_interp) &
          bind(c)
       import
       implicit none
       type(c_ptr), intent(in) :: mf(*), cmf(*)
       type(c_ptr), value :: cgeom, fgeom
       type(c_ptr), intent(in) :: lo_bc(*), hi_bc(*)
       type(c_funptr), intent(in) :: cfill(*), ffill(*)
       type(c_funptr), value :: pre_interp, post_interp
       real(amrex_real), value :: time
       integer, value :: scomp, dcomp, ncomp, rr, interp
     end subroutine amrex_fi_fillcoarsepatch_faces
  end interface

contains

  subroutine amrex_fillpatch_single (mf, told, mfold, tnew, mfnew, geom, fill_physbc, &
       &                             time, scomp, dcomp, ncomp)
    type(amrex_multifab), intent(inout) :: mf
    type(amrex_multifab), intent(in   ) :: mfold, mfnew
    integer, intent(in) :: scomp, dcomp, ncomp
    real(amrex_real), intent(in) :: told, tnew, time
    type(amrex_geometry), intent(in) :: geom
    procedure(amrex_physbc_proc) :: fill_physbc

    real(amrex_real) :: teps
    real(amrex_real) :: stime(2)
    type(c_ptr) :: smf(2)
    integer :: ns

    teps = 1.e-4_amrex_real * abs(tnew - told)
    if (abs(time-tnew) .le. teps) then
       ns = 1
       smf  (1) = mfnew%p
       stime(1) =  tnew
    else if (abs(time-told) .le. teps) then
       ns = 1
       smf  (1) = mfold%p
       stime(1) =  told
    else
       ns = 2
       smf  (1) = mfold%p
       smf  (2) = mfnew%p
       stime(1) =  told
       stime(2) =  tnew
    end if

    ! scomp-1 and dcomp-1 because of Fortran index starts with 1
    call amrex_fi_fillpatch_single(mf%p, time, smf, stime, ns, scomp-1, dcomp-1, ncomp, geom%p, &
         &                         c_funloc(fill_physbc))

  end subroutine amrex_fillpatch_single

  SUBROUTINE amrex_fillpatch_two &
    ( MF, &
      OldTimeCrse, OldMFCrse, NewTimeCrse, NewMFCrse, &
      GeomCrse, FillPhysBCCrse, &
      OldTimeFine, OldMFFine, NewTimeFine, NewMFFine, &
      GeomFine, FillPhysBCFine, &
      Time, sComp, dComp, nComp, RefRatio, interp, LoBC, HiBC, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: &
      OldMFCrse, NewMFCrse, OldMFFine, NewMFFine
    INTEGER             , INTENT(in)    :: &
      sComp, dComp, nComp, RefRatio, interp
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,sComp+nComp-1), HiBC(amrex_spacedim,sComp+nComp-1)
    REAL(amrex_real)    , INTENT(in)    :: &
      OldTimeCrse, NewTimeCrse, OldTimeFine, NewTimeFine, Time
    TYPE(amrex_geometry), INTENT(in)    :: &
      GeomCrse, GeomFine
    PROCEDURE(amrex_physbc_proc)        :: &
      FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: &
      pre_interp, post_interp

    REAL(amrex_real) :: teps
    REAL(amrex_real) :: CrseTime(2), FineTime(2)
    TYPE(c_ptr)      :: pCrseMF (2), pFineMF (2)
    TYPE(c_ptr)      :: pLoBC(sComp+nComp-1), pHiBC(sComp+nComp-1)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: nCrse, nFine, iComp

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       nCrse       = 1
       pCrseMF (1) = NewMFCrse   % p
       CrseTime(1) = NewTimeCrse

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       nCrse       = 1
       pCrseMF (1) = OldMFCrse   % p
       CrseTime(1) = OldTimeCrse

    ELSE

       nCrse       = 2
       pCrseMF (1) = OldMFCrse   % p
       pCrseMF (2) = NewMFCrse   % p
       CrseTime(1) = OldTimeCrse
       CrseTime(2) = NewTimeCrse

    END IF

    ! Fine level
    teps = 1.0e-4_amrex_real * ABS( NewTimeFine - OldTimeFine )
    IF( ABS( Time - NewTimeFine ) .LE. teps )THEN

       nFine       = 1
       pFineMF (1) = NewMFFine   % p
       FineTime(1) = NewTimeFine

    ELSE IF( ABS( Time - OldTimeFine ) .LE. teps )THEN

       nFine       = 1
       pFineMF (1) = OldMFFine   % p
       FineTime(1) = OldTimeFine

    ELSE

       nFine       = 2
       pFineMF (1) = OldMFFine   % p
       pFineMF (2) = NewMFFine   % p
       FineTime(1) = OldTimeFine
       FineTime(2) = NewTimeFine

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillpatch_two &
           ( MF % p, Time, pCrseMF, CrseTime, nCrse, pFineMF, FineTime, nFine, &
             sComp-1, dComp-1, nComp, GeomCrse % p, GeomFine % p, &
             C_FUNLOC( FillPhysBCCrse ), C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC,&
             pre_interp_ptr, post_interp_ptr )

  END SUBROUTINE amrex_fillpatch_two

  SUBROUTINE amrex_fillpatch_dgconservative_two &
    ( MF, MF_G, &
      OldTimeCrse, OldMFCrse, OldMFCrse_G, NewTimeCrse, NewMFCrse, NewMFCrse_G, &
      GeomCrse, FillPhysBCCrse, &
      OldTimeFine, OldMFFine, OldMFFine_G, NewTimeFine, NewMFFine, NewMFFine_G, &
      GeomFine, FillPhysBCFine, &
      Time, sComp, dComp, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: &
      OldMFCrse, OldMFCrse_G, NewMFCrse, NewMFCrse_G, &
      OldMFFine, OldMFFine_G, NewMFFine, NewMFFine_G, MF_G
    INTEGER             , INTENT(in)    :: &
      sComp, dComp, nComp, RefRatio, interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,sComp+nComp-1), HiBC(amrex_spacedim,sComp+nComp-1)
    REAL(amrex_real)    , INTENT(in)    :: &
      OldTimeCrse, NewTimeCrse, OldTimeFine, NewTimeFine, Time
    TYPE(amrex_geometry), INTENT(in)    :: &
      GeomCrse, GeomFine
    TYPE(c_ptr)         , INTENT(in)    :: &
      vpCoarseToFineProjectionMatrix
    PROCEDURE(amrex_physbc_proc)        :: &
      FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: &
      pre_interp, post_interp

    REAL(amrex_real) :: teps
    REAL(amrex_real) :: CrseTime (2), FineTime (2)
    TYPE(c_ptr)      :: pCrseMF  (2), pFineMF  (2)
    TYPE(c_ptr)      :: pCrseMF_G(2), pFineMF_G(2)
    TYPE(c_ptr)      :: pLoBC(sComp+nComp-1), pHiBC(sComp+nComp-1)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: nCrse, nFine, iComp

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       nCrse        = 1
       pCrseMF  (1) = NewMFCrse   % p
       pCrseMF_G(1) = NewMFCrse_G % p
       CrseTime (1) = NewTimeCrse

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       nCrse        = 1
       pCrseMF  (1) = OldMFCrse   % p
       pCrseMF_G(1) = OldMFCrse_G % p
       CrseTime (1) = OldTimeCrse

    ELSE

       nCrse        = 2
       pCrseMF  (1) = OldMFCrse   % p
       pCrseMF  (2) = NewMFCrse   % p
       pCrseMF_G(1) = OldMFCrse_G % p
       pCrseMF_G(2) = NewMFCrse_G % p
       CrseTime (1) = OldTimeCrse
       CrseTime (2) = NewTimeCrse

    END IF

    ! Fine level
    teps = 1.0e-4_amrex_real * ABS( NewTimeFine - OldTimeFine )
    IF( ABS( Time - NewTimeFine ) .LE. teps )THEN

       nFine        = 1
       pFineMF  (1) = NewMFFine   % p
       pFineMF_G(1) = NewMFFine_G % p
       FineTime (1) = NewTimeFine

    ELSE IF( ABS( Time - OldTimeFine ) .LE. teps )THEN

       nFine        = 1
       pFineMF  (1) = OldMFFine   % p
       pFineMF_G(1) = OldMFFine_G % p
       FineTime (1) = OldTimeFine

    ELSE

       nFine        = 2
       pFineMF  (1) = OldMFFine   % p
       pFineMF  (2) = NewMFFine   % p
       pFineMF_G(1) = OldMFFine_G % p
       pFineMF_G(2) = NewMFFine_G % p
       FineTime (1) = OldTimeFine
       FineTime (2) = NewTimeFine

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillpatch_dgconservative_two &
           ( MF % p, MF_G % p, Time, &
             pCrseMF, pCrseMF_G, CrseTime, nCrse, &
             pFineMF, pFineMF_G, FineTime, nFine, &
             sComp-1, dComp-1, nComp, GeomCrse % p, GeomFine % p, &
             C_FUNLOC( FillPhysBCCrse ), C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr )

  END SUBROUTINE amrex_fillpatch_dgconservative_two

  SUBROUTINE amrex_fillpatch_dgpointwise_two &
    ( MF, &
      OldTimeCrse, OldMFCrse, NewTimeCrse, NewMFCrse, &
      GeomCrse, FillPhysBCCrse, &
      OldTimeFine, OldMFFine, NewTimeFine, NewMFFine, &
      GeomFine, FillPhysBCFine, &
      Time, sComp, dComp, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: &
      OldMFCrse, NewMFCrse, &
      OldMFFine, NewMFFine
    INTEGER             , INTENT(in)    :: &
      sComp, dComp, nComp, RefRatio, interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,sComp+nComp-1), HiBC(amrex_spacedim,sComp+nComp-1)
    REAL(amrex_real)    , INTENT(in)    :: &
      OldTimeCrse, NewTimeCrse, OldTimeFine, NewTimeFine, Time
    TYPE(amrex_geometry), INTENT(in)    :: &
      GeomCrse, GeomFine
    TYPE(c_ptr)         , INTENT(in)    :: &
      vpCoarseToFineProjectionMatrix
    PROCEDURE(amrex_physbc_proc)        :: &
      FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: &
      pre_interp, post_interp

    REAL(amrex_real) :: teps
    REAL(amrex_real) :: CrseTime (2), FineTime (2)
    TYPE(c_ptr)      :: pCrseMF  (2), pFineMF  (2)
    TYPE(c_ptr)      :: pLoBC(sComp+nComp-1), pHiBC(sComp+nComp-1)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: nCrse, nFine, iComp

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       nCrse       = 1
       pCrseMF (1) = NewMFCrse % p
       CrseTime(1) = NewTimeCrse

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       nCrse       = 1
       pCrseMF (1) = OldMFCrse % p
       CrseTime(1) = OldTimeCrse

    ELSE

       nCrse       = 2
       pCrseMF (1) = OldMFCrse  % p
       pCrseMF (2) = NewMFCrse  % p
       CrseTime(1) = OldTimeCrse
       CrseTime(2) = NewTimeCrse

    END IF

    ! Fine level
    teps = 1.0e-4_amrex_real * ABS( NewTimeFine - OldTimeFine )
    IF( ABS( Time - NewTimeFine ) .LE. teps )THEN

       nFine       = 1
       pFineMF (1) = NewMFFine % p
       FineTime(1) = NewTimeFine

    ELSE IF( ABS( Time - OldTimeFine ) .LE. teps )THEN

       nFine       = 1
       pFineMF (1) = OldMFFine % p
       FineTime(1) = OldTimeFine

    ELSE

       nFine       = 2
       pFineMF (1) = OldMFFine % p
       pFineMF (2) = NewMFFine % p
       FineTime(1) = OldTimeFine
       FineTime(2) = NewTimeFine

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillpatch_dgpointwise_two &
           ( MF % p, Time, &
             pCrseMF, CrseTime, nCrse, &
             pFineMF, FineTime, nFine, &
             sComp-1, dComp-1, nComp, GeomCrse % p, GeomFine % p, &
             C_FUNLOC( FillPhysBCCrse ), C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr )

  END SUBROUTINE amrex_fillpatch_dgpointwise_two

  subroutine amrex_fillpatch_two_faces(mf, told_c, mfold_c, tnew_c, mfnew_c, geom_c, fill_physbc_cx, &
#if (AMREX_SPACEDIM > 1)
        &                                   fill_physbc_cy, &
#if (AMREX_SPACEDIM > 2)
        &                                   fill_physbc_cz, &
#endif
#endif
        &                                   told_f, mfold_f, tnew_f, mfnew_f, geom_f, fill_physbc_fx, &
#if (AMREX_SPACEDIM > 1)
        &                                   fill_physbc_fy, &
#if (AMREX_SPACEDIM > 2)
        &                                   fill_physbc_fz, &
#endif
#endif
        &                               time, scomp, dcomp, ncomp, rr, interp, lo_bc, hi_bc, &
        &                               pre_interp, post_interp)
    type(amrex_multifab), intent(inout) :: mf(amrex_spacedim)
    type(amrex_multifab), intent(in   ) :: mfold_c(amrex_spacedim), mfnew_c(amrex_spacedim)
    type(amrex_multifab), intent(in   ) :: mfold_f(amrex_spacedim), mfnew_f(amrex_spacedim)
    integer, intent(in) :: scomp, dcomp, ncomp, rr, interp
    !                 (BC dir        , comp        , MF)
    integer, dimension(amrex_spacedim,scomp+ncomp-1,amrex_spacedim), target, intent(in) :: lo_bc, hi_bc
    real(amrex_real), intent(in) :: told_c, tnew_c, told_f, tnew_f, time
    type(amrex_geometry), intent(in) :: geom_c, geom_f
    procedure(amrex_physbc_proc) :: fill_physbc_cx, fill_physbc_fx
#if (AMREX_SPACEDIM > 1)
    procedure(amrex_physbc_proc) :: fill_physbc_cy, fill_physbc_fy
#endif
#if (AMREX_SPACEDIM > 2)
    procedure(amrex_physbc_proc) :: fill_physbc_cz, fill_physbc_fz
#endif
    procedure(amrex_interp_hook_arr_proc), optional :: pre_interp
    procedure(amrex_interp_hook_arr_proc), optional :: post_interp

    real(amrex_real) :: teps
    real(amrex_real) :: c_time(2), f_time(2)
    type(c_ptr) :: faces(amrex_spacedim)
    type(c_ptr) :: c_mf(amrex_spacedim*2), f_mf(amrex_spacedim*2)
    type(c_funptr) :: cfill(amrex_spacedim), ffill(amrex_spacedim)
    type(c_ptr) :: lo_bc_ptr(amrex_spacedim*(scomp+ncomp-1))
    type(c_ptr) :: hi_bc_ptr(amrex_spacedim*(scomp+ncomp-1))
    type(c_funptr) :: pre_interp_ptr, post_interp_ptr
    integer :: ncrse, nfine, dim, mfid, nc, i

    cfill(1) = c_funloc(fill_physbc_cx)
    ffill(1) = c_funloc(fill_physbc_fx)
#if AMREX_SPACEDIM>=2
    cfill(2) = c_funloc(fill_physbc_cy)
    ffill(2) = c_funloc(fill_physbc_fy)
#if AMREX_SPACEDIM>=3
    cfill(3) = c_funloc(fill_physbc_cz)
    ffill(3) = c_funloc(fill_physbc_fz)
#endif
#endif

    do dim = 1, amrex_spacedim
       faces(dim) = mf(dim)%p
    end do

    ! coarse level
    teps = 1.e-4_amrex_real * abs(tnew_c - told_c)
    if (abs(time-tnew_c) .le. teps) then
       ncrse= 1
       c_time(1) =  tnew_c
       do dim = 1, amrex_spacedim
          c_mf(dim) = mfnew_c(dim)%p
       end do
    else if (abs(time-told_c) .le. teps) then
       ncrse= 1
       c_time(1) =  told_c
       do dim = 1, amrex_spacedim
          c_mf(dim) = mfold_c(dim)%p
       end do
    else
       ncrse= 2
       c_time(1) =  told_c
       c_time(2) =  tnew_c
       do dim = 1, amrex_spacedim
          c_mf(dim)                = mfold_c(dim)%p
          c_mf(dim+amrex_spacedim) = mfnew_c(dim)%p
       end do
    end if

    ! fine level
    teps = 1.e-4_amrex_real * abs(tnew_f - told_f)
    if (abs(time-tnew_f) .le. teps) then
       nfine= 1
       f_time(1) =  tnew_f
       do dim = 1, amrex_spacedim
          f_mf(dim) = mfnew_f(dim)%p
       enddo
    else if (abs(time-told_f) .le. teps) then
       nfine= 1
       f_time(1) =  told_f
       do dim = 1, amrex_spacedim
          f_mf(dim) = mfold_f(dim)%p
       end do
    else
       nfine= 2
       f_time(1) =  told_f
       f_time(2) =  tnew_f
       do dim = 1, amrex_spacedim
          f_mf(dim)                = mfold_f(dim)%p
          f_mf(dim+amrex_spacedim) = mfnew_f(dim)%p
       end do
    end if

    ! lo_bc & hi_bc: (BC dir, comp, MF)
    nc = scomp+ncomp-1
    do mfid = 1, amrex_spacedim
       do i = 1, scomp-1
          lo_bc_ptr((mfid-1)*nc + i) = c_null_ptr
          hi_bc_ptr((mfid-1)*nc + i) = c_null_ptr
       end do
       do i = scomp, nc
          lo_bc_ptr((mfid-1)*nc + i) = c_loc(lo_bc(1,i,mfid))
          hi_bc_ptr((mfid-1)*nc + i) = c_loc(hi_bc(1,i,mfid))
       end do
    end do

    pre_interp_ptr = c_null_funptr
    if (present(pre_interp)) pre_interp_ptr = c_funloc(pre_interp)
    post_interp_ptr = c_null_funptr
    if (present(post_interp)) post_interp_ptr = c_funloc(post_interp)

    ! scomp-1 and dcomp-1 because of Fortran index starts with 1
    call amrex_fi_fillpatch_two_faces(faces, time, c_mf, c_time, ncrse,&
         &                                         f_mf, f_time, nfine,&
         &                            scomp-1, dcomp-1, ncomp,         &
         &                            geom_c%p, geom_f%p,              &
         &                            cfill, ffill,                    &
         &                            rr, interp, lo_bc_ptr, hi_bc_ptr,&
         &                            pre_interp_ptr, post_interp_ptr)
  end subroutine amrex_fillpatch_two_faces


  SUBROUTINE amrex_fillcoarsepatch_default &
    ( MF, &
      OldTimeCrse, OldMFCrse, NewTimeCrse, NewMFCrse, &
      CrseGeom, FillPhysBCCrse, FineGeom, FillPhysBCFine, &
      Time, sComp, dComp, nComp, RefRatio, interp, LoBC, HiBC, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: OldMFCrse, NewMFCrse
    INTEGER             , INTENT(in)    :: sComp, dComp, nComp, RefRatio, interp
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,sComp+nComp-1), HiBC(amrex_spacedim,sComp+nComp-1)
    REAL(amrex_real)    , INTENT(in)    :: OldTimeCrse, NewTimeCrse, Time
    TYPE(amrex_geometry), INTENT(in)    :: CrseGeom, FineGeom
    PROCEDURE(amrex_physbc_proc)        :: FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: pre_interp
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: post_interp

    REAL(amrex_real) :: teps
    TYPE(c_ptr)      :: pCrseMF
    TYPE(c_ptr)      :: pLoBC(sComp+nComp-1), pHiBC(sComp+nComp-1)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: iComp

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       pCrseMF = NewMFCrse % p

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       pCrseMF = OldMFCrse % p

    ELSE

       pCrseMF = NewMFCrse % p

       CALL amrex_abort( "amrex_fillcoarsepatch: how did this happen?" )

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillcoarsepatch &
           ( MF % p, Time, pCrseMF, sComp-1, dComp-1, nComp, &
             CrseGeom % p, FineGeom % p, &
             C_FUNLOC( FillPhysBCCrse ), &
             C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             pre_interp_ptr, post_interp_ptr)

  END SUBROUTINE amrex_fillcoarsepatch_default


  SUBROUTINE amrex_fillcoarsepatch_dgconservative &
    ( MF, MF_G, &
      OldTimeCrse, OldMFCrse, OldMFCrse_G, &
      NewTimeCrse, NewMFCrse, NewMFCrse_G, &
      CrseGeom, FillPhysBCCrse, FineGeom, FillPhysBCFine, &
      Time, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF, MF_G
    TYPE(amrex_multifab), INTENT(in)    :: OldMFCrse, OldMFCrse_G, &
                                           NewMFCrse, NewMFCrse_G
    INTEGER             , INTENT(in)    :: nComp, RefRatio, &
                                           interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,nComp), HiBC(amrex_spacedim,nComp)
    REAL(amrex_real)    , INTENT(in)    :: OldTimeCrse, NewTimeCrse, Time
    TYPE(amrex_geometry), INTENT(in)    :: CrseGeom, FineGeom
    PROCEDURE(amrex_physbc_proc)        :: FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: pre_interp
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: post_interp
    TYPE(c_ptr)         , INTENT(in)    :: vpCoarseToFineProjectionMatrix

    REAL(amrex_real) :: teps
    TYPE(c_ptr)      :: pCrseMF, pCrseMF_G
    TYPE(c_ptr)      :: pLoBC(nComp), pHiBC(nComp)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: iComp

    INTEGER :: sComp, dComp

    sComp = 1
    dComp = 1

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       pCrseMF   = NewMFCrse   % p
       pCrseMF_G = NewMFCrse_G % p

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       pCrseMF   = OldMFCrse   % p
       pCrseMF_G = OldMFCrse_G % p

    ELSE

       pCrseMF = NewMFCrse % p

       CALL amrex_abort( "amrex_fillcoarsepatch: how did this happen?" )

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillcoarsepatch_dgconservative &
           ( MF % p, MF_G % p, Time, pCrseMF, pCrseMF_G, &
             sComp-1, dComp-1, nComp, &
             CrseGeom % p, FineGeom % p, &
             C_FUNLOC( FillPhysBCCrse ), &
             C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr)

  END SUBROUTINE amrex_fillcoarsepatch_dgconservative


  SUBROUTINE amrex_fillcoarsepatch_dgpointwise &
    ( MF, &
      OldTimeCrse, OldMFCrse, &
      NewTimeCrse, NewMFCrse, &
      CrseGeom, FillPhysBCCrse, FineGeom, FillPhysBCFine, &
      Time, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: OldMFCrse, &
                                           NewMFCrse
    INTEGER             , INTENT(in)    :: nComp, RefRatio, &
                                           interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,nComp), HiBC(amrex_spacedim,nComp)
    REAL(amrex_real)    , INTENT(in)    :: OldTimeCrse, NewTimeCrse, Time
    TYPE(amrex_geometry), INTENT(in)    :: CrseGeom, FineGeom
    PROCEDURE(amrex_physbc_proc)        :: FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: pre_interp
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: post_interp
    TYPE(c_ptr)         , INTENT(in)    :: vpCoarseToFineProjectionMatrix

    REAL(amrex_real) :: teps
    TYPE(c_ptr)      :: pCrseMF
    TYPE(c_ptr)      :: pLoBC(nComp), pHiBC(nComp)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: iComp

    INTEGER :: sComp, dComp

    sComp = 1
    dComp = 1

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       pCrseMF = NewMFCrse % p

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       pCrseMF = OldMFCrse % p

    ELSE

       pCrseMF = NewMFCrse % p

       CALL amrex_abort( "amrex_fillcoarsepatch: how did this happen?" )

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillcoarsepatch_dgpointwise &
           ( MF % p, Time, pCrseMF, &
             sComp-1, dComp-1, nComp, &
             CrseGeom % p, FineGeom % p, &
             C_FUNLOC( FillPhysBCCrse ), &
             C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr)

  END SUBROUTINE amrex_fillcoarsepatch_dgpointwise


  subroutine amrex_fillcoarsepatch_faces (mf, told_c, mfold_c, tnew_c, mfnew_c, &
       &                                  geom_c, fill_physbc_cx, &
#if (AMREX_SPACEDIM > 1)
       &                                  fill_physbc_cy, &
#if (AMREX_SPACEDIM > 2)
       &                                  fill_physbc_cz, &
#endif
#endif
       &                                  geom_f, fill_physbc_fx, &
#if (AMREX_SPACEDIM > 1)
       &                                  fill_physbc_fy, &
#if (AMREX_SPACEDIM > 2)
       &                                  fill_physbc_fz, &
#endif
#endif
       &                                  time, scomp, dcomp, ncomp, rr, interp, lo_bc, hi_bc, &
       &                                  pre_interp, post_interp)
    type(amrex_multifab), intent(inout) :: mf(amrex_spacedim)
    type(amrex_multifab), intent(in   ) :: mfold_c(amrex_spacedim), mfnew_c(amrex_spacedim)
    integer, intent(in) :: scomp, dcomp, ncomp, rr, interp
    !                 (BC dir        , comp        , MF)
    integer, dimension(amrex_spacedim,scomp+ncomp-1,amrex_spacedim), target, intent(in) :: lo_bc, hi_bc
    real(amrex_real), intent(in) :: told_c, tnew_c, time
    type(amrex_geometry), intent(in) :: geom_c, geom_f
    procedure(amrex_physbc_proc) :: fill_physbc_cx, fill_physbc_fx
#if (AMREX_SPACEDIM > 1)
    procedure(amrex_physbc_proc) :: fill_physbc_cy, fill_physbc_fy
#endif
#if (AMREX_SPACEDIM > 2)
    procedure(amrex_physbc_proc) :: fill_physbc_cz, fill_physbc_fz
#endif
    procedure(amrex_interp_hook_arr_proc), optional :: pre_interp
    procedure(amrex_interp_hook_arr_proc), optional :: post_interp

    real(amrex_real) :: teps
    type(c_ptr) :: faces(amrex_spacedim)
    type(c_ptr) :: c_mf(amrex_spacedim)
    type(c_funptr) :: cfill(amrex_spacedim), ffill(amrex_spacedim)
    type(c_ptr) :: lo_bc_ptr(amrex_spacedim*(scomp+ncomp-1)), hi_bc_ptr(amrex_spacedim*(scomp+ncomp-1))
    type(c_funptr) :: pre_interp_ptr, post_interp_ptr
    integer :: i, nc, dim, mfid

    cfill(1) = c_funloc(fill_physbc_cx)
    ffill(1) = c_funloc(fill_physbc_fx)
#if (AMREX_SPACEDIM >= 2)
    cfill(2) = c_funloc(fill_physbc_cy)
    ffill(2) = c_funloc(fill_physbc_fy)
#if (AMREX_SPACEDIM >= 3)
    cfill(3) = c_funloc(fill_physbc_cz)
    ffill(3) = c_funloc(fill_physbc_fz)
#endif
#endif

    do dim = 1, amrex_spacedim
       faces(dim) = mf(dim)%p
    end do

    ! coarse level
    teps = 1.e-4_amrex_real * abs(tnew_c - told_c)
    if (abs(time-tnew_c) .le. teps) then
       do dim = 1, amrex_spacedim
          c_mf(dim) = mfnew_c(dim)%p
       end do
    else if (abs(time-told_c) .le. teps) then
       do dim = 1, amrex_spacedim
          c_mf(dim) = mfold_c(dim)%p
       end do
    else
       call amrex_abort("amrex_fillcoarsepatch_faces: how did this happen?")
    end if

    ! lo_bc & hi_bc: (BC dir, comp, MF)
    nc = scomp+ncomp-1
    do mfid = 1, amrex_spacedim
       do i = 1, scomp-1
          lo_bc_ptr((mfid-1)*nc + i) = c_null_ptr
          hi_bc_ptr((mfid-1)*nc + i) = c_null_ptr
       end do
       do i = scomp, nc
          lo_bc_ptr((mfid-1)*nc + i) = c_loc(lo_bc(1,i,mfid))
          hi_bc_ptr((mfid-1)*nc + i) = c_loc(hi_bc(1,i,mfid))
       end do
    end do

    pre_interp_ptr = c_null_funptr
    if (present(pre_interp)) pre_interp_ptr = c_funloc(pre_interp)
    post_interp_ptr = c_null_funptr
    if (present(post_interp)) post_interp_ptr = c_funloc(post_interp)

    ! scomp-1 and dcomp-1 because of Fortran index starts with 1
    call amrex_fi_fillcoarsepatch_faces(faces, time, c_mf, scomp-1, dcomp-1, ncomp, &
         &                              geom_c%p, geom_f%p,              &
         &                              cfill, ffill,                    &
         &                              rr, interp, lo_bc_ptr, hi_bc_ptr,&
         &                              pre_interp_ptr, post_interp_ptr)
  end subroutine amrex_fillcoarsepatch_faces

end module amrex_fillpatch_module
