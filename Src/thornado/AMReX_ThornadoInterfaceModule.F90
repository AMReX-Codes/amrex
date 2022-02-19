MODULE AMReX_ThornadoInterfaceModule

  USE ISO_C_BINDING

  USE amrex_fort_module, ONLY: &
    DP => amrex_real

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex_InitializeMeshRefinement_Thornado
  PUBLIC :: amrex_FinalizeMeshRefinement_Thornado

  INTERFACE ThornadoInterface

    SUBROUTINE amrex_fi_initializemeshrefinement_thornado &
      ( nNodes, ProjectionMatrix, WeightsX_q ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int), INTENT(in) :: nNodes(*)
       REAL(DP)      , INTENT(in) :: ProjectionMatrix(*)
       REAL(DP)      , INTENT(in) :: WeightsX_q(*)
    END SUBROUTINE amrex_fi_initializemeshrefinement_thornado

    SUBROUTINE amrex_fi_finalizemeshrefinement_thornado() BIND(c)
       IMPORT
       IMPLICIT NONE
    END SUBROUTINE amrex_fi_finalizemeshrefinement_thornado

  END INTERFACE ThornadoInterface


CONTAINS


  SUBROUTINE amrex_InitializeMeshRefinement_Thornado &
    ( nNodes, ProjectionMatrix, WeightsX_q )
    INTEGER , INTENT(in) :: nNodes(*)
    REAL(DP), INTENT(in) :: ProjectionMatrix(*)
    REAL(DP), INTENT(in) :: WeightsX_q(*)

    CALL amrex_fi_initializemeshrefinement_thornado &
           ( nNodes, ProjectionMatrix, WeightsX_q )
  END SUBROUTINE amrex_InitializeMeshRefinement_Thornado


  SUBROUTINE amrex_FinalizeMeshRefinement_Thornado
    CALL amrex_fi_finalizemeshrefinement_thornado
  END SUBROUTINE amrex_FinalizeMeshRefinement_Thornado


END MODULE AMReX_ThornadoInterfaceModule
