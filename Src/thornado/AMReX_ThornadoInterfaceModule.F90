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
      ( nNodes, ProjectionMatrix, WeightsX1, WeightsX2, WeightsX3, &
        LX_X1_Refined_Packed, &
        LX_X2_Refined_Packed, &
        LX_X3_Refined_Packed ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int), INTENT(in) :: nNodes(*)
       REAL(DP)      , INTENT(in) :: ProjectionMatrix(*)
       REAL(DP)      , INTENT(in) :: WeightsX1(*), WeightsX2(*), WeightsX3(*)
       REAL(DP)      , INTENT(in) :: LX_X1_Refined_Packed(*), &
                                     LX_X2_Refined_Packed(*), &
                                     LX_X3_Refined_Packed(*)
    END SUBROUTINE amrex_fi_initializemeshrefinement_thornado

    SUBROUTINE amrex_fi_finalizemeshrefinement_thornado() BIND(c)
       IMPORT
       IMPLICIT NONE
    END SUBROUTINE amrex_fi_finalizemeshrefinement_thornado

  END INTERFACE ThornadoInterface


CONTAINS


  SUBROUTINE amrex_InitializeMeshRefinement_Thornado &
    ( nNodes, ProjectionMatrix, WeightsX1, WeightsX2, WeightsX3, &
      LX_X1_Refined_Packed, LX_X2_Refined_Packed, LX_X3_Refined_Packed )
    INTEGER , INTENT(in) :: nNodes(*)
    REAL(DP), INTENT(in) :: ProjectionMatrix(*)
    REAL(DP), INTENT(in) :: WeightsX1(*), WeightsX2(*), WeightsX3(*)
    REAL(DP), INTENT(in) :: LX_X1_Refined_Packed(*), &
                            LX_X2_Refined_Packed(*), &
                            LX_X3_Refined_Packed(*)

    CALL amrex_fi_initializemeshrefinement_thornado &
           ( nNodes, ProjectionMatrix, &
             WeightsX1, WeightsX2, WeightsX3, &
             LX_X1_Refined_Packed, &
             LX_X2_Refined_Packed, &
             LX_X3_Refined_Packed )
  END SUBROUTINE amrex_InitializeMeshRefinement_Thornado


  SUBROUTINE amrex_FinalizeMeshRefinement_Thornado
    CALL amrex_fi_finalizemeshrefinement_thornado
  END SUBROUTINE amrex_FinalizeMeshRefinement_Thornado


END MODULE AMReX_ThornadoInterfaceModule
