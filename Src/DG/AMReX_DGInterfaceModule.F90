MODULE AMReX_DGInterfaceModule

  USE ISO_C_BINDING

  USE amrex_fort_module, ONLY: &
    DP => amrex_real

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex_InitializeMeshRefinement_DG
  PUBLIC :: amrex_FinalizeMeshRefinement_DG

  INTERFACE DGInterface

    SUBROUTINE amrex_fi_initializemeshrefinement_dg &
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
    END SUBROUTINE amrex_fi_initializemeshrefinement_dg

    SUBROUTINE amrex_fi_finalizemeshrefinement_dg() BIND(c)
       IMPORT
       IMPLICIT NONE
    END SUBROUTINE amrex_fi_finalizemeshrefinement_dg

  END INTERFACE DGInterface


CONTAINS


  SUBROUTINE amrex_InitializeMeshRefinement_DG &
    ( nNodes, ProjectionMatrix, WeightsX1, WeightsX2, WeightsX3, &
      LX_X1_Refined_Packed, LX_X2_Refined_Packed, LX_X3_Refined_Packed )
    INTEGER , INTENT(in) :: nNodes(*)
    REAL(DP), INTENT(in) :: ProjectionMatrix(*)
    REAL(DP), INTENT(in) :: WeightsX1(*), WeightsX2(*), WeightsX3(*)
    REAL(DP), INTENT(in) :: LX_X1_Refined_Packed(*), &
                            LX_X2_Refined_Packed(*), &
                            LX_X3_Refined_Packed(*)

    CALL amrex_fi_initializemeshrefinement_dg &
           ( nNodes, ProjectionMatrix, &
             WeightsX1, WeightsX2, WeightsX3, &
             LX_X1_Refined_Packed, &
             LX_X2_Refined_Packed, &
             LX_X3_Refined_Packed )
  END SUBROUTINE amrex_InitializeMeshRefinement_DG


  SUBROUTINE amrex_FinalizeMeshRefinement_DG
    CALL amrex_fi_finalizemeshrefinement_dg
  END SUBROUTINE amrex_FinalizeMeshRefinement_DG


END MODULE AMReX_DGInterfaceModule
