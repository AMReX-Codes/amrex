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
        LX_X3_Refined_Packed, &
        LX_X1_Up_1D, LX_X1_Dn_1D, &
        LX_X2_Up_1D, LX_X2_Dn_1D, &
        LX_X3_Up_1D, LX_X3_Dn_1D, iGF_SqrtGm ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int), INTENT(in) :: nNodes(*)
       REAL(DP)      , INTENT(in) :: ProjectionMatrix(*)
       REAL(DP)      , INTENT(in) :: WeightsX1(*), WeightsX2(*), WeightsX3(*)
       REAL(DP)      , INTENT(in) :: LX_X1_Refined_Packed(*), &
                                     LX_X2_Refined_Packed(*), &
                                     LX_X3_Refined_Packed(*), &
                                     LX_X1_Up_1D(*), LX_X1_Dn_1D(*), &
                                     LX_X2_Up_1D(*), LX_X2_Dn_1D(*), &
                                     LX_X3_Up_1D(*), LX_X3_Dn_1D(*)
       INTEGER(c_int), INTENT(in), VALUE :: iGF_SqrtGm
    END SUBROUTINE amrex_fi_initializemeshrefinement_dg

    SUBROUTINE amrex_fi_finalizemeshrefinement_dg() BIND(c)
       IMPORT
       IMPLICIT NONE
    END SUBROUTINE amrex_fi_finalizemeshrefinement_dg

  END INTERFACE DGInterface


CONTAINS


  SUBROUTINE amrex_InitializeMeshRefinement_DG &
    ( nNodes, ProjectionMatrix, WeightsX1, WeightsX2, WeightsX3, &
      LX_X1_Refined_Packed, LX_X2_Refined_Packed, LX_X3_Refined_Packed, &
      LX_X1_Up_1D, LX_X1_Dn_1D, &
      LX_X2_Up_1D, LX_X2_Dn_1D, &
      LX_X3_Up_1D, LX_X3_Dn_1D, iGF_SqrtGm )
    INTEGER , INTENT(in) :: nNodes(*), iGF_SqrtGm
    REAL(DP), INTENT(in) :: ProjectionMatrix(*)
    REAL(DP), INTENT(in) :: WeightsX1(*), WeightsX2(*), WeightsX3(*)
    REAL(DP), INTENT(in) :: LX_X1_Refined_Packed(*), &
                            LX_X2_Refined_Packed(*), &
                            LX_X3_Refined_Packed(*), &
                            LX_X1_Up_1D(*), LX_X1_Dn_1D(*), &
                            LX_X2_Up_1D(*), LX_X2_Dn_1D(*), &
                            LX_X3_Up_1D(*), LX_X3_Dn_1D(*)

    CALL amrex_fi_initializemeshrefinement_dg &
           ( nNodes, ProjectionMatrix, &
             WeightsX1, WeightsX2, WeightsX3, &
             LX_X1_Refined_Packed, &
             LX_X2_Refined_Packed, &
             LX_X3_Refined_Packed, &
             LX_X1_Up_1D, LX_X1_Dn_1D, &
             LX_X2_Up_1D, LX_X2_Dn_1D, &
             LX_X3_Up_1D, LX_X3_Dn_1D, iGF_SqrtGm  )
  END SUBROUTINE amrex_InitializeMeshRefinement_DG


  SUBROUTINE amrex_FinalizeMeshRefinement_DG
    CALL amrex_fi_finalizemeshrefinement_dg
  END SUBROUTINE amrex_FinalizeMeshRefinement_DG


END MODULE AMReX_DGInterfaceModule
