target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE
   AMReX_EBFArrayBox.H
   AMReX_EBFArrayBox.cpp
   AMReX_EBFabFactory.H
   AMReX_EBFabFactory.cpp
   AMReX_EBMultiFabUtil.H
   AMReX_EBMultiFabUtil_C.H
   AMReX_EBMultiFabUtil_${DIM}D_C.H
   AMReX_EBMultiFabUtil.cpp
   AMReX_EBCellFlag.H
   AMReX_EBCellFlag.cpp
   AMReX_EBDataCollection.H
   AMReX_EBDataCollection.cpp
   AMReX_MultiCutFab.H
   AMReX_MultiCutFab.cpp
   AMReX_EBSupport.H
   AMReX_EBInterpolater.H
   AMReX_EBInterpolater.cpp
   AMReX_EBFluxRegister.H
   AMReX_EBFluxRegister_C.H
   AMReX_EBFluxRegister_${DIM}D_C.H
   AMReX_EBFluxRegister.cpp
   AMReX_EBAmrUtil.H
   AMReX_EBAmrUtil.cpp
   AMReX_EB_utils.H
   AMReX_EB_utils.cpp
   AMReX_algoim.H
   AMReX_algoim_K.H
   AMReX_algoim.cpp
   AMReX_EB2_IF_AllRegular.H
   AMReX_EB2_IF_Box.H
   AMReX_EB2_IF_Cylinder.H
   AMReX_EB2_IF_Ellipsoid.H
   AMReX_EB2_IF_Plane.H
   AMReX_EB2_IF_Sphere.H
   AMReX_EB2_IF_Torus.H
   AMReX_distFcnElement.H
   AMReX_EB2_IF_Spline.H
   AMReX_EB2_IF_Polynomial.H
   AMReX_EB2_IF_Complement.H
   AMReX_EB2_IF_Intersection.H
   AMReX_EB2_IF_Lathe.H
   AMReX_EB2_IF_Rotation.H
   AMReX_EB2_IF_Scale.H
   AMReX_EB2_IF_Translation.H
   AMReX_EB2_IF_Union.H
   AMReX_EB2_IF_Extrusion.H
   AMReX_EB2_IF_Difference.H
   AMReX_EB2_IF.H
   AMReX_EB2_IF_Base.H
   AMReX_distFcnElement.cpp
   AMReX_EB2_GeometryShop.H
   AMReX_EB2.H
   AMReX_EB2_IndexSpaceI.H
   AMReX_EB2_Level.H
   AMReX_EB2_Graph.H
   AMReX_EB2_MultiGFab.H
   AMReX_EB2.cpp
   AMReX_EB2_Level.cpp
   AMReX_EB2_MultiGFab.cpp
   AMReX_EB2_${DIM}D_C.cpp
   AMReX_EB2_C.H AMReX_EB2_${DIM}D_C.H
   )

if (ENABLE_FORTRAN)
   target_sources(amrex
      PRIVATE
      AMReX_ebcellflag_mod.F90
      AMReX_EBFluxRegister_nd.F90
      AMReX_compute_normals.F90
      AMReX_EB_F.H
      AMReX_WriteEBSurface.cpp
      AMReX_WriteEBSurface.H
      AMReX_WriteEB_F.H
      AMReX_eb_to_pvd.F90
      # LS
      AMReX_EB_geometry.F90
      AMReX_EB_bc_fill_nd.F90
      AMReX_EB_Tagging.F90
      AMReX_EB_levelset_F.F90
      AMReX_EB_levelset.H
      AMReX_EB_levelset.cpp
      )

   if (DIM EQUAL 3)
      target_sources(amrex
         PRIVATE
         AMReX_EB_LSCoreBase.H
         AMReX_EB_LSCoreBase.cpp
         AMReX_EB_LSCore.H
         AMReX_EB_LSCoreI.H
         AMReX_EB_LSCore_F.H
         )
  endif ()

endif ()
