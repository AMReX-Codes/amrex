add_amrex_define( AMREX_USE_EB NO_LEGACY IF AMReX_EB )

target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE
   AMReX_EBFArrayBox.H
   AMReX_EBFArrayBox.cpp
   AMReX_EBFabFactory.H
   AMReX_EBFabFactory.cpp
   AMReX_EBMultiFabUtil.H
   AMReX_EBMultiFabUtil_C.H
   AMReX_EBMultiFabUtil_${AMReX_SPACEDIM}D_C.H
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
   AMReX_EBMFInterpolater.H
   AMReX_EBMFInterpolater.cpp
   AMReX_EBFluxRegister.H
   AMReX_EBFluxRegister_C.H
   AMReX_EBFluxRegister_${AMReX_SPACEDIM}D_C.H
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
   AMReX_EB2_IF_Parser.H
   AMReX_EB2_IF.H
   AMReX_EB2_IF_Base.H
   AMReX_EB2_IF_DevicePtr.H
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
   AMReX_EB2_${AMReX_SPACEDIM}D_C.cpp
   AMReX_EB2_ND_C.cpp
   AMReX_EB2_C.H
   AMReX_EB2_${AMReX_SPACEDIM}D_C.H
   AMReX_EB_STL_utils.H
   AMReX_EB_STL_utils.cpp
   AMReX_EB_chkpt_file.H
   AMReX_EB_chkpt_file.cpp
   AMReX_EB_triGeomOps_K.H
   AMReX_EB2_Level_STL.H
   AMReX_EB2_Level_STL.cpp
   AMReX_EB2_IndexSpace_STL.H
   AMReX_EB2_IndexSpace_STL.cpp
   AMReX_EB2_Level_chkpt_file.H
   AMReX_EB2_Level_chkpt_file.cpp
   AMReX_EB2_IndexSpace_chkpt_file.H
   AMReX_EB2_IndexSpace_chkpt_file.cpp
   )

if (AMReX_SPACEDIM EQUAL 3)
   target_sources(amrex
   PRIVATE
   AMReX_WriteEBSurface.cpp
   AMReX_WriteEBSurface.H
   AMReX_EBToPVD.H
   AMReX_EBToPVD.cpp
   )
endif ()

if (AMReX_FORTRAN)
   target_sources(amrex
      PRIVATE
      AMReX_ebcellflag_mod.F90
      AMReX_EBFluxRegister_nd.F90
      )
endif ()
