set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS AMReX_BaseIF.H  AMReX_GeometryShop.H   AMReX_Moments.H AMReX_GeomIntersectUtils.H )
list ( APPEND CXXSRC     AMReX_GeometryShop.cpp AMReX_Moments.cpp AMReX_GeomIntersectUtils.cpp )

list ( APPEND ALLHEADERS AMReX_LSquares.H   AMReX_SphereIF.H    AMReX_AnisotropicIF.H	AMReX_PlaneIF.H	AMReX_AnisotropicDxPlaneIF.H )
list ( APPEND ALLHEADERS AMReX_IrregNode.H   AMReX_LoHiSide.H )
   
list ( APPEND CXXSRC     AMReX_LSquares.cpp AMReX_SphereIF.cpp  AMReX_PlaneIF.cpp )
list ( APPEND CXXSRC     AMReX_IrregNode.cpp AMReX_LoHiSide.cpp )

list ( APPEND ALLHEADERS AMReX_IntVectSet.H     AMReX_FlatPlateGeom.H    AMReX_EBCellFAB.H    AMReX_FaceIndex.H     AMReX_VolIndex.H )
list ( APPEND CXXSRC     AMReX_IntVectSet.cpp  AMReX_FlatPlateGeom.cpp  AMReX_EBCellFAB.cpp  AMReX_FaceIndex.cpp   AMReX_VolIndex.cpp )

list ( APPEND ALLHEADERS AMReX_Stencils.H   AMReX_GraphNode.H     AMReX_EBGraph.H )  
list ( APPEND CXXSRC     AMReX_Stencils.cpp AMReX_GraphNode.cpp   AMReX_EBGraph.cpp )

list ( APPEND ALLHEADERS AMReX_VoFIterator.H   AMReX_FaceIterator.H )
list ( APPEND CXXSRC     AMReX_VoFIterator.cpp AMReX_FaceIterator.cpp )

list ( APPEND ALLHEADERS AMReX_BaseIVFAB.H AMReX_BaseIVFABI.H AMReX_BaseIVFactory.H )
list ( APPEND ALLHEADERS AMReX_BaseIFFAB.H AMReX_BaseIFFABI.H AMReX_EBData.H  AMReX_EBDataFactory.H )
list ( APPEND ALLHEADERS AMReX_BaseEBFaceFABI.H )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFAB.H )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFABI.H )
list ( APPEND ALLHEADERS AMReX_EBCellFactory.H )

list ( APPEND CXXSRC     AMReX_EBData.cpp )

list ( APPEND ALLHEADERS AMReX_EBISBox.H   AMReX_PolyGeom.H )
list ( APPEND CXXSRC     AMReX_EBISBox.cpp AMReX_PolyGeom.cpp )

list ( APPEND ALLHEADERS AMReX_EBIndexSpace.H   AMReX_EBISLayout.H   AMReX_EBISLevel.H   AMReX_GeometryService.H )
list ( APPEND CXXSRC     AMReX_EBIndexSpace.cpp AMReX_EBISLayout.cpp AMReX_EBISLevel.cpp AMReX_GeometryService.cpp )


list ( APPEND ALLHEADERS AMReX_EBLevelGrid.H  AMReX_AllRegularService.H )
list ( APPEND CXXSRC     AMReX_EBLevelGrid.cpp  AMReX_AllRegularService.cpp )

list ( APPEND ALLHEADERS AMReX_EBDebugOut.H AMReX_EBDebugDump.H AMReX_EBFaceFAB.H AMReX_BaseEBFaceFAB.H AMReX_EBArith.H )
list ( APPEND CXXSRC     AMReX_EBDebugOut.cpp AMReX_EBFaceFAB.cpp AMReX_EBArith.cpp )

list ( APPEND ALLHEADERS AMReX_RedistStencil.H  AMReX_EBLevelRedist.H AMReX_EBFluxFAB.H AMReX_IrregFAB.H )
list ( APPEND ALLHEADERS AMReX_EBFluxFactory.H   AMReX_IrregFABFactory.H )
list ( APPEND CXXSRC     AMReX_RedistStencil.cpp AMReX_EBLevelRedist.cpp AMReX_EBFluxFAB.cpp )
list ( APPEND CXXSRC     AMReX_IrregFAB.cpp AMReX_EBFluxFactory.cpp AMReX_IrregFABFactory.cpp )

list ( APPEND ALLHEADERS AMReX_EBNormalizeByVolumeFraction.H   AMReX_EBLoHiCenter.H   AMReX_FabArrayIO.H )
list ( APPEND CXXSRC     AMReX_EBNormalizeByVolumeFraction.cpp AMReX_EBLoHiCenter.cpp AMReX_FabArrayIO.cpp )

list ( APPEND ALLHEADERS AMReX_UnionIF.H   AMReX_IntersectionIF.H   AMReX_ComplementIF.H   AMReX_LatheIF.H   AMReX_TransformIF.H AMReX_ExtrudeIF.H )
list ( APPEND CXXSRC     AMReX_UnionIF.cpp AMReX_IntersectionIF.cpp AMReX_ComplementIF.cpp AMReX_LatheIF.cpp AMReX_TransformIF.cpp AMReX_ExtrudeIF.cpp )

list ( APPEND ALLHEADERS AMReX_PolynomialIF.H   )
list ( APPEND CXXSRC     AMReX_PolynomialIF.cpp )

list ( APPEND ALLHEADERS AMReX_STLAsciiReader.H    AMReX_STLExplorer.H )
list ( APPEND ALLHEADERS AMReX_STLMesh.H   AMReX_STLBox.H )
list ( APPEND ALLHEADERS AMReX_STLIF.H	     AMReX_STLUtil.H   AMReX_STLReader.H )

list ( APPEND CXXSRC     AMReX_STLAsciiReader.cpp  AMReX_STLExplorer.cpp )
list ( APPEND CXXSRC     AMReX_STLMesh.cpp AMReX_STLBox.cpp	  AMReX_STLIF.cpp    AMReX_STLUtil.cpp )

list ( APPEND ALLHEADERS AMReX_CellEdge.H   AMReX_KDTree.H   AMReX_PXStuff.H )
list ( APPEND CXXSRC     AMReX_CellEdge.cpp AMReX_KDTree.cpp AMReX_PXStuff.cpp )

# 
#  Collect sources
# 
set ( ALLSRC  ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR} )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})
