set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND ALLHEADERS AMReX_AggStencil.H         AMReX_BaseIVFABI.H       AMReX_EBFluxFactory.H )
list ( APPEND ALLHEADERS AMReX_FlatPlateGeom.H      AMReX_PlaneIF.H )
list ( APPEND ALLHEADERS AMReX_AggStencilI.H        AMReX_BaseIVFactory.H    AMReX_EBGraph.H )
list ( APPEND ALLHEADERS AMReX_GeometryService.H    AMReX_PolyGeom.H )
list ( APPEND ALLHEADERS AMReX_AllRegularService.H  AMReX_ComplementIF.H     AMReX_EBIndexSpace.H )
list ( APPEND ALLHEADERS AMReX_GeometryShop.H       AMReX_PolynomialIF.H )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFAB.H      AMReX_EBArith.H )        AMReX_EBISBox.H )
list ( APPEND ALLHEADERS AMReX_GraphNode.H          AMReX_RedistStencil.H )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFABI.H     AMReX_EBCellFAB.H        AMReX_EBISLayout.H )
list ( APPEND ALLHEADERS AMReX_IntersectionIF.H     AMReX_SlabService.H )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFactory.H  AMReX_EBCellFactory.H    AMReX_EBISLevel.H )
list ( APPEND ALLHEADERS AMReX_IntVectSet.H         AMReX_SphereIF.H )
list ( APPEND ALLHEADERS AMReX_BaseEBFaceFAB.H      AMReX_EBDataFactory.H    AMReX_EBLevelGrid.H )
list ( APPEND ALLHEADERS AMReX_IrregFABFactory.H    AMReX_Stencils.H )
list ( APPEND ALLHEADERS AMReX_BaseEBFaceFABI.H     AMReX_EBData.H           AMReX_EBLevelRedist.H )
list ( APPEND ALLHEADERS AMReX_IrregFAB.H           AMReX_TransformIF.H )
list ( APPEND ALLHEADERS AMReX_BaseIFFAB.H          AMReX_EBDataVarMacros.H  AMReX_EBLoHiCenter.H )
list ( APPEND ALLHEADERS AMReX_IrregNode.H          AMReX_UnionIF.H )
list ( APPEND ALLHEADERS AMReX_BaseIFFABI.H         AMReX_EBDebugDump.H      AMReX_EBNormalizeByVolumeFraction.H )
list ( APPEND ALLHEADERS AMReX_LatheIF.H            AMReX_VoFIterator.H )
list ( APPEND ALLHEADERS AMReX_BaseIFFactory.H      AMReX_EBDebugOut.H       AMReX_FabArrayIO.H )
list ( APPEND ALLHEADERS AMReX_LoHiSide.H           AMReX_VolIndex.H )
list ( APPEND ALLHEADERS AMReX_BaseIF.H             AMReX_EBFaceFAB.H        AMReX_FaceIndex.H )
list ( APPEND ALLHEADERS AMReX_LSquares.H )
list ( APPEND ALLHEADERS AMReX_BaseIVFAB.H          AMReX_EBFluxFAB.H        AMReX_FaceIterator.H )
list ( APPEND ALLHEADERS AMReX_Moments.H )

list ( APPEND CXXSRC AMReX_AllRegularService.cpp  AMReX_EBGraph.cpp                       )
list ( APPEND CXXSRC AMReX_FabArrayIO.cpp         AMReX_IrregFAB.cpp         AMReX_PolynomialIF.cpp )
list ( APPEND CXXSRC AMReX_ComplementIF.cpp       AMReX_EBIndexSpace.cpp                  )
list ( APPEND CXXSRC AMReX_FaceIndex.cpp          AMReX_IrregFABFactory.cpp  AMReX_RedistStencil.cpp )
list ( APPEND CXXSRC AMReX_EBArith.cpp            AMReX_EBISBox.cpp                       )
list ( APPEND CXXSRC AMReX_FaceIterator.cpp       AMReX_IrregNode.cpp        AMReX_SlabService.cpp )
list ( APPEND CXXSRC AMReX_EBCellFAB.cpp          AMReX_EBISLayout.cpp                    )
list ( APPEND CXXSRC AMReX_FlatPlateGeom.cpp      AMReX_LatheIF.cpp          AMReX_SphereIF.cpp )
list ( APPEND CXXSRC AMReX_EBData.cpp             AMReX_EBISLevel.cpp                     )
list ( APPEND CXXSRC AMReX_GeometryService.cpp    AMReX_LoHiSide.cpp         AMReX_Stencils.cpp )
list ( APPEND CXXSRC AMReX_EBDebugOut.cpp         AMReX_EBLevelGrid.cpp                   )
list ( APPEND CXXSRC AMReX_GeometryShop.cpp       AMReX_LSquares.cpp         AMReX_TransformIF.cpp )
list ( APPEND CXXSRC AMReX_EBFaceFAB.cpp          AMReX_EBLevelRedist.cpp                 )
list ( APPEND CXXSRC AMReX_GraphNode.cpp          AMReX_Moments.cpp          AMReX_UnionIF.cpp )
list ( APPEND CXXSRC AMReX_EBFluxFAB.cpp          AMReX_EBLoHiCenter.cpp                  )
list ( APPEND CXXSRC AMReX_IntersectionIF.cpp     AMReX_PlaneIF.cpp          AMReX_VoFIterator.cpp )
list ( APPEND CXXSRC AMReX_EBFluxFactory.cpp      AMReX_EBNormalizeByVolumeFraction.cpp   )
list ( APPEND CXXSRC AMReX_IntVectSet.cpp         AMReX_PolyGeom.cpp         AMReX_VolIndex.cpp )

#
# This is an option to display arrays in GDB. It is for debugging
# AMReX and currently not supported by cmake.# 
# if(USE_ARRAYVIEW EQUAL TRUE)
#   set(CXX_source_files ${CXX_source_files} AMReX_DatasetClient.cpp)
#   set(CXX_header_files ${CXX_header_files} AMReX_DatasetClient.H)
# endif()

# 
#  Collect sources
# 
set ( ALLSRC  ${CXXSRC} ${F90SRC} ${F77SRC} )

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
