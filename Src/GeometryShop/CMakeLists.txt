set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND ALLHEADERS AMReX_AggStencil.H	      AMReX_CoordinateSystem.H		 AMReX_EBIndexSpace.H		      AMReX_IFData.H	       AMReX_MinimalCCCMImplem.H     AMReX_STLIF.H               )
list ( APPEND ALLHEADERS AMReX_AggStencilI.H	      AMReX_CoordinateSystemImplem.H	 AMReX_EBLevelGrid.H		      AMReX_IFDataImplem.H     AMReX_MomentIterator.H	     AMReX_STLMesh.H             )
list ( APPEND ALLHEADERS AMReX_AllRegularService.H    AMReX_CutCellMoments.H		 AMReX_EBLevelRedist.H		      AMReX_IFSlicer.H	       AMReX_MomentIteratorImplem.H  AMReX_STLReader.H           )
list ( APPEND ALLHEADERS AMReX_AnisotropicDxPlaneIF.H AMReX_CutCellMomentsImplem.H	 AMReX_EBLoHiCenter.H		      AMReX_IFSlicerImplem.H   AMReX_Moments.H		     AMReX_STLUtil.H             )
list ( APPEND ALLHEADERS AMReX_AnisotropicIF.H	      AMReX_DivNormalRefinement.H	 AMReX_EBNormalizeByVolumeFraction.H  AMReX_IndexTM.H	       AMReX_MonomialPowers.H	              )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFAB.H	      AMReX_DivNormalRefinementImplem.H  AMReX_EB_TYPEDEFS.H		      AMReX_IndexTMI.H	       AMReX_MultiIndex.H	     AMReX_SmoothAbsoluteValue.H )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFABI.H	      AMReX_EBArith.H			 AMReX_EllipsoidIF.H		      AMReX_IndexedMoments.H   AMReX_MultiIndexImplem.H      AMReX_SmoothIntersection.H  )
list ( APPEND ALLHEADERS AMReX_BaseEBCellFactory.H    AMReX_EBCellFAB.H			 AMReX_ExtrudeIF.H		      AMReX_IntVectSet.H       AMReX_NoRefinement.H	     AMReX_SmoothUnion.H         )
list ( APPEND ALLHEADERS AMReX_BaseEBFaceFAB.H	      AMReX_EBCellFactory.H		 AMReX_FabArrayIO.H		      AMReX_IntersectionIF.H   AMReX_NoRefinementImplem.H    AMReX_SphereIF.H            )
list ( APPEND ALLHEADERS AMReX_BaseEBFaceFABI.H	      AMReX_EBData.H			 AMReX_FaceIndex.H		      AMReX_IrregFAB.H	       AMReX_NormalDerivative.H      AMReX_Stencils.H            )
list ( APPEND ALLHEADERS AMReX_BaseIF.H		      AMReX_EBDataFactory.H		 AMReX_FaceIterator.H		      AMReX_IrregFABFactory.H  AMReX_NormalDerivativeNew.H   AMReX_TransformIF.H         )
list ( APPEND ALLHEADERS AMReX_BaseIFFAB.H	      AMReX_EBDataVarMacros.H		 AMReX_Factorial.H		      AMReX_IrregNode.H        AMReX_Notation.H		     AMReX_UnionIF.H             )
list ( APPEND ALLHEADERS AMReX_BaseIFFABI.H	      AMReX_EBDebugDump.H		 AMReX_FixedRefinement.H	      AMReX_KDStruct.H	       AMReX_PXStuff.H		     AMReX_VoFIterator.H         )
list ( APPEND ALLHEADERS AMReX_BaseIFFactory.H	      AMReX_EBDebugOut.H		 AMReX_FixedRefinementImplem.H	      AMReX_KDTree.H	       AMReX_PlaneIF.H		     AMReX_VolIndex.H            )
list ( APPEND ALLHEADERS AMReX_BaseIVFAB.H	      AMReX_EBFaceFAB.H			 AMReX_FlatPlateGeom.H		      AMReX_LSProblem.H        AMReX_PolyGeom.H		     AMReX_WrappedGShop.H        )
list ( APPEND ALLHEADERS AMReX_BaseIVFABI.H	      AMReX_EBFluxFAB.H			 AMReX_GenericArithmetic.H	      AMReX_LSProblemImplem.H  AMReX_PolynomialIF.H	     AMReX_ZCylinder.H           )
list ( APPEND ALLHEADERS AMReX_BaseIVFactory.H	      AMReX_EBFluxFactory.H		 AMReX_GenericArithmeticI.H	      AMReX_LSquares.H	       AMReX_RedistStencil.H                                     )
list ( APPEND ALLHEADERS AMReX_CH_EBIS_ORDER.H	      AMReX_EBGraph.H			 AMReX_GeomIntersectUtils.H	      AMReX_LatheIF.H	       AMReX_RefinementCriterion.H                               )
list ( APPEND ALLHEADERS AMReX_CellEdge.H	      AMReX_EBISBox.H			 AMReX_GeometryService.H	      AMReX_LoHiSide.H	       AMReX_STLAsciiReader.H                                    )
list ( APPEND ALLHEADERS AMReX_ComplementIF.H	      AMReX_EBISLayout.H		 AMReX_GeometryShop.H		      AMReX_MetaPrograms.H     AMReX_STLBox.H                                            )
list ( APPEND ALLHEADERS AMReX_ConstrainedLS.H	      AMReX_EBISLevel.H			 AMReX_GraphNode.H		      AMReX_MinimalCCCM.H      AMReX_STLExplorer.H                                       )    



list ( APPEND CXXSRC AMReX_AllRegularService.cpp     AMReX_EBFluxFactory.cpp		    AMReX_ExtrudeIF.cpp		  AMReX_IntersectionIF.cpp    AMReX_PXStuff.cpp		AMReX_SmoothAbsoluteValue.cpp)
list ( APPEND CXXSRC AMReX_BaseIFFAB.cpp	     AMReX_EBGraph.cpp			    AMReX_FabArrayIO.cpp	  AMReX_IrregFAB.cpp	      AMReX_PlaneIF.cpp		AMReX_SmoothIntersection.cpp )
list ( APPEND CXXSRC AMReX_CellEdge.cpp	             AMReX_EBISBox.cpp			    AMReX_FaceIndex.cpp		  AMReX_IrregFABFactory.cpp   AMReX_PolyGeom.cpp	AMReX_SmoothUnion.cpp                )
list ( APPEND CXXSRC AMReX_ComplementIF.cpp	     AMReX_EBISLayout.cpp		    AMReX_FaceIterator.cpp	  AMReX_IrregNode.cpp	      AMReX_PolynomialIF.cpp	AMReX_SphereIF.cpp           )
list ( APPEND CXXSRC AMReX_ConstrainedLS.cpp         AMReX_EBISLevel.cpp		    AMReX_FlatPlateGeom.cpp	  AMReX_KDTree.cpp	      AMReX_RedistStencil.cpp	AMReX_Stencils.cpp                   )
list ( APPEND CXXSRC AMReX_CutCellMoments.cpp        AMReX_EBIndexSpace.cpp		    AMReX_GeomIntersectUtils.cpp  AMReX_LSProblem.cpp	      AMReX_STLAsciiReader.cpp	AMReX_TransformIF.cpp        )
list ( APPEND CXXSRC AMReX_EBArith.cpp	             AMReX_EBLevelGrid.cpp		    AMReX_GeometryService.cpp	  AMReX_LSquares.cpp	      AMReX_STLBox.cpp		AMReX_UnionIF.cpp                    )
list ( APPEND CXXSRC AMReX_EBCellFAB.cpp	     AMReX_EBLevelRedist.cpp		    AMReX_GeometryShop.cpp	  AMReX_LatheIF.cpp	      AMReX_STLExplorer.cpp	AMReX_VoFIterator.cpp        )
list ( APPEND CXXSRC AMReX_EBData.cpp	             AMReX_EBLoHiCenter.cpp		    AMReX_GraphNode.cpp		  AMReX_LoHiSide.cpp	      AMReX_STLIF.cpp		AMReX_VolIndex.cpp                   )
list ( APPEND CXXSRC AMReX_EBDebugOut.cpp	     AMReX_EBNormalizeByVolumeFraction.cpp  AMReX_IFData.cpp		  AMReX_MinimalCCCM.cpp       AMReX_STLMesh.cpp		AMReX_WrappedGShop.cpp       )
list ( APPEND CXXSRC AMReX_EBFaceFAB.cpp	     AMReX_Ellipsoid.cpp		    AMReX_IFSlicer.cpp		  AMReX_Moments.cpp	      AMReX_STLUtil.cpp		AMReX_ZCylinder.cpp          )
list ( APPEND CXXSRC AMReX_EBFluxFAB.cpp	     AMReX_EllipsoidIF.cpp		    AMReX_IntVectSet.cpp	  AMReX_NormalDerivative.cpp                                    )    

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
