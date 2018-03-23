
#ifndef _EB_TYPEDEFS_H_
#define _EB_TYPEDEFS_H_


#include "AMReX_CH_EBIS_ORDER.H"
#include "AMReX_IndexedMoments.H"
#include "AMReX_MomentIterator.H"


///
typedef IndexedMoments<BL_SPACEDIM  , CH_EBIS_ORDER> IndMomSpaceDim;

///
typedef IndexedMoments<BL_SPACEDIM-1, CH_EBIS_ORDER> IndMomSDMinOne;

///
typedef IndexTM<int, BL_SPACEDIM>    IvSpaceDim;

///
typedef IndexTM<int, BL_SPACEDIM-1>  IvSDMinOne;

///
typedef IndexTM<Real,BL_SPACEDIM>    RvSpaceDim;

///
typedef IndexTM<Real,BL_SPACEDIM-1>  RvSDMinOne;

///
typedef MomentIterator<BL_SPACEDIM-1, CH_EBIS_ORDER> MomItSDMinOne;

///
typedef MomentIterator<BL_SPACEDIM  , CH_EBIS_ORDER> MomItSpaceDim;




#endif
