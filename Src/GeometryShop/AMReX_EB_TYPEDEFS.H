#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

///meant to be overridden at compile time if necessary
#ifndef _EB_TYPEDEFS_H_
#define _EB_TYPEDEFS_H_


#include "AMReX_CH_EBIS_ORDER.H"
#include "AMReX_IndexedMoments.H"
#include "AMReX_MomentIterator.H"


///
typedef IndexedMoments<SpaceDim  , CH_EBIS_ORDER> IndMomSpaceDim;

///
typedef IndexedMoments<SpaceDim-1, CH_EBIS_ORDER> IndMomSDMinOne;

///
typedef IndexTM<int, SpaceDim>    IvSpaceDim;

///
typedef IndexTM<int, SpaceDim-1>  IvSDMinOne;

///
typedef IndexTM<Real,SpaceDim>    RvSpaceDim;

///
typedef IndexTM<Real,SpaceDim-1>  RvSDMinOne;

///
typedef MomentIterator<SpaceDim-1, CH_EBIS_ORDER> MomItSDMinOne;

///
typedef MomentIterator<SpaceDim  , CH_EBIS_ORDER> MomItSpaceDim;




#endif
