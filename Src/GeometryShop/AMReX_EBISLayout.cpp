/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "AMReX_EBISLayout.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBArith.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBDataFactory.H"

namespace amrex
{
  //so i can pass the ebgraph to ebdatafactory
  void null_deleter_fab_ebg(FabArray<EBGraph>* a_ptr)
  {
  }

/****************/
  void
  EBISLayoutImplem::define(const Box               & a_domain,
                           const BoxArray          & a_grids,
                           const int               & a_nghost,
                           const FabArray<EBGraph> & a_graph,
                           const FabArray<EBData>  & a_data)
  {
    BL_PROFILE("EBISLayoutImplem::define");

    m_domain = a_domain;
    m_nghost = a_nghost;
    m_dblInputDom = a_grids;
    m_fineLevels.resize(0);
    m_coarLevels.resize(0);
    m_maxCoarseningRatio = 2;
    m_maxRefinementRatio = 1;//ug--face refinement means you have to have to do this once.


    DistributionMapping dm(a_grids);
    m_ebGraph.define(a_grids, dm, 1, m_nghost);

    std::shared_ptr<FabArray<EBGraph> > graphptr(m_ebGraph, &null_deleter_fab_ebg);
    EBDataFactory ebdatafact(graphptr);
    m_ebData .define(a_grids, dm, 1, m_nghost, MFInfo(), ebdatafact);

    m_ebGraph.copy(a_graph, 0, 0, 1);
    m_ebData .copy(a_data , 0, 0, 1);
    m_ebGraph.FillBoundary();
    m_ebData .FillBoundary();

    m_defined = true;
  }

/****************/
  const EBISBox&
  EBISLayoutImplem::operator[](const MFIter & a_dit) const
  {
    EBISBox retval(m_ebGraph[a_dit], m_ebData[a_dit]);
    return retval;
  }
/****************/
  EBISLayoutImplem::EBISLayoutImplem()
  {
    m_maxCoarseningRatio = 1;
    m_maxRefinementRatio = 1;
    m_defined = false;
  }
/****************/
  EBISLayoutImplem::~EBISLayoutImplem()
  {

  }
/****************/
  VolIndex
  EBISLayoutImplem::coarsen(const VolIndex & a_vof,
                            const int      & a_ratio,
                            const MFIter   & a_mfi) const
  {
    BL_ASSERT(a_ratio > 0);
    BL_ASSERT(a_ratio <= m_maxCoarseningRatio);
    BL_ASSERT(a_ratio%2 == 0);

    //for ratio of 2, just use ebisbox
    EBISBox ebisBoxFine(m_ebGraph[a_mfi], m_ebData[a_mfi]);
    VolIndex coarVoF = ebisBoxFine.coarsen(a_vof);
    //for ratio > 2, chase its tail
    int icoarlev = 0;
    for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      VolIndex fineVoF = coarVoF;
      const EBISLayout& ebisl = m_coarLevels[icoarlev];
      EBISBox  ebisBox = ebisl[a_mfi];
      coarVoF = ebisBox.coarsen(fineVoF);
      icoarlev++;
    }
    return coarVoF;
  }
/****************/
  FaceIndex
  EBISLayoutImplem::coarsen(const FaceIndex & a_face,
                            const int       & a_ratio,
                            const MFIter    & a_mfi) const
  {
    BL_ASSERT(a_ratio > 0);
    BL_ASSERT(a_ratio <= m_maxCoarseningRatio);
    BL_ASSERT(a_ratio%2 == 0);

    //for ratio of 2, just use ebisbox
    EBISBox ebisBoxFine = (*this)[a_mfi];
    FaceIndex coarFace = ebisBoxFine.coarsen(a_face);
    //for ratio > 2, chase its tail
    int icoarlev = 0;
    for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      FaceIndex fineFace = coarFace;
      const EBISLayout& ebisl = m_coarLevels[icoarlev];
      EBISBox ebisBox = ebisl[a_mfi];
      coarFace = ebisBox.coarsen(fineFace);
      icoarlev++;
    }
    return coarFace;
  }
/****************/
  st::vector<VolIndex>
  EBISLayoutImplem::refine(const VolIndex  & a_vof,
                           const int       & a_ratio,
                           const MFIter    & a_mfi) const
  {
    BL_ASSERT(a_ratio > 0);
    BL_ASSERT(a_ratio <= m_maxRefinementRatio);
    BL_ASSERT(a_ratio%2 == 0);

    //for ratio of 2, just use ebisbox
    EBISBox ebisBoxCoar = (*this)[a_mfi];
    std::vector<VolIndex> fineVoFs = ebisBoxCoar.refine(a_vof);
    //for ratio > 2, chase its tail
    int ifinelev = 0;
    for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      std::vector<VolIndex> coarVoFs = fineVoFs;

      const EBISLayout& ebisl = m_fineLevels[ifinelev];
      EBISBox ebisBox = ebisl[a_mfi];
      fineVoFs.resize(0);
      for (int ivof = 0; ivof < coarVoFs.size(); ivof++)
      {
        fineVoFs.append(ebisBox.refine(coarVoFs[ivof]));
      }
      ifinelev++;
    }
    return fineVoFs;

  }
/****************/
  std::vector<FaceIndex>
  EBISLayoutImplem::refine(const FaceIndex & a_face,
                           const int       & a_ratio,
                           const MFIter    & a_mfi) const
  {
    BL_ASSERT(a_ratio > 0);
    BL_ASSERT(a_ratio <= m_maxRefinementRatio);
    BL_ASSERT(a_ratio%2 == 0);

    //for ratio of 2, just use ebisbox
    EBISBox ebisBoxCoar2 =         (*this)[a_mfi];
    EBISBox ebisBoxFine2 = m_fineLevels[0][a_mfi];
    std::vector<FaceIndex> fineFaces = ebisBoxCoar2.refine(a_face,ebisBoxFine2);
    //for ratio > 2, chase its tail
    int ifinelev = 0;
    for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      std::vector<FaceIndex> coarFaces = fineFaces;

      EBISBox ebisBoxCoar = m_fineLevels[ifinelev    ][a_mfi];
      EBISBox ebisBoxFine = m_fineLevels[ifinelev + 1][a_mfi];
      fineFaces.resize(0);
      for (int iface = 0; iface < coarFaces.size(); iface++)
      {
        fineFaces.append(ebisBoxCoar.refine(coarFaces[iface],ebisBoxFine));
      }
      ifinelev++;
    }
    return fineFaces;

  }
/****************/
  void
  EBISLayoutImplem::setMaxRefinementRatio(const int& a_maxRefine, const EBIndexSpace* ebisPtr)
  {
    BL_ASSERT(a_maxRefine % 2 == 0);
    BL_ASSERT(a_maxRefine > 0);
    if (a_maxRefine <= m_maxRefinementRatio)
    {
      return;
    }

    m_maxRefinementRatio = a_maxRefine;
    //figure out how many levels i will need
    int nlevels = 0;
    for (int irat = 2; irat  <= a_maxRefine; irat *= 2)
    {
      nlevels++;
    }

    m_fineLevels.resize(nlevels);
    int irat = 2;
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
      Box fineDomain = m_domain;
      BoxArray fineBLDomain = m_dblInputDom;
      fineDomain.  refine(irat);
      fineBLDomain.refine(irat);
      ebisPtr->fillEBISLayout(m_fineLevels[ilev], fineBLDomain, fineDomain, m_nghost);
      irat *= 2;
    }
  }
/****************/
  void
  EBISLayoutImplem::setMaxCoarseningRatio(const int&                a_maxCoarsen,
                                          const EBIndexSpace* const a_ebisPtr)
  {
    BL_ASSERT(a_maxCoarsen % 2 == 0);
    BL_ASSERT(a_maxCoarsen > 0);
    if (a_maxCoarsen <= m_maxCoarseningRatio)
    {
      return;
    }

    m_maxCoarseningRatio = a_maxCoarsen;
    //figure out how many levels i will need
    int nlevels = 0;
    for (int irat = 4; irat  <= a_maxCoarsen; irat *= 2)
    {
      nlevels++;
    }

    m_coarLevels.resize(nlevels);
    int irat = 2;
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
      Box  coarDomain  = m_domain;
      BoxArray  coarBLDomain = m_dblInputDom;
      coarDomain.  coarsen(irat);
      coarBLDomain.coarsen(irat);

      a_ebisPtr->fillEBISLayout(m_coarLevels[ilev], coarBLDomain, coarDomain, m_nghost);
      irat *= 2;
    }
  }
/****************/
  int
  EBISLayoutImplem::getMaxCoarseningRatio() const
  {
    return m_maxCoarseningRatio;
  }
/****************/
  int
  EBISLayoutImplem::getMaxRefinementRatio() const
  {
    return m_maxRefinementRatio;
  }
/****************/
}
