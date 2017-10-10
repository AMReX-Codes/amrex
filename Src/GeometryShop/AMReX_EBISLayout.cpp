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
#include "AMReX_VoFIterator.H"
#include "AMReX_EBDataFactory.H"
#include "AMReX_parstream.H"

namespace amrex
{
  static const IntVect   ebisl_debiv(D_DECL(994,213,7));
  static const IntVect   ebisl_debivlo(D_DECL(29,12, 0));
  static const IntVect   ebisl_debivhi(D_DECL(29,13, 0));
  static const VolIndex  ebisl_debvoflo(ebisl_debivlo, 0);
  static const VolIndex  ebisl_debvofhi(ebisl_debivhi, 0);
  static const FaceIndex ebisl_debface( ebisl_debvoflo, ebisl_debvofhi);

  void EBISL_checkGraph(const BoxArray          & a_grids,
                        const DistributionMapping & a_dm,
                        const FabArray<EBGraph> & a_graph,
                        const string & a_identifier)
  {
    for(MFIter mfi(a_grids, a_dm); mfi.isValid(); ++mfi)
    {
      const EBGraph & graph = a_graph[mfi];
      Box region     = graph.getRegion();
      Box fullRegion = graph.getFullRegion();
      if(!fullRegion.contains(region))
      {
        amrex::AllPrint() << "ebislayout:" << a_identifier;
        amrex::AllPrint() << ", region = " << region << ", fullRegion = " << fullRegion  << endl;
      }
//      const Box& region = a_grids[mfi];
//      const EBGraph & graph = a_graph[mfi];
//      if(region.contains(ebisl_debiv))
//      {
//        amrex::AllPrint() << "ebislayout:" << a_identifier;
//        int ireg = 0; int icov = 0;
//        if(graph.isRegular(ebisl_debiv))
//        {
//          ireg = 1;
//        }
//        if(graph.isCovered(ebisl_debiv))
//        {
//          icov = 1;
//        }
//        amrex::AllPrint() << ", ireg = " << ireg << ", icov = " << icov << endl;
//      }
    }
  }
  void EBISL_checkData(const BoxArray          & a_grids,
                        const DistributionMapping & a_dm,
                        const FabArray<EBData> & a_data,
                        const string & a_identifier)
  {
    for(MFIter mfi(a_grids, a_dm); mfi.isValid(); ++mfi)
    {
      const Box& region   = a_grids[mfi];
      const EBData & data = a_data[mfi];
      const BaseIFFAB<Real>& fab = data.getFaceData(1);
      if(fab.hasFace(ebisl_debface))
      {
        Real value = fab(ebisl_debface, 0);
        pout() << "ebislayout:" << a_identifier;
        pout() << ", region = " << region << ", facedat(" << ebisl_debface << ", 0)=" << value << endl;
      }
    }
  }
  /****************/
  void
  EBISLayoutImplem::define(const Box               & a_domain,
                           const BoxArray          & a_grids,
                           const DistributionMapping & a_dm,
                           const int               & a_nghost,
                           const FabArray<EBGraph> & a_graph,
                           const FabArray<EBData>  & a_data)
  {
    BL_PROFILE("EBISLayoutImplem::define");
    //pout() << "in ebislayoutimplem::define with nghost = " << a_nghost << endl;
    m_dm = a_dm;
    m_domain = a_domain;
    m_nghost = a_nghost;
    m_dblInputDom = a_grids;
    m_fineLevels.resize(0);
    m_coarLevels.resize(0);
    m_maxCoarseningRatio = 2;
    m_maxRefinementRatio = 1;//ug--face refinement means you have to have to do this once.
    int dstGhostData = a_nghost;
    int dstGhostGraph = a_nghost+2; //because of irregular faces at box boundaries
    int srcGhost = 0;

//begin debug
//    BoxArray inpgrids = a_graph.boxArray();
//    DistributionMapping inpdm = a_graph.DistributionMap();
//    EBISL_checkGraph(inpgrids, inpdm, a_graph, string(" input graph"));
//    BoxArray inpgrids = a_data.boxArray();
//    DistributionMapping inpdm = a_data.DistributionMap();
//    EBISL_checkData(inpgrids, inpdm, a_data, string(" input to ebisl::define"));
// end debug

    m_ebGraph = shared_ptr<FabArray<EBGraph> >(new FabArray<EBGraph>(a_grids, a_dm, 1, dstGhostGraph,
                                                                     MFInfo(),DefaultFabFactory<EBGraph>()));
    
    m_ebGraph->copy(a_graph, 0, 0, 1, srcGhost, dstGhostGraph);

    EBDataFactory ebdatafact(m_ebGraph);
    m_ebData  = shared_ptr<FabArray<EBData > >(new FabArray<EBData>(a_grids, a_dm, 1, m_nghost, MFInfo(), ebdatafact));
      
    BL_PROFILE_VAR("EBISLayout_copy_ebdata",copy_data);
    m_ebData ->copy(a_data , 0, 0, 1, srcGhost, dstGhostData);
    BL_PROFILE_VAR_STOP(copy_data);

//begin debug
//    EBISL_checkGraph(a_grids, a_dm, *m_ebGraph, string(" my graph after copy"));
//    EBISL_checkData(a_grids, a_dm, *m_ebData, string(" output from ebisl::define"));
// end debug
      
    m_defined = true;
  }
      
  /****************/
  EBISBox
  EBISLayoutImplem::operator[](const MFIter & a_dit) const
  {
    EBISBox retval((*m_ebGraph)[a_dit], (*m_ebData)[a_dit]);
    return retval;
  }
  /****************/
  EBISBox
  EBISLayoutImplem::operator[](const int& a_boxindex) const
  {
    EBISBox retval((*m_ebGraph)[a_boxindex], (*m_ebData)[a_boxindex]);
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
    EBISBox ebisBoxFine((*m_ebGraph)[a_mfi], (*m_ebData)[a_mfi]);
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
  Vector<VolIndex>
  EBISLayoutImplem::refine(const VolIndex  & a_vof,
                           const int       & a_ratio,
                           const MFIter    & a_mfi) const
  {
    BL_ASSERT(a_ratio > 0);
    BL_ASSERT(a_ratio <= m_maxRefinementRatio);
    BL_ASSERT(a_ratio%2 == 0);
      
    //for ratio of 2, just use ebisbox
    EBISBox ebisBoxCoar = (*this)[a_mfi];
    Vector<VolIndex> fineVoFs = ebisBoxCoar.refine(a_vof);
    //for ratio > 2, chase its tail
    int ifinelev = 0;
    for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      Vector<VolIndex> coarVoFs = fineVoFs;
      
      const EBISLayout& ebisl = m_fineLevels[ifinelev];
      EBISBox ebisBox = ebisl[a_mfi];
      fineVoFs.resize(0);
      for (int ivof = 0; ivof < coarVoFs.size(); ivof++)
      {
        Vector<VolIndex> newvofs = ebisBox.refine(coarVoFs[ivof]);
        fineVoFs.insert(fineVoFs.end(), newvofs.begin(), newvofs.end());
      }
      ifinelev++;
    }
    return fineVoFs;
      
  }
  /****************/
  Vector<FaceIndex>
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
    Vector<FaceIndex> fineFaces = ebisBoxCoar2.refine(a_face,ebisBoxFine2);
    //for ratio > 2, chase its tail
    int ifinelev = 0;
    for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      Vector<FaceIndex> coarFaces = fineFaces;
      
      EBISBox ebisBoxCoar = m_fineLevels[ifinelev    ][a_mfi];
      EBISBox ebisBoxFine = m_fineLevels[ifinelev + 1][a_mfi];
      fineFaces.resize(0);
      for (int iface = 0; iface < coarFaces.size(); iface++)
      {
        Vector<FaceIndex> newfaces = ebisBoxCoar.refine(coarFaces[iface],ebisBoxFine);
        fineFaces.insert(fineFaces.end(), newfaces.begin(), newfaces.end());
      }
      ifinelev++;
    }
    return fineFaces;
      
  }
  /****************/
  void
  EBISLayoutImplem::setMaxRefinementRatio(const int& a_maxRefine)
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
    const EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int irat = 2;
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
      Box fineDomain = m_domain;
      BoxArray fineBLDomain = m_dblInputDom;
      fineDomain.  refine(irat);
      fineBLDomain.refine(irat);
      ebisPtr->fillEBISLayout(m_fineLevels[ilev],  fineBLDomain, m_dm, fineDomain, m_nghost);
      irat *= 2;
    }
  }
  /****************/
  void
  EBISLayoutImplem::setMaxCoarseningRatio(const int&                a_maxCoarsen)
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
    const EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
      Box  coarDomain  = m_domain;
      BoxArray  coarBLDomain = m_dblInputDom;
      coarDomain.  coarsen(irat);
      coarBLDomain.coarsen(irat);
      
      ebisPtr->fillEBISLayout(m_coarLevels[ilev], coarBLDomain, m_dm, coarDomain, m_nghost);
      irat *= 2;
    }
  }
}
