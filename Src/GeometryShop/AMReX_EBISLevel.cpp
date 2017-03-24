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



#include "AMReX_Print.H"
#include "AMReX_EBISLevel.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_IrregNode.H"
#include "AMReX_PolyGeom.H"
#include "AMReX_FaceIterator.H"


namespace amrex
{
  EBIndexSpace* AMReX_EBIS::s_instance = NULL;
  ///
  EBIndexSpace* 
  AMReX_EBIS::
  instance()
  {
    if (s_instance == NULL)
    {
      s_instance = new EBIndexSpace();
    }

    return  s_instance;
  }



  EBISLevel::EBISLevel(const Box             & a_domain,
                       const RealVect        & a_origin,
                       const Real            & a_dx,
                       const int             & a_nCellMax,
                       const GeometryService & a_geoserver)
  {
    // this is the method called by EBIndexSpace::buildFirstLevel
    BL_PROFILE("EBISLevel::EBISLevel_geoserver_domain");
    m_domain = a_domain;
    m_dx = a_dx;
    m_origin = a_origin;
    m_nCellMax = a_nCellMax;

    defineFromGeometryService(a_geoserver);


    if(a_geoserver.canGenerateMultiCells())
    {
      fixRegularNextToMultiValued();
    }
  }
  ///
  void
  EBISLevel::defineFromGeometryService(const GeometryService & a_geoserver)
  {
    
    m_grids.define(m_domain);
    m_grids.maxSize(m_nCellMax);
    DistributionMapping dm(m_grids);
    m_graph.define(m_grids, dm, 1, 1);
    m_data .define(m_grids, dm, 1, 1);
    for (MFIter mfi(m_graph); mfi.isValid(); ++mfi)
    {
      Box valid  = mfi.validbox();
      Box ghostRegion = valid;
      ghostRegion.grow(1);
      ghostRegion &= m_domain;

      EBGraph& ebgraph = m_graph[mfi];
      EBData& ebdata   = m_data [mfi];
      GeometryService::InOut inout = a_geoserver.InsideOutside(ghostRegion, m_domain, m_origin, m_dx);

      if (inout == GeometryService::Regular)
      {
        ebgraph.setToAllRegular();
      }
      else if (inout == GeometryService::Covered)
      {
        ebgraph.setToAllCovered();
      }
      else
      {
        BaseFab<int>             regIrregCovered;
        std::vector<IrregNode>   nodes;

        a_geoserver.fillGraph(regIrregCovered, nodes, valid,
                              ghostRegion, m_domain,
                              m_origin, m_dx);

        ebgraph.buildGraph(regIrregCovered, nodes, ghostRegion, m_domain);
        ebdata.define(ebgraph, nodes, valid, ghostRegion);
      }
    }
    m_graph.FillBoundary();
    m_data. FillBoundary();
  }
  ///
  EBISLevel::
  EBISLevel(EBISLevel             & a_fineEBIS,
            const GeometryService & a_geoserver)
  { // method used by EBIndexSpace::buildNextLevel
    BL_PROFILE("EBISLevel::EBISLevel_fineEBIS");

    m_domain = a_fineEBIS.m_domain;
    m_domain.coarsen(2);
    m_dx = 2.*a_fineEBIS.m_dx;
    m_origin = a_fineEBIS.m_origin;

    //create coarsened vofs from fine.
    //create coarse faces from fine
    coarsenVoFsAndFaces(a_fineEBIS);

    //fix the regular next to the multivalued cells
    //to be full irregular cells
    fixRegularNextToMultiValued();
  }

  //steps to coarsen an ebislevel:
  //1. coarsen vofs
  //1a. make a layout over refine(mydbl,2) and copy
  //    fine layout into it
  //1b.do connectivity bizbaz to make my vofs, volfrac, vof->fineVofs
  //2. make faces doing connectivity jive
  //2.23 make coarse geometric stuff (centroids and all that) from fine
  //3. make coarse layout from coarsen(finelayout). and copy my data into it
  //   to make fine->coarserVoF
  // (so finer ebislevel does change in this function)
  void 
  EBISLevel::
  coarsenVoFsAndFaces(EBISLevel& a_fineEBIS)
  {
    BL_PROFILE("EBISLevel::coarsenVoFsAndFaces");

    BoxArray gridsReCo = m_grids;
    gridsReCo.refine(2);

    DistributionMapping dmco(m_grids);
    DistributionMapping dmfc(gridsReCo);
    //need two because of coarsen faces
    m_graph.define(m_grids, dmco, 1, 2);
    FabArray<EBGraph> ebgraphReCo(gridsReCo, dmfc, 1, 3);

    int srcGhost =0; int dstGhost = 3;
    ebgraphReCo.copy(a_fineEBIS.m_graph, 0, 0, 1, srcGhost, dstGhost);
    ///first deal with the graph
    for (MFIter mfi(m_graph); mfi.isValid(); ++mfi)
    {
      const EBGraph& fineEBGraph = ebgraphReCo[mfi];
      EBGraph      & coarEBGraph = m_graph[mfi];
      const Box    & coarRegion  = mfi.validbox();
      coarEBGraph.coarsenVoFs(fineEBGraph, coarRegion);
    }
    m_graph.FillBoundary();
    for (MFIter mfi(m_graph); mfi.isValid(); ++mfi)
    {
      const EBGraph& fineEBGraph = ebgraphReCo[mfi];
      EBGraph      & coarEBGraph = m_graph[mfi];
      const Box    & coarRegion  = mfi.validbox();
      coarEBGraph.coarsenFaces(coarRegion, coarRegion);
      coarEBGraph.fixFineToCoarse(fineEBGraph);
    }
    m_graph.FillBoundary();

    //now deal with the data
    std::shared_ptr<FabArray<EBGraph> > graphptrCoar(    m_graph, &null_deleter_fab_ebg);
    std::shared_ptr<FabArray<EBGraph> > graphptrReCo(ebgraphReCo, &null_deleter_fab_ebg);
    EBDataFactory ebdfCoar(graphptrCoar);
    EBDataFactory ebdfReCo(graphptrReCo);
    FabArray<EBData> ebdataReCo;
    int nghostData = 1;

    m_data    .define(a_grids  , dm, 1, 0         , MFInfo(), ebdfCoar);
    ebdataReCo.define(gridsReCo, dm, 1, nghostData, MFInfo(), ebdfReCo);
    dstGhost = 1;
    ebdataReCo.copy(a_fineEBIS.m_data, 0, 0, 1, srcGhost, dstGhost);    
    for (MFIter mfi(m_graph); mfi.isValid(); ++mfi)
    {
      const EBGraph& fineEBGraph = ebgraphReCo[mfi];
      const EBGraph& coarEBGraph =     m_graph[mfi];
      const EBData & fineEBData  =  ebdataReCo[mfi];
      EBData       & coarEBData  =      m_data[mfi];

      m_data[mfi].coarsenVoFs (fineEBData, fineEBGraph, coarEBGraph, m_grids[mfi]);
      m_data[mfi].coarsenFaces(fineEBData, fineEBGraph, coarEBGraph, m_grids[mfi]);
    }
  }

  //now fix the multivalued next to regular thing for the graph and the data
  //the oldgraph/newgraph thing is necessary because the graphs are
  //reference counted and they have to be kept consistent with the data
  void 
  EBISLevel::
  fixRegularNextToMultiValued()
  {
    BL_PROFILE("EBISLevel::fixRegularNextToMultiValued");
    bool vofsAdded = false;
    for (MFIter mfi(m_graph); mfi.isValid(); ++mfi)
    {
      IntVectSet vofsToChange;
      BL_PROFILE("EBISLevel::fixRegularNextToMultiValued_loop1");
      m_graph[mfi].getRegNextToMultiValued(vofsToChange,
                                           m_grids[mfi]);

      m_graph[dit()].addFullIrregularVoFs(vofsToChange);


      m_data[dit() ].addFullIrregularVoFs(vofsToChange,
                                          m_graph[mfi]);

    }
    m_data .FillBoundary();
    m_graph.FillBoundary();
  }

  ///
  void 
  EBISLevel::
  fillEBISLayout(EBISLayout&              a_ebisLayout,
                 const DisjointBoxLayout& a_grids,
                 const int&               a_nghost) const
  {
    CH_assert(a_nghost >= 0);
  
    //a_ebisLayout.define(m_domain, a_grids, a_nghost, m_graph, m_data);
    //return; // caching disabled for now.... ugh.  bvs
    a_ebisLayout.define(m_domain, a_grids, a_nghost, m_graph, m_data);
  }


}

