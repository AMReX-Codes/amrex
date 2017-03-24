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
      if (a_fixRegularNextToMultiValued)
      {
        fixRegularNextToMultiValued();
      }
    }
  }
  ///
  void
  EBISLevel::defineFromGeometryService(const GeometryService & a_geoserver)
  {
    
    m_grids.define(m_domain);
    m_grids.maxSize(maxgrid);
    DistributionMapping dm(m_grids);
    m_graph.define(m_grids, dm, 1, 1);
    m_data .define(m_grids, dm, 1, 1);
    for (MFIter mfi(m_graph), mfi.isValid(); ++mfi);
    {
      Box valid  = mfi.validbox();
      Box ghostRegion = valid;
      ghostRegion.grow(1);
      ghostRegion &= a_domain;

      EBGraph& ebgraph = a_graph[dit()];
      EBGraph& ebdata  = a_graph[dit()];
      GeometryService::InOut inout = a_geoserver.InsideOutside(ghostRegion, m_domain, m_origin, m_dx);

      EBGraph& ebgraph= m_graph[mfi];
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
        std::vector<IrregNode>&  nodes;

        a_geoserver.fillGraph(regIrregCovered, nodes, validRegion,
                              ghostRegion, m_domain,
                              m_origin, m_dx);

        ebgraph.buildGraph(regIrregCovered, nodes, ghostRegion, m_domain);
      }
      ebdata.define(ebgraph, nodes, validBox, ghostRegion);
    }
    ebdata.FillBoundary();
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
    coarsenVoFs(a_fineEBIS);

    //overallMemoryUsage();
    //create coarse faces from fine
    coarsenFaces(a_fineEBIS);
    //overallMemoryUsage();
    //fix the regular next to the multivalued cells
    //to be full irregular cells
    //  dumpDebug(string("EBIS::before FRNTM"));
    if (a_fixRegularNextToMultiValued)
    {
      fixRegularNextToMultiValued();
    }
    //  dumpDebug(string("EBIS::after FRNTM"));

    //overallMemoryUsage();
    // fix the fine->coarseVoF thing.
    fixFineToCoarse(a_fineEBIS);
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
  coarsenVoFs(EBISLevel& a_fineEBIS)
  {
    BL_PROFILE("EBISLevel::coarsenVoFs");

    BoxArray gridsReCo = m_grids;
    gridsReCo.refine(2);

    DistributionMapping dmfc(gridsReCo);
    m_graph.define(m_grids, dm, 1, 1);
    FabArray<EBGraph> ebgraphReCo(gridsReCo, dmfc, 1, 1);

    int srcGhost =0; int dstGhost = 1;
    ebgraphReCo.copy(a_fineEBIS.m_graph, 0, 0, 1, srcGhost, dstGhost);

    for (MFIter mfi(m_graph), mfi.isValid(); ++mfi);
    {
      const EBGraph& fineEBGraph = ebgraphReCo[mfi];
      const Box& coarRegion      = m_grids[mfi];
      EBGraph& coarEBGraph       = m_graph[mfi];
      coarEBGraph.coarsenVoFs(fineEBGraph, coarRegion);
    }

    std::shared_ptr<FabArray<EBGraph> > graphptrCoar(    m_graph, &null_deleter_fab_ebg);
    std::shared_ptr<FabArray<EBGraph> > graphptrReCo(ebgraphReCo, &null_deleter_fab_ebg);
    EBDataFactory ebdfCoar(graphptrCoar);
    EBDataFactory ebdfReCo(graphptrReCo);
    FabArray<EBData> ebdataReCo
    m_data    .define(a_grids  , dm, 1, m_nghost, MFInfo(), ebdfCoar);
    ebdataReCo.define(gridsReCo, dm, 1, m_nghost, MFInfo(), ebdfReCo);
    ebdataReCo.copy(a_fineEBIS.m_data, 0, 0, 1, srcGhost, dstGhost);    
    a_fineEBIS.m_data.copyTo(interv, ebdataReCo, interv);

    for (MFIter mfi(m_graph), mfi.isValid(); ++mfi);
    {
      const EBGraph& fineEBGraph = ebgraphReCo[mfi];
      const EBGraph& coarEBGraph =     m_graph[mfi];
      const EBData & fineEBData  =  ebdataReCo[mfi];
      EBData       & coarEBData  =      m_data[mfi];

      m_data[mfi].coarsenVoFs(fineEBData, fineEBGraph, coarEBGraph, m_grids[mfi]);
    }
  }

  //now fix the multivalued next to regular thing for the graph and the data
  //the oldgraph/newgraph thing is necessary because the graphs are
  //reference counted and they have to be kept consistent with the data
  void EBISLevel::fixRegularNextToMultiValued()
  {
    BL_PROFILE("EBISLevel::fixRegularNextToMultiValued");

    EBGraphFactory graphfact(m_domain);
    LayoutData<IntVectSet> vofsToChange(m_grids);
    LevelData<EBGraph> oldGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
    Interval interv(0,0);

    m_graph.copyTo(interv, oldGhostGraph, interv);

    for (MFIter mfi(m_graph), mfi.isValid(); ++mfi);
    {
      BL_PROFILE("EBISLevel::fixRegularNextToMultiValued_loop1");
      m_graph[mfi].getRegNextToMultiValued(vofsToChange[dit()],
                                             oldGhostGraph[dit()]);

      m_graph[dit()].addFullIrregularVoFs(vofsToChange[ dit()],
                                          oldGhostGraph[dit()]);
    }

    EBDataFactory datafact;
    LevelData<EBGraph> newGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
    LevelData<EBData>  newGhostData(m_grids, 1, IntVect::Unit, datafact);

    m_graph.copyTo(interv, newGhostGraph, interv);

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      BL_PROFILE("EBISLevel::fixRegularNextToMultiValued_loop2");

      Box localBox = m_grids.get(dit());
      localBox.grow(1);
      localBox &= m_domain;
      newGhostData[dit()].defineVoFData(oldGhostGraph[dit()],  localBox);
      newGhostData[dit()].defineFaceData(oldGhostGraph[dit()], localBox);
    }

    m_data.copyTo(interv,  newGhostData,  interv);

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      BL_PROFILE("EBISLevel::fixRegularNextToMultiValued_loop3");

      m_data[dit() ].addFullIrregularVoFs(vofsToChange[dit()],
                                          newGhostGraph[dit()],
                                          newGhostData[dit()].getVolData(),
                                          oldGhostGraph[dit()]);
    }
  }

  void EBISLevel::fixFineToCoarse(EBISLevel& a_fineEBIS)
  {
    BL_PROFILE("EBISLevel::fixFineToCoarse");
    // make a coarse layout from the fine layout so that we
    //can fix the fine->coarseVoF thing.
    DisjointBoxLayout coarFromFineDBL;

    coarsen(coarFromFineDBL, a_fineEBIS.m_graph.getBoxes(), 2);
    coarFromFineDBL.close();
    EBGraphFactory ebgraphfact(m_domain);
    LevelData<EBGraph> coarFromFineEBGraph(coarFromFineDBL,1, IntVect::Zero, ebgraphfact);
    Interval interv(0,0);
    m_graph.copyTo(interv, coarFromFineEBGraph, interv);

    for (DataIterator dit = a_fineEBIS.m_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBGraph& fineEBGraph       =  a_fineEBIS.m_graph[dit()];
      const EBGraph& coarEBGraph = coarFromFineEBGraph[dit()];

      coarEBGraph.fixFineToCoarse(fineEBGraph);
    }

  }

  void EBISLevel::coarsenFaces(EBISLevel& a_fineEBIS)
  {
    BL_PROFILE("EBISLevel::coarsenFaces");
    //now make a fine ebislayout with two ghost cell
    //on the same mapping as m_dbl
    //so that i can do the vofs and faces of this level.
    //this one will have the fine from coarse stuff fixed
    DisjointBoxLayout fineFromCoarDBL;
    refine(fineFromCoarDBL, m_grids, 2);
    fineFromCoarDBL.close();

    //no need for ghost cells here
    EBGraphFactory ebgraphfactfine(a_fineEBIS.m_domain);
    EBGraphFactory ebgraphfactcoar(m_domain);
    LevelData<EBGraph> fineEBGraphGhostLD(fineFromCoarDBL,1,3*IntVect::Unit, ebgraphfactfine);
    Interval interv(0,0);
    a_fineEBIS.m_graph.copyTo(interv, fineEBGraphGhostLD, interv);
    LevelData<EBGraph> coarEBGraphGhostLD(m_grids,        1,  IntVect::Unit, ebgraphfactcoar);
    m_graph.copyTo(           interv, coarEBGraphGhostLD, interv);

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& fineEBGraphGhost = fineEBGraphGhostLD[dit()];
      const EBGraph& coarEBGraphGhost = coarEBGraphGhostLD[dit()];
      EBGraph& coarEBGraph = m_graph[dit()];
      coarEBGraph.coarsenFaces(coarEBGraphGhost, fineEBGraphGhost);
    }
    //redefine coarebghostgraphld so i can use the faces for the ebdata
    coarEBGraphGhostLD.define(m_grids, 1,  IntVect::Unit, ebgraphfactcoar);
    m_graph.copyTo(interv, coarEBGraphGhostLD, interv);

    EBDataFactory ebdatafact;
    LevelData<EBData> fineEBDataGhostLD(fineFromCoarDBL,1, 2*IntVect::Unit, ebdatafact);
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box localBox = grow(fineFromCoarDBL.get(dit()), 2);
      localBox &= a_fineEBIS.m_domain;
      fineEBDataGhostLD[dit()].defineVoFData(fineEBGraphGhostLD[dit()], localBox);;
      fineEBDataGhostLD[dit()].defineFaceData(fineEBGraphGhostLD[dit()], localBox);
    }
    a_fineEBIS.m_data.copyTo(interv, fineEBDataGhostLD, interv);

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBData&   fineEBData      = fineEBDataGhostLD[dit()];
      const EBGraph& fineEBGraphGhost = fineEBGraphGhostLD[dit()];
      const EBGraph& coarEBGraphGhost = coarEBGraphGhostLD[dit()];

      EBData& coarEBData   = m_data[dit()];
      coarEBData.coarsenFaces(fineEBData,  fineEBGraphGhost, coarEBGraphGhost, m_grids.get(dit()));
    }

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

