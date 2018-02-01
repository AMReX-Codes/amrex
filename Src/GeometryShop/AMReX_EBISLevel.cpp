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
#include "AMReX_EBDataFactory.H"
#include "AMReX_FaceIterator.H"
#include "AMReX_FabArrayIO.H"
#include "AMReX_Utility.H"
#include "AMReX_ParmParse.H"
#include "AMReX_VectorIO.H"


namespace amrex
{
  static std::string eb_surface_filename;

  static const IntVect   ebl_debiv(D_DECL(994,213,7));
  static const IntVect   ebl_debivlo(D_DECL(190,15,0));
  static const IntVect   ebl_debivhi(D_DECL(191,15,0));
  static const VolIndex  ebl_debvoflo(ebl_debivlo, 0);
  static const VolIndex  ebl_debvofhi(ebl_debivhi, 0);
  static const FaceIndex ebl_debface(ebl_debvoflo, ebl_debvofhi);

  void EBISLevel_checkGraph(const BoxArray          & a_grids,
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
        amrex::AllPrint() << "ebislevel:" << a_identifier;
        amrex::AllPrint() << ", region = " << region << ", fullRegion = " << fullRegion  << endl;
      }
//      const Box& region = a_grids[mfi];
//      if(region.contains(ebl_debiv))
//      {
//        amrex::AllPrint() << "ebislevel:" << a_identifier;
//        int ireg = 0; int icov = 0;
//        if(graph.isRegular(ebl_debiv))
//        {
//          ireg = 1;
//        }
//        if(graph.isCovered(ebl_debiv))
//        {
//          icov = 1;
//        }
//        amrex::AllPrint() << ", ireg = " << ireg << ", icov = " << icov << endl;
//      }
    }
  }
  void EBISLevel_checkData(const BoxArray            & a_grids,
                           const DistributionMapping & a_dm,
                           const FabArray<EBData>    & a_data,
                           const string              & a_identifier)
  {
    for(MFIter mfi(a_grids, a_dm); mfi.isValid(); ++mfi)
    {
      const Box& region = a_grids[mfi];
      const EBData & data = a_data[mfi];
      int ibox = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        const BaseIFFAB<Real>& facedat = data.getFaceData(idir);
        if(facedat.hasFace(ebl_debface))
        {
          Real areafrac = facedat(ebl_debface, 0);
          amrex::AllPrint() << "ebislevel:" << a_identifier;
          amrex::AllPrint() << ", ibox = "<< ibox << ", valid = " << region << ", areaFrac( " << ebl_debface << ") = " << areafrac << endl;
        }
      }
      ibox++;
    }
  }
  void 
  EBISLevel::
  write(const string& a_dirname) const
  {
    //this creates the directory of all the stuff
    UtilCreateCleanDirectory(a_dirname, true);
    writeHeader(a_dirname);
    string graphdirname = a_dirname + "/_graph";
    string  datadirname = a_dirname + "/_data";
//    UtilCreateDirectoryDestructive(graphdirname, true); done in functions below
//    UtilCreateDirectoryDestructive( datadirname, true);done in functions below
    FabArrayIO<EBGraph>::write(m_graph, graphdirname);
    FabArrayIO<EBData >::write(m_data ,  datadirname);
  }


  void 
  EBISLevel::
  read(const string& a_dirname)
  {
    readHeader(a_dirname);
    string graphdirname = a_dirname + "/_graph";
    string  datadirname = a_dirname + "/_data";

    FabArrayIO<EBGraph>::read(m_graph, graphdirname);
    FabArrayIO<EBData >::read(m_data ,  datadirname);
  }
  void 
  EBISLevel::
  writeHeader(const string& a_dirname) const
  {
    std::ofstream headerfile;
    string filename = a_dirname + string("/headerfile");
    headerfile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    headerfile << m_nCellMax << endl;
    headerfile << m_domain  << endl;
    headerfile << m_origin  << endl;
    headerfile << m_dx  << endl;

    //why box array has to be weird, who knows?
    //putting this last because it probably leaves the is in the wrong
    //place
    m_grids.writeOn(headerfile);

    headerfile.flush();
    headerfile.close();
  }


  void 
  EBISLevel::
  readHeader(const string& a_dirname)
  {
    std::ifstream headerfile;
    string filename = a_dirname + string("/headerfile");
    headerfile.open(filename.c_str(), std::ios::in);
    headerfile >> m_nCellMax;
    headerfile >> m_domain;
    headerfile >> m_origin;
    headerfile >> m_dx;

    //why box array has to be weird, who knows?
    //putting this last because it probably leaves the is in the wrong
    //place
    readBoxArray(m_grids, headerfile, false);
 
    headerfile.close();
  }
    
  void null_deleter_fab_ebg(FabArray<EBGraph>* a_ptr)
  {
  }

  EBIndexSpace* AMReX_EBIS::s_instance = nullptr;
  ///
  EBIndexSpace* 
  AMReX_EBIS::
  instance()
  {
    if (s_instance == nullptr)
    {
      s_instance = new EBIndexSpace();
    }

    return  s_instance;
  }

  void
  AMReX_EBIS::
  reset()
  {
    delete s_instance;
    s_instance = nullptr;
  }

  EBISLevel::EBISLevel(const Box             & a_domain,
                       const RealVect        & a_origin,
                       const Real            & a_dx,
                       const int             & a_nCellMax,
                       const GeometryService & a_geoserver)
      : m_nCellMax (a_nCellMax),
        m_domain   (a_domain),
        m_origin   (a_origin),
        m_dx       (a_dx)
  {
    // this is the method called by EBIndexSpace::buildFirstLevel
    BL_PROFILE("EBISLevel::EBISLevel()");

    ParmParse pp("ebis");
    m_build_eb_surface = false;

    int n_name = pp.countval("eb_surface_filename");
    if (n_name > 0) {
      m_build_eb_surface = true;
      pp.get("eb_surface_filename",eb_surface_filename);
    }

    defineFromGeometryService(a_geoserver);

    if(a_geoserver.canGenerateMultiCells())
    {
      fixRegularNextToMultiValued();
    }
  }

  void
  EBISLevel::buildEBSurface(const GeometryService & a_geoserver)
  {

#if AMREX_SPACEDIM==2
    Vector<Vector<Segment>> surfaceFragmentsG;
#elif AMREX_SPACEDIM==3
    Vector<Vector<Triangle>> surfaceFragmentsG;
#else
    amrex::Abort("buildEBSurface not applicable to 1D");
#endif

    for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
    {
      auto& eb_graph = m_graph[mfi];
      auto ivsirreg = eb_graph.getIrregCells(eb_graph.getRegion());
      auto& intersections = m_intersections[mfi];

#if AMREX_SPACEDIM==2
      list<Segment> segments;

      std::vector<IntVect> pt(4);
      std::vector<RealVect> x(4);
      std::vector<bool> v(4);

      surfaceFragmentsG.push_back(Vector<Segment>());
      auto& segVec = surfaceFragmentsG.back();

      for(IVSIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
      {
        pt[0] = ivsit();
        pt[1] = pt[0] + BASISV(0);
        pt[2] = pt[1] + BASISV(1);
        pt[3] = pt[0] + BASISV(1);

        for (int i=0; i<4; ++i)
        {
          for (int idir=0; idir<SpaceDim; ++idir)
          {
            x[i][idir] = m_origin[idir] + pt[i][idir]*m_dx;
          }
          v[i] = a_geoserver.pointOutside(x[i]);
        }

        int segCase = 0;
        if (v[0]) segCase |= 1;
        if (v[1]) segCase |= 2;
        if (v[2]) segCase |= 4;
        if (v[3]) segCase |= 8;

        switch (segCase)
        {
        case 1:
        case 14:
          segVec.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                   intersections.find(Edge(pt[0],pt[3]))));
          break;
        case 2:
        case 13:
          segVec.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                   intersections.find(Edge(pt[1],pt[2]))));
          break;
        case 3:
        case 12:
          segVec.push_back(Segment(intersections.find(Edge(pt[1],pt[2])),
                                   intersections.find(Edge(pt[0],pt[3]))));
          break;
        case 4:
        case 11:
          segVec.push_back(Segment(intersections.find(Edge(pt[1],pt[2])),
                                   intersections.find(Edge(pt[2],pt[3]))));
          break;
        case 6:
        case 9:
          segVec.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                   intersections.find(Edge(pt[2],pt[3]))));
          break;
        case 7:
        case 8:
          segVec.push_back(Segment(intersections.find(Edge(pt[2],pt[3])),
                                   intersections.find(Edge(pt[0],pt[3]))));
          break;
        case 5:
        case 10:
          segVec.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                   intersections.find(Edge(pt[1],pt[2]))));
          segVec.push_back(Segment(intersections.find(Edge(pt[2],pt[3])),
                                   intersections.find(Edge(pt[0],pt[3]))));
          break;
        }

      }
#else
      std::vector<IntVect> pt(8);
      std::vector<RealVect> x(8);
      std::vector<bool> v(8);

      surfaceFragmentsG.push_back(Vector<Triangle>());

      for(IVSIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
      {
        pt[0] = ivsit();
        pt[1] = pt[0] + BASISV(0);
        pt[2] = pt[1] + BASISV(1);
        pt[3] = pt[0] + BASISV(1);
        pt[4] = pt[0] + BASISV(2);
        pt[5] = pt[4] + BASISV(0);
        pt[6] = pt[5] + BASISV(1);
        pt[7] = pt[4] + BASISV(1);

        for (int i=0; i<8; ++i)
        {
          for (int idir=0; idir<SpaceDim; ++idir)
          {
            x[i][idir] = m_origin[idir] + pt[i][idir]*m_dx;
          }
          v[i] = a_geoserver.pointOutside(x[i]);
        }

        int cubeindex = 0;
        if (v[0]) cubeindex |= 1;
        if (v[1]) cubeindex |= 2;
        if (v[2]) cubeindex |= 4;
        if (v[3]) cubeindex |= 8;
        if (v[4]) cubeindex |= 16;
        if (v[5]) cubeindex |= 32;
        if (v[6]) cubeindex |= 64;
        if (v[7]) cubeindex |= 128;
      
        /* Cube is entirely in/out of the surface */
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] == 0) break;

        std::array<NodeMapIt,12> vertlist;

        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 1)
          vertlist[0]  = intersections.find(Edge(pt[0],pt[1]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 2)
          vertlist[1]  = intersections.find(Edge(pt[1],pt[2]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 4)
          vertlist[2]  = intersections.find(Edge(pt[2],pt[3]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 8)
          vertlist[3]  = intersections.find(Edge(pt[3],pt[0]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 16)
          vertlist[4]  = intersections.find(Edge(pt[4],pt[5]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 32)
          vertlist[5]  = intersections.find(Edge(pt[5],pt[6]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 64)
          vertlist[6]  = intersections.find(Edge(pt[6],pt[7]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 128)
          vertlist[7]  = intersections.find(Edge(pt[7],pt[4]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 256)
          vertlist[8]  = intersections.find(Edge(pt[0],pt[4]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 512)
          vertlist[9]  = intersections.find(Edge(pt[1],pt[5]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 1024)
          vertlist[10]  = intersections.find(Edge(pt[2],pt[6]));
        if (GeomIntersectUtils::MarchingCubesEdgeTable[cubeindex] & 2048)
          vertlist[11]  = intersections.find(Edge(pt[3],pt[7]));

        /* Create the triangles */
        int nTriang = 0;
        for (int i=0;GeomIntersectUtils::MarchingCubesTriTable[cubeindex][i]!=-1;i+=3)
          nTriang++;

        for (int j=0; j<nTriang; ++j)
        {
          auto& triangleVec = surfaceFragmentsG.back();
          int j3 = 3*j;
          triangleVec.push_back(Triangle(vertlist[GeomIntersectUtils::MarchingCubesTriTable[cubeindex][j3  ]],
                                         vertlist[GeomIntersectUtils::MarchingCubesTriTable[cubeindex][j3+1]],
                                         vertlist[GeomIntersectUtils::MarchingCubesTriTable[cubeindex][j3+2]]));
        }
      }
#endif
    }

    // Write contour fragments to surface file
    amrex::Print() << "EBISLevel: Writing EB surface to: " << eb_surface_filename << '\n';

    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(eb_surface_filename, 0755))
            amrex::CreateDirectoryFailed(eb_surface_filename);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();


    for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
    {
      auto& intersections = m_intersections[mfi];
      auto& fragVec = surfaceFragmentsG[mfi.LocalIndex()];
      if (fragVec.size() > 0)
      {
        Vector<NodeMapIt> orderedNodes(intersections.size());
        for (NodeMapIt it=intersections.begin(); it!=intersections.end(); ++it)
        {
          orderedNodes[it->first.ID] = it;
        }

        int SizeOfEdgeData = 2 * AMREX_SPACEDIM;
        int SizeOfNodeData = AMREX_SPACEDIM;
        int SizeOfFragData = AMREX_SPACEDIM;

        Vector<int>  flattenedEdges(SizeOfEdgeData * orderedNodes.size());
        Vector<Real> flattenedNodes(SizeOfNodeData * orderedNodes.size());
        Vector<long> flattenedFrags(SizeOfFragData * fragVec.size());

        int*  feptr = flattenedEdges.dataPtr();
        Real* fnptr = flattenedNodes.dataPtr();
        long* ffptr = flattenedFrags.dataPtr();

        for (int i=0; i<orderedNodes.size(); ++i)
        {
          auto eoffset = i*SizeOfEdgeData;
          auto noffset = i*SizeOfNodeData;
          const auto& edge = orderedNodes[i]->first;
          const auto& node = orderedNodes[i]->second;
          AMREX_ASSERT(node.size() == SizeOfNodeData);
          auto& L = edge.IV_l;
          auto& R = edge.IV_r;
          for (int d=0; d<AMREX_SPACEDIM; ++d)
          {
            feptr[eoffset+d               ] = L[d];
            feptr[eoffset+d+AMREX_SPACEDIM] = R[d];
            fnptr[noffset+d] = node[d];
          }
        }
        for (int i=0; i<fragVec.size(); ++i)
        {
          auto foffset = i*SizeOfFragData;
          auto frag = fragVec[i];
          AMREX_ASSERT(frag.size() == SizeOfFragData);
          for (int d=0; d<SizeOfFragData; ++d)
          {
            ffptr[foffset+d] = frag[d]->first.ID;
          }
        }

        auto FullDataPath = eb_surface_filename;
        if (!FullDataPath.empty() && FullDataPath[FullDataPath.size()-1] != '/')
          FullDataPath += '/';
        FullDataPath += "Data";

        auto nGrid = m_intersections.size();
        auto nDigits = std::log10(nGrid) + 1;
        FullDataPath += Concatenate("_",mfi.index(),nDigits);

        std::ofstream ofs(FullDataPath.c_str(),std::ios::binary);
        writeIntData( flattenedEdges.dataPtr(), flattenedEdges.size(), ofs);
        writeRealData(flattenedNodes.dataPtr(), flattenedNodes.size(), ofs);
        writeLongData(flattenedFrags.dataPtr(), flattenedFrags.size(), ofs);

        ofs.close();
      }
    }

    int nGrids = m_intersections.size();
    int ioProc = ParallelDescriptor::IOProcessorNumber();

    Vector<long> nNodes(nGrids,0);
    Vector<long> nElts(nGrids,0);

    for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
    {
      nNodes[mfi.index()] = m_intersections[mfi].size();
      nElts[mfi.index()] = surfaceFragmentsG[mfi.LocalIndex()].size();
    }

    ParallelDescriptor::ReduceLongSum(nNodes.dataPtr(), nNodes.size(), ioProc);
    ParallelDescriptor::ReduceLongSum(nElts.dataPtr(),  nElts.size(),  ioProc);

    std::string FullHeaderPath = eb_surface_filename;
    if (ParallelDescriptor::IOProcessor())
    {
      //std::string FullHeaderPath = eb_surface_filename;
      if (!FullHeaderPath.empty() && FullHeaderPath[FullHeaderPath.size()-1] != '/')
        FullHeaderPath += '/';
      FullHeaderPath += "Header";

      // Write header info
      std::ofstream ofs(FullHeaderPath.c_str());
      ofs << "EBsurfaceFormat-V1" << '\n';
      ofs << AMREX_SPACEDIM << '\n'; // number of nodes per element
      ofs << AMREX_SPACEDIM << '\n'; // number of components of floating point data at each node

      long nDataFiles = 0;
      for (int i=0; i<nGrids; ++i) {
        if (nElts[i] != 0) nDataFiles++;
      }
      ofs << nDataFiles << '\n';
      for (int i=0; i<nGrids; ++i) {
        if (nElts[i] != 0) ofs << i<< " ";
      }
      ofs << '\n';
      for (int i=0; i<nGrids; ++i) {
        if (nElts[i] != 0) ofs << nNodes[i] << " ";
      }
      ofs << '\n';
      for (int i=0; i<nGrids; ++i) {
        if (nElts[i] != 0) ofs << nElts[i] << " ";
      }
      ofs << '\n';
      ofs << FPC::NativeIntDescriptor() << '\n';
      ofs << FPC::NativeRealDescriptor() << '\n';
      ofs << FPC::NativeLongDescriptor() << '\n';
      ofs.close();
    }
  }

  ///
  void
  EBISLevel::defineFromGeometryService(const GeometryService & a_geoserver)
  {

    BL_PROFILE("EBISLevel::defineFromGeometryService()");   
 
    m_grids.define(m_domain);
    m_grids.maxSize(m_nCellMax);
    m_dm.define(m_grids);
    int ngrowGraph =2;
    int ngrowData =0;
    m_graph.define(m_grids, m_dm, 1, ngrowGraph, MFInfo(), DefaultFabFactory<EBGraph>());

    m_intersections.define(m_grids, m_dm);

    LayoutData<Vector<IrregNode> > allNodes(m_grids, m_dm);


    for (MFIter mfi(m_grids, m_dm); mfi.isValid(); ++mfi)
    {
      const Box& valid  = mfi.validbox();
      Box ghostRegion = valid;
      ghostRegion.grow(ngrowGraph);
      Box ghostRegionInt = ghostRegion;
      ghostRegionInt &= m_domain;

      EBGraph& ebgraph = m_graph[mfi];
      GeometryService::InOut inout = a_geoserver.InsideOutside(ghostRegionInt, m_domain, m_origin, m_dx);
      ebgraph.setDomain(m_domain);
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
        Vector<IrregNode>&   nodes = allNodes[mfi];
        NodeMap&        intersects = m_intersections[mfi];

        a_geoserver.fillGraph(regIrregCovered, nodes, intersects, valid,
                              ghostRegionInt, m_domain,
                              m_origin, m_dx);

        if (!m_build_eb_surface)
        {
          intersects.clear();
        }

        ebgraph.buildGraph(regIrregCovered, nodes, ghostRegion, m_domain);

      }
    }

    if (m_build_eb_surface)
    {
      buildEBSurface(a_geoserver);

      for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
      {
        m_intersections[mfi].clear();
      }
    }

    m_graph.FillBoundary();

//begin debug
//    EBISLevel_checkGraph(m_grids, m_dm, m_graph, string("after initial build"));
// end debug

    std::shared_ptr<FabArray<EBGraph> > graphptr(&m_graph, &null_deleter_fab_ebg);
    EBDataFactory ebdf(graphptr);

    m_data.define(m_grids, m_dm, 1, ngrowData, MFInfo(), ebdf);

    for (MFIter mfi(m_grids, m_dm); mfi.isValid(); ++mfi)
    {
      const Box& valid  = mfi.validbox();
      Box ghostRegion = valid;
      ghostRegion.grow(ngrowData);
      ghostRegion &= m_domain;

      const EBGraph& ebgraph = m_graph[mfi];
      EBData& ebdata   = m_data [mfi];
      if (ebgraph.isAllRegular() || ebgraph.isAllCovered())
      {
        
        ebdata.define(ebgraph,  ghostRegion);
      }
      else
      {
        const Vector<IrregNode>&   nodes = allNodes[mfi];
        ebdata.define(ebgraph, nodes, valid, ghostRegion);
      }
    }

//begin debug
//    EBISLevel_checkData(m_grids, dm, m_data, string("after initial build"));
// end debug

  }
  ///
  EBISLevel::
  EBISLevel(EBISLevel             & a_fineEBIS,
            const GeometryService & a_geoserver)
      : m_nCellMax (a_fineEBIS.m_nCellMax),
        m_domain   (amrex::coarsen(a_fineEBIS.m_domain, 2)),
        m_origin   (a_fineEBIS.m_origin),
        m_dx       (2.*a_fineEBIS.m_dx)
  { // method used by EBIndexSpace::buildNextLevel
    BL_PROFILE("EBISLevel::EBISLevel_fineEBIS");

    m_grids.define(m_domain);
    m_grids.maxSize(m_nCellMax);

    m_dm.define(m_grids);

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

    int nghostGraph = 2;
    int nghostData  = 0;
    int srcGhost = 0;
  
    //need two because of coarsen faces
    m_graph.define(m_grids, m_dm, 1, nghostGraph, MFInfo(), DefaultFabFactory<EBGraph>());
    FabArray<EBGraph> ebgraphReCo(gridsReCo, m_dm, 1, nghostGraph+1,
                                  MFInfo(), DefaultFabFactory<EBGraph>());
    //pout() << "ebislevel::coarsenvofsandfaces: doing ebgraph copy" << endl;

//begin debug
//    EBISLevel_checkGraph(a_fineEBIS.m_grids, a_fineEBIS.m_dm, a_fineEBIS.m_graph, string(" source graph for copy"));
//end debug

    ebgraphReCo.copy(a_fineEBIS.m_graph, 0, 0, 1, srcGhost, nghostGraph+1);

//begin debug
//    EBISLevel_checkGraph( gridsReCo, m_dm, ebgraphReCo, string(" ebgraphReCo after copy "));
//end debug

    ///first deal with the graph
    //pout() << "ebislevel::coarsenvofsandfaces: doing coarsenvofs " << endl;
    for (MFIter mfi(m_grids, m_dm); mfi.isValid(); ++mfi)
    {
      EBGraph      & fineEBGraph = ebgraphReCo[mfi];
      EBGraph      & coarEBGraph = m_graph[mfi];
      const Box    & coarRegion  = mfi.validbox();
      //Box coarRegion2 = m_grids[mfi];
      coarEBGraph.coarsenVoFs(fineEBGraph, coarRegion);
    }

    //pout() << "ebislevel::coarsenvofsandfaces: doing finetocoarse " << endl;
    for (MFIter mfi(m_grids, m_dm); mfi.isValid(); ++mfi)
    {
      EBGraph      & fineEBGraph = ebgraphReCo[mfi];
      EBGraph      & coarEBGraph = m_graph[mfi];
      const Box    & coarRegion  = mfi.validbox();
      //pout() << "coarsen faces for box" << coarRegion << endl;
      coarEBGraph.coarsenFaces(fineEBGraph, coarRegion);
      //pout() << "fixFineToCoarse for box" << coarRegion << endl;
      coarEBGraph.fixFineToCoarse(fineEBGraph, coarRegion);
    }
    //after fixing up fine to coarse, copy info back
    //pout() << "ebislevel::doing copy back " << endl;
    a_fineEBIS.m_graph.copy(ebgraphReCo, 0, 0, 1, 0, 0);
    a_fineEBIS.m_graph.FillBoundary();

    //pout() << "out of copyback and making new holders" << endl;

    //now deal with the data
    std::shared_ptr<FabArray<EBGraph> > graphptrCoar(&    m_graph, &null_deleter_fab_ebg);
    std::shared_ptr<FabArray<EBGraph> > graphptrReCo(&ebgraphReCo, &null_deleter_fab_ebg);
    EBDataFactory ebdfCoar(graphptrCoar);
    EBDataFactory ebdfReCo(graphptrReCo);
    FabArray<EBData> ebdataReCo;

    //pout() << "making m_data" << endl;
    m_data    .define(m_grids  , m_dm, 1, nghostData, MFInfo(), ebdfCoar);

    //pout() << "making ebdataReCo" << endl;
    ebdataReCo.define(gridsReCo, m_dm, 1, 0, MFInfo(), ebdfReCo);

//begin debug
//   EBISLevel_checkData(a_fineEBIS.m_grids, a_fineEBIS.m_dm, a_fineEBIS.m_data, string(" source data for copy"));
//end debug

    //pout() << "doing ebdatareco copy" << endl;
    ebdataReCo.copy(a_fineEBIS.m_data, 0, 0, 1, 0, 0);    

//begin debug
    //EBISLevel_checkData(gridsReCo, m_dm, ebdataReCo, string(" ebdataReCo after copy "));
//end debug

    //pout() << "coarsening data" << endl;
    for (MFIter mfi(m_grids, m_dm); mfi.isValid(); ++mfi)
    {
      const EBGraph& fineEBGraph = ebgraphReCo[mfi];
      const EBGraph& coarEBGraph =     m_graph[mfi];
      const EBData & fineEBData  =  ebdataReCo[mfi];
      const Box& valid = mfi.validbox();
      m_data[mfi].coarsenVoFs (fineEBData, fineEBGraph, coarEBGraph, valid);
      m_data[mfi].coarsenFaces(fineEBData, fineEBGraph, coarEBGraph, valid);
    }
//    m_data.FillBoundary();
  }

  //now fix the multivalued next to regular thing for the graph and the data
  //the oldgraph/newgraph thing is necessary because the graphs are
  //reference counted and they have to be kept consistent with the data
  void 
  EBISLevel::
  fixRegularNextToMultiValued()
  {
    BL_PROFILE("EBISLevel::fixRegularNextToMultiValued");

    for (MFIter mfi(m_grids, m_dm); mfi.isValid(); ++mfi)
    {
      IntVectSet vofsToChange;
      const Box& valid = mfi.validbox();
      m_graph[mfi].getRegNextToMultiValued(vofsToChange, valid);
      m_graph[mfi].addFullIrregularVoFs(vofsToChange);
      m_data[ mfi].addFullIrregularVoFs(vofsToChange, valid);

    }
    m_data .FillBoundary();
    m_graph.FillBoundary();
  }

  ///
  void 
  EBISLevel::
  fillEBISLayout(EBISLayout     & a_ebisLayout,
                 const BoxArray & a_grids,
                 const DistributionMapping & a_dm,
                 const int      & a_nghost) const
  {
    BL_ASSERT(a_nghost >= 0);
  
    //a_ebisLayout.define(m_domain, a_grids, a_nghost, m_graph, m_data);
    //return; // caching disabled for now.... ugh.  bvs
    a_ebisLayout.define(m_domain, a_grids, a_dm, a_nghost, m_graph, m_data);
  }


}

