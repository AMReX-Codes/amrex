#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "parstream.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "IndexTM.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "AMRIO.H"

#include "EBCFCopy.H"
#include "EBISLevel.H"
#include "EBIndexSpace.H"
#include "EBGraphFactory.H"
#include "EBDataFactory.H"
#include "EBISLayout.H"
#include "VoFIterator.H"
#include "IrregNode.H"
#include "AllRegularService.H"
#include "PolyGeom.H"
#include "EBLevelDataOps.H"
#include "FaceIterator.H"
#include "NamespaceHeader.H"


EBIndexSpace* Chombo_EBIS::s_instance = NULL;
bool          Chombo_EBIS::s_aliased  = false;
int EBISLevel::s_ebislGhost = 6;
EBIndexSpace* Chombo_EBIS::instance()
{
  if ((!s_aliased) && (s_instance == NULL))
  {
    s_instance = new EBIndexSpace();
  }

  return  s_instance;
}
////
void 
EBISLevel::
checkGraph() const
{
#ifndef NDEBUG
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & graph  = m_graph[dit()];
      const Box     & grid   = m_grids[dit()];
      IntVectSet ivsgrid(grid);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(FaceIterator faceit(ivsgrid, graph, idir, FaceStop::SurroundingNoBoundary); faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();
              for(SideIterator sit; sit.ok(); ++sit)
                {
                  const IntVect& iv = face.gridIndex(sit());
                  if(graph.getRegion().contains(iv) && graph.isCovered(iv))
                    {
                      pout() << "cell " << iv << "is both covered and part of a face" << endl;
                      MayDay::Error("inconsistent graph description");
                    }
                }
            }
        }
    }
#endif
}
////
void Chombo_EBIS::alias(const EBIndexSpace* a_input)
{
  s_instance = (EBIndexSpace*)(a_input);
  s_aliased  = true;
}

Real EBISLevel::s_tolerance = 1.0e-12;
bool EBISLevel::s_verbose   = false;
bool EBISLevel::s_recursive = false;

long long EBISLevel::numVoFsOnProc() const
{
  CH_TIME("EBISLevel::numVoFsOnProc");
  long long retval = 0;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& ebgraph = m_graph[dit()];
      long long numVoFsBox = ebgraph.numVoFs(m_grids.get(dit()));
      retval += numVoFsBox;
    }

  return retval;
}
///
bool isPowerOfTwo (int x)
{
  while (((x % 2) == 0) && (x > 1)) /* While x is even and > 1 */
    {
      x /= 2;
    }
 return (x == 1);
}
///
void EBISLevel::makeLoads(Vector<unsigned long long>&       a_loads,
                          Vector<Box>&                      a_boxes,
                          const Box&                        a_region,
                          const ProblemDomain&              a_domain,
                          const GeometryService&            a_geoserver,
                          const RealVect&                   a_origin,
                          const Real&                       a_dx,
                          const int                         a_ncellmax,
                          const EBIndexSpace* const         a_ebisPtr)
{
#ifdef CH_MPI
  if(EBIndexSpace::s_useMemoryLoadBalance)
    {
      pout() << "using memory for load balance" << endl;
      for (int i = 0; i < a_boxes.size(); i++)
        {
          if (a_boxes[i].ixType() ==  IndexType::TheNodeType())
            {
              a_boxes[i].convert(IndexType::TheCellType());
            }
        }
      Vector<int> procs;
      LoadBalance(procs, a_boxes);
      DisjointBoxLayout dbl(a_boxes, procs);
      DataIterator dit = dbl.dataIterator();

      dit.enablePeak();
      dit.clearPeak();
      EBGraphFactory graphfact(a_domain);
      EBDataFactory   datafact;
      LevelData<EBGraph> graph(dbl, 1, IntVect::Unit, graphfact);
      LayoutData<Vector<IrregNode> > allNodes(dbl);
      defineGraphFromGeo(graph, allNodes, a_geoserver, dbl,      
                         a_domain,  a_origin, a_dx, false);

      LevelData<EBData>   data(dbl, 1, IntVect::Zero,  datafact);
      for (dit.reset(); dit.ok(); ++dit)
        {
          data[dit()].defineVoFData( graph[dit()], dbl.get(dit()));
          data[dit()].defineFaceData(graph[dit()], dbl.get(dit()));
        }
      dit.disablePeak();
      dit.mergePeak();
      a_loads = dit.getPeak();
    }
  else
#endif
    {
      pout() << "using old ebis load balance" << endl;
      for (int i = 0; i < a_boxes.size(); i++)
        {
          if (a_boxes[i].ixType() ==  IndexType::TheNodeType())
            {
              a_boxes[i].convert(IndexType::TheCellType());
              a_loads[i] = 8;
            }
          else
            {
              a_loads[i] = 1;
            }
        }
    }
}
///
void EBISLevel::makeBoxes(Vector<Box>&               a_boxes,
                          Vector<unsigned long long>&              a_loads,
                          const Box&                 a_region,
                          const ProblemDomain&       a_domain,
                          const GeometryService&     a_geoserver,
                          const RealVect&            a_origin,
                          const Real&                a_dx,
                          const int                  a_ncellmax,
                          const EBIndexSpace* const  a_ebisPtr)
{
  bool allPowersOfTwo = true;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      bool powerOfTwoThisDir = isPowerOfTwo(a_domain.size(idir));
      allPowersOfTwo = allPowersOfTwo && powerOfTwoThisDir;
    }
  if(allPowersOfTwo && s_recursive) //the recursive makeboxes really only likes powers of two
    {
      pout() << "EBISLevel::makeBoxes -- doing recursive" << endl;
      std::list<Box> boxes;
      makeBoxes(boxes, a_region, a_domain, a_geoserver, a_origin, a_dx, a_ncellmax);
      a_boxes.resize(boxes.size());
      a_loads.resize(boxes.size());
      std::list<Box>::iterator it = boxes.begin();
      for (int i = 0; i < a_boxes.size(); ++i, ++it)
        {
          a_boxes[i]=*it;
        }
      mortonOrdering(a_boxes);
    }
  else
    {
      if (EBIndexSpace::s_MFSingleBox)
        {
          pout() << "EBISLevel::makeBoxes -- doing single box" << endl;
          a_boxes.resize(1);
          a_loads.resize(1);
          a_boxes[0] = a_region;
          a_loads[0] = 1;
          return;
        }
      pout() << "EBISLevel::makeBoxes -- doing domain split" << endl;
      domainSplit(a_domain, a_boxes, a_ncellmax, 1);
      mortonOrdering(a_boxes);
      a_loads.resize(a_boxes.size(), 1);
    }
  makeLoads(a_loads, a_boxes, a_region, a_domain, a_geoserver, a_origin, a_dx, a_ncellmax, a_ebisPtr);
}

void EBISLevel::makeBoxes(std::list<Box>&        a_boxes,
                          const Box&             a_region,
                          const ProblemDomain&   a_domain,
                          const GeometryService& a_geoserver,
                          const RealVect&        a_origin,
                          const Real&            a_dx,
                          const int              a_ncellmax)
{
  int longdir;
  int length = a_region.longside(longdir);

  if (length > a_ncellmax)
    {
      int n = length/2;
      //CH_assert(n*2==length);
      Box low(a_region), high;
      high = low.chop(longdir, a_region.smallEnd(longdir)+n);
      makeBoxes(a_boxes, low,  a_domain,
                a_geoserver, a_origin, a_dx, a_ncellmax);
      makeBoxes(a_boxes, high, a_domain,
                a_geoserver, a_origin, a_dx, a_ncellmax);
    }
  else
    {
      if (a_geoserver.InsideOutside(a_region, a_domain, a_origin, a_dx) == GeometryService::Irregular)
        {
          Box n = a_region;
          n.convert(IndexType::TheNodeType());
          a_boxes.push_back(n);
        }
      else
        {
          a_boxes.push_back(a_region);
        }
    }
}

#ifdef CH_USE_HDF5
EBISLevel::EBISLevel(HDF5Handle& a_handle)
{
  CH_TIME("EBISLevel::EBISLevel_hdf5");
  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;


  HDF5HeaderData header;
  header.readFromFile(a_handle);
  m_origin = header.m_realvect["EBIS_origin"];
  m_domain = header.m_box     ["EBIS_domain"];
  m_dx =     header.m_real    ["EBIS_dx"]    ;
  m_tolerance = m_dx*1E-4;
  m_level = 0;

  //read in the grids
  Vector<Box> boxes;
  read(a_handle,boxes);
  Vector<int> procAssign;
  LoadBalance(procAssign, boxes);
  m_grids.define(boxes, procAssign);//this should use m_domain for periodic...
  EBGraphFactory graphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, graphfact);

  //read the graph  in from the file
  std::string graphName("EBIS_graph");
  int eekflag = read(a_handle, m_graph, graphName, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }

  //need a ghosted layout so that the data can be defined properly
  LevelData<EBGraph> ghostGraph(m_grids, 1, IntVect::Unit, graphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, ghostGraph, interv);

  //now the data for the graph
  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit()].defineVoFData(ghostGraph[dit()], m_grids.get(dit()));
      m_data[dit()].defineFaceData(ghostGraph[dit()], m_grids.get(dit()));
    }
  //read the data  in from the file
  std::string  dataName("EBIS_data");

  eekflag = read(a_handle, m_data ,  dataName, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
#if 0
  pout() << "EBISLevel::EBISLevel 1 - m_grids - m_dx: " << m_dx << endl;
  pout() << "--------" << endl;
  pout() << m_grids.boxArray().size() << endl;
  pout() << "--------" << endl;
  pout() << endl;
#endif
}

void EBISLevel::write(HDF5Handle& a_handle) const
{
  CH_TIME("EBISLevel::write");
  HDF5HeaderData header;
  //this naming stuff kinda depends on the fact
  //that we are only outputting the finest level.
  //we could be slick and incorporate the
  //level number in there if we wanted.
  header.m_int["num_levels"] = 1;
  header.m_int["num_components"] = 1;
  header.m_string["component_0"] = "phi0";
  header.m_realvect["EBIS_origin"] = m_origin;
  header.m_box     ["EBIS_domain"] = m_domain.domainBox();
  header.m_real    ["EBIS_dx"]     = m_dx;
  header.writeToFile(a_handle);
  //write the grids to the file
  CH_XD::write(a_handle, m_grids);

  std::string graphName("EBIS_graph");
  std::string  dataName("EBIS_data");
  int eekflag = CH_XD::write(a_handle, m_graph, graphName);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }
  eekflag = CH_XD::write(a_handle, m_data ,  dataName);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
}
#endif
///
void
EBISLevel::defineGraphFromGeo(LevelData<EBGraph>             & a_graph,
                              LayoutData<Vector<IrregNode> > & a_allNodes,
                              const GeometryService          & a_geoserver,
                              const DisjointBoxLayout        & a_grids,
                              const ProblemDomain            & a_domain,
                              const RealVect                 & a_origin,
                              const Real                     & a_dx,
                              const bool                     & a_distributedData)
{
  //define the graph stuff
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box region = a_grids.get(dit());
      region.grow(1);
      Box ghostRegion = grow(region,1);
      ghostRegion &= a_domain;
      region &= a_domain;

      EBGraph& ebgraph = a_graph[dit()];
      GeometryService::InOut inout;
      if (!a_distributedData)
        {
          inout = a_geoserver.InsideOutside(region, a_domain, a_origin, a_dx);
        }
      else
        {
          inout = a_geoserver.InsideOutside(region, a_domain, a_origin, a_dx, dit());
        }
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
          BaseFab<int>       regIrregCovered(ghostRegion, 1);
          Vector<IrregNode>&  nodes = a_allNodes[dit()];

          if (!a_distributedData)
            {
              a_geoserver.fillGraph(regIrregCovered, nodes, region,
                                    ghostRegion, a_domain,
                                    a_origin, a_dx);
            }
          else
            {
              a_geoserver.fillGraph(regIrregCovered, nodes, region,
                                    ghostRegion, a_domain,
                                    a_origin, a_dx, dit());
            }
          ebgraph.buildGraph(regIrregCovered, nodes, region, a_domain);

        }
    }
}

EBISLevel::EBISLevel(const ProblemDomain   & a_domain,
                     const RealVect        & a_origin,
                     const Real            & a_dx,
                     const GeometryService & a_geoserver,
                     const EBIndexSpace    * const a_ebisPtr,
                     const bool            & a_distributedData,
                     const bool            & a_fixRegularNextToMultiValued,
                     const int             & a_level)
{
  // this is the method called by EBIndexSpace::buildFirstLevel
  CH_TIME("EBISLevel::EBISLevel_geoserver_domain");
  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;

  m_domain = a_domain;
  m_dx = a_dx;
  m_tolerance = a_dx*1E-4;
  m_origin = a_origin;

  m_level = a_level;

  if (!a_distributedData)
    { // this is the original code

      //divide up the domain into a layout
      Vector<Box> vbox;
      Vector<unsigned long long> irregCount;
      {
        CH_TIME("EBISLevel::EBISLevel_makeboxes");
        makeBoxes(vbox,
                  irregCount,
                  a_domain.domainBox(),
                  a_domain,
                  a_geoserver,
                  a_origin,
                  a_dx,
                  a_ebisPtr->getNCellMax(),
                  a_ebisPtr);
      }

      // pout()<<vbox<<"\n\n";
      //load balance the boxes
      Vector<int> procAssign;
      UnLongLongLoadBalance(procAssign, irregCount, vbox);
      //   pout()<<irregCount<<std::endl;
      //   pout()<<procAssign<<std::endl;
      m_grids.define(vbox, procAssign,a_domain);//this should use a_domain for periodic

    }
  else
    {
      // permit the geometry service to construct a layout
      int nCellMax = a_ebisPtr->getNCellMax();
      (const_cast<GeometryService*>(&a_geoserver))->makeGrids(a_domain, m_grids, nCellMax, 15);
    }

  RealVect dx2D;
  for (int i = 0; i < SpaceDim; i++)
    {
      dx2D[i]=a_dx;
    }

  (const_cast<GeometryService*>(&a_geoserver))->postMakeBoxLayout(m_grids,dx2D);
  LayoutData<Vector<IrregNode> > allNodes(m_grids);

  EBGraphFactory graphfact(a_domain);
  m_graph.define(m_grids, 1, IntVect::Unit, graphfact);

  defineGraphFromGeo(m_graph, allNodes, a_geoserver, m_grids,
                     m_domain,m_origin, m_dx, 
                     a_distributedData);

  checkGraph();

  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit()].define(m_graph[dit()], allNodes[dit()], m_grids.get(dit()));

    }

  if(a_geoserver.canGenerateMultiCells())
    {
      if (a_fixRegularNextToMultiValued)
        {
          fixRegularNextToMultiValued();
        }
    }
}

//now fix the multivalued next to regular thing for the graph and the data
//the oldgraph/newgraph thing is necessary because the graphs are
//reference counted and they have to be kept consistent with the data
void EBISLevel::fixRegularNextToMultiValued()
{
  CH_TIME("EBISLevel::fixRegularNextToMultiValued");

  EBGraphFactory graphfact(m_domain);
  LayoutData<IntVectSet> vofsToChange(m_grids);
  LevelData<EBGraph> oldGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
  Interval interv(0,0);

  m_graph.copyTo(interv, oldGhostGraph, interv);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      CH_TIME("EBISLevel::fixRegularNextToMultiValued_loop1");
      m_graph[dit()].getRegNextToMultiValued(vofsToChange[dit()],
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
      CH_TIME("EBISLevel::fixRegularNextToMultiValued_loop2");

      Box localBox = m_grids.get(dit());
      localBox.grow(1);
      localBox &= m_domain;
      newGhostData[dit()].defineVoFData(oldGhostGraph[dit()],  localBox);
      newGhostData[dit()].defineFaceData(oldGhostGraph[dit()], localBox);
    }

  m_data.copyTo(interv,  newGhostData,  interv);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      CH_TIME("EBISLevel::fixRegularNextToMultiValued_loop3");

      m_data[dit() ].addFullIrregularVoFs(vofsToChange[dit()],
                                          newGhostGraph[dit()],
                                          newGhostData[dit()].getVolData(),
                                          oldGhostGraph[dit()]);
    }
}

//checks to see the vofs are in the correct cells.
//checks to see that the faces are over the correct cells
//checks that volume fractions, area fractions are positive
//bail out with MayDay::Error if anything fails
void EBISLevel::sanityCheck(const EBIndexSpace* const a_ebisPtr)
{
#if 0
  pout() << "EBISLevel::sanityCheck" << endl;
  CH_TIME("EBISLevel::sanityCheck");
  EBISLayout ghostLayout;

  a_ebisPtr->fillEBISLayout(ghostLayout, m_grids, m_domain, 1);

  Real maxcentval = 0.5+ s_tolerance;
  for (DataIterator dit = m_grids.dataIterator();  dit.ok(); ++dit)
    {
      const Box& thisBox = m_grids.get(dit());
      const EBISBox& ebisBox = ghostLayout[dit()];
      for (BoxIterator bit(thisBox); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          for (int ivof = 0; ivof < vofs.size(); ivof++)
            {
              const VolIndex& vof = vofs[ivof];
              if (vof.gridIndex() != iv)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has grid index = " << vof.gridIndex() << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 1");
                }
              if (vof.cellIndex() < 0)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has negative cell index = " << vof.cellIndex() << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 2");
                }
              Real volFrac = ebisBox.volFrac(vof);
              if (volFrac < 0.0)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has invalid volume fraction = " << volFrac << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 5");
                }
              RealVect volCentroid = ebisBox.centroid(vof);
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Real volcentdir = volCentroid[idir];
                  if (volFrac > s_tolerance)
                    {
                      if (volcentdir > maxcentval || volcentdir < -maxcentval)
                        {
                          pout() << "EBISLevel::sanityCheck: Error" << endl;
                          pout() << "VoF at Intvect = " << iv
                                 << " has invalid vol centroid = " << volcentdir
                                 << " at direction "<< idir << endl;
                          MayDay::Error("EBISLevel::sanityCheck Error 51");
                        }
                    }
                }
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
                      IntVect iv2 = iv + sign(sit())*BASISV(idir);
                      ////check for regular next to covered and multivalued next to regular
                      if (m_domain.contains(iv2))
                        {
                          if (ebisBox.isRegular(iv))
                            {
                              if (ebisBox.isCovered(iv2))
                                {
                                  pout() << iv << " is regular and " <<  iv2 << " is covered" << endl;
                                  MayDay::Error("EBISLevel::sanityCheck error 420 ");
                                }
                              else
                                {
                                  Vector<VolIndex> otherVoFs = ebisBox.getVoFs(iv2);
                                  if (otherVoFs.size() > 1)
                                    {
                                      pout() << iv << " is regular and " <<  iv2 << " is multivalued" << endl;
                                      MayDay::Error("EBISLevel::sanityCheck error 420.2 ");
                                    }
                                }
                            }
                        }
                      IntVect ivlo, ivhi;
                      if (sit() == Side::Lo)
                        {
                          ivlo = iv2;
                          ivhi = iv;
                        }
                      else
                        {
                          ivlo = iv;
                          ivhi = iv2;
                        }
                      for (int iface = 0; iface < faces.size(); iface++)
                        {
                          const FaceIndex& face = faces[iface];
                          if (face.gridIndex(Side::Lo) != ivlo)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "face at IntVects = " << ivlo << "  " << ivhi
                                     << "has low IntVect  = " << face.gridIndex(Side::Lo)
                                     << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 3");
                            }
                          if (face.gridIndex(Side::Hi) != ivhi)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "face at IntVects = " << ivlo << "  " << ivhi
                                     << "has high IntVect = " << face.gridIndex(Side::Hi)
                                     << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 4");
                            }
                          Real areaFrac = ebisBox.areaFrac(face);
                          if (areaFrac  < 0.0)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "VoF at Intvect = " << iv
                                     << "has invalid area fraction = " << areaFrac << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 51");
                            }
                          if (areaFrac  >  s_tolerance)
                            {
                              RealVect faceCentroid = ebisBox.centroid(face);
                              for (int idir = 0; idir < SpaceDim; idir++)
                                {
                                  Real facecentdir = faceCentroid[idir];
                                  if (facecentdir > maxcentval || facecentdir < -maxcentval)
                                    {
                                      pout() << "EBISLevel::sanityCheck: Error" << endl;
                                      pout() << "VoF at Intvect = " << iv
                                             << " has invalid face centroid = " << facecentdir
                                             << " at direction "<< idir << endl;
                                      MayDay::Error("EBISLevel::sanityCheck Error 51");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif
}

EBISLevel::EBISLevel()
{
  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;

  m_level = 0;

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
void EBISLevel::dumpDebug(const string& a_string)
{
  if (m_domain.domainBox() == EBGraphImplem::s_doDebug)
    {
      pout() << a_string <<  ": ";
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          if (m_grids[dit()].contains(EBGraphImplem::s_ivDebug))
            {
              pout() << "EBIS1: " << EBGraphImplem::s_ivDebug;
              if (m_graph[dit()].isRegular(EBGraphImplem::s_ivDebug))
                {
                  pout() << " is regular" << endl;
                }
              else if (m_graph[dit()].isCovered(EBGraphImplem::s_ivDebug))
                {
                  pout() << " is covered" << endl;
                }
              else
                {
                  pout() << " is irregular" << endl;
                }
            }
        }

    }
}

void EBISLevel::coarsenVoFs(EBISLevel& a_fineEBIS)
{
  CH_TIME("EBISLevel::coarsenVoFs");

  //so that i can do the vofs and faces of this level.
  DisjointBoxLayout fineFromCoarDBL;
  refine(fineFromCoarDBL, m_grids, 2);
  fineFromCoarDBL.close();

  //no need for ghost cells here except to define the face data
  //you need the graph to be one bigger
  EBGraphFactory ebgraphfactFine(a_fineEBIS.m_domain);
  LevelData<EBGraph> fineFromCoarEBGraph(fineFromCoarDBL,1, IntVect::Unit, ebgraphfactFine);

  Interval interv(0,0);
  a_fineEBIS.m_graph.copyTo(interv, fineFromCoarEBGraph, interv);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& fineEBGraph = fineFromCoarEBGraph[dit()];
      const Box& coarRegion      = m_grids.get(dit());
      EBGraph& coarEBGraph = m_graph[dit()];
      coarEBGraph.coarsenVoFs(fineEBGraph, coarRegion);
    }

  EBGraphFactory ebgraphfactCoar(m_domain);
  LevelData<EBGraph> coarGhostEBGraph(m_grids,1, IntVect::Unit, ebgraphfactCoar);
  m_graph.copyTo(interv, coarGhostEBGraph, interv);

  //dumpDebug(string("EBIS::coarsenVoFs"));

  EBDataFactory ebdatafact;
  LevelData<EBData> fineFromCoarEBData(fineFromCoarDBL,1, IntVect::Zero, ebdatafact);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      CH_TIME("EBISLevel::coarsenVoFs_defineData");

      const Box& localBox = fineFromCoarDBL.get(dit());
      fineFromCoarEBData[dit()].defineVoFData(fineFromCoarEBGraph[dit()],  localBox);
      fineFromCoarEBData[dit()].defineFaceData(fineFromCoarEBGraph[dit()], localBox);
    }

  a_fineEBIS.m_data.copyTo(interv, fineFromCoarEBData, interv);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& fineEBGraph =  fineFromCoarEBGraph[dit()];
      const EBData& fineEBData = fineFromCoarEBData[dit()];
      const EBGraph& coarEBGraph = coarGhostEBGraph[dit()];

      m_data[dit()].coarsenVoFs(fineEBData, fineEBGraph, coarEBGraph, m_grids.get(dit()));
    }
}

void EBISLevel::fixFineToCoarse(EBISLevel& a_fineEBIS)
{
  CH_TIME("EBISLevel::fixFineToCoarse");
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
  CH_TIME("EBISLevel::coarsenFaces");
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

EBISLevel::EBISLevel(EBISLevel             & a_fineEBIS,
                     const GeometryService & a_geoserver,
                     const EBIndexSpace    * const a_ebisPtr,
                     const bool            & a_distributedData,
                     const bool            & a_fixRegularNextToMultiValued)
{ // method used by EBIndexSpace::buildNextLevel
  CH_TIME("EBISLevel::EBISLevel_fineEBIS");

  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;

  m_domain = coarsen(a_fineEBIS.m_domain,2);
  m_dx = 2.*a_fineEBIS.m_dx;
  m_tolerance = 2.*a_fineEBIS.m_tolerance;
  m_origin = a_fineEBIS.m_origin;

  m_level = a_fineEBIS.m_level + 1;

  if (!a_distributedData)
    { // this is the original method

      Vector<Box> vbox;
      Vector<unsigned long long> irregCount;
      {
        CH_TIME("EBISLevel::EBISLevel_fineEBIS_makeboxes 2");
        makeBoxes(vbox, irregCount, m_domain.domainBox(), m_domain, a_geoserver,
                  m_origin, m_dx, a_ebisPtr->getNCellMax(), a_ebisPtr);
      }

      //pout()<<vbox<<"\n\n";
      //load balance the boxes
      Vector<int> procAssign;
      UnLongLongLoadBalance(procAssign, irregCount, vbox);

      //pout()<<procAssign<<std::endl;
      //define the layout.  this includes the domain and box stuff
      m_grids.define(vbox, procAssign);//this should use m_domain for periodic

    }
  else
    {
      // permit the geometry service to construct a layout
      int nCellMax = a_ebisPtr->getNCellMax();
      (const_cast<GeometryService*>(&a_geoserver))->makeGrids(m_domain, m_grids, nCellMax, 15);
    }

  EBGraphFactory ebgraphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, ebgraphfact);

  EBDataFactory ebdatafact;
  m_data.define(m_grids, 1, IntVect::Zero, ebdatafact);

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
  checkGraph();
#if 0
  pout() << "EBISLevel::EBISLevel 4 - m_grids - m_dx: " << m_dx << endl;
  pout() << "--------" << endl;
  pout() << m_grids.boxArray().size() << endl;
  pout() << "--------" << endl;
  pout() << endl;
#endif
}

EBISLevel::~EBISLevel()
{
}

void EBISLevel::fillEBISLayout(EBISLayout&              a_ebisLayout,
                               const DisjointBoxLayout& a_grids,
                               const int&               a_nghost) const
{
  CH_assert(a_nghost >= 0);
  
  //a_ebisLayout.define(m_domain, a_grids, a_nghost, m_graph, m_data);
  //return; // caching disabled for now.... ugh.  bvs

  EBISLayout& l = m_cache[a_grids];
  if (!l.isDefined() || (a_nghost > l.getGhost()))
    {
      CH_TIME("ebisllevel::fillebislayout cache miss");
      //int thisghost = Max(s_ebislGhost, a_nghost);
      int thisghost = a_nghost;
      l.define(m_domain, a_grids, thisghost, m_graph, m_data);
      m_cacheMisses++;
      m_cacheStale++;
      //pout()<<"a_nghost:"<<a_nghost;
    }
  else
    {
      CH_TIME("cache_hit");
      m_cacheHits++;
    }
  a_ebisLayout = l;// refcount is at least 2 now.
  if (m_cacheStale == 1)
    {
      refreshCache();
      m_cacheStale = 0;
    }
  //int isize = m_cache.size();
  //  pout()<<"fillebisl::m_level:"<<m_level<<" m_cache.size():"<<isize<<" m_cacheHits:"<<m_cacheHits<<" m_cacheMisses:"<<m_cacheMisses<<"\n";
}

void EBISLevel::dumpCache() const
{
  pout()<<std::endl;
  dmap::iterator d = m_cache.begin();
  while (d != m_cache.end())
    {
      d++;
	  
    }

} 
void EBISLevel::refreshCache() const
{
  dmap::iterator d = m_cache.begin();
  while (d != m_cache.end())
    {
      if(d->second.refCount() ==1)
	{
	  m_cache.erase(d++);
	}
      else
	{
	  d++;
	}
    }

 // int s=m_cache.size();
 // pout()<<" m_level:"<<m_level<<" m_cache.size():"<<s<<" m_cacheHits:"<<m_cacheHits<<" m_cacheMisses:"<<m_cacheMisses<<"\n";
}

void EBISLevel::clearMultiBoundaries()
{
  CH_TIME("EBISLevel::clearMultiBoundaries");
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      EBData& data = m_data[dit];
      data.clearMultiBoundaries();
    }
}

void EBISLevel::setBoundaryPhase(int phase)
{
  CH_TIME("EBISLevel::setBoundaryPhase");
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      EBData& data = m_data[dit];
      data.setBoundaryPhase(phase);
    }
}

//throughout this routine A refers to objects belonging to this
//ebindexspace,  B refers to the other one.
void setOtherVoFSingleValued(VolData&        a_vol,
                             const VolIndex& a_sourceVoF,
                             int&            a_otherPhase,
                             EBData&         a_otherFluidData,
                             bool&           a_sourceIvContainedInOtherFluid,
                             bool&           a_skipCell)
{
  VolIndex vother = a_sourceVoF;
  a_skipCell = false;
  bool fixBoundaryData = false;
  VolData refVol;
  //if a cell is full and has a unit boundary area, it means that the irregular
  //boundary lives on a cell boundary.  We have to link to a cell over.
  // In this case, in the presence of multi valued cells, to quote Yeats:
  //Mere anarchy is loosed upon the world,
  //The blood-dimmed tide is loosed, and everywhere
  //The ceremony of innocence is drowned.
  if ((a_vol.m_volFrac == 1) && (a_vol.m_averageFace.m_bndryArea == 1))
    {
      a_skipCell = true;
      //figure out which cell we are moving to. Then we use cellIndex = 0
      //because the cells involved are single valued.
      int whichWayToGo = -23;
      int sign = 1;
      bool found = false;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real normDir = a_vol.m_averageFace.m_normal[idir];
          if (Abs(normDir) == 1)
            {
              found = true;
              whichWayToGo = idir;
              sign = -1;
              if (normDir < 0)
                {
                  sign = 1;
                }
            }
        }
      if (!found)
        {
          MayDay::Error("EBIndexSpace::setOtherVoFSingleValued - Unit normal direction not found");
        }
      IntVect ivOther = sign*BASISV(whichWayToGo);
      ivOther += a_sourceVoF.gridIndex();
      vother = VolIndex(ivOther, 0);
    }
  else if ((a_vol.m_volFrac == 0) &&
           (a_vol.m_averageFace.m_bndryArea == 1) &&
           (a_sourceIvContainedInOtherFluid))
    {
      a_skipCell = true;
      refVol = a_otherFluidData.getVolData()(a_sourceVoF, 0);
      // use -1*normal of the VoF with this IV in the other fluid
      RealVect normal = -refVol.m_averageFace.m_normal;
      //figure out which cell we are moving to. Then we use cellIndex = 0
      //because the cells involved are single valued.
      int whichWayToGo = -23;
      int sign = 1;
      bool found = false;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real normDir = normal[idir];
          if (Abs(normDir) == 1)
            {
              found = true;
              whichWayToGo = idir;
              sign = 1;
              if (normDir < 0)
                {
                  sign = -1;
                }
            }
        }
      if (found)
        {
          IntVect ivOther = sign*BASISV(whichWayToGo);
          ivOther += a_sourceVoF.gridIndex();
          vother = VolIndex(ivOther, 0);
          fixBoundaryData = true;
        }
      else
        {
          MayDay::Warning("EBIndexSpace::setOtherVoFSingleValued - Unit normal direction not found; EB probably intersects cell corner");
        }
    }
  a_vol.m_phaseFaces[0].m_volIndex   = vother;
  a_vol.m_phaseFaces[0].m_bndryPhase = a_otherPhase;
  if (fixBoundaryData)
    {
      Real eps = 1e-6;
      BoundaryData&    boundaryData0 = a_vol.m_phaseFaces[0];
      BoundaryData& avgBoundaryData0 = a_vol.m_averageFace;
      BoundaryData&    boundaryData1 = refVol.m_averageFace;
      if ((boundaryData0.m_normal.vectorLength() < eps) &&
          (boundaryData1.m_normal.vectorLength() > 1-eps) &&
          (boundaryData1.m_normal.vectorLength() < 1+eps))
        {
          // fix boundary at phaseFace
          boundaryData0.m_normal = -1.0*(boundaryData1.m_normal);
          boundaryData0.m_bndryArea = boundaryData1.m_bndryArea;
          boundaryData0.m_bndryCentroid = boundaryData1.m_bndryCentroid;
          // make average face reflect new phaseFace data
          avgBoundaryData0.m_normal = boundaryData0.m_normal;
          avgBoundaryData0.m_bndryArea = boundaryData0.m_bndryArea;
          avgBoundaryData0.m_bndryCentroid = boundaryData0.m_bndryCentroid;
        }
    }
}

void testAndFixBoundaryData(VolData&       a_volData0,
                            const VolData& a_volData1)
{
  Real eps = 1e-6;
  BoundaryData&    boundaryData0 = a_volData0.m_phaseFaces[0];
  BoundaryData& avgBoundaryData0 = a_volData0.m_averageFace;
  const BoundaryData&    boundaryData1 = a_volData1.m_phaseFaces[0];
  //  const BoundaryData& avgBoundaryData1 = a_volData1.m_averageFace;
  if ((boundaryData0.m_bndryArea > 0) &&
      (boundaryData0.m_normal.vectorLength() < eps) &&
      (boundaryData1.m_normal.vectorLength() > 1-eps) &&
      (boundaryData1.m_normal.vectorLength() < 1+eps))
    {
      // fix boundary at phaseFace
      boundaryData0.m_normal = -1.0*(boundaryData1.m_normal);
      boundaryData0.m_bndryArea = boundaryData1.m_bndryArea;
      boundaryData0.m_bndryCentroid = boundaryData1.m_bndryCentroid;
      // make average face reflect new phaseFace data
      avgBoundaryData0.m_normal = boundaryData0.m_normal;
      avgBoundaryData0.m_bndryArea = boundaryData0.m_bndryArea;
      avgBoundaryData0.m_bndryCentroid = boundaryData0.m_bndryCentroid;
    }
}

void EBISLevel::reconcileIrreg(EBISLevel& a_otherPhase)
{
  CH_TIME("EBISLevel::reconcileIrreg");

  //both layouts are made independently so will not work
  //with the same data iterator but they ought to because they
  //are both breaking up the domain the same way.   This is yet
  //another subtle thing thing that will break when the
  //number of fluids is greater than two.

  //define ghosted coarse graph so i can use the faces for the ebdata
  EBGraphFactory ebgraphfactcoarA(m_domain);
  LevelData<EBGraph> ghostedEBGraphLDCoarA(m_grids, 1, IntVect::Unit,
                                           ebgraphfactcoarA);
  Interval interv(0, 0);
  m_graph.copyTo(interv, ghostedEBGraphLDCoarA, interv);

  EBGraphFactory ebgraphfactcoarB(a_otherPhase.m_domain);
  LevelData<EBGraph> ghostedEBGraphLDCoarB(a_otherPhase.m_grids,
                                           1,
                                           IntVect::Unit,
                                           ebgraphfactcoarB);
  a_otherPhase.m_graph.copyTo(interv, ghostedEBGraphLDCoarB, interv);

  DataIterator dita = m_grids.dataIterator();
  DataIterator ditb= a_otherPhase.m_grids.dataIterator();
  dita.begin(); ditb.begin();
  for ( ; dita.ok(); ++dita, ++ditb)
    {
      EBGraph& ebgrapCoarA =              m_graph[dita];
      EBGraph& ebgrapCoarB = a_otherPhase.m_graph[ditb];
      EBData&  ebdataCoarA =              m_data[dita];
      EBData&  ebdataCoarB = a_otherPhase.m_data[ditb];
      Box region = ebgrapCoarA.getRegion();

      // ghosted graphs for filling faces for ebdata
      EBGraph& ghostedEBGraphCoarA = ghostedEBGraphLDCoarA[dita];
      EBGraph& ghostedEBGraphCoarB = ghostedEBGraphLDCoarB[ditb];

      CH_assert(ebgrapCoarB.getRegion() == region);
      IntVectSet seta = ebgrapCoarA.getIrregCells(region);
      IntVectSet setb = ebgrapCoarB.getIrregCells(region);

      IntVectSet abDifference = seta - setb;
      IntVectSet baDifference = setb - seta;
      //this should only happen when one side changed a regular to an
      //irregular.   This means that the other has to change a covered
      //to an irregular
      ebgrapCoarB.addEmptyIrregularVoFs(abDifference);
      ebgrapCoarA.addEmptyIrregularVoFs(baDifference);
      // repeat add for ghosted graphs
      ghostedEBGraphCoarB.addEmptyIrregularVoFs(abDifference);
      ghostedEBGraphCoarA.addEmptyIrregularVoFs(baDifference);
      // use ghosted graphs to add to ebdata
      ebdataCoarB.addEmptyIrregularVoFs(abDifference, ghostedEBGraphCoarB);
      ebdataCoarA.addEmptyIrregularVoFs(baDifference, ghostedEBGraphCoarA);
    } //end loop over grids
}

void EBISLevel::levelStitch(EBISLevel&       a_otherPhase,
                            const EBISLevel* a_finePtrA,
                            const EBISLevel* a_finePtrB)
{
  CH_TIME("EBISLevel::levelStitch");

  //I have no idea what to do if only one of the inputs is null.
  //either both or neither makes sense
  CH_assert(((a_finePtrA != NULL) && (a_finePtrB != NULL)) ||
            ((a_finePtrA == NULL) && (a_finePtrB == NULL)));
  EBISLayout ebislFineA, ebislFineB;
  DisjointBoxLayout dblFineA, dblFineB;
  if (a_finePtrA != NULL)
    {
      int nghost = 0; //should not need any as this is all about connectivity within a coarse cell
      refine(dblFineA,              m_grids, 2);
      refine(dblFineB, a_otherPhase.m_grids, 2);
      a_finePtrA->fillEBISLayout(ebislFineA, dblFineA, nghost);
      a_finePtrB->fillEBISLayout(ebislFineB, dblFineB, nghost);
    }

  //both layouts are made independently so will not work
  //with the same data iterator but they ought to because they
  //are both breaking up the domain the same way.   This is yet
  //another subtle thing thing that will break when the
  //number of fluids is greater than two.
  DataIterator dita = m_grids.dataIterator();
  DataIterator ditb= a_otherPhase.m_grids.dataIterator();
  dita.begin(); ditb.begin();
  for ( ; dita.ok(); ++dita, ++ditb)
    {
      EBGraph& ebgrapCoarA =              m_graph[dita];
      EBGraph& ebgrapCoarB = a_otherPhase.m_graph[ditb];
      EBData&  ebdataCoarA =              m_data[dita];
      EBData&  ebdataCoarB = a_otherPhase.m_data[ditb];
      // Box region = ebgrapCoarA.getRegion();
      Box region = m_grids[dita];

      // CH_assert(ebgrapCoarB.getRegion() == region);
      CH_assert(a_otherPhase.m_grids[ditb] == region);
      IntVectSet seta = ebgrapCoarA.getIrregCells(region);
      IntVectSet setb = ebgrapCoarB.getIrregCells(region);

      //different EBIS's can (correctly) disagree about which cells are multivalued
      IntVectSet setMultiA = ebgrapCoarA.getMultiCells(region);
      IntVectSet setMultiB = ebgrapCoarB.getMultiCells(region);
      IntVectSet setMulti = setMultiA | setMultiB;
      int aphase =              m_phase;
      int bphase = a_otherPhase.m_phase;

      IntVectSet setANoMulti = seta - setMulti;
      IntVectSet setBNoMulti = setb - setMulti;
      {
        //EBIndexSpace does the right thing with all the geometric information
        //in the case of two single valued vofs.   All that needs to be set is
        //the phase on the other side of the irregular face.
        for (IVSIterator ita(setANoMulti); ita.ok(); ++ita)
          {
            VolIndex v(ita(), 0); //0 because we know this whole set is single valued

            VolData& volDatA = ebdataCoarA.getVolData()(v,0);
            volDatA.m_phaseFaces.resize(1);
            //this sets boundary area and normal and centroid to
            //whatever the ebindexspace set it to.
            volDatA.m_phaseFaces[0]=volDatA.m_averageFace;
            bool skip = false;
            bool vContainedInOtherIVS = setBNoMulti.contains(v.gridIndex());
            setOtherVoFSingleValued(volDatA, v, bphase, ebdataCoarB, vContainedInOtherIVS, skip);
            if (skip && setMultiB.contains(volDatA.m_phaseFaces[0].m_volIndex.gridIndex()))
              {
                MayDay::Error("coordinate face boundary also a multi valued one");
              }
          }//end loop over cells  where both phases are singlevalued

        for (IVSIterator itb(setBNoMulti); itb.ok(); ++itb)
          {
            VolIndex v(itb(), 0); //0 because we know this whole set is single valued

            VolData& volDatB = ebdataCoarB.getVolData()(v,0);
            volDatB.m_phaseFaces.resize(1);
            //this sets boundary area and normal and centroid to
            //whatever the ebindexspace set it to.
            volDatB.m_phaseFaces[0]=volDatB.m_averageFace;
            bool skip = false;
            bool vContainedInOtherIVS = setANoMulti.contains(v.gridIndex());
            setOtherVoFSingleValued(volDatB, v, aphase, ebdataCoarA, vContainedInOtherIVS, skip);
            if (skip && setMultiA.contains(volDatB.m_phaseFaces[0].m_volIndex.gridIndex()))
              {
                MayDay::Error("coordinate face boundary also a multi valued one");
              }
          }//end loop over cells where both phases are singlevalued

        // now fix incorrect boundary data, where possible
        // first, for the phase a volData
//         for (IVSIterator ita(setANoMulti); ita.ok(); ++ita)
        // hack to get around mismatching sets
        IntVectSet abIntersectionNoMulti = setANoMulti & setBNoMulti;
        for (IVSIterator ita(abIntersectionNoMulti); ita.ok(); ++ita)
          {
            VolIndex v(ita(), 0);

            VolData& volDatA = ebdataCoarA.getVolData()(v,0);
            VolIndex vother = volDatA.m_phaseFaces[0].m_volIndex;
            // need to fix the normal
            VolData& volDatB = ebdataCoarB.getVolData()(vother,0);
            // only send the first boundary data, since single valued cell
            testAndFixBoundaryData(volDatA, volDatB);
          }
        // next, for the phase b volData
//         for (IVSIterator itb(setBNoMulti); itb.ok(); ++itb)
        // hack to get around mismatching sets
        for (IVSIterator itb(abIntersectionNoMulti); itb.ok(); ++itb)
          {
            VolIndex v(itb(), 0); //0 because we know this whole set is single valued

            VolData& volDatB = ebdataCoarB.getVolData()(v,0);
            VolIndex vother = volDatB.m_phaseFaces[0].m_volIndex;
            // need to fix the normal
            VolData& volDatA = ebdataCoarA.getVolData()(vother,0);
            // only send the first boundary data, since single valued cell
            testAndFixBoundaryData(volDatB, volDatA);
          }
      }
      if (!setMulti.isEmpty())
        {
          //now have to do the union of the each set
          //of multivalued cells.   either phase being multivalued can
          //make either ebindex space get the wrong answer for
          //the geometric information
          const EBISBox& ebisbxFineA = ebislFineA[dita()];
          const EBISBox& ebisbxFineB = ebislFineB[ditb()];

          IVSIterator it(setMulti); //see above derivation.
          for (it.begin(); it.ok(); ++it)
            {
              cellStitch(ebdataCoarA, ebdataCoarB,
                         ebgrapCoarA, ebgrapCoarB,
                         ebisbxFineA, ebisbxFineB,
                         it(), aphase, bphase);
            }//end loop over multivalued cells
        }
    } //end loop over grids
}

void getFinerBoundaryData(Vector<Vector<BoundaryData> > & a_bndryDataFineA,
                          Vector<VolIndex>&               a_otherCellIndex,
                          const VolIndex&                 a_coarVoFA,
                          const EBGraph&                  a_ebgrapCoarA,
                          const EBGraph&                  a_ebgrapCoarB,
                          const EBISBox&                  a_ebisbxFineA)
{
  //divdes fine boundary datas a by which coarse vof b they are connected to.
  // to which they are connected?

  const IntVect ivCoar = a_coarVoFA.gridIndex();

  Vector<VolIndex>          coarVoFsB = a_ebgrapCoarB.getVoFs(ivCoar);
  Vector<Vector<VolIndex> > fineVoFsB(coarVoFsB.size());
  a_bndryDataFineA.resize(coarVoFsB.size());
  a_otherCellIndex.resize(coarVoFsB.size());
  for (int ib = 0; ib < coarVoFsB.size(); ib++)
    {
      fineVoFsB[ib] = a_ebgrapCoarB.refine(coarVoFsB[ib]);
      a_otherCellIndex[ib] = coarVoFsB[ib];
    }

  Vector<VolIndex> allFinerVoFsA = a_ebgrapCoarA.refine(a_coarVoFA);
  for (int ib = 0; ib < coarVoFsB.size(); ib++)
    {
      a_bndryDataFineA[ib].resize(0);
      for (int ifineB = 0; ifineB < fineVoFsB[ib].size(); ifineB++)
        {
          const VolIndex& vofFineB = fineVoFsB[ib][ifineB];
          for (int ifineA = 0; ifineA < allFinerVoFsA.size(); ifineA++)
            {
              const VolIndex& vofFineA = allFinerVoFsA[ifineA];
              if (vofFineA.gridIndex() == vofFineB.gridIndex())
                {
                  if (a_ebisbxFineA.isIrregular(vofFineA.gridIndex()))
                    {
                      const VolData& voldat = a_ebisbxFineA.getEBData().getVolData()(vofFineA, 0);
                      const Vector<BoundaryData>& boundaryDat = voldat.m_phaseFaces;
                      for (int iface = 0; iface < boundaryDat.size(); iface++)
                        {
                          if (boundaryDat[iface].m_volIndex.cellIndex() == vofFineB.cellIndex())
                            {
                              a_bndryDataFineA[ib].push_back(boundaryDat[iface]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void coarsenBoundaryInfo(EBData&                                a_ebdataCoar,
                         const EBISBox&                         a_ebisbxFine,
                         const Vector< Vector<BoundaryData > >& a_finerBoundaryData,
                         const VolIndex&                        a_vof,
                         const int&                             a_otherPhase,
                         const Vector<VolIndex>&                a_otherCellIndex)
{
  //this coarsening factor for boundary area fractions can be thought of this way
  //(1) take fine area and multiply by dxf^{SD-1}
  //(2) add up areas (now we have the real coarse area).
  //(3) divide by dxc^{SD-1} = dxf{SD-1}*2{SD-1}
  Real faceCoarsenFactor = D_TERM(1.0, * 2.0,* 2.0);

  //the point of this routine is to fix the boundary area and normal
  //stuff in this object
  VolData& coarData = a_ebdataCoar.getVolData()(a_vof, 0);
  coarData.m_phaseFaces.resize(a_finerBoundaryData.size());
  for (int ib = 0; ib < a_finerBoundaryData.size(); ib++)
    {
      //initialization is important  in both of these things.
      //since i am taking the average
      RealVect aveNormal  = RealVect::Zero;
      RealVect aveCentro  = RealVect::Zero;
      Real sumArea=0;
      Real aveArea=0;
      int numfaces=0;
      for (int jb = 0; jb < a_finerBoundaryData[ib].size(); jb++)
        {
          const BoundaryData& boundaryDat = a_finerBoundaryData[ib][jb];
          const Real    & fineArea = boundaryDat.m_bndryArea;
          const RealVect& fineNorm = boundaryDat.m_normal;
          const RealVect& fineCent = boundaryDat.m_bndryCentroid;
          sumArea += fineArea;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              aveNormal[idir] += fineArea*fineNorm[idir];
              aveCentro[idir] += fineArea*fineCent[idir];
            }
          numfaces++;
        }

      //we are taking an average-weighted average
      //of each normal and centroid components so we need to divide by the area
      if (sumArea > 1.0e-12)
        {
          aveNormal /= sumArea;
          aveCentro /= sumArea;
        }
      Real sumSq;
      PolyGeom::unifyVector(aveNormal, sumSq);
      //this takes into account the fact that boundary areas are normalized
      //by dx^(SD-1).  See note above.
      aveArea  = sumArea/faceCoarsenFactor;

      coarData.m_phaseFaces[ib].m_normal            = aveNormal;
      coarData.m_phaseFaces[ib].m_bndryCentroid     = aveCentro;
      coarData.m_phaseFaces[ib].m_bndryArea         = aveArea;
      coarData.m_phaseFaces[ib].m_bndryPhase        = a_otherPhase;
      coarData.m_phaseFaces[ib].m_volIndex          = a_otherCellIndex[ib];
    }

}

void EBISLevel::cellStitch(EBData&        a_ebdataCoarA,
                           EBData&        a_ebdataCoarB,
                           const EBGraph& a_ebgrapCoarA,
                           const EBGraph& a_ebgrapCoarB,
                           const EBISBox& a_ebisbxFineA,
                           const EBISBox& a_ebisbxFineB,
                           const IntVect& a_iv,
                           const int&     a_phaseA,
                           const int&     a_phaseB)
{
  //pout() << "entering cellStitch" << endl;
  Vector<VolIndex> coarVoFsA = a_ebgrapCoarA.getVoFs(a_iv);
  Vector<VolIndex> coarVoFsB = a_ebgrapCoarB.getVoFs(a_iv);
  for (int ivof = 0; ivof < coarVoFsA.size(); ivof++)
    {
      const VolIndex& vof = coarVoFsA[ivof];
      // each element in finerVoFs is a list of finer vofs whose
      // opposing number on the other side of the irregular boundary
      // are connected to each other.   This will correspond with
      // the Vector<VolIndex> of the coarse vofs in the *other* phase.
      Vector<Vector< BoundaryData> > finerBndryData;
      Vector<VolIndex>  cellIndexB;
      getFinerBoundaryData(finerBndryData, cellIndexB, vof, a_ebgrapCoarA, a_ebgrapCoarB, a_ebisbxFineA);


      //once the above list is made, we can correctly coarsen the boundary
      //information using only the information in this phase.
      coarsenBoundaryInfo(a_ebdataCoarA, a_ebisbxFineA, finerBndryData, vof, a_phaseB, cellIndexB);
    }

  for (int ivof = 0; ivof < coarVoFsB.size(); ivof++)
    {
      const VolIndex& vof = coarVoFsB[ivof];
      // each element in finerVoFs is a list of finer vofs whose
      // opposing number on the other side of the irregular boundary
      // are connected to each other.   This will correspond with
      // the Vector<VolIndex> of the coarse vofs in the *other* phase.
      Vector<Vector< BoundaryData> > finerBndryData;
      Vector<VolIndex>  cellIndexA;
      getFinerBoundaryData(finerBndryData, cellIndexA, vof, a_ebgrapCoarB, a_ebgrapCoarA, a_ebisbxFineB);


      //once the above list is made, we can correctly coarsen the boundary
      //information using only the information in this phase.
      coarsenBoundaryInfo(a_ebdataCoarB, a_ebisbxFineB, finerBndryData, vof, a_phaseA, cellIndexA);
    }
}

#ifdef CH_USE_HDF5
void EBISLevel::write(HDF5Handle& a_handle,
                      const int& a_levelNumber) const
{
  CH_TIME("EBISLevel::write_with_lev_number");
  HDF5HeaderData header;
  char levelcstr[256];
  sprintf(levelcstr, "%d", a_levelNumber);
  string levelstring(levelcstr);

  header.m_realvect["EBIS_origin"] = m_origin;

  string domstring = string("EBIS_domain_lev_") + levelstring;
  header.m_box[domstring] = m_domain.domainBox();
  
  string dxstring = string("EBIS_dx_lev_") + levelstring;
  header.m_real[dxstring]     = m_dx;

  //  pout() << "writing handle  for level " << a_levelNumber << endl;
  header.writeToFile(a_handle);

  //write the grids to the file
  string boxstring = string("EBIS_boxes_lev_") + levelstring;
  CH_XD::write(a_handle, m_grids, boxstring);

  string graphstring = string("EBIS_graph_lev_") + levelstring;
  string datastring = string("EBIS_data_lev_") + levelstring;

  //  pout() << "writing graph for level " << a_levelNumber << endl;
  int eekflag = CH_XD::write(a_handle, m_graph, graphstring);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }
  //  pout() << "writing data  for level " << a_levelNumber << endl;
  eekflag = CH_XD::write(a_handle, m_data ,  datastring);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
}

EBISLevel::EBISLevel(HDF5Handle& a_handle, 
                     const int& a_levelNumber)
{
  CH_TIME("EBISLevel::EBISLevel_hdf5_with_lev_number");

  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;
  char levelcstr[256];
  sprintf(levelcstr, "%d", a_levelNumber);
  string levelstring(levelcstr);

  HDF5HeaderData header;
  header.readFromFile(a_handle);
  m_origin = header.m_realvect["EBIS_origin"];

  string domstring = string("EBIS_domain_lev_") + levelstring;
  m_domain = header.m_box[domstring];

  string dxstring = string("EBIS_dx_lev_") + levelstring;
  m_dx =     header.m_real[dxstring]    ;

  m_tolerance = m_dx*1E-4;

  //read in the grids
  Vector<Box> boxes;
  string boxstring = string("EBIS_boxes_lev_") + levelstring;
  read(a_handle,boxes, boxstring);
  Vector<int> procAssign;
  LoadBalance(procAssign, boxes);
  m_grids.define(boxes, procAssign);//this should use m_domain for periodic...
  EBGraphFactory graphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, graphfact);


  string graphstring = string("EBIS_graph_lev_") + levelstring;
  int eekflag = read(a_handle, m_graph, graphstring, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }

  //need a ghosted layout so that the data can be defined properly
  LevelData<EBGraph> ghostGraph(m_grids, 1, IntVect::Unit, graphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, ghostGraph, interv);

  //now the data for the graph
  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit()].defineVoFData(ghostGraph[dit()], m_grids.get(dit()));
      m_data[dit()].defineFaceData(ghostGraph[dit()], m_grids.get(dit()));
    }
  //read the data  in from the file
  string  datastring  = string("EBIS_data_lev_") + levelstring;
  eekflag = read(a_handle, m_data ,  datastring, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
#if 0
  pout() << "EBISLevel::EBISLevel 5 - m_grids - m_dx: " << m_dx << endl;
  pout() << "--------" << endl;
  pout() << m_grids.boxArray().size() << endl;
  pout() << "--------" << endl;
  pout() << endl;
#endif
}
#endif

const ProblemDomain& EBISLevel::getDomain() const
{
  return m_domain;
}

const Real& EBISLevel::getDX() const
{
  return m_dx;
}

const RealVect& EBISLevel::getOrigin() const
{
  return m_origin;
}

const DisjointBoxLayout& EBISLevel::getGrids() const
{
  return m_grids;
}

// All boxes containing an irregular cell
DisjointBoxLayout EBISLevel::getIrregGrids(const ProblemDomain& a_domain) const
{
  DisjointBoxLayout irregGrids;
  DataIterator dit =   m_graph.dataIterator();
  Vector<Box> localBoxes;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular())
        {
          localBoxes.push_back(b);

        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, 0);

  broadcast(allBoxes, 0);

  Vector<Box> boxes;
  for (int i = 0; i < allBoxes.size(); i++)
    {
      boxes.append(allBoxes[i]);
    }
  Vector<int> procs;
  LoadBalance(procs, boxes);
  irregGrids = DisjointBoxLayout(boxes, procs, a_domain);

  return irregGrids;
}

// All boxes with only irregular or regular cells
DisjointBoxLayout EBISLevel::getFlowGrids(const ProblemDomain& a_domain) const
{
  DisjointBoxLayout flowGrids;
  DataIterator dit =   m_graph.dataIterator();
  Vector<Box> localBoxes;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular() || graph.isAllRegular())
        {
          localBoxes.push_back(b);

        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, 0);

  broadcast(allBoxes, 0);
  Vector<Box> boxes;
  for (int i = 0; i < allBoxes.size(); i++)
    {
      boxes.append(allBoxes[i]);
    }
  Vector<int> procs;
  LoadBalance(procs, boxes);

  flowGrids = DisjointBoxLayout(boxes, procs, a_domain);
  return flowGrids;
}

// All boxes with an irregular or covered cell
DisjointBoxLayout EBISLevel::getCoveredGrids(const ProblemDomain& a_domain) const
{
  DisjointBoxLayout coveredGrids;
  DataIterator dit = m_graph.dataIterator();
  Vector<Box> localBoxes;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular() || graph.isAllCovered())
        {
          localBoxes.push_back(b);
        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, 0);

  broadcast(allBoxes, 0);

  Vector<Box> boxes;
  for (int i = 0; i < allBoxes.size(); i++)
    {
      boxes.append(allBoxes[i]);
    }
  Vector<int> procs;
  LoadBalance(procs, boxes);
  coveredGrids = DisjointBoxLayout(boxes, procs, a_domain);

  return coveredGrids;
}

IntVectSet EBISLevel::irregCells() const
{
  DataIterator dit =   m_graph.dataIterator();
  IntVectSet rtn;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular())
        {
          rtn |= graph.getIrregCells(b);
        }
    }
  return rtn;
}

void EBISLevel::getGraphSummary(long long & a_irrVoFs,
                                long long & a_arcs,
                                long long & a_multiVoFs,
                                long long & a_zeroVoFs,
                                long long & a_zeroVoFsArcs)
{
  CH_TIME("EBISLevel::getGraphSummary");

  a_irrVoFs      = 0;
  a_arcs         = 0;
  a_multiVoFs    = 0;
  a_zeroVoFs     = 0;
  a_zeroVoFsArcs = 0;

  for (DataIterator dit(m_grids); dit.ok(); ++dit)
  {
    const Box & curBox = m_grids[dit()];

    EBGraph & curGraph = m_graph[dit()];
    const EBData  & curData  = m_data[dit()];

    if (curGraph.m_implem->m_tag == EBGraphImplem::HasIrregular)
    {

      BaseFab<GraphNode> & curGraphFAB  = curGraph.m_implem->m_graph;
      const BaseIVFAB<VolData> & curDataIVFAB = curData.m_implem->m_volData;

      for (BoxIterator bit(curBox); bit.ok(); ++bit)
      {
        const IntVect   & curIV = bit();
        GraphNode & curGraphNode = curGraphFAB(curIV,0);

        if (curGraphNode.isIrregular())
        {
          int vofs = curGraphNode.m_cellList->size();

          a_irrVoFs += vofs;

          if (vofs > 1)
          {
            a_multiVoFs += vofs;
          }

          for (int vof=0; vof < vofs; vof++)
          {
            VolIndex volIndex(curIV,vof);

            Real volFrac = curDataIVFAB(volIndex,0).m_volFrac;

            if (volFrac == 0.0)
            {
              a_zeroVoFs++;
            }

            const GraphNodeImplem & curGraphNodeImplem = (*curGraphNode.m_cellList)[vof];
            int localVoFArcs = 0;

            for (int iside = 0; iside < 2*SpaceDim; iside++)
            {
              int arcs = curGraphNodeImplem.m_arc[iside].size();
              localVoFArcs += arcs;

            }

            if (volFrac == 0.0)
            {
              a_zeroVoFsArcs += localVoFArcs;

              if (vofs == 1 && localVoFArcs == 0)
              {
                delete curGraphNode.m_cellList;
                curGraphNode.m_cellList = 0;
                a_irrVoFs--;
              }
              else
              {
                a_arcs += localVoFArcs;
              }
            }
            else
            {
              a_arcs += localVoFArcs;
            }
          }
        }
      }
    }
  }

  a_irrVoFs      = EBLevelDataOps::parallelSum(a_irrVoFs);
  a_arcs         = EBLevelDataOps::parallelSum(a_arcs);
  a_multiVoFs    = EBLevelDataOps::parallelSum(a_multiVoFs);
  a_zeroVoFs     = EBLevelDataOps::parallelSum(a_zeroVoFs);
  a_zeroVoFsArcs = EBLevelDataOps::parallelSum(a_zeroVoFsArcs);
}

void EBISLevel::printGraphSummary(char const * a_prefix)
{
  CH_TIME("EBISLevel::printGraphSummary");

  long long irrVoFs;
  long long arcs;
  long long multiVoFs;
  long long zeroVoFs;
  long long zeroVoFsArcs;

  getGraphSummary(irrVoFs,arcs,multiVoFs,zeroVoFs,zeroVoFsArcs);

  pout() << a_prefix << "irreg VoFs    : " << irrVoFs      << endl;
  pout() << a_prefix << "arcs          : " << arcs         << endl;
  pout() << a_prefix << "multi VoFs    : " << multiVoFs    << endl;
  pout() << a_prefix << "zero VoFs     : " << zeroVoFs     << endl;
  pout() << a_prefix << "zero VoFs arcs: " << zeroVoFsArcs << endl;
}

#include "NamespaceFooter.H"
