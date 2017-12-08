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
    UtilCreateDirectoryDestructive(a_dirname, true);
    writeHeader(a_dirname);
    string graphdirname = a_dirname + "/_graph";
    string  datadirname = a_dirname + "/_data";
    UtilCreateDirectoryDestructive(graphdirname, true);
    UtilCreateDirectoryDestructive( datadirname, true);
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
    BL_PROFILE("EBISLevel::EBISLevel_geoserver_domain");

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

#if AMREX_SPACEDIM==2
  static
  std::list<Segment>::iterator FindMySeg(std::list<Segment>& segs, const NodeMapIt& idx)
  {
    for (std::list<Segment>::iterator it=segs.begin(); it!=segs.end(); ++it)
    {
      AMREX_ASSERT(it->size() == 2);
      if ( ((*it)[0] == idx) || ((*it)[1] == idx) )
	return it;
    }
    return segs.end();
  }
#endif
  
  void
  EBISLevel::buildEBSurface(const GeometryService & a_geoserver)
  {

    Vector<Vector<Triangle>> triangles;
    Vector<list<list<Segment>>> contourLines;

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
          segments.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                     intersections.find(Edge(pt[0],pt[3]))));
          break;
        case 2:
        case 13:
          segments.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                     intersections.find(Edge(pt[1],pt[2]))));
          break;
        case 3:
        case 12:
          segments.push_back(Segment(intersections.find(Edge(pt[1],pt[2])),
                                     intersections.find(Edge(pt[0],pt[3]))));
          break;
        case 4:
        case 11:
          segments.push_back(Segment(intersections.find(Edge(pt[1],pt[2])),
                                     intersections.find(Edge(pt[2],pt[3]))));
          break;
        case 6:
        case 9:
          segments.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                     intersections.find(Edge(pt[2],pt[3]))));
          break;
        case 7:
        case 8:
          segments.push_back(Segment(intersections.find(Edge(pt[2],pt[3])),
                                     intersections.find(Edge(pt[0],pt[3]))));
          break;
        case 5:
        case 10:
          segments.push_back(Segment(intersections.find(Edge(pt[0],pt[1])),
                                     intersections.find(Edge(pt[1],pt[2]))));
          segments.push_back(Segment(intersections.find(Edge(pt[2],pt[3])),
                                     intersections.find(Edge(pt[0],pt[3]))));
          break;
        }
      }

      contourLines.push_back(list<list<Segment>>());
      auto& cLines = contourLines.back();
      if (segments.size() > 0)
      {
        cLines.push_back(list<Segment>());
        cLines.back().push_back(segments.front());
        segments.pop_front();

        auto idx = cLines.back().back()[1];
        while (segments.begin() != segments.end())
        {
          auto segIt = FindMySeg(segments,idx);
          if (segIt != segments.end())
          {
            const auto& idx_l = (*segIt)[0];
            const auto& idx_r = (*segIt)[1];

            if ( idx_l == idx )
            {
              idx = idx_r;
              cLines.back().push_back(*segIt);
            }
            else
            {
              idx = idx_l;
              cLines.back().push_back(Segment(idx_r,idx_l));
            }

            segments.erase(segIt);
          }
          else
          {
            cLines.push_back(list<Segment>());
            cLines.back().push_back(segments.front());
            segments.pop_front();

            idx = cLines.back().back()[1];
          }
        }

        // Connect up the line segments as much as possible
        bool changed;
        do
        {
          changed = false;
          for (std::list<list<Segment>>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
          {
            if (!it->empty())
            {
              const auto& idx_l = it->front()[0];
              const auto& idx_r = it->back()[1];
              for (std::list<list<Segment>>::iterator it1 = cLines.begin(); it1!=cLines.end(); ++it1)
              {
                if (!it1->empty() && it!=it1)
                {
                  if (idx_r == it1->front()[0])
                  {
                    it->splice(it->end(),*it1);
                    changed = true;
                  }
                  else if (idx_r == it1->back()[1])
                  {
                    it1->reverse();
                    for (list<Segment>::iterator it2=it1->begin(); it2!=it1->end(); ++it2)
                    {
                      const auto tmp = (*it2)[0];
                      (*it2)[0] = (*it2)[1];
                      (*it2)[1] = tmp;
                    }
                    it->splice(it->end(),*it1);
                    changed = true;
                  }
                  else if (idx_l == it1->front()[0])
                  {
                    it1->reverse();
                    for (list<Segment>::iterator it2=it1->begin(); it2!=it1->end(); ++it2)
                    {
                      const auto tmp = (*it2)[0];
                      (*it2)[0] = (*it2)[1];
                      (*it2)[1] = tmp;
                    }
                    it->splice(it->begin(),*it1);
                    changed = true;
                  }
                }
              }
            }
          }
        } while(changed);

        // Clear out empty placeholders for lines we connected up to others.
        for (std::list<list<Segment>>::iterator it = cLines.begin(); it!=cLines.end();)
        {
          if (it->empty())
            cLines.erase(it++);
          else
            it++;
        }
      }

#else
      static int edgeTable[256]={
        0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
        0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
        0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
        0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
        0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
        0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
        0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
        0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
        0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
        0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
        0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
        0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
        0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
        0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
        0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
        0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
        0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
        0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
        0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
        0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
        0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
        0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
        0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
        0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
        0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
        0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
        0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
        0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
        0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

      static int triTable[256][16] =
        {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
         {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
         {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
         {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
         {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
         {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
         {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
         {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
         {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
         {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
         {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
         {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
         {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
         {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
         {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
         {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
         {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
         {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
         {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
         {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
         {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
         {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
         {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
         {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
         {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
         {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
         {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
         {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
         {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
         {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
         {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
         {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
         {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
         {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
         {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
         {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
         {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
         {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
         {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
         {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
         {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
         {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
         {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
         {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
         {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
         {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
         {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
         {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
         {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
         {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
         {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
         {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
         {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
         {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
         {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
         {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
         {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
         {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
         {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
         {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
         {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
         {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
         {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
         {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
         {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
         {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
         {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
         {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
         {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
         {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
         {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
         {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
         {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
         {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
         {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
         {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
         {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
         {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
         {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
         {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
         {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
         {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
         {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
         {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
         {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
         {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
         {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
         {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
         {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
         {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
         {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
         {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
         {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
         {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
         {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
         {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
         {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
         {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
         {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
         {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
         {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
         {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
         {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
         {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
         {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
         {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
         {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
         {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
         {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
         {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
         {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
         {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
         {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
         {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
         {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
         {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
         {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
         {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
         {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
         {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
         {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
         {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
         {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
         {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
         {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
         {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
         {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
         {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
         {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
         {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
         {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
         {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
         {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
         {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
         {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
         {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
         {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
         {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
         {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
         {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
         {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
         {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
         {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
         {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
         {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
         {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
         {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
         {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
         {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
         {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
         {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
         {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
         {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
         {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
         {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
         {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

      std::vector<IntVect> pt(8);
      std::vector<RealVect> x(8);
      std::vector<bool> v(8);

      triangles.push_back(Vector<Triangle>());

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
        if (edgeTable[cubeindex] == 0) break;

        std::array<NodeMapIt,12> vertlist;

        if (edgeTable[cubeindex] & 1)
          vertlist[0]  = intersections.find(Edge(pt[0],pt[1]));
        if (edgeTable[cubeindex] & 2)
          vertlist[1]  = intersections.find(Edge(pt[1],pt[2]));
        if (edgeTable[cubeindex] & 4)
          vertlist[2]  = intersections.find(Edge(pt[2],pt[3]));
        if (edgeTable[cubeindex] & 8)
          vertlist[3]  = intersections.find(Edge(pt[3],pt[0]));
        if (edgeTable[cubeindex] & 16)
          vertlist[4]  = intersections.find(Edge(pt[4],pt[5]));
        if (edgeTable[cubeindex] & 32)
          vertlist[5]  = intersections.find(Edge(pt[5],pt[6]));
        if (edgeTable[cubeindex] & 64)
          vertlist[6]  = intersections.find(Edge(pt[6],pt[7]));
        if (edgeTable[cubeindex] & 128)
          vertlist[7]  = intersections.find(Edge(pt[7],pt[4]));
        if (edgeTable[cubeindex] & 256)
          vertlist[8]  = intersections.find(Edge(pt[0],pt[4]));
        if (edgeTable[cubeindex] & 512)
          vertlist[9]  = intersections.find(Edge(pt[1],pt[5]));
        if (edgeTable[cubeindex] & 1024)
          vertlist[10]  = intersections.find(Edge(pt[2],pt[6]));
        if (edgeTable[cubeindex] & 2048)
          vertlist[11]  = intersections.find(Edge(pt[3],pt[7]));

        /* Create the triangles */
        int nTriang = 0;
        for (int i=0;triTable[cubeindex][i]!=-1;i+=3)
          nTriang++;

        for (int j=0; j<nTriang; ++j)
        {
          auto& triangleVec = triangles.back();
          int j3 = 3*j;
          triangleVec.push_back(Triangle(vertlist[triTable[cubeindex][j3  ]],
                                         vertlist[triTable[cubeindex][j3+1]],
                                         vertlist[triTable[cubeindex][j3+2]]));
        }
      }

#endif
    }

#if AMREX_SPACEDIM==2
    for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
    {
      auto& intersections = m_intersections[mfi];
      auto& cLines = contourLines[mfi.LocalIndex()];

      if (cLines.size() > 0)
      {
        std::string FullDataPath = eb_surface_filename;
        if (!FullDataPath.empty() && FullDataPath[FullDataPath.size()-1] != '/')
          FullDataPath += '/';
        FullDataPath += "Data";

        auto nGrid = m_intersections.size();
        auto nDigits = std::log10(nGrid) + 1;
        FullDataPath += Concatenate("_",mfi.index(),nDigits);

        std::ofstream ofs(FullDataPath.c_str());
        ofs << cLines.size() << '\n';
        for (std::list<list<Segment>>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
        {
          ofs << it->size() + 1 << '\n';
          auto slit=it->begin();
          const auto& pt0 = (*slit)[0]->second;
          ofs << pt0[0] << " " << pt0[1] << '\n';
          for (  ; slit!=it->end(); ++slit)
          {
            const auto& pt = (*slit)[1]->second;
            ofs << pt[0] << " " << pt[1] << '\n';
          }
        }
        ofs.close();
      }
    }
#else
    for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
    {
      auto& intersections = m_intersections[mfi];
      auto& fragVec = triangles[mfi.LocalIndex()];
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

    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(eb_surface_filename, 0755))
            amrex::CreateDirectoryFailed(eb_surface_filename);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    amrex::Print() << "EBISLevel: Writing EB surface to: " << eb_surface_filename << '\n';
    int nGrids = m_intersections.size();
    int ioProc = ParallelDescriptor::IOProcessorNumber();

    Vector<long> nNodes(nGrids,0);
    Vector<long> nElts(nGrids,0);

    for (MFIter mfi(m_intersections); mfi.isValid(); ++mfi)
    {
      nNodes[mfi.index()] = m_intersections[mfi].size();
      nElts[mfi.index()] = triangles[mfi.LocalIndex()].size();
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
#endif
  }

  ///
  void
  EBISLevel::defineFromGeometryService(const GeometryService & a_geoserver)
  {
    
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

