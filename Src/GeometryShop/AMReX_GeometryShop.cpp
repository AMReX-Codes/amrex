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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <list>

#include "AMReX_GeometryShop.H"
#include "AMReX_RealVect.H"
#include "AMReX.H"
#include "AMReX_Print.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_BoxIterator.H"


namespace amrex
{

//  static const IntVect   gs_debiv(D_DECL(994,213,7));
//  static const IntVect   gs_debivlo(D_DECL(190,15,0));
//  static const IntVect   gs_debivhi(D_DECL(191,15,0));
  
  bool GeometryShop::isRegularEveryPoint(const Box&           a_region,
                                         const Box&           a_domain,
                                         const RealVect&      a_origin,
                                         const Real&          a_dx) const
  {
//    BL_PROFILE("GeometryShop::isRegularEveryPoint");

    // All corner indices for the current box
    Box allCorners(a_region);
    allCorners.surroundingNodes();

    RealVect physCorner;
    BoxIterator bit(allCorners);
    // If every corner is inside, the box is regular

    for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
      {
        // Current corner
        const IntVect& corner = bit();

        // Compute physical coordinate of corner


        for (int idir = 0; idir < SpaceDim; ++idir)
        {
          physCorner[idir] = a_dx*corner[idir] + a_origin[idir];
        }

        // If the implicit function value is positive then the current
        // corner is outside the domain
        Real functionValue = m_implicitFunction->value(physCorner);

        if (functionValue > 0.0 )
        {
          return false;
        }
      }
      bit.reset();
      ++bit;
    }


    return true;
  }

  bool GeometryShop::isCoveredEveryPoint(const Box&           a_region,
                                         const Box&           a_domain,
                                         const RealVect&      a_origin,
                                         const Real&          a_dx) const
  {
//    BL_PROFILE("GeometryShop::isCoveredEveryPoint");

    // All corner indices for the current box
    Box allCorners(a_region);
    allCorners.surroundingNodes();

    RealVect physCorner;
    BoxIterator bit(allCorners);

    // If every corner is inside, the box is regular
    for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
      {
        // Current corner
        IntVect corner = bit();

        // Compute physical coordinate of corner

        for (int idir = 0; idir < SpaceDim; ++idir)
        {
          physCorner[idir] = a_dx*corner[idir] + a_origin[idir];
        }

        // If the implicit function value is positive then the current
        // corner is outside the domain
        Real functionValue = m_implicitFunction->value(physCorner);

        if (functionValue < 0.0 )
        {


          return false;
        }
      }
      bit.reset();
      ++bit;
    }

    return true;
  }
  GeometryShop::GeometryShop(const BaseIF& a_localGeom,
                             int           a_verbosity,
                             Real          a_thrshdVoF)
  {
    m_implicitFunction.reset(a_localGeom.newImplicitFunction());

    m_verbosity = a_verbosity;

    Real arg1 = 10.0;
    Real arg2 = -m_verbosity;

    m_threshold = 1.0e-15*pow(arg1, arg2);

    m_thrshdVoF = a_thrshdVoF;

  }

  /**********************************************/
  GeometryService::InOut 
  GeometryShop::
  InsideOutside(const Box&           a_region,
                const Box&           a_domain,
                const RealVect&      a_origin,
                const Real&          a_dx) const
  {
    GeometryShop::InOut rtn;
//begin debug
//    bool debugc = (a_region.contains(gs_debivlo) || a_region.contains(gs_debivhi));
//end debug
    if(isRegularEveryPoint(a_region, a_domain, a_origin, a_dx))
    {
      rtn = GeometryShop::Regular;
//begin debug
//      if(debugc)
//      {
//        amrex::AllPrint() << "geometryshop::insideoutside:"<< gs_debiv << " in an all regular box" << endl;
//      }
//end debug
    }
    else if(isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx))
    {
      rtn = GeometryShop::Covered;
//begin debug
//      if(debugc)
//      {
//        amrex::AllPrint() << "geometryshop::insideoutside:"<< gs_debiv << " in an all covered box" << endl;
//      }
//end debug
    }
    else
    {
      rtn = GeometryShop::Irregular;
    }
    return rtn;

  }



  struct Edge
  {
    Edge(const IntVect& lhs, const IntVect& rhs)
      {
	ID = -1;
	if (lhs < rhs) {
	  IV_l = lhs;
	  IV_r = rhs;
	} else {
	  IV_l = rhs;
	  IV_r = lhs;
	}
      }

    bool operator== (const Edge& rhs) const
      {
	return IV_l==rhs.IV_l && IV_r==rhs.IV_r;
      }
    bool operator< (const Edge& rhs) const
      {
        if (IV_l==rhs.IV_l) {
	  return IV_r < rhs.IV_r;
        }
	return IV_l < rhs.IV_l;
      }

    IntVect IV_l,IV_r;
    size_t ID;
  };

  std::ostream& operator<<(std::ostream& os, const Edge& e)
  {
    os << "EDGE: " << e.IV_l << " : " << e.IV_r << " ID=" << e.ID;
    return os;
  }

  typedef std::map<Edge, RealVect> PMap;

  typedef PMap::iterator PMapIt;

  struct PMapItCompare
  {
    bool operator() (const PMapIt& lhs, const PMapIt& rhs) const {
      const RealVect& l = lhs->second;
      const RealVect& r = rhs->second;
      const void* vl=&l;
      const void* vr=&r;
      return vl<vr;
    }
  };

  struct Triangle
  {
    Triangle() : p(3), mArea(-1) {}
    Real Area();
    const PMapIt& operator[] (int n) const { return p[n]; }
    PMapIt& operator[] (int n) { return p[n]; }
    static int xComp,yComp,zComp;
    Array<PMapIt> p;
    
  private:
    void my_area();
    Real mArea;
  };

  int Triangle::xComp = 0;
  int Triangle::yComp = 1;
  int Triangle::zComp = 2;

  Real
  Triangle::Area()
  {
    if (mArea<0)
      my_area();
    return mArea;
  }

  void
  Triangle::my_area()
  {
#ifndef NDEBUG
    for (int i=0; i<p.size(); ++i)
    {
      BL_ASSERT(xComp>=0 && xComp<(*p[i]).second.size());
      BL_ASSERT(yComp>=0 && yComp<(*p[i]).second.size());
      BL_ASSERT(zComp>=0 && zComp<(*p[i]).second.size());
    }
#endif

    const RealVect& p0 = (*p[0]).second;
    const RealVect& p1 = (*p[1]).second;
    const RealVect& p2 = (*p[2]).second;
    mArea = 0.5*sqrt(
      pow(  ( p1[yComp] - p0[yComp])*(p2[zComp]-p0[zComp]) 
	    -(p1[zComp] - p0[zComp])*(p2[yComp]-p0[yComp]), 2)
      
      + pow(( p1[zComp] - p0[zComp])*(p2[xComp]-p0[xComp]) 
	    -(p1[xComp] - p0[xComp])*(p2[zComp]-p0[zComp]), 2)
      
      + pow(( p1[xComp] - p0[xComp])*(p2[yComp]-p0[yComp]) 
	    -(p1[yComp] - p0[yComp])*(p2[xComp]-p0[xComp]), 2));
  }    

  struct Segment
  {
    Segment() {}
    Segment(const PMapIt& p1, const PMapIt& p2) : p({p1, p2}) {}
    const PMapIt& ID_l() const {return p[0];}
    const PMapIt& ID_r() const {return p[1];}
    const PMapIt& operator[] (int n) const { return p[n]; }
    PMapIt& operator[] (int n) { return p[n]; }
    void flip ();
    std::array<PMapIt,2> p;
  };

  void
  Segment::flip()
  {
    PMapIt ptmp = p[0];
    p[0] = p[1];
    p[1] = ptmp;
  }

  list<Segment>::iterator FindMySeg(list<Segment>& segs, const PMapIt& idx)
  {
    for (list<Segment>::iterator it=segs.begin(); it!=segs.end(); ++it)
    {
      if ( ((*it).ID_l() == idx) || ((*it).ID_r() == idx) )
	return it;
    }
    return segs.end();
  }

  
  PMap vertCache;
  list<Segment> segments;

  pair<bool,PMapIt> SetIntersect(const IntVect& p1,const IntVect& p2,const RealVect& intersect,PMap& vertCache)
  {
    PMapIt fwd,rev;
    Edge edge(p1,p2);
    fwd = vertCache.find(edge);
    if (fwd == vertCache.end())
    {
      rev = vertCache.find(Edge(p2,p1));
        
      if (rev == vertCache.end())
      {
	edge.ID = vertCache.size(); // Number this edge
	std::pair<Edge,RealVect> ent(edge,intersect);
	std::pair<PMapIt,bool> it = vertCache.insert(ent);
	BL_ASSERT(it.second);
	return make_pair(true,it.first);
      }
      else
      {
	return make_pair(false,rev);
      }
    }
    return make_pair(false,fwd);
  }

  /*********************************************/
  void
  GeometryShop::fillGraph(BaseFab<int>        & a_regIrregCovered,
                          Vector<IrregNode>   & a_nodes,
                          const Box           & a_validRegion,
                          const Box           & a_ghostRegion,
                          const Box           & a_domain,
                          const RealVect      & a_origin,
                          const Real          & a_dx) const
  {
    AMREX_ASSERT(a_domain.contains(a_ghostRegion));
    a_nodes.resize(0);
    a_regIrregCovered.resize(a_ghostRegion, 1);


    Real thrshd = m_thrshdVoF;

    IntVectSet ivsirreg; // wz. no need to be a set
    IntVectSet ivsdrop ;
    long int numCovered=0, numReg=0, numIrreg=0;


    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
      {
        const IntVect iv =bit();

        Box miniBox(iv, iv);
        GeometryShop::InOut inout = InsideOutside(miniBox, a_domain, a_origin, a_dx);

        if (inout == GeometryShop::Covered)
          {
            // set covered cells to -1
            a_regIrregCovered(iv, 0) = -1;
            numCovered++;
          }
        else if (inout == GeometryShop::Regular)
          {
            // set regular cells to 1
            a_regIrregCovered(iv, 0) =  1;
            numReg++;
          }
        else
          {
            // set irregular cells to 0
            a_regIrregCovered(iv, 0) =  0;
            if (a_validRegion.contains(iv))
              {
                ivsirreg |= iv;
                numIrreg++;
              }
          }
      }

    //if a regular is next to a  covered, change to irregular with correct arcs and so on.
    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
      {
        const IntVect iv =bit();

      if(a_regIrregCovered(iv, 0) == -1)
        {
          fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);
        }
      }

    // now loop through irregular cells and make nodes for each  one.
    for(IVSIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
      {
        const IntVect& iv = ivsit();

        Real     volFrac, bndryArea;
        RealVect normal, volCentroid, bndryCentroid;
        std::array<Vector<int>,SpaceDim> loArc;
        std::array<Vector<int>,SpaceDim> hiArc;
        std::array<Vector<Real>,SpaceDim> loAreaFrac;
        std::array<Vector<Real>,SpaceDim> hiAreaFrac;
        std::array<Vector<RealVect>,SpaceDim> loFaceCentroid;
        std::array<Vector<RealVect>,SpaceDim> hiFaceCentroid;
        computeVoFInternals(volFrac,
                            loArc,
                            hiArc,
                            loAreaFrac,
                            hiAreaFrac,
                            bndryArea,
                            normal,
                            volCentroid,
                            bndryCentroid,
                            loFaceCentroid,
                            hiFaceCentroid,
                            a_regIrregCovered,
                            a_domain,
                            a_origin,
                            a_dx,
                            iv);


        if (thrshd > 0. && volFrac < thrshd)
          {
            ivsdrop |= iv;
            a_regIrregCovered(iv, 0) = -1;
          }//CP record these nodes to be removed
        else
          {
            IrregNode newNode;
            newNode.m_cell          = iv;
            newNode.m_volFrac       = volFrac;
            newNode.m_cellIndex     = 0;
            newNode.m_volCentroid   = volCentroid;
            newNode.m_bndryCentroid = bndryCentroid;

            for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
              {
                int loNodeInd = newNode.index(faceDir, Side::Lo);
                int hiNodeInd = newNode.index(faceDir, Side::Hi);
                newNode.m_arc[loNodeInd]          = loArc[faceDir];
                newNode.m_arc[hiNodeInd]          = hiArc[faceDir];
                newNode.m_areaFrac[loNodeInd]     = loAreaFrac[faceDir];
                newNode.m_areaFrac[hiNodeInd]     = hiAreaFrac[faceDir];
                newNode.m_faceCentroid[loNodeInd] = loFaceCentroid[faceDir];
                newNode.m_faceCentroid[hiNodeInd] = hiFaceCentroid[faceDir];
              }
            a_nodes.push_back(newNode);
          }

      } // end loop over cut cells in the box

#if BL_SPACEDIM==2
    Vector<IntVect> pt(4);
    Vector<RealVect> x(4);
    Vector<Real> v(4);
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
	  x[i][idir] = a_origin[idir] + pt[i][idir]*a_dx;
	}
	v[i] = m_implicitFunction->value(x[i]);
      }

      int segCase = 0;
      if (v[0] < 0) segCase |= 1;
      if (v[1] < 0) segCase |= 2;
      if (v[2] < 0) segCase |= 4;
      if (v[3] < 0) segCase |= 8;
      
      switch (segCase)
      {
      case 1:
      case 14:
	segments.push_back(Segment(vertCache.find(Edge(pt[0],pt[1])),
				   vertCache.find(Edge(pt[0],pt[3]))));
	break;
      case 2:
      case 13:
	segments.push_back(Segment(vertCache.find(Edge(pt[0],pt[1])),
				   vertCache.find(Edge(pt[1],pt[2]))));
	break;
      case 3:
      case 12:
	segments.push_back(Segment(vertCache.find(Edge(pt[1],pt[2])),
				   vertCache.find(Edge(pt[0],pt[3]))));
	break;
      case 4:
      case 11:
	segments.push_back(Segment(vertCache.find(Edge(pt[1],pt[2])),
				   vertCache.find(Edge(pt[2],pt[3]))));
       break;
      case 6:
      case 9:
	segments.push_back(Segment(vertCache.find(Edge(pt[0],pt[1])),
				   vertCache.find(Edge(pt[2],pt[3]))));
	break;
      case 7:
      case 8:
	segments.push_back(Segment(vertCache.find(Edge(pt[2],pt[3])),
				   vertCache.find(Edge(pt[0],pt[3]))));
	break;
      case 5:
      case 10:
	segments.push_back(Segment(vertCache.find(Edge(pt[0],pt[1])),
				   vertCache.find(Edge(pt[1],pt[2]))));
	segments.push_back(Segment(vertCache.find(Edge(pt[2],pt[3])),
				   vertCache.find(Edge(pt[0],pt[3]))));
	break;
      }
    }

    list<list<Segment>> cLines;
    cLines.push_back(list<Segment>());
    cLines.back().push_back(segments.front());

    PMapIt idx = segments.front().ID_r();
    segments.pop_front();
    while (segments.begin() != segments.end())
    {
      list<Segment>::iterator segIt = FindMySeg(segments,idx);
      if (segIt != segments.end())
      {
	const PMapIt& idx_l = (*segIt).ID_l();
	const PMapIt& idx_r = (*segIt).ID_r();

	if ( idx_l == idx )
	{
	  idx = idx_r;
	  cLines.back().push_back(*segIt);
	}
	else
	{
	  idx = idx_l;
	  Segment nseg(*segIt); nseg.flip();
	  cLines.back().push_back(nseg);
	}

	segments.erase(segIt);
      }
      else
      {
	cLines.push_back(list<Segment>());
	idx = segments.front().ID_r();
	segments.pop_front();
      }
    }

    // Connect up the line fragments as much as possible
    bool changed;
    do
    {
      changed = false;
      for (std::list<list<Segment>>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
      {
	if (!it->empty())
	{
	  const PMapIt& idx_l = it->front().ID_l();
	  const PMapIt& idx_r = it->back().ID_r();
	  for (std::list<list<Segment>>::iterator it1 = cLines.begin(); it1!=cLines.end(); ++it1)
	  {
	    if (!it1->empty() && it!=it1)
	    {
	      if (idx_r == it1->front().ID_l())
	      {
		it->splice(it->end(),*it1);
		changed = true;
	      }
	      else if (idx_r == it1->back().ID_r())
	      {
		it1->reverse();
		for (list<Segment>::iterator it2=it1->begin(); it2!=it1->end(); ++it2)
		  it2->flip();
		it->splice(it->end(),*it1);
		changed = true;
	      }
	      else if (idx_l == it1->front().ID_l())
	      {
		it1->reverse();
		for (list<Segment>::iterator it2=it1->begin(); it2!=it1->end(); ++it2)
		  it2->flip();
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

    std::ofstream ofs("contourLines.dat");
    for (std::list<list<Segment>>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
    {
      list<Segment>::const_iterator slit=it->begin();
      const RealVect& pt0 = slit->ID_l()->second;
      ofs << pt0[0] << " " << pt0[1] << std::endl;
      for (  ; slit!=it->end(); ++slit)
      {
	const RealVect& pt = slit->ID_r()->second;
	ofs << pt[0] << " " << pt[1] << std::endl;
      }
    }
    ofs.close();

    amrex::Print() << "************** number of cLines: " << cLines.size() << std::endl;
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

    Vector<IntVect> pt(8);
    Vector<RealVect> x(8);
    Vector<Real> v(8);
    Vector<Triangle> triangles; // result

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
	  x[i][idir] = a_origin[idir] + pt[i][idir]*a_dx;
	}
	v[i] = m_implicitFunction->value(x[i]);
      }

      int cubeindex = 0;
      if (v[0] < 0) cubeindex |= 1;
      if (v[1] < 0) cubeindex |= 2;
      if (v[2] < 0) cubeindex |= 4;
      if (v[3] < 0) cubeindex |= 8;
      if (v[4] < 0) cubeindex |= 16;
      if (v[5] < 0) cubeindex |= 32;
      if (v[6] < 0) cubeindex |= 64;
      if (v[7] < 0) cubeindex |= 128;
      
      /* Cube is entirely in/out of the surface */
      if (edgeTable[cubeindex] == 0) break;

      std::array<PMapIt,12> vertlist;

      if (edgeTable[cubeindex] & 1)
        vertlist[0]  = vertCache.find(Edge(pt[0],pt[1]));
      if (edgeTable[cubeindex] & 2)
        vertlist[1]  = vertCache.find(Edge(pt[1],pt[2]));
      if (edgeTable[cubeindex] & 4)
        vertlist[2]  = vertCache.find(Edge(pt[2],pt[3]));
      if (edgeTable[cubeindex] & 8)
        vertlist[3]  = vertCache.find(Edge(pt[3],pt[0]));
      if (edgeTable[cubeindex] & 16)
        vertlist[4]  = vertCache.find(Edge(pt[4],pt[5]));
      if (edgeTable[cubeindex] & 32)
        vertlist[5]  = vertCache.find(Edge(pt[5],pt[6]));
      if (edgeTable[cubeindex] & 64)
        vertlist[6]  = vertCache.find(Edge(pt[6],pt[7]));
      if (edgeTable[cubeindex] & 128)
        vertlist[7]  = vertCache.find(Edge(pt[7],pt[4]));
      if (edgeTable[cubeindex] & 256)
        vertlist[8]  = vertCache.find(Edge(pt[0],pt[4]));
      if (edgeTable[cubeindex] & 512)
        vertlist[9]  = vertCache.find(Edge(pt[1],pt[5]));
      if (edgeTable[cubeindex] & 1024)
        vertlist[10]  = vertCache.find(Edge(pt[2],pt[6]));
      if (edgeTable[cubeindex] & 2048)
        vertlist[11]  = vertCache.find(Edge(pt[3],pt[7]));

      /* Create the triangles */
      int nTriang = 0;
      for (int i=0;triTable[cubeindex][i]!=-1;i+=3)
        nTriang++;

      for (int j=0; j<nTriang; ++j)
      {
        triangles.push_back(Triangle());
        Triangle& tri = triangles.back();
        int j3 = 3*j;
        tri[0] = vertlist[triTable[cubeindex][j3  ]];
        tri[1] = vertlist[triTable[cubeindex][j3+1]];
        tri[2] = vertlist[triTable[cubeindex][j3+2]];
      }
    }
    amrex::Print() << "************** number of triangles: " << triangles.size() << std::endl;

    Vector<PMapIt> orderedNodes(vertCache.size());
    for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it)
    {
      orderedNodes[it->first.ID] = it;
    }

    std::ofstream ofs("contourSurfaces.dat");
    ofs << orderedNodes.size() << " " << triangles.size() << std::endl;
    for (int i=0; i<orderedNodes.size(); ++i)
    {
      const RealVect& vec = orderedNodes[i]->second;
      ofs << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
    }

    for (int i=0; i<triangles.size(); ++i)
    {
      const Triangle& tri = triangles[i];
      for (int j=0; j<3; ++j)
      {
        ofs << tri[j]->first.ID + 1 << " ";
      }
      ofs << std::endl;
    }
    ofs.close();

#endif

    amrex::Print() << "************** number of intersections: " << vertCache.size() << std::endl;

    // CP: fix sweep that removes cells with volFrac less than a certain threshold
    for(IVSIterator ivsit(ivsdrop); ivsit.ok(); ++ivsit)
      {
        const IntVect& iv = ivsit();
  
        for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
              {
                int isign = sign(sit());
                IntVect otherIV = iv + isign*BASISV(faceDir);
                if (a_validRegion.contains(otherIV))
                  {
                    if (a_regIrregCovered(otherIV,0) == 0)
                      {
                        // i am in the case where the other cell
                        // is also irregular.   I just made a previously
                        // irregular cell covered so I have to check to
                        // see if it had any faces pointed this way.
                        int inode = -1;
                        bool found = false;
                        for (int ivec = 0; ivec < a_nodes.size() && ! found; ivec++)
                          {
                            if (a_nodes[ivec].m_cell == otherIV)
                              {
                                inode = ivec;
                                found = true;
                              }
                          }
                        if (!found && a_validRegion.contains(otherIV))
                          {
                            amrex::Abort("something wrong in our logic");
                          }
                        if (found)
                          {
                            int arcindex = a_nodes[inode].index(faceDir, flip(sit()));
                            a_nodes[inode].m_arc[         arcindex].resize(0);
                            a_nodes[inode].m_areaFrac[    arcindex].resize(0);
                            a_nodes[inode].m_faceCentroid[arcindex].resize(0);
                          }
                      }
                  }//valid region
              }//sit
          }//facedir
        
        //also need to fix regular cells next to new covered cell
        fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);

      }//ivsdrop
    if(m_verbosity > 2)
    {
      amrex::AllPrint() << "numIrreg  = " << numIrreg << "\n";
      amrex::AllPrint() << "number of nodes  = " << a_nodes.size() << "\n";
    }
  }
  /*************/
  void
  GeometryShop::
  fixRegularCellsNextToCovered(Vector<IrregNode>    & a_nodes, 
                               BaseFab<int>        & a_regIrregCovered,
                               const Box           & a_validRegion,
                               const Box           & a_domain,
                               const IntVect       & a_iv,
                               const Real          & a_dx) const

  {
    Box grownBox(a_iv, a_iv);
    grownBox.grow(1);
    grownBox  &= a_domain;
    IntVectSet ivstocheck(grownBox);
    ivstocheck -= a_iv;
    Box ghostRegion = a_regIrregCovered.box();
    //first check neighbors in each direction.  
    //If any of these are regular, they are replaced 
    //by irregular cells with a boundary face facing the covered cell.
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            int ishift = sign(sit());
            IntVect ivshift = a_iv + ishift*BASISV(idir);
            ivstocheck -= ivshift;
            int bfvalshift = -1;
            if(ghostRegion.contains(ivshift))
              {
                bfvalshift = a_regIrregCovered(ivshift, 0);
              }
            if(bfvalshift  == 1)
              {
                a_regIrregCovered(ivshift, 0) =  0;

                if(a_validRegion.contains(ivshift))
                  {
                    IrregNode newNode;
                    getFullNodeWithCoveredFace(newNode, 
                                               a_regIrregCovered,
                                               ivshift, 
                                               a_domain);
                    a_nodes.push_back(newNode);
                  }

              }
          }
      }
    //next we loop through the remaining cells (corner cells in 2d, corner and edge cells in 3D)
    //if any of these are regular, we change them to irregular 
    for(IVSIterator ivsit(ivstocheck); ivsit.ok(); ++ivsit)
      {
        const IntVect& iv = ivsit();
        if(ghostRegion.contains(iv))
          {
            if(a_regIrregCovered(iv, 0) == 1)
              {
                a_regIrregCovered(iv, 0) = 0;
                IrregNode newNode;
                newNode.makeRegular(iv, a_dx, a_domain);
                a_nodes.push_back(newNode);
              }
          }
      }
  }
  /**********************************************/
  void
  GeometryShop::
  getFullNodeWithCoveredFace(IrregNode            & a_newNode, 
                             const BaseFab<int>   & a_regIrregCovered,
                             const IntVect        & a_iv,
                             const Box            & a_domain) const
  {

    a_newNode.m_cell          = a_iv;
    a_newNode.m_volFrac       = 1.0;
    a_newNode.m_cellIndex     = 0;
    a_newNode.m_volCentroid   = RealVect::Zero;
    //set regular cell values then fix up
    a_newNode.m_bndryCentroid = RealVect::Zero;
    int coveredDir;
    Side::LoHiSide coveredSide;
    bool found = false;

    for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            int ishift = sign(sit());
            IntVect ivshift = a_iv + ishift*BASISV(faceDir);
            Vector<int> arc;
            Vector<Real> areaFrac;
            Vector<RealVect> faceCentroid;
            if(!a_domain.contains(ivshift))
              {
                // boundary arcs always -1
                arc.resize(1,-1);
                areaFrac.resize(1, 1.0);
                faceCentroid.resize(1, RealVect::Zero);
              }
            else if (a_regIrregCovered(ivshift, 0) >= 0)
              {
                //irregular cell or regular cell
                //compute vof internals returns something special if 
                //connected to a regular cell but EBGraph treats both the  same.
                //it just  knows that the cell index of a regular cell is 0
                arc.resize(1,0);
                areaFrac.resize(1, 1.0);
                faceCentroid.resize(1, RealVect::Zero);
              }
            else if (a_regIrregCovered(ivshift, 0) < 0)
              {
                found = true;
                coveredDir= faceDir;
                coveredSide = sit();
                // covered face!
                arc.resize(0);
                areaFrac.resize(0);
                faceCentroid.resize(0);
              }
            else
              {
                amrex::Error("logic error");
              }
          
            int nodeInd = a_newNode.index(faceDir, sit());
            a_newNode.m_arc[nodeInd]          = arc;
            a_newNode.m_areaFrac[nodeInd]     = areaFrac;
            a_newNode.m_faceCentroid[nodeInd] = faceCentroid;
          }
      }
    //fix boundary centroid
    if(found)
      {
        int centsign = sign(coveredSide);
        a_newNode.m_bndryCentroid[coveredDir] =  centsign*0.5;
      }
  }
  /*********************************************/
  void
  GeometryShop::computeVoFInternals(Real&               a_volFrac,
                                    std::array<Vector<int>,SpaceDim>& a_loArc,
                                    std::array<Vector<int>,SpaceDim>& a_hiArc,
                                    std::array<Vector<Real>,SpaceDim>& a_loAreaFrac,
                                    std::array<Vector<Real>,SpaceDim>& a_hiAreaFrac,
                                    Real&                a_bndryArea,
                                    RealVect&            a_normal,
                                    RealVect&            a_volCentroid,
                                    RealVect&            a_bndryCentroid,
                                    std::array<Vector<RealVect>,SpaceDim>& a_loFaceCentroid,
                                    std::array<Vector<RealVect>,SpaceDim>& a_hiFaceCentroid,
                                    const BaseFab<int>&  a_regIrregCovered,
                                    const Box&           a_domain,
                                    const RealVect&      a_origin,
                                    const Real&          a_dx,
                                    const IntVect&       a_iv) const
  {

    // need maxDx to properly scale a_bndryArea
    Real maxDx = a_dx;


    if (SpaceDim == 2)
      {
        // In 2D a vof is a faceMo = edgeMo[4],boundary length and normal vector
        faceMo Face;
        edgeMo edges[4];

        bool faceCovered;
        bool faceRegular;
        bool faceDontKnow;
        int faceNormal = 2;

        // get edgeType and intersection points
        edgeData2D(edges,
                   faceCovered,
                   faceRegular,
                   faceDontKnow,
                   a_dx,
                   a_iv,
                   a_domain,
                   a_origin);

        AMREX_ASSERT(faceRegular || faceCovered || faceDontKnow);
        AMREX_ASSERT((!(faceRegular && faceCovered)) && (!(faceRegular && faceDontKnow)) && (!(faceDontKnow && faceCovered)));
        // define the faceMo
        Face.define(edges,faceNormal,faceCovered,faceRegular,faceDontKnow);

        Moments geom;

        // answer0,answer1 are vectors whose components contain geometric information
        Vector<Real> answer0;
        Vector<Real> answer1;

        int order = 0;
        answer0 = geom.momentCalc2D(order,Face);

        // extract the info from answer0 and answer1

        // volfrac
        a_volFrac = answer0[0];

        order = 1;
        answer1 = geom.momentCalc2D(order,Face);

        // centroid
        for (int idir=0; idir<SpaceDim;++idir)
          {
            a_volCentroid[idir] = answer1[SpaceDim-1-idir];
          }

        if (a_volFrac <= 0.0)
          {
            a_volCentroid = RealVect::Zero;
          }
        else
          {
            a_volCentroid /= a_volFrac;
          }

        // normal
        Real normalVec[SpaceDim];
        Face.getNormal(normalVec);
        for (int idir = 0;idir < SpaceDim;++idir)
          {
            a_normal[idir] = normalVec[idir];
          }

        for (int idir = 0;idir < SpaceDim;++idir)
          {
            // (nx,ny)->(nxdy,nydx)
            a_normal[idir] = normalVec[idir]*a_dx;
          }
        Real anisBd = 0.0;
        for (int idir = 0;idir < SpaceDim;++idir)
          {
            anisBd += a_normal[idir]*a_normal[idir];
          }
        anisBd = sqrt(anisBd);
        if (anisBd !=0.0)
          {
            a_normal /= anisBd;
          }

        // compute bndryArea and bndryCentroid
        a_bndryArea = Face.getBdLength();
        if (a_bndryArea <= 0.0)
          {
            a_bndryCentroid = RealVect::Zero;
            a_bndryArea = 0.0;
          }
        else
          {
            for (int idir = 0;idir < SpaceDim; ++idir)
              {
                a_bndryCentroid[idir] = answer0[SpaceDim-idir]/a_bndryArea;
              }
            a_bndryArea *= anisBd;
            a_bndryArea /= maxDx;
          }

        for (int edgeNormal = 0; edgeNormal < SpaceDim; ++edgeNormal)
          {
            // loside
            bool coveredLo = edges[edgeNormal*2].isCovered();
            IntVect loiv = a_iv;
            loiv[edgeNormal] -= 1;
            if(a_domain.contains(loiv))
              {
                coveredLo= coveredLo || (a_regIrregCovered(loiv) < 0);
              }
            if (coveredLo)
              {
                a_loArc[edgeNormal].resize(0);
                a_loFaceCentroid[edgeNormal].resize(0);
                a_loAreaFrac[edgeNormal].resize(0);
              }
            else if (!coveredLo)
              {
                a_loArc[edgeNormal].resize(1);
                a_loFaceCentroid[edgeNormal].resize(1);
                a_loAreaFrac[edgeNormal].resize(1);
                IntVect otherIV = a_iv;
                otherIV[edgeNormal] -= 1;
                if (a_domain.contains(otherIV))
                  {
                    //stuff that uses this does not care whether the arcs go
                    //to a regular or irregular cell--I am taking away that distinction
                    a_loArc[edgeNormal][0]= 0;
                  }
                else
                  {
                    // boundary arcs always -1
                    a_loArc[edgeNormal][0] = -1;
                  }
                a_loFaceCentroid[edgeNormal][0] = edges[edgeNormal*2].getEdgeCentroid();
                a_loAreaFrac[edgeNormal][0] = edges[edgeNormal*2].getEdgeLength();
              }

            // hiside
            bool coveredHi = edges[edgeNormal*2+1].isCovered();
            IntVect hiiv = a_iv;
            hiiv[edgeNormal] += 1;
            if(a_domain.contains(hiiv))
              {
                coveredHi = coveredHi || (a_regIrregCovered(hiiv, 0) < 0);
              }
            if (coveredHi)
              {
                a_hiArc[edgeNormal].resize(0);
                a_hiFaceCentroid[edgeNormal].resize(0);
                a_hiAreaFrac[edgeNormal].resize(0);
              }
            else if (!coveredHi)
              {
                a_hiArc[edgeNormal].resize(1);
                a_hiFaceCentroid[edgeNormal].resize(1);
                a_hiAreaFrac[edgeNormal].resize(1);
                IntVect otherIV = a_iv;
                otherIV[edgeNormal] += 1;
                if (a_domain.contains(otherIV))
                  {
                    //stuff that uses this does not care whether the arcs go
                    //to a regular or irregular cell--I am taking away that distinction.
                    //It only ever made sense in the context of a geometry gneration package that could
                    //generate multi-valued cells
                    a_hiArc[edgeNormal][0] = 0;
                  }
                else
                  {
                    // boundary arcs always -1
                    a_hiArc[edgeNormal][0] = -1;
                  }
                a_hiFaceCentroid[edgeNormal][0] = edges[edgeNormal*2+1].getEdgeCentroid();
                a_hiAreaFrac[edgeNormal][0] = edges[edgeNormal*2+1].getEdgeLength();

              }

          }
      }

    if (SpaceDim==3)
      {
        // 1) using the intvect, build up the classes in Moments: edgeMO,
        //    faceMo and finally vofMo
        // 2) check for covered or regular faces
        // 3) call  momentCalc3D
        // 4) keep track of what the output means and fill the variables
        //    requested

        faceMo Faces[6];
        int index = -1;

        for (int faceNormal = 0;faceNormal < SpaceDim;++faceNormal)
          {
            for (int hiLoFace = 0;hiLoFace < 2;++hiLoFace)
              {
                index += 1;
                edgeMo edges[4];
                bool faceCovered;
                bool faceRegular;
                bool faceDontKnow;
                edgeData3D(edges,
                           faceCovered,
                           faceRegular,
                           faceDontKnow,
                           hiLoFace,
                           faceNormal,
                           a_dx,
                           a_iv,
                           a_domain,
                           a_origin);

                AMREX_ASSERT(faceRegular || faceCovered || faceDontKnow);
                AMREX_ASSERT((!(faceRegular && faceCovered)) && (!(faceRegular && faceDontKnow)) && (!(faceDontKnow && faceCovered)));
                Faces[index].define(edges,faceNormal,faceCovered,faceRegular,faceDontKnow);

                // if the face is covered we will deal with it later
                if (!Faces[index].isCovered())
                  {
                    Moments geom;
                    // answer0 and answer1 have all the geometric facts
                    Vector<Real> answer0;
                    Vector<Real> answer1;

                    int order = 0;
                    answer0 = geom.momentCalc2D(order,Faces[index]);

                    // area of this face
                    Real area=answer0[0];

                    order = 1;
                    answer1 = geom.momentCalc2D(order,Faces[index]);

                    // first we get the centroid of the face
                    RealVect faceCentroid;
                    for (int idir = 0; idir<SpaceDim;++idir)
                      {
                        faceCentroid[idir] = answer1[SpaceDim-1-idir];
                      }

                    if (area > 0.0)
                      {
                        faceCentroid /= area;
                      }
                    else
                      {
                        faceCentroid = RealVect::Zero;
                      }

                    // record these facts in member data
                    Faces[index].setFaceCentroid(faceCentroid);
                    Faces[index].setFaceArea(area);
                  }

                else if (Faces[index].isCovered())
                  {
                    RealVect faceCentroid = RealVect::Zero;
                    Faces[index].setFaceCentroid(faceCentroid);
                    Real area = 0.0;
                    Faces[index].setFaceArea(area);
                  }
              }
          }
#if 1
        // iterate over faces recalculating face area
        // face order is xLo,xHi,yLo,yHi,zLo,zHi
        for (int iFace = 0; iFace < 2*SpaceDim; ++iFace)
          {
            // recalculate face area
            faceMo& face = Faces[iFace];

            // collect exactly two irregular edges or mayday.
            Vector<RealVect> crossingPt;

            Vector<int> cPtHiLo;
            Vector<int> edgeHiLo;
            Vector<int> edgeDir;

            if (!(face.isCovered()) && !(face.isRegular()))
              {
                for ( int iEdge = 0; iEdge < 4; ++iEdge)
                  {
                    const edgeMo& curEdge = face.retrieveEdge(iEdge);

                    if (curEdge.dontKnow())
                      {
                        // loPt of edge
                        RealVect loPt = curEdge.getLo();

                        // hiPt of edge
                        RealVect hiPt = curEdge.getHi();

                        // direction that varies over edge
                        int direction = curEdge.direction();
                        bool intersectLo = curEdge.getIntersectLo();

                        // for irregular or regular edges at most one pt away from corner
                        if (intersectLo)
                          {
                            crossingPt.push_back(loPt);
                            cPtHiLo.push_back(-1);

                            int hilo = iEdge % 2;
                            edgeHiLo.push_back(hilo);

                            edgeDir.push_back(direction);
                          }
                        else // forced by dontknow
                          {
                            crossingPt.push_back(hiPt);
                            cPtHiLo.push_back(1);

                            int hilo = iEdge % 2;
                            edgeHiLo.push_back(hilo);

                            edgeDir.push_back(direction);
                          }
                      }
                  }
              }

            if (crossingPt.size() == 2)
              {
                // get midpoint of line connecting intersection points
                RealVect midPt = crossingPt[0];
                midPt += crossingPt[1];
                midPt *= 0.5;

                // faceNormal
                int faceNormal = Faces[iFace].getFaceNormal();

                // directions over which the face varies
                int dir1 = (faceNormal + 1) % SpaceDim;
                int dir2 = (faceNormal + 2) % SpaceDim;

                // find max of (deltaDir1,deltaDir2)
                Real deltaDir1 = std::abs(crossingPt[0][dir1] - crossingPt[1][dir1]);
                Real deltaDir2 = std::abs(crossingPt[0][dir2] - crossingPt[1][dir2]);

                // maxDir will be the direction of integration (independent variable)
                // minDir will the direction in which the integrand varies(dependent variable)
                int maxDir;
                int minDir;
                if (deltaDir1 > deltaDir2)
                  {
                    maxDir = dir1;
                    minDir = dir2;
                  }
                else
                  {
                    maxDir = dir2;
                    minDir = dir1;
                  }
                // flip area?
                bool complementArea;

                AMREX_ASSERT (cPtHiLo.size() == 2);

                // loEdge-loEdge
                if (edgeHiLo[0] == 0 && edgeHiLo[1] == 0)
                  {
                    // both crossingPts must be Hi or both must be Lo
                    AMREX_ASSERT((cPtHiLo[0] == 1 && cPtHiLo[1] == 1) ||
                           (cPtHiLo[0] == -1 && cPtHiLo[1] == -1));

                    // prismArea gives triangle

                    // if the cPtHiLo[0]= hiPt then one wants triangle
                    if (cPtHiLo[0] == 1)
                      {
                        complementArea = false;
                      }
                    else
                      {
                        complementArea = true;
                      }
                  }

                // hiEdge-hiEdge
                else if (edgeHiLo[0] == 1 && edgeHiLo[1] == 1)
                  {
                    // both crossingPts must be Hi or both must be Lo
                    AMREX_ASSERT((cPtHiLo[0] == 1 && cPtHiLo[1] == 1) ||
                           (cPtHiLo[0] == -1 && cPtHiLo[1] == -1));
                    // prismArea gives the trapezoid

                    // cPtHiLo[0] == Lo => one wants the triangle area
                    if (cPtHiLo[0] == -1)
                      {
                        complementArea = true;
                      }
                    else
                      {
                        complementArea = false;
                      }
                  }

                // hiEdge-loEdge
                else if (edgeHiLo[0] == 1 && edgeHiLo[1] == 0)
                  {
                    // cpPtHiLo must be the same or the opposite of edgeHiLo
                    AMREX_ASSERT((cPtHiLo[0] == 1 && cPtHiLo[1] == -1) ||
                           (cPtHiLo[0] == -1 && cPtHiLo[1] == 1));
                    // maxDir > minDir =>prismArea gives trapezoid
                    // maxDir < minDir =>prismArea gives triangle

                    // (cPtHiLo[1] == 1) => one wants trapezoid
                    if (cPtHiLo[1] == 1)
                      {
                        if (maxDir < minDir)
                          {
                            complementArea = true;
                          }
                        else
                          {
                            complementArea = false;
                          }
                      }
                    else
                      // (cPtHiLo[1] == -1) => one wants triangle
                      {
                        if (maxDir < minDir)
                          {
                            complementArea = false;
                          }
                        else
                          {
                            complementArea = true;
                          }
                      }
                  }

                // loEdge-hiEdge
                else if (edgeHiLo[0] == 0 && edgeHiLo[1] == 1)
                  {
                    // triangle + triangle complement or two trapezoids comprise this case

                    // two trapezoids
                    if (cPtHiLo[1] == 1 && cPtHiLo[0] == 1)
                      {
                        AMREX_ASSERT(edgeDir[0] == edgeDir[1]);
                        complementArea = false;
                      }
                    else if (cPtHiLo[1] == -1 && cPtHiLo[0] == -1)
                      {
                        AMREX_ASSERT(edgeDir[0] == edgeDir[1]);
                        complementArea = true;
                      }
                    // triangle + triangle complement
                    // if the cPtHiLo[0]= loPt then one wants triangle
                    else if (cPtHiLo[0] == -1 && cPtHiLo[1] == 1 )
                      {
                        AMREX_ASSERT(edgeDir[0] != edgeDir[1]);
                        // if maxDir < minDir prismArea gives trapezoid
                        // if maxDir > minDir prismArea gives triangle
                        if (maxDir < minDir)
                          {
                            complementArea = true;
                          }
                        else
                          {
                            complementArea = false;
                          }
                      }

                    // cPtHiLo[0]= hiPt => one wants trapezoid
                    else if (cPtHiLo[0] == 1 && cPtHiLo[1] == -1)
                      {
                        if (maxDir < minDir)
                          {
                            complementArea = false;
                          }
                        else
                          {
                            complementArea = true;
                          }
                      }
                  }

                else
                  {
                    amrex::Error("cPtDir or mindir or maxDir not set correctly");
                  }

                // segLo is the lo end of segment within face[iEdge] for Brent Rootfinder
                RealVect segLo = midPt;
                segLo[minDir] = -0.5;

                // segHi is the hi end of segment within face[iEdge] for Brent Rootfinder
                RealVect segHi = midPt;
                segHi[minDir] = 0.5;

                // put segLo and segHi in physical coordinates
                RealVect physSegLo;
                RealVect physSegHi;
                RealVect physMidPt;
                for (int idir = 0; idir < SpaceDim; ++idir)
                  {
                    physSegLo[idir] = a_dx*(segLo[idir] + a_iv[idir] + 0.5) + a_origin[idir];
                    physSegHi[idir] = a_dx*(segHi[idir] + a_iv[idir] + 0.5) + a_origin[idir];
                    physMidPt[idir] = a_dx*(midPt[idir] + a_iv[idir] + 0.5) + a_origin[idir];
                  }

                // find upDir
                // std::pair<int,Side::LoHiSide> upDir;

                // physIntercept is along the segment[physSegLo,physSegHi]
                // this segment passes through midPt with direction minDir
                Real physIntercept = 0;
                bool dropOrder = false;



                Real fLo = m_implicitFunction->value(physSegLo);
                Real fHi = m_implicitFunction->value(physSegHi);

                // This guards against the "root must be bracketed" error
                // by dropping order
                if (fLo*fHi > 0.0)
                  {
                    dropOrder = true;
                  }
                else
                  {
                    physIntercept = BrentRootFinder(physSegLo, physSegHi, minDir);
                  }

                if(!dropOrder)
                {
                  // put physIntercept into relative coordinates
                  Real intercept = physIntercept - a_origin[minDir];
                  intercept  /= a_dx;
                  intercept -= (a_iv[minDir]+0.5);

                  // push_back third pt onto crossingPt
                  crossingPt.push_back(midPt);
                  crossingPt[2][minDir] = intercept;

                  // integrate w.r.t xVec using Prismoidal Rule
                  RealVect xVec;
                  RealVect yVec;

                  // the order of (xVec,yVec) will be sorted out in PrismoidalAreaCalc
                  xVec[0] = crossingPt[0][maxDir];
                  xVec[1] = crossingPt[2][maxDir];
                  xVec[2] = crossingPt[1][maxDir];

                  yVec[0] = crossingPt[0][minDir];
                  yVec[1] = crossingPt[2][minDir];
                  yVec[2] = crossingPt[1][minDir];

                  // Prismoidal's rule
                  Real area = PrismoidalAreaCalc(xVec,yVec);

                  // Only use area if it is valid
                  if (area >= 0.0 && area <= 1.0)
                  {
                    // assign area to this value or (1 - this value)
                    if (complementArea)
                    {
                      area = 1.0 - area;
                    }

                    Faces[iFace].setFaceArea(area);
                  }
                }
              }
          }
#endif
        // fill in some arguments of computeVofInternals for the faces
        for (int faceNormal = 0;faceNormal < SpaceDim;++faceNormal)
          {
            bool coveredLo = Faces[faceNormal*2].isCovered();
            IntVect loiv = a_iv;
            loiv[faceNormal] -= 1;
            if(a_domain.contains(loiv))
              {
                coveredLo = coveredLo || (a_regIrregCovered(loiv) < 0);
              }
            if (coveredLo) 
              {
                a_loArc[faceNormal].resize(0);
                a_loFaceCentroid[faceNormal].resize(0);
                a_loAreaFrac[faceNormal].resize(0);
              }
            else if (!coveredLo)
              {
                a_loArc[faceNormal].resize(1);
                a_loFaceCentroid[faceNormal].resize(1);
                a_loAreaFrac[faceNormal].resize(1);
                IntVect otherIV = a_iv;
                otherIV[faceNormal] -= 1;

                if (a_domain.contains(otherIV))
                  {
                    //stuff that uses this does not care whether the arcs go
                    //to a regular or irregular cell--I am taking away that distinction.
                    //It only ever made sense in the context of a geometry gneration package that could
                    //generate multi-valued cells
                    a_loArc[faceNormal][0] = 0;
                  }
                else
                  {
                    // boundary arcs always -1
                    a_loArc[faceNormal][0] = -1;
                  }

                a_loFaceCentroid[faceNormal][0] = Faces[faceNormal*2].getFaceCentroid();

                a_loAreaFrac[faceNormal][0] = Faces[faceNormal*2].getFaceArea();

              }
            else
              {
                amrex::Error("is it coveredLo?");
              }
            bool coveredHi = Faces[faceNormal*2+1].isCovered();
            IntVect hiiv = a_iv;
            hiiv[faceNormal] += 1;
            if(a_domain.contains(hiiv))
              {
                coveredHi = coveredHi ||  (a_regIrregCovered(hiiv, 0) < 0);
              }
            if (coveredHi)
              {
                a_hiArc[faceNormal].resize(0);
                a_hiFaceCentroid[faceNormal].resize(0);
                a_hiAreaFrac[faceNormal].resize(0);
              }
            else if (!coveredHi)
              {
                a_hiArc[faceNormal].resize(1);
                a_hiFaceCentroid[faceNormal].resize(1);
                a_hiAreaFrac[faceNormal].resize(1);
                IntVect otherIV = a_iv;
                otherIV[faceNormal] += 1;
                if (a_domain.contains(otherIV))
                  {
                    //stuff that uses this does not care whether the arcs go
                    //to a regular or irregular cell--I am taking away that distinction.
                    //It only ever made sense in the context of a geometry gneration package that could
                    //generate multi-valued cells
                    a_hiArc[faceNormal][0] = 0;
                  }
                else if (!a_domain.contains(otherIV))
                  {
                    // boundaryArcs always -1
                    a_hiArc[faceNormal][0] = -1;
                  }
                a_hiFaceCentroid[faceNormal][0] = Faces[faceNormal*2+1].getFaceCentroid();
                a_hiAreaFrac[faceNormal][0] = Faces[faceNormal*2+1].getFaceArea();
              }
            else
              {
                amrex::Error("is it coveredHi?");
              }
          }

        // We have enough face data to construct the vof
        vofMo Vof;
        Vof.define(Faces);

        Moments geom;
        int order = 0;
        Vector<Real> answer0;
        answer0 = geom.momentCalc3D(order,Vof);

        a_volFrac = answer0[0];

        Vector<Real> answer1;
        order = 1;
        answer1 = geom.momentCalc3D(order,Vof);

        for (int idir=0; idir<SpaceDim;++idir)
          {
            a_volCentroid[idir] = answer1[SpaceDim-1-idir];
          }

        if (a_volFrac == 0.0)
          {
            a_volCentroid = RealVect::Zero;
          }
        else
          {
            a_volCentroid /= a_volFrac;

          }

        Real normalVec[SpaceDim];
        Vof.getNormal(normalVec);
        for (int idir = 0;idir < SpaceDim;++idir)
          {
            a_normal[idir] = normalVec[idir];
          }

        for (int idir = 0;idir < SpaceDim;++idir)
          {
            // (nx,ny,nz)->(nxdydz,nydxdz,nzdxdy)
            a_normal[idir] = normalVec[idir]*a_dx*a_dx;
          }
        Real anisBd = 0.0;
        for (int idir = 0;idir < SpaceDim;++idir)
          {
            anisBd += a_normal[idir]*a_normal[idir];
          }
        anisBd = sqrt(anisBd);
        if (anisBd !=0.0)
          {
            a_normal /= anisBd;
          }

        a_bndryArea = Vof.getBdArea();
        if (a_bndryArea > 0.0)
          {
            for (int idir = 0;idir < SpaceDim;++idir)
              {
                a_bndryCentroid[idir] = answer0[SpaceDim-idir]/a_bndryArea;
              }
            a_bndryArea *= anisBd;
            a_bndryArea /= (maxDx*maxDx);
          }
        else
          {
            a_bndryCentroid = RealVect::Zero;
            a_bndryArea = 0.0;
          }
      }

    // clipping
    // only report adjustments when discrepancy is above the threshold
    //    bool thisVofClipped = false;
    Real discrepancy = 0.0;
    Real volDiscrepancy = 0.0;
    // volFrac out of bounds
    if (a_volFrac < 0.0)
      {
        volDiscrepancy = std::abs(a_volFrac);
        char message[1024];
        if (SpaceDim ==2)
          {
            sprintf(message,"vol fraction (%e) out of bounds. Clipping: (%d,%d)",
                    a_volFrac,a_iv[0],a_iv[1]);
          }
        else if (SpaceDim == 3)
          {
            sprintf(message,"vol frac (%e) out of bounds. Clipping: (%d,%d,%d)",
                    a_volFrac,a_iv[0],a_iv[1],a_iv[2]);
          }
        else
          {
            sprintf(message,"SpaceDim not 2 or 3");
          }

        if (volDiscrepancy > m_threshold)
          {
            amrex::AllPrint() << std::string(message) << "\n";
          }
        // do the clipping
        //thisVofClipped = true;
        a_volFrac = 0.0;
      }

    if (a_volFrac > 1.0)
      {
        volDiscrepancy = std::abs(1.0 - a_volFrac);
        char message[1024];
        if (SpaceDim ==2)
          {
            sprintf(message,"vol fraction (%e) out of bounds. Clipping: (%d,%d)",
                    a_volFrac,a_iv[0],a_iv[1]);
          }
        else if (SpaceDim == 3)
          {
            sprintf(message,"vol frac (%e) out of bounds. Clipping: (%d,%d,%d)",
                    a_volFrac,a_iv[0],a_iv[1],a_iv[2]);
          }
        else
          {
            sprintf(message,"SpaceDim not 2 or 3");
          }
        if (volDiscrepancy>m_threshold)
          {
            amrex::AllPrint() << message << "\n";
          }
        // do the clipping
        //thisVofClipped = true;
        a_volFrac = 1.0;
      }

    // area frac out of bounds
    for (int idir = 0; idir<SpaceDim; ++idir)
      {
        for (int num = 0; num < a_loAreaFrac[idir].size();num ++)
          {
            // lo frac too high
            if (a_loAreaFrac[idir][num] > 1.0)
              {
                discrepancy = std::abs(1 - a_loAreaFrac[idir][num]);
                char message[1024];
                if (SpaceDim ==2)
                  {
                    sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d)",
                            a_loAreaFrac[idir][num],a_iv[0],a_iv[1]);
                  }
                else if (SpaceDim == 3)
                  {
                    sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                            a_loAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                  }
                else
                  {
                    sprintf(message,"SpaceDim not 2 or 3");
                  }
                if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                  {
                    amrex::AllPrint() << message << "\n";
                  }

                // do the clipping
                //thisVofClipped = true;
                a_loAreaFrac[idir][num] = 1.0;
              }
            // lo frac too low
            if (a_loAreaFrac[idir][num] < 0.0)
              {
                discrepancy = std::abs(a_loAreaFrac[idir][num]);
                char message[1024];
                if (SpaceDim ==2)
                  {
                    sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d)",
                            a_loAreaFrac[idir][num],a_iv[0],a_iv[1]);
                  }
                else if (SpaceDim == 3)
                  {
                    sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                            a_loAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                  }
                else
                  {
                    sprintf(message,"SpaceDim not 2 or 3");
                  }
                if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                  {
                    amrex::AllPrint() << message << "\n";
                  }

                // do the clipping
                //thisVofClipped = true;
                a_loAreaFrac[idir][num] = 0.0;
              }

          }
      }

    for (int idir = 0; idir<SpaceDim; ++idir)
      {
        for (int num = 0; num < a_hiAreaFrac[idir].size();num ++)
          {
            // hi frac too high
            if (a_hiAreaFrac[idir][num] > 1.0)
              {
                discrepancy = std::abs(1 - a_hiAreaFrac[idir][num]);
                char message[1024];
                if (SpaceDim ==2)
                  {
                    sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d)",
                            a_hiAreaFrac[idir][num],a_iv[0],a_iv[1]);
                  }
                else if (SpaceDim == 3)
                  {
                    sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                            a_hiAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                  }
                else
                  {
                    sprintf(message,"SpaceDim not 2 or 3");
                  }
                if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                  {
                    amrex::AllPrint() << message << "\n";
                  }

                // do the clipping
                //thisVofClipped = true;
                a_hiAreaFrac[idir][num] = 1.0;
              }
            // hi frac too low
            if (a_hiAreaFrac[idir][num] < 0.0)
              {
                discrepancy = std::abs(a_hiAreaFrac[idir][num]);
                char message[1024];
                if (SpaceDim ==2)
                  {
                    sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d)",
                            a_hiAreaFrac[idir][num],a_iv[0],a_iv[1]);
                  }
                else if (SpaceDim == 3)
                  {
                    sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                            a_hiAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                  }
                else
                  {
                    sprintf(message,"SpaceDim not 2 or 3");
                  }
                if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                  {
                    amrex::AllPrint() << message << "\n";
                  }

                // do the clipping
                //thisVofClipped = true;
                a_hiAreaFrac[idir][num] = 0.0;
              }

          }

      }

    // bndry area out of bounds
    if (a_bndryArea < 0.0)
      {
        discrepancy = std::abs(a_bndryArea);
        char message[1024];
        if (SpaceDim ==2)
          {
            sprintf(message,"boundary area fraction (%e) out of bounds. Clipping: (%d,%d)",
                    a_bndryArea,a_iv[0],a_iv[1]);
          }
        else if (SpaceDim == 3)
          {
            sprintf(message,"boundary area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                    a_bndryArea,a_iv[0],a_iv[1],a_iv[2]);
          }
        else
          {
            sprintf(message,"SpaceDim not 2 or 3");
          }
        if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
          {
            amrex::AllPrint() << message << "\n";
          }

        // do the clipping
        //thisVofClipped = true;
        a_bndryArea = 0.0;
      }

    if (a_bndryArea > sqrt(2.0))
      {
        discrepancy = std::abs(sqrt(2.0) - a_bndryArea);
        char message[1024];
        if (SpaceDim ==2)
          {
            sprintf(message,"boundary area fraction (%e) out of bounds. Clipping: (%d,%d)",
                    a_bndryArea,a_iv[0],a_iv[1]);
          }
        else if (SpaceDim == 3)
          {
            sprintf(message,"boundary area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                    a_bndryArea,a_iv[0],a_iv[1],a_iv[2]);
          }
        else
          {
            sprintf(message,"SpaceDim not 2 or 3");
          }
        if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
          {
            amrex::AllPrint() << message << "\n";
          }
        // do the clipping
        //thisVofClipped = true;
        a_bndryArea = sqrt(2.0);
      }

    // volCentroid out of bounds
    for (int idir =0;idir<SpaceDim;++idir)
      {
        if (a_volCentroid[idir] > 0.5)
          {
            discrepancy = std::abs(0.5 - a_volCentroid[idir]);
            char message[1024];
            if (SpaceDim ==2)
              {
                sprintf(message,"volCentroid (%e) out of bounds. Clipping: (%d,%d)",
                        a_volCentroid[idir],a_iv[0],a_iv[1]);
              }
            else if (SpaceDim == 3)
              {
                sprintf(message,"volCentroid(%e) out of bounds. Clipping: (%d,%d,%d)",
                        a_volCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
              }
            else
              {
                sprintf(message,"SpaceDim not 2 or 3");
              }
            if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
              {
                amrex::AllPrint() << message << "\n";

              }
            // do the clipping
            //thisVofClipped = true;
            a_volCentroid[idir] = 0.5;
          }
        if (a_volCentroid[idir] < -0.5)
          {
            discrepancy = std::abs(-0.5 - a_volCentroid[idir]);
            char message[1024];
            if (SpaceDim ==2)
              {
                sprintf(message,"volCentroid (%e) out of bounds. Clipping: (%d,%d)",
                        a_volCentroid[idir],a_iv[0],a_iv[1]);
              }
            else if (SpaceDim == 3)
              {
                sprintf(message,"volCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                        a_volCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
              }
            else
              {
                sprintf(message,"SpaceDim not 2 or 3");
              }
            if (discrepancy > m_threshold&& volDiscrepancy>m_threshold)
              {
                amrex::AllPrint() << message << "\n";
              }
            // do the clipping
            //thisVofClipped = true;
            a_volCentroid[idir] = -0.5;
          }

        // boundary centroid out of bounds
        if (a_bndryCentroid[idir] > 0.5)
          {
            discrepancy = std::abs(0.5 - a_bndryCentroid[idir]);
            char message[1024];
            if (SpaceDim ==2)
              {
                sprintf(message,"bndryCentroid (%e) out of bounds. Clipping: (%d,%d)",
                        a_bndryCentroid[idir],a_iv[0],a_iv[1]);
              }
            else if (SpaceDim == 3)
              {
                sprintf(message,"boundary Centroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                        a_bndryCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
              }
            else
              {
                sprintf(message,"SpaceDim not 2 or 3");
              }
            if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
              {
                amrex::AllPrint() << message << "\n";
              }
            // do the clipping
            //thisVofClipped = true;
            a_bndryCentroid[idir] = 0.5;
          }

        if (a_bndryCentroid[idir] < -0.5)
          {
            discrepancy = std::abs(-0.5 - a_bndryCentroid[idir]);
            char message[1024];
            if (SpaceDim ==2)
              {
                sprintf(message,"bndryCentroid (%e) out of bounds. Clipping: (%d,%d)",
                        a_bndryCentroid[idir],a_iv[0],a_iv[1]);
              }
            else if (SpaceDim == 3)
              {
                sprintf(message,"bndryCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                        a_bndryCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
              }
            else
              {
                sprintf(message,"SpaceDim not 2 or 3");
              }
            if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
              {
                amrex::AllPrint() << message << "\n";
              }
            // do the clipping
            //thisVofClipped = true;
            a_bndryCentroid[idir] = -0.5;
          }
      }
    // loFaceCentroid out of bounds
    for (int idir = 0; idir<SpaceDim; ++idir)
      {
        for (int num = 0; num < a_loFaceCentroid[idir].size();num ++)
          {
            for (int jdir = 0; jdir<SpaceDim; ++jdir)
              {
                if (a_loFaceCentroid[idir][num][jdir] > 0.5)
                  {
                    discrepancy = std::abs(0.5 - a_loFaceCentroid[idir][num][jdir]);
                    char message[1024];
                    if (SpaceDim ==2)
                      {
                        sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                      }
                    else if (SpaceDim == 3)
                      {
                        sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                      }
                    else
                      {
                        sprintf(message,"SpaceDim not 2 or 3");
                      }
                    if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
                      {
                        amrex::AllPrint() << message << "\n";

                      }
                    // do the clipping
                    //thisVofClipped = true;
                    a_loFaceCentroid[idir][num][jdir] = 0.5;
                  }
                if (a_loFaceCentroid[idir][num][jdir] < -0.5)
                  {
                    discrepancy = std::abs(-0.5 - a_loFaceCentroid[idir][num][jdir]);
                    char message[1024];
                    if (SpaceDim ==2)
                      {
                        sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                      }
                    else if (SpaceDim == 3)
                      {
                        sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                      }
                    else
                      {
                        sprintf(message,"SpaceDim not 2 or 3");
                      }
                    if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
                      {
                        amrex::AllPrint() << message << "\n";

                      }
                    // do the clipping
                    //thisVofClipped = true;
                    a_loFaceCentroid[idir][num][jdir] = -0.5;
                  }
              }
          }
      }
    // hiFaceCentroid out of bounds
    for (int idir = 0; idir<SpaceDim; ++idir)
      {
        for (int num = 0; num < a_hiFaceCentroid[idir].size();num ++)
          {
            for (int jdir = 0; jdir<SpaceDim; ++jdir)
              {
                if (a_hiFaceCentroid[idir][num][jdir] > 0.5)
                  {
                    discrepancy = std::abs(0.5 - a_hiFaceCentroid[idir][num][jdir]);
                    char message[1024];
                    if (SpaceDim ==2)
                      {
                        sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                      }
                    else if (SpaceDim == 3)
                      {
                        sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                      }
                    else
                      {
                        sprintf(message,"SpaceDim not 2 or 3");
                      }
                    if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
                      {
                        amrex::AllPrint() << message << "\n";
                      }
                    // do the clipping
                    //thisVofClipped = true;
                    a_hiFaceCentroid[idir][num][jdir] = 0.5;
                  }
                if (a_hiFaceCentroid[idir][num][jdir] < -0.5)
                  {
                    discrepancy = std::abs(0.5 - a_hiFaceCentroid[idir][num][jdir]);
                    char message[1024];
                    if (SpaceDim ==2)
                      {
                        sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                      }
                    else if (SpaceDim == 3)
                      {
                        sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                      }
                    else
                      {
                        sprintf(message,"SpaceDim not 2 or 3");
                      }
                    if (discrepancy > m_threshold && volDiscrepancy > m_threshold)
                      {
                        amrex::AllPrint() << message << "\n";
                      }
                    // do the clipping
                    //thisVofClipped = true;
                    a_hiFaceCentroid[idir][num][jdir] = -0.5;
                  }
              }
          }
      }

  }


  void GeometryShop::edgeData3D(edgeMo a_edges[4],
                                bool& a_faceCovered,
                                bool& a_faceRegular,
                                bool& a_faceDontKnow,
                                const int a_hiLoFace,
                                const int a_faceNormal,
                                const Real& a_dx,
                                const IntVect& a_iv,
                                const Box& a_domain,
                                const RealVect& a_origin) const
  {
    a_faceRegular = true;
    a_faceCovered = true;
    a_faceDontKnow = false;

    int index = -1;

    // edge order is lexigraphic xLo,xHi,yLo,yHi
    for (int dom = 0; dom < 3; ++dom)
      {
        if (dom != a_faceNormal)
          {
            for (int lohi = 0; lohi < 2; ++lohi)
              {
                // which edge 0,1,2, or 3 in lexigraphic order is given by index
                index += 1;
                // range is the direction along which the edge varies
                int range = 3 - a_faceNormal - dom;

                RealVect LoPt;
                bool LoPtChanged = false;
                Real funcLo;

                RealVect HiPt;
                bool HiPtChanged = false;
                Real funcHi;

                // put LoPt in physical coordinates
                LoPt[a_faceNormal] = a_origin[a_faceNormal]+
                  (a_iv[a_faceNormal]+a_hiLoFace)*a_dx;
                LoPt[dom] = a_origin[dom] + (a_iv[dom]+lohi)*a_dx;
                LoPt[range] = a_origin[range] + (a_iv[range])*a_dx;

                // put HiPt in physical coordinates
                HiPt[a_faceNormal] = a_origin[a_faceNormal] +
                  (a_iv[a_faceNormal]+a_hiLoFace)*a_dx;
                HiPt[dom] = a_origin[dom] + (a_iv[dom]+lohi)*a_dx;
                HiPt[range] = a_origin[range] + (a_iv[range]+1)*a_dx;

                // find the midpoint
                RealVect MidPt = LoPt;
                MidPt += HiPt;
                MidPt /= 2.0;

                Real signHi;
                Real signLo;

                // placeholders for edgeType
                bool covered  = false;
                bool regular  = false;
                bool dontKnow = false;

                //                RealVect interceptPt = RealVect::Zero;

                funcHi = m_implicitFunction->value(HiPt);
                funcLo = m_implicitFunction->value(LoPt);

                // For level set data negative -> in the fluid
                //                    positive -> out of the fluid
                signHi = -funcHi;
                signLo = -funcLo;

                edgeType(covered,regular,dontKnow,signHi,signLo);

                // now we know the boolean values so we can set the edge Hi  and Lo pts
                if (covered)
                  {
                    a_faceRegular=false;

                    LoPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_dx;
                    LoPtChanged = true;

                    HiPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_dx;
                    HiPtChanged = true;
                  }
                else if (regular)
                  {
                    a_faceCovered = false;
                  }
                else if (dontKnow)
                  {
                    a_faceRegular = false;
                    a_faceCovered = false;
                    a_faceDontKnow = true;

                    Real intercept;
                  
                    // find where the surface intersects the edge
                    intercept = BrentRootFinder(LoPt, HiPt, range);

                    if (funcHi >= 0 && funcLo*funcHi <= 0)
                      {
                        HiPt[range] = intercept;
                        HiPtChanged = true;
                      }
                    else if (funcLo >= 0 && funcLo*funcHi <= 0)
                      {
                        LoPt[range] = intercept;
                        LoPtChanged = true;
                      }
                    else
                      {
                        amrex::Error("Bogus intersection calculated");
                      }

                    IntVect p1 = a_iv + a_hiLoFace*BASISV(a_faceNormal) + lohi*BASISV(dom);
                    IntVect p2 = p1   + BASISV(range);
                    RealVect iPoint;
                    iPoint[range] = intercept;
                    pair<bool,PMapIt> bpit = SetIntersect(p1,p2,iPoint,vertCache);
                    if (bpit.first)
                    {
                      const PMapIt& pmit = bpit.second;
                      IntVect diff = p2 - p1;
                      RealVect& intersect = pmit->second;

                      intersect[a_faceNormal] = a_origin[a_faceNormal]+
                        (a_iv[a_faceNormal]+a_hiLoFace)*a_dx;
                      intersect[dom] = a_origin[dom] + (a_iv[dom]+lohi)*a_dx;
                      intersect[range] = intercept;
                    }
                  }

                // put LoPt and HiPt in local coordinates
                if (a_hiLoFace == 0)
                  {
                    LoPt[a_faceNormal] = -0.5;
                    HiPt[a_faceNormal] = -0.5;
                  }
                else
                  {
                    LoPt[a_faceNormal] =  0.5;
                    HiPt[a_faceNormal] =  0.5;
                  }

                if (lohi == 0)
                  {
                    LoPt[dom] = -0.5;
                    HiPt[dom] = -0.5;
                  }
                else
                  {
                    LoPt[dom] =  0.5;
                    HiPt[dom] =  0.5;
                  }

                if (LoPtChanged)
                  {
                    LoPt[range] -= a_origin[range];
                    LoPt[range] /= a_dx;
                    LoPt[range] -= (a_iv[range] + 0.5);
                  }
                else
                  {
                    LoPt[range] = -0.5;
                  }

                if (HiPtChanged)
                  {
                    HiPt[range] -= a_origin[range];
                    HiPt[range] /= a_dx;
                    HiPt[range] -= (a_iv[range] + 0.5);
                  }
                else
                  {
                    HiPt[range] =  0.5;
                  }

                AMREX_ASSERT((!(regular && covered)) && (!(regular && dontKnow)) && (!(dontKnow && covered)));
                AMREX_ASSERT(regular || covered || (!(LoPtChanged && HiPtChanged)));
                AMREX_ASSERT(regular || covered || dontKnow);
                AMREX_ASSERT(regular || covered || LoPtChanged || HiPtChanged);
                bool intersectLo = LoPtChanged;
                edgeMo Edge;
                // range means the coordinate direction that varies over the length of the edge
                Edge.define(LoPt,HiPt,intersectLo, range,covered,regular,dontKnow);
                a_edges[index] = Edge;
              }
          }
      }
    return;
  }
  void GeometryShop::edgeData2D(edgeMo a_edges[4],
                                bool& a_faceCovered,
                                bool& a_faceRegular,
                                bool& a_faceDontKnow,
                                const Real& a_dx,
                                const IntVect& a_iv,
                                const Box& a_domain,
                                const RealVect& a_origin) const
  {
    // index counts which edge:xLo=0,xHi=1,yLo=2,yHi=3
    int index = -1;

    a_faceRegular = true;
    a_faceCovered = true;
    a_faceDontKnow = false;

    // domain means the direction normal to the edge
    for (int domain = 0; domain < 2;++domain)
      {
        for (int lohi = 0; lohi < 2;++lohi)
          {
            index += 1;

            // range is the direction along the edge
            int range = 1-domain;

            // Express HiPt. LoPt and MidPt in physical coordinates
            RealVect LoPt;
            bool LoPtChanged = false;
            LoPt[domain] = a_origin[domain] + (a_iv[domain]+lohi)*a_dx;
            LoPt[range]  = a_origin[range]  + (a_iv[range])*a_dx;

            RealVect HiPt;
            bool HiPtChanged = false;
            HiPt[domain] = a_origin[domain] + (a_iv[domain]+lohi)*a_dx;
            HiPt[range] =  a_origin[range] + (a_iv[range]+1)*a_dx;

            RealVect MidPt = HiPt;
            MidPt += LoPt;
            MidPt /= 2.0;

            // decide which type of edge
            bool covered;
            bool regular;
            bool dontKnow;

            // function value
            Real funcHi = m_implicitFunction->value(HiPt);
            Real funcLo = m_implicitFunction->value(LoPt);

            // the sign of signHi and signLo determine edgetype
            Real signHi;
            Real signLo;

            // For level set data negative -> in the fluid
            //                    positive -> out of the fluid
            signHi = -funcHi;
            signLo = -funcLo;

            edgeType(covered,regular,dontKnow,signHi,signLo);

            // Given edgeType, set the edge Hi  and Lo pts
            if (covered)
              {
                a_faceRegular=false;

                LoPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_dx;
                LoPtChanged = true;

                HiPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_dx;
                HiPtChanged = true;
              }
            else if (regular)
              {
                a_faceCovered = false;
              }
            else if (dontKnow)
              {
                a_faceRegular  = false;
                a_faceCovered  = false;
                a_faceDontKnow = true;

                // find where the surface intersects the edge
                Real intercept;

                intercept = BrentRootFinder(LoPt, HiPt, range);

                // choose the midpoint for an ill-conditioned problem
                if (intercept<LoPt[range] || intercept>HiPt[range])
                  {
                    amrex::AllPrint()<<"GeometryShop::edgeData: Ill-conditioned edge data"<<"\n";
                    intercept = (LoPt[range]+HiPt[range])/2.0;
                  }

                if (funcHi >= 0 && funcLo*funcHi <= 0)
                  {
                    HiPt[range] = intercept;
                    HiPtChanged = true;
                  }
                else if (funcLo >= 0 && funcLo*funcHi <= 0)
                  {
                    LoPt[range] = intercept;
                    LoPtChanged = true;
                  }
                else
                  {
                    amrex::Error("Bogus intersection calculated");
                  }

		IntVect p1 = a_iv + lohi*BASISV(domain);
		IntVect p2 = p1   + BASISV(range);
		RealVect iPoint;
		iPoint[range] = intercept;
		pair<bool,PMapIt> bpit = SetIntersect(p1,p2,iPoint,vertCache);
		if (bpit.first)
		{
		  const PMapIt& pmit = bpit.second;
		  IntVect diff = p2 - p1;
		  RealVect& intersect = pmit->second;
                  intersect[domain] = a_origin[domain] + p1[domain]*a_dx;
                  intersect[range] = intercept;
		}
              }

            // express the answer relative to dx and cell-center
            if (lohi == 0)
              {
                LoPt[domain] = -0.5;
                HiPt[domain] = -0.5;
              }
            else
              {
                LoPt[domain] =  0.5;
                HiPt[domain] =  0.5;
              }

            if (LoPtChanged)
              {
                LoPt[range] -= a_origin[range];
                LoPt[range] /= a_dx;
                LoPt[range] -= (a_iv[range] + 0.5);
              }
            else
              {
                LoPt[range] = -0.5;
              }

            if (HiPtChanged)
              {
                HiPt[range] -= a_origin[range];
                HiPt[range] /= a_dx;
                HiPt[range] -= (a_iv[range] + 0.5);
              }
            else
              {
                HiPt[range] =  0.5;
              }

            AMREX_ASSERT(regular || covered || (!(LoPtChanged && HiPtChanged)));
            AMREX_ASSERT(regular || covered || dontKnow);
            AMREX_ASSERT((!(regular && covered)) && (!(regular && dontKnow)) && (!(dontKnow && covered)));
            AMREX_ASSERT(regular || covered || LoPtChanged || HiPtChanged);
            // default is something invalid
            bool intersectLo = LoPtChanged;

            // define this edge
            // Note we have some irregular edges of 0 length
            a_edges[index].define(LoPt,HiPt,intersectLo, range,covered,regular,dontKnow);
          }
      }
  }

  void GeometryShop::edgeType(bool& a_covered,
                              bool& a_regular,
                              bool& a_dontKnow,
                              Real& a_signHi,
                              Real& a_signLo) const
  {
    // if signHi and signLo are both positive
    if (a_signHi > 0.0 && a_signLo > 0.0)
      {
        a_covered  = false;
        a_regular  = true;
        a_dontKnow = false;
      }

    // if signHi and signLo are both negative
    else if (a_signHi <= 0.0 && a_signLo <= 0.0)
      {
        a_covered  = true;
        a_regular  = false;
        a_dontKnow = false;
      }

    // if signHi or signLo both are zero
    else if (a_signHi == 0.0 && a_signLo == 0.0)
      {
        a_covered  = true;
        a_regular  = false;
        a_dontKnow = false;
      }

    // otherwise signLo*signHi <= 0
    // in this case we will look for an intersection point
    else
      {
        a_covered  = false;
        a_regular  = false;
        a_dontKnow = true;
      }

    return;
  }

  //  The following is an implementation of "Brent's Method"
  //    for one-dimensional root finding. Pseudo-code for this
  //    algorithm can be found on p. 253 of "Numerical Recipes"
  //    ISBN 0-521-30811-9
  Real GeometryShop::BrentRootFinder(const RealVect& a_x1,
                                     const RealVect& a_x2,
                                     const int&      a_range) const
  {
    const Real tol = 1.0e-12;

    //  Max allowed iterations and floating point precision
    const unsigned int  MAXITER = 100;
    const Real      EPS   = 3.0e-15;

    unsigned int i;
    RealVect aPt;
    RealVect bPt;
    Real c, fa, fb, fc;
    Real d, e;
    Real tol1, xm;
    Real p, q, r, s;

    aPt = a_x1;
    bPt = a_x2;

    fa = -m_implicitFunction->value(aPt);
    fb = -m_implicitFunction->value(bPt);

    //  Init these to be safe
    c = d = e = 0.0;

    if (fb*fa > 0)
      {
        amrex::AllPrint() << "fa " << fa << " fb " << fb << "\n";
        amrex::Error("GeometryShop::BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
      }

    fc = fb;

    for (i = 0; i < MAXITER; i++)
      {
        if (fb*fc > 0)
          {
            //  Rename a, b, c and adjust bounding interval d
            c = aPt[a_range];
            fc  = fa;
            d = bPt[a_range] - aPt[a_range];
            e = d;
          }

        if (std::abs(fc) < std::abs(fb))
          {
            aPt[a_range] = bPt[a_range];
            bPt[a_range] = c;
            c = aPt[a_range];
            fa  = fb;
            fb  = fc;
            fc  = fa;
          }

        //  Convergence check
        tol1  = 2.0 * EPS * std::abs(bPt[a_range]) + 0.5 * tol;
        xm    = 0.5 * (c - bPt[a_range]);

        if (std::abs(xm) <= tol1 || fb == 0.0)
          {
            break;
          }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
          {
            //  Attempt inverse quadratic interpolation
            s = fb / fa;
            if (aPt[a_range] == c)
              {
                p = 2.0 * xm * s;
                q = 1.0 - s;
              }
            else
              {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q-r) - (bPt[a_range]-aPt[a_range]) * (r-1.0));
                q = (q-1.0) * (r-1.0) * (s-1.0);
              }

            //  Check whether in bounds
            if (p > 0) q = -q;

            p = std::abs(p);

            if (2.0 * p < std::min(3.0*xm*q-std::abs(tol1*q), std::abs(e*q)))
              {
                //  Accept interpolation
                e = d;
                d = p / q;
              }
            else
              {
                //  Interpolation failed, use bisection
                d = xm;
                e = d;
              }
          }
        else
          {
            //  Bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
          }

        //  Move last best guess to a
        aPt[a_range] = bPt[a_range];
        fa  = fb;

        //  Evaluate new trial root
        if (std::abs(d) > tol1)
          {
            bPt[a_range] = bPt[a_range] + d;
          }
        else
          {
            if (xm < 0) bPt[a_range] = bPt[a_range] - tol1;
            else        bPt[a_range] = bPt[a_range] + tol1;
          }

        fb = -m_implicitFunction->value(bPt);
      }

    if (i >= MAXITER)
      {
        std::cerr  << "BrentRootFinder: exceeding maximum iterations: "
                   << MAXITER << "\n";
      }
    //  //  Keep statistics
    //     statCount++;
    //     statSum += i;
    //     statSum2  += i*i;

    return bPt[a_range];
  }

  // Coordinates are in the box: [-0.5,0.5] x [-0.5,0.5]
  // It is assumed that the spacing in x is uniform
  Real GeometryShop::PrismoidalAreaCalc(RealVect& a_xVec,
                                        RealVect& a_yVec) const
  {
    // The area of the parabola determined by the three input points
    Real retval;

    // See if parabolic approximation, a*x^2 + b*x + c, stays in bounds
    Real aScale;  // a times 2h^2
    Real bScale;  // b times 2h
    Real cScale;  // c

    aScale = a_yVec[2] - 2*a_yVec[1] + a_yVec[0];
    bScale = a_yVec[2] - a_yVec[0];
    cScale = a_yVec[1];

    if ((std::abs(bScale) > 2*std::abs(aScale)) ||
        (std::abs(8*aScale*cScale - bScale*bScale) <= std::abs(4*aScale)))
      {
        // Compute the area of the parabola

        // Integrate w.r.t. x, y is a function of x
        // (x_vec[1],y_vec[1]) is the middle point

        // Learn whether x_vec[0] or x_vec[2] is smaller (with origin = cell center)
        Real largeX = a_xVec[2];
        Real smallX = a_xVec[0];

        // Y's that go with the largeX and smallX
        Real largeX_Y = a_yVec[2];
        Real smallX_Y = a_yVec[0];

        // Swap large for small, if necessary
        if (a_xVec[0] > a_xVec[2])
          {
            largeX = a_xVec[0];
            smallX = a_xVec[2];

            largeX_Y = a_yVec[0];
            smallX_Y = a_yVec[2];
          }

        // Find h
        Real h = 0.5 * (largeX - smallX);
        AMREX_ASSERT (h >= 0.0);

        // Prismoidal rule, which is exact on cubics.
        retval = (h/3.0) * (a_yVec[0] + 4.0 * a_yVec[1] + a_yVec[2]);

        // yVec relative to cell center. Hence add 0.5, but
        // adding 0.5 to each y in the line above adds (h/3) * 3 = h)
        retval += h;

        // Unaccounted for rectangle at the high end of the face.
        retval += (0.5 - largeX) * (largeX_Y + 0.5);

        // Similarly, the unaccounted for rectangle may be at the small end.
        retval += (smallX + 0.5) * (smallX_Y + 0.5);
      }
    else
      {
        // If the extreme of the parabola occurs between a_xVec[0] and a_xVec[2]
        // and the lies outside the current cell then return a negative area so
        // the linear approximation is used
        retval = -1.0;
      }

    return retval;
  }

}//namespace amrex
