#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <algorithm>
using std::sort;

#include "CH_Timer.H"
#include "RefCountedPtr.H"

#include "STLExplorer.H"
#include "KDTree.H"

#include "NamespaceHeader.H"

#include "STLUtil.H"
using namespace STLUtil;

/*
 * Class to explore an STL mesh
 * has member functions that
 * 1) build a K-D tree of the vertices
 * 2) search for nearest triangles given
 *  a) a point (return 1 triangle)
 *  b) a cell/volume (return all triangles with some part inside the volume)
 * 3) find intersection of a triangle with a single edge
 * 4) find intersection(s) of the whole mesh and an edge
 */

/// Constructor - just need a mesh
STLExplorer::STLExplorer(RefCountedPtr<STLMesh> a_stlmesh)
{
  CH_TIME("STLExplorer::STLExplorer");

  m_msh = a_stlmesh;
  // m_sb = NULL;
  m_freestlbox = false; // let the user free it if they created it
  m_printdebug = false;
  m_vertmap.resize( m_msh->vertices.vertex.size() );
  m_trimap.resize( m_msh->triangles.corners.size() );
  if (m_printdebug)
    m_msh->PrintMesh();
}

/// builds cellToTriangles - connectivity between box and stlmesh
void STLExplorer::Explore(RefCountedPtr<STLBox> a_sb)
{
  CH_TIME("STLExplorer::Explore_STLBox");

  // warning message if this box is going to be destroyed
  if (m_freestlbox)
  {
    pout() << endl;
    pout() << "STLExplorer::Explore Warning, the input STLBox will be deleted upon destruction.  This is because a new STLBox was created in your first call to Explore()." << endl;
  }

  m_sb = a_sb;

  DoExplore();
}

void STLExplorer::Explore(const Box&           a_region,
                          const ProblemDomain& a_domain,
                          const RealVect&      a_origin,
                          const RealVect&      a_dx)
{
  CH_TIME("STLExplorer::Explore_Box");

  // if m_sb is null, create a new STLBox.
  if (!m_sb)
  {
    RefCountedPtr<STLBox> temp(new STLBox(m_msh,a_region,a_domain,a_origin,a_dx));
    m_sb = temp;
    m_freestlbox = true; // we created the stlbox, so we will delete it as well
  }
  else
    m_sb->SetMeshBox(m_msh,a_region,a_domain,a_origin,a_dx); // if m_sb already exists, just over-write

  DoExplore();
}

void STLExplorer::DoExplore()
{
  CH_TIME("STLExplorer::DoExplore");

  //pout() << "Exploring with STLBox: box = " << m_sb->m_region << ", domain = " << m_sb->m_domain << ", orig = " << m_sb->m_origin << ", dx = " << m_sb->m_dx << "\n";
  
  // exploration divided among internal functions
  FindCellsOnVertices();

  FindCellsOnEdges();

  FindCellsInTriangles();

  RemoveCellsOutsideDomain();

  FindCellEdgesOnBondary();

  BuildKDTree();
}


STLExplorer::~STLExplorer()
{
  CH_TIME("STLExplorer::~STLExplorer");

  // get rid of the KDTree
  Vector<pair<IntVect,bool> >* pdata;
  KDGetGlobalData( m_ptree , (void **) &pdata );
  delete pdata;

  int KDError = KDFree( m_ptree );
  if (KDError!=0)
  {
    pout() << endl;
    pout() << "KDFree returned an error" << endl;
  }
  KDTreeFinalize();
}

/// return the point of intersection between an edge and the mesh
/// and whether the two nodes are inside or outside
/// a_intersectPt is bogus if both nodes are inside or both outside
void STLExplorer::GetCellEdgeIntersection(const CellEdge& a_celledge,
                                          RealVect&       a_intersectPt,
                                          bool&           a_isNode0Inside,
                                          bool&           a_isNode1Inside)
{
  CH_TIME("STLExplorer::GetCellEdgeIntersection");

  if (m_sb->m_edgemap.size()==0)
  {
    a_intersectPt = RealVect::Zero;
    a_isNode0Inside = false; // cannot determine if edge is inside or outside
    a_isNode1Inside = false;
    // pout() << "STLExplorer: Warning, no cell edges found that intersect the STLMesh" << endl;
    return;
  }
    
  //if (!m_sb->m_region.contains(a_celledge.m_node0) || !m_sb->m_region.contains(a_celledge.m_node1))
  //{
  //  pout() << "STLExplorer: warning, you are asking for an edge that is not in the box used to construct this explorer: " << a_celledge.m_node0 << "-->" << a_celledge.m_node1 << "\n";
  //}

  //FindEdgeInOut(a_celledge,a_isNode0Inside,a_isNode1Inside);
  FindEdgeInOutWithKDTree(a_celledge,a_isNode0Inside,a_isNode1Inside);
  if (a_isNode0Inside != a_isNode1Inside) // found edge on a boundary
  {
    EdgeMapIt it = m_sb->m_edgemap.find(a_celledge);

    if (it==m_sb->m_edgemap.end())
    {
      // pout() << "STLExplorer: inconsistent data, FindEdgeInOut thinks the edge is on the boundary, but the edge is not loaded into m_edgemap: " << a_celledge.m_node0 << "-->" << a_celledge.m_node1 << endl;
      //MayDay::Abort("STLExplorer: inconsistent data, FindEdgeInOut thinks the edge is on the boundary, but the edge is not loaded into m_edgemap");
      a_intersectPt = IVToRV(a_celledge.m_node0, m_sb->m_origin, m_sb->m_dx); // makes a little more sense than random data in memory
    }
    else
    {
      // return data from edgemap
      a_intersectPt = it->second;
    }
  }
  else
    a_intersectPt = RealVect::Zero; // default

}

/// return whether or not a point is inside or outside the domain
void STLExplorer::GetPointInOut(const IntVect&  a_point,
                                bool&           a_inout)
{
  CH_TIME("STLExplorer::GetPointInOut");

  //if (!m_sb->m_region.contains(a_point))
  //{
  //  pout() << "STLExplorer: warning, you are asking for a node that is not in the box used to construct this explorer: " << a_point << "\n";
  //}
  
  int KDError;

  // first convert a_point to arrays of Real's (for KDTree)
  RealVect RVpoint = IVToRV(a_point, m_sb->m_origin, m_sb->m_dx);
  Real pointArr[SpaceDim];
  for (int i=0; i<SpaceDim; i++)
    pointArr[i] = RVpoint[i];

  // set up output - result array
  Real closestArr[SpaceDim];
  for (int i=0; i<SpaceDim; i++)
    closestArr[i] = INFINITY;
  // pointer to data
  pair<IntVect,bool>* dat = NULL;

  // do search
  KDError = KDNearestNeighbor( m_ptree , pointArr , closestArr , (void**) &dat , NULL , 0 );
  if (KDError!=0)
  {
    pout() << endl;
    pout() << "KDNearestNeighbor returned an error" << endl;
  }

  RealVect closestRV;
  for (int i=0; i<SpaceDim; i++)
    closestRV[i] = closestArr[i];
  if (m_printdebug)
  {
    pout() << "In KDTree: searched for " << a_point;
    pout() << " found " << closestRV;
    pout() << " with IntVect " << dat->first << " and boolean " << dat->second << endl;
  }

  // set output - a point must be inside if the closest point in the tree is inside,
  // then the query point should be inside as well
  a_inout = dat->second;
}

/*
 * Helper functions to generate STL mesh <--> Chombo connectivity
 */

// sets map<IntVect,TriInCell> m_sb->m_cellmap (adds on cells as we find them)
// and  Vector<Intvect> m_vertmap
void STLExplorer::FindCellsOnVertices()
{
  CH_TIME("STLExplorer::FindCellsOnVertices");

  pair<IntVect,TriInCell> tmprecord;
  tmprecord.second.vertices.resize(1);
  for (int ivert = 0; ivert < m_msh->vertices.vertex.size(); ivert++)
  {
    // find IntVect of each point by essentially floor(x/dx)
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      tmprecord.first[idir] = (int) floor( (m_msh->vertices.vertex[ivert][idir]-m_sb->m_origin[idir]) / m_sb->m_dx[idir] );
    }
    /*
     * Don't want to do this yet -- what if no nodes are in the box?  still relevant!
    // check if in the box at all??
    if (!m_sb->m_region.contains(tmprecord.first))
    {
      m_vertmap[ivert] = IntVect(-1,-1,-1); // dummy value for outside of box
      continue;
    }
    */

    tmprecord.second.vertices[0]=ivert;

    // and store away in m_sb->m_cellmap
    pair<CellMapIt, bool> rvCell = m_sb->m_cellmap.insert(tmprecord);
    if (!rvCell.second)
    {
      // this cell already exists, instead add to the vector of vertices
      rvCell.first->second.vertices.push_back(ivert);
    }
    // and in m_vertmap
    m_vertmap[ivert] = tmprecord.first;

  }
  if (m_printdebug)
  {
    PMap(m_sb->m_cellmap);
    pout() << "m_vertmap: ";
    PVec(m_vertmap);
    pout() << endl;
  }

}

// goes along triangle edges and finds cells (sets left & right triangles)
// sets m_sb->m_cellmap (inserts cells) and m_trimap
void STLExplorer::FindCellsOnEdges()
{
  CH_TIME("STLExplorer::FindCellsOnEdges");

  // pick an edge
  for (int iedge = 0; iedge < m_msh->edges.edge.size(); iedge++)
  {
    // get indices of left and right triangles
    Vector<int> triLR = m_msh->connect.edgeToTriangle[iedge];

    // get cell at beginning of edge
    IntVect curcell = m_vertmap[ m_msh->edges.edge[iedge][0] ];
    /*
     * Don't check this yet... wait until we've found all relevant cells
    if (curcell<IntVect::Zero) // outside domain
      continue;
    */

    // insert triangles into curcell
    CellMapIt tmp = m_sb->m_cellmap.find(curcell);
    if (triLR[0]>=0) // don't include if negative (i.e. no triangle)
    {
      tmp->second.triangles.push_back(triLR[0]);
      m_trimap[ triLR[0] ].push_back(curcell);
    }
    if (triLR[1]>=0) // don't include if negative (i.e. no triangle)
    {
      tmp->second.triangles.push_back(triLR[1]);
      m_trimap[ triLR[1] ].push_back(curcell);
    }

    RealVect alpha; // distance along edge
    // first vertex (should we include origin here?)
    RealVect v0 = m_msh->vertices.vertex[ m_msh->edges.edge[iedge][0] ];
    // vector from v0 to v1 (should we include origin here?)
    RealVect d = m_msh->vertices.vertex[ m_msh->edges.edge[iedge][1] ] - v0;
    // signs of d - determines which octant (as if v0 was origin) we should search
    IntVect signs = RVSign(d);
    // possible edges of curcell that the edge will hit: x = dx*(i or i+1) and y = dy*(j or j+1) etc
    IntVect dcells = (signs+1)/2;

    if (m_printdebug)
    {
      pout() << "  Debug: RVSign(d)=" << RVSign(d) << ", signs=" << signs << endl;
      pout() << "  Along edge from " << v0 << " to " << d+v0 << ":" << endl;
      pout() << curcell << " --> ";
    }

    // start marching along edge
    while (true)
    {
      // get distance along edge when we get to next cell
      // (should we include origin here? yes!)
      for (int idir=0; idir<SpaceDim; idir++)
      {
        if (Abs(d[idir]) < m_msh->tol)
          alpha[idir] = INFINITY;
        else
          alpha[idir] = (m_sb->m_origin[idir] + m_sb->m_dx[idir] * (Real) (curcell[idir]+dcells[idir]) - v0[idir]) / d[idir];

        if ((alpha[idir] < m_msh->tol) || signs[idir]==0 || alpha[idir]==-0)
          alpha[idir] = INFINITY;
      }

      int newDir = alpha.minDir(false); // false just means regular minDir (no Abs())

      // check if we've gone too far or wrong direction
      if (alpha[newDir] > 1.0)
        break; // stop if we get beyond v1
      if (std::isinf(alpha[newDir]))
      {
        // actually print out what alpha was before we set it to INFINITY
        RealVect alphareal;
        for (int idir=0; idir<SpaceDim; idir++)
          alphareal[idir] = (m_sb->m_origin[idir] + m_sb->m_dx[idir] * (Real) (curcell[idir]+dcells[idir]) - v0[idir]) / d[idir];
        // pout() << "STLExplorer: bad alpha at iedge = " << iedge << endl;
        // pout() << "alpha=" << alphareal << ", d=" << d << endl;
        // pout() << "cell=" << curcell << ", v0=" << v0 << endl;
        //MayDay::Abort("STLExplorer: in FindCellOnEdge, cannot have negative distance along edge");
      }

      // move to new location
      curcell[newDir] += signs[newDir]; // actual motion follows signs along specified direction

      if (m_printdebug)
      {
        pout() << curcell << " --> ";
      }

      // add triangles to new cell
      pair<IntVect,TriInCell> currecord; currecord.first = curcell;
      pair<CellMapIt, bool> rvCell = m_sb->m_cellmap.insert(currecord);
      if (triLR[0]>=0) // don't include if negative (i.e. no triangle)
      {
        rvCell.first->second.triangles.push_back(triLR[0]);
        m_trimap[ triLR[0] ].push_back(curcell);
      }
      if (triLR[1]>=0) // don't include if negative (i.e. no triangle)
      {
        rvCell.first->second.triangles.push_back(triLR[1]);
        m_trimap[ triLR[1] ].push_back(curcell);
      }

    }
    if (m_printdebug)
    {
      pout() << " Done" << endl;
    }
    /*
    pout() << "edge " << iedge << " from ";
    PRV(m_msh->vertices.vertex[ m_msh->edges.edge[iedge][0] ]);
    pout() << " to ";
    PRV(m_msh->vertices.vertex[ m_msh->edges.edge[iedge][1] ]);
    pout() << ": \n";
    pout() << "  L:"; PVec(m_trimap[ triLR[0] ]);
    pout() << "\n";
    pout() << "  R:"; PVec(m_trimap[ triLR[1] ]);
    pout() << "\n";
    */
  }
}

// fills in any cells completely contained by the triangle
// sets m_sb->m_cellmap (inserts cells) and m_trimap
void STLExplorer::FindCellsInTriangles()
{
  CH_TIMERS("STLExplorer::FindCellsInTriangles");
  CH_TIMER("sort cells",tsort);
  CH_TIMER("fill line",tline);
  CH_TIMER("insert cells",tinsert);
  CH_TIMER("remove duplicate cells",tunique);

  if (SpaceDim>3)
    MayDay::Abort("STLExplorer: FindCellsInTriangles not implemented for greater than 3D");

  vector<IntVect> cells;
  for (int itri = 0; itri < m_msh->triangles.corners.size(); itri++)
  {
    // sort the cells by components of their intvects
    cells = m_trimap[itri].stdVector();

    /*
    //debug stuff
    RealVect node0 = m_msh->vertices.vertex[ m_msh->triangles.corners[itri][0] ];
    RealVect node1 = m_msh->vertices.vertex[ m_msh->triangles.corners[itri][1] ];
    RealVect node2 = m_msh->vertices.vertex[ m_msh->triangles.corners[itri][2] ];
    pout() << "tris(:,:," << itri+1 << ") = [" << node0[0] << "," << node0[1] << "," << node0[2] << ";" \
         << node1[0] << "," << node1[1] << "," << node1[2] << ";" \
         << node2[0] << "," << node2[1] << "," << node2[2] << "];\n";
         */

    for (int idir = 0; idir < SpaceDim; idir++)
    //int idir=0; // just 1 iteration
    {

      CH_START(tsort);
      // sort along direction idir first
      IVCompareSWO comparator(idir);
      sort ( cells.begin() , cells.end() , comparator );
      CH_STOP(tsort);

      CH_START(tunique);
      // remove duplicates
      vector<IntVect>::iterator it = unique( cells.begin(), cells.end() );
      cells.resize( it - cells.begin() );
      CH_STOP(tunique);

      //debug stuff
      //Vector<IntVect> tmpprint(cells);
      //pout() << "  Filling in cells in " << idir << " direction: \n";
      //PVec(tmpprint); pout() << "\n";

      int i=0;
      IntVect tmp;
      while (true)
      {

        if (i>=cells.size())
          break;

        if (i==0)
        {
          i+=1; continue;
        }

        CH_START(tline);

        // only fill-in sequential points along a single cartesian direction
        // go to next slice of cells
        if (SpaceDim>0 && (cells[i][idir] != cells[i-1][idir]))
        {
          // fill in missing line(s) or 
          // fill in cells when we jump position from one line to the next 
          // don't do this only if i'th cell is adjacent to (i-1)'th, only 
          // in a new line (so idir value is incremented, idir+1 is the same)
          // note, i is incremented in FillInCellLine
          if (SpaceDim==2 && \
              !( cells[i][idir] == cells[i-1][idir]+1 && \
                 cells[i][(idir+1)%SpaceDim] == cells[i-1][(idir+1)%SpaceDim] ))
            FillInCellLine(cells,i,itri,idir,(idir+1)%SpaceDim);
          i += 1;
          CH_STOP(tline);
          continue;
        }

        // go to next line of cells
        if (SpaceDim>1 && (cells[i][(idir+1)%SpaceDim] != cells[i-1][(idir+1)%SpaceDim]))
        {
          // fill in missing line(s) or 
          // fill in cells when we jump position from one line to the next 
          // don't do this only if i'th cell is adjacent to (i-1)'th, only 
          // in a new line (so idir+1 value is incremented, idir+2 is the same)
          // note, i is incremented in FillInCellLine
          if (SpaceDim==3 && \
              !( cells[i][(idir+1)%SpaceDim] == cells[i-1][(idir+1)%SpaceDim]+1 && \
                 cells[i][(idir+2)%SpaceDim] == cells[i-1][(idir+2)%SpaceDim] ))
            FillInCellLine(cells,i,itri,(idir+1)%SpaceDim,(idir+2)%SpaceDim);
          i += 1;
          CH_STOP(tline);
          continue;
        }

        CH_STOP(tline);

        CH_START(tinsert);

        int dirfill = (idir+SpaceDim-1)%SpaceDim;

        // fill in this line of cells
        if (cells[i][dirfill] != cells[i-1][dirfill]+1)
        {
          // if both x and y values match with the previous cell,
          // and the z values do not go up by 1,
          // fill in the intervening cells so that they do
          tmp = cells[i-1];
          int ninsert = cells[i][dirfill] - cells[i-1][dirfill] - 1;
          for (int k = 0; k < ninsert; k++)
          {
            tmp[dirfill] += 1; // increment the z-value
            // and insert before i'th value, push old i'th value out to i+k
            cells.insert( cells.begin()+i+k , tmp );
            // add to cell map as well
            TriInCell tmptic; tmptic.triangles.push_back(itri);
            pair<CellMapIt, bool> rvCell = m_sb->m_cellmap.insert( make_pair(tmp,tmptic) );
            if (!rvCell.second)
              rvCell.first->second.triangles.push_back(itri);
          }

          i += ninsert; // reset iterator for inserted elements

        }

        CH_STOP(tinsert);

        i+=1; // regular increment

      }
    }

    m_trimap[itri] = cells; // put back in m_trimap
  }

  CH_START(tunique);

  // now remove duplicate triangle entries - fix up m_sb->m_cellmap
  vector<int> tmptri;
  for (CellMapIt it = m_sb->m_cellmap.begin(); it != m_sb->m_cellmap.end(); it++)
  {
    it->second.triangles.sort();
    tmptri = it->second.triangles.stdVector();
    vector<int>::iterator it2 = unique(tmptri.begin(),tmptri.end());
    tmptri.resize( it2 - tmptri.begin() );
    it->second.triangles = tmptri; // re-assign
  }

  CH_STOP(tunique);
    
  if (m_printdebug)
  {
    int itri = 102; IntVect iv(D_DECL(33,32,13));
    CellMapIt it = m_sb->m_cellmap.find(iv); 
    pout() << " After fill-in:"; PVec(m_trimap[itri]);
    pout() << endl;
    pout() << " Found IntVect " << iv << "? " << (it!=m_sb->m_cellmap.end()) << endl;
  }

  if (m_printdebug)
  {
    pout() << "Filled in cells." << endl;
    //PMap(m_sb->m_cellmap);
    //pout() << " now triangle index -> cells:\n";
    //PVec(m_trimap);
    //pout() << "\n";
  }

}

// removes all cells that lie outside of the Box from consideration
// modifies m_sb->m_cellmap only (call before setting m_sb->m_edgemap and m_sb->m_nodemap)
void STLExplorer::RemoveCellsOutsideDomain()
{
  CH_TIME("STLExplorer::RemoveCellsOutsideDomain");

  // tricky to erase from a map... see qualapps.blogspot.com for this sol'n
  CellMapIt it = m_sb->m_cellmap.begin();
  while ( it != m_sb->m_cellmap.end() )
  {
    if (!m_sb->m_region.contains(it->first))
    {
      CellMapIt next = it;
      ++next;

      // if the current cell is not in the box, erase it
      // note, we are not modifying m_vertmap or m_trimap!!
      // they are only temporary data holders and may contain IntVects
      // that lie outside of m_box
      m_sb->m_cellmap.erase(it);

      it = next;
    }
    else
    {
      ++it;
    }
  }

}

// stores the set of cell edges on the boundary
// builds m_sb->m_edgemap and m_sb->m_nodemap (from m_sb->m_cellmap)
void STLExplorer::FindCellEdgesOnBondary()
{
  CH_TIMERS("STLExplorer::FindCellEdgesOnBoundary");
  CH_TIMER("check edge",tcheck);
  CH_TIMER("compute intersections",tinters);
  CH_TIMER("insert into nodemap and edgemap",tinsert);

  // list of increments to the cell IntVect to get all of the edges in the cell
  Vector<IntVect> edgeinc(0);
  Vector<int> edgedir(0); // list of directions of each edge in edgeinc

//  if (SpaceDim==3)
//  {
#if CH_SPACEDIM == 3
    edgeinc.push_back(IntVect(0,0,0)); edgedir.push_back(0); // along x-dir
    edgeinc.push_back(IntVect(0,0,1)); edgedir.push_back(0);
    edgeinc.push_back(IntVect(0,1,0)); edgedir.push_back(0);
    edgeinc.push_back(IntVect(0,1,1)); edgedir.push_back(0);
    edgeinc.push_back(IntVect(0,0,0)); edgedir.push_back(1); // along y-dir
    edgeinc.push_back(IntVect(0,0,1)); edgedir.push_back(1);
    edgeinc.push_back(IntVect(1,0,0)); edgedir.push_back(1);
    edgeinc.push_back(IntVect(1,0,1)); edgedir.push_back(1);
    edgeinc.push_back(IntVect(0,0,0)); edgedir.push_back(2); // along z-dir
    edgeinc.push_back(IntVect(0,1,0)); edgedir.push_back(2);
    edgeinc.push_back(IntVect(1,0,0)); edgedir.push_back(2);
    edgeinc.push_back(IntVect(1,1,0)); edgedir.push_back(2);
//  }
//  else if (SpaceDim==2)
//  {
#endif
#if CH_SPACEDIM == 2
    edgeinc.push_back(IntVect(0,0)); edgedir.push_back(0); // along x-dir
    edgeinc.push_back(IntVect(0,1)); edgedir.push_back(0);
    edgeinc.push_back(IntVect(0,0)); edgedir.push_back(1); // along y-dir
    edgeinc.push_back(IntVect(1,0)); edgedir.push_back(1);
#endif
//  }
//  else
//  {
#if CH_SPACEDIM > 3
    MayDay::Abort("STLEXplorer::FindCellEdgesOnBoundary only implemented for 2D and 3D");
#endif
//  }

  //pout() << "STLExplorer::FindCellEdgesOnBoundary: looking for edges in " << m_sb->m_cellmap.size() << " cells...\n";

  Vector<int> tris(0);
  Vector<RealVect> pts(0); 
  CellEdge curedge(m_sb->m_cellmap.begin()->first,0); // initialize
  int nedgeMultipleIntersections = 0;
  int nnodeOnTriangle = 0;
  int nedgeOnTriangle = 0;
  //int printtri = 0; // debug stuff
  for (CellMapIt it = m_sb->m_cellmap.begin(); it != m_sb->m_cellmap.end(); it++)
  {
    /*
    //debug stuff
    if (it->first == IntVect(33,32,13))
    {
      pout() << "Found odd cell...";
      printtri=1;
    }
    else
      printtri=0;
    RealVect c0 = m_msh->vertices.vertex[ m_msh->triangles.corners[ it->second.triangles[0] ][0] ];
    RealVect c1 = m_msh->vertices.vertex[ m_msh->triangles.corners[ it->second.triangles[0] ][1] ];
    RealVect c2 = m_msh->vertices.vertex[ m_msh->triangles.corners[ it->second.triangles[0] ][2] ];
    if (printtri>=1)
    {
      printtri=2;
      pout() << "tri = [ " << c0[0] << "," << c0[1] << "," << c0[2] << ";" \
           << c1[0] << "," << c1[1] << "," << c1[2] << ";" \
           << c2[0] << "," << c2[1] << "," << c2[2] << "];\n";
    }
    */

    for (int iedgel = 0; iedgel < edgeinc.size(); iedgel++)
    {
      CH_START(tcheck);

      // put enough info in curedge (note, we haven't updated cell or lohi)
      curedge.m_dir = edgedir[iedgel];
      curedge.m_node0 = it->first + edgeinc[iedgel]; // iterate over all edges in cell
      curedge.m_node1 = curedge.m_node0; curedge.m_node1[curedge.m_dir]++;

      // don't do anything if we already have the edge/nodes
      EdgeMapIt ittmp = m_sb->m_edgemap.find(curedge);
      
      CH_STOP(tcheck);

      if (ittmp != m_sb->m_edgemap.end())
        continue;

      CH_START(tinters);

      // check intersection of all triangles in the cell with the edge
      tris.resize(0); pts.resize(0);
      for (int itri = 0; itri < it->second.triangles.size(); itri++)
      {
        RealVect inter = FindPlaneLineIntersection(curedge,it->second.triangles[itri]);
        
        /*
        //debug stuff
        bool t0 = IsPointInTriangle(inter,it->second.triangles[itri]);
        bool t1 = IsPointOnCellEdge(inter,curedge);
        Vector<int> t2 = it->second.triangles; // triangles associated with this cell
        //RealVect c0 = m_msh->vertices.vertex[ m_msh->triangles.corners[ it->second.triangles[itri] ][0] ];
        //RealVect c1 = m_msh->vertices.vertex[ m_msh->triangles.corners[ it->second.triangles[itri] ][1] ];
        //RealVect c2 = m_msh->vertices.vertex[ m_msh->triangles.corners[ it->second.triangles[itri] ][2] ];
        RealVect node0 = IVToRV(curedge.m_node0, m_sb->m_origin, m_sb->m_dx);
        RealVect node1 = IVToRV(curedge.m_node1, m_sb->m_origin, m_sb->m_dx);
        if (printtri>=2 && itri==0)
        {
          pout() << "edges(:,:," << iedgel+1 << ") = [" << node0[0] << "," << node0[1] << "," << node0[2] << ";" << node1[0] << "," << node1[1] << "," << node1[2] << "];\n";
          if (iedgel >= edgeinc.size())
            printtri=0;
        }
        */

        if (IsPointInTriangle(inter,it->second.triangles[itri]) && IsPointOnCellEdge(inter,curedge))
        {
          // some checks for poorly conditioned geometry
          if (IsPointInTriangle(IVToRV(curedge.m_node0, m_sb->m_origin, m_sb->m_dx),it->second.triangles[itri]))
            nnodeOnTriangle++;
            
          if (IsPointInTriangle(IVToRV(curedge.m_node1, m_sb->m_origin, m_sb->m_dx),it->second.triangles[itri]))
            nnodeOnTriangle++;

          // this is actually a good check because we already know that inter is in the triangle
          // note, we will not catch cell edges that lie in the plane of a triangle but where 
          // the midpoint of the edge is not in the triangle (e.g. the edge cuts into the tip of a 
          // triangle at 10% and 20% along the edge)
          if (Abs(m_msh->triangles.normal[it->second.triangles[itri]][curedge.m_dir]) < 1.0e-8)
            nedgeOnTriangle++;
          
          tris.push_back(it->second.triangles[itri]);
          pts.push_back(inter);
        }
      }
      
      CH_STOP(tinters);

      // now tris is a list of triangles that intersect this edge
      // if none, just skip the edge
      // first check if we found a point on an STLMesh edge -
      // i.e. we found 2 triangles that are adjacent

      /*
       * This was used to deal (a little bit) with edges that have 
       * more than one intersection point.  Now, we just take the median
       * of the points (odd # of intersections) or declare zero 
       * intersections (even # of intersections), to be consistent with
       * Chombo libraries
       */
      /*
      if (tris.size()==2 && (pts[0]-pts[1]).vectorLength() < m_msh->tol)
      {
        // if the intersection points are the same
        tris.resize(1);
      }

      if (tris.size()==2)
      {
        // if the two triangles share an edge.  Hopefully this will never be called,
        // since the previous block is much faster and easier!
        // the edge could have a flipped entry in edgeToTriangle
        Vector<int> tris2(2); tris2[0]=tris[1]; tris2[1]=tris[0];
        // search through all edges to see if these triangles share an edge
        int iedge;
        for (iedge = 0; iedge < m_msh->connect.edgeToTriangle.size(); iedge++)
        {
          if (m_msh->connect.edgeToTriangle[iedge].stdVector() == tris.stdVector() || \
              m_msh->connect.edgeToTriangle[iedge].stdVector() == tris2.stdVector())
          {
            pout() << "STLExplorer: Building boundary nodes: Warning, two intersection points were found greater than STLMesh::tol apart, but an STL edge between to adjacent triangles was found!\n";
            tris.resize(1); // we can safely ignore one of the triangles
            break;
          }
        }
      }

      if (tris.size()==1)
      */

      if ( tris.size() % 2 == 1 ) // odd number of intersections
      {

        int itri; // index of relevant intersection (tris and pts)
        if (tris.size()==1)
          itri = 0;
        else
        {

          // check if all of the intersection points we found are really the same point
          bool areAllPtsEqual=true;
          for (int i=1; i<pts.size(); i++)
            areAllPtsEqual = areAllPtsEqual && ((pts[i]-pts[0]).vectorLength() < m_msh->tol);

          if (areAllPtsEqual)
            itri = 0;
          else
          {
            nedgeMultipleIntersections++;
            // force to 1 intersection by taking median intersection point
            RealVect node0 = IVToRV(curedge.m_node0, m_sb->m_origin, m_sb->m_dx);

            // sort based on dist. along edge, find median triangle index
            vector<pair<Real,int> > edgeintmap; edgeintmap.resize(tris.size());
            for (int i=0; i<tris.size(); i++)
              edgeintmap[i] = make_pair((pts[i]-node0).vectorLength() , i);
            RealIntCompare comparator;
            sort(edgeintmap.begin(), edgeintmap.end(), comparator);
            // watch out, this is integer division, but it's what we want - the median entry
            itri = edgeintmap[ edgeintmap.size()/2 ].second;

            //map<Real,int> edgeintmap; // distance along edge -> triangle index
            //for (int i=0; i<tris.size(); i++)
            //  edgeintmap.insert(make_pair( (pts[i]-node0).vectorLength() , tris[i] ));
            //// watch out, this is integer division, but it's what we want - the median entry
            //itri = ( edgeintmap.begin() + ((int) edgeintmap.size())/2 )->second; 
          }
        }

        CH_START(tinsert);

        // good condition, only 1 intersection point
        // now find out which way is inside/outside
        bool isNode1Inside = WhichNodeIsInside(curedge,tris[itri]);

        // put node0 in the m_sb->m_nodemap
        pair<NodeMapIt, bool>  rvNode = m_sb->m_nodemap.insert(make_pair( curedge.m_node0 ,  !isNode1Inside ));
        if (!rvNode.second && (rvNode.first->second != (!isNode1Inside) ))
        {
          // node was not inserted and the existing node has the opposite in/out'ness
          pout() << endl;
          pout() << "STLExplorer: Building boundary nodes: Warning, a node is specified as inside for one edge but outside for another!" << endl;
          pout() << curedge.m_node0 << " at " << IVToRV(curedge.m_node0, m_sb->m_origin, m_sb->m_dx) << endl; 
        }

        // put node1 in the m_sb->m_nodemap
        rvNode = m_sb->m_nodemap.insert(make_pair( curedge.m_node1 , isNode1Inside ));
        if (!rvNode.second && (rvNode.first->second != (isNode1Inside) ))
        {
          // node was not inserted and the existing node has the opposite in/out'ness
          pout() << endl;
          pout() << "STLExplorer: Building boundary nodes: Warning, a node is specified as inside for one edge but outside for another!\n";
          pout() << curedge.m_node1 << " at " << IVToRV(curedge.m_node1, m_sb->m_origin, m_sb->m_dx) << endl;
        }

        // and store away the intersection point between the edge and the triangle
        m_sb->m_edgemap.insert(make_pair(curedge,pts[itri]));

        /*
        //pout() << "Inserting edge and nodes: " << curedge.m_node0 << "-->" << curedge.m_node1 << "\n";
        RealVect node0 = IVToRV(curedge.m_node0, m_sb->m_origin, m_sb->m_dx);
        RealVect node1 = IVToRV(curedge.m_node1, m_sb->m_origin, m_sb->m_dx);
        pout() << "edges(:,:,end+1) = [" << node0[0] << "," << node0[1] << "," << node0[2] << ";" \
             << node1[0] << "," << node1[1] << "," << node1[2] << "];\n";
        pout() << "% edge: " << curedge.m_node0 << "-->" << curedge.m_node1 << ", node " << isNode1Inside << " inside domain\n";
        */
        CH_STOP(tinsert);

      }
      else if (tris.size()>1)
      { 
        // check if all of the intersection points we found are really the same point
        bool areAllPtsEqual=true;
        for (int i=1; i<pts.size(); i++)
          areAllPtsEqual = areAllPtsEqual && ((pts[i]-pts[0]).vectorLength() < m_msh->tol);

        if (!areAllPtsEqual)
        {
          // even number of intersections, more than 2
          nedgeMultipleIntersections++;
        }
      }

      /*
      if (tris.size()>1)
      {
        pout() << "STLExplorer: Building boundary nodes: Warning, a Chombo edge has multiple intersections with the embedded boundary.  We cannot deal with this yet.  Please make a finer mesh!\n";
        // print intersection points
        pout() << " Intersection points: "; PVec(pts); pout() << "\n";
      }
      */

    }
  }
  if (m_printdebug)
    PMap(m_sb->m_nodemap);
  if (nedgeMultipleIntersections>0)
  {
    pout() << endl;
    pout() << "STLExplorer: Building boundary nodes: Warning, there were " << nedgeMultipleIntersections << " edges that had " << endl << "  multiple intersection points.  You should consider using a finer mesh!" << endl;
  }
  if (nnodeOnTriangle>0)
  {
    pout() << endl;
    pout() << "STLExplorer: Building boundary nodes: Warning, there were approximately " << nnodeOnTriangle << " nodes that lie in the plane of a triangle. " << endl << " This may give unexpected results.  Try perturbing the geometry." << endl; 
  }
  if (nedgeOnTriangle>0)
  {
    pout() << endl;
    pout() << "STLExplorer: Building boundary nodes: Warning, there were " << nedgeOnTriangle << " edges that lie in the plane of a triangle. " << endl << " This may give unexpected results.  Try perturbing the geometry." << endl; 
  }

}

typedef pair<RealVect,pair<IntVect,bool>*> BUILDNODE;

void STLExplorer::BuildKDTree()
{
  CH_TIME("STLExplorer::BuildKDTree");

  // build KDTree of nodes in m_nodemap for later searching
  int KDError;
  m_ptree = (KDTree *) KDCreate( SpaceDim , &KDError );
  if (KDError!=0)
  {
    pout() << endl;
    pout() << "KDCreate returned an error" << endl;
  }

  vector<pair<IntVect,bool> >* data = new vector<pair<IntVect,bool> >(m_sb->m_nodemap.size());

  KDSetGlobalData( m_ptree , (void *) data );

  // don't bother inserting if there are no edges
  if (m_sb->m_edgemap.size()==0)
    return;

#if 1
  vector<BUILDNODE> allNodes(m_sb->m_nodemap.size());

  int i = 0;
  for (NodeMapIt it = m_sb->m_nodemap.begin(); it != m_sb->m_nodemap.end(); it++)
  {
    RealVect nodeRV;
    nodeRV = IVToRV(it->first, m_sb->m_origin, m_sb->m_dx);

    allNodes[i].first = nodeRV;

    (*data)[i].first = it->first;
    (*data)[i].second = it->second;
    
    allNodes[i].second = &(*data)[i];
    i++;
  }

  STLExplorer::RecursiveKDTreeInsert(allNodes,0,allNodes.size()-1,0);
#else
  // insert each node and it's associated boolean
  Real nodeArr[SpaceDim];
  RealVect nodeRV;
  int i = 0;
  for (NodeMapIt it = m_sb->m_nodemap.begin(); it != m_sb->m_nodemap.end(); it++)
  {
    nodeRV = IVToRV(it->first, m_sb->m_origin, m_sb->m_dx);
    for (int j=0; j<SpaceDim; j++)
      nodeArr[j] = nodeRV[j];
    (*data)[i].first = it->first;
    (*data)[i].second = it->second;
    KDError = KDInsert( m_ptree , nodeArr , &(*data)[i] );
    if (KDError!=0)
    {
      pout() << endl;
      pout() << "KDInsert returned an error" << endl;
    }
    i++;
  }
#endif

  KDTreeStatistics(m_ptree);
  pout() << endl;
  // KDTreePrint(m_ptree);
}

class SortNodes
{
  public:
    SortNodes(const int & a_dim)
    { 
      m_dim = a_dim;
    }

    bool operator()(BUILDNODE m_a,
                    BUILDNODE m_b)
    {
      const RealVect &posA = m_a.first;
      const RealVect &posB = m_b.first;

      if (posA[m_dim] < posB[m_dim])
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  protected:
    int m_dim;
};

void STLExplorer::RecursiveKDTreeInsert(vector<BUILDNODE> &allNodes,
                                        const int         &nstart,
                                        const int         &nend,
                                        const int         &depth)
{
  if (nstart <= nend)
  {
    int dim = depth % SpaceDim;

    vector<BUILDNODE>::iterator istart = allNodes.begin() + nstart  ;
    vector<BUILDNODE>::iterator iend   = allNodes.begin() + nend + 1;

#if 0
    pout() << "RecursiveKDTreeInsert" << endl;
    pout() << "  nstart = " << nstart << endl;
    pout() << "  nend   = " << nend   << endl;
    pout() << "  npivot = " << npivot << endl;
    pout() << "  depth  = " << depth  << endl;
    pout() << "  dim    = " << dim    << endl;

    pout() << "  input vector:" << endl;
    for (vector<BUILDNODE>::iterator i = istart; i != iend; i++)
    {
      pout() << "    " << i->first << endl;
    }
#endif

    struct SortNodes sortNodes(dim);

    sort(istart,iend,sortNodes);

#if 0
    pout() << "  sorted vector:" << endl;
    for (vector<BUILDNODE>::iterator i = istart; i != iend; i++)
    {
      pout() << "    " << i->first << endl;
    }
    pout() << endl;
#endif

    int npivot;

    for (npivot = (nend + nstart)/2; npivot < nend; npivot++)
    {
      if (allNodes[npivot].first[dim] != allNodes[npivot+1].first[dim])
      {
        break;
      }
    }

    Real nodeArr[SpaceDim];

    for (int i = 0; i < SpaceDim; i++)
    {
      nodeArr[i] = allNodes[npivot].first[i];
    }

    int KDError;

#if 0
    pout() << "  inserting: " << allNodes[npivot].first << endl;
    pout() << endl;
#endif

    KDError = KDInsert(m_ptree,nodeArr,allNodes[npivot].second);

    if (KDError!=0)
    {
      pout() << "KDInsert returned an error" << endl;
    }

#if 0
    KDTreePrint(m_ptree);
    pout() << endl;
#endif

    RecursiveKDTreeInsert(allNodes,nstart  ,npivot-1,depth+1);
    RecursiveKDTreeInsert(allNodes,npivot+1,nend    ,depth+1);
  }
}

/*
 * Functions for low-level computations
 */

// given a triangle and an edge, find the intersection
RealVect STLExplorer::FindPlaneLineIntersection(const CellEdge& celledge,
                                                const int&      triangle)
{
  RealVect normal  = m_msh->triangles.normal[triangle];
  RealVect corner0 = m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][0] ];
  RealVect node0   = IVToRV(celledge.m_node0, m_sb->m_origin, m_sb->m_dx);

  /*
   * Eqn of the plane of the triangle is
   * dot( n , (x - x_t) ) = 0
   * where n is the normal and x_t is a point on the plane (here corner0)
   * So to find intersection along the edge, we fix y and z (say) and find
   * the x that solves the above equation.  In this case, we would have
   * x = x_t - ( n_y*(y_n-y_t) + n_z*(z_n-z_t) )/n_x
   * where []_n are the fixed y and z values of the two nodes
   * the formula is suitably modified for edges pointing in the y or z directions
   */
  Real ndotdeltax = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    if (idir != celledge.m_dir) // do dot product: normal*(x_node - x_corner), except for dir
      ndotdeltax += normal[idir] * (node0[idir] - corner0[idir]);
  }

  // start at node and move over to plane
  RealVect intersectPt = node0;
  intersectPt[ celledge.m_dir ] = corner0[ celledge.m_dir ] - ndotdeltax/normal[ celledge.m_dir ];

  // bad conditioning if normal[dir]~0 (edge lies in plane of triangle)
  if ( Abs(normal[celledge.m_dir]) < 1.0e-8)
  {
    //pout() << "STLExplorer::FindPlaneLineIntersection: Warning, poorly conditioned for triangle " << triangle << " and edge " << celledge.m_node0 << "-->" << celledge.m_node1 << "\n";
    // just output midpoint of edge
    intersectPt[celledge.m_dir] = node0[celledge.m_dir] + 0.5*m_sb->m_dx[celledge.m_dir]; 
  }
  return intersectPt;

}

// returns whether node0 or node1 of the cell edge is inside the domain
bool STLExplorer::WhichNodeIsInside(const CellEdge& celledge,
                                    const int&      triangle)
{
  // vector from 0 to 1
  RealVect vedge = IVToRV(celledge.m_node1, m_sb->m_origin, m_sb->m_dx) - IVToRV(celledge.m_node0, m_sb->m_origin, m_sb->m_dx);

  /*
  if (m_printdebug && celledge.m_node0 == IntVect(D_DECL(4,7,0)) && celledge.m_node1 == IntVect(D_DECL(5,7,0)))
  {
    pout() << "  WhichNodeIsIniside:" << endl;
    pout() << "   node0 = " << celledge.m_node0 << ", node1 = " << celledge.m_node1 << endl;
    pout() << "   edge vect = " << vedge << ", normal = " << (m_msh->triangles.normal[triangle]) << endl;
  }
  */

  if (vedge.dotProduct(m_msh->triangles.normal[triangle])<=0)
    return 0;
  else
    return 1;
  /*{
    pout() << "Warning, unknown which node on this edge is inside/outside of the domain\n";
    pout() << " Two endpoints are: "; PRV(node0); pout() << " and "; PRV(node1); pout() << "\n";
    return 0;
  }*/
}

// returns whether node0 and node1 
void STLExplorer::FindEdgeInOut(const CellEdge& celledge,
                                bool&           isNode0Inside,
                                bool&           isNode1Inside)
{
  CH_TIME("STLExplorer::FindEdgeInOut");

  // temporary map for distance between a point and a Chombo node
  map<Real,IntVect> distmap0;
  map<Real,IntVect> distmap1;
  RealVect d;

  RealVect node0 = IVToRV(celledge.m_node0, m_sb->m_origin, m_sb->m_dx);
  RealVect node1 = IVToRV(celledge.m_node1, m_sb->m_origin, m_sb->m_dx);

  // check if node0|1 is inside/outside by finding the closest node
  // from the list of nodes adjacent to the boundary
  //
  // This should be done with a KD tree.  For now, brute force it
  //
  // first compute distance to every node on boundary
  for (NodeMapIt it = m_sb->m_nodemap.begin(); it != m_sb->m_nodemap.end(); it++)
  {
    d = IVToRV(it->first, m_sb->m_origin, m_sb->m_dx) - node0;
    distmap0.insert(make_pair(d.vectorLength(),it->first));
    d = IVToRV(it->first, m_sb->m_origin, m_sb->m_dx) - node1;
    distmap1.insert(make_pair(d.vectorLength(),it->first));
  }

  // retrieve inside/outside information for the closest nodes
  NodeMapIt it0 = m_sb->m_nodemap.find( distmap0.begin()->second );
  isNode0Inside = it0->second;
  NodeMapIt it1 = m_sb->m_nodemap.find( distmap1.begin()->second );
  isNode1Inside = it1->second;

}

// same as above, but uses a KDTree instead of exhaustive search
void STLExplorer::FindEdgeInOutWithKDTree(const CellEdge& celledge,
                                          bool&           isNode0Inside,
                                          bool&           isNode1Inside)
{
  CH_TIME("STLExplorer::FindEdgeInOutWithKDTree");

  int KDError;

  // first convert celledge.m_node0|1 to arrays of Real's (for KDTree)
  RealVect node0 = IVToRV(celledge.m_node0, m_sb->m_origin, m_sb->m_dx);
  RealVect node1 = IVToRV(celledge.m_node1, m_sb->m_origin, m_sb->m_dx);
  Real node0Arr[SpaceDim];
  Real node1Arr[SpaceDim];
  for (int i=0; i<SpaceDim; i++)
  {
    node0Arr[i] = node0[i];
    node1Arr[i] = node1[i];
  }

  // set up output - result array
  Real closestArr[SpaceDim];
  for (int i=0; i<SpaceDim; i++)
    closestArr[i] = INFINITY;
  // two pointers to bools
  pair<IntVect,bool>* inout0 = NULL;
  pair<IntVect,bool>* inout1 = NULL;
  //bool* inout0 = NULL;
  //bool* inout1 = NULL;
  // and distance between points
  //Real dist;

  // do search for node0
  KDError = KDNearestNeighbor( m_ptree , node0Arr , closestArr , (void**) &inout0 , NULL , 0 );
  if (KDError!=0)
    pout() << "KDNearestNeighbor returned an error" << endl;

  RealVect closestRV;
  for (int i=0; i<SpaceDim; i++)
    closestRV[i] = closestArr[i];
  if (m_printdebug)
  {
    pout() << "In KDTree: searched for " << node0;
    pout() << " found " << closestRV;
    pout() << " with IntVect " << inout0->first << " and boolean " << inout0->second << endl;
  }

  // do search for node1
  KDError = KDNearestNeighbor( m_ptree , node1Arr , closestArr , (void**) &inout1 , NULL , 0 );
  if (KDError!=0)
    pout() << "KDNearestNeighbor returned an error" << endl;

  for (int i=0; i<SpaceDim; i++)
    closestRV[i] = closestArr[i];
  if (m_printdebug)
  {
    pout() << "In KDTree: searched for " << node1;
    pout() << " found " << closestRV;
    pout() << " with IntVect " << inout1->first << " and boolean " << inout1->second << endl;
  }

  isNode0Inside = inout0->second;
  isNode1Inside = inout1->second;

}

// return true if point lies within the triangle (up to the tolerance mesh.tol??)
bool STLExplorer::IsPointInTriangle(const RealVect& point,
                                    const int& triangle)
{
  if (SpaceDim==3)
  {

    // stolen from blackpawn.com
    // compute vectors
    RealVect v0 = m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][1] ] - \
                  m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][0] ];
    RealVect v1 = m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][2] ] - \
                  m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][0] ];
    RealVect v2 = point - \
                  m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][0] ];

    // compute dot products
    Real dot00 = v0.dotProduct(v0);
    Real dot01 = v0.dotProduct(v1);
    Real dot02 = v0.dotProduct(v2);
    Real dot11 = v1.dotProduct(v1);
    Real dot12 = v1.dotProduct(v2);

    // compute barycentric coords
    Real u = (dot11 * dot02 - dot01 * dot12) / (dot00 * dot11 - dot01 * dot01);
    Real v = (dot00 * dot12 - dot01 * dot02) / (dot00 * dot11 - dot01 * dot01);

    // poor conditioning if denominator~0, i.e. v0~v1 so skinny isocelese triangle
    //if ( Abs(dot00 * dot11 - dot01 * dot01) < 1.0e-5 )
    //{
    //  pout() << "STLExplorer::IsPointInTriangle: Warning, poorly conditioned for triangle " << triangle << "\n";
    //}

    // make sure the point is in the plane of the triangle as well
    RealVect projection = m_msh->vertices.vertex[ m_msh->triangles.corners[triangle][0] ] + u*v0 + v*v1;
    Real dnormal = (point-projection).vectorLength();

    Real newtol = sqrt(m_msh->tol);

    return (u >= -newtol) && (v >= -newtol) && (u + v < 1+newtol) && (dnormal < newtol);
  }
  else
  {
    MayDay::Abort("Point in triangle not implemented for anything but 3D");
    return 0;
  }
}

// return true if the point lies on the edge (up to the tolerance mesh.tol)
bool STLExplorer::IsPointOnCellEdge(const RealVect& point,
                                    const CellEdge& celledge)
{
  /*
  // for each direction, must be between node0 and node1 of edge
  bool condition = true;
  RealVect node0 = IVToRV(celledge.m_node0, m_sb->m_origin, m_sb->m_dx);
  RealVect node1 = IVToRV(celledge.m_node1, m_sb->m_origin, m_sb->m_dx);
  RealVect cent = (node0 + node1)/2;
  node0 = (node1 - node0);
  Real length = node0.vectorLength();
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    condition = condition && (Abs(point[idir]-cent[idir]) < (length/2 + m_msh->tol));
  }
  return condition;
  */

  // another try, should be easier
  bool condition = true;
  RealVect node0 = IVToRV(celledge.m_node0, m_sb->m_origin, m_sb->m_dx);
  RealVect node1 = IVToRV(celledge.m_node1, m_sb->m_origin, m_sb->m_dx);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    if (idir == celledge.m_dir)
    {
      // between endpoints, assuming node0[idir]<node1[idir], which it should be for a CellEdge
      condition = condition && ( point[idir]>(node0[idir]-m_msh->tol) && point[idir]<(node1[idir]+m_msh->tol) ); 
    }
    else
    {
      // within a small radius around node0/node1 (should be the same in this direction)
      condition = condition && ( Abs(point[idir]-node0[idir]) < m_msh->tol );
      if (Abs(node0[idir]-node1[idir])>m_msh->tol)
        pout() << "STLExplorer::IsPointOnCellEdge: Bad input, CellEdge has non-adjacent nodes" << endl;
    }
  }

  return condition;
}
 
// note, i is an input and is incremented in FillInCellLine
void STLExplorer::FillInCellLine(vector<IntVect>& cells,
                                 int&             i,
                                 const int&       itri,
                                 const int&       idir0,
                                 const int&       idir1)
{
  // fill in a gap between cells[i] and cells[i-1]
  // and make sure to increment i afterwards to point to the original i'th cell 
  // idir0 is the direction for the line jump (idir0 tells you which line you are on)
  // idir1 is the direction along a single line
  IntVect c0 = cells[i-1];
  IntVect c1 = cells[i];
  int origsize = cells.size();

  // well just do something stupid, take all cells in the rectangle formed by c0 and c1
  // should probably make a line between c0 and c1 and take all of the cells
  // along that line, e.g. first the cells taken by going right/up then by going up/right
  int lo0 = c0[idir0]; // in direction 0, we know that c0 < c1
  int hi0 = c1[idir0];
  int lo1 = Min(c0[idir1],c1[idir1]); // in direction 1, not sure which is larger/smaller
  int hi1 = Max(c0[idir1],c1[idir1]);

  int ninrect = (hi0-lo0+1) * (hi1-lo1+1); // number of cells in the rectangle (number we will add plus 2)
  
  IntVect tmpcell; // cell to add. start off with bottom left (lo) corner
  tmpcell[idir0] = lo0; 
  tmpcell[idir1] = lo1;
  if (SpaceDim==3)
  {
    tmpcell[(idir1+1)%SpaceDim] = c0[(idir1+1)%SpaceDim];
    // quick sanity check
    if (c0[(idir1+1)%SpaceDim] != c1[(idir1+1)%SpaceDim] || idir1 != (idir0+1)%SpaceDim)
      pout() << "STLExplorer::FillInCellLine: Warning, inconsistent inputs." << endl;
  }
  
  int ninserted = 0;
  for (int j = lo0; j < hi0+1; j++)
  {
    for (int k = lo1; k < hi1+1; k++)
    {
      tmpcell[idir0] = j;
      tmpcell[idir1] = k;
      if (tmpcell != c0 && tmpcell != c1)
      {
        cells.insert( cells.begin()+i+ninserted , tmpcell ); // insert new cells  

        // add to cell map as well
        TriInCell tmptic; tmptic.triangles.push_back(itri);
        pair<CellMapIt, bool> rvCell = m_sb->m_cellmap.insert( make_pair(tmpcell,tmptic) );
        if (!rvCell.second)
          rvCell.first->second.triangles.push_back(itri);

        ninserted++;
      }
    }
  }
  i+=ninserted; // increment i

  // sanity check
  if (cells.size()-origsize != ninrect-2 || ninserted != ninrect-2)
  {
    pout() << "STLExplorer::FillInCellLine: Warning, wrong number of cells to add." << endl;
  }
}

// Note, this is only a temporary data holder
// and may contain invalid data (IntVects outside m_box)
void STLExplorer::GetVertMap(Vector<IntVect>** a_vertmap)
{
  *a_vertmap = &m_vertmap;
}

// Note, this is only a temporary data holder
// and may contain invalid data (IntVects outside m_box)
void STLExplorer::GetTriMap(Vector<Vector<IntVect> >** a_trimap)
{
  *a_trimap = &m_trimap;
}

void STLExplorer::GetSTLBox(RefCountedPtr<STLBox>& a_sb)
{
  a_sb = m_sb;
}

#include "NamespaceFooter.H"
