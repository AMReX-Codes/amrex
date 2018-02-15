
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

#include "AMReX_STLAsciiReader.H"
#include <sstream>
#include "AMReX_parstream.H"
#include "AMReX_STLUtil.H"

namespace amrex
{
  using namespace STLUtil;

/*
 * Reads ascii STL files and generates a mesh
 * see http://www.cplusplus.com/doc/tutorial/files/
 */

/// Constructor - read from standard input
  STLAsciiReader::STLAsciiReader()
  {
    m_header = NULL;

    ReadData(cin,0);
  }

/// Constructor - read from file name
  STLAsciiReader::STLAsciiReader(const string& a_filename)
  {
    m_header = NULL;

    ifstream curFile;
    curFile.open(a_filename.c_str(),ios::in);
    if (!curFile.good() || !curFile.is_open())
    {
      amrex::Abort("STLAsciiReader - unable to open file");
    }

    ReadData(curFile,0);

    curFile.close();
  }

/// Destructor
  STLAsciiReader::~STLAsciiReader()
  {
    delete m_header;
  }

/// Return pointer to header string
  string* STLAsciiReader::GetHeader() const
  {
    return m_header;
  }


/// Return number of elements
  void STLAsciiReader::GetNtri(int& a_ntri) const
  {
    a_ntri = m_ntri;
  }

/// Return whether number of elements from header matches file
  void STLAsciiReader::GetNtriMatch(bool& a_ntriMatch) const
  {
    a_ntriMatch = m_ntriMatch;
  }

/// Return pointer to the mesh
  RefCountedPtr<STLMesh> STLAsciiReader::GetMesh() const
  {
    return m_stlmesh;
  }

  void STLAsciiReader::ReadData(istream&   a_file,
                                const int offset)
  {

    BL_PROFILE("STLAsciiReader::ReadData");

    RefCountedPtr<STLMesh> temp(new STLMesh());
    m_stlmesh = temp;
    m_header = new string;
    m_stlmesh->tol = 1.0e-10; // maybe do something fancier later, for now just constant

    string line; // temp variable for the whole line
    string word; // temp variable for text words (need to get rid of before extracting numbers)
    istringstream iss;

    // read first line "solid mySolid"
    getline(a_file,line);
    iss.str(line);
    iss >> word; // extract "solid"
    iss >> *m_header; // put the rest of the line in m_header


    // start big loop
    RealVect normal;
    Vector<RealVect> verts(3);
    int itri = 0;
    while (a_file.good())
    {

      // CH_START(tmesh);

      // each triangle has 7 lines
      getline(a_file,line); // "facet normal # # #"
      iss.str(line); iss.clear(); iss.seekg(0,ios::beg);
      iss >> word; // if this is a triangle, will be "facet"; if end of the file, will be "endsolid"

      if (!word.compare("endsolid"))
        break;

      iss >> word; // get rid of "facet" and "normal"

      iss >> normal[0] >> normal[1] >> normal[2]; // read normal vector

      getline(a_file,line); // "outer loop" (ignore)

      getline(a_file,line); // "vertex # # #"
      iss.str(line); iss.clear(); iss.seekg(0,ios::beg);
      iss >> word; // get rid of "vertex"
      iss >> verts[0][0] >> verts[0][1] >> verts[0][2]; // read vertex 0

      getline(a_file,line); // "vertex # # #"
      iss.str(line); iss.clear(); iss.seekg(0,ios::beg);
      iss >> word; // get rid of "vertex"
      iss >> verts[1][0] >> verts[1][1] >> verts[1][2]; // read vertex 1

      getline(a_file,line); // "vertex # # #"
      iss.str(line); iss.clear(); iss.seekg(0,ios::beg);
      iss >> word; // get rid of "vertex"
      iss >> verts[2][0] >> verts[2][1] >> verts[2][2]; // read vertex 2

      getline(a_file,line); // "endloop" (ignore)
      getline(a_file,line); // "endfacet" (ignore)

      // CH_STOP(tmesh);
      // CH_START(tdat);

      // now update mesh

      // check for degenerate triangles
      if ((verts[0]-verts[1]).vectorLength() < m_stlmesh->tol || \
          (verts[1]-verts[2]).vectorLength() < m_stlmesh->tol || \
          (verts[2]-verts[0]).vectorLength() < m_stlmesh->tol)
      {
        pout() << "STLAsciiReader: Building mesh: Warning, encountered degenerate triangle\n";
        pout() << " itri = " << itri << " corners = "; PRV(verts[0]); PRV(verts[1]); PRV(verts[2]); pout() << "\n";
      }

      // initialize triangle
      m_stlmesh->triangles.corners.resize(itri+1);
      m_stlmesh->triangles.corners[itri].resize(3);
      m_stlmesh->triangles.normal.resize(itri+1);

      // insert normal
      m_stlmesh->triangles.normal[itri] = normal;

      for (int i = 0; i < SpaceDim; i++)
        m_stlmesh->triangles.corners[itri][i]=-1; // initialize to -1

      // see if vertices exist already
      bool condition;
      for (int ivertg = 0; ivertg < m_stlmesh->vertices.vertex.size(); ivertg++)
      {
        for (int ivertl = 0; ivertl < SpaceDim; ivertl++)
        {
          // if all dimensions match within tol, then it's the same point
          condition = true;
          for (int j = 0; j < SpaceDim; j++)
            condition = condition && Abs(verts[ivertl][j]-m_stlmesh->vertices.vertex[ivertg][j])<m_stlmesh->tol;
          if (condition)
            m_stlmesh->triangles.corners[itri][ivertl] = ivertg;
        }
      }


      // if vertices don't exist, add them
      for (int ivertl = 0; ivertl < SpaceDim; ivertl++)
      {
        if (m_stlmesh->triangles.corners[itri][ivertl]==-1)
        {
          m_stlmesh->vertices.vertex.push_back(verts[ivertl]);
          m_stlmesh->triangles.corners[itri][ivertl] = m_stlmesh->vertices.vertex.size()-1;
          // and initialize a place in connect.vertexToTriangle
          m_stlmesh->connect.vertexToTriangle.resize( m_stlmesh->vertices.vertex.size() );
        }
      }

      // see if edges exist already, if not add them
      Vector<int> tmpedge(2);
      Vector<int> tmpedg2(2);
      for (int iedgel = 0; iedgel < SpaceDim; iedgel ++)
      {
        // create edges
        tmpedge[0] = m_stlmesh->triangles.corners[itri][ (iedgel)   % 3 ];
        tmpedge[1] = m_stlmesh->triangles.corners[itri][ (iedgel+1) % 3 ];
        tmpedg2[0]=tmpedge[1]; tmpedg2[1]=tmpedge[0]; // swap nodes
        bool foundedge = false; // flag to see if we need to add a new edge

        for (int iedgeg = 0; iedgeg < m_stlmesh->edges.edge.size(); iedgeg++)
        {
          // if all indices are the same, it's the same edge
          if ( m_stlmesh->edges.edge[iedgeg].stdVector() == tmpedge.stdVector() || \
               m_stlmesh->edges.edge[iedgeg].stdVector() == tmpedg2.stdVector() )
          {
            foundedge=true;
            if (m_stlmesh->connect.edgeToTriangle[iedgeg][0]==-1)
              m_stlmesh->connect.edgeToTriangle[iedgeg][0] = itri;
            else if (m_stlmesh->connect.edgeToTriangle[iedgeg][1]==-1)
              m_stlmesh->connect.edgeToTriangle[iedgeg][1] = itri;
            else
            {
              //MayDay::Abort("STLAsciiReader: Building mesh: edge has more than two triangles connected to it");
              pout() << "STLAsciiReader: Building mesh: edge has more than two triangles connected\n";
              printf(" iedgel=%i, iedgeg=%i, itri=%i, connected tri1=%i, connected tri2=%i\n",iedgel,iedgeg,itri,m_stlmesh->connect.edgeToTriangle[iedgeg][0],m_stlmesh->connect.edgeToTriangle[iedgeg][1]);
              //pout() << " old left="; PRV(m_stlmesh->vertices.vertex[ tmpedge[0] ]);
              //pout() << " right="; PRV(m_stlmesh->vertices.vertex[ tmpedge[1] ]); pout() << "\n";
              //pout() << " current left="; PRV(verts[ (iedgel) % 3 ]);
              //pout() << " right="; PRV(verts[ (iedgel+1) % 3 ]); pout() << "\n";
              pout() << " this triangle: "; PRV(verts[0]); PRV(verts[1]); PRV(verts[2]); pout() << "\n";
              int itri1 = m_stlmesh->connect.edgeToTriangle[iedgeg][0];
              pout() << " old triangle1: "; PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][0] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][1] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][2] ]); pout() << "\n";
              itri1 = m_stlmesh->connect.edgeToTriangle[iedgeg][1];
              pout() << " old triangle2: "; PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][0] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][1] ]);
              PRV(m_stlmesh->vertices.vertex[ m_stlmesh->triangles.corners[itri1][2] ]); pout() << "\n";
              //m_stlmesh->PrintMesh(); pout() << "\n";
            }
          }
        }

        // add edge to list and update connectivity
        if (!foundedge)
        {
          m_stlmesh->edges.edge.push_back(tmpedge);
          Vector<int> tmpe2t(2);
          tmpe2t[0]=itri; tmpe2t[1]=-1;
          m_stlmesh->connect.edgeToTriangle.push_back(tmpe2t);
        }
      }

      // add vertex to triangle connectivity
      m_stlmesh->connect.vertexToTriangle[ m_stlmesh->triangles.corners[itri][0] ].push_back(itri);
      m_stlmesh->connect.vertexToTriangle[ m_stlmesh->triangles.corners[itri][1] ].push_back(itri);
      m_stlmesh->connect.vertexToTriangle[ m_stlmesh->triangles.corners[itri][2] ].push_back(itri);


      itri++; // move to next triangle
    }

    m_ntri = itri; // set number of triangles
    m_ntriMatch = m_stlmesh->triangles.corners.size()==m_ntri;
    //m_stlmesh->PrintMesh();

  }
}
