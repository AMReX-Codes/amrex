#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "STLMesh.H"
#include <math.h>
using namespace std;

#include "STLUtil.H"

#include "NamespaceHeader.H"

using namespace STLUtil;

// print info about the mesh
void STLMesh::PrintMesh()
{

  pout() << "STL mesh with " << triangles.corners.size() << " elements\n";

  pout() << "Vertices (total " << vertices.vertex.size() << ")";
  PVec(vertices.vertex); pout() << "\n";

  pout() << "Triangles (total " << triangles.corners.size() << ")";
  PVec(triangles.corners);pout() << "\n";

  pout() << "Edges (total " << edges.edge.size() << ")";
  PVec(edges.edge);pout() << "\n";

  pout() << "Edge to triangle";
  PVec(connect.edgeToTriangle);pout() << "\n";

  pout() << "Vertex to triangle";
  PVec(connect.vertexToTriangle);pout() << "\n";

}

void STLMesh::Transform(const Real     scale,
                        const RealVect translate,
                        const Real     theta,
                        RealVect       axis)
{

  if (SpaceDim>3)
    MayDay::Abort("STL Mesh: cannot do transformations in more than 3D");

  if (SpaceDim==3 && Abs(axis.vectorLength()-1.0) < tol)
    axis /= axis.vectorLength(); // normalize

  // build rotation matrix
  Real c = cos(theta);
  Real s = sin(theta);
  Vector<RealVect> rot(SpaceDim);
  if (SpaceDim==3)
  {
    // see wikipedia article on rotation matrix, angle-and-axis formula
    // column 1
    rot[0][0] = axis[0]*axis[0]*(1-c)+c;
    rot[0][1] = axis[1]*axis[0]*(1-c)+axis[2]*s;
    rot[0][2] = axis[2]*axis[0]*(1-c)-axis[1]*s;
    // column 2
    rot[1][0] = axis[0]*axis[1]*(1-c)-axis[2]*s;
    rot[1][1] = axis[1]*axis[1]*(1-c)+c;
    rot[1][2] = axis[2]*axis[1]*(1-c)+axis[0]*s;
    // column 3
    rot[2][0] = axis[0]*axis[2]*(1-c)+axis[1]*s;
    rot[2][1] = axis[1]*axis[2]*(1-c)-axis[0]*s;
    rot[2][2] = axis[2]*axis[2]*(1-c)+c;
  }
  else if (SpaceDim==2)
  {
    // ignore axis - assume z-axis, rotate by theta ccw
    rot[0][0] = c;
    rot[0][1] = s;

    rot[1][0] = -s;
    rot[1][1] = c;
  }

  // do actual transformation of vertices
  RealVect tmpv;
  for (int i=0; i<vertices.vertex.size(); i++)
  {
    vertices.vertex[i] = scale*vertices.vertex[i] + translate;
    tmpv = RealVect::Zero;
    for (int idir=0; idir<SpaceDim; idir++)
      tmpv += vertices.vertex[i][idir]*rot[idir]; // column-by-column matrix mult.
    vertices.vertex[i] = tmpv;
  }

  // transform normal vectors
  for (int i=0; i<triangles.normal.size(); i++)
  {
    // just rotate the normal vectors
    tmpv = RealVect::Zero;
    for (int idir=0; idir<SpaceDim; idir++)
      tmpv += triangles.normal[i][idir]*rot[idir]; // column-by-column matrix mult.
    triangles.normal[i] = tmpv;
  }

}

#include "NamespaceFooter.H"

