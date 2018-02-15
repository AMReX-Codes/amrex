#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STLMESH_H_
#define _STLMESH_H_

#include "RealVect.H"

#include "NamespaceHeader.H"

/*
 * This data structure class holds all of the data
 * pertaining to the STL surface mesh, including the
 * connectivity data generated while reading.  It also
 * has functions to print and perform affine transformations
 * on the mesh.
 * The data here is all public and filled in by an STLReader
 */

/// The whole mesh
class STLMesh
{
public:

  /// Vertices
  struct vertlist {

    // master list of vertices
    Vector<RealVect> vertex;

  };

  /// Edges
  struct edgelist {

    // list of edges (indices into list of vertices)
    // each edge has 2 vertices
    Vector<Vector<int> > edge;

  };

  /// Triangles
  struct trilist {

    // list of triangles
    // each triangle has 3 vertices
    Vector<Vector<int> > corners;
    // and a normal vector
    Vector<RealVect> normal;

  };

  /// Connectivity
  struct conninfo {

    // triangles connected to each vertex (indices)
    // vertexToTriangle.length()==vertex.length() (first index into vertices, then triangles)
    Vector<Vector<int> > vertexToTriangle;
    // each edge has a left and right triangle
    // (left as you go from vertex[0] to vertex[1]
    // edgeToTriangle.length()==edge.length() (first index into edges, then triangles)
    Vector<Vector<int> > edgeToTriangle;

  };

  // vertices
  vertlist vertices;
  // edges
  edgelist edges;
  // triangles
  trilist triangles;
  // connectivity
  conninfo connect;

  // tolerance for collapsing nodes and such
  Real tol;

  void PrintMesh(); // print info about the mesh

  // transform stuff
  void Transform(const Real     scale,
                 const RealVect translate,
                 const Real     theta,
                 RealVect       axis);

};

#include "NamespaceFooter.H"
#endif

