
#ifndef AMREX_IFDATAIMPLEM_H_
#define AMREX_IFDATAIMPLEM_H_


#include <iostream>
#include <iomanip>

//#include "AMReX_MultiIndex.H"
//#include "AMReX_NormalDerivative.H"
#include "AMReX_NormalDerivativeNew.H"
#include "AMReX_MomentIterator.H"
//#include "AMReX_ParmParse.H"
//#include "AMReX.H"




template <int dim> void IFData<dim>::setNormalDerivatives()

template <int dim> void IFData<dim>::makeCornerSigns()

template <int dim> void IFData<dim>::findIntersectionPts()
template <int dim> void IFData<dim>::defineLocalCoords()

template <int dim> bool IFData<dim>::isConnected(int          & a_edgeDir,
                                                 const Vertex & a_vertex1,
                                                 const Vertex & a_vertex2)

template <int dim> void IFData<dim>::makeEdgeKey(const int    & a_edgeDir,
                                                 const Vertex & a_vertex0,
                                                 const Vertex & a_vertex1)

template <int dim> Real IFData<dim>::rootFinder(const EdgeIndex& a_thisEdge)

template <int dim> Real IFData<dim>::BrentRootFinder(const RvDim & a_loPt,
                                                     const RvDim & a_hiPt,
                                                     const int   & a_edgeDir) const
//                                                     const Real  & a_smallestRoot,
//                                                     const Real  & a_biggestRoot) const

template <int dim> void IFData<dim>::checkIntersection(bool       & a_hiOn,
                                                       bool       & a_loOn,
                                                       const Real & a_pt) const

template <int dim> void IFData<dim>::remakeCornerSigns()


// equals operator
template <int dim> void IFData<dim>::operator=(const IFData & a_IFData)


#endif
