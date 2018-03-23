
#ifndef _MULTIINDEX_H_
#define _MULTIINDEX_H_


#include <map>
using std::map;

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_IndexTM.H"
#include "AMReX_Notation.H"
#include "AMReX_IFData.H"

// Generate and return a vector of all multi-indices whose magnitude (the sum
// of their entries) equals "a_magnitude"
template <int dim> void generateMultiIndices(Vector<IndexTM<int,dim> > &      a_indices,
                                             const int                      & a_magnitude);


#include "AMReX_MultiIndexImplem.H"

#endif
