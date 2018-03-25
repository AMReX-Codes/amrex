
#ifndef ____MOMENTITERATORImplem_H_
#define ____MOMENTITERATORImplem_H_

#include "AMReX_IndexedMoments.H"


template<int D, int P>
MomentIterator<D, P>::
MomentIterator()
{
  m_indicies = &(IndexedMoments<D,P>::getMonomialPowers());
  m_current = 0;
}


#endif
