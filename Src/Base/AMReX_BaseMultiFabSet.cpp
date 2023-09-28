#include<AMReX_BaseMultiFabSet.H>

namespace amrex {

template<class MF>
void
BaseMultiFabSet<MF>::initializeSet(int setSize) {
    m_nSet = setSize;
    // m_mfptr_array = new MF*[setSize];
    // m_ixtype_set = new IndexType[setSize];
    // m_ngrow_set = new IntVect[setSize];

    m_mfptr_set(setSize);
    m_ixtype_set(setSize);
    m_ngrow_set(setSize);
}


}