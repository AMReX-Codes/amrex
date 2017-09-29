#ifndef AMREX_EB_DATA_COLLECTION_H_
#define AMREX_EB_DATA_COLLECTION_H_

#include <AMReX_Geometry.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBSupport.H>

namespace amrex {

template <class T> class FabArray;
class MultiFab;

class EBDataCollection
{
public:

    EBDataCollection (const Geometry& a_geom, const BoxArray& a_ba, const DistributionMapping& a_dm,
                      int a_ngrow, EBSupport a_support);

    ~EBDataCollection ();

    EBDataCollection (const EBDataCollection&) = delete;
    EBDataCollection (EBDataCollection&&) = delete;
    EBDataCollection& operator= (const EBDataCollection&) = delete;
    EBDataCollection& operator= (EBDataCollection&&) = delete;

    const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const;

private:

    int m_ngrow;
    EBSupport m_support;
    Geometry m_geom;

    // have to use pointer to break include loop
    FabArray<EBCellFlagFab>* m_cellflags;
};

}

#endif
