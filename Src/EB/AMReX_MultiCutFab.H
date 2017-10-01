#ifndef AMREX_MULTICUTFAB_H_
#define AMREX_MULTICUTFAB_H_

#include <AMReX_MultiFab.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

class MultiCutFab
{
public:

    MultiCutFab ();

    MultiCutFab (const BoxArray& ba, const DistributionMapping& dm,
                 int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags);

    ~MultiCutFab ();

    MultiCutFab (MultiCutFab&& rhs) noexcept = default;
    
    MultiCutFab (const MultiCutFab& rhs) = delete;
    MultiCutFab& operator= (const MultiCutFab& rhs) = delete;
    MultiCutFab& operator= (MultiCutFab&& rhs) = delete;

    void define (const BoxArray& ba, const DistributionMapping& dm,
                 int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags);

    const FArrayBox& operator[] (const MFIter& mfi) const;

private:

    MultiFab m_data;
    const FabArray<EBCellFlagFab>* m_cellflags;
};

}

#endif
