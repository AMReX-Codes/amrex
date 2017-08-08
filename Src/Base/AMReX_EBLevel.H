#ifndef AMREX_EBLEVEL_H
#define AMREX_EBLEVEL_H

#include <AMReX_EBLevelGrid.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

class EBLevel
    : public EBLevelGrid
{
public:
    EBLevel ();
    ~EBLevel ();

    EBLevel (const EBLevel& rhs) = default;

    EBLevel (EBLevel&& rhs) = delete;
    EBLevel& operator= (const EBLevel& rhs) = delete;
    EBLevel& operator= (EBLevel&& rhs) = delete;

    EBLevel (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng);

    const FabArray<EBCellFlagFab>& CellFlags () const { return *m_cellflags; }
    const EBISLayout& getEBISL               () const { return *m_ebisl; }

private:
    std::shared_ptr<FabArray<EBCellFlagFab> > m_cellflags;
    std::shared_ptr<EBISLayout>               m_ebisl;
};

class MultiFab;
Box getLevelDomain (const MultiFab& mf);

const EBLevel& getEBLevel (const MultiFab& mf);

const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab (const MultiFab& mf);

}

#endif
