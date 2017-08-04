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

    const FabArray<EBFlagFab>& Flags () const { return *m_flags; }

private:
    std::shared_ptr<FabArray<EBFlagFab> > m_flags;
};

}

#endif
