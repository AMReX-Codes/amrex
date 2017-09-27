#ifndef AMREX_EBTOWER_H_
#define AMREX_EBTOWER_H_

#include <AMReX_FabArray.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

class EBTower
{
public:

    static void Build ();
    static void Destroy ();

private:

    EBTower ();
    ~EBTower ();

    static EBTower* m_instance;

    Array<Box> m_domains;

    Array<BoxArray> m_covered_ba;
    Array<BoxArray> m_irregular_ba;

    Array<FabArray<EBCellFlagFab> > m_cellflags;
};

}

#endif
