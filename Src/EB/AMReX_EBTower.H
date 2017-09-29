#ifndef AMREX_EBTOWER_H_
#define AMREX_EBTOWER_H_

#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiFab.H>

namespace amrex {

class EBLevelGrid;

class EBTower
{
public:

    static void Build ();
    static void Destroy ();
    static const EBTower* get () { return m_instance; }

    static void fillEBCellFlag (FabArray<EBCellFlagFab>& a_flag, const Geometry& a_geom);
    static void fillVolFrac (MultiFab& a_volfrac, const Geometry& a_geom);

private:

    EBTower ();
    ~EBTower ();

    void initCellFlags (int lev, const EBLevelGrid& eblg);
    void initVolFrac (int lev, const EBLevelGrid& eblg);

    int getIndex (const Box& domain) const;

    static EBTower* m_instance;

    Array<Box> m_domains;

    Array<BoxArray> m_covered_ba;
    Array<BoxArray> m_irregular_ba;

    Array<FabArray<EBCellFlagFab> > m_cellflags;
    Array<MultiFab> m_volfrac;
};

}

#endif
