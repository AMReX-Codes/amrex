#ifndef AMREX_EBTOWER_H_
#define AMREX_EBTOWER_H_

#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>

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
    static void fillBndryCent (MultiCutFab& a_bndrycent, const Geometry& a_geom);
    static void fillFaceGeometry (std::array<MultiCutFab*,AMREX_SPACEDIM>& a_areafrac,
                                  std::array<MultiCutFab*,AMREX_SPACEDIM>& a_facecent,
                                  const Geometry& a_geom);

private:

    EBTower ();
    ~EBTower ();

    void initCellFlags (int lev, const EBLevelGrid& eblg);
    void initVolFrac (int lev, const EBLevelGrid& eblg);
    void initBndryCent (int lev, const EBLevelGrid& eblg);
    void initFaceGeometry (int lev, const EBLevelGrid& eblg);

    int getIndex (const Box& domain) const;

    static EBTower* m_instance;

    Array<Box> m_domains;

    Array<BoxArray> m_covered_ba;
    Array<BoxArray> m_irregular_ba;

    Array<FabArray<EBCellFlagFab>> m_cellflags;
    Array<MultiFab> m_volfrac;
    Array<MultiCutFab> m_bndrycent;
    Array<std::array<MultiCutFab,AMREX_SPACEDIM>> m_areafrac;
    Array<std::array<MultiCutFab,AMREX_SPACEDIM>> m_facecent;
};

}

#endif
