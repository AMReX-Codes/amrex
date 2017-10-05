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

    Vector<Box> m_domains;

    Vector<BoxArray> m_covered_ba;
    Vector<BoxArray> m_irregular_ba;

    Vector<FabArray<EBCellFlagFab>> m_cellflags;
    Vector<MultiFab> m_volfrac;
    Vector<MultiCutFab> m_bndrycent;
    Vector<std::array<MultiCutFab,AMREX_SPACEDIM>> m_areafrac;
    Vector<std::array<MultiCutFab,AMREX_SPACEDIM>> m_facecent;
};

}

#endif
