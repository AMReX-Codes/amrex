#ifndef AMREX_EB_DATA_COLLECTION_H_
#define AMREX_EB_DATA_COLLECTION_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBSupport.H>
#include <AMReX_Array.H>

namespace amrex {

template <class T> class FabArray;
class MultiFab;
class iMultiFab;
class MultiCutFab;
namespace EB2 { class Level; }

class EBDataCollection
{
public:

    EBDataCollection (const EB2::Level& a_level, const Geometry& a_geom,
                      const BoxArray& a_ba, const DistributionMapping& a_dm,
                      Vector<int>  a_ngrow, EBSupport a_support);

    ~EBDataCollection ();

    EBDataCollection (const EBDataCollection&) = delete;
    EBDataCollection (EBDataCollection&&) = delete;
    EBDataCollection& operator= (const EBDataCollection&) = delete;
    EBDataCollection& operator= (EBDataCollection&&) = delete;

    [[nodiscard]] const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const;
    [[nodiscard]] const MultiFab& getLevelSet () const;
    [[nodiscard]] const MultiFab& getVolFrac () const;
    [[nodiscard]] const MultiCutFab& getCentroid () const;
    [[nodiscard]] const MultiCutFab& getBndryCent () const;
    [[nodiscard]] const MultiCutFab& getBndryArea () const;
    [[nodiscard]] const MultiCutFab& getBndryNormal () const;
    [[nodiscard]] Array<const MultiCutFab*, AMREX_SPACEDIM> getAreaFrac () const;
    [[nodiscard]] Array<const MultiCutFab*, AMREX_SPACEDIM> getFaceCent () const;
    [[nodiscard]] Array<const MultiCutFab*, AMREX_SPACEDIM> getEdgeCent () const;
    [[nodiscard]] const iMultiFab* getCutCellMask () const;

    // public for cuda
    void extendDataOutsideDomain (IntVect const& level_ng);

private:

    Vector<int> m_ngrow;
    EBSupport m_support;
    Geometry m_geom;

    // have to use pointer to break include loop

    // EBSupport::basic
    FabArray<EBCellFlagFab>* m_cellflags = nullptr;
    MultiFab* m_levelset = nullptr;

    // EBSupport::volume
    MultiFab* m_volfrac = nullptr;
    MultiCutFab* m_centroid = nullptr;

    // EBSupport::full
    MultiCutFab* m_bndrycent = nullptr;
    MultiCutFab* m_bndryarea = nullptr;
    MultiCutFab* m_bndrynorm = nullptr;
    Array<MultiCutFab*,AMREX_SPACEDIM> m_areafrac {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    Array<MultiCutFab*,AMREX_SPACEDIM> m_facecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    Array<MultiCutFab*,AMREX_SPACEDIM> m_edgecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};

    // for levels created by addRegularCoarseLevels only
    iMultiFab* m_cutcellmask = nullptr;
};

}

#endif
