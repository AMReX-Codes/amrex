#ifndef AMREX_EB_DATA_COLLECTION_H_
#define AMREX_EB_DATA_COLLECTION_H_

#include <AMReX_Geometry.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBSupport.H>
#include <AMReX_Array.H>

namespace amrex {

template <class T> class FabArray;
class MultiFab;
class MultiCutFab;
namespace EB2 { class Level; }

class EBDataCollection
{
public:

    EBDataCollection (const EB2::Level& a_level, const Geometry& a_geom,
                      const BoxArray& a_ba, const DistributionMapping& a_dm,
                      const Vector<int>& a_ngrow, EBSupport a_support);

    ~EBDataCollection ();

    EBDataCollection (const EBDataCollection&) = delete;
    EBDataCollection (EBDataCollection&&) = delete;
    EBDataCollection& operator= (const EBDataCollection&) = delete;
    EBDataCollection& operator= (EBDataCollection&&) = delete;

    const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const;
    const MultiFab& getVolFrac () const;
    const MultiCutFab& getCentroid () const;
    const MultiCutFab& getBndryCent () const;
    const MultiCutFab& getBndryArea () const;
    const MultiCutFab& getBndryNormal () const;
    Array<const MultiCutFab*, AMREX_SPACEDIM> getAreaFrac () const;
    Array<const MultiCutFab*, AMREX_SPACEDIM> getFaceCent () const;

private:

    Vector<int> m_ngrow;
    EBSupport m_support;
    Geometry m_geom;

    // have to use pointer to break include loop

    // EBSupport::basic
    FabArray<EBCellFlagFab>* m_cellflags = nullptr;

    // EBSupport::volume
    MultiFab* m_volfrac = nullptr;
    MultiCutFab* m_centroid = nullptr;

    // EBSupport::full
    MultiCutFab* m_bndrycent = nullptr;
    MultiCutFab* m_bndryarea = nullptr;
    MultiCutFab* m_bndrynorm = nullptr;
    Array<MultiCutFab*,AMREX_SPACEDIM> m_areafrac {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    Array<MultiCutFab*,AMREX_SPACEDIM> m_facecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
};

}

#endif
