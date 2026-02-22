#ifndef AMREX_EB_DATA_COLLECTION_H_
#define AMREX_EB_DATA_COLLECTION_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBSupport.H>
#include <AMReX_Array.H>

/**
 * \file AMReX_EBDataCollection.H
 *
 * Manages EB geometry data (flags, centroids, etc.) extracted from EB2::Level.
 */

namespace amrex {

class EBFArrayBoxFactory;
template <class T> class FabArray;
class MultiFab;
class iMultiFab;
class MultiCutFab;
namespace EB2 { class Level; }

//! Aggregates the EB data required to build EB-aware MultiFabs/FabFactories.
class EBDataCollection
{
public:

    friend class EBFArrayBoxFactory;

    /**
     * \brief Allocate EB data components with the requested support level.
     *
     * \param a_level   EB2 level providing data.
     * \param a_geom    Geometry used for filling.
     * \param a_ba      Grid layout.
     * \param a_dm      Distribution mapping for data ownership.
     * \param a_ngrow   Grow cells per component (basic/volume/full).
     * \param a_support Requested EBSupport tier (basic/volume/full).
     */
    EBDataCollection (const EB2::Level& a_level, const Geometry& a_geom,
                      const BoxArray& a_ba, const DistributionMapping& a_dm,
                      Vector<int>  a_ngrow, EBSupport a_support);

    ~EBDataCollection ();

    EBDataCollection (const EBDataCollection&) = delete;
    EBDataCollection (EBDataCollection&&) = delete;
    EBDataCollection& operator= (const EBDataCollection&) = delete;
    EBDataCollection& operator= (EBDataCollection&&) = delete;

    //! Access the EBCellFlag data (always available at least for EBSupport::basic).
    [[nodiscard]] const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const;
    //! Signed-distance field of the EB (basic support and above).
    [[nodiscard]] const MultiFab& getLevelSet () const;
    //! Volume fraction (EBSupport::volume/full).
    [[nodiscard]] const MultiFab& getVolFrac () const;
    //! Cell-centroid MultiCutFab.
    [[nodiscard]] const MultiCutFab& getCentroid () const;
    //! Boundary centroid MultiCutFab.
    [[nodiscard]] const MultiCutFab& getBndryCent () const;
    //! Boundary area MultiCutFab.
    [[nodiscard]] const MultiCutFab& getBndryArea () const;
    //! Boundary normal MultiCutFab.
    [[nodiscard]] const MultiCutFab& getBndryNormal () const;
    //! Area fractions for faces in each coordinate direction.
    [[nodiscard]] Array<const MultiCutFab*, AMREX_SPACEDIM> getAreaFrac () const;
    //! Face centroid MultiCutFabs.
    [[nodiscard]] Array<const MultiCutFab*, AMREX_SPACEDIM> getFaceCent () const;
    //! Edge centroid MultiCutFabs.
    [[nodiscard]] Array<const MultiCutFab*, AMREX_SPACEDIM> getEdgeCent () const;
    //! Optional cut-cell mask (only for addRegularCoarseLevels).
    [[nodiscard]] const iMultiFab* getCutCellMask () const;

    // public for cuda
    //! Extend ghost data beyond the domain to satisfy EB support requirements.
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
