#ifndef AMREX_EBFABFACTORY_H_
#define AMREX_EBFABFACTORY_H_
#include <AMReX_Config.H>

#include <AMReX_FabFactory.H>

#include <AMReX_EBData.H>
#include <AMReX_EBDataCollection.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBSupport.H>
#include <AMReX_Array.H>
#include <AMReX_MFIter.H>

/**
 * \file AMReX_EBFabFactory.H
 *
 * FabFactory that produces FArrayBoxes aware of embedded-boundary metadata.
 */

namespace amrex
{

namespace EB2 {
    class Level;
    class IndexSpace;
}

/// \ingroup amrex_eb
//! FabFactory specialization that exposes EB geometry alongside state FABs.
class EBFArrayBoxFactory
    : public FabFactory<FArrayBox>
{
public:

    /**
     * \brief Build an EB-aware factory backed by EB2 level data.
     *
     * \param a_level   EB level to sample.
     * \param a_geom    Geometry describing the level.
     * \param a_ba      Grid layout.
     * \param a_dm      Distribution mapping for FAB ownership.
     * \param a_ngrow   Grow-cell requirements (basic/volume/full).
     * \param a_support Requested EB support level.
     */
    EBFArrayBoxFactory (const EB2::Level& a_level, const Geometry& a_geom,
                        const BoxArray& a_ba, const DistributionMapping& a_dm,
                        const Vector<int>& a_ngrow, EBSupport a_support);

    /**
     * \brief Build an EB factory for face-centered grids.
     *
     * \param a_level   EB level to sample.
     * \param a_geom    Geometry describing the level.
     * \param a_ba      Grid layout.
     * \param a_dm      Distribution mapping for FAB ownership.
     * \param a_ngrow   Grow-cell requirements (basic/volume/full).
     * \param a_support Requested EB support level.
     * \param face_dir  Face direction (0=x, 1=y, 2=z).
     */
    EBFArrayBoxFactory (const EB2::Level& a_level, const Geometry& a_geom,
                        const BoxArray& a_ba, const DistributionMapping& a_dm,
                        const Vector<int>& a_ngrow, EBSupport a_support, int face_dir);
    ~EBFArrayBoxFactory () override = default;

    EBFArrayBoxFactory (const EBFArrayBoxFactory&) = default;
    EBFArrayBoxFactory (EBFArrayBoxFactory&&) noexcept = default;

    EBFArrayBoxFactory () = delete;
    EBFArrayBoxFactory& operator= (const EBFArrayBoxFactory&) = delete;
    EBFArrayBoxFactory& operator= (EBFArrayBoxFactory&&) = delete;

    //! Allocate an EB-aware FAB for `box` with `ncomps` components.
    AMREX_NODISCARD
    FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const final;

    //! Create an alias view into an existing EB-aware FAB.
    AMREX_NODISCARD
    FArrayBox* create_alias (FArrayBox const& rhs, int scomp, int ncomp) const final;

    //! Destroy a FAB produced by this factory.
    void destroy (FArrayBox* fab) const final;

    //! Produce a deep copy of this factory (excluding EB level ownership).
    AMREX_NODISCARD
    EBFArrayBoxFactory* clone () const final;

    //! EB cell flags for all boxes.
    [[nodiscard]] const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const noexcept
        { return m_ebdc->getMultiEBCellFlagFab(); }

    //! Level set field.
    [[nodiscard]] const MultiFab& getLevelSet () const noexcept { return m_ebdc->getLevelSet(); }

    //! Volume fraction.
    [[nodiscard]] const MultiFab& getVolFrac () const noexcept { return m_ebdc->getVolFrac(); }

    //! Cell centroid.
    [[nodiscard]] const MultiCutFab& getCentroid () const noexcept { return m_ebdc->getCentroid(); }

    //! Boundary centroid.
    [[nodiscard]] const MultiCutFab& getBndryCent () const noexcept { return m_ebdc->getBndryCent(); }

    //! Boundary normals.
    [[nodiscard]] const MultiCutFab& getBndryNormal () const noexcept { return m_ebdc->getBndryNormal(); }

    //! Boundary area.
    [[nodiscard]] const MultiCutFab& getBndryArea () const noexcept { return m_ebdc->getBndryArea(); }

    //! Area fractions along each axis.
    [[nodiscard]] Array<const MultiCutFab*,AMREX_SPACEDIM> getAreaFrac () const noexcept {
        return m_ebdc->getAreaFrac();
    }

    //! Face centroids along each axis. For y-face-centroid, the two
    //! components are for x and then z-directions, respectively.
    [[nodiscard]] Array<const MultiCutFab*,AMREX_SPACEDIM> getFaceCent () const noexcept {
        return m_ebdc->getFaceCent();
    }

    //! Edge centroids.
    [[nodiscard]] Array<const MultiCutFab*,AMREX_SPACEDIM> getEdgeCent () const noexcept {
        return m_ebdc->getEdgeCent();
    }

    //! \return True if the EB is regular everywhere.
    [[nodiscard]] bool isAllRegular () const noexcept;

    //! \return EB level pointer.
    [[nodiscard]] EB2::Level const* getEBLevel () const noexcept { return m_parent; }
    //! \return EB IndexSpace pointer.
    [[nodiscard]] EB2::IndexSpace const* getEBIndexSpace () const noexcept;
    //! Maximum coarse level supported by the EB hierarchy.
    [[nodiscard]] int maxCoarseningLevel () const noexcept;

    //! Distribution mapping describing the FABs this factory allocates.
    [[nodiscard]] const DistributionMapping& DistributionMap () const noexcept;
    //! BoxArray describing the FAB layout.
    [[nodiscard]] const BoxArray& boxArray () const noexcept;
    //! Geometry inherited from the EB level.
    [[nodiscard]] const Geometry& Geom () const noexcept { return m_geom; }
    //! Face direction if this is an FC factory (-1 for CC mode).
    [[nodiscard]] int faceDir () const noexcept { return m_face_dir; }

    //! True if EB geometric data beyond flags is available.
    [[nodiscard]] bool hasEBInfo() const noexcept;

    //! Returns nullptr unless this level is built by EB2::addRegularCoarseLevels.
    //! One should use getMultiEBCellFlagFab for normal levels.
    [[nodiscard]] iMultiFab const* getCutCellMask () const noexcept { return m_ebdc->getCutCellMask(); }

    //! Convenience wrapper returning EBData views for the FAB referenced by \p mfi.
    [[nodiscard]] EBData getEBData (MFIter const& mfi) const noexcept;

    //! Return GPU-accessible array views for all FABs owned by this factory.
    [[nodiscard]] EBDataArrays getEBDataArrays () const noexcept;

private:

    EBSupport m_support;
    Geometry m_geom;
    std::shared_ptr<EBDataCollection> m_ebdc;
    EB2::Level const* m_parent = nullptr;
    Gpu::DeviceVector<Array4<Real const>> m_eb_data;
    int m_face_dir = -1;  //!< Face direction (0=x, 1=y, 2=z) or -1 for cell-centered
};

/// \ingroup amrex_eb
//! Convenience overload that builds the EB level internally from Geometry + IndexSpace::top().
std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const Geometry& a_geom,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support);

/// \ingroup amrex_eb
//! Convenience overload that reuses an already-built EB level.
std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const EB2::Level*,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support);

/// \ingroup amrex_eb
//! Convenience overload that builds a level from a provided IndexSpace.
std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const EB2::IndexSpace*, const Geometry& a_geom,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support);

}

#endif
