#ifndef AMREX_EBFABFACTORY_H_
#define AMREX_EBFABFACTORY_H_
#include <AMReX_Config.H>

#include <AMReX_FabFactory.H>

#include <AMReX_EBDataCollection.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBSupport.H>
#include <AMReX_Array.H>

namespace amrex
{

namespace EB2 {
    class Level;
    class IndexSpace;
}

class EBFArrayBoxFactory
    : public FabFactory<FArrayBox>
{
public:

    EBFArrayBoxFactory (const EB2::Level& a_level, const Geometry& a_geom,
                        const BoxArray& a_ba, const DistributionMapping& a_dm,
                        const Vector<int>& a_ngrow, EBSupport a_support);
    ~EBFArrayBoxFactory () override = default;

    EBFArrayBoxFactory (const EBFArrayBoxFactory&) = default;
    EBFArrayBoxFactory (EBFArrayBoxFactory&&) noexcept = default;

    EBFArrayBoxFactory () = delete;
    EBFArrayBoxFactory& operator= (const EBFArrayBoxFactory&) = delete;
    EBFArrayBoxFactory& operator= (EBFArrayBoxFactory&&) = delete;

    AMREX_NODISCARD
    FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const final;

    AMREX_NODISCARD
    FArrayBox* create_alias (FArrayBox const& rhs, int scomp, int ncomp) const final;

    void destroy (FArrayBox* fab) const final;

    AMREX_NODISCARD
    EBFArrayBoxFactory* clone () const final;

    [[nodiscard]] const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const noexcept
        { return m_ebdc->getMultiEBCellFlagFab(); }

    [[nodiscard]] const MultiFab& getLevelSet () const noexcept { return m_ebdc->getLevelSet(); }

    [[nodiscard]] const MultiFab& getVolFrac () const noexcept { return m_ebdc->getVolFrac(); }

    [[nodiscard]] const MultiCutFab& getCentroid () const noexcept { return m_ebdc->getCentroid(); }

    [[nodiscard]] const MultiCutFab& getBndryCent () const noexcept { return m_ebdc->getBndryCent(); }

    [[nodiscard]] const MultiCutFab& getBndryNormal () const noexcept { return m_ebdc->getBndryNormal(); }

    [[nodiscard]] const MultiCutFab& getBndryArea () const noexcept { return m_ebdc->getBndryArea(); }

    [[nodiscard]] Array<const MultiCutFab*,AMREX_SPACEDIM> getAreaFrac () const noexcept {
        return m_ebdc->getAreaFrac();
    }

    // For y-face-centroid, the two components are for x and then z-directions, respectively.
    [[nodiscard]] Array<const MultiCutFab*,AMREX_SPACEDIM> getFaceCent () const noexcept {
        return m_ebdc->getFaceCent();
    }

    [[nodiscard]] Array<const MultiCutFab*,AMREX_SPACEDIM> getEdgeCent () const noexcept {
        return m_ebdc->getEdgeCent();
    }

    [[nodiscard]] bool isAllRegular () const noexcept;

    [[nodiscard]] EB2::Level const* getEBLevel () const noexcept { return m_parent; }
    [[nodiscard]] EB2::IndexSpace const* getEBIndexSpace () const noexcept;
    [[nodiscard]] int maxCoarseningLevel () const noexcept;

    [[nodiscard]] const DistributionMapping& DistributionMap () const noexcept;
    [[nodiscard]] const BoxArray& boxArray () const noexcept;
    [[nodiscard]] const Geometry& Geom () const noexcept { return m_geom; }

private:

    EBSupport m_support;
    Geometry m_geom;
    std::shared_ptr<EBDataCollection> m_ebdc;
    EB2::Level const* m_parent = nullptr;
};

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const Geometry& a_geom,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support);

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const EB2::Level*,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support);

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const EB2::IndexSpace*, const Geometry& a_geom,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support);

}

#endif
