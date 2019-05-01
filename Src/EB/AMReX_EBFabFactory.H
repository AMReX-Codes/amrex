#ifndef AMREX_EBFABFACTORY_H_
#define AMREX_EBFABFACTORY_H_

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
    virtual ~EBFArrayBoxFactory () = default;

    EBFArrayBoxFactory (const EBFArrayBoxFactory&) = default;
    EBFArrayBoxFactory (EBFArrayBoxFactory&&) noexcept = default;

    EBFArrayBoxFactory () = delete;
    EBFArrayBoxFactory& operator= (const EBFArrayBoxFactory&) = delete;
    EBFArrayBoxFactory& operator= (EBFArrayBoxFactory&&) = delete;

    virtual FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const final;

    virtual void destroy (FArrayBox* fab) const final;

#ifdef AMREX_USE_GPU
    virtual FArrayBox* createDeviceAlias (FArrayBox const& src) const final;
    virtual void destroyDeviceAlias (FArrayBox* fab) const final;
#endif

    virtual EBFArrayBoxFactory* clone () const final;

    const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const noexcept
        { return m_ebdc->getMultiEBCellFlagFab(); }

    const MultiFab& getVolFrac () const noexcept { return m_ebdc->getVolFrac(); }

    const MultiCutFab& getCentroid () const noexcept { return m_ebdc->getCentroid(); }

    const MultiCutFab& getBndryCent () const noexcept { return m_ebdc->getBndryCent(); }

    const MultiCutFab& getBndryNormal () const noexcept { return m_ebdc->getBndryNormal(); }

    const MultiCutFab& getBndryArea () const noexcept { return m_ebdc->getBndryArea(); }

    Array<const MultiCutFab*,AMREX_SPACEDIM> getAreaFrac () const noexcept {
        return m_ebdc->getAreaFrac();
    }

    Array<const MultiCutFab*,AMREX_SPACEDIM> getFaceCent () const noexcept {
        return m_ebdc->getFaceCent();
    }

    EB2::Level const* getEBLevel () const noexcept { return m_parent; }
    EB2::IndexSpace const* getEBIndexSpace () const noexcept;
    int maxCoarseningLevel () const noexcept;

    const DistributionMapping& DistributionMap () const noexcept;
    const BoxArray& boxArray () const noexcept;

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
