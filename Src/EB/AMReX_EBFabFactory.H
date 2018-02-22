#ifndef AMREX_EBFABFACTORY_H_
#define AMREX_EBFABFACTORY_H_

#include <AMReX_FabFactory.H>

#include <AMReX_EBDataCollection.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBSupport.H>

namespace amrex
{

class EBFArrayBoxFactory
    : public FabFactory<FArrayBox>
{
public:

    EBFArrayBoxFactory (const Geometry& a_geom, const BoxArray& a_ba, const DistributionMapping& a_dm,
                        const Vector<int>& a_ngrow, EBSupport a_support);
    virtual ~EBFArrayBoxFactory () = default;

    EBFArrayBoxFactory (const EBFArrayBoxFactory&) = default;
    EBFArrayBoxFactory (EBFArrayBoxFactory&&) noexcept = default;

    EBFArrayBoxFactory () = delete;
    EBFArrayBoxFactory& operator= (const EBFArrayBoxFactory&) = delete;
    EBFArrayBoxFactory& operator= (EBFArrayBoxFactory&&) = delete;

    virtual FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override;

    virtual EBFArrayBoxFactory* clone () const override;

    virtual Box getDomain () const override { return m_geom.Domain(); }

    const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const
        { return m_ebdc->getMultiEBCellFlagFab(); }

    const MultiFab& getVolFrac () const { return m_ebdc->getVolFrac(); }

//    const MultiCutFab& getCentroid () const { return m_ebdc->getCentroid(); }

    const MultiCutFab& getBndryCent () const { return m_ebdc->getBndryCent(); }

    std::array<const MultiCutFab*,AMREX_SPACEDIM> getAreaFrac () const {
        return m_ebdc->getAreaFrac();
    }

    std::array<const MultiCutFab*,AMREX_SPACEDIM> getFaceCent () const {
        return m_ebdc->getFaceCent();
    }

private:

    EBSupport m_support;
    Geometry m_geom;
    std::shared_ptr<EBDataCollection> m_ebdc;

};

}

#endif

