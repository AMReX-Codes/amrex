#ifndef AMREX_EBFABFACTORY_H_
#define AMREX_EBFABFACTORY_H_

#include <AMReX_FabFactory.H>

#include <AMReX_EBLevel.H>
#include <AMReX_EBDataCollection.H>
#include <AMReX_Geometry.H>

namespace amrex
{
class EBFArrayBoxFactory
    : public FabFactory<FArrayBox>
{
public:

    EBFArrayBoxFactory (const Geometry& a_geom, const EBLevel& a_eblevel);
    virtual ~EBFArrayBoxFactory () = default;

    EBFArrayBoxFactory (const EBFArrayBoxFactory&) = default;
    EBFArrayBoxFactory (EBFArrayBoxFactory&&) noexcept = default;

    EBFArrayBoxFactory () = delete;
    EBFArrayBoxFactory& operator= (const EBFArrayBoxFactory&) = delete;
    EBFArrayBoxFactory& operator= (EBFArrayBoxFactory&&) = delete;

    const EBLevel& getEBLevel () const { return m_eblevel; }
    const EBISLayout& getEBISL () const { return m_eblevel.getEBISL(); }

    virtual FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override;

    virtual EBFArrayBoxFactory* clone () const override;

private:

    Geometry m_geom;
    std::shared_ptr<EBDataCollection> m_ebdc;
    EBLevel m_eblevel;
};

}

#endif

