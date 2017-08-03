#ifndef AMREX_EBFABFACTORY_H_
#define AMREX_EBFABFACTORY_H_

#include <AMReX_FabFactory.H>

#include <AMReX_EBLevel.H>

namespace amrex
{
class EBFArrayBoxFactory
    : public FabFactory<FArrayBox>
{
public:

    EBFArrayBoxFactory (const EBLevel& a_eblevel);
    virtual ~EBFArrayBoxFactory () = default;

    EBFArrayBoxFactory (const EBFArrayBoxFactory&) = default;
    EBFArrayBoxFactory (EBFArrayBoxFactory&&) noexcept = default;

    EBFArrayBoxFactory () = delete;
    EBFArrayBoxFactory& operator= (const EBFArrayBoxFactory&) = delete;
    EBFArrayBoxFactory& operator= (EBFArrayBoxFactory&&) = delete;

    const EBLevel& getEBLevel () const { return m_eblevel; }
    const EBISLayout& getEBISLayout () const { return m_ebisl; }

    virtual FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override;

    virtual EBFArrayBoxFactory* clone () const override;

private:
    EBLevel m_eblevel;
    EBISLayout m_ebisl;
};

}

#endif

