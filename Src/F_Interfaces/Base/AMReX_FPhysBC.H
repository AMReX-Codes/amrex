#ifndef AMREX_FPHYSBC_H_
#define AMREX_FPHYSBC_H_

#include <AMReX_PhysBCFunct.H>

namespace amrex
{

class FPhysBC
    : public PhysBCFunctBase
{
public:
    using fill_physbc_funptr_t = void(*)(MultiFab* mf, int scomp, int ncomp, Real time);

    FPhysBC (fill_physbc_funptr_t fill) : fill_physbc(fill) {}

    virtual ~FPhysBC () {}
    virtual void FillBoundary (MultiFab& mf, int scomp, int ncomp, Real time) override;
    using amrex::PhysBCFunctBase::FillBoundary;

private:
    fill_physbc_funptr_t fill_physbc = nullptr;
};

}
#endif
