#ifndef AMREX_FPHYSBC_H_
#define AMREX_FPHYSBC_H_

#include <AMReX_PhysBCFunct.H>

namespace amrex
{

class FPhysBC
    : public PhysBCFunctBase
{
public:
    FPhysBC () {}
    virtual ~FPhysBC () {}
    virtual void FillBoundary (MultiFab& mf, int scomp, int ncomp, Real time) override;
    using amrex::PhysBCFunctBase::FillBoundary;

    using fill_physbc_funptr_t = void(*)(MultiFab* mf, int scomp, int ncomp, Real time);
    fill_physbc_funptr_t fill_physbc = nullptr;
};

}
#endif
