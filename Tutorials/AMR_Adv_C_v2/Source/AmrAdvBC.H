#ifndef AmrAdvBC_H_
#define AmrAdvBC_H_

#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. And FillBoundary does nothing.

class AmrAdvPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    AmrAdvPhysBC () {}
    virtual ~AmrAdvPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};

#endif
