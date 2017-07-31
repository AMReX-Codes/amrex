#ifndef AmrAdvPhysBC_H_
#define AmrAdvPhysBC_H_

#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. And FillBoundary does nothing.

class AmrCoreAdvPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    AmrCoreAdvPhysBC () {}
    virtual ~AmrCoreAdvPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};

#endif
