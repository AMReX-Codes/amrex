#ifndef AmrAdvBC_H_
#define AmrAdvBC_H_

#include <PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. And FillBoundary does nothing.

class AmrAdvPhysBC
    : public PhysBCFunctBase
{
public:
    AmrAdvPhysBC () {}
    virtual ~AmrAdvPhysBC () {}
    virtual void FillBoundary (MultiFab& mf, int, int, Real time) override { }
    using PhysBCFunctBase::FillBoundary;
};

#endif
