#ifndef AMReX_AmrInSituBridge_H
#define AMReX_AmrInSituBridge_H

#include <AMReX_InSituBridge.H>

namespace amrex
{

class Amr;

/// Bridge for simulations that use amrex::Amr
/// inherits intialize and finalize from amrex::InSituBridge
class AmrInSituBridge : public InSituBridge
{
public:
    AmrInSituBridge() {}
    ~AmrInSituBridge() {}

    AmrInSituBridge(const AmrInSituBridge&) = delete;
    void operator=(const AmrInSituBridge&) = delete;

    // invoke the in situ analysis with data from an Amr
    // instance
    int update(amrex::Amr *dataSource);
};

}

#endif
