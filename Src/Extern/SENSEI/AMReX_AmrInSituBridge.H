#ifndef AMReX_AmrInSituBridge_H
#define AMReX_AmrInSituBridge_H
#include <AMReX_Config.H>

#include <AMReX_InSituBridge.H>

namespace amrex
{

class Amr;

/// Contains the bridge code for simulations that use amrex::Amr
class AmrInSituBridge : public InSituBridge
{
public:
    AmrInSituBridge() {}
    ~AmrInSituBridge() {}

    AmrInSituBridge(const AmrInSituBridge&) = delete;
    void operator=(const AmrInSituBridge&) = delete;

    /// invoke the in situ analysis with data from an amrex::Amr instance
    int update(amrex::Amr *dataSource);
};

}

#endif
