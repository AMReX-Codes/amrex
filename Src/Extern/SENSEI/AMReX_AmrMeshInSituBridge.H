#ifndef AMReX_AmrMeshInSituBridge_H
#define AMReX_AmrMeshInSituBridge_H
#include <AMReX_Config.H>

#include <AMReX_InSituBridge.H>
#include <AMReX_Vector.H>

namespace amrex
{
class AmrMesh;
class MultiFab;

/// SENSEI bridge for code simulations that use amrex::AmrMesh/Core
class AmrMeshInSituBridge : public InSituBridge
{
public:
    AmrMeshInSituBridge() {}
    ~AmrMeshInSituBridge() {}

    AmrMeshInSituBridge(const AmrMeshInSituBridge&) = delete;
    void operator=(const AmrMeshInSituBridge&) = delete;

    /// invoke the in situ analysis with data from an AmrMesh instance.
    int update(unsigned int step, double time, amrex::AmrMesh *mesh,
        const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
        const std::vector<std::vector<std::string>> &names);
};

}

#endif
