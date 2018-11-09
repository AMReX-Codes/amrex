
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the desciptor set up in CNS::variableSetUp.

void cns_bcfill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
#if AMREX_USE_GPU
    bool run_on_gpu = Gpu::inLaunchRegion();
#else
    bool run_on_gpu = false;
#endif

    if (run_on_gpu) {
        // Without EXT_DIR (e.g., inflow), we can pass a nullptr
        GpuBndryFuncFab(nullptr)(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
    } else {
        CpuBndryFuncFab(nullptr)(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
    }
}
