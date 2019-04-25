
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

struct CnsFillExtDir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
        }
};

namespace {
    static CnsFillExtDir cns_fill_ext_dir;
    static GpuBndryFuncFab<CnsFillExtDir> gpu_bndry_func(cns_fill_ext_dir);
    static CpuBndryFuncFab cpu_bndry_func(nullptr); // Without EXT_DIR (e.g., inflow), we can pass a nullptr
}

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
#ifdef AMREX_USE_GPU
    bool run_on_gpu = Gpu::inLaunchRegion();
#else
    bool run_on_gpu = false;
#endif

    if (run_on_gpu) {
        AMREX_ASSERT(Gpu::isGpuPtr(&data));
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
    } else {
        cpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
    }
}
