#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

struct NullFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& /*iv*/, Array4<Real> const& /*dest*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     GeometryData const& /*geom*/, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
            // no physical boundaries to fill because it is all periodic
        }
};

void nullfill (Box const& bx, FArrayBox& data,
               const int dcomp, const int numcomp,
               Geometry const& geom, const Real time,
               const Vector<BCRec>& bcr, const int bcomp,
               const int scomp)
{
    GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill{});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}
