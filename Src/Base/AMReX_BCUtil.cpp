#include <AMReX_BCUtil.H>
#include <AMReX_PhysBCFunct.H>

namespace amrex
{

namespace {

void dummy_cpu_fill_extdir (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    // do something for external Dirichlet (BCType::ext_dir) if there are
}

struct dummy_gpu_fill_extdir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir) if there are
        }
};

}

void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc)
{
    if (geom.isAllPeriodic()) return;
    if (phi.nGrow() == 0) return;

    AMREX_ALWAYS_ASSERT(phi.ixType().cellCentered());

    if (Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<dummy_gpu_fill_extdir> gpu_bndry_func(dummy_gpu_fill_extdir{});
        PhysBCFunct<GpuBndryFuncFab<dummy_gpu_fill_extdir> > physbcf
            (geom, bc, gpu_bndry_func);
        physbcf.FillBoundary(phi, 0, phi.nComp(), 0.0, 0);
    }
    else
    {
        CpuBndryFuncFab cpu_bndry_func(dummy_cpu_fill_extdir);;
        PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
        physbcf.FillBoundary(phi, 0, phi.nComp(), 0.0, 0);
    }
}

}
