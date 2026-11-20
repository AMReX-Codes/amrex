#include "WMLES.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

struct WmlesFillExtDir
{
    Real* inflow_state = nullptr;

    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real /*time*/,
                     const BCRec* bcr, const int bcomp,
                     const int /*orig_comp*/) const
        {
            const Box& domain_box = geom.Domain();

            const BCRec& bc = bcr[bcomp+0];

            AMREX_D_TERM(int i = iv[0];,
                         int j = iv[1];,
                         int k = iv[2];);

            // Inflow at low-x
            if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
            {
                for (int nc = 0; nc < numcomp; ++nc)
                   dest(i,j,k,dcomp+nc) = inflow_state[dcomp+nc];
            }
        }
};

void wmles_bcfill (Box const& bx, FArrayBox& data,
                  const int dcomp, const int numcomp,
                  Geometry const& geom, const Real time,
                  const Vector<BCRec>& bcr, const int bcomp,
                  const int scomp)
{
    GpuBndryFuncFab<WmlesFillExtDir> gpu_bndry_func(WmlesFillExtDir{WMLES::h_prob_parm->inflow_state});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}
