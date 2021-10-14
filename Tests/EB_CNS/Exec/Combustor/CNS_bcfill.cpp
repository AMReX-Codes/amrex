#include "CNS.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

struct CnsFillExtDir
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

            // The combustor problem is hard-wired for inflow at low-z right now
            if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
            {
                for (int nc = 0; nc < numcomp; ++nc)
                   dest(i,j,k,dcomp+nc) = inflow_state[dcomp+nc];
            }
        }
};

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the descriptor set up in CNS::variableSetUp.

void cns_bcfill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
    GpuBndryFuncFab<CnsFillExtDir> gpu_bndry_func(CnsFillExtDir{CNS::h_prob_parm->inflow_state});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}
