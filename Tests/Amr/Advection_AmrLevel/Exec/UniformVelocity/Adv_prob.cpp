#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuContainers.H>

#include "AmrLevelAdv.H"

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex::Real* /*problo*/,
                         const amrex::Real* /*probhi*/)
    {
        // Read the prob block from the input file using ParmParse
        amrex::ParmParse pp("prob");
        amrex::Vector<amrex::Real> adv_vel(AMREX_SPACEDIM);
        pp.queryarr("adv_vel", adv_vel);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            AmrLevelAdv::h_prob_parm->adv_vel[i] = adv_vel[i];
        }

        // Transfer the problem-specific data to the GPU
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         AmrLevelAdv::h_prob_parm, AmrLevelAdv::h_prob_parm+1,
                         AmrLevelAdv::d_prob_parm);
    }
}
