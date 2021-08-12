#ifndef Prob_H_
#define Prob_H_

#include <AMReX_REAL.H>

#include "Adv_F.H"
#include "Prob_Parm.H"

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        // Read the prob block from the input file using ParmParse
        amrex::ParmParse pp("prob");
        pp.query("adv_vel", h_prob_parm->adv_vel);

        // Transfer the problem-specific data to the GPU
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_prob_parm, h_prob_parm,
                         d_prob_parm);
    }
}
#endif
