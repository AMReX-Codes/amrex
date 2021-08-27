#include <AMReX_REAL.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex::Real* /*problo*/,
                         const amrex::Real* /*probhi*/)
    {
        // Nothing needs to be done here,
        // since there are no extra inputs to be read from probin file
    }
}
