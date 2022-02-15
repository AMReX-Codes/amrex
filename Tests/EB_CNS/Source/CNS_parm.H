#ifndef CNS_PARM_H_
#define CNS_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct Parm
{
    amrex::Real eos_gamma = 1.4;
    amrex::Real eos_mu = 28.97;  // mean molecular weight

    amrex::Real cv;
    amrex::Real cp;

    amrex::Real Pr  = 0.72;     // Prandtl number
    amrex::Real C_S = 1.458e-5; // constant in Sutherland's law
    amrex::Real T_S = 110.4;    // Sutherland temperature

    amrex::Real smallr = 1.e-19;
    amrex::Real smallp = 1.e-10;

    amrex::Real const_visc_mu = -1.0;
    amrex::Real const_visc_ki = -1.0;
    amrex::Real const_lambda  = -1.0;

    static constexpr int level_mask_interior   = 0; // valid cells
    static constexpr int level_mask_covered    = 1; // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2; // ghost cells not covered
    static constexpr int level_mask_physbnd    = 3; // outside domain

    void Initialize ();
};

#endif
