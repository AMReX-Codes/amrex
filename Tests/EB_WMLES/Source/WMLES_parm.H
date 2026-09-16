#ifndef WMLES_PARM_H_
#define WMLES_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct Parm
{
    amrex::Real eos_gamma = 1.4;
    amrex::Real eos_mu = 28.97;  // mean molecular weight

    amrex::Real cv;
    amrex::Real cp;

    amrex::Real Pr  = 0.72;     // Prandtl number (molecular)
    amrex::Real C_S = 1.458e-5; // constant in Sutherland's law
    amrex::Real T_S = 110.4;    // Sutherland temperature

    amrex::Real smallr = 1.e-19;
    amrex::Real smallp = 1.e-10;

    amrex::Real const_visc_mu = -1.0;
    amrex::Real const_visc_ki = -1.0;
    amrex::Real const_lambda  = -1.0;

    // ---------- WMLES / SGS model parameters ----------
    amrex::Real Cs_smag = 0.18;       // Smagorinsky constant
    amrex::Real Pr_t    = 0.90;       // turbulent Prandtl number
    amrex::Real sgs_filter_width_ratio = 1.0; // filter width factor relative to grid

    // ---------- Wall model parameters ----------
    int  wall_model_enable       = 0;  // 0=off, 1=equilibrium, 2=ODE-based
    int  wall_model_sample_npts  = 3;  // number of sampling points along wall-normal
    amrex::Real wall_model_sample_h = 0.1; // sampling height factor (relative to dx)
    int  wall_model_update_interval = 1;  // update wall stress every N time steps
    amrex::Real wall_model_kappa = 0.41;  // von Karman constant
    amrex::Real wall_model_B     = 5.2;   // log-law intercept
    amrex::Real wall_model_yplus_target = 30.0; // target y+ for sampling

    static constexpr int level_mask_interior   = 0; // valid cells
    static constexpr int level_mask_covered    = 1; // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2; // ghost cells not covered
    static constexpr int level_mask_physbnd    = 3; // outside domain

    void Initialize ();
};

#endif
