#ifndef _Gravity_F_H_
#define _Gravity_F_H_
#include <AMReX_BLFort.H>

BL_FORT_PROC_DECL(FORT_AVGDOWN_PHI, fort_avgdown_phi)
    (BL_FORT_FAB_ARG(crse_fab),
     const BL_FORT_FAB_ARG(fine_fab),
     const int ovlo[], const int ovhi[], const int rat[]);

BL_FORT_PROC_DECL(FORT_EDGE_INTERP, fort_edge_interp)
    (const int* flo, const int* fhi,
     const int* nc, const int* refRatio, const int* dir,
     BL_FORT_FAB_ARG(fine));

BL_FORT_PROC_DECL(FORT_PC_EDGE_INTERP, fort_pc_edge_interp)
    (const int* lo, const int* hi,
     const int* nc, const int* refRatio, const int* dir,
     const BL_FORT_FAB_ARG(crse),
     BL_FORT_FAB_ARG(fine));

BL_FORT_PROC_DECL(FORT_GET_FLUXES, fort_get_fluxes)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(phi_cc),
     const BL_FORT_FAB_ARG(xgrad),
     const BL_FORT_FAB_ARG(ygrad),
     const BL_FORT_FAB_ARG(zgrad),
     const amrex::Real* dx);

BL_FORT_PROC_DECL(FORT_SET_HOMOG_BCS, fort_set_homog_bcs)
    (const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     BL_FORT_FAB_ARG(phi),
     const amrex::Real* dx);

BL_FORT_PROC_DECL(FORT_ADD_MONOPOLE_BCS, fort_add_monopole_bcs)
    (const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const int* nparticles,
     const amrex::Real* part_locs, const amrex::Real* part_mass,
     BL_FORT_FAB_ARG(phi), const amrex::Real* dx);

#ifdef CGRAV
BL_FORT_PROC_DECL(FORT_PRESCRIBE_GRAV,fort_prescribe_grav)
    (const int lo[], const int hi[],
     const amrex::Real* dx,
     const BL_FORT_FAB_ARG(S),
     const amrex::Real* problo,
     const int* addToSelfGrav);
#endif
#endif

