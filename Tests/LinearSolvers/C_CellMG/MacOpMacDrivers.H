
#ifndef _MACOPMACDRIVERS_H_

//====================================================================
// This header file contains declarations for the MacOperator mac projection
// driver routines
//====================================================================


// the driver functions for the MAC solves
void mac_level_driver( Amr* parent, const MacBndry &mac_bndry, const BCRec& phys_bc,
                       const BoxArray &grids,
                       int use_cg_solve, int level, int Density,
                       const Real *dx, Real dt,
                       Real mac_tol, Real mac_abs_tol, Real rhs_scale,
                       MultiFab *area,  MultiFab &volume,
                       MultiFab &S,     MultiFab &Rhs,
                       MultiFab *u_mac, MultiFab *mac_phi );

// a driver function for computing a sync MAC solve
void mac_sync_driver( Amr* parent, const MacBndry &mac_bndry, const BCRec& phys_bc,
                      const BoxArray &grids,
                      int use_cg_solve, int level, 
                      const Real *dx, Real dt,
                      Real mac_sync_tol, Real mac_abs_tol, Real rhs_scale,
                      MultiFab *area,  MultiFab &volume,
                      MultiFab &Rhs,   MultiFab *rho_half,
                      MultiFab *mac_sync_phi );


// apply the mac pressure gradient to a velocity field
// init, means that velocities are initialized here
void mac_vel_update( int init,
                     FArrayBox &uvel,
                     FArrayBox &vvel,
#if (BL_SPACEDIM == 3 )
                     FArrayBox &wvel,
#endif
                     const FArrayBox &phi,
                     const FArrayBox *rho, int rho_comp,  
                     const Box &grd, int level, int n,
                     const Real *dx, Real scale );

#endif
