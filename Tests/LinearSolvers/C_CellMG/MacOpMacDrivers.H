
#ifndef _MACOPMACDRIVERS_H_

//====================================================================
// This header file contains declarations for the MacOperator mac projection
// driver routines
//====================================================================


// the driver functions for the MAC solves
void mac_level_driver( amrex::AmrCore* parent, const amrex::MacBndry &mac_bndry, const amrex::BCRec& phys_bc,
                       const amrex::BoxArray &grids,
                       int use_cg_solve, int level, int Density,
                       const amrex::Real *dx, amrex::Real dt,
                       amrex::Real mac_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
                       amrex::MultiFab *area,  amrex::MultiFab &volume,
                       amrex::MultiFab &S,     amrex::MultiFab &Rhs,
                       amrex::MultiFab *u_mac, amrex::MultiFab *mac_phi );

// a driver function for computing a sync MAC solve
void mac_sync_driver( amrex::AmrCore* parent, const amrex::MacBndry &mac_bndry, const amrex::BCRec& phys_bc,
                      const amrex::BoxArray &grids,
                      int use_cg_solve, int level, 
                      const amrex::Real *dx, amrex::Real dt,
                      amrex::Real mac_sync_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
                      amrex::MultiFab *area,  amrex::MultiFab &volume,
                      amrex::MultiFab &Rhs,   amrex::MultiFab *rho_half,
                      amrex::MultiFab *mac_sync_phi );


// apply the mac pressure gradient to a velocity field
// init, means that velocities are initialized here
void mac_vel_update( int init,
                     amrex::FArrayBox &uvel,
                     amrex::FArrayBox &vvel,
#if (BL_SPACEDIM == 3 )
                     amrex::FArrayBox &wvel,
#endif
                     const amrex::FArrayBox &phi,
                     const amrex::FArrayBox *rho, int rho_comp,  
                     const amrex::Box &grd, int level, int n,
                     const amrex::Real *dx, amrex::Real scale );

#endif
