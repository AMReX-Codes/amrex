#ifndef NODAL_PROJECTION_H
#define NODAL_PROJECTION_H

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>

//
// Solves
//
//     div(sigma*grad(phi)) = div(vel) + S_nd + S_cc
//
// and  performs the projection
//
//     vel = vel - ( sigma / alpha ) * grad(phi)
//
// vel, sigma, alpha, and S_cc are cell-centered variables, while
// phi and S_nd are nodal-centered variables.
//
// By default alpha is assumed to be 1. Use setAlpha to change the default.
//
// The user can provide a custom RHS instead of letting NodalProjector to compute
// one. Use setCustomRHS to provide a custom RHS.
//
namespace amrex {

class NodalProjector
{

public:

    NodalProjector ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                     const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                     const Vector<Geometry>&                      a_geom,
                     const LPInfo&                                a_lpinfo,
                     const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                     const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    NodalProjector ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                     const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                     const Vector<Geometry>&                      a_geom,
                     const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                     const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} )
    : NodalProjector(a_vel, a_sigma, a_geom, LPInfo(), a_S_cc, a_S_nd ) {}


    void project ( amrex::Real a_rtol = 1.0e-11, amrex::Real a_atol = 1.0e-14 );
    void project ( const Vector<amrex::MultiFab*>& a_phi, amrex::Real a_rtol = 1.0e-11,
                   amrex::Real a_atol = 1.0e-14 );

    Vector< const MultiFab* > getGradPhi () const {return GetVecOfConstPtrs(m_fluxes);}
    Vector< const MultiFab* > getPhi     () const {return GetVecOfConstPtrs(m_phi);}

    void computeRHS ( const amrex::Vector<amrex::MultiFab*>&       a_rhs,
                      const amrex::Vector<amrex::MultiFab*>&       a_vel,
                      const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                      const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    void setAlpha     (const amrex::Vector<const amrex::MultiFab*> a_alpha)
        {m_alpha=a_alpha;m_has_alpha=true;}
    void setCustomRHS (const amrex::Vector<const amrex::MultiFab*> a_rhs);


    // Methods to set verbosity
    void setVerbose (int  v) noexcept { m_verbose = v; }


    // Set domain BC
    void setDomainBC ( std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                       std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi );

    // Methods to get underlying objects
    // Use these to modify properties of MLMG and linear operator
    MLNodeLaplacian& getLinOp () noexcept { return *m_linop; }
    MLMG&            getMLMG  () noexcept { return *m_mlmg;  }

    // Methods to set MF for sync
    void setSyncResidualFine (MultiFab* a_sync_resid_fine) {m_sync_resid_fine=a_sync_resid_fine;}
    void setSyncResidualCrse (MultiFab* a_sync_resid_crse, IntVect a_ref_ratio, BoxArray a_fine_grids )
        {m_sync_resid_crse=a_sync_resid_crse; m_ref_ratio=a_ref_ratio; m_fine_grids=a_fine_grids;}

private:

    void setOptions ();
    void setCoarseBoundaryVelocityForSync ();
    void computeSyncResidual ();
    void averageDown (const amrex::Vector<amrex::MultiFab*> a_var);

    bool m_has_rhs   = false;
    bool m_has_alpha = false;
    bool m_need_bcs  = true;

    // Verbosity
    int  m_verbose        = 0;

    // Geometry
    Vector<Geometry>               m_geom;

    // EB factory if any
#ifdef AMREX_USE_EB
    Vector<EBFArrayBoxFactory const *>  m_ebfactory;
#endif

    // Cell-centered data
    Vector<amrex::MultiFab*>        m_vel;
    Vector<amrex::MultiFab>         m_fluxes;
    Vector<const amrex::MultiFab*>  m_alpha;
    Vector<amrex::MultiFab*>        m_S_cc;
    Vector<const amrex::MultiFab*>  m_sigma;

    // Node-centered data
    Vector<amrex::MultiFab>         m_phi;
    Vector<amrex::MultiFab>         m_rhs;
    Vector<const amrex::MultiFab*>  m_S_nd;

    // Linear operator
    std::unique_ptr< MLNodeLaplacian > m_linop;

    // Solver
    std::unique_ptr< MLMG > m_mlmg;

     // Boundary conditions
    std::array<LinOpBCType,AMREX_SPACEDIM>  m_bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM>  m_bc_hi;

    // Members for syncronization
    IntVect   m_ref_ratio;
    BoxArray  m_fine_grids;  // Grid at crse level + 1
    MultiFab* m_sync_resid_crse = nullptr;
    MultiFab* m_sync_resid_fine = nullptr;

    void printInfo ();
};

}

#endif
