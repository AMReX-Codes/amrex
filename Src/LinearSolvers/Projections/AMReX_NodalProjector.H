#ifndef NODAL_PROJECTION_H
#define NODAL_PROJECTION_H

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>

namespace amrex {

class NodalProjector
{

public:


    NodalProjector () {};

    NodalProjector ( const Vector<Geometry>&                a_geom,
                     const Vector<BoxArray>&                a_grids,
                     const Vector<DistributionMapping>&     a_dmap,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi,
                     const LPInfo&                          a_lpinfo );

    NodalProjector ( const Vector<Geometry>&                a_geom,
                     const Vector<BoxArray>&                a_grids,
                     const Vector<DistributionMapping>&     a_dmap,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi )
        : NodalProjector(a_geom, a_grids, a_dmap, a_bc_lo, a_bc_hi, LPInfo()) {};

#ifdef AMREX_USE_EB
    NodalProjector ( const Vector<Geometry>&                a_geom,
                     const Vector<BoxArray>&                a_grids,
                     const Vector<DistributionMapping>&     a_dmap,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi,
                     Vector<EBFArrayBoxFactory const *>     a_ebfactory,
                     const LPInfo&                          a_lpinfo );

    NodalProjector ( const Vector<Geometry>&                a_geom,
                     const Vector<BoxArray>&                a_grids,
                     const Vector<DistributionMapping>&     a_dmap,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                     std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi,
                     Vector<EBFArrayBoxFactory const *>     a_ebfactory)
        : NodalProjector(a_geom, a_grids, a_dmap, a_bc_lo, a_bc_hi, a_ebfactory, LPInfo()) {};
#endif


    void define ( const Vector<Geometry>&                a_geom,
                  const Vector<BoxArray>&                a_grids,
                  const Vector<DistributionMapping>&     a_dmap,
                  std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                  std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi,
                  const LPInfo&                          a_lpinfo);

#ifdef AMREX_USE_EB
    void define ( const Vector<Geometry>&                a_geom,
                  const Vector<BoxArray>&                a_grids,
                  const Vector<DistributionMapping>&     a_dmap,
                  std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                  std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi,
                  Vector<EBFArrayBoxFactory const *>     a_ebfactory,
                  const LPInfo&                          a_lpinfo );
#endif

    void project ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                   const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                   const amrex::Vector<amrex::MultiFab*>        a_S_cc = {},
                   const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    void project2 ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                    const amrex::Vector<const amrex::MultiFab*>& a_alpha,
                    const amrex::Vector<const amrex::MultiFab*>& a_beta,
                    const amrex::Vector<const amrex::MultiFab*>& a_rhs );

    Vector< const MultiFab* > getPhi ()     const {return GetVecOfConstPtrs(m_phi);};
    Vector< const MultiFab* > getGradPhi () const {return GetVecOfConstPtrs(m_fluxes);};

    void computeRHS ( const amrex::Vector<amrex::MultiFab*>&       a_rhs,
                      const amrex::Vector<amrex::MultiFab*>&       a_vel,
                      const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                      const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

private:

    bool m_ok = false;

    // Verbosity
    int  m_verbose = 0;
    int  m_mg_verbose = 0;
    int  m_mg_cg_verbose = 0;

    // Control MLMG behavior
    int  m_mg_maxiter = 100;
    int  m_mg_cg_maxiter = 100;
    Real m_mg_rtol = 1.0e-11;
    Real m_mg_atol = 1.0e-14;
    std::string m_bottom_solver_type = "bicgcg";


    Vector<Geometry>                    m_geom;
    Vector<BoxArray>                    m_grids;
    Vector<DistributionMapping>         m_dmap;
    LPInfo                              m_lpinfo;

#ifdef AMREX_USE_EB
    Vector<EBFArrayBoxFactory const *>  m_ebfactory;
#endif

    // Boundary conditions
    std::array<LinOpBCType,AMREX_SPACEDIM>  m_bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM>  m_bc_hi;

    // Cell-centered data
    Vector< std::unique_ptr< MultiFab > >  m_fluxes;

    // Node-centered data
    Vector< std::unique_ptr< MultiFab > >  m_phi;
    Vector< std::unique_ptr< MultiFab > >  m_rhs;

    // Linear operator
    std::unique_ptr< MLNodeLaplacian > m_matrix;

    // Solver
    std::unique_ptr< MLMG > m_solver;

    void readParameters ();

    void printInfo ();

    void setup ();

};

}

#endif
