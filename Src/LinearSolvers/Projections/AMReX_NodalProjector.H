#ifndef NODAL_PROJECTION_H
#define NODAL_PROJECTION_H

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>

class NodalProjector
{

public:


    NodalProjector () {};

    NodalProjector ( const amrex::Vector<amrex::Geometry>&              a_geom,
                      const amrex::Vector<amrex::BoxArray>&             a_grids,
                      const amrex::Vector<amrex::DistributionMapping>&  a_dmap,
                      std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_lo,
                      std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_hi,
                      amrex::Vector<amrex::EBFArrayBoxFactory const *>  a_ebfactory );

    void define ( const amrex::Vector<amrex::Geometry>&             a_geom,
                  const amrex::Vector<amrex::BoxArray>&             a_grids,
                  const amrex::Vector<amrex::DistributionMapping>&  a_dmap,
                  std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_lo,
                  std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_hi,
                  amrex::Vector<amrex::EBFArrayBoxFactory const *>  a_ebfactory );

    void project ( amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                   const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_sigma,
                   const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_S_cc = {},
                   const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_S_nd = {} );

    void project2 (  amrex::Vector< std::unique_ptr< amrex::MultiFab > >&       a_vel,
                     const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_alpha,
                     const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_beta,
                     const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_rhs );

    amrex::Vector< const amrex::MultiFab* > getPhi ()     const {return GetVecOfConstPtrs(m_phi);};
    amrex::Vector< const amrex::MultiFab* > getGradPhi () const {return GetVecOfConstPtrs(m_fluxes);};

    void computeRHS (  amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_rhs,
                       const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                       const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_S_cc = {},
                       const amrex::Vector< std::unique_ptr< amrex::MultiFab > >& a_S_nd = {} );

private:

    bool m_ok = false;

    // Verbosity
    int  m_verbose = 0;
    int  m_mg_verbose = 0;
    int  m_mg_cg_verbose = 0;

    // Control MLMG behavior
    int  m_mg_maxiter = 100;
    int  m_mg_cg_maxiter = 100;
    amrex::Real m_mg_rtol = 1.0e-11;
    amrex::Real m_mg_atol = 1.0e-14;
    std::string m_bottom_solver_type = "bicgcg";

    // Max coarsening level
    int  m_mg_max_coarsening_level = 100;


    amrex::Vector<amrex::Geometry>                    m_geom;
    amrex::Vector<amrex::BoxArray>                    m_grids;
    amrex::Vector<amrex::DistributionMapping>         m_dmap;
    amrex::Vector<amrex::EBFArrayBoxFactory const *>  m_ebfactory;


    // Boundary conditions
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_lo;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_hi;

    // Cell-centered data
    amrex::Vector< std::unique_ptr< amrex::MultiFab > >  m_fluxes;

    // Node-centered data
    amrex::Vector< std::unique_ptr< amrex::MultiFab > >  m_phi;
    amrex::Vector< std::unique_ptr< amrex::MultiFab > >  m_rhs;

    // Linear operator
    std::unique_ptr< amrex::MLNodeLaplacian > m_matrix;

    // Solver
    std::unique_ptr< amrex::MLMG > m_solver;

    void readParameters ();

    void printInfo ();

    void setup ();

};


#endif
