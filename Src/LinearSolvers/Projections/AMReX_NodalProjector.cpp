#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_NodalProjector.H>
#include <AMReX_ParmParse.H>

namespace amrex {

NodalProjector::NodalProjector ( const amrex::Vector<amrex::Geometry>&               a_geom,
                                   const amrex::Vector<amrex::BoxArray>&             a_grids,
                                   const amrex::Vector<amrex::DistributionMapping>&  a_dmap,
                                   std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_lo,
                                   std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_hi )
{
    define(a_geom, a_grids, a_dmap, a_bc_lo, a_bc_hi);
}

#ifdef AMREX_USE_EB
NodalProjector::NodalProjector ( const amrex::Vector<amrex::Geometry>&               a_geom,
                                   const amrex::Vector<amrex::BoxArray>&             a_grids,
                                   const amrex::Vector<amrex::DistributionMapping>&  a_dmap,
                                   std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_lo,
                                   std::array<amrex::LinOpBCType,AMREX_SPACEDIM>     a_bc_hi,
                                   amrex::Vector<amrex::EBFArrayBoxFactory const *>  a_ebfactory )
{
    define(a_geom, a_grids, a_dmap, a_bc_lo, a_bc_hi, a_ebfactory);
}
#endif

void
NodalProjector::define ( const  amrex::Vector<amrex::Geometry>&              a_geom,
                         const  amrex::Vector<amrex::BoxArray>&              a_grids,
                         const  amrex::Vector<amrex::DistributionMapping>&   a_dmap,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM>       a_bc_lo,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM>       a_bc_hi)
{

    m_geom      = a_geom;
    m_grids     = a_grids;
    m_dmap      = a_dmap;
    m_bc_lo     = a_bc_lo;
    m_bc_hi     = a_bc_hi;

    int nlev( m_grids.size() );

    // Resize member data
    m_phi.resize(nlev);
    m_fluxes.resize(nlev);
    m_rhs.resize(nlev);

    // Allocate member data
    int ng(1);      // We use 1 ghost node only -- it should be enough

    for (int lev(0); lev < nlev; ++lev )
    {
        // Cell-centered data
        m_fluxes[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], AMREX_SPACEDIM, ng));

        // Node-centered data
	const auto& ba_nd = m_grids[lev].surroundingNodes();
        m_phi[lev].reset(new MultiFab(ba_nd, m_dmap[lev], 1, ng));
        m_rhs[lev].reset(new MultiFab(ba_nd, m_dmap[lev], 1, ng));
    }

    // Get inputs from ParmParse
    readParameters();

    // object is ready
    m_ok = true;

    // First setup
    setup();
}


#ifdef AMREX_USE_EB
void
NodalProjector::define ( const  amrex::Vector<amrex::Geometry>&                   a_geom,
                         const  amrex::Vector<amrex::BoxArray>&                   a_grids,
                         const  amrex::Vector<amrex::DistributionMapping>&        a_dmap,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM>            a_bc_lo,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM>            a_bc_hi,
                         amrex::Vector<amrex::EBFArrayBoxFactory const *>         a_ebfactory )
{

    m_geom      = a_geom;
    m_grids     = a_grids;
    m_dmap      = a_dmap;
    m_bc_lo     = a_bc_lo;
    m_bc_hi     = a_bc_hi;
    m_ebfactory = a_ebfactory;

    int nlev( m_grids.size() );

    // Resize member data
    m_phi.resize(nlev);
    m_fluxes.resize(nlev);
    m_rhs.resize(nlev);

    // Allocate member data
    int ng(1);      // We use 1 ghost node only -- it should be enough

    for (int lev(0); lev < nlev; ++lev )
    {
        // Cell-centered data
        m_fluxes[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], AMREX_SPACEDIM, ng, MFInfo(), *m_ebfactory[lev]));

        // Node-centered data
	const auto& ba_nd = m_grids[lev].surroundingNodes();
        m_phi[lev].reset(new MultiFab(ba_nd, m_dmap[lev], 1, ng, MFInfo(), *m_ebfactory[lev]));
        m_rhs[lev].reset(new MultiFab(ba_nd, m_dmap[lev], 1, ng, MFInfo(), *m_ebfactory[lev]));
    }

    // Get inputs from ParmParse
    readParameters();

    // object is ready
    m_ok = true;

    // First setup
    setup();
}
#endif

//
// Perform projection:
//
//     vel = vel - sigma*grad(phi)
//
//  where phi is the solution of
//
//   div( sigma * grad(phi) ) = div(vel) + S_cc + S_nd
//
//  where vel, sigma, S_cc are cell-centered variables
//  and phi and S_nd are nodal variables
//
//  grad(phi) is node-centered.
//
void
NodalProjector::project ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                          const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                          const amrex::Vector<amrex::MultiFab*>        a_S_cc,
                          const amrex::Vector<const amrex::MultiFab*>& a_S_nd )

{
    AMREX_ALWAYS_ASSERT(m_ok);
    BL_PROFILE("NodalProjector::project");

    if (m_verbose > 0)
        amrex::Print() << "Nodal Projection:" << std::endl;

    // Setup solver -- ALWAYS do this because matrix may change
    setup();

    // Compute RHS
    computeRHS( GetVecOfPtrs(m_rhs), a_vel, a_S_cc, a_S_nd );

    // Print diagnostics
    if (m_verbose > 0)
    {
        amrex::Print() << " >> Before projection:" << std::endl;
        printInfo();
    }

    // Set matrix coefficients
    for (int lev(0); lev < a_sigma.size(); ++lev)
        m_matrix -> setSigma(lev, *a_sigma[lev]);

    // Solve
    m_solver -> solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), m_mg_rtol, m_mg_atol );

    // Get fluxes -- fluxes = - sigma*grad(phi)
    m_solver -> getFluxes( GetVecOfPtrs(m_fluxes) );

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        // vel = vel + fluxes = vel - sigma * grad(phi)
        MultiFab::Add( *a_vel[lev], *m_fluxes[lev], 0, 0, AMREX_SPACEDIM, 0);

        // set m_fluxes = -fluxes/sigma = grad(phi)
        m_fluxes[lev] -> mult(- 1.0, m_fluxes[lev]->nGrow() );
        for (int n(0); n < AMREX_SPACEDIM; ++n)
            MultiFab::Divide(*m_fluxes[lev], *a_sigma[lev], 0, n, 1, m_fluxes[lev]->nGrow() );

        // Fill boundaries and apply scale factor to phi
        m_phi[lev] -> FillBoundary( m_geom[lev].periodicity());

    }

    // Print diagnostics
    if (m_verbose > 0)
    {
        computeRHS( GetVecOfPtrs(m_rhs), a_vel, a_S_cc, a_S_nd );

        amrex::Print() << " >> After projection:" << std::endl;
        printInfo();
    }
}


//
// Perform projection:
//
//     vel = vel - (1/beta)*grad(phi)
//
//  where phi is the solution of
//
//   div( (alpha/beta) * grad(phi) ) = RHS
//
//  where vel, alpha and beta are cell-centered variables
//  RHS is cell-centered
//
//  grad(phi) is node-centered.
//
void
NodalProjector::project2 ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                           const amrex::Vector<const amrex::MultiFab*>& a_alpha,
                           const amrex::Vector<const amrex::MultiFab*>& a_beta,
                           const amrex::Vector<const amrex::MultiFab*>& a_rhs )

{
    AMREX_ALWAYS_ASSERT(m_ok);
    BL_PROFILE("NodalProjector::project");

    if (m_verbose > 0)
        amrex::Print() << "Nodal Projection:" << std::endl;

    // Setup solver -- ALWAYS do this because matrix may change
    setup();

    // Print diagnostics
    if (m_verbose > 0)
    {
        amrex::Print() << " >> Before projection:" << std::endl;
        printInfo();
    }

    // Set matrix coefficients
    Vector< std::unique_ptr<MultiFab> > sigma(m_phi.size());
    for (int lev(0); lev < sigma.size(); ++lev)
    {
#ifdef AMREX_USE_EB
        sigma[lev].reset(new MultiFab( m_grids[lev], m_dmap[lev], 3, 1 , MFInfo(), *m_ebfactory[lev] ) );
#else
        sigma[lev].reset(new MultiFab( m_grids[lev], m_dmap[lev], 3, 1 ) );
#endif
        MultiFab::Copy(*sigma[lev],*a_alpha[lev],0,0,1,0);
        MultiFab::Divide(*sigma[lev],*a_beta[lev],0,0,1,0);

        m_matrix -> setSigma(lev, *sigma[lev]);
    }

    // Solve
    m_solver -> solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(a_rhs), m_mg_rtol, m_mg_atol );

    // Get fluxes -- fluxes = -  (alpha/beta) * grad(phi)
    m_solver -> getFluxes( GetVecOfPtrs(m_fluxes) );

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        // fluxes -> fluxes/alpha = -grad(phi)/beta
        MultiFab::Divide( *m_fluxes[lev], *a_alpha[lev], 0, 0, 1, 0 );
        MultiFab::Divide( *m_fluxes[lev], *a_alpha[lev], 0, 1, 1, 0 );
        MultiFab::Divide( *m_fluxes[lev], *a_alpha[lev], 0, 2, 1, 0 );

        // vel = vel + fluxes = vel - grad(phi) / beta
        MultiFab::Add( *a_vel[lev], *m_fluxes[lev], 0, 0, AMREX_SPACEDIM, 0);

        // set m_fluxes = -fluxes*beta = grad(phi)
        m_fluxes[lev] -> mult(- 1.0, m_fluxes[lev]->nGrow() );
        for (int n(0); n < AMREX_SPACEDIM; ++n)
            MultiFab::Multiply(*m_fluxes[lev], *a_beta[lev], 0, n, 1, m_fluxes[lev]->nGrow() );

        // Fill boundaries and apply scale factor to phi
        m_phi[lev] -> FillBoundary( m_geom[lev].periodicity());

    }
}

//
// Read from input file
//
void
NodalProjector::readParameters ()
{
    ParmParse pp("projection");
    pp.query( "verbose"                , m_verbose );
    pp.query( "mg_verbose"             , m_mg_verbose );
    pp.query( "mg_cg_verbose"          , m_mg_cg_verbose );
    pp.query( "mg_maxiter"             , m_mg_maxiter );
    pp.query( "mg_cg_maxiter"          , m_mg_cg_maxiter );
    pp.query( "mg_rtol"                , m_mg_rtol );
    pp.query( "mg_atol"                , m_mg_atol );
    pp.query( "mg_max_coarsening_level", m_mg_max_coarsening_level );
    pp.query( "bottom_solver_type"     , m_bottom_solver_type );
}



//
// Setup object before solve
//
void
NodalProjector::setup ()
{
    BL_PROFILE("NodalProjector::setup");
    AMREX_ALWAYS_ASSERT(m_ok);

    // Initialize all variables
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        m_phi[lev] -> setVal(0.0);
        m_fluxes[lev] -> setVal(0.0);
        m_rhs[lev] ->  setVal(0.0);
    }

    //
    // Setup Matrix
    //
    LPInfo                       info;
    info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
#ifdef AMREX_USE_EB
    m_matrix.reset(new MLNodeLaplacian(m_geom, m_grids, m_dmap, info, m_ebfactory));
#else
    m_matrix.reset(new MLNodeLaplacian(m_geom, m_grids, m_dmap, info));
#endif

    m_matrix->setGaussSeidel(true);
    m_matrix->setHarmonicAverage(false);
    m_matrix->setDomainBC(m_bc_lo, m_bc_hi);

    //
    // Setup solver
    //
    m_solver.reset(new MLMG(*m_matrix));

    m_solver->setMaxIter(m_mg_maxiter);
    m_solver->setVerbose(m_mg_verbose);
    m_solver->setCGVerbose(m_mg_cg_verbose);
    m_solver->setCGMaxIter(m_mg_cg_maxiter);

    if (m_bottom_solver_type == "smoother")
    {
        m_solver->setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (m_bottom_solver_type == "bicg")
    {
        m_solver->setBottomSolver(MLMG::BottomSolver::bicgstab);
    }
    else if (m_bottom_solver_type == "cg")
    {
        m_solver->setBottomSolver(MLMG::BottomSolver::cg);
    }
    else if (m_bottom_solver_type == "bicgcg")
    {
        m_solver->setBottomSolver(MLMG::BottomSolver::bicgcg);
    }
    else if (m_bottom_solver_type == "cgbicg")
    {
        m_solver->setBottomSolver(MLMG::BottomSolver::cgbicg);
    }
    else if (m_bottom_solver_type == "hypre")
    {
#ifdef AMREX_USE_HYPRE
        m_solver->setBottomSolver(MLMG::BottomSolver::hypre);
#else
        amrex::Abort("AMReX was not built with HYPRE support");
#endif
    }

}


//
// Compute RHS: div(u) + S_nd + S_cc
//
void
NodalProjector::computeRHS (  const amrex::Vector<amrex::MultiFab*>&       a_rhs,
                              const amrex::Vector<amrex::MultiFab*>&       a_vel,
                              const amrex::Vector<amrex::MultiFab*>&       a_S_cc,
                              const amrex::Vector<const amrex::MultiFab*>& a_S_nd )
{
    AMREX_ALWAYS_ASSERT(m_ok);
    AMREX_ALWAYS_ASSERT(a_rhs.size()==m_phi.size());
    AMREX_ALWAYS_ASSERT(a_vel.size()==m_phi.size());
    AMREX_ALWAYS_ASSERT((a_S_cc.size()==0) || (a_S_cc.size()==m_phi.size()) );
    AMREX_ALWAYS_ASSERT((a_S_nd.size()==0) || (a_S_nd.size()==m_phi.size()) );

    BL_PROFILE("NodalProjector::computeRHS");

    bool has_S_cc(!a_S_cc.empty());
    bool has_S_nd(!a_S_nd.empty());

    // Check the type of grids used
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_rhs[lev]->ixType().nodeCentered(), "a_rhs is not node centered");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_vel[lev]->ixType().cellCentered(), "a_vel is not cell centered");

        if (has_S_cc)
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_S_cc[lev]->ixType().cellCentered(),"a_S_cc is not cell centered");

        if (has_S_nd)
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_S_nd[lev]->ixType().nodeCentered(),"a_S_nd is not node centered");
    }

    m_matrix -> compRHS( a_rhs, a_vel, a_S_nd, a_S_cc );
}


void
NodalProjector::printInfo ()
{
    for (int lev(0); lev < m_rhs.size(); ++lev)
    {
        amrex::Print() << "  * On lev " << lev
                       << " max(abs(divu)) = "
                       << m_rhs[lev]->norm0(0,0,false,true)
                       << std::endl;
    }
}

}
