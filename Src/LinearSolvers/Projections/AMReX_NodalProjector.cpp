#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_NodalProjector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {


//
// Add defaults for boundary conditions???
//
NodalProjector::NodalProjector ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                                 const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                                 const amrex::Vector<amrex::Geometry>&        a_geom,
                                 const LPInfo&                                a_lpinfo,
                                 const amrex::Vector<amrex::MultiFab*>&       a_S_cc,
                                 const amrex::Vector<const amrex::MultiFab*>& a_S_nd )
    : m_geom(a_geom),
      m_vel(a_vel),
      m_S_cc(a_S_cc),
      m_sigma(a_sigma),
      m_S_nd(a_S_nd)
{
    int nlevs = a_vel.size();

    Vector<BoxArray> ba(nlevs);
    Vector<DistributionMapping> dm(nlevs);
    for (int lev = 0; lev < nlevs; ++lev)
    {
        ba[lev] = a_vel[lev]->boxArray();
        dm[lev] = a_vel[lev]->DistributionMap();
    }

    // Resize member data
    m_phi.resize(nlevs);
    m_fluxes.resize(nlevs);
    m_rhs.resize(nlevs);


#ifdef AMREX_USE_EB
    bool has_eb = a_vel[0] -> hasEBFabFactory();
    if (has_eb)
    {
        m_ebfactory.resize(nlevs,nullptr);
        for (int lev = 0; lev < nlevs; ++lev )
        {
            m_ebfactory[lev] = dynamic_cast<EBFArrayBoxFactory const*>(&(a_vel[lev]->Factory()));

            // Cell-centered data
            m_fluxes[lev].define(ba[lev], dm[lev], AMREX_SPACEDIM, 0, MFInfo(), a_vel[lev]->Factory());

            // Node-centered data
            auto tmp = ba[lev];
            const auto& ba_nd = tmp.surroundingNodes();
            m_phi[lev].define(ba_nd, dm[lev], 1, 1, MFInfo(), a_vel[lev]->Factory());
            m_rhs[lev].define(ba_nd, dm[lev], 1, 0, MFInfo(), a_vel[lev]->Factory());
        }
    }
    else
#endif
    {
        for (int lev = 0; lev < nlevs; ++lev )
        {
            // Cell-centered data
            m_fluxes[lev].define(ba[lev], dm[lev], AMREX_SPACEDIM, 0);

            // Node-centered data
            BoxArray tmp = ba[lev];
            const auto& ba_nd = tmp.surroundingNodes();
            m_phi[lev].define(ba_nd, dm[lev], 1, 1);
            m_rhs[lev].define(ba_nd, dm[lev], 1, 0);
        }
    }

    // Initialize all variables
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        m_phi[lev].setVal(0.0);
        m_fluxes[lev].setVal(0.0);
        m_rhs[lev].setVal(0.0);
    }

    //
    // Setup linear operator
    //
#ifdef AMREX_USE_EB
    m_linop.reset(new MLNodeLaplacian(m_geom, ba, dm, a_lpinfo, m_ebfactory));
#else
    m_linop.reset(new MLNodeLaplacian(m_geom, ba, dm, a_lpinfo));
#endif

    m_linop->setGaussSeidel(true);
    m_linop->setHarmonicAverage(false);

    //
    // Setup solver
    //
    m_mlmg.reset(new MLMG(*m_linop));

    setOptions();
}


//
// Set options by using default values and values read in input file
//
void
NodalProjector::setOptions ()
{

    // Default values
    int          bottom_verbose(0);
    int          maxiter(100);
    int          bottom_maxiter(100);
    Real         bottom_rtol(1.0e-4);
    Real         bottom_atol(-1.0);
    std::string  bottom_solver("bicgcg");

    int          num_pre_smooth (2);
    int          num_post_smooth(2);

    // Read from input file
    ParmParse pp("nodal_proj");
    pp.query( "verbose"       , m_verbose );
    pp.query( "bottom_verbose", bottom_verbose );
    pp.query( "maxiter"       , maxiter );
    pp.query( "bottom_maxiter", bottom_maxiter );
    pp.query( "bottom_rtol"   , bottom_rtol );
    pp.query( "bottom_atol"   , bottom_atol );
    pp.query( "bottom_solver" , bottom_solver );

    pp.query( "num_pre_smooth"  , num_pre_smooth );
    pp.query( "num_post_smooth" , num_post_smooth );

    // Set default/input values
    m_mlmg->setVerbose(m_verbose);
    m_mlmg->setBottomVerbose(bottom_verbose);
    m_mlmg->setMaxIter(maxiter);
    m_mlmg->setBottomMaxIter(bottom_maxiter);
    m_mlmg->setBottomTolerance(bottom_rtol);
    m_mlmg->setBottomToleranceAbs(bottom_atol);

    m_mlmg->setPreSmooth(num_pre_smooth);
    m_mlmg->setPostSmooth(num_post_smooth);

    if (bottom_solver == "smoother")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (bottom_solver == "bicg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::bicgstab);
    }
    else if (bottom_solver == "cg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::cg);
    }
    else if (bottom_solver == "bicgcg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::bicgcg);
    }
    else if (bottom_solver == "cgbicg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::cgbicg);
    }
#ifdef AMREX_USE_HYPRE
    else if (bottom_solver == "hypre")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::hypre);
    }
#endif
}

void
NodalProjector::setDomainBC ( std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                              std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi )
{
    m_bc_lo=a_bc_lo;
    m_bc_hi=a_bc_hi;
    m_linop->setDomainBC(m_bc_lo,m_bc_hi);
    m_need_bcs = false;
}

void
NodalProjector::setCustomRHS (const amrex::Vector<const amrex::MultiFab*> a_rhs)
{
    AMREX_ALWAYS_ASSERT(m_rhs.size()==a_rhs.size());

    for (int lev=0; lev < m_rhs.size(); ++lev)
    {
        MultiFab::Copy(m_rhs[lev], *a_rhs[lev], 0, 0, 1, 0);
    }

    m_has_rhs = true;
}


void
NodalProjector::project ( Real a_rtol, Real a_atol )
{
    BL_PROFILE("NodalProjector::project");
    AMREX_ALWAYS_ASSERT(!m_need_bcs);

    if (m_verbose > 0)
        amrex::Print() << "Nodal Projection:" << std::endl;

    //
    // Average fine grid velocity values down to the coarse grid
    // By doing this operation at this time we:
    //
    //   1) fill regions covered by a finer grid with some "valid" data, i.e. not NaNs
    //      for example. This is required internally by MLMG.
    //   2) make the velocity field ready for projection on the entirety of each level.
    //
    averageDown(m_vel);

    // Set matrix coefficients
    for (int lev = 0; lev < m_sigma.size(); ++lev)
    {
        m_linop -> setSigma(lev, *m_sigma[lev]);
    }

    // Compute RHS if necessary
    if (!m_has_rhs)
    {
        computeRHS( GetVecOfPtrs(m_rhs), m_vel, m_S_cc, m_S_nd );
    }

    // Print diagnostics
    if (m_verbose > 0)
    {
        amrex::Print() << " >> Before projection:" << std::endl;
        printInfo();
        amrex::Print() << std::endl;
    }

    // Solve
    // phi comes out already averaged-down and ready to be used by caller if needed
    m_mlmg -> solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), a_rtol, a_atol );

    // Get fluxes -- fluxes = -  (alpha/beta) * grad(phi)
    m_mlmg -> getFluxes( GetVecOfPtrs(m_fluxes) );

    // At this time, the fluxes are "correct" only on regions not covered by finer grids.
    // We average the fluxes down so that they are "correct" everywhere in each level.
    // This is necessary because the caller can access grad(phi) and may use it for
    // computations involving the whole level.
    averageDown(GetVecOfPtrs(m_fluxes));

    // Compute sync residual BEFORE performing projection
    computeSyncResidual();

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        if (m_has_alpha)
        {
            // fluxes -> fluxes/alpha = -grad(phi)/beta
            for (int n = 0; n < AMREX_SPACEDIM; ++n)
            {
                MultiFab::Divide( m_fluxes[lev], *m_alpha[lev], 0, n, 1, 0 );
            }
        }

        //
        // vel = vel + fluxes = vel - grad(phi) / beta
        //
        // Since we already averaged-down the velocity field and -grad(phi),
        // we perform the projection by simply adding the two of them.
        // In virtue of the linearity of the operations involved, this is equivalent
        // to averaging down the velocity only once AFTER summing -grad(phi)
        //
        MultiFab::Add( *m_vel[lev], m_fluxes[lev], 0, 0, AMREX_SPACEDIM, 0);

        // set m_fluxes = grad(phi)
        m_fluxes[lev].mult(-1.0);
        for (int n = 0; n < AMREX_SPACEDIM; ++n)
        {
            if (m_has_alpha)
            {
                MultiFab::Multiply(m_fluxes[lev], *m_alpha[lev], 0, n, 1, 0);
            }
            MultiFab::Divide(m_fluxes[lev], *m_sigma[lev], 0, n, 1, 0);
        }

    }

    // Print diagnostics
    if ( (m_verbose > 0) && (!m_has_rhs))
    {
        computeRHS( GetVecOfPtrs(m_rhs), m_vel, m_S_cc, m_S_nd );
        amrex::Print() << " >> After projection:" << std::endl;
        printInfo();
        amrex::Print() << std::endl;
    }

}

void
NodalProjector::project ( const Vector<MultiFab*>& a_phi, Real a_rtol, Real a_atol )
{
    AMREX_ALWAYS_ASSERT(a_phi.size()==m_phi.size());

    for (int lev=0; lev < m_phi.size(); ++lev )
    {
        MultiFab::Copy(m_phi[lev],*a_phi[lev],0,0,1,m_phi[lev].nGrow());
    }

    project(a_rtol, a_atol);

    for (int lev=0; lev < m_phi.size(); ++lev )
    {
        MultiFab::Copy(*a_phi[lev],m_phi[lev],0,0,1,m_phi[lev].nGrow());
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
    AMREX_ALWAYS_ASSERT(!m_need_bcs); // This is needed to use linop
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

        if (has_S_cc && (a_S_cc[lev]!=nullptr) )
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_S_cc[lev]->ixType().cellCentered(),"a_S_cc is not cell centered");

        if (has_S_nd && (a_S_nd[lev]!=nullptr) )
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_S_nd[lev]->ixType().nodeCentered(),"a_S_nd is not node centered");
    }

    m_linop -> compRHS( a_rhs, a_vel, a_S_nd, a_S_cc );
}


void
NodalProjector::printInfo ()
{
    for (int lev(0); lev < m_rhs.size(); ++lev)
    {
        amrex::Print() << "  * On lev " << lev
                       << " max(abs(rhs)) = "
                       << m_rhs[lev].norm0(0,0,false,true)
                       << std::endl;
    }
}


void
NodalProjector::computeSyncResidual ()
{
    BL_PROFILE("NodalProjector::computeSyncResidual");

    if ( (m_sync_resid_fine != nullptr) || (m_sync_resid_crse != nullptr) )
    {
        int c_lev = 0;

        setCoarseBoundaryVelocityForSync();

        if (m_sync_resid_fine != nullptr)
        {
            MultiFab* rhptr = nullptr;
            if (!m_S_cc.empty())
                rhptr = m_S_cc[c_lev];
            m_linop->compSyncResidualFine(*m_sync_resid_fine, m_phi[c_lev], *m_vel[c_lev], rhptr);
        }

        if (m_sync_resid_crse != nullptr)
        {
            MultiFab* rhptr = nullptr;
            if (!m_S_cc.empty())
                rhptr = m_S_cc[c_lev];

            // this requires sigma to have 2 ghost cells (valid at -2)
            m_linop->compSyncResidualCoarse(*m_sync_resid_crse, m_phi[c_lev], *m_vel[c_lev],
                                           rhptr, m_fine_grids, m_ref_ratio);
        }


    }

}


// Set velocity in ghost cells to zero except for inflow
// 1) At non-inflow faces, the normal component of velocity will be completely zero'd
// 2) The normal velocity at corners -- even periodic corners -- just outside inflow faces will be zero'd
void
NodalProjector::setCoarseBoundaryVelocityForSync ()
{

    const BoxArray& grids = m_vel[0]->boxArray();
    const Box& domainBox  = m_geom[0].Domain();

    for (int idir=0; idir < AMREX_SPACEDIM; idir++)
    {
        if (m_bc_lo[idir] != LinOpBCType::inflow && m_bc_hi[idir] != LinOpBCType::inflow)
        {
            m_vel[0]->setBndry(0.0, idir, 1);
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_vel[0]); mfi.isValid(); ++mfi)
            {
                int i = mfi.index();
                FArrayBox& v_fab = (*m_vel[0])[mfi];

                const Box& reg = grids[i];
                const Box& bxg1 = amrex::grow(reg,1);
                BoxList bxlist(reg);

                if (m_bc_lo[idir] == LinOpBCType::inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir))
                {
                    Box bx;                // bx is the region we *protect* from zero'ing
                    bx = amrex::adjCellLo(reg, idir);
                    bxlist.push_back(bx);
                }

                if (m_bc_hi[idir] == LinOpBCType::inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
                    Box bx;                // bx is the region we *protect* from zero'ing
                    bx = amrex::adjCellHi(reg, idir);
                    bxlist.push_back(bx);
                }

                BoxList bxlist2 = amrex::complementIn(bxg1, bxlist);

                for (BoxList::iterator it=bxlist2.begin(); it != bxlist2.end(); ++it)
                {
                    Box ovlp = *it & v_fab.box();
                    if (ovlp.ok())
                    {
                        v_fab.setVal<RunOn::Host>(0.0, ovlp, idir, 1);
                    }
                }
            }
        }
    }

}


void
NodalProjector::averageDown (const amrex::Vector<amrex::MultiFab*> a_var)
{
    // If not cartesian, we should average down by using volume weighting
    // We check that coord sys is Cartesian only for coarsest level and assume
    // geom for all other levels are Cartesian as well
    AMREX_ALWAYS_ASSERT(m_geom[0].IsCartesian());

    int f_lev = a_var.size()-1;
    int c_lev = 0;

    for (int lev = f_lev-1; lev >= c_lev; --lev)
    {
        IntVect rr   = m_geom[lev+1].Domain().size() / m_geom[lev].Domain().size();

#ifdef AMREX_USE_EB
        EB_average_down(*a_var[lev+1], *a_var[lev], 0, a_var[lev]->nComp(), rr);
#else
        average_down(*a_var[lev+1], *a_var[lev], 0, a_var[lev]->nComp(), rr);
#endif

    }

}


}
