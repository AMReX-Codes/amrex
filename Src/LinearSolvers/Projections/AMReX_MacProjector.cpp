#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MacProjector.H>
#include <AMReX_ParmParse.H>

namespace amrex {

MacProjector::MacProjector(
    const Vector<Geometry>& a_geom,
    MLMG::Location a_umac_loc,
    MLMG::Location a_beta_loc,
    MLMG::Location a_phi_loc,
    MLMG::Location a_divu_loc)
    : m_geom(a_geom),
      m_needs_level_bcs(a_geom.size(),true),
      m_umac_loc(a_umac_loc),
      m_beta_loc(a_beta_loc),
      m_phi_loc(a_phi_loc),
      m_divu_loc(a_divu_loc)
{
    amrex::ignore_unused(m_divu_loc, m_beta_loc, m_phi_loc, m_umac_loc);
}

MacProjector::MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                            MLMG::Location a_umac_loc,
                            const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                            MLMG::Location a_beta_loc,
                            MLMG::Location  a_phi_loc,
                            const Vector<Geometry>& a_geom,
                            const LPInfo& a_lpinfo,
                            const Vector<MultiFab const*>& a_divu,
                            MLMG::Location a_divu_loc,
                            const Vector<iMultiFab const*>& a_overset_mask)
    : m_umac(a_umac),
      m_geom(a_geom),
      m_needs_level_bcs(a_geom.size(),true),
      m_umac_loc(a_umac_loc),
      m_beta_loc(a_beta_loc),
      m_phi_loc(a_phi_loc),
      m_divu_loc(a_divu_loc)
{
    amrex::ignore_unused(m_divu_loc, m_beta_loc, m_phi_loc, m_umac_loc);
    initProjector(a_lpinfo, a_beta, a_overset_mask);
    setDivU(a_divu);
}

void MacProjector::initProjector (
    const LPInfo& a_lpinfo,
    const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
    const Vector<iMultiFab const*>& a_overset_mask)
{
    const int nlevs = a_beta.size();
    Vector<BoxArray> ba(nlevs);
    Vector<DistributionMapping> dm(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        ba[ilev] = amrex::convert(
            a_beta[ilev][0]->boxArray(), IntVect::TheZeroVector());
        dm[ilev] = a_beta[ilev][0]->DistributionMap();
    }

    m_rhs.resize(nlevs);
    m_phi.resize(nlevs);
    m_fluxes.resize(nlevs);
    m_divu.resize(nlevs);

#ifdef AMREX_USE_EB
    bool has_eb = a_beta[0][0]->hasEBFabFactory();
    if (has_eb) {
        m_eb_factory.resize(nlevs, nullptr);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_eb_factory[ilev] = dynamic_cast<EBFArrayBoxFactory const*>(
                &(a_beta[ilev][0]->Factory()));
            m_rhs[ilev].define(
                ba[ilev], dm[ilev], 1, 0, MFInfo(), a_beta[ilev][0]->Factory());
            m_phi[ilev].define(
                ba[ilev], dm[ilev], 1, 1, MFInfo(), a_beta[ilev][0]->Factory());
            m_rhs[ilev].setVal(0.0);
            m_phi[ilev].setVal(0.0);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_fluxes[ilev][idim].define(
                    amrex::convert(ba[ilev], IntVect::TheDimensionVector(idim)),
                    dm[ilev], 1, 0, MFInfo(), a_beta[ilev][0]->Factory());
            }
        }

        m_eb_abeclap = std::make_unique<MLEBABecLap>(m_geom, ba, dm, a_lpinfo, m_eb_factory);
        m_linop = m_eb_abeclap.get();

        if (m_phi_loc == MLMG::Location::CellCentroid)
            m_eb_abeclap->setPhiOnCentroid();

        m_eb_abeclap->setScalars(0.0, 1.0);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_eb_abeclap->setBCoeffs(ilev, a_beta[ilev], m_beta_loc);
        }
    } else
#endif
    {
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_rhs[ilev].define(ba[ilev], dm[ilev], 1, 0);
            m_phi[ilev].define(ba[ilev], dm[ilev], 1, 1);
            m_rhs[ilev].setVal(0.0);
            m_phi[ilev].setVal(0.0);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_fluxes[ilev][idim].define(
                    amrex::convert(ba[ilev], IntVect::TheDimensionVector(idim)),
                    dm[ilev], 1, 0);
            }
        }

        if (a_overset_mask.empty()) {
            m_abeclap = std::make_unique<MLABecLaplacian>(m_geom, ba, dm, a_lpinfo);
        } else {
            m_abeclap = std::make_unique<MLABecLaplacian>(m_geom, ba, dm, a_overset_mask, a_lpinfo);
        }

        m_linop = m_abeclap.get();

        m_abeclap->setScalars(0.0, 1.0);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_abeclap->setBCoeffs(ilev, a_beta[ilev]);
        }
    }

    m_mlmg = std::make_unique<MLMG>(*m_linop);

    setOptions();

    m_needs_init = false;
}

void MacProjector::updateBeta (
    const Vector<Array<MultiFab const*, AMREX_SPACEDIM>>& a_beta)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_linop != nullptr,
        "MacProjector::updateBeta: initProjector must be called before calling this method");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_poisson == nullptr,
        "MacProjector::updateBeta: should not be called for constant beta");

    const int nlevs = a_beta.size();
#ifdef AMREX_USE_EB
    const bool has_eb = a_beta[0][0]->hasEBFabFactory();
    if (has_eb) {
        for (int ilev=0; ilev < nlevs; ++ilev)
            m_eb_abeclap->setBCoeffs(ilev, a_beta[ilev], m_beta_loc);
    } else
#endif
    {
        for (int ilev=0; ilev < nlevs; ++ilev)
            m_abeclap->setBCoeffs(ilev, a_beta[ilev]);
    }
}

void MacProjector::setUMAC(
    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& a_umac)
{
    m_umac = a_umac;
}

void MacProjector::setDivU(const Vector<MultiFab const*>& a_divu)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_linop != nullptr,
        "MacProjector::setDivU: initProjector must be called before calling this method");

    for (int ilev = 0, N = a_divu.size(); ilev < N; ++ilev) {
        if (a_divu[ilev]) {
            if (!m_divu[ilev].ok()) {
#ifdef AMREX_USE_EB
                m_divu[ilev].define(
                    a_divu[ilev]->boxArray(),
                    a_divu[ilev]->DistributionMap(),
                    1,0,MFInfo(), a_divu[ilev]->Factory());
#else
                m_divu[ilev].define(
                    a_divu[ilev]->boxArray(),
                    a_divu[ilev]->DistributionMap(), 1, 0);
#endif
            }
            MultiFab::Copy(m_divu[ilev], *a_divu[ilev], 0, 0, 1, 0);
        }
    }
}

void
MacProjector::setDomainBC (const Array<LinOpBCType,AMREX_SPACEDIM>& lobc,
                           const Array<LinOpBCType,AMREX_SPACEDIM>& hibc)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_linop != nullptr,
        "MacProjector::setDomainBC: initProjector must be called before calling this method");
    m_linop->setDomainBC(lobc, hibc);
    m_needs_domain_bcs = false;
}


void
MacProjector::setLevelBC (int amrlev, const MultiFab* levelbcdata)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_needs_domain_bcs,
                                     "setDomainBC must be called before setLevelBC");
    m_linop->setLevelBC(amrlev, levelbcdata);
    m_needs_level_bcs[amrlev] = false;
}



void
MacProjector::project (Real reltol, Real atol)
{
    const int nlevs = m_rhs.size();

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        if (m_needs_level_bcs[ilev]) {
            m_linop->setLevelBC(ilev, nullptr);
            m_needs_level_bcs[ilev] = false;
        }
    }

    averageDownVelocity();

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            u[idim] = m_umac[ilev][idim];
        }
        MultiFab divu(m_rhs[ilev].boxArray(), m_rhs[ilev].DistributionMap(),
                      1, 0, MFInfo(), m_rhs[ilev].Factory());
#ifdef AMREX_USE_EB
        if (m_umac_loc != MLMG::Location::FaceCentroid)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_umac[ilev][idim]->nGrow() > 0,
                                                 "MacProjector: with EB, umac must have at least one ghost cell if not already_on_centroid");
                m_umac[ilev][idim]->FillBoundary(m_geom[ilev].periodicity());
            }
        }

        EB_computeDivergence(divu, u, m_geom[ilev], (m_umac_loc == MLMG::Location::FaceCentroid));
#else
        computeDivergence(divu, u, m_geom[ilev]);
#endif

        // For mlabeclaplacian, we solve -del dot (beta grad phi) = rhs
        //   and set up RHS as (m_divu - divu), where m_divu is a user-provided source term
        // For mlpoisson, we solve `del dot grad phi = rhs/(-const_beta)`
        //   and set up RHS as (m_divu - divu)*(-1/const_beta)
        AMREX_ASSERT(m_poisson == nullptr || m_const_beta != Real(0.0));
        MultiFab::Copy(m_rhs[ilev], divu, 0, 0, 1, 0);
        m_rhs[ilev].mult(m_poisson ? Real(1.0)/m_const_beta : Real(-1.0));

        if (m_divu[ilev].ok())
        {
            MultiFab::Saxpy(m_rhs[ilev], m_poisson ? Real(-1.0)/m_const_beta : Real(1.0),
                            m_divu[ilev], 0, 0, 1, 0);
        }

        // Always reset initial phi to be zero. This is needed to handle the
        // situation where the MacProjector is being reused.
        m_phi[ilev].setVal(0.0);
    }

    m_mlmg->solve(amrex::GetVecOfPtrs(m_phi), amrex::GetVecOfConstPtrs(m_rhs), reltol, atol);

    m_mlmg->getFluxes(amrex::GetVecOfArrOfPtrs(m_fluxes), m_umac_loc);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (m_poisson) {
                MultiFab::Saxpy(*m_umac[ilev][idim], m_const_beta, m_fluxes[ilev][idim], 0,0,1,0);
            } else {
                MultiFab::Add(*m_umac[ilev][idim], m_fluxes[ilev][idim], 0, 0, 1, 0);
            }
#ifdef AMREX_USE_EB
            EB_set_covered_faces(m_umac[ilev], 0.0);
#endif
        }
    }

    averageDownVelocity();
}

void
MacProjector::project (const Vector<MultiFab*>& phi_inout, Real reltol, Real atol)
{
    const int nlevs = m_rhs.size();
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        MultiFab::Copy(m_phi[ilev], *phi_inout[ilev], 0, 0, 1, 0);
    }

    project(reltol, atol);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        MultiFab::Copy(*phi_inout[ilev], m_phi[ilev], 0, 0, 1, 0);
    }
}

void
MacProjector::getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                         const Vector<MultiFab*>& a_sol, MLMG::Location a_loc) const
{
    m_linop->getFluxes(a_flux, a_sol, a_loc);
    if (m_poisson) {
        for (auto const& mfarr : a_flux) {
            for (auto const& mfp : mfarr) {
                mfp->mult(m_const_beta);
            }
        }
    }
}

//
// Set options by using default values and values read in input file
//
void
MacProjector::setOptions ()
{
    // Default values
    int          maxorder(3);
    int          bottom_verbose(0);
    int          maxiter(200);
    int          bottom_maxiter(200);
    Real         bottom_rtol(1.0e-4_rt);
    Real         bottom_atol(-1.0_rt);
    std::string  bottom_solver("bicg");

    int num_pre_smooth(2);
    int num_post_smooth(2);

    // Read from input file
    ParmParse pp("mac_proj");
    pp.query( "verbose"       , m_verbose );
    pp.query( "maxorder"      , maxorder );
    pp.query( "bottom_verbose", bottom_verbose );
    pp.query( "maxiter"       , maxiter );
    pp.query( "bottom_maxiter", bottom_maxiter );
    pp.query( "bottom_rtol"   , bottom_rtol );
    pp.query( "bottom_atol"   , bottom_atol );
    pp.query( "bottom_solver" , bottom_solver );

    pp.query( "num_pre_smooth"  , num_pre_smooth );
    pp.query( "num_post_smooth" , num_post_smooth );

    // Set default/input values
    m_linop->setMaxOrder(maxorder);
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
    else if (bottom_solver == "hypre")
    {
#ifdef AMREX_USE_HYPRE
        m_mlmg->setBottomSolver(MLMG::BottomSolver::hypre);
#else
        amrex::Abort("AMReX was not built with HYPRE support");
#endif
    }
}

void
MacProjector::averageDownVelocity ()
{
    int finest_level = m_umac.size() - 1;


    for (int lev = finest_level; lev > 0; --lev)
    {

        IntVect rr  = m_geom[lev].Domain().size() / m_geom[lev-1].Domain().size();

#ifdef AMREX_USE_EB
        EB_average_down_faces(GetArrOfConstPtrs(m_umac[lev]),
                              m_umac[lev-1],
                              rr, m_geom[lev-1]);
#else
        average_down_faces(GetArrOfConstPtrs(m_umac[lev]),
                           m_umac[lev-1],
                           rr, m_geom[lev-1]);
#endif
    }
}

#ifndef AMREX_USE_EB
void MacProjector::initProjector (Vector<BoxArray> const& a_grids,
                                  Vector<DistributionMapping> const& a_dmap,
                                  const LPInfo& a_lpinfo, Real const a_const_beta,
                                  const Vector<iMultiFab const*>& a_overset_mask)
{
    m_const_beta = a_const_beta;

    const int nlevs = a_grids.size();
    Vector<BoxArray> ba(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        ba[ilev] = amrex::convert(a_grids[ilev], IntVect::TheZeroVector());
    }
    auto const& dm = a_dmap;

    m_rhs.resize(nlevs);
    m_phi.resize(nlevs);
    m_fluxes.resize(nlevs);
    m_divu.resize(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        m_rhs[ilev].define(ba[ilev], dm[ilev], 1, 0);
        m_phi[ilev].define(ba[ilev], dm[ilev], 1, 1);
        m_rhs[ilev].setVal(0.0);
        m_phi[ilev].setVal(0.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_fluxes[ilev][idim].define(
                amrex::convert(ba[ilev], IntVect::TheDimensionVector(idim)),
                dm[ilev], 1, 0);
        }
    }

    if (a_overset_mask.empty()) {
        m_poisson = std::make_unique<MLPoisson>(m_geom, ba, dm, a_lpinfo);
    } else {
        m_poisson = std::make_unique<MLPoisson>(m_geom, ba, dm, a_overset_mask, a_lpinfo);
    }

    m_linop = m_poisson.get();

    m_mlmg = std::make_unique<MLMG>(*m_linop);

    setOptions();

    m_needs_init = false;
}

MacProjector::MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                            const Real a_const_beta,
                            const Vector<Geometry>& a_geom,
                            const LPInfo& a_lpinfo,
                            const Vector<iMultiFab const*>& a_overset_mask,
                            const Vector<MultiFab const*>& a_divu)
    : m_const_beta(a_const_beta),
      m_umac(a_umac),
      m_geom(a_geom),
      m_umac_loc(MLMG::Location::FaceCenter),
      m_beta_loc(MLMG::Location::FaceCenter),
      m_phi_loc(MLMG::Location::CellCenter),
      m_divu_loc(MLMG::Location::CellCenter)
{
    const int nlevs = a_umac.size();
    Vector<BoxArray> ba(nlevs);
    Vector<DistributionMapping> dm(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        ba[ilev] = a_umac[ilev][0]->boxArray();
        dm[ilev] = a_umac[ilev][0]->DistributionMap();
    }
    initProjector(ba, dm, a_lpinfo, a_const_beta, a_overset_mask);
    setDivU(a_divu);
}

void MacProjector::updateBeta (Real a_const_beta)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_linop != nullptr,
        "MacProjector::updateBeta: initProjector must be called before calling this method");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_poisson != nullptr,
        "MacProjector::updateBeta: should not be called for variable beta");

    m_const_beta = a_const_beta;
}

#endif

}
