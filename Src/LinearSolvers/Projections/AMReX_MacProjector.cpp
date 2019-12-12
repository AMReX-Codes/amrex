#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MacProjector.H>
#include <AMReX_ParmParse.H>

namespace amrex {

MacProjector::MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                            const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                            const Vector<Geometry>& a_geom,
                            const LPInfo& a_lpinfo,
                            const Vector<MultiFab const*>& a_divu)
    : m_umac(a_umac),
      m_geom(a_geom)
{
    int nlevs = a_umac.size();
    Vector<BoxArray> ba(nlevs);
    Vector<DistributionMapping> dm(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        ba[ilev] = amrex::convert(a_umac[ilev][0]->boxArray(), IntVect::TheZeroVector());
        dm[ilev] = a_umac[ilev][0]->DistributionMap();
    }

    m_rhs.resize(nlevs);
    m_phi.resize(nlevs);
    m_fluxes.resize(nlevs);

#ifdef AMREX_USE_EB
    bool has_eb = a_umac[0][0]->hasEBFabFactory();
    if (has_eb)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_umac[0][0]->nGrow() > 0,
                                         "MacProjector: with EB, umac must have at least one ghost cell");

        m_eb_factory.resize(nlevs,nullptr);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_eb_factory[ilev] = dynamic_cast<EBFArrayBoxFactory const*>(&(a_umac[ilev][0]->Factory()));
            m_rhs[ilev].define(ba[ilev],dm[ilev],1,0,MFInfo(),a_umac[ilev][0]->Factory());
            m_phi[ilev].define(ba[ilev],dm[ilev],1,1,MFInfo(),a_umac[ilev][0]->Factory());
            m_rhs[ilev].setVal(0.0);
            m_phi[ilev].setVal(0.0);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_fluxes[ilev][idim].define(amrex::convert(ba[ilev],IntVect::TheDimensionVector(idim)),
                                            dm[ilev],1,0,MFInfo(),a_umac[ilev][0]->Factory());
            }
        }

        m_eb_abeclap.reset(new MLEBABecLap(a_geom, ba, dm, a_lpinfo, m_eb_factory));
        m_linop = m_eb_abeclap.get();

        m_eb_abeclap->setScalars(0.0, 1.0);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_eb_abeclap->setBCoeffs(ilev, a_beta[ilev]);
        }
    }
    else
#endif
    {
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_rhs[ilev].define(ba[ilev],dm[ilev],1,0);
            m_phi[ilev].define(ba[ilev],dm[ilev],1,1);
            m_rhs[ilev].setVal(0.0);
            m_phi[ilev].setVal(0.0);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_fluxes[ilev][idim].define(amrex::convert(ba[ilev],IntVect::TheDimensionVector(idim)),
                                            dm[ilev],1,0);
            }
        }

        m_abeclap.reset(new MLABecLaplacian(a_geom, ba, dm, a_lpinfo));
        m_linop = m_abeclap.get();

        m_abeclap->setScalars(0.0, 1.0);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            m_abeclap->setBCoeffs(ilev, a_beta[ilev]);
        }
    }

    for (int ilev = 0, N = a_divu.size(); ilev < N; ++ilev) {
        if (a_divu[ilev]) {
            MultiFab::Copy(m_rhs[ilev], *a_divu[ilev], 0, 0, 1, 0);
        }
    }

    m_mlmg.reset(new MLMG(*m_linop));

    setOptions();
}

void
MacProjector::setDomainBC (const Array<LinOpBCType,AMREX_SPACEDIM>& lobc,
                           const Array<LinOpBCType,AMREX_SPACEDIM>& hibc)
{
    m_linop->setDomainBC(lobc, hibc);
    for (int ilev = 0, N = m_umac.size(); ilev < N; ++ilev) {
        m_linop->setLevelBC(ilev, nullptr);
    }

    m_needs_domain_bcs = false;
}


void
MacProjector::setLevelBC (int amrlev, const MultiFab* levelbcdata)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_needs_domain_bcs,
                                     "setDomainBC must be called before setLevelBC");
    m_linop->setLevelBC(amrlev, levelbcdata);
}



void
MacProjector::project (Real reltol, Real atol, MLMG::Location loc)
{
    const int nlevs = m_rhs.size();

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            u[idim] = m_umac[ilev][idim];
        }
        MultiFab divu(m_rhs[ilev].boxArray(), m_rhs[ilev].DistributionMap(),
                      1, 0, MFInfo(), m_rhs[ilev].Factory());
#ifdef AMREX_USE_EB
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_umac[ilev][idim]->FillBoundary(m_geom[ilev].periodicity());
        }
        bool already_on_centroid = false;
        if (loc == MLMG::Location::FaceCentroid) already_on_centroid = true;
        EB_computeDivergence(divu, u, m_geom[ilev], already_on_centroid);
#else
        computeDivergence(divu, u, m_geom[ilev]);
#endif
        MultiFab::Subtract(m_rhs[ilev], divu, 0, 0, 1, 0);
    }

    m_mlmg->solve(amrex::GetVecOfPtrs(m_phi), amrex::GetVecOfConstPtrs(m_rhs), reltol, atol);

    m_mlmg->getFluxes(amrex::GetVecOfArrOfPtrs(m_fluxes), loc);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            MultiFab::Add(*m_umac[ilev][idim], m_fluxes[ilev][idim], 0, 0, 1, 0);
#ifdef AMREX_USE_EB
            EB_set_covered_faces(m_umac[ilev], 0.0);
#endif
        }
    }
}

void
MacProjector::project (const Vector<MultiFab*>& phi_inout, Real reltol, Real atol, MLMG::Location loc)
{

    const int nlevs = m_rhs.size();

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            u[idim] = m_umac[ilev][idim];
        }
        MultiFab divu(m_rhs[ilev].boxArray(), m_rhs[ilev].DistributionMap(),
                      1, 0, MFInfo(), m_rhs[ilev].Factory());
#ifdef AMREX_USE_EB
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_umac[ilev][idim]->FillBoundary(m_geom[ilev].periodicity());
        }
        bool already_on_centroid = false;
        if (loc == MLMG::Location::FaceCentroid) already_on_centroid = true;
        EB_computeDivergence(divu, u, m_geom[ilev], already_on_centroid);
#else
        computeDivergence(divu, u, m_geom[ilev]);
#endif
        MultiFab::Subtract(m_rhs[ilev], divu, 0, 0, 1, 0);

        MultiFab::Copy(m_phi[ilev], *phi_inout[ilev], 0, 0, 1, 0);
    }

    m_mlmg->solve(amrex::GetVecOfPtrs(m_phi), amrex::GetVecOfConstPtrs(m_rhs), reltol, atol);

    m_mlmg->getFluxes(amrex::GetVecOfArrOfPtrs(m_fluxes), loc);


    for (int ilev = 0; ilev < nlevs; ++ilev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            MultiFab::Add(*m_umac[ilev][idim], m_fluxes[ilev][idim], 0, 0, 1, 0);
#ifdef AMREX_USE_EB
            EB_set_covered_faces(m_umac[ilev], 0.0);
#endif
        }
    }
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        MultiFab::Copy(*phi_inout[ilev], m_phi[ilev], 0, 0, 1, 0);
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
    Real         bottom_rtol(1.0e-4);
    Real         bottom_atol(-1.0);
    std::string  bottom_solver("bicg");

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

    // Set default/input values
    m_linop->setMaxOrder(maxorder);
    m_mlmg->setVerbose(m_verbose);
    m_mlmg->setBottomVerbose(bottom_verbose);
    m_mlmg->setMaxIter(maxiter);
    m_mlmg->setBottomMaxIter(bottom_maxiter);
    m_mlmg->setBottomTolerance(bottom_rtol);
    m_mlmg->setBottomToleranceAbs(bottom_atol);

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

}
