#include <FMultiGrid.H>

FMultiGrid::FMultiGrid (const Geometry            & geom, 
			const DistributionMapping & dmap,
			const BoxArray            & ba,
			int                         baselevel,
			IntVect                     crse_ratio)
    :
    m_nlevels(1),
    m_baselevel(baselevel),
    m_crse_ratio(crse_ratio),
    m_ncomp(1),
    m_stencil(CC_CROSS_STENCIL),
    m_verbose(0),
    m_geom(PArray<Geometry           >(m_nlevels)),
    m_dmap(PArray<DistributionMapping>(m_nlevels)),
    m_ba  (PArray<BoxArray           >(m_nlevels)),
    m_coeffs_set(0),
    m_bc_set(0),
    m_bndry(0),
    m_mgt_solver(0)
{
    m_geom.set(0, &geom);
    m_dmap.set(0, &dmap);
    m_ba  .set(0, &ba  );

    if (m_baselevel > 0 && m_crse_ratio == IntVect::TheZeroVector()) 
	BoxLib::Abort("FMultiGrid: must set crse_ratio if baselevel > 0");
}

FMultiGrid::FMultiGrid (const std::vector<Geometry>            & geom, 
			const std::vector<DistributionMapping> & dmap,
			const std::vector<BoxArray>            & ba,
			int                                      baselevel,
			IntVect                                  crse_ratio)
    :
    m_nlevels(ba.size()),
    m_baselevel(baselevel),
    m_crse_ratio(crse_ratio),
    m_ncomp(1),
    m_stencil(CC_CROSS_STENCIL),
    m_verbose(0),
    m_geom(PArray<Geometry           >(m_nlevels)),
    m_dmap(PArray<DistributionMapping>(m_nlevels)),
    m_ba  (PArray<BoxArray           >(m_nlevels)),
    m_coeffs_set(0),
    m_bc_set(0),
    m_bndry(0),
    m_mgt_solver(0)
{
    for (int lev=0; lev < m_nlevels; ++lev) {
	m_geom.set(lev, &geom[lev]);
	m_dmap.set(lev, &dmap[lev]);
	m_ba  .set(lev, &ba  [lev]);
    }

    if (m_baselevel > 0 && m_crse_ratio == IntVect::TheZeroVector()) 
	BoxLib::Abort("FMultiGrid: must set crse_ratio if baselevel > 0");
}

FMultiGrid::FMultiGrid (const PArray<Geometry>            & geom, 
					  const PArray<DistributionMapping> & dmap,
					  const PArray<BoxArray>            & ba,
					  int                                 baselevel,
			                  IntVect                             crse_ratio)
    :
    m_nlevels(ba.size()),
    m_baselevel(baselevel),
    m_crse_ratio(crse_ratio),
    m_ncomp(1),
    m_stencil(CC_CROSS_STENCIL),
    m_verbose(0),
    m_geom(PArray<Geometry           >(m_nlevels)),
    m_dmap(PArray<DistributionMapping>(m_nlevels)),
    m_ba  (PArray<BoxArray           >(m_nlevels)),
    m_coeffs_set(0),
    m_bc_set(0),
    m_bndry(0),
    m_mgt_solver(0)
{
    for (int lev=0; lev < m_nlevels; ++lev) {
	m_geom.set(lev, &geom[lev]);
	m_dmap.set(lev, &dmap[lev]);
	m_ba  .set(lev, &ba  [lev]);
    }

    if (m_baselevel > 0 && m_crse_ratio == IntVect::TheZeroVector()) 
	BoxLib::Abort("FMultiGrid: must set crse_ratio if baselevel > 0");
}

FMultiGrid::~FMultiGrid ()
{
    delete m_bndry;
    delete m_mgt_solver;
}

void
FMultiGrid::set_bc (int * mg_bc)
{
    BL_ASSERT(m_bc_set == 0);

    // The values of phys_bc & ref_ratio do not matter
    // because we are going to use mg_bc.

    IntVect ref_ratio = IntVect::TheZeroVector();
    Array<int> lo_bc(BL_SPACEDIM, 0);
    Array<int> hi_bc(BL_SPACEDIM, 0);
    BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());

    m_bndry = new MacBndry(m_ba[0], m_ncomp, m_geom[0]);
    m_bndry->setHomogValues(phys_bc, ref_ratio);

    for (int i=0; i < 2*BL_SPACEDIM; ++i) {
	m_mg_bc[i] = mg_bc[i];
    }

    m_bc_set = 1;
}

void
FMultiGrid::set_bc (int           * mg_bc,
		    const MultiFab& phi)
{
    BL_ASSERT(m_bc_set == 0);
    BL_ASSERT(phi.nComp() == m_ncomp);

    // The values of phys_bc & ref_ratio do not matter
    // because we are going to use mg_bc.

    Array<int> lo_bc(BL_SPACEDIM, 0);
    Array<int> hi_bc(BL_SPACEDIM, 0);
    BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());

    m_bndry = new MacBndry(m_ba[0], m_ncomp, m_geom[0]);
    m_bndry->setBndryValues(phi, 0, 0, m_ncomp, phys_bc); 

    for (int i=0; i < 2*BL_SPACEDIM; ++i) {
	m_mg_bc[i] = mg_bc[i];
    }

    m_bc_set = 1;
}

void
FMultiGrid::set_bc (int           * mg_bc,
		    const MultiFab& crse_phi,
		    const MultiFab& phi)
{
    if (m_crse_ratio == IntVect::TheZeroVector())
	BoxLib::Abort("FMultiGrid: calling wrong set_bc function");

    BL_ASSERT(m_bc_set == 0);
    BL_ASSERT(phi.nComp() == m_ncomp);
    BL_ASSERT(crse_phi.nComp() == m_ncomp);

    // The values of phys_bc & ref_ratio do not matter
    // because we are going to use mg_bc.

    Array<int> lo_bc(BL_SPACEDIM, 0);
    Array<int> hi_bc(BL_SPACEDIM, 0);
    BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());

    const int in_rad     = 0;
    const int out_rad    = 1;
    const int extent_rad = 2;

    BoxArray crse_boxes(phi.boxArray());
    crse_boxes.coarsen(m_crse_ratio);

    BndryRegister crse_br(crse_boxes, in_rad, out_rad, extent_rad, m_ncomp);
    crse_br.copyFrom(crse_phi, crse_phi.nGrow(), 0, 0, m_ncomp);

    m_bndry = new MacBndry(m_ba[0], m_ncomp, m_geom[0]);
    m_bndry->setBndryValues(crse_br, 0, phi, 0, 0, m_ncomp, m_crse_ratio, 
			    phys_bc);

    for (int i=0; i < 2*BL_SPACEDIM; ++i) {
	m_mg_bc[i] = mg_bc[i];
    }

    m_bc_set = 1;
}

void 
FMultiGrid::set_abeclap_coeffs (Real alpha,
				Real beta,
				const PArray<MultiFab>& b)
{
    if (m_bc_set == 0) 
	BoxLib::Abort("MultiFab: set_bc must be called before set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    Array<PArray<MultiFab> > bb(1);
    bb[0].resize(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	bb[0].set(i, &b[i]);
    }

    m_mgt_solver->set_abeclap_coeffs(alpha, beta, bb, xa, xb);

    m_coeffs_set = 1;
}

void 
FMultiGrid::set_abeclap_coeffs (const MultiFab& a,
				Real beta,
				const PArray<MultiFab>& b)
{
    if (m_bc_set == 0) 
	BoxLib::Abort("MultiFab: set_bc must be called before set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    PArray<MultiFab> aa(1);
    aa.set(0, &a);

    Array<PArray<MultiFab> > bb(1);
    bb[0].resize(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	bb[0].set(i, &b[i]);
    }

    m_mgt_solver->set_abeclap_coeffs(aa, beta, bb, xa, xb);

    m_coeffs_set = 1;
}

void 
FMultiGrid::set_gravity_coeffs (const PArray<MultiFab>& b)
{
    if (m_bc_set == 0) 
	BoxLib::Abort("MultiFab: set_bc must be called calling set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    Array<PArray<MultiFab> > bb(1);
    bb[0].resize(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	bb[0].set(i, &b[i]);
    }

    m_mgt_solver->set_gravity_coefficients(bb, xa, xb);

    m_coeffs_set = 1;
}

void 
FMultiGrid::set_abeclap_coeffs (Real alpha,
				Real beta,
				const Array<PArray<MultiFab> >& b)
{
    if (m_bc_set == 0) 
	BoxLib::Abort("MultiFab: set_bc must be called before set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    m_mgt_solver->set_abeclap_coeffs(alpha, beta, b, xa, xb);

    m_coeffs_set = 1;
}

void 
FMultiGrid::set_abeclap_coeffs (const PArray<MultiFab>& a,
				Real beta,
				const Array<PArray<MultiFab> >& b)
{
    if (m_bc_set == 0) 
	BoxLib::Abort("MultiFab: set_bc must be called before set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    m_mgt_solver->set_abeclap_coeffs(a, beta, b, xa, xb);

    m_coeffs_set = 1;
}

void 
FMultiGrid::set_gravity_coeffs (const Array< PArray<MultiFab> >& b)
{
    if (m_bc_set == 0) 
	BoxLib::Abort("MultiFab: set_bc must be called calling set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    m_mgt_solver->set_gravity_coefficients(b, xa, xb);

    m_coeffs_set = 1;
}

void 
FMultiGrid::set_const_gravity_coeffs ()
{
    if (m_bc_set == 0) 
	BoxLib::Abort("FMultiGrid: set_bc must be called before set_*_coeffs is called");

    BL_ASSERT(m_coeffs_set == 0);

    build_mgt_solver();

    Array< Array<Real> > xa, xb;
    make_xa_xb(xa, xb);

    m_mgt_solver->set_const_gravity_coeffs(xa, xb);

    m_coeffs_set = 1;
}

Real
FMultiGrid::solve (MultiFab& phi,
		   MultiFab& rhs,
		   Real rel_tol, Real abs_tol,
		   int always_use_bnorm,
		   int need_grad_phi,
		   int verbose)
{
    if (m_bc_set == 0)
	BoxLib::Abort("FMultiGrid::solve: set_bc not called");

    if (m_coeffs_set == 0)
	BoxLib::Abort("FMultiGrid::solve: coeffs not set");

    if (m_nlevels > 1)
	BoxLib::Abort("FMultiGrid::solve: inconsistent # of levels");

    MultiFab* phi_p[1];
    MultiFab* rhs_p[1];
    phi_p[0] = &phi;
    rhs_p[0] = &rhs;

    Real final_resnorm;
    m_mgt_solver->solve(phi_p, rhs_p, *m_bndry, rel_tol, abs_tol, 
			always_use_bnorm, final_resnorm, need_grad_phi);
    return final_resnorm;
}

Real
FMultiGrid::solve (PArray<MultiFab>& phi,
		   PArray<MultiFab>& rhs,
		   Real rel_tol, Real abs_tol,
		   int always_use_bnorm,
		   int need_grad_phi,
		   int verbose)
{
    if (m_bc_set == 0)
	BoxLib::Abort("FMultiGrid::solve: set_bc not called");

    if (m_coeffs_set == 0)
	BoxLib::Abort("FMultiGrid::solve: coeffs not set");

    MultiFab* phi_p[m_nlevels];
    MultiFab* rhs_p[m_nlevels];
    for (int lev=0; lev < m_nlevels; ++lev) {
      phi_p[lev] = &phi[lev];
      rhs_p[lev] = &rhs[lev];
    }    

    Real final_resnorm;
    m_mgt_solver->solve(phi_p, rhs_p, *m_bndry, rel_tol, abs_tol, 
			always_use_bnorm, final_resnorm, need_grad_phi);
    return final_resnorm;
}

void
FMultiGrid::get_grad_phi (int amr_lev, PArray<MultiFab>& grad_phi)
{
    int lev = amr_lev - m_baselevel;
    const Real* dx = m_geom[lev].CellSize();
    m_mgt_solver->get_fluxes(lev, grad_phi, dx);
}

void 
FMultiGrid::compute_residual (MultiFab & phi,
			      MultiFab & rhs,
			      MultiFab & res)
{
    if (m_bc_set == 0)
	BoxLib::Abort("FMultiGrid::solve: set_bc not called");

    if (m_coeffs_set == 0)
	BoxLib::Abort("FMultiGrid::solve: coeffs not set");

    if (m_nlevels > 1)
	BoxLib::Abort("FMultiGrid::solve: inconsistent # of levels");

    MultiFab* phi_p[1];
    MultiFab* rhs_p[1];
    MultiFab* res_p[1];
    phi_p[0] = &phi;
    rhs_p[0] = &rhs;
    res_p[0] = &res;

    m_mgt_solver->compute_residual(phi_p, rhs_p, res_p, *m_bndry);
}

void 
FMultiGrid::compute_residual (PArray<MultiFab> & phi,
			      PArray<MultiFab> & rhs,
			      PArray<MultiFab> & res)
{
    if (m_bc_set == 0)
	BoxLib::Abort("FMultiGrid::solve: set_bc not called");

    if (m_coeffs_set == 0)
	BoxLib::Abort("FMultiGrid::solve: coeffs not set");

    MultiFab* phi_p[m_nlevels];
    MultiFab* rhs_p[m_nlevels];
    MultiFab* res_p[m_nlevels];
    for (int lev=0; lev < m_nlevels; ++lev) {
      phi_p[lev] = &phi[lev];
      rhs_p[lev] = &rhs[lev];
      res_p[lev] = &res[lev];
    }    

    m_mgt_solver->compute_residual(phi_p, rhs_p, res_p, *m_bndry);
}

void
FMultiGrid::build_mgt_solver ()
{
    BL_ASSERT(m_mgt_solver == 0);
    BL_ASSERT(m_bc_set != 0);

    std::vector<Geometry>            geom(m_nlevels);
    std::vector<DistributionMapping> dmap(m_nlevels);
    std::vector<BoxArray>              ba(m_nlevels);
    for (int lev=0; lev < m_nlevels; ++lev) {
	geom[lev] = m_geom[lev];
	dmap[lev] = m_dmap[lev];
	ba  [lev] = m_ba  [lev];
    }

    bool nodal = false;
    int  nc = 0;
    m_mgt_solver = new MGT_Solver(geom, m_mg_bc, ba, dmap, nodal, m_stencil,
				  nodal, nc, m_ncomp, m_verbose);
}

void
FMultiGrid::make_xa_xb (Array < Array<Real> > & xa,
				 Array < Array<Real> > & xb)
{
    BL_ASSERT( m_baselevel >= 0 );
    BL_ASSERT( m_baselevel == 0 || m_crse_ratio != IntVect::TheZeroVector() );

    xa.resize(m_nlevels);
    xb.resize(m_nlevels);

    for (int lev=0; lev < m_nlevels; ++lev) {
	xa[lev].resize(BL_SPACEDIM);
	xb[lev].resize(BL_SPACEDIM);
	if (lev + m_baselevel == 0) {
	    // For level 0, the boundary lives exactly on the faces
	    for (int n=0; n<BL_SPACEDIM; n++) {
		xa[lev][n] = 0.0;
		xb[lev][n] = 0.0;
	    }
	} else if (lev == 0) {
	    const Real* dx = m_geom[0].CellSize();
	    for (int n=0; n<BL_SPACEDIM; n++) {
		xa[lev][n] = 0.5 * m_crse_ratio[n] * dx[n];
		xb[lev][n] = 0.5 * m_crse_ratio[n] * dx[n];
	    }	    
	} else {
	    const Real* dx_crse = m_geom[lev-1].CellSize();
	    for (int n=0; n<BL_SPACEDIM; n++) {
		xa[lev][n] = 0.5 * dx_crse[n];
		xb[lev][n] = 0.5 * dx_crse[n];
	    }
	}
    }
}


