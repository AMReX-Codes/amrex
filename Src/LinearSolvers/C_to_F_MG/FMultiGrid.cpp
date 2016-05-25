#include <FMultiGrid.H>

FMultiGrid::FMultiGrid (const Geometry & geom, 
			int              baselevel,
			IntVect          crse_ratio)
    :
    m_nlevels(1),
    m_baselevel(baselevel),
    m_crse_ratio(crse_ratio),
    m_stencil(CC_CROSS_STENCIL),
    m_maxorder(0),
    m_verbose(0),
    m_geom(m_nlevels),
    m_bndry(0),
    m_mgt_solver(0)
{
    m_geom[0] = geom;

    if (m_baselevel > 0 && m_crse_ratio == IntVect::TheZeroVector()) 
	BoxLib::Abort("FMultiGrid: must set crse_ratio if baselevel > 0");
}

FMultiGrid::FMultiGrid (const std::vector<Geometry> & geom, 
			int                           baselevel,
			IntVect                       crse_ratio)
    :
    m_nlevels(geom.size()),
    m_baselevel(baselevel),
    m_crse_ratio(crse_ratio),
    m_stencil(CC_CROSS_STENCIL),
    m_maxorder(0),
    m_verbose(0),
    m_geom(m_nlevels),
    m_bndry(0),
    m_mgt_solver(0)
{
    m_geom = geom;

    if (m_baselevel > 0 && m_crse_ratio == IntVect::TheZeroVector()) 
	BoxLib::Abort("FMultiGrid: must set crse_ratio if baselevel > 0");
}

FMultiGrid::FMultiGrid (const PArray<Geometry> & geom, 
			int                      baselevel,
			IntVect                  crse_ratio)
    :
    m_nlevels(geom.size()),
    m_baselevel(baselevel),
    m_crse_ratio(crse_ratio),
    m_stencil(CC_CROSS_STENCIL),
    m_maxorder(0),
    m_verbose(0),
    m_geom(m_nlevels),
    m_bndry(0),
    m_mgt_solver(0)
{
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
	m_geom[ilev] = geom[ilev];
    }

    if (m_baselevel > 0 && m_crse_ratio == IntVect::TheZeroVector()) 
	BoxLib::Abort("FMultiGrid: must set crse_ratio if baselevel > 0");
}

void
FMultiGrid::set_bc (int * mg_bc)
{
    BL_ASSERT(!m_bc.initilized);
    m_bc = Boundary(mg_bc);
}

void
FMultiGrid::set_bc (int     * mg_bc,
		    MultiFab& phi)
{
    BL_ASSERT(!m_bc.initilized);

    m_bc = Boundary(mg_bc, 0, &phi);
}

void
FMultiGrid::set_bc (int     * mg_bc,
		    MultiFab& crse_phi,
		    MultiFab& phi)
{
    BL_ASSERT(m_crse_ratio != IntVect::TheZeroVector());
    BL_ASSERT(!m_bc.initilized);

    m_bc = Boundary(mg_bc, &crse_phi, &phi);
}

void
FMultiGrid::set_bc (const MacBndry& mac_bndry)
{
    BL_ASSERT(m_bndry == 0);
    m_bndry = const_cast<MacBndry*>(&mac_bndry);
}

void 
FMultiGrid::set_const_gravity_coeffs ()
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq);

    m_coeff.eq_type = const_gravity_eq;
}

void 
FMultiGrid::set_gravity_coeffs (PArray<MultiFab>& b)
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq);
    BL_ASSERT(m_nlevels == 1);
    BL_ASSERT(b.size() == BL_SPACEDIM);

    m_coeff.eq_type = gravity_eq;
    m_coeff.coeffs_set = true;

    Copy(m_coeff.b, b);
}

void 
FMultiGrid::set_gravity_coeffs (Array< PArray<MultiFab> >& b)
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq);
    BL_ASSERT(b.size() == m_nlevels);
    BL_ASSERT(b[0].size() == BL_SPACEDIM);

    m_coeff.eq_type = gravity_eq;
    m_coeff.coeffs_set = true;

    Copy(m_coeff.b, b);
}

void 
FMultiGrid::set_mac_coeffs (PArray<MultiFab>& b)
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq);
    BL_ASSERT(m_nlevels == 1);
    BL_ASSERT(b.size() == BL_SPACEDIM);

    m_coeff.eq_type = macproj_eq;
    m_coeff.coeffs_set = true;

    Copy(m_coeff.b, b);
}

void
FMultiGrid::set_scalars (Real alpha, Real beta)
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq || m_coeff.eq_type == general_eq);
    BL_ASSERT(!m_coeff.scalars_set);
    
    m_coeff.eq_type = general_eq;
    m_coeff.scalars_set = true;

    m_coeff.alpha = alpha;
    m_coeff.beta = beta;
}

void
FMultiGrid::set_coefficients (const MultiFab& a, PArray<MultiFab> & b)
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq || m_coeff.eq_type == general_eq);
    BL_ASSERT(!m_coeff.coeffs_set);
    BL_ASSERT(m_nlevels == 1);
    
    m_coeff.eq_type = general_eq;
    m_coeff.coeffs_set   = true;

    Copy(m_coeff.a, a);
    Copy(m_coeff.b, b);
}

void
FMultiGrid::set_coefficients (PArray<MultiFab>& a, Array<PArray<MultiFab> > & b)
{
    BL_ASSERT(m_coeff.eq_type == invalid_eq || m_coeff.eq_type == general_eq);
    BL_ASSERT(!m_coeff.coeffs_set);
    BL_ASSERT(m_nlevels == a.size() && m_nlevels == b.size());
    BL_ASSERT(BL_SPACEDIM == b[0].size());
    
    m_coeff.eq_type = general_eq;
    m_coeff.coeffs_set   = true;

    Copy(m_coeff.a, a);
    Copy(m_coeff.b, b);
}

Real
FMultiGrid::solve (MultiFab& phi,
		   MultiFab& rhs,
		   Real rel_tol, Real abs_tol,
		   int always_use_bnorm,
		   int need_grad_phi,
		   int verbose)
{
    PArray<MultiFab> phi_p;
    PArray<MultiFab> rhs_p;
    Copy(phi_p, phi);
    Copy(rhs_p, rhs);
    return solve(phi_p, rhs_p, rel_tol, abs_tol,
		 always_use_bnorm, need_grad_phi, verbose);
}

Real
FMultiGrid::solve (PArray<MultiFab>& phi,
		   PArray<MultiFab>& rhs,
		   Real rel_tol, Real abs_tol,
		   int always_use_bnorm ,
		   int need_grad_phi,
		   int verbose)
{
    BL_ASSERT(  m_bc.initilized || m_bndry != 0);
    BL_ASSERT(!(m_bc.initilized && m_bndry != 0));
    BL_ASSERT(m_coeff.eq_type != invalid_eq);
    BL_ASSERT(m_mgt_solver == 0);

    MultiFab* phi_p[m_nlevels];
    MultiFab* rhs_p[m_nlevels];
    for (int ilev=0; ilev < m_nlevels; ++ilev) 
    {
	phi_p[ilev] = &phi[ilev];
	rhs_p[ilev] = &rhs[ilev];
    }    

    init_mgt_solver(phi);

    Real final_resnorm;
    m_mgt_solver->solve(phi_p, rhs_p, *m_bndry, rel_tol, abs_tol, 
			always_use_bnorm, final_resnorm, need_grad_phi);
    return final_resnorm;
}

void
FMultiGrid::get_fluxes (PArray<MultiFab>& grad_phi, int ilev)
{
    BL_ASSERT(ilev < m_nlevels);

    const Real* dx = m_geom[ilev].CellSize();
    m_mgt_solver->get_fluxes(ilev, grad_phi, dx);
}

void
FMultiGrid::get_fluxes (Array<PArray<MultiFab> >& grad_phi)
{
    for (int ilev = 0; ilev < m_nlevels; ++ilev) 
    {	
	const Real* dx = m_geom[ilev].CellSize();
	m_mgt_solver->get_fluxes(ilev, grad_phi[ilev], dx);
    }
}

void
FMultiGrid::get_fluxes (PArray<PArray<MultiFab> >& grad_phi)
{
    for (int ilev = 0; ilev < m_nlevels; ++ilev) 
    {
	const Real* dx = m_geom[ilev].CellSize();
	m_mgt_solver->get_fluxes(ilev, grad_phi[ilev], dx);
    }
}

void 
FMultiGrid::compute_residual (MultiFab & phi,
			      MultiFab & rhs,
			      MultiFab & res)
{
    PArray<MultiFab> phi_p;
    PArray<MultiFab> rhs_p;
    PArray<MultiFab> res_p;
    Copy(phi_p, phi);
    Copy(rhs_p, rhs);
    Copy(res_p, res);
    compute_residual(phi_p, rhs_p, res_p);
}

void
FMultiGrid::compute_residual (PArray<MultiFab> & phi,
			      PArray<MultiFab> & rhs,
			      PArray<MultiFab> & res)
{
    BL_ASSERT(  m_bc.initilized || m_bndry != 0);
    BL_ASSERT(!(m_bc.initilized && m_bndry != 0));
    BL_ASSERT(m_coeff.eq_type != invalid_eq);
    BL_ASSERT(m_mgt_solver == 0);

    MultiFab* phi_p[m_nlevels];
    MultiFab* rhs_p[m_nlevels];
    MultiFab* res_p[m_nlevels];
    for (int ilev=0; ilev < m_nlevels; ++ilev) 
    {
	phi_p[ilev] = &phi[ilev];
	rhs_p[ilev] = &rhs[ilev];
	res_p[ilev] = &res[ilev];
    }    

    init_mgt_solver(phi);

    m_mgt_solver->compute_residual(phi_p, rhs_p, res_p, *m_bndry);
}

void 
FMultiGrid::applyop (MultiFab & phi,
		     MultiFab & res)
{
    PArray<MultiFab> phi_p;
    PArray<MultiFab> res_p;
    Copy(phi_p, phi);
    Copy(res_p, res);
    applyop(phi_p, res_p);
}

void
FMultiGrid::applyop (PArray<MultiFab> & phi,
		     PArray<MultiFab> & res)
{
    BL_ASSERT(  m_bc.initilized || m_bndry != 0);
    BL_ASSERT(!(m_bc.initilized && m_bndry != 0));
    BL_ASSERT(m_coeff.eq_type != invalid_eq);
    BL_ASSERT(m_mgt_solver == 0);
    BL_ASSERT(m_bndry == 0);

    MultiFab* phi_p[m_nlevels];
    MultiFab* res_p[m_nlevels];
    for (int ilev=0; ilev < m_nlevels; ++ilev) 
    {
	phi_p[ilev] = &phi[ilev];
	res_p[ilev] = &res[ilev];
    }    

    init_mgt_solver(phi);

    m_mgt_solver->applyop(phi_p, res_p, *m_bndry);
}

void
FMultiGrid::Copy (PArray<MultiFab>& dst, const MultiFab& src)
{
    dst.resize(1);
    dst.set(0, &src);
}

void
FMultiGrid::Copy (PArray<MultiFab>& dst, PArray<MultiFab>& src)
{
    int nlevels = src.size();
    dst.resize(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
	dst.set(ilev, &src[ilev]);
    }
}

void 
FMultiGrid::Copy (Array<PArray<MultiFab> >& dst, PArray<MultiFab>& src)
{
    int ndim = src.size();
    dst.resize(1);
    dst[0].resize(ndim);
    for(int idim = 0; idim < ndim; ++idim) {
	dst[0].set(idim, &src[idim]);
    }
}

void 
FMultiGrid::Copy (Array<PArray<MultiFab> >& dst, Array<PArray<MultiFab> >& src)
{
    int nlevels = src.size();
    int ndim = src[0].size();
    dst.resize(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
	dst[ilev].resize(ndim);
	for(int idim = 0; idim < ndim; ++idim) {
	    dst[ilev].set(idim, &src[ilev][idim]);
	}
    }
}

void
FMultiGrid::Boundary::set_bndry_values (MacBndry& bndry, IntVect crse_ratio)
{
    // The values of phys_bc & ref_ratio do not matter
    // because we are not going to use those parts of MacBndry.
    IntVect ref_ratio = IntVect::TheZeroVector();
    Array<int> lo_bc(BL_SPACEDIM, 0);
    Array<int> hi_bc(BL_SPACEDIM, 0);
    BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());

    if (crse_phi == 0 && phi == 0) 
    {
	bndry.setHomogValues(phys_bc, ref_ratio);
    }
    else if (crse_phi == 0 && phi != 0)
    {
	bndry.setBndryValues(*phi, 0, 0, phi->nComp(), phys_bc); 
    }
    else if (crse_phi != 0 && phi != 0) 
    {
	BL_ASSERT(crse_ratio != IntVect::TheZeroVector());

	const int ncomp      = phi->nComp();
	const int in_rad     = 0;
	const int out_rad    = 1;
	const int extent_rad = 2;

	BoxArray crse_boxes(phi->boxArray());
	crse_boxes.coarsen(crse_ratio);

	BndryRegister crse_br(crse_boxes, in_rad, out_rad, extent_rad, ncomp);
	crse_br.copyFrom(*crse_phi, crse_phi->nGrow(), 0, 0, ncomp);

	bndry.setBndryValues(crse_br, 0, *phi, 0, 0, ncomp, crse_ratio, phys_bc);
    }
    else
    {
	BoxLib::Abort("FMultiGrid::Boundary::build_bndry: How did we get here?");
    }
}

void
FMultiGrid::ABecCoeff::set_coeffs (MGT_Solver & mgt_solver, FMultiGrid& fmg)
{
    BL_ASSERT( fmg.m_baselevel >= 0 );
    BL_ASSERT( fmg.m_baselevel == 0 || fmg.m_crse_ratio != IntVect::TheZeroVector() );

    Array< Array<Real> > xa(fmg.m_nlevels);
    Array< Array<Real> > xb(fmg.m_nlevels);

    for (int lev=0; lev < fmg.m_nlevels; ++lev) {
	xa[lev].resize(BL_SPACEDIM);
	xb[lev].resize(BL_SPACEDIM);
	if (lev + fmg.m_baselevel == 0) {
	    // For level 0, the boundary lives exactly on the faces
	    for (int n=0; n<BL_SPACEDIM; n++) {
		xa[lev][n] = 0.0;
		xb[lev][n] = 0.0;
	    }
	} else if (lev == 0) {
	    const Real* dx = fmg.m_geom[0].CellSize();
	    for (int n=0; n<BL_SPACEDIM; n++) {
		xa[lev][n] = 0.5 * fmg.m_crse_ratio[n] * dx[n];
		xb[lev][n] = 0.5 * fmg.m_crse_ratio[n] * dx[n];
	    }	    
	} else {
	    const Real* dx_crse = fmg.m_geom[lev-1].CellSize();
	    for (int n=0; n<BL_SPACEDIM; n++) {
		xa[lev][n] = 0.5 * dx_crse[n];
		xb[lev][n] = 0.5 * dx_crse[n];
	    }
	}
    }

    switch (eq_type) {
        case const_gravity_eq:
	{
	    mgt_solver.set_const_gravity_coeffs(xa, xb);
	    break;
	}
        case (gravity_eq):
	{
	    BL_ASSERT(coeffs_set);
	    mgt_solver.set_gravity_coefficients(b, xa, xb);
	    break;
	}
        case (macproj_eq):
	{
	    BL_ASSERT(coeffs_set);
	    mgt_solver.set_mac_coefficients(b, xa, xb);
	    break;
	}
        case (general_eq):
	{
	    BL_ASSERT(scalars_set && coeffs_set);
	    mgt_solver.set_abeclap_coeffs(alpha, a, beta, b, xa, xb);
	    break;
	}
        default:
	{
	    BoxLib::Abort("FMultiGrid::ABecCoeff::set_coeffs: How did we get here?");
	}
    }
}

void
FMultiGrid::init_mgt_solver (PArray<MultiFab>& phi)
{
    BL_ASSERT(  m_bc.initilized || m_bndry != 0);
    BL_ASSERT(!(m_bc.initilized && m_bndry != 0));
    BL_ASSERT(m_coeff.eq_type != invalid_eq);
    BL_ASSERT(m_mgt_solver == 0);

    int ncomp = phi[0].nComp();

    std::vector<DistributionMapping> dmap(m_nlevels);
    std::vector<BoxArray> ba(m_nlevels);

    for (int ilev = 0; ilev < m_nlevels; ++ilev) 
    {
	dmap[ilev] = phi[ilev].DistributionMap();
	ba  [ilev] = phi[ilev].boxArray();
    }

    bool nodal = false;
    int  nc = 0;
    Array<int> mg_bc;
    if (m_bc.initilized) 
    {
	mg_bc = Array<int>(m_bc.get_mg_bc(), 2*BL_SPACEDIM);
    } 
    else 
    {
	mg_bc.resize(2*BL_SPACEDIM);
	
	for ( int i = 0; i < BL_SPACEDIM; ++i )
        { 
            if ( m_geom[0].isPeriodic(i) )
            {
                mg_bc[i*2 + 0] = 0;
                mg_bc[i*2 + 1] = 0;
            }
            else
            {
                mg_bc[i*2 + 0] = m_bndry->phys_bc_lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
                mg_bc[i*2 + 1] = m_bndry->phys_bc_hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
            }
        }
    }

    RAII_mgt_solver.resize(1, PArrayManage);
    RAII_mgt_solver.set(0, new MGT_Solver (m_geom, mg_bc.dataPtr() , ba, dmap, nodal,
					   m_stencil, nodal, nc, ncomp, m_verbose));
    m_mgt_solver = &(RAII_mgt_solver[0]);

    if (m_maxorder > 0) m_mgt_solver->set_maxorder(m_maxorder);

    if (m_bndry == 0) {
	RAII_bndry.resize(1, PArrayManage);
	RAII_bndry.set(0, new MacBndry (ba[0], ncomp, m_geom[0]));
	m_bndry = &(RAII_bndry[0]);

	m_bc.set_bndry_values(*m_bndry, m_crse_ratio);
    }

    m_coeff.set_coeffs(*m_mgt_solver, *this);
}
