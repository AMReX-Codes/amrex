
#include "hg_multi.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGRES      hgres_
#  define   FORT_HGRES_TERRAIN      hgres_terrain_
#  define   FORT_HGRESU     hgresu_
#  define   FORT_HGRLX      hgrlx_
#  define   FORT_HGRLX_TERRAIN      hgrlx_terrain_
#  define   FORT_HGRLXU     hgrlxu_
#  define   FORT_HGRLXL     hgrlxl_
#  define   FORT_HGRLNF     hgrlnf_
#  define   FORT_HGRLNF_TERRAIN     hgrlnf_terrain
#  define   FORT_HGRLNB     hgrlnb_
#  define   FORT_HGCG       hgcg_
#  define   FORT_HGCG1      hgcg1_
#  define   FORT_HGCG2      hgcg2_
#  define   FORT_HGIP       hgip_
#else
#  define   FORT_HGRES      HGRES
#  define   FORT_HGRES_TERRAIN      HGRES_TERRAIN
#  define   FORT_HGRESU     HGRESU
#  define   FORT_HGRLX      HGRLX
#  define   FORT_HGRLX_TERRAIN      HGRLX_TERRAIN
#  define   FORT_HGRLXU     HGRLXU
#  define   FORT_HGRLXL     HGRLXL
#  define   FORT_HGRLNF     HGRLNF
#  define   FORT_HGRLNF_TERRAIN     HGRLNF_TERRAIN
#  define   FORT_HGRLNB     HGRLNB
#  define   FORT_HGCG       HGCG
#  define   FORT_HGCG1      HGCG1
#  define   FORT_HGCG2      HGCG2
#  define   FORT_HGIP       HGIP
#endif

extern "C" 
{
    
#if (BL_SPACEDIM == 1)
#error not relevant
#else
    void FORT_HGRES_TERRAIN(Real*, intS, const Real*, intS, const Real*, intS, Real*, intS, Real*, intS, intS);
    void FORT_HGRLX_TERRAIN(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
    void FORT_HGRLNF_TERRAIN(Real*, intS, Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS, intS, const int&, const int&);
    void FORT_HGRESU(Real*, intS, const Real*, const Real*, const Real*, Real*, intS);
#if (defined (HG_CROSS_STENCIL) && BL_SPACEDIM==3)
    void FORT_HGRES(Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS);
    void FORT_HGRLX(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
    void FORT_HGRLXU(Real*, Real*, Real*, Real*, intS, Real*, intS);
    void FORT_HGRLXL(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS, intS, const int&);
    void FORT_HGRLNF(Real*, intS, Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS, intS, const int&, const int&);
#    else
    void FORT_HGRES(Real*, intS, const Real*, intS, const Real*, intS, CRealPS, intS, intS, CRealRS, const int&, const int&);
    void FORT_HGRLX(Real*, intS, const Real*, intS, CRealPS, intS, const Real*, intS, intS, CRealRS, const int&, const int&);
    void FORT_HGRLXL(Real*, intS, Real*, intS, RealPS, intS, Real*, intS, intS, intS, RealRS, const int&, const int&, const int&);
    void FORT_HGRLNF(Real*, intS, Real*, intS, Real*, intS, RealPS, intS, Real*, intS, intS, intS, RealRS, const int&, const int&, const int&, const int&);
#    endif
    void FORT_HGRLNB(Real*, intS, Real*, intS, intS, const int&, const int&);
    
    void FORT_HGCG1(Real*, Real*, Real*, Real*, Real*, Real*, Real*, intS, const Real&, Real&);
    void FORT_HGCG2(Real*, Real*, intS, const Real&);
    void FORT_HGIP(Real*, Real*, Real*, intS, Real&);

#endif
}

void holy_grail_amr_multigrid::level_residual(MultiFab& r, MultiFab& s, MultiFab& d, int mglev, bool iclear)
{
    assert(r.boxArray() == s.boxArray());
    assert(r.boxArray() == d.boxArray());
    assert(mglev >= 0);
    fill_borders(d, lev_interface[mglev], mg_boundary, -1, m_hg_terrain);
    
    if(m_hg_terrain)
    {
	
	for (MultiFabIterator r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    DependentMultiFabIterator s_dmfi(r_mfi, s);
	    DependentMultiFabIterator d_dmfi(r_mfi, d);
	    DependentMultiFabIterator c_dmfi(r_mfi, cen[mglev]);
	    DependentMultiFabIterator sg_dmfi(r_mfi, sigma[mglev]);
	    const Box& rbox = r_mfi->box();
	    const Box& sbox = s_dmfi->box();
	    const Box& dbox = d_dmfi->box();
	    const Box& cenbox = c_dmfi->box();
	    const Box& sigbox = sg_dmfi->box();
	    const Box& freg = lev_interface[mglev].part_fine(r_mfi.index());
	    FORT_HGRES_TERRAIN(r_mfi->dataPtr(), DIMLIST(rbox),
		s_dmfi->dataPtr(), DIMLIST(sbox),
		d_dmfi->dataPtr(), DIMLIST(dbox),
		sg_dmfi->dataPtr(), DIMLIST(sigbox),
		c_dmfi->dataPtr(), DIMLIST(cenbox),
		DIMLIST(freg));
	}
	
	if (iclear) 
	{
	    clear_part_interface(r, lev_interface[mglev]);
	}
    }
    else if (m_hg_cross_stencil)
    {
	for (MultiFabIterator r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    DependentMultiFabIterator s_dmfi(r_mfi, s);
	    DependentMultiFabIterator d_dmfi(r_mfi, d);
	    DependentMultiFabIterator sn_dmfi(r_mfi, sigma_node[mglev]);
	    DependentMultiFabIterator m_dmfi(r_mfi, mask[mglev]);
	    const Box& rbox = r_mfi->box();
	    const Box& freg = lev_interface[mglev].part_fine(r_mfi.index());
	    FORT_HGRESU(r_mfi->dataPtr(), DIMLIST(rbox),
		s_dmfi->dataPtr(), d_dmfi->dataPtr(), sn_dmfi->dataPtr(), m_dmfi->dataPtr(), DIMLIST(freg));
	}
	
    }
    else
    {
#if BL_SPACEDIM != 3
	Real hx = h[mglev][0];
	Real hy = h[mglev][1];
#  if (BL_SPACEDIM == 3)
	Real hz = h[mglev][2];
#  endif
	
	for (MultiFabIterator r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    DependentMultiFabIterator s_dmfi(r_mfi, s);
	    DependentMultiFabIterator d_dmfi(r_mfi, s);
	    const Box& rbox = r_mfi->box();
	    const Box& sbox = s_dmfi->box();
	    const Box& dbox = d_dmfi->box();
	    const Box& freg = lev_interface[mglev].part_fine(r_mfi.index());
	    DependentMultiFabIterator sg_dmfi(r_mfi, sigma[mglev]);
	    DependentMultiFabIterator s0_dmfi(r_mfi, sigma_nd[0][mglev]);
	    DependentMultiFabIterator s1_dmfi(r_mfi, sigma_nd[1][mglev]);
	    DependentMultiFabIterator s2_dmfi(r_mfi, sigma_nd[2][mglev]);
	    // this branch is the only one that can be reached here
	    const Box& sigbox = sg_dmfi->box();
#if	(BL_SPACEDIM==2)
	    FORT_HGRES(r_mfi->dataPtr(), DIMLIST(rbox),
		s_dmfi->dataPtr(), DIMLIST(sbox),
		d_dmfi->dataPtr(), DIMLIST(dbox),
		sigma_nd[0][mglev][igrid].dataPtr(),
		s1_dmfi->dataPtr(), DIMLIST(sigbox),
		DIMLIST(freg), hx, hy,
		IsRZ(), mg_domain[mglev].bigEnd(0) + 1
		);
#else
	    FORT_HGRES(r_mfi->dataPtr(), DIMLIST(rbox),
		s_dmfi->dataPtr(), DIMLIST(sbox),
		d_dmfi->dataPtr(), DIMLIST(dbox),
		s0_dmfi->dataPtr(),
		s1_dmfi->dataPtr(),
		s2_dmfi->dataPtr(), DIMLIST(sigbox),
		DIMLIST(freg), hx, hy, hz
		);
#endif
	}
	
	if (iclear) 
	{
	    clear_part_interface(r, lev_interface[mglev]);
	}
#endif
    }
}

void holy_grail_amr_multigrid::relax(int mglev, int i1, bool is_zero)
{
    Real hx = h[mglev][0];
    Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
    Real hz = h[mglev][2];
#endif
    
    Box tdom = mg_domain[mglev];
    tdom.convert(IntVect::TheNodeVector());
    
    for (int icount = 0; icount < i1; icount++) 
    {
	
	if (smoother_mode == 0 || smoother_mode == 1 || line_solve_dim == -1) 
	{
	    
	    if ( is_zero == false )
		fill_borders(corr[mglev], lev_interface[mglev], mg_boundary, -1, m_hg_terrain);
	    else
		is_zero = false;
	    for (MultiFabIterator r_mfi(resid[mglev]); r_mfi.isValid(); ++r_mfi)
	    {
		DependentMultiFabIterator c_dmfi(r_mfi, corr[mglev]);
		DependentMultiFabIterator cn_dmfi(r_mfi, cen[mglev]);
		DependentMultiFabIterator sn_dmfi(r_mfi, sigma_node[mglev]);
		DependentMultiFabIterator m_dmfi(r_mfi, mask[mglev]);
		const Box& sbox = r_mfi->box();
		const Box& freg = lev_interface[mglev].part_fine(r_mfi.index());
		if (line_solve_dim == -1) 
		{
		    // Gauss-Seidel section:
		    if(m_hg_terrain)
		    {
			DependentMultiFabIterator sg_dmfi(r_mfi, sigma[mglev]);
			const Box& fbox = c_dmfi->box();
			const Box& cenbox = cn_dmfi->box();
			const Box& sigbox = sg_dmfi->box();
			FORT_HGRLX(c_dmfi->dataPtr(), DIMLIST(fbox),
			    r_mfi->dataPtr(), DIMLIST(sbox),
			    sg_dmfi->dataPtr(), DIMLIST(sigbox),
			    cn_dmfi->dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg));
		    } 
		    else if (m_hg_cross_stencil)
		    {
			FORT_HGRLXU(c_dmfi->dataPtr(),
			    r_mfi->dataPtr(),
			    sn_dmfi->dataPtr(),
			    cn_dmfi->dataPtr(), DIMLIST(sbox),
			    m_dmfi->dataPtr(),
			    DIMLIST(freg));
		    }
		    else
		    {
#if	BL_SPACEDIM != 3
			const Box& fbox = c_dmfi->box();
			const Box& cenbox = cn_dmfi->box();
			const Box& sigbox = sigma[mglev][igrid].box();
#if (BL_SPACEDIM==2)
			FORT_HGRLX(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			    resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			    sigma_nd[0][mglev][igrid].dataPtr(),
			    sigma_nd[1][mglev][igrid].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), hx, hy,
			    IsRZ(), mg_domain[mglev].bigEnd(0) + 1
			    );
#else
			FORT_HGRLX(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			    resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			    sigma_nd[0][mglev][igrid].dataPtr(),
			    sigma_nd[1][mglev][igrid].dataPtr(),
			    sigma_nd[2][mglev][igrid].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), hx, hy, hz
			    );
#endif
#endif
		    }
		}
		else 
		{
		    // Grid-by-grid line solve section:
		    if(m_hg_terrain)
			BoxLib::Error("Terrain line solves not implemented");
		    const Box& fbox = c_dmfi->box();
		    const Box& cenbox = cn_dmfi->box();
		    if(m_hg_cross_stencil)
		    {
			const Box& sigbox = sn_dmfi->box();
			FORT_HGRLXL(c_dmfi->dataPtr(), DIMLIST(fbox),
			    r_mfi->dataPtr(), DIMLIST(sbox),
			    sn_dmfi->dataPtr(), DIMLIST(sigbox),
			    cn_dmfi->dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom), line_solve_dim);
		    }
		    else
		    {
#if BL_SPACEDIM!=3
			const Box& sigbox = sigma[mglev][igrid].box();
#if (BL_SPACEDIM==2)
			FORT_HGRLXL(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			    resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			    sigma_nd[0][mglev][igrid].dataPtr(),
			    sigma_nd[1][mglev][igrid].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom), hx, hy,
			    IsRZ(), mg_domain[mglev].bigEnd(0) + 1, line_solve_dim
			    );
#else
			FORT_HGRLXL(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			    resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			    sigma_nd[0][mglev][igrid].dataPtr(),
			    sigma_nd[1][mglev][igrid].dataPtr(),
			    sigma_nd[2][mglev][igrid].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom), hx, hy, hz
			    );
#endif
#  endif
		    }
		}
      }
      sync_borders(corr[mglev], lev_interface[mglev], mg_boundary);
    }
    else 
    {
	// Full-level line solve section:
	if (line_order.length() == 0) 
	{
	    build_line_order(line_solve_dim);
	}
	int lev = lev_min;
	while (ml_index[lev] < mglev)
	    lev++;
	
	for (int ipass = 0; ipass <= 1; ipass++) 
	{
	    if (is_zero == false)
		fill_borders(corr[mglev], lev_interface[mglev], mg_boundary, -1, m_hg_terrain);
	    else
		is_zero = false;
	    
	    // Forward solve:
	    // PARALLEL TODO
	    for (int i = 0; i < mg_mesh[mglev].length(); i++) 
	    {
		
		// Do grids in order along line_solve_dim:
		int igrid = line_order[lev][i];
		const Box& fbox = corr[mglev][igrid].box();
		const Box& sbox = resid[mglev][igrid].box();
		const Box& wbox = work[mglev][igrid].box();
		const Box& cenbox = cen[mglev][igrid].box();
		const Box& freg = corr[mglev].box(igrid);
		if(m_hg_terrain)
		{
		    const Box& sigbox = sigma[mglev][igrid].box();
		    FORT_HGRLNF(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			sigma[mglev][igrid].dataPtr(), DIMLIST(sigbox),
			cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			DIMLIST(freg), DIMLIST(tdom), line_solve_dim, ipass);
		}
		else if (m_hg_cross_stencil)
		{
		    const Box& sigbox = sigma_node[mglev][igrid].box();
		    FORT_HGRLNF(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			sigma_node[mglev][igrid].dataPtr(), DIMLIST(sigbox),
			cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			DIMLIST(freg), DIMLIST(tdom), line_solve_dim, ipass);
		}
		else
		{
#  if	BL_SPACEDIM!=3
		    const Box& sigbox = sigma[mglev][igrid].box();
#if BL_SPACEDIM==2
		    FORT_HGRLNF(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			sigma_nd[0][mglev][igrid].dataPtr(),
			sigma_nd[1][mglev][igrid].dataPtr(), DIMLIST(sigbox),
			cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			DIMLIST(freg), DIMLIST(tdom), hx, hy,
			IsRZ(), mg_domain[mglev].bigEnd(0) + 1,
			line_solve_dim, ipass
			);
#else
		    FORT_HGRLNF(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			sigma_nd[0][mglev][igrid].dataPtr(),
			sigma_nd[1][mglev][igrid].dataPtr(),
			sigma_nd[2][mglev][igrid].dataPtr(), DIMLIST(sigbox),
			cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			DIMLIST(freg), DIMLIST(tdom), hx, hy, hz
			);
#endif
#  endif
		}		
		// Copy work arrays to following grids:
		for (ListIterator<int> j(line_after[lev][igrid]); j; j++) 
		{
		    Box b = (freg & corr[mglev].box(j()));
		    internal_copy(corr[mglev], j(), igrid, b);
		    internal_copy(work[mglev], j(), igrid, b);
		}
	    }
	    
	    // Back substitution:
	    // PARALLEL TODO
	    for (int i = mg_mesh[mglev].length() - 1; i >= 0; i--) 
	    {
		
		// Do grids in reverse order along line_solve_dim:
		int igrid = line_order[lev][i];
		const Box& freg = corr[mglev].box(igrid);
		
		// Copy solution array from following grids:
		for (ListIterator<int> j(line_after[lev][igrid]); j; j++) 
		{
		    Box b = (freg & corr[mglev].box(j()));
		    internal_copy(corr[mglev], igrid, j(), b);
		}
		
		const Box& fbox = corr[mglev][igrid].box();
		const Box& wbox = work[mglev][igrid].box();
		FORT_HGRLNB(corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
		    work[mglev][igrid].dataPtr(), DIMLIST(wbox),
		    DIMLIST(freg), line_solve_dim, ipass);
	    }
	}
    }
  }
}

void holy_grail_amr_multigrid::build_line_order(int lsd)
{
    line_order.resize(lev_max + 1);
    line_after.resize(lev_max + 1);
    
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	int mglev = ml_index[lev], ngrids = mg_mesh[mglev].length();
	
	line_order[lev].resize(ngrids);
	line_after[lev].resize(ngrids);
	
	for (int igrid = 0; igrid < ngrids; igrid++) 
	{
	    line_order[lev].set(igrid, igrid);
	    
	    // bubble sort, replace with something faster if necessary:
	    for (int i = igrid; i > 0; i--) 
	    {
		if (ml_mesh[lev][line_order[lev][i]].smallEnd(lsd) <
		    ml_mesh[lev][line_order[lev][i-1]].smallEnd(lsd)) 
		{
		    int tmp              = line_order[lev][i-1];
		    line_order[lev][i-1] = line_order[lev][i];
		    line_order[lev][i]   = tmp;
		}
		else 
		{
		    break;
		}
	    }
	    
	    for (int i = 0; i < ngrids; i++) 
	    {
		if (bdryLo(ml_mesh[lev][i], lsd).intersects
		    (bdryHi(ml_mesh[lev][igrid], lsd))) 
		{
		    line_after[lev][igrid].append(i);
		}
	    }
	}
    }
}

void holy_grail_amr_multigrid::cgsolve(int mglev)
{
    assert(mglev == 0);
    
    MultiFab& r = cgwork[0];
    MultiFab& p = cgwork[1];
    MultiFab& z = cgwork[2];
    MultiFab& x = cgwork[3];
    MultiFab& w = cgwork[4];
    MultiFab& c = cgwork[5];
    MultiFab& zero_array = cgwork[6];
    MultiFab& ipmask = cgwork[7];
    
    
    Real alpha, rho;
    
    // x (corr[0]) should be all 0.0 at this point
    for(MultiFabIterator r_mfi(r); r_mfi.isValid(); ++r_mfi)
    {
	DependentMultiFabIterator r_dmfi(r_mfi, resid[mglev]);
	r_mfi->copy(r_dmfi());
	r_mfi->negate();
    }
    //r.copy(resid[mglev]);
    //r.negate();
    
    if (singular) 
    {
	// singular systems are very sensitive to solvability
	w.setVal(1.0);
	alpha = inner_product(r, w) / mg_domain[mglev].volume();
	r.plus(-alpha, 0);
    }
    
    rho = 0.0;
    for(MultiFabIterator r_mfi(r); r_mfi.isValid(); ++r_mfi)
    {
	DependentMultiFabIterator z_dmfi(r_mfi, z);
	DependentMultiFabIterator c_dmfi(r_mfi, c);
	DependentMultiFabIterator i_dmfi(r_mfi, ipmask);
	DependentMultiFabIterator p_dmfi(r_mfi, p);
	z_dmfi->copy(r_mfi());
	z_dmfi->mult(c_dmfi());
	const Box& reg = p_dmfi->box();
	FORT_HGIP(z_dmfi->dataPtr(), r_mfi->dataPtr(), i_dmfi->dataPtr(), DIMLIST(reg), rho);
	p_dmfi->copy(z_dmfi());
    }
    ParallelDescriptor::ReduceRealSum(rho);
    Real tol = 1.e-3 * rho;
    
    int i = 0;
    while (tol > 0.0) 
    {
	i++;
	if (i > 250 && pcode >= 2)
	    cout << "Conjugate-gradient iteration failed to converge" << endl;
	Real rho_old = rho;
	// safe to set the clear flag to 0 here---bogus values make it
	// into r but are cleared from z by the mask in c
	level_residual(w, zero_array, p, 0, false);
	alpha = 0.0;
	for(MultiFabIterator p_mfi(p); p_mfi.isValid(); ++p_mfi)
	{
	    DependentMultiFabIterator w_dmfi(p_mfi, w);
	    DependentMultiFabIterator i_dmfi(p_mfi, ipmask);
	    const Box& reg = p_mfi->box();
	    FORT_HGIP(p_mfi->dataPtr(), w_dmfi->dataPtr(), i_dmfi->dataPtr(), DIMLIST(reg), alpha);
	    
	}
	ParallelDescriptor::ReduceRealSum(alpha);
	alpha = rho / alpha;
	rho = 0.0;
	for(MultiFabIterator r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    DependentMultiFabIterator p_dmfi(r_mfi, p);
	    DependentMultiFabIterator z_dmfi(r_mfi, z);
	    DependentMultiFabIterator x_dmfi(r_mfi, x);
	    DependentMultiFabIterator w_dmfi(r_mfi, w);
	    DependentMultiFabIterator c_dmfi(r_mfi, c);
	    DependentMultiFabIterator i_dmfi(r_mfi, ipmask);
	    const Box& reg = p_dmfi->box();
	    FORT_HGCG1(r_mfi->dataPtr(), p_dmfi->dataPtr(), z_dmfi->dataPtr(), x_dmfi->dataPtr(),
		w_dmfi->dataPtr(), c_dmfi->dataPtr(), i_dmfi->dataPtr(), DIMLIST(reg), alpha, rho);
	}
	ParallelDescriptor::ReduceRealSum(rho);
	if (pcode >= 3)
	    cout << i << " " << rho << endl;
	if (rho <= tol || i > 250)
	    break;
	alpha = rho / rho_old;
	for(MultiFabIterator p_mfi(p); p_mfi.isValid(); ++p_mfi)
	{
	    DependentMultiFabIterator z_dmfi(p_mfi, z);
	    const Box& reg = p_mfi->box();
	    FORT_HGCG2(p_mfi->dataPtr(), z_dmfi->dataPtr(), DIMLIST(reg), alpha);
	}
    }
    
    if (pcode >= 2)
	cout << i << " iterations required for conjugate-gradient" << endl;
}
