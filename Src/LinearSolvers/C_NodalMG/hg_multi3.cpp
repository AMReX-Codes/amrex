
#include "hg_multi.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define   FORT_HGRES		hgres_
#define   FORT_HGRES_TERRAIN	hgres_terrain_
#define   FORT_HGRES_FULL	hgres_full_
#define   FORT_HGRESU		hgresu_
#define   FORT_HGRESUR		hgresur_
#define   FORT_HGRLX		hgrlx_
#define   FORT_HGRLX_TERRAIN	hgrlx_terrain_
#define   FORT_HGRLX_FULL	hgrlx_full_
#define   FORT_HGRLXU		hgrlxu_
#define   FORT_HGRLXUR		hgrlxur_
#define   FORT_HGRLXL		hgrlxl_
#define   FORT_HGRLXL_FULL	hgrlxl_full_
#define   FORT_HGRLNF		hgrlnf_
#define   FORT_HGRLNF_TERRAIN	hgrlnf_terrain_
#define   FORT_HGRLNF_FULL	hgrlnf_full_
#define   FORT_HGRLNB		hgrlnb_
#define   FORT_HGCG		hgcg_
#define   FORT_HGCG1		hgcg1_
#define   FORT_HGCG2		hgcg2_
#define   FORT_HGIP		hgip_
#elif defined( BL_FORT_USE_UPPERCASE )
#define   FORT_HGRES		HGRES
#define   FORT_HGRES_TERRAIN    HGRES_TERRAIN
#define   FORT_HGRES_FULL	HGRES_FULL
#define   FORT_HGRESU		HGRESU
#define   FORT_HGRESUR		HGRESUR
#define   FORT_HGRLX		HGRLX
#define   FORT_HGRLX_TERRAIN    HGRLX_TERRAIN
#define   FORT_HGRLX_FULL	HGRLX_FULL
#define   FORT_HGRLXU		HGRLXU
#define   FORT_HGRLXUR		HGRLXUR
#define   FORT_HGRLXL		HGRLXL
#define   FORT_HGRLXL_FULL	HGRLXL_FULL
#define   FORT_HGRLNF		HGRLNF
#define   FORT_HGRLNF_TERRAIN   HGRLNF_TERRAIN
#define   FORT_HGRLNF_FULL	HGRLNF_FULL
#define   FORT_HGRLNB		HGRLNB
#define   FORT_HGCG		HGCG
#define   FORT_HGCG1		HGCG1
#define   FORT_HGCG2		HGCG2
#define   FORT_HGIP		HGIP
#elif defined( BL_FORT_USE_LOWERCASE )
#define   FORT_HGRES		hgres
#define   FORT_HGRES_TERRAIN    hgres_terrain
#define   FORT_HGRES_FULL	hgres_full
#define   FORT_HGRESU		hgresu
#define   FORT_HGRESUR		hgresur
#define   FORT_HGRLX		hgrlx
#define   FORT_HGRLX_TERRAIN    hgrlx_terrain
#define   FORT_HGRLX_FULL	hgrlx_full
#define   FORT_HGRLXU		hgrlxu
#define   FORT_HGRLXUR		hgrlxur
#define   FORT_HGRLXL		hgrlxl
#define   FORT_HGRLXL_FULL	hgrlxl_full
#define   FORT_HGRLNF		hgrlnf
#define   FORT_HGRLNF_TERRAIN   hgrlnf_terrain
#define   FORT_HGRLNF_FULL	hgrlnf_full
#define   FORT_HGRLNB		hgrlnb
#define   FORT_HGCG		hgcg
#define   FORT_HGCG1		hgcg1
#define   FORT_HGCG2		hgcg2
#define   FORT_HGIP		hgip
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{

#if (BL_SPACEDIM == 1)
#error not relevant
#endif
    void FORT_HGRES_TERRAIN (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRLX_TERRAIN (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRLNF_TERRAIN(Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS, intS,
			     const int*, const int*);
#if (BL_SPACEDIM==3)
    void FORT_HGRES_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRESU        (Real*, intS, const Real*, const Real*,
			     const Real*, const Real*, intS, const int*);
    void FORT_HGRESUR       (Real*, intS, const Real*, const Real*,
			     const Real*, intS, const int*);
    void FORT_HGRLX         (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS, intS);
    void FORT_HGRLX_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS, intS);
    void FORT_HGRLXU        (Real*, const Real*, const Real*,
			     const Real*, intS, const Real*, intS, const int*);
    void FORT_HGRLXUR        (Real*, const Real*, const Real*,
			     const Real*, intS, intS, const int*);
    void FORT_HGRLXL        (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS, intS,
			     intS, const int*);
    void FORT_HGRLNF        (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     const Real*, intS,
			     const Real*, intS, intS, intS,
			     const int*, const int*);
    void FORT_HGRLNF_FULL   (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS, intS,
			     const int*, const int*);
#else
    void FORT_HGRES         (Real*, intS, const Real*, intS,
			     const Real*, intS, CRealPS, intS, intS,
			     CRealPS, const int* );
    void FORT_HGRES_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRESU        (Real*, intS, const Real*, const Real*,
			     const Real*, const Real*, intS, const int*);
    void FORT_HGRESUR        (Real*, intS, const Real*, const Real*,
			      const Real*, intS, const int*);
    void FORT_HGRLX         (Real*, intS, const Real*, intS,
			     CRealPS, intS, const Real*, intS, intS,
			     CRealPS, const int*, const int*);
    void FORT_HGRLX_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRLXU        (Real*, const Real*, const Real*,
			     const Real*, intS, const Real*, intS, const int*);
    void FORT_HGRLXUR        (Real*, const Real*, const Real*,
			     const Real*, intS, intS, const int*);
    void FORT_HGRLXL        (Real*, intS, const Real*, intS,
			     CRealPS, intS, Real*, intS, intS, intS,
			     CRealPS, const int*, const int*, const int*);
    void FORT_HGRLXL_FULL   (Real*, intS, const Real*, intS,
			     CRealPS, intS, Real*, intS, intS, intS,
			     CRealPS, const int*, const int*, const int*);
    void FORT_HGRLNF        (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     CRealPS, intS, Real*, intS, intS, intS,
			     CRealPS, const int*,
			     const int*, const int*, const int*);
    void FORT_HGRLNF_FULL   (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     CRealPS, intS, Real*, intS, intS, intS,
			     CRealPS, const int*, const int*, const int*,
			     const int*);
#endif
    void FORT_HGRLNB        (Real*, intS, Real*, intS, intS,
			     const int*, const int*);

    void FORT_HGCG1         (Real*, const Real*, Real*, Real*,
			     const Real*, const Real*,
			     const Real*, intS, const Real&, Real&);
    void FORT_HGCG2         (Real*, const Real*, intS, const Real&);
    void FORT_HGIP          (const Real*, const Real*,
			     const Real*, intS, Real&);

}

void
holy_grail_amr_multigrid::level_residual (MultiFab& r,
                                          MultiFab& s,
                                          MultiFab& d,
                                          int       mglev,
                                          bool      iclear,
                                          int       for_fill_sync_reg)
{
    BL_ASSERT(mglev >= 0);
    BL_ASSERT(r.boxArray() == s.boxArray());
    BL_ASSERT(r.boxArray() == d.boxArray());

    HG_TEST_NORM(d, "level_residual a");
    fill_borders(d, lev_interface[mglev],
		 mg_boundary, -1, is_dense(m_stencil));
    HG_TEST_NORM(d, "level_residual a1");
    HG_TEST_NORM(s, "level_residual");
    HG_TEST_NORM(r, "level_residual");

    if ( m_stencil == terrain || m_stencil == full )
    {
	HG_TEST_NORM(sigma[mglev], "level_residual");
	for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    const Box& rbox = r[r_mfi].box();
	    const Box& sbox = s[r_mfi].box();
	    const Box& dbox = d[r_mfi].box();
	    const Box& cenbox = cen[mglev][r_mfi].box();
	    const Box& sigbox = sigma[mglev][r_mfi].box();
            Box freg = (for_fill_sync_reg > 0) ?
		BoxLib::surroundingNodes(mg_mesh[mglev][r_mfi.index()]) :
		Box(lev_interface[mglev].part_fine(r_mfi.index()));
	    if ( m_stencil == terrain )
	    {
		FORT_HGRES_TERRAIN(r[r_mfi].dataPtr(), DIMLIST(rbox),
				   s[r_mfi].dataPtr(), DIMLIST(sbox),
				   d[r_mfi].dataPtr(), DIMLIST(dbox),
				   sigma[mglev][r_mfi].dataPtr(), DIMLIST(sigbox),
				   cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
				   DIMLIST(freg));
	    }
	    else if ( m_stencil == full )
	    {
		FORT_HGRES_FULL(r[r_mfi].dataPtr(), DIMLIST(rbox),
				s[r_mfi].dataPtr(), DIMLIST(sbox),
				d[r_mfi].dataPtr(), DIMLIST(dbox),
				sigma[mglev][r_mfi].dataPtr(), DIMLIST(sigbox),
				cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
				DIMLIST(freg));
	    }
	}
	if (iclear)
	{
	    clear_part_interface(r, lev_interface[mglev]);
	}
    }
    else if (m_stencil == cross)
    {
        const int isRZ = getCoordSys();
	HG_TEST_NORM(sigma_node[mglev], "level_residual");
	HG_TEST_NORM(mask[mglev], "level_residual");
	for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    const Box& rbox = r[r_mfi].box();
            Box freg = (for_fill_sync_reg > 0) ?
		BoxLib::surroundingNodes(mg_mesh[mglev][r_mfi.index()]) :
		Box(lev_interface[mglev].part_fine(r_mfi.index()));

#if 1
	    FORT_HGRESU(r[r_mfi].dataPtr(), DIMLIST(rbox),
			s[r_mfi].dataPtr(), d[r_mfi].dataPtr(),
                        sigma_node[mglev][r_mfi].dataPtr(), mask[mglev][r_mfi].dataPtr(), DIMLIST(freg),
                        &isRZ);
#else
	    FORT_HGRESUR(r[r_mfi].dataPtr(), DIMLIST(rbox),
			s[r_mfi].dataPtr(), d[r_mfi].dataPtr(),
                        sigma_node[mglev][r_mfi].dataPtr(), DIMLIST(freg),
                        &isRZ);
#endif
	}
    }
    HG_TEST_NORM(r, "level_residual: out");
}

void
holy_grail_amr_multigrid::relax (int  mglev,
                                 int  i1,
                                 bool is_zero)
{
    Box tdom = mg_domain[mglev];
    tdom.convert(IntVect::TheNodeVector());

    HG_TEST_NORM(corr[mglev], "relax corr a");
    HG_TEST_NORM(resid[mglev],  "relax resid  a");
    HG_TEST_NORM(cen[mglev],  "relax cen  a");
    HG_TEST_NORM(sigma[mglev],  "relax sigma  a");
    HG_DEBUG_OUT( "relax: i1 = " << i1 << "is_zero = " << is_zero << std::endl );
    for (int icount = 0; icount < i1; icount++)
    {
	HG_DEBUG_OUT( "icount = " << icount << std::endl );
	if (smoother_mode == 0 || smoother_mode == 1 || line_solve_dim == -1)
	{
	    if (is_zero == false)
		fill_borders(corr[mglev], lev_interface[mglev],
			     mg_boundary, -1, is_dense(m_stencil));
	    else
		is_zero = false;
	    HG_TEST_NORM(corr[mglev], "relax corr b");
	    for (MFIter r_mfi(resid[mglev]); r_mfi.isValid(); ++r_mfi)
	    {
		const Box& sbox = resid[mglev][r_mfi].box();
		const Box& freg = lev_interface[mglev].part_fine(r_mfi.index());
		if (line_solve_dim == -1)
		{
                    //
		    // Gauss-Seidel section:
                    //
		    if (m_stencil == terrain || m_stencil == full )
		    {
			const Box& fbox = corr[mglev][r_mfi].box();
			const Box& cenbox = cen[mglev][r_mfi].box();
			const Box& sigbox = sigma[mglev][r_mfi].box();
			if ( m_stencil == terrain )
			{
			    FORT_HGRLX_TERRAIN(
				corr[mglev][r_mfi].dataPtr(), DIMLIST(fbox),
				resid[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
				sigma[mglev][r_mfi].dataPtr(), DIMLIST(sigbox),
				cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
				DIMLIST(freg));
			}
			else
			{
			    FORT_HGRLX_FULL(
				corr[mglev][r_mfi].dataPtr(), DIMLIST(fbox),
				resid[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
				sigma[mglev][r_mfi].dataPtr(), DIMLIST(sigbox),
				cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
				DIMLIST(freg));
			}
		    }
		    else if (m_stencil == cross)
		    {
			const int isRZ = getCoordSys();
#if 1
			FORT_HGRLXU(corr[mglev][r_mfi].dataPtr(),
                                    resid[mglev][r_mfi].dataPtr(),
                                    sigma_node[mglev][r_mfi].dataPtr(),
                                    cen[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
                                    mask[mglev][r_mfi].dataPtr(),
                                    DIMLIST(freg),
				    &isRZ);
#else
			FORT_HGRLXUR(corr[mglev][r_mfi].dataPtr(),
                                    resid[mglev][r_mfi].dataPtr(),
                                    sigma_node[mglev][r_mfi].dataPtr(),
                                    cen[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
                                    DIMLIST(freg),
				    &isRZ);
#endif
		    }
		}
		else
		{
		    BoxLib::Abort( "holy_grail_amr_multigrid::relax():"
				   "line solves not implemented" );
                    //
		    // Grid-by-grid line solve section:
                    //
#if HG_CAN_DO_LINE_SOLVES
		    if(m_stencil == terrain || m_stencil == full )
		    {
			BoxLib::Abort( "holy_grail_amr_multigrid::relax():"
				       "Terrain line solves not implemented" );
		    }
		    const Box& fbox = corr[mglev][r_mfi].box();
		    const Box& cenbox = cen[mglev][r_mfi].box();
		    if (m_stencil == cross)
		    {
#if BL_SPACEDIM==2
			const Real hxyz[BL_SPACEDIM] = { D_DECL(h[mglev][0],
								h[mglev][1],
								h[mglev][2]) };
			const Box& sigbox = sigma_nd[0][mglev][r_mfi].box();
			const int isRZ = getCoordSys();
			const int imax = mg_domain[mglev].bigEnd(0) + 1;
			FORT_HGRLXL(
			    corr[mglev][r_mfi].dataPtr(), DIMLIST(fbox),
			    resid[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
			    D_DECL(sigma_nd[0][mglev][r_mfi].dataPtr(),
				   sigma_nd[1][mglev][r_mfi].dataPtr(),
				   sigma_nd[2][mglev][r_mfi].dataPtr()), DIMLIST(sigbox),
			    cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom),
			    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]),
			    &isRZ, &imax, &line_solve_dim);
#else
                        const Box& sigbox = sigma_node[mglev][r_mfi].box();
                        FORT_HGRLXL(
			    corr[mglev][r_mfi].dataPtr(), DIMLIST(fbox),
			    resid[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
			    sigma_node[mglev][r_mfi].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom), &line_solve_dim);
#endif
		    }
		    else
		    {
#if BL_SPACEDIM!=3
			const Real hxyz[BL_SPACEDIM] = { D_DECL(h[mglev][0],
								h[mglev][1],
								h[mglev][2]) };
			const Box& sigbox = sigma_nd[0][mglev][r_mfi].box();
			const int isRZ = getCoordSys();
			const int imax = mg_domain[mglev].bigEnd(0) + 1;
#if (BL_SPACEDIM==2)
			FORT_HGRLXL_FULL(
			    corr[mglev][r_mfi].dataPtr(), DIMLIST(fbox),
			    resid[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
			    D_DECL(sigma_nd[0][mglev][r_mfi].dataPtr(),
				   sigma_nd[1][mglev][r_mfi].dataPtr(),
				   sigma_nd[2][mglev][r_mfi].dataPtr()), DIMLIST(sigbox),
			    cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom),
			    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]),
			    &isRZ, &imax, &line_solve_dim);
#else
			FORT_HGRLXL(
			    corr[mglev][r_mfi].dataPtr(), DIMLIST(fbox),
			    resid[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
			    D_DECL(sigma_nd[0][mglev][r_mfi].dataPtr(),
				   sigma_nd[1][mglev][r_mfi].dataPtr(),
				   sigma_nd[2][mglev][r_mfi].dataPtr()), DIMLIST(sigbox),
			    cen[mglev][r_mfi].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg), DIMLIST(tdom),
			    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]));
#endif
#endif
		    }
#endif /* HG_CAN_DO_LINE_SOLVES */
		}
      }
      HG_TEST_NORM(corr[mglev], "relax corr b1");
      sync_borders(corr[mglev], lev_interface[mglev], mg_boundary);
      HG_TEST_NORM(corr[mglev], "relax corr b2");
    }
    else
    {
	BoxLib::Abort( "holy_grail_amr_multigrid::relax():"
		       "Line Solves aren't parallelized" );
#ifdef HG_CAN_DO_LINE_SOLVES
        //
	// Full-level line solve section:
        //
	if (line_order.length() == 0)
	{
	    build_line_order(line_solve_dim);
	}
	int lev = lev_min;
	while (ml_index[lev] < mglev)
	{
	    lev++;
	}

	for (int ipass = 0; ipass <= 1; ipass++)
	{
	    if (is_zero == false)
		fill_borders(corr[mglev], lev_interface[mglev],
			     mg_boundary, -1, is_dense(m_stencil));
	    else
		is_zero = false;
	    //
	    // Forward solve:
            //
	    task_list tl;
	    for (int i = 0; i < mg_mesh[mglev].length(); i++)
	    {
		//
		// Do grids in order along line_solve_dim:
                //
		int igrid = line_order[lev][i];
		const Box& fbox = corr[mglev][igrid].box();
		const Box& sbox = resid[mglev][igrid].box();
		const Box& wbox = work[mglev][igrid].box();
		const Box& cenbox = cen[mglev][igrid].box();
		const Box& freg = corr[mglev].box(igrid);
		if (m_stencil == terrain || m_stencil == full )
		{
		    const Box& sigbox = sigma[mglev][igrid].box();
		    if ( m_stencil == terrain )
		    {
			FORT_HGRLNF_TERRAIN(
			    corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			    resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			    work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			    sigma[mglev][igrid].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg),
			    DIMLIST(tdom), &line_solve_dim, &ipass);
		    }
		    else if ( m_stencil == full )
		    {
			FORT_HGRLNF_FULL(
			    corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			    resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			    work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			    sigma[mglev][igrid].dataPtr(), DIMLIST(sigbox),
			    cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			    DIMLIST(freg),
			    DIMLIST(tdom), &line_solve_dim, &ipass);
		    }
		}
		else if (m_stencil == cross)
		{
#if BL_SPACEDIM == 2
		    const Real hxyz[BL_SPACEDIM] = { D_DECL(h[mglev][0],
							    h[mglev][1],
							    h[mglev][2]) };
		    const int isRZ = getCoordSys();
		    const int imax = mg_domain[mglev].bigEnd(0) + 1;
		    const Box& sigbox = sigma_nd[0][mglev][igrid].box();
		    FORT_HGRLNF(
			corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			D_DECL(
			    sigma_nd[0][mglev][igrid].dataPtr(),
			    sigma_nd[1][mglev][igrid].dataPtr(),
			    sigma_nd[2][mglev][igrid].dataPtr()), DIMLIST(sigbox),
			cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			DIMLIST(freg), DIMLIST(tdom),
			D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]),
			&isRZ, &imax, &line_solve_dim, &ipass);
#else
		    const Box& sigbox = sigma_node[mglev][igrid].box();
		    FORT_HGRLNF_TERRAIN(
			corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
			resid[mglev][igrid].dataPtr(), DIMLIST(sbox),
			work[mglev][igrid].dataPtr(), DIMLIST(wbox),
			sigma_node[mglev][igrid].dataPtr(), DIMLIST(sigbox),
			cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
			DIMLIST(freg), DIMLIST(tdom), &line_solve_dim, &ipass);
#endif
		}
                //
		// Copy work arrays to following grids:
                //
		for (ListIterator<int> j(line_after[lev][igrid]); j; j++)
		{
		    const int jgrid = j();
		    const Box b = (freg & corr[mglev].box(jgrid));
		    tl.add_task(
			new task_copy_link(corr[mglev], jgrid, igrid, b, 0));
		    tl.add_task(
			new task_copy_link(work[mglev], jgrid, igrid, b, 0));
		}
	    }
	    tl.execute("holy_grail_amr_multigrid::relax(1)");
            //
	    // Back substitution:
            //
	    for (int i = mg_mesh[mglev].length() - 1; i >= 0; i--)
	    {
		//
		// Do grids in reverse order along line_solve_dim:
                //
		const int igrid = line_order[lev][i];
		const Box& freg = corr[mglev].box(igrid);
		//
		// Copy solution array from following grids:
                //
		for (ListIterator<int> j(line_after[lev][igrid]); j; j++)
		{
		    const int jgrid = j();
		    const Box b = (freg & corr[mglev].box(jgrid));
		    tl.add_task(
			new task_copy_link(corr[mglev], igrid, jgrid, b, 0));
		}

		const Box& fbox = corr[mglev][igrid].box();
		const Box& wbox = work[mglev][igrid].box();
		FORT_HGRLNB(
		    corr[mglev][igrid].dataPtr(), DIMLIST(fbox),
		    work[mglev][igrid].dataPtr(), DIMLIST(wbox),
		    DIMLIST(freg), &line_solve_dim, &ipass);
	    }
	    tl.execute("holy_grail_amr_multigrid::relax(2)");
	}
#endif /* HG_CAN_DO_LINE_SOLVES */
   }
  }
    HG_TEST_NORM(corr[mglev], "relax a1");
}

void
holy_grail_amr_multigrid::build_line_order (int lsd)
{
    line_order.resize(lev_max + 1);
    line_after.resize(lev_max + 1);

    for (int lev = lev_min; lev <= lev_max; lev++)
    {
	int mglev = ml_index[lev], ngrids = mg_mesh[mglev].size();

	line_order[lev].resize(ngrids);
	line_after[lev].resize(ngrids);

	for (int igrid = 0; igrid < ngrids; igrid++)
	{
	    line_order[lev].set(igrid, igrid);
	    //
	    // bubble sort, replace with something faster if necessary:
            //
	    for (int i = igrid; i > 0; i--)
	    {
		if (ml_mesh[lev][line_order[lev][i]].smallEnd(lsd)
		    < ml_mesh[lev][line_order[lev][i-1]].smallEnd(lsd))
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
		if (BoxLib::bdryLo(ml_mesh[lev][i], lsd).intersects(BoxLib::bdryHi(ml_mesh[lev][igrid], lsd)))
		{
		    line_after[lev][igrid].push_back(i);
		}
	    }
	}
    }
}

void
holy_grail_amr_multigrid::cgsolve (int mglev)
{
    BL_ASSERT(mglev == 0);

    MultiFab& r = cgwork[0];
    MultiFab& p = cgwork[1];
    MultiFab& z = cgwork[2];
    MultiFab& x = cgwork[3];
    MultiFab& w = cgwork[4];
    MultiFab& c = cgwork[5];
    MultiFab& zero_array = cgwork[6];
    MultiFab& ipmask = cgwork[7];
    //
    // x (corr[0]) should be all 0.0 at this point
    //
    for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
    {
	r[r_mfi].copy(resid[mglev][r_mfi]);
	r[r_mfi].negate();
    }

    if (singular)
    {
        //
	// Singular systems are very sensitive to solvability
        //
	w.setVal(1.0);
	Real aa = inner_product(r, w) / mg_domain[mglev].volume();
	r.plus(-aa, 0);
    }

    Real rho = 0.0;
    for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
    {
	z[r_mfi].copy(r[r_mfi]);
	z[r_mfi].mult(c[r_mfi]);
	const Box& reg = p[r_mfi].box();
	FORT_HGIP(z[r_mfi].dataPtr(), r[r_mfi].dataPtr(),
                  ipmask[r_mfi].dataPtr(), DIMLIST(reg), rho);
	p[r_mfi].copy(z[r_mfi]);
    }
    ParallelDescriptor::ReduceRealSum(rho);
    if ( pcode >= 3 && ParallelDescriptor::IOProcessor() )
    {
	std::cout << "HG: cgsolve rho = " << rho << std::endl;
    }

    const Real tol = HG::cgsolve_tolfact * rho;

    int i = 0;
    while (tol > 0.0)
    {
	if ( ++i > HG::cgsolve_maxiter && ParallelDescriptor::IOProcessor() )
	{
	    BoxLib::Error( "cgsolve:"
			   "Conjugate-gradient iteration failed to converge" );
	}
	Real rho_old = rho;
        //
	// safe to set the clear flag to 0 here---bogus values make it
	// into r but are cleared from z by the mask in c
        //
	level_residual(w, zero_array, p, 0, false, 0);
	Real alpha = 0.0;
	for (MFIter p_mfi(p); p_mfi.isValid(); ++p_mfi)
	{
	    const Box& reg = p[p_mfi].box();
	    FORT_HGIP(p[p_mfi].dataPtr(), w[p_mfi].dataPtr(),
                      ipmask[p_mfi].dataPtr(), DIMLIST(reg), alpha);
	}
	ParallelDescriptor::ReduceRealSum(alpha);
	alpha = rho / alpha;
	rho = 0.0;
	for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    const Box& reg = p[r_mfi].box();
	    FORT_HGCG1(r[r_mfi].dataPtr(),
                       p[r_mfi].dataPtr(),
                       z[r_mfi].dataPtr(),
		       x[r_mfi].dataPtr(),
                       w[r_mfi].dataPtr(),
                       c[r_mfi].dataPtr(),
		       ipmask[r_mfi].dataPtr(),
		       DIMLIST(reg), alpha, rho);
	}
	ParallelDescriptor::ReduceRealSum(rho);
	if (pcode >= 3  && ParallelDescriptor::IOProcessor())
	{
	    std::cout << "HG: cgsolve iter(" << i << ") rho=" << rho << std::endl;
	}
	if (rho <= tol)
	    break;
	alpha = rho / rho_old;
	for (MFIter p_mfi(p); p_mfi.isValid(); ++p_mfi)
	{
	    const Box& reg = p[p_mfi].box();
	    FORT_HGCG2(p[p_mfi].dataPtr(),z[p_mfi].dataPtr(),DIMLIST(reg),alpha);
	}
    }

    if (pcode >= 2  && ParallelDescriptor::IOProcessor())
    {
	std::cout << "HG: "
                  << i << " iterations required for conjugate-gradient" << std::endl;
    }
}
