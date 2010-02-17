#include <iostream>

#include "hg_multi.H"

#include <Profiler.H>

#if defined( BL_FORT_USE_UNDERSCORE )
#define   FORT_HGRES_TERRAIN	hgres_terrain_
#define   FORT_HGRES_FULL	hgres_full_
#define   FORT_HGRES_CROSS	hgres_cross_
#define   FORT_HGRLX_TERRAIN	hgrlx_terrain_
#define   FORT_HGRLX_FULL	hgrlx_full_
#define   FORT_HGRLX		hgrlx_
#define   FORT_HGRESU		hgresu_
#define   FORT_HGRLXU		hgrlxu_
#elif defined( BL_FORT_USE_UPPERCASE )
#define   FORT_HGRES_TERRAIN    HGRES_TERRAIN
#define   FORT_HGRES_FULL	HGRES_FULL
#define   FORT_HGRES_CROSS	HGRES_CROSS
#define   FORT_HGRLX		HGRLX
#define   FORT_HGRLX_TERRAIN    HGRLX_TERRAIN
#define   FORT_HGRLX_FULL	HGRLX_FULL
#define   FORT_HGRLXL		HGRLXL
#define   FORT_HGRLXL_FULL	HGRLXL_FULL
#define   FORT_HGRESU		HGRESU
#define   FORT_HGRLXU		HGRLXU
#elif defined( BL_FORT_USE_LOWERCASE )
#define   FORT_HGRES_TERRAIN    hgres_terrain
#define   FORT_HGRES_FULL	hgres_full
#define   FORT_HGRES_CROSS	hgres_cross
#define   FORT_HGRLX		hgrlx
#define   FORT_HGRLX_TERRAIN    hgrlx_terrain
#define   FORT_HGRLX_FULL	hgrlx_full
#define   FORT_HGRESU		hgresu
#define   FORT_HGRLXU		hgrlxu
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_HGRES_TERRAIN (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRLX_TERRAIN (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     const Real*, intS, intS);
#if (BL_SPACEDIM==3)
    void FORT_HGRESU        (Real*, intS, const Real*, const Real*,
			     const Real*, const Real*, intS, const int*);

    void FORT_HGRLXU        (Real*, const Real*, const Real*,
			     const Real*, intS, const Real*, intS, const int*);

    void FORT_HGRES_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRES_CROSS   (Real*, intS, const Real*, const Real*,
			     const Real*, intS, const int*);
    void FORT_HGRLX_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS, intS);
    void FORT_HGRLX          (Real*, const Real*, const Real*,
			     const Real*, intS, intS, const int*);
    void FORT_HGRLXL        (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS, intS,
			     intS, const int*);
#else
    void FORT_HGRES_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS, const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRES_CROSS    (Real*, intS, const Real*, const Real*,
			      const Real*, intS, const int*);
    void FORT_HGRLX_FULL    (Real*, intS, const Real*, intS,
			     const Real*, intS,
			     const Real*, intS, intS);
    void FORT_HGRLX          (Real*, const Real*, const Real*,
			     const Real*, intS, intS, const int*);
    void FORT_HGRLXL        (Real*, intS, const Real*, intS,
			     CRealPS, intS, Real*, intS, intS, intS,
			     CRealPS, const int*, const int*, const int*);
    void FORT_HGRLXL_FULL   (Real*, intS, const Real*, intS,
			     CRealPS, intS, Real*, intS, intS, intS,
			     CRealPS, const int*, const int*, const int*);
#endif
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

#if 0
	    FORT_HGRES_CROSS(r[r_mfi].dataPtr(), DIMLIST(rbox),
			     s[r_mfi].dataPtr(), d[r_mfi].dataPtr(),
                             sigma_node[mglev][r_mfi].dataPtr(), DIMLIST(freg),
                             &isRZ);
#else
	    FORT_HGRESU(r[r_mfi].dataPtr(), DIMLIST(rbox),
			s[r_mfi].dataPtr(), d[r_mfi].dataPtr(),
                        sigma_node[mglev][r_mfi].dataPtr(), mask[mglev][r_mfi].dataPtr(), DIMLIST(freg),
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::relax()");

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

#if 0
			FORT_HGRLX(corr[mglev][r_mfi].dataPtr(),
                                  resid[mglev][r_mfi].dataPtr(),
                                  sigma_node[mglev][r_mfi].dataPtr(),
                                  cen[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
                                  DIMLIST(freg),&isRZ);
#else
			FORT_HGRLXU(corr[mglev][r_mfi].dataPtr(),
                                    resid[mglev][r_mfi].dataPtr(),
                                    sigma_node[mglev][r_mfi].dataPtr(),
                                    cen[mglev][r_mfi].dataPtr(), DIMLIST(sbox),
                                    mask[mglev][r_mfi].dataPtr(),
                                    DIMLIST(freg),
				    &isRZ);
#endif
		    }
		}
		else
		{
		    BoxLib::Abort( "holy_grail_amr_multigrid::relax():"
				   "line solves not implemented" );
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

