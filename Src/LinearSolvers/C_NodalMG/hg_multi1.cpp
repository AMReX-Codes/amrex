
#include "hg_multi.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define   FORT_HGSRST		hgsrst_
#define   FORT_HGSCON		hgscon_
#define   FORT_HGCEN		hgcen_
#define   FORT_HGCEN_TERRAIN	hgcen_terrain_
#define   FORT_HGCEN_FULL	hgcen_full_
#define   FORT_HGINTS		hgints_
#define   FORT_HGINTS_DENSE	hgints_dense_
#define   FORT_FACRST1		acrst1_
#define   FORT_FANRST2		anrst2_
#elif defined( BL_FORT_USE_UPPERCASE )
#define   FORT_HGSRST		HGSRST
#define   FORT_HGSCON		HGSCON
#define   FORT_HGCEN		HGCEN
#define   FORT_HGCEN_TERRAIN	HGCEN_TERRAIN
#define   FORT_HGCEN_FULL	HGCEN_FULL
#define   FORT_HGINTS		HGINTS
#define   FORT_HGINTS_DENSE	HGINTS_DENSE
#define   FORT_FACRST1		ACRST1
#define   FORT_FANRST2		ANRST2
#elif defined( BL_FORT_USE_LOWERCASE )
#define   FORT_HGSRST		hgsrst
#define   FORT_HGSCON		hgscon
#define   FORT_HGCEN		hgcen
#define   FORT_HGCEN_TERRAIN	hgcen_terrain
#define   FORT_HGCEN_FULL	hgcen_full
#define   FORT_HGINTS		hgints
#define   FORT_HGINTS_DENSE	hgints_dense
#define   FORT_FACRST1		acrst1
#define   FORT_FANRST2		anrst2
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_FACRST1       (Real*, intS, intS, const Real*, intS, intRS,
			     const int&, const int*, const int*, const int*);
    void FORT_FANRST2       (Real*, intS, intS, const Real*, intS, intRS,
			     const int&, const int*, const int*, const int*);

    void FORT_HGSRST        (RealPS, intS, intS, CRealPS, intS, intRS);
    void FORT_HGINTS        (Real*, intS, intS, Real*, intS,
			     const Real*, intS, intS, intRS);
    void FORT_HGINTS_DENSE  (Real*, intS, intS, CRealPS, intS,
			     const Real*, intS, intS, intRS);
    void FORT_HGCEN         (Real*, intS, Real*, intS, intS, const int*);
    void FORT_HGCEN_FULL    (Real*, intS, Real*, intS, intS);
    void FORT_HGCEN_TERRAIN (Real*, intS, Real*, intS, intS);
    void FORT_HGSCON        (Real*, intS, RealPS, intS, intS, CRealPS);
}

class task_interpolate_patch : public task
{
public:
    task_interpolate_patch (task_list&              tl_,
                            MultiFab&               dmf_,
                            int                     dgrid_,
                            const Box&              dbx_,
                            const MultiFab&         smf_,
                            const IntVect&          rat_,
                            const amr_interpolator* interp_,
                            const level_interface&  lev_interface_);
    virtual bool ready ();
    virtual ~task_interpolate_patch ();

private:

    task::task_proxy        tf;
    MultiFab&               dmf;
    const int               dgrid;
    const Box               dbx;
    const MultiFab&         smf;
    const IntVect           rat;
    const amr_interpolator* interp;
    const level_interface&  lev_interface;
};

task_interpolate_patch::task_interpolate_patch (task_list&              tl_,
                                                MultiFab&               dmf_,
                                                int                     dgrid_,
                                                const Box&              dbx_,
                                                const MultiFab&         smf_,
                                                const IntVect&          rat_,
                                                const amr_interpolator* interp_,
                                                const level_interface&  lev_interface_)
    :
    task(tl_),
    dmf(dmf_),
    dgrid(dgrid_),
    dbx(dbx_),
    smf(smf_),
    rat(rat_),
    interp(interp_),
    lev_interface(lev_interface_),
    tf(0)
{
    BL_ASSERT(dbx.sameType(dmf.box(dgrid)));

    tf = m_task_list.add_task(new task_fill_patch(m_task_list,
                                                  dmf,
                                                  dgrid,
                                                  interp->box(dbx, rat),
                                                  smf,
                                                  lev_interface,
                                                  0,
                                                  -1,
                                                  -1));
    depend_on(tf);
}

bool
task_interpolate_patch::ready ()
{
    if (is_local(dmf, dgrid))
    {
        BL_ASSERT(is_started());
        BL_ASSERT(!tf.null());
	task_fab* tff = dynamic_cast<task_fab*>(tf.get());
        BL_ASSERT(tff != 0);
	interp->fill(dmf[dgrid], dbx, tff->fab(), tff->fab().box(), rat);
    }
    return true;
}

task_interpolate_patch::~task_interpolate_patch ()
{
    delete interp;
}

class holy_grail_interpolator_dense
    :
    public bilinear_interpolator
{
public:

    holy_grail_interpolator_dense(Real* Sigptr[BL_SPACEDIM],
				  const Box& Sigbox)
	: sigbox(Sigbox)
    {
        D_TERM(sigptr[0] = Sigptr[0];,
               sigptr[1] = Sigptr[1];,
               sigptr[2] = Sigptr[2];);
    }
    virtual void fill (FArrayBox& patch,
                       const Box& region,
                       const FArrayBox& cgr,
                       const Box& cb,
                       const IntVect& rat) const;
protected:

    Real*     sigptr[BL_SPACEDIM];
    const Box sigbox;
};

class holy_grail_interpolator
    :
    public bilinear_interpolator
{
public:
    holy_grail_interpolator(Real *Sigptr,
			    const Box& Sigbox)
	: sigptr(Sigptr),
	  sigbox(Sigbox)
    {}
    virtual void fill (FArrayBox&       patch,
                       const Box&       region,
                       const FArrayBox& cgr,
                       const Box&       cb,
                       const IntVect&   rat) const;
protected:

    Real*     sigptr;
    const Box sigbox;
};

class holy_grail_sigma_restrictor
    :
    public cell_average_restrictor
{
public:
    explicit holy_grail_sigma_restrictor(holy_grail_amr_multigrid::stencil m_hg_stencil_)
	:
        cell_average_restrictor(0),
        m_hg_stencil(m_hg_stencil_)
    {}
    virtual void fill (FArrayBox&       patch,
                       const Box&       region,
                       const FArrayBox& fgr,
                       const IntVect&   rat) const;
private:

    holy_grail_amr_multigrid::stencil m_hg_stencil;
};

void
holy_grail_amr_multigrid::set_line_solve_dimension (int dim)
{
  if (dim != -1)
  {
      BoxLib::Abort(
          "holy_grail_amr_multigrid::holy_grail_amr_multigrid():"
          "LineSolves not supported in parallel" );
  }
  line_solve_dim = dim;
}

void
holy_grail_amr_multigrid::set_smoother_mode (int mode)
{
    smoother_mode = mode;
}

bool
holy_grail_amr_multigrid::is_dense (stencil sval)
{
    return sval == terrain || sval == full;
}

holy_grail_amr_multigrid::holy_grail_amr_multigrid(const Array<BoxArray>& Mesh,
						   const Array<IntVect>&  Gen_ratio,
						   const Box&             fdomain,
						   int                    Lev_min_min,
						   int                    Lev_min_max,
						   int                    Lev_max_max,
						   const amr_fluid_boundary& Boundary,
						   stencil                stencil_,
						   int                    Pcode)
    : amr_multigrid(Mesh,
                    Gen_ratio,
                    Lev_min_min,
                    Lev_min_max,
                    Lev_max_max,
                    Boundary.pressure(),
                    Pcode),
      boundary(Boundary),
      smoother_mode(2),
      line_solve_dim(-1),
      m_stencil(stencil_)
{
    build_mesh(fdomain);
}

void
holy_grail_amr_multigrid::alloc_hg_multi (PArray<MultiFab>& Dest,
                                          PArray<MultiFab>& Source,
                                          PArray<MultiFab>& Coarse_source,
                                          PArray<MultiFab>& Sigma,
                                          Real              H[],
                                          int               Lev_min,
                                          int               Lev_max,
                                          int               for_fill_sync_reg)
{
    BL_ASSERT(Dest.size() > Lev_max);
    BL_ASSERT(Dest[Lev_min].nGrow() == 1);

    if (Source.size())
    {
	source_owned = false;
	amr_multigrid::alloc_amr_multi(Dest, Source, Coarse_source, Lev_min, Lev_max);
    }
    else
    {
	source_owned = true;
	PArray<MultiFab> Src;
	Src.resize(Lev_max + 1);
	for (int lev = Lev_min; lev <= Lev_max; lev++)
	{
	    const BoxArray& mesh = Dest[lev].boxArray();
	    Src.set(lev, new MultiFab(mesh, 1, Dest[Lev_min].nGrow()));
	    Src[lev].setVal(0.0);
	}
	amr_multigrid::alloc_amr_multi(Dest, Src,    Coarse_source, Lev_min, Lev_max);
    }

    h = new Real[mglev_max + 1][BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	h[mglev_max][i] = H[i];
	for (int mglev = mglev_max - 1; mglev >= 0; mglev--)
	{
	    int rat =
		mg_domain[mglev+1].length(i) / mg_domain[mglev].length(i);
	    h[mglev][i] = rat * h[mglev+1][i];
	}
    }

    build_sigma(Sigma, for_fill_sync_reg);

    if (for_fill_sync_reg > 0) return;

    alloc_sync_caches();

    int ib = dest[lev_min].nGrow();
    const BoxArray& mesh0 = corr[0].boxArray();

    cgwork.resize(8);
    cgwork.set(0, new MultiFab(mesh0, 1, ib));
    cgwork.set(1, new MultiFab(mesh0, 1, ib));
    cgwork.set(2, new MultiFab(mesh0, 1, ib));
    cgwork.set(3, &corr[0]);
    cgwork.set(4, &work[0]);
    cgwork.set(5, &cen[0]);
    cgwork.set(6, new MultiFab(mesh0, 1, ib));
    cgwork[6].setVal(0.0);
    cgwork.set(7, new MultiFab(mesh0, 1, ib));

    BL_ASSERT(cgwork[3].nGrow() == ib &&
	cgwork[4].nGrow() == ib &&
	cgwork[5].nGrow() == ib);

    for (MFIter g_mfi(cgwork[7]); g_mfi.isValid(); ++g_mfi)
    {
	FArrayBox& gtmp = cgwork[7][g_mfi];
	const Box& valid = g_mfi.validbox();
	gtmp.setVal(0.0);
	gtmp.setVal(1.0, valid, 0);
	Box b = BoxLib::bdryLo(valid, 1);
	gtmp.setVal(0.5, b, 0);
	b = BoxLib::bdryHi(valid, 1);
	gtmp.setVal(0.5, b, 0);
	b = BoxLib::bdryLo(valid, 0);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, BoxLib::bdryLo(b, 1), 0);
	gtmp.setVal(0.25, BoxLib::bdryHi(b, 1), 0);
	b = BoxLib::bdryHi(valid, 0);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, BoxLib::bdryLo(b, 1), 0);
	gtmp.setVal(0.25, BoxLib::bdryHi(b, 1), 0);
#if (BL_SPACEDIM == 3)
	b = BoxLib::bdryLo(valid, 2);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, BoxLib::bdryLo(b, 0), 0);
	gtmp.setVal(0.25, BoxLib::bdryHi(b, 0), 0);
	Box bb = BoxLib::bdryLo(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, BoxLib::bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, BoxLib::bdryHi(bb, 0), 0);
	bb = BoxLib::bdryHi(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, BoxLib::bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, BoxLib::bdryHi(bb, 0), 0);
	b = BoxLib::bdryHi(valid, 2);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, BoxLib::bdryLo(b, 0), 0);
	gtmp.setVal(0.25, BoxLib::bdryHi(b, 0), 0);
	bb = BoxLib::bdryLo(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, BoxLib::bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, BoxLib::bdryHi(bb, 0), 0);
	bb = BoxLib::bdryHi(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, BoxLib::bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, BoxLib::bdryHi(bb, 0), 0);
#endif
    }

    singular = false;
    if (mg_boundary->singular())
    {
	long sng = 0;
	for (int i = 0; i < mg_mesh[0].size(); i++)
	{
	    sng += mg_mesh[0][i].numPts();
	}
	singular = (sng == mg_domain[0].numPts());
    }

    // if (m_stencil == terrain)	// FIXME
    if ( is_dense(m_stencil) )
	integrate = 1;
}


void
holy_grail_sigma_restrictor::fill (FArrayBox&       patch,
				   const Box&       region,
				   const FArrayBox& fgr,
				   const IntVect&   rat) const
{
    BL_ASSERT(patch.box().cellCentered());
    BL_ASSERT(rat[0] == 2 && rat[1] == 2 ||
	rat[0] == 2 && rat[1] == 1 ||
	rat[0] == 1 && rat[1] == 2);

    if (m_hg_stencil == holy_grail_amr_multigrid::terrain)
    {
	FORT_HGSRST(
	    D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
	    DIMLIST(patch.box()),
	    DIMLIST(region),
	    D_DECL( fgr.dataPtr(0), fgr.dataPtr(1), fgr.dataPtr(2)),
	    DIMLIST(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]));
	const int integ = 0;
	FORT_FACRST1(patch.dataPtr(BL_SPACEDIM),
                     DIMLIST(patch.box()),
                     DIMLIST(region),
                     fgr.dataPtr(BL_SPACEDIM),
                     DIMLIST(fgr.box()),
                     D_DECL(rat[0], rat[1], rat[2]), 1, &integ, 0, 0);

#if (BL_SPACEDIM == 2)
	patch.mult(Real(rat[1]) / rat[0], region, 0, 1);
	patch.mult(Real(rat[0]) / rat[1], region, 1, 1);
	// component 2 remains unchanged
#else
	FORT_FACRST1(patch.dataPtr(BL_SPACEDIM+1),
                     DIMLIST(patch.box()),
                     DIMLIST(region),
                     fgr.dataPtr(BL_SPACEDIM+1),
                     DIMLIST(fgr.box()),
                     D_DECL(rat[0], rat[1], rat[2]), 1, &integ, 0, 0);
	patch.mult(Real(rat[1] * rat[2]) / rat[0], region, 0, 1);
	patch.mult(Real(rat[0] * rat[2]) / rat[1], region, 1, 1);
	patch.mult(Real(rat[0] * rat[1]) / rat[2], region, 2, 1);
	patch.mult(Real(rat[1]),                   region, 3, 1);
	patch.mult(Real(rat[0]),                   region, 4, 1);
#endif
    }
    else if (m_hg_stencil ==  holy_grail_amr_multigrid::full)
    {
	FORT_HGSRST(
	    D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
	    DIMLIST(patch.box()),
	    DIMLIST(region),
	    D_DECL( fgr.dataPtr(0), fgr.dataPtr(1), fgr.dataPtr(2)),
	    DIMLIST(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]));
#if (BL_SPACEDIM == 2)
	patch.mult(Real(rat[1]) / rat[0], region, 0, 1);
	patch.mult(Real(rat[0]) / rat[1], region, 1, 1);
#else
	patch.mult(Real(rat[1] * rat[2]) / rat[0], region, 0, 1);
	patch.mult(Real(rat[0] * rat[2]) / rat[1], region, 1, 1);
	patch.mult(Real(rat[0] * rat[1]) / rat[2], region, 2, 1);
#endif
    }
    else
    {
	if (fgr.nComp() == 1)
	{
	    FORT_HGSRST(
		D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
		DIMLIST(patch.box()),
		DIMLIST(region),
		D_DECL(fgr.dataPtr(), fgr.dataPtr(), fgr.dataPtr()),
		DIMLIST(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
	}
	else
	{
	    FORT_HGSRST(
		D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
		DIMLIST(patch.box()),
		DIMLIST(region),
		D_DECL(fgr.dataPtr(0), fgr.dataPtr(1), fgr.dataPtr(2)),
		DIMLIST(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
	}
    }
}

void
holy_grail_amr_multigrid::build_sigma (PArray<MultiFab>& Sigma,
                                       int for_fill_sync_reg)
{
    if (m_stencil == terrain || m_stencil == full )
    {
	//
	// For terrain stencils we have as many sigma arrays passed as
	// arguments and used at the lev_interface as we build for internal
	// multigrid purposes.  This simplifies handling as we do not
	// need to maintain separate arrays for different purposes.
	//
	const int ncomp = m_stencil == terrain ?
	    2 * BL_SPACEDIM - 1 :
	    BL_SPACEDIM;

	sigma.resize(mglev_max+1);

	for (int mglev = 0; mglev <= mglev_max; mglev++)
	{
            if (for_fill_sync_reg == 0 || mglev == mglev_max)
            {
		sigma.set(mglev, new MultiFab(mg_mesh[mglev], ncomp, 1));
		sigma[mglev].setVal(1.0e50);
		int lev;
		if ((lev = get_amr_level(mglev)) >= 0)
		{
		    if ( Sigma[lev].nComp() != 1 && 
			 Sigma[lev].nComp() != ncomp)
		    {
		      BoxLib::Error("Sigma has wrong number of components");
		    }
		    for (MFIter s_mfi(Sigma[lev]); s_mfi.isValid(); ++s_mfi)
		    {
			Box copyon(s_mfi.validbox());
                        if (for_fill_sync_reg != 0) copyon.grow(1);
		        if (Sigma[lev].nComp() == 1)
		        {
			    for (int i = 0; i < ncomp; ++i)
		 	    {
				sigma[mglev][s_mfi].copy(Sigma[lev][s_mfi], copyon, 0, copyon, i, 1);
			    }
			}
		        else if (Sigma[lev].nComp() == ncomp)
			    sigma[mglev][s_mfi].copy(Sigma[lev][s_mfi], copyon, 0, copyon, 0, ncomp);
		    }
		}
	    }
	}

        if (for_fill_sync_reg == 0)
	{
	    for (int mglev = mglev_max; mglev > 0; mglev--)
	    {
		IntVect rat = mg_domain[mglev].length()
		    / mg_domain[mglev-1].length();
		restrict_level(
		    sigma[mglev-1], sigma[mglev], rat,
		    holy_grail_sigma_restrictor(m_stencil),
		    default_level_interface, 0);
	    }
	    for (int mglev = 0; mglev <= mglev_max; mglev++)
	    {
		// FIXME, not terrain sigma?
		fill_borders(sigma[mglev], lev_interface[mglev],
			     boundary.terrain_sigma(),
			     -1, is_dense(m_stencil));
		HG_TEST_NORM( sigma[mglev], "build_sigma bb");
	    }
        }
        else if (for_fill_sync_reg == 1)
        {
	    // FIXME, not terrain sigma?
            boundary.terrain_sigma()->
  	      fill_sync_reg_borders(sigma[mglev_max], 
				    lev_interface[mglev_max],-1);
        }
    }
    else
    {
	//
	// Intended functionality:  sigma_split exists only at coarser levels,
	// since only after coarsening sigma is different in different directions.
	// sigma exists at all levels, and is intended for use on fine grids
	// and at lev_interface points, where all components are the same.  To
	// save storage it is aliased to the first component of sigma_split
	// on all but the finest level.
	//
	// sigma_split replaced by sigma_nd in more recent version, used
	// only as a local variable here
	//
	PArray<MultiFab> sigma_split;
        if (for_fill_sync_reg == 0)
	{
	    sigma_split.resize(mglev_max);
	    for (int mglev = 0; mglev < mglev_max; mglev++)
	    {
		sigma_split.set(mglev,
				new MultiFab(mg_mesh[mglev], BL_SPACEDIM, 1));
	    }
	}
	sigma.resize(mglev_max + 1);
	sigma.set(mglev_max, new MultiFab(mg_mesh[mglev_max], 1, 1));
	//
	// Level project:
	// Any value can fill values in the border cells that fill_borders
	// will not touch---those touching coarser grids.  The values in these
	// cells will only be seen by the interpolation, and the quantity
	// being interpolated will always be zero at these edges, but we
	// should insure that no NaN's or other garbage is there that could
	// cause a floating point fault.
	//
	// Sync project:
	// Ghost values will be seen by multilevel interpolation, so put
	// a huge value in ghost cells so that coarse-fine lev_interface
	// interpolation will linear, as finite element derivation requires.
	//
        if (for_fill_sync_reg == 0)
	{
	    for (int mglev = 0; mglev < mglev_max; mglev++)
	    {
		MultiFab& target = sigma_split[mglev];
		target.setVal(1.0e50);
		int lev;
		if ((lev = get_amr_level(mglev)) >= 0)
		{
		    MultiFab& s_comp = Sigma[lev];
		    for (int i = 0; i < BL_SPACEDIM; i++)
		    {
			for (MFIter s_mfi(s_comp); s_mfi.isValid(); ++s_mfi)
			{
			    target[s_mfi].copy(s_comp[s_mfi], s_mfi.validbox(),
                                               0, target.boxArray()[s_mfi.index()], i, 1);
			}
		    }
		}
	    }
	}
	sigma[mglev_max].setVal(1.0e50);
	HG_TEST_NORM(sigma[mglev_max], "build_sigma aa1");
	HG_TEST_NORM(Sigma[  lev_max], "build_sigma aa10");

	for (MFIter S_mfi(Sigma[lev_max]); S_mfi.isValid(); ++S_mfi)
	{
   	     if (for_fill_sync_reg == 0)
             {
		sigma[mglev_max][S_mfi].copy(
		    Sigma[lev_max][S_mfi], mg_mesh[mglev_max][S_mfi.index()], 0,
		    mg_mesh[mglev_max][S_mfi.index()], 0, 1);
  	     }
             else
  	     {
		sigma[mglev_max][S_mfi].copy(
		    Sigma[lev_max][S_mfi], BoxLib::grow(mg_mesh[mglev_max][S_mfi.index()], 1), 0,
		    BoxLib::grow(mg_mesh[mglev_max][S_mfi.index()], 1), 0, 1);
	     }
	}

	HG_TEST_NORM(sigma[mglev_max], "build_sigma aa2");
        if (for_fill_sync_reg == 0)
	{
	    if (mglev_max > 0)
	    {
		IntVect rat = mg_domain[mglev_max].length()
		    / mg_domain[mglev_max-1].length();
		restrict_level(
		    sigma_split[mglev_max-1], sigma[mglev_max], rat,
		    holy_grail_sigma_restrictor(m_stencil),
		    default_level_interface, 0);
	    }

	    fill_borders(sigma[mglev_max], lev_interface[mglev_max],
			 boundary.scalar(),
			 -1, is_dense(m_stencil));

	    for (int mglev = mglev_max - 1; mglev > 0; mglev--)
	    {
		IntVect rat = mg_domain[mglev].length()
		    / mg_domain[mglev-1].length();
		restrict_level(
		    sigma_split[mglev-1], sigma_split[mglev], rat,
		    holy_grail_sigma_restrictor(m_stencil),
		    default_level_interface, 0);
	    }

	    for (int mglev = 0; mglev < mglev_max; mglev++)
	    {
		HG_TEST_NORM(sigma_split[mglev], "build_sigma 0");
		fill_borders(sigma_split[mglev], lev_interface[mglev],
			     boundary.scalar(),
			     -1, is_dense(m_stencil));
		HG_TEST_NORM(sigma_split[mglev], "build_sigma");
	    }
        }
        else if (for_fill_sync_reg == 1)
        {
  	    boundary.scalar()->fill_sync_reg_borders(sigma[mglev_max],
  					             lev_interface[mglev_max],
                                                     -1);
        }

	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    sigma_nd[i].resize(mglev_max + 1);
	}

        if (for_fill_sync_reg == 0)
	{
	    for (int mglev = 0; mglev < mglev_max; mglev++)
	    {
		MultiFab& s = sigma_split[mglev];
		for (int i = 0; i < BL_SPACEDIM; i++)
		{
		    sigma_nd[i].set(mglev, new MultiFab(mg_mesh[mglev], 1, 1));
		    MultiFab& d = sigma_nd[i][mglev];
		    for (MFIter s_mfi(s); s_mfi.isValid(); ++s_mfi)
		    {
			d[s_mfi].copy(s[s_mfi], i, 0);
		    }
		}
		delete sigma_split.remove(mglev);
		sigma.set(mglev, &sigma_nd[0][mglev]);
	    }
	}
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    sigma_nd[i].set(mglev_max, &sigma[mglev_max]);
	}

	if (m_stencil == cross)
	{
	    sigma_node.resize(mglev_max + 1);
	    for (int mglev = 0; mglev <= mglev_max; mglev++)
	    {
	      if (for_fill_sync_reg == 0 || mglev == mglev_max)
	      {
		    BoxArray mesh = mg_mesh[mglev];
		    mesh.convert(IndexType(IntVect::TheNodeVector()));
		    sigma_node.set(mglev, new MultiFab(mesh, BL_SPACEDIM, 1));
		    sigma_node[mglev].setVal(1.0e50);
	      }
	    }

	    for (int mglev = 0; mglev <= mglev_max; mglev++)
	    {
		for ( int i = 0; i < BL_SPACEDIM; ++i )
		{
		    HG_TEST_NORM(sigma_nd[i][mglev], "build_sigma pre hgscon");
		}
		if (for_fill_sync_reg == 0 || mglev == mglev_max)
		{
		    const Real hxyz[BL_SPACEDIM] = { D_DECL(h[mglev][0],
							    h[mglev][1],
							    h[mglev][2]) };
		    for (MFIter s_mfi(sigma[mglev]); s_mfi.isValid(); ++s_mfi)
		    {
			const Box& scbox = sigma[mglev][s_mfi].box();
			const Box& snbox = sigma_node[mglev][s_mfi].box();
			Box reg = (for_fill_sync_reg > 0) ?
			    BoxLib::surroundingNodes(mg_mesh[mglev][s_mfi.index()]) :
			    Box(lev_interface[mglev].part_fine(s_mfi.index()));
			FORT_HGSCON(
			    sigma_node[mglev][s_mfi].dataPtr(),
			    DIMLIST(snbox),
			    D_DECL(sigma_nd[0][mglev][s_mfi].dataPtr(),
				   sigma_nd[1][mglev][s_mfi].dataPtr(),
				   sigma_nd[2][mglev][s_mfi].dataPtr()), DIMLIST(scbox),
			    DIMLIST(reg),
			    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]));
		    }
		    if (mglev < mglev_max)
		    {
			sigma_nd[0].remove(mglev);
			for (int i = 1; i < BL_SPACEDIM; i++)
			{
			    delete sigma_nd[i].remove(mglev);
			}
		    }
		    else
		    {
			for (int i = 0; i < BL_SPACEDIM; i++)
			{
			    sigma_nd[i].remove(mglev);
			}
		    }
		}
		HG_TEST_NORM(sigma_node[mglev], "build_sigma hgscon");
	    }
	}
    }

    cen.resize(mglev_max + 1);
    for (int mglev = 0; mglev <= mglev_max; mglev++)
    {
        if (for_fill_sync_reg == 0 || mglev == mglev_max)
        {
            cen.set(mglev, new MultiFab(corr[mglev].boxArray(), 1,
                    dest[lev_min].nGrow()));
            MultiFab& ctmp = cen[mglev];
            ctmp.setVal(0.0);

            if ( m_stencil == terrain || m_stencil == full)
            {
                for (MFIter c_mfi(ctmp); c_mfi.isValid(); ++c_mfi)
                {
                    const Box& cenbox = ctmp[c_mfi].box();
                    const Box& reg = 
                                 lev_interface[mglev].part_fine(c_mfi.index());
                    const Box& sigbox = sigma[mglev][c_mfi].box();
                    if ( m_stencil == terrain )
                    {
                        FORT_HGCEN_TERRAIN(ctmp[c_mfi].dataPtr(), DIMLIST(cenbox),
                                           sigma[mglev][c_mfi].dataPtr(), DIMLIST(sigbox),
                                           DIMLIST(reg));
                    }
                    else
                    {
                        FORT_HGCEN_FULL(ctmp[c_mfi].dataPtr(), DIMLIST(cenbox),
                                        sigma[mglev][c_mfi].dataPtr(), DIMLIST(sigbox),
                                        DIMLIST(reg));
                    }
                }
            }
            else
            {
                HG_TEST_NORM(sigma_node[mglev], "buildsigma");
                for (MFIter c_mfi(cen[mglev]); c_mfi.isValid(); ++c_mfi)
                {
                    const Box& cenbox = cen[mglev][c_mfi].box();
                    const Box& reg = 
                                  lev_interface[mglev].part_fine(c_mfi.index());
                    const Box& sigbox = sigma_node[mglev][c_mfi].box();
                    const int isRZ = getCoordSys();
                    FORT_HGCEN(cen[mglev][c_mfi].dataPtr(), DIMLIST(cenbox),
                               sigma_node[mglev][c_mfi].dataPtr(), DIMLIST(sigbox),
                               DIMLIST(reg), &isRZ);
                }
            }
            HG_TEST_NORM(ctmp, "buildsigma");
            clear_part_interface(ctmp, lev_interface[mglev]);
        }
    }
    
    if (m_stencil == cross)
    {
	mask.resize(mglev_max + 1);
        if (for_fill_sync_reg == 0)
        {
	    for (int mglev = 0; mglev <= mglev_max; mglev++)
	    {
		mask.set(mglev, new MultiFab(corr[mglev].boxArray(),
				             1, dest[lev_min].nGrow()));
		mask[mglev].setVal(0.0);
		for (MFIter m_mfi(mask[mglev]); m_mfi.isValid(); ++m_mfi)
		{
		    mask[mglev][m_mfi].setVal(1.0,lev_interface[mglev].part_fine(m_mfi.index()), 0);
		}
		HG_TEST_NORM(mask[mglev], "buildsigma 2");
		clear_part_interface(mask[mglev], lev_interface[mglev]);
	    }
        }
        else
        {
	    mask.set(mglev_max, new MultiFab(corr[mglev_max].boxArray(),
					     1, dest[lev_min].nGrow()));
            mask[mglev_max].setVal(0.0);
	    for (MFIter m_mfi(mask[mglev_max]); m_mfi.isValid(); ++m_mfi)
            {
                Box subbox(m_mfi.validbox());
		mask[mglev_max][m_mfi].setVal(1.0, subbox, 0);
            }
        }
    }
}


void
holy_grail_amr_multigrid::clear_hg_multi ()
{
    line_order.clear();
    line_after.clear();

    delete_sync_caches();

    delete cgwork.remove(0);
    delete cgwork.remove(1);
    delete cgwork.remove(2);
    cgwork.remove(3);
    cgwork.remove(4);
    cgwork.remove(5);
    delete cgwork.remove(6);
    delete cgwork.remove(7);

    if (m_stencil == terrain || m_stencil == full )
    {
	for (int mglev = 0; mglev <= mglev_max; mglev++)
	{
	    delete sigma.remove(mglev);
	}
    }
    else
    {
	for (int mglev = 0; mglev <= mglev_max; mglev++)
	{
	    delete sigma.remove(mglev);
	    if (sigma_node.size() && sigma_node.defined(mglev))
	    {
		delete sigma_node.remove(mglev);
	    }
	}
    }

    for (int mglev = 0; mglev <= mglev_max; mglev++)
    {
	delete cen.remove(mglev);
	if (mask.size() && mask.defined(mglev))
	{
	    delete mask.remove(mglev);
	}
    }

    delete [] h;
    if (source_owned)
    {
	for (int lev = lev_min; lev <= lev_max; lev++)
	{
	    if (source.defined(lev))
	    {
		delete source.remove(lev);
	    }
	}
    }

    amr_multigrid::clear_amr_multi();
}

void
holy_grail_amr_multigrid::sync_resid_clear ()
{
    line_order.clear();
    line_after.clear();

    delete sigma.remove(mglev_max);
    delete cen.remove(mglev_max);
    if (m_stencil == cross )
    {
	delete sigma_node.remove(mglev_max);
	delete mask.remove(mglev_max);
    }
    delete [] h;
    if ( source_owned )
    {
	for (int lev = lev_min; lev <= lev_max; lev++)
	{
	    if (source.defined(lev))
	    {
		delete source.remove(lev);
	    }
	}
    }

    amr_multigrid::clear_amr_multi();
}

bool
holy_grail_amr_multigrid::can_coarsen (const BoxArray& mesh,
                                       const Box&      domain) const
{
    int retval = 1;
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	retval &= ((domain.smallEnd(i)&1) == 0);
	retval &= ((domain.bigEnd(i)&1)   == 1);
	retval &= (domain.length(i) >= 4);
	for (int igrid = 0; igrid < mesh.size(); igrid++)
	{
	    retval &= (
		(mesh[igrid].smallEnd(i)&1) == 0 &&
		(mesh[igrid].bigEnd(i)&1)   == 1 &&
		(mesh[igrid].length(i) >= 4)
		);
	}
    }
    return retval != 0;
}

void
holy_grail_amr_multigrid::sync_interfaces ()
{
    for (int lev = lev_min+1; lev <= lev_max; lev++)
    {
	int mglev = ml_index[lev];
	int mgc = ml_index[lev-1];
	IntVect rat = mg_domain[mglev].length() / mg_domain[mgc].length();
	MultiFab& target = dest[lev];
	task_list tl;
	for (int iface = 0;
	     iface < lev_interface[mglev].nboxes(level_interface::FACEDIM);
	     iface++)
	{
            //
	    // Find a fine grid touching this face.
            //
	    int igrid = lev_interface[mglev].grid(level_interface::FACEDIM,
						  iface, 0);
	    if (igrid < 0)
		igrid = lev_interface[mglev].grid(level_interface::FACEDIM,
						  iface, 1);
	    const unsigned int geo =
		lev_interface[mglev].geo(level_interface::FACEDIM, iface);
            //
	    // Reject fine-fine interfaces and those without an interior fine grid
            //
	    const Box& nbox =
		lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	    if (   geo == level_interface::ALL
		|| igrid < 0
		|| lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
		continue;
	    tl.add_task(
		new task_interpolate_patch(tl, target, igrid, nbox,
					   dest[lev-1], rat,
					   new bilinear_interpolator(),
					   lev_interface[mgc]));
	}
	tl.execute("holy_grail_amr_multigrid::sync_interfaces");
    }
}

void
holy_grail_amr_multigrid::sync_periodic_interfaces ()
{
    for (int lev = lev_min+1; lev <= lev_max; lev++)
    {
	int mglev = ml_index[lev];
	int mgc = ml_index[lev-1];
	IntVect rat = mg_domain[mglev].length() / mg_domain[mgc].length();
	Box idomain = mg_domain[mglev];
	idomain.convert(type(dest[lev])).grow(-1);
	MultiFab& target = dest[lev];
	task_list tl;
	for (int iface = 0;
	     iface < lev_interface[mglev].nboxes(level_interface::FACEDIM);
	     iface++)
	{
            //
	    // Find a fine grid touching this face.
            //
	    int igrid =
		lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	    if (igrid < 0)
		igrid =
		    lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	    const unsigned int geo =
		lev_interface[mglev].geo(level_interface::FACEDIM, iface);
            //
	    // Use only exterior coarse-fine faces with an interior fine grid.
            //
	    const Box& nbox =
		lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	    if ( geo == level_interface::ALL
		 || igrid < 0
		 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
		continue;
	    if (idomain.intersects(nbox))
                continue;
	    tl.add_task(
		new task_interpolate_patch(tl, target, igrid, nbox,
					   dest[lev-1], rat,
					   new bilinear_interpolator(),
					   lev_interface[mgc]));
	}
	tl.execute("holy_grail_amr_multigrid::sync_periodic_interfaces");
    }
}

void
holy_grail_amr_multigrid::mg_restrict_level (int lto,
                                             int lfrom)
{
    IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
    HG_TEST_NORM( resid[lto], "mg_restrict_level: resid in" );
    HG_TEST_NORM( work[lfrom], "mg_restrict_level: work in" );
    if (get_amr_level(lto) >= 0)
    {
	restrict_level(
	    resid[lto], work[lfrom], rat,
	    bilinear_restrictor((integrate==0)?0:1, is_dense(m_stencil)),
	    lev_interface[lfrom], mg_boundary);
    }
    else
    {
	mg_restrict(lto, lfrom);
    }
    HG_TEST_NORM( resid[lto], "mg_restrict_level: resid out" );
    HG_TEST_NORM( work[lfrom], "mg_restrict_level: work out" );
}

void
holy_grail_amr_multigrid::mg_restrict (int lto,
                                       int lfrom)
{
    HG_TEST_NORM( work[lfrom], "mg_restrict 1");
    HG_TEST_NORM( resid[lto], "mg_restrict 11");
    fill_borders(work[lfrom], lev_interface[lfrom],
		 mg_boundary, -1, is_dense(m_stencil));
    const IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
    for (MFIter w_mfi(work[lfrom]); w_mfi.isValid(); ++w_mfi)
    {
	const Box& fbox = work[lfrom][w_mfi].box();
	const Box& cbox = resid[lto][w_mfi].box();
	const Box& creg = lev_interface[lto].part_fine(w_mfi.index());
	FORT_FANRST2(resid[lto][w_mfi].dataPtr(), DIMLIST(cbox), DIMLIST(creg),
                     work[lfrom][w_mfi].dataPtr(), DIMLIST(fbox),
                     D_DECL(rat[0], rat[1], rat[2]), 1, &integrate, 0, 0);
    }
    clear_part_interface(resid[lto], lev_interface[lto]);
    HG_TEST_NORM( resid[lto], "mg_restrict 21");
}

void
holy_grail_interpolator_dense::fill (FArrayBox&       patch,
				     const Box&       region,
				     const FArrayBox& cgr,
				     const Box&       cb,
				     const IntVect&   rat) const
{
    FORT_HGINTS_DENSE(
	patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	D_DECL(sigptr[0], sigptr[1], sigptr[2]),
	DIMLIST(sigbox),
	cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb),
	D_DECL(rat[0], rat[1], rat[2]));

}

void
holy_grail_interpolator::fill (FArrayBox&       patch,
			       const Box&       region,
			       const FArrayBox& cgr,
			       const Box&       cb,
			       const IntVect&   rat) const
{
    FORT_HGINTS(
	patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	sigptr, DIMLIST(sigbox),
	cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb),
	D_DECL(rat[0], rat[1], rat[2]));

}

void
holy_grail_amr_multigrid::mg_interpolate_level (int lto,
                                                int lfrom)
{
    if (get_amr_level(lfrom) >= 0)
    {
        //
	// general version---attempt to use special stencils for multilevel.
        //
	const int ltmp = lfrom + 1;
	MultiFab& target = work[ltmp];
	const IntVect rat =
	    mg_domain[ltmp].length() / mg_domain[lfrom].length();
	HG_TEST_NORM( corr[lfrom], "mg_interpolate_level");
	task_list tl;
	for (int igrid = 0; igrid < target.size(); igrid++)
	{
	    amr_interpolator* hgi;
	    if ( is_dense(m_stencil) )
	    {
		Real* sigptr[BL_SPACEDIM] = { D_DECL(0, 0, 0) };
		if ( is_local(sigma[ltmp], igrid) )
		{
		    for (int i = 0; i < BL_SPACEDIM; i++)
		    {
			sigptr[i] = sigma[ltmp][igrid].dataPtr(i);
		    }
		}
		const Box sigbox =
		    BoxLib::grow(sigma[ltmp].box(igrid),
			 sigma[ltmp].nGrow());
		BL_ASSERT( is_remote(sigma[ltmp], igrid)
			|| sigbox == sigma[ltmp][igrid].box());
		hgi = new holy_grail_interpolator_dense(sigptr, sigbox);
	    }
	    else
	    {
		Real* sigptr = 0;
		if ( is_local(sigma_node[ltmp], igrid ) )
		{
		    sigptr = sigma_node[ltmp][igrid].dataPtr();
		}
		const Box sigbox =
		    BoxLib::grow(sigma_node[ltmp].box(igrid),
			 sigma_node[ltmp].nGrow());
		BL_ASSERT(    is_remote(sigma_node[ltmp], igrid)
			   || sigbox == sigma_node[ltmp][igrid].box());
		// const Box& sigbox = sigma_node[ltmp][igrid].box();
		hgi = new holy_grail_interpolator(sigptr, sigbox);
	    }
	    tl.add_task(
		new task_interpolate_patch(tl, target, igrid, target.box(igrid),
					   corr[lfrom], rat,
					   hgi,
					   lev_interface[lfrom]));
	}
	tl.execute("holy_grail_amr_multigrid::mg_interpolate_level");
	HG_TEST_NORM( target, "mg_interpolate_level a");
	if (lto > ltmp)
	{
	    corr[ltmp].copy(target);
	    mg_interpolate_level(lto, ltmp);
	}
    }
    else
    {
        //
	// Multigrid interpolation, grids known to match up
	// special stencil needed for multigrid convergence.
        //
	const IntVect rat = mg_domain[lto].length()
	    / mg_domain[lfrom].length();
	for (MFIter w_mfi(work[lto]); w_mfi.isValid(); ++w_mfi)
	{
	    const Box& fbox = work[lto][w_mfi].box();
	    const Box& freg = w_mfi.validbox();
	    const Box& cbox = corr[lfrom][w_mfi].box();
	    const Box& creg = corr[lfrom].boxArray()[w_mfi.index()];
	    if ( m_stencil == terrain || m_stencil == full )
	    {
		const Box& sigbox = sigma[lto][w_mfi].box();
		FORT_HGINTS_DENSE(
		    work[lto][w_mfi].dataPtr(), DIMLIST(fbox), DIMLIST(freg),
		    D_DECL(sigma[lto][w_mfi].dataPtr(0),
			   sigma[lto][w_mfi].dataPtr(1),
			   sigma[lto][w_mfi].dataPtr(2)), DIMLIST(sigbox),
		    corr[lfrom][w_mfi].dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		    D_DECL(rat[0], rat[1], rat[2]));
	    }
	    else
	    {
		const Box& sigbox = sigma_node[lto][w_mfi].box();
		FORT_HGINTS(
		    work[lto][w_mfi].dataPtr(), DIMLIST(fbox), DIMLIST(freg),
		    sigma_node[lto][w_mfi].dataPtr(), DIMLIST(sigbox),
		    corr[lfrom][w_mfi].dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		    D_DECL(rat[0], rat[1], rat[2]));
	    }
	}
	HG_TEST_NORM( work[lto], "mg_interpolate_level");
    }
}
