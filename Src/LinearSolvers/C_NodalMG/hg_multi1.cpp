
#include "hg_multi.H"

#ifdef BL_FORT_USE_UNDERSCORE
#define   FORT_HGSRST		hgsrst_
#define   FORT_HGSCON		hgscon_
#define   FORT_HGCEN		hgcen_
#define   FORT_HGCEN_FULL	hgcen_full_
#define   FORT_HGCEN_NO_SIGMA   hgcen_no_sigma_
#define   FORT_HGCEN_TERRAIN    hgcen_terrain_
#define   FORT_HGINTS		hgints_
#define   FORT_HGINTS_NO_SIGMA  hgints_no_sigma_
#define   FORT_FACRST1		acrst1_
#define   FORT_FANRST2		anrst2_
#else
#define   FORT_HGSRST		HGSRST
#define   FORT_HGSCON		HGSCON
#define   FORT_HGCEN		HGCEN
#define   FORT_HGCEN_FULL	HGCEN_FULL
#define   FORT_HGCEN_NO_SIGMA   HGCEN_NO_SIGMA
#define   FORT_HGCEN_TERRAIN	HGCEN_TERRAIN
#define   FORT_HGINTS		HGINTS
#define   FORT_HGINTS_NO_SIGMA  HGINTS_NO_SIGMA
#define   FORT_FACRST1		ACRST1
#define   FORT_FANRST2		ANRST2
#endif

extern "C" 
{
    void FORT_FACRST1        (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);
    void FORT_FANRST2        (Real*, intS, intS, const Real*, intS, intRS, const int&, const int*, const int*, const int*);

    void FORT_HGSRST         (RealPS, intS, intS, CRealPS, intS, intRS);
    void FORT_HGINTS         (Real*, intS, intS, Real*, intS, const Real*, intS, intS, intRS);
    void FORT_HGCEN_TERRAIN  (Real*, intS, Real*, intS, intS);
#if BL_SPACEDIM == 2
    void FORT_HGCEN_FULL     (Real*, intS, CRealPS, intS, intS, CRealPS, const int*, const int*);
#endif
    void FORT_HGCEN_NO_SIGMA (Real*, intS, RealPS, intS, intS, CRealPS);
    void FORT_HGINTS_NO_SIGMA(Real*, intS, intS, CRealPS, intS, const Real*, intS, intS, intRS);
    void FORT_HGSCON         (Real*, intS, RealPS, intS, intS, CRealPS);
    void FORT_HGCEN          (Real*, intS, Real*, intS, intS);
}

class task_interpolate_patch : public task
{
public:
    task_interpolate_patch(MultiFab& dmf_, int dgrid_, const Box& dbx_,
	const MultiFab& smf_, const IntVect& rat_,
	const amr_interpolator_class& interp_, const level_interface& lev_interface_)
	: dmf(dmf_), dgrid(dgrid_), dbx(dbx_),
	smf(smf_), rat(rat_), interp(interp_), lev_interface(lev_interface_)
    {
	assert(dbx.sameType(dmf[dgrid].box()));
	const Box cb = interp.box(dbx, rat);
	tf = new task_fill_patch( cb, smf, lev_interface, 0, -1, -1);
    }
    virtual bool ready()
    {
	BoxLib::Abort( "FIXME task_interpolate_patch::ready" );
	if ( tf->ready() )
	{
	    interp.fill(dmf[dgrid], dbx, tf->fab(), tf->fab().box(), rat);
	    return true;
	}
	return false;
    }
    virtual ~task_interpolate_patch()
    {
	delete tf;
    }
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	return tf->init(sno, comm);
    }
private:
    task_fab* tf;
    MultiFab& dmf;
    const int dgrid;
    const Box dbx;
    const MultiFab& smf;
    const IntVect rat;
    const amr_interpolator_class& interp;
    const level_interface& lev_interface;
};

void holy_grail_amr_multigrid::alloc(PArray<MultiFab>& Dest, PArray<MultiFab>& Source, PArray<MultiFab>& Coarse_source, PArray<MultiFab>& Sigma, Real H[], int Lev_min, int Lev_max)
{
    assert(Dest.length() > Lev_max);
    assert(Dest[Lev_min].nGrow() == 1);
    
    if (Source.ready()) 
    {
	source_owned = false;
	amr_multigrid::alloc(Dest, Source, Coarse_source, Lev_min, Lev_max);
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
	amr_multigrid::alloc(Dest, Src,    Coarse_source, Lev_min, Lev_max);
    }
    
    h = new Real[mglev_max + 1][BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	h[mglev_max][i] = H[i];
	for (int mglev = mglev_max - 1; mglev >= 0; mglev--) 
	{
	    int rat = mg_domain[mglev+1].length(i) / mg_domain[mglev].length(i);
	    h[mglev][i] = rat * h[mglev+1][i];
	}
    }
    
    
    build_sigma(Sigma);
    
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
    
    assert(cgwork[3].nGrow() == ib &&
	cgwork[4].nGrow() == ib &&
	cgwork[5].nGrow() == ib);
    
    for (MultiFabIterator g_mfi(cgwork[7]); g_mfi.isValid(); ++g_mfi)
    {
	FArrayBox& gtmp = g_mfi();
	const Box& valid = g_mfi.validbox();
	gtmp.setVal(0.0);
	gtmp.setVal(1.0, valid, 0);
	Box b = bdryLo(valid, 1);
	gtmp.setVal(0.5, b, 0);
	b = bdryHi(valid, 1);
	gtmp.setVal(0.5, b, 0);
	b = bdryLo(valid, 0);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, bdryLo(b, 1), 0);
	gtmp.setVal(0.25, bdryHi(b, 1), 0);
	b = bdryHi(valid, 0);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, bdryLo(b, 1), 0);
	gtmp.setVal(0.25, bdryHi(b, 1), 0);
#if (BL_SPACEDIM == 3)
	b = bdryLo(valid, 2);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, bdryLo(b, 0), 0);
	gtmp.setVal(0.25, bdryHi(b, 0), 0);
	Box bb = bdryLo(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, bdryHi(bb, 0), 0);
	bb = bdryHi(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, bdryHi(bb, 0), 0);
	b = bdryHi(valid, 2);
	gtmp.setVal(0.5, b, 0);
	gtmp.setVal(0.25, bdryLo(b, 0), 0);
	gtmp.setVal(0.25, bdryHi(b, 0), 0);
	bb = bdryLo(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, bdryHi(bb, 0), 0);
	bb = bdryHi(b, 1);
	gtmp.setVal(0.25, bb, 0);
	gtmp.setVal(0.125, bdryLo(bb, 0), 0);
	gtmp.setVal(0.125, bdryHi(bb, 0), 0);
#endif
    }
    
    singular = false;
    if (mg_boundary->singular()) 
    {
	long sng = 0;
	for (int i = 0; i < mg_mesh[0].length(); i++) 
	{
	    sng += mg_mesh[0][i].numPts();
	}
	singular = (sng == mg_domain[0].numPts());
    }
    
    if (m_hg_terrain)
	integrate = 1;
}


void holy_grail_sigma_restrictor_class::fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const
{
    assert(patch.box().cellCentered());
    assert(rat[0] == 2 && rat[1] == 2 ||
	rat[0] == 2 && rat[1] == 1 ||
	rat[0] == 1 && rat[1] == 2);
    
    if (m_hg_terrain)
    {
	FORT_HGSRST(D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
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
	patch.mult((Real) rat[1] / rat[0], region, 0, 1);
	patch.mult((Real) rat[0] / rat[1], region, 1, 1);
	// component 2 remains unchanged
#else
	FORT_FACRST1(patch.dataPtr(BL_SPACEDIM+1),
	    DIMLIST(patch.box()),
	    DIMLIST(region),
	    fgr.dataPtr(BL_SPACEDIM+1),
	    DIMLIST(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]), 1, &integ, 0, 0);
	patch.mult((Real) rat[1] * rat[2] / rat[0], region, 0, 1);
	patch.mult((Real) rat[0] * rat[2] / rat[1], region, 1, 1);
	patch.mult((Real) rat[0] * rat[1] / rat[2], region, 2, 1);
	patch.mult((Real) rat[1],                   region, 3, 1);
	patch.mult((Real) rat[0],                   region, 4, 1);
#endif
	
    }
    else
    {
	
	if (fgr.nComp() == 1) 
	{
	    FORT_HGSRST(D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
		DIMLIST(patch.box()),
		DIMLIST(region),
		D_DECL(fgr.dataPtr(), fgr.dataPtr(), fgr.dataPtr()),
		DIMLIST(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
	}
	else 
	{
	    FORT_HGSRST(D_DECL(patch.dataPtr(0), patch.dataPtr(1), patch.dataPtr(2)),
		DIMLIST(patch.box()),
		DIMLIST(region),
		D_DECL(fgr.dataPtr(0), fgr.dataPtr(1), fgr.dataPtr(2)),
		DIMLIST(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
	}
	
    }
}

void holy_grail_amr_multigrid::build_sigma(PArray<MultiFab>& Sigma)
{
    
    if (m_hg_terrain)
    {
	
	// For terrain stencils we have as many sigma arrays passed as
	// arguments and used at the lev_interface as we build for internal
	// multigrid purposes.  This simplifies handling as we do not
	// need to maintain separate arrays for different purposes.
	
	const int ncomp = 2 * BL_SPACEDIM - 1;
	
	sigma.resize(mglev_max+1);
	
	for (int mglev = 0; mglev <= mglev_max; mglev++) 
	{
	    sigma.set(mglev, new MultiFab(mg_mesh[mglev], ncomp, 1));
	    MultiFab& target = sigma[mglev];
	    target.setVal(1.e20);
	    int lev;
	    if ((lev = get_amr_level(mglev)) >= 0) 
	    {
		MultiFab& s_in = Sigma[lev];
		for (MultiFabIterator s_mfi(s_in); s_mfi.isValid(); ++s_mfi)
		{
		    DependentMultiFabIterator t_dmfi(s_mfi, target);
		    t_dmfi->copy(s_mfi(), s_mfi.validbox(), 0, t_dmfi.validbox(), 0, ncomp);
		}
	    }
	}
	
	for (int mglev = mglev_max; mglev > 0; mglev--) 
	{
	    IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
	    restrict_level(sigma[mglev-1], sigma[mglev], rat, holy_grail_sigma_restrictor_class(), level_interface(), 0);
	}
	for (int mglev = 0; mglev <= mglev_max; mglev++) 
	{
	    fill_borders(sigma[mglev], lev_interface[mglev], boundary.terrain_sigma(), -1, m_hg_terrain);
	}
	
    }
    else
    {
	
	// Intended functionality:  sigma_split exists only at coarser levels,
	// since only after coarsening is sigma different in different directions.
	// sigma exists at all levels, and is intended for use on fine grids
	// and at lev_interface points, where all components are the same.  To
	// save storage it is aliased to the first component of sigma_split
	// on all but the finest level.
	
	// sigma_split replaced by sigma_nd in more recent version, used
	// only as a local variable here
	
	PArray<MultiFab> sigma_split;
	sigma_split.resize(mglev_max);
	for (int mglev = 0; mglev < mglev_max; mglev++) 
	{
	    sigma_split.set(mglev, new MultiFab(mg_mesh[mglev], BL_SPACEDIM, 1));
	}
	sigma.resize(mglev_max + 1);
	sigma.set(mglev_max, new MultiFab(mg_mesh[mglev_max], 1, 1));
	
	// Level project:
	// Any value can fill values in the border cells that fill_borders
	// will not touch---those touching coarser grids.  The values in these
	// cells will only be seen by the interpolation, and the quantity
	// being interpolated will always be zero at these edges, but we
	// should insure that no NaN's or other garbage is there that could
	// cause a floating point fault.
	
	// Sync project:
	// Ghost values will be seen by multilevel interpolation, so put
	// a huge value in ghost cells so that coarse-fine lev_interface
	// interpolation will linear, as finite element derivation requires.
	
	for (int mglev = 0; mglev < mglev_max; mglev++) 
	{
	    MultiFab& target = sigma_split[mglev];
	    target.setVal(1.e20);
	    int lev;
	    if ((lev = get_amr_level(mglev)) >= 0) 
	    {
		MultiFab& s_comp = Sigma[lev];
		for (int i = 0; i < BL_SPACEDIM; i++) 
		{
		    for (MultiFabIterator s_mfi(s_comp); s_mfi.isValid(); ++s_mfi)
		    {
			DependentMultiFabIterator t_dmfi(s_mfi, target);
			t_dmfi->copy(s_mfi(), s_mfi.validbox(), 0, t_dmfi.validbox(), i, 1);
		    }
		}
	    }
	}
	
	sigma[mglev_max].setVal(1.e20);
	for (MultiFabIterator S_mfi(Sigma[lev_max]); S_mfi.isValid(); ++S_mfi)
	{
	    DependentMultiFabIterator s_dmfi(S_mfi, sigma[mglev_max]);
	    s_dmfi->copy(S_mfi(), mg_mesh[mglev_max][S_mfi.index()], 0, mg_mesh[mglev_max][S_mfi.index()], 0, 1);
	}
	
	if (mglev_max > 0) 
	{
	    IntVect rat = mg_domain[mglev_max].length() / mg_domain[mglev_max-1].length();
	    restrict_level(sigma_split[mglev_max-1], sigma[mglev_max], rat, holy_grail_sigma_restrictor_class(), level_interface(), 0);
	}
	fill_borders(sigma[mglev_max], lev_interface[mglev_max], boundary.scalar(), -1, m_hg_terrain);
	for (int mglev = mglev_max - 1; mglev > 0; mglev--) 
	{
	    IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
	    restrict_level(sigma_split[mglev-1], sigma_split[mglev], rat, holy_grail_sigma_restrictor_class(), level_interface(), 0);
	}
	for (int mglev = 0; mglev < mglev_max; mglev++) 
	{
	    HG_TEST_NORM(sigma_split[mglev], "build_sigma 0");
	    fill_borders(sigma_split[mglev], lev_interface[mglev], boundary.scalar(), -1, m_hg_terrain);
	    HG_TEST_NORM(sigma_split[mglev], "build_sigma");
	}
	
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    sigma_nd[i].resize(mglev_max + 1);
	}
	
	for (int mglev = 0; mglev < mglev_max; mglev++) 
	{
	    MultiFab& s = sigma_split[mglev];
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
		sigma_nd[i].set(mglev, new MultiFab(mg_mesh[mglev], 1, 1));
		MultiFab& d = sigma_nd[i][mglev];
		for (MultiFabIterator s_mfi(s); s_mfi.isValid(); ++s_mfi)
		{
		    DependentMultiFabIterator d_dmfi(s_mfi, d);
		    d_dmfi->copy(s_mfi(), i, 0);
		}
	    }
	    delete sigma_split.remove(mglev);
	    sigma.set(mglev, &sigma_nd[0][mglev]);
	}
	
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    sigma_nd[i].set(mglev_max, &sigma[mglev_max]);
	}
	
	if (m_hg_cross_stencil)
	{
#if	BL_SPACEDIM==3
	    sigma_node.resize(mglev_max + 1);
	    for (int mglev = 0; mglev <= mglev_max; mglev++) 
	    {
		BoxArray mesh = mg_mesh[mglev];
		mesh.convert(IndexType(IntVect::TheNodeVector()));
		sigma_node.set(mglev, new MultiFab(mesh, BL_SPACEDIM, 1));
		sigma_node[mglev].setVal(1.e20);
	    }
	    
	    for (int mglev = 0; mglev <= mglev_max; mglev++) 
	    {
		const Real hxyz[BL_SPACEDIM] = { D_DECL(h[mglev][0], h[mglev][1], h[mglev][2]) };
		for (MultiFabIterator s_mfi(sigma[mglev]); s_mfi.isValid(); ++s_mfi)
		{
		    DependentMultiFabIterator sn_dmfi(s_mfi, sigma_node[mglev]);
		    DependentMultiFabIterator D_DECL(sn0(s_mfi, sigma_nd[0][mglev]),
			sn1(s_mfi, sigma_nd[1][mglev]),
			sn2(s_mfi, sigma_nd[2][mglev]));
		    const Box& scbox = s_mfi->box();
		    const Box& snbox = sn_dmfi->box();
		    const Box& reg = lev_interface[mglev].part_fine(s_mfi.index());
		    FORT_HGSCON(sn_dmfi->dataPtr(),
			DIMLIST(snbox),
			D_DECL(sn0->dataPtr(), sn1->dataPtr(), sn2->dataPtr()),
			DIMLIST(scbox),
			DIMLIST(reg),
			D_DECL(&hxyz[0], &hxyz[1], &hxyz[2])
			);
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
#endif
	}
	
    }
    
    cen.resize(mglev_max + 1);
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	cen.set(mglev, new MultiFab(corr[mglev].boxArray(), 1, dest[lev_min].nGrow()));
	MultiFab& ctmp = cen[mglev];
	ctmp.setVal(0.0);
	
	if (m_hg_terrain)
	{
	    
	    for (MultiFabIterator c_mfi(ctmp); c_mfi.isValid(); ++c_mfi)
	    {
		DependentMultiFabIterator s_dmfi(c_mfi, sigma[mglev]);
		const Box& cenbox = c_mfi->box();
		const Box& reg = lev_interface[mglev].part_fine(c_mfi.index());
		const Box& sigbox = s_dmfi->box();
		FORT_HGCEN_TERRAIN(c_mfi->dataPtr(), DIMLIST(cenbox), s_dmfi->dataPtr(), DIMLIST(sigbox), DIMLIST(reg));
	    }
	}
	else if (m_hg_full_stencil)
	{
#if BL_SPACEDIM != 3
	    const Real hxyz[BL_SPACEDIM] = { D_DECL( h[mglev][0], h[mglev][1], h[mglev][2] ) };
	    for (MultiFabIterator c_mfi(cen[mglev]); c_mfi.isValid(); ++c_mfi)
	    {
		const Box& cenbox = c_mfi->box();
		const Box& reg = lev_interface[mglev].part_fine(c_mfi.index());
		DependentMultiFabIterator D_DECL(sn0(c_mfi, sigma_nd[0][mglev]),
		    sn1(c_mfi, sigma_nd[1][mglev]),
		    sn2(c_mfi, sigma_nd[2][mglev]));
		const Box& sigbox = sn0->box();
#if (BL_SPACEDIM == 2)
		const int imax = mg_domain[mglev].bigEnd(0) + 1;
		const int isRZ = IsRZ();
		FORT_HGCEN_FULL(c_mfi->dataPtr(), DIMLIST(cenbox),
		    D_DECL(sn0->dataPtr(), sn1->dataPtr(), sn2->dataPtr()), DIMLIST(sigbox),
		    DIMLIST(reg),
		    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]), &isRZ, &imax);
#else
		FORT_HGCEN_NO_SIGMA(c_mfi->dataPtr(), DIMLIST(cenbox),
		    D_DECL(sn0->dataPtr(), sn1->dataPtr(),sn2->dataPtr()),DIMLIST(sigbox),
		    DIMLIST(reg),
		    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]));
#endif
	    }
#endif
	}
	else
	{
	    HG_TEST_NORM(sigma_node[mglev],"buildsigma");
	    for (MultiFabIterator c_mfi(cen[mglev]); c_mfi.isValid(); ++c_mfi)
	    {
		const Box& cenbox = c_mfi->box();
		const Box& reg = lev_interface[mglev].part_fine(c_mfi.index());
#if BL_SPACEDIM == 2
		const Real hxyz[BL_SPACEDIM] = { D_DECL( h[mglev][0], h[mglev][1], h[mglev][2] ) };
		DependentMultiFabIterator D_DECL(sn0(c_mfi, sigma_nd[0][mglev]),
		    sn1(c_mfi, sigma_nd[1][mglev]),
		    sn2(c_mfi, sigma_nd[2][mglev]));
		const Box& sigbox = sn0->box();
		FORT_HGCEN_NO_SIGMA(c_mfi->dataPtr(), DIMLIST(cenbox),
		    D_DECL(sn0->dataPtr(),sn1->dataPtr(),sn2->dataPtr()),DIMLIST(sigbox),
		    DIMLIST(reg),
		    D_DECL(&hxyz[0], &hxyz[1], &hxyz[2]));	    
#else
		DependentMultiFabIterator sn_dmfi(c_mfi, sigma_node[mglev]);
		const Box& sigbox = sn_dmfi->box();
		FORT_HGCEN(c_mfi->dataPtr(), DIMLIST(cenbox),
		    sn_dmfi->dataPtr(), DIMLIST(sigbox),
		    DIMLIST(reg));
#endif
	    }
	}	
	HG_TEST_NORM(ctmp, "buildsigma");
	clear_part_interface(ctmp, lev_interface[mglev]);
    }
    
    if (m_hg_cross_stencil)
    {
#if BL_SPACEDIM==3
	mask.resize(mglev_max + 1);
	for (int mglev = 0; mglev <= mglev_max; mglev++) 
	{
	    mask.set(mglev, new MultiFab(corr[mglev].boxArray(), 1, dest[lev_min].nGrow()));
	    MultiFab& mtmp = mask[mglev];
	    mtmp.setVal(0.0);
	    for (MultiFabIterator m_mfi(mtmp); m_mfi.isValid(); ++m_mfi)
	    {
		m_mfi->setVal(1.0, lev_interface[mglev].part_fine(m_mfi.index()), 0);
	    }
	    HG_TEST_NORM(mtmp, "buildsigma 2");
	    clear_part_interface(mtmp, lev_interface[mglev]);
	}
#endif
    }
}

void holy_grail_amr_multigrid::clear()
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
    
    if (m_hg_terrain)
    {
	for (int mglev = 0; mglev <= mglev_max; mglev++) 
	{
	    delete sigma.remove(mglev);
	}
    }
    else if (m_hg_full_stencil)
    {
	
	delete sigma.remove(mglev_max);
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    sigma_nd[i].remove(mglev_max);
	}
	for (int mglev = 0; mglev < mglev_max; mglev++) 
	{
	    sigma.remove(mglev);
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
		delete sigma_nd[i].remove(mglev);
	    }
	}
    }
    else
    {
	for (int mglev = 0; mglev <= mglev_max; mglev++) 
	{
	    delete sigma.remove(mglev);
	    if (sigma_node.ready() && sigma_node.defined(mglev))
		delete sigma_node.remove(mglev);
#if BL_SPACEDIM == 2
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
#endif	
	}
    }
//#endif
    
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	delete cen.remove(mglev);
	if (mask.ready() && mask.defined(mglev))
	    delete mask.remove(mglev);
    }
    
    delete [] h;
    if (source_owned) 
    {
	for (int lev = lev_min; lev <= lev_max; lev++) 
	{
	    if (source.defined(lev)) delete source.remove(lev);
	}
    }
    
    amr_multigrid::clear();
}

bool holy_grail_amr_multigrid::can_coarsen(const BoxArray& mesh, const Box& domain) const
{
    int retval = 1;
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	retval &= ((domain.smallEnd(i)&1) == 0);
	retval &= ((domain.bigEnd(i)&1)   == 1);
	retval &= (domain.length(i) >= 4);
	for (int igrid = 0; igrid < mesh.length(); igrid++) 
	{
	    retval &= ((mesh[igrid].smallEnd(i)&1) == 0 &&
		(mesh[igrid].bigEnd(i)&1)   == 1 &&
		(mesh[igrid].length(i) >= 4));
	}
    }
    return retval != 0;
}

void holy_grail_amr_multigrid::sync_interfaces()
{
    for (int lev = lev_min+1; lev <= lev_max; lev++) 
    {
	int mglev = ml_index[lev];
	int mgc = ml_index[lev-1];
	IntVect rat = mg_domain[mglev].length() / mg_domain[mgc].length();
	MultiFab& target = dest[lev];
	task_list tl;
	for (int iface = 0; iface < lev_interface[mglev].nboxes(level_interface::FACEDIM); iface++) 
	{
	    // find a fine grid touching this face
	    int igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	    if (igrid < 0)
		igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	    const unsigned int geo = lev_interface[mglev].geo(level_interface::FACEDIM, iface);
	    // reject fine-fine interfaces and those without an interior fine grid
	    const Box& nbox = lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
		continue;
	    tl.add_task(
		new task_interpolate_patch(target, igrid, nbox, dest[lev-1], rat, bilinear_interpolator_class(), lev_interface[mgc])
		);
	}
	tl.execute();
    }
}

void holy_grail_amr_multigrid::sync_periodic_interfaces()
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
	for (int iface = 0; iface < lev_interface[mglev].nboxes(level_interface::FACEDIM); iface++) 
	{
	    // find a fine grid touching this face
	    int igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	    if (igrid < 0)
		igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	    const unsigned int geo = lev_interface[mglev].geo(level_interface::FACEDIM, iface);
	    // use only exterior coarse-fine faces with an interior fine grid
	    const Box& nbox = lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
		continue;
	    if ( idomain.intersects(nbox) ) continue;
	    tl.add_task(
		new task_interpolate_patch(target ,igrid, nbox, dest[lev-1], rat, bilinear_interpolator_class(), lev_interface[mgc])
	    );
	}
	tl.execute();
    }
}

void holy_grail_amr_multigrid::mg_restrict_level(int lto, int lfrom)
{
    IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
    if (get_amr_level(lto) >= 0) 
    {
	restrict_level(resid[lto], work[lfrom], rat, bilinear_restrictor_class( (integrate==0)?0:1, m_hg_terrain), lev_interface[lfrom], mg_boundary);
    }
    else 
    {
	mg_restrict(lto, lfrom);
    }
}

void holy_grail_amr_multigrid::mg_restrict(int lto, int lfrom)
{
    fill_borders(work[lfrom], lev_interface[lfrom], mg_boundary, -1, m_hg_terrain);
    const IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
    for (MultiFabIterator w_mfi(work[lfrom]); w_mfi.isValid(); ++w_mfi)
    {
	DependentMultiFabIterator r_dmfi(w_mfi, resid[lto]);
	const Box& fbox = w_mfi->box();
	const Box& cbox = r_dmfi->box();
	const Box& creg = lev_interface[lto].part_fine(w_mfi.index());
	FORT_FANRST2(r_dmfi->dataPtr(), DIMLIST(cbox), DIMLIST(creg),
	    w_mfi->dataPtr(), DIMLIST(fbox),
	    D_DECL(rat[0], rat[1], rat[2]), 1, &integrate, 0, 0);
    }
    clear_part_interface(resid[lto], lev_interface[lto]);
}

void holy_grail_interpolator_class_not_cross::fill(FArrayBox& patch, const Box& region, const FArrayBox& cgr, const Box& cb, const IntVect& rat) const
{
    FORT_HGINTS_NO_SIGMA(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	D_DECL(sigptr[0], sigptr[1], sigptr[2]),
	DIMLIST(sigbox),
	cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb),
	D_DECL(rat[0], rat[1], rat[2]));
    
}

void holy_grail_interpolator_class::fill(FArrayBox& patch,
					 const Box& region,
					 const FArrayBox& cgr,
					 const Box& cb,
					 const IntVect& rat) const
{
    FORT_HGINTS(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	sigptr, DIMLIST(sigbox),
	cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb),
	D_DECL(rat[0], rat[1], rat[2]));
    
}

void holy_grail_amr_multigrid::mg_interpolate_level(int lto, int lfrom)
{
    if (get_amr_level(lfrom) >= 0) 
    {
	// general version---attempt to use special stencils for multilevel
	const int ltmp = lfrom + 1;
	MultiFab& target = work[ltmp];
	const IntVect rat = mg_domain[ltmp].length() / mg_domain[lfrom].length();
	task_list tl;
	for (int igrid = 0; igrid < target.length(); igrid++) 
	{
	    amr_interpolator_class* hgi;
	    if (m_hg_terrain)
	    {
		Real* sigptr[BL_SPACEDIM];
		for (int i = 0; i < BL_SPACEDIM; i++) 
		{
		    sigptr[i] = sigma[ltmp][igrid].dataPtr(i);
		}
		const Box& sigbox = sigma[ltmp][igrid].box();
		hgi = new holy_grail_interpolator_class_not_cross(sigptr, sigbox);
	    }
	    else if (m_hg_full_stencil)
	    {
		Real* sigptr[BL_SPACEDIM];
		for (int i = 0; i < BL_SPACEDIM; i++) 
		{
		    sigptr[i] = sigma_nd[i][ltmp][igrid].dataPtr();
		}
		const Box& sigbox = sigma[ltmp][igrid].box();
		hgi = new holy_grail_interpolator_class_not_cross(sigptr, sigbox);
	    }
	    else
	    {
#if BL_SPACEDIM != 3
		Real* sigptr[BL_SPACEDIM];
		for (int i = 0; i < BL_SPACEDIM; i++) 
		{
		    sigptr[i] = sigma_nd[i][ltmp][igrid].dataPtr();
		}
		const Box& sigbox = sigma_nd[0][ltmp][igrid].box();
		hgi = new holy_grail_interpolator_class_not_cross(sigptr, sigbox);
#else
		Real* sigptr = sigma_node[ltmp][igrid].dataPtr();
		const Box& sigbox = sigma_node[ltmp][igrid].box();
		hgi = new holy_grail_interpolator_class(sigptr, sigbox);
#endif
	    }
	    // BUG --  must keep a clone of the interpolator
	    tl.add_task(
		new task_interpolate_patch(target ,igrid, target.box(igrid), corr[lfrom], rat, *hgi, lev_interface[lfrom])
		);
	    delete hgi;
	}
	tl.execute();
	if (lto > ltmp) 
	{
	    corr[ltmp].copy(target);
	    mg_interpolate_level(lto, ltmp);
	}
    }
    else 
    {
	// multigrid interpolation, grids known to match up
	// special stencil needed for multigrid convergence
	const IntVect rat = mg_domain[lto].length() / mg_domain[lfrom].length();
	for (MultiFabIterator w_mfi(work[lto]); w_mfi.isValid(); ++w_mfi)
	{
	    DependentMultiFabIterator c_dmfi(w_mfi, corr[lfrom]);
	    const Box& fbox = w_mfi->box();
	    const Box& freg = w_mfi.validbox();
	    const Box& cbox = c_dmfi->box();
	    const Box& creg = c_dmfi.validbox();
	    if (m_hg_terrain)
	    {
		DependentMultiFabIterator sn_dmfi(w_mfi, sigma[lto]);
		const Box& sigbox = sn_dmfi->box();
		FORT_HGINTS_NO_SIGMA(w_mfi->dataPtr(), DIMLIST(fbox), DIMLIST(freg),
		    D_DECL(sn_dmfi->dataPtr(0), sn_dmfi->dataPtr(1), sn_dmfi->dataPtr(2)), 
		    DIMLIST(sigbox),
		    c_dmfi->dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		    D_DECL(rat[0], rat[1], rat[2]));
	    }
	    else if (m_hg_full_stencil)
	    {
#if (BL_SPACEDIM != 3)
		DependentMultiFabIterator D_DECL(s0_dmfi(w_mfi, sigma_nd[0][lto]),
						 s1_dmfi(w_mfi, sigma_nd[1][lto]),
						 s2_dmfi(w_mfi, sigma_nd[2][lto]));
		const Box& sigbox = s0_dmfi->box();
		FORT_HGINTS_NO_SIGMA(w_mfi->dataPtr(), DIMLIST(fbox), DIMLIST(freg),
		    D_DECL(s0_dmfi->dataPtr(), s1_dmfi->dataPtr(), s2_dmfi->dataPtr()),
		    DIMLIST(sigbox),
		    c_dmfi->dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		    D_DECL(rat[0], rat[1], rat[2]));
#endif
	    }
	    else
	    {
#if (BL_SPACEDIM != 3)
		DependentMultiFabIterator D_DECL(s0_dmfi(w_mfi, sigma_nd[0][lto]),
						 s1_dmfi(w_mfi, sigma_nd[1][lto]),
						 s2_dmfi(w_mfi, sigma_nd[2][lto]));
		const Box& sigbox = s0_dmfi->box();
		FORT_HGINTS_NO_SIGMA(w_mfi->dataPtr(), DIMLIST(fbox), DIMLIST(freg),
					  D_DECL(s0_dmfi->dataPtr(),
						 s1_dmfi->dataPtr(),
						 s2_dmfi->dataPtr()), DIMLIST(sigbox),
					  c_dmfi->dataPtr(), DIMLIST(cbox), DIMLIST(creg),
					  D_DECL(rat[0], rat[1], rat[2]));
#else
		DependentMultiFabIterator sn_dmfi(w_mfi, sigma_node[lto]);
		const Box& sigbox = sn_dmfi->box();
		FORT_HGINTS(w_mfi->dataPtr(), DIMLIST(fbox), DIMLIST(freg),
		    sn_dmfi->dataPtr(),
		    DIMLIST(sigbox),
		    c_dmfi->dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		    D_DECL(rat[0], rat[1], rat[2]));
#endif
	    }
	}
    }
}
