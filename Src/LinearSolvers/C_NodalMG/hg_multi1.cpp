
#include "hg_multi.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGSRST     hgsrst_
#  define   FORT_HGSCON     hgscon_
#  define   FORT_HGCEN      hgcen_
#  define   FORT_HGINTS     hgints_
#  define   FORT_FACRST1    acrst1_
#  define   FORT_FANRST2    anrst2_
#  define   FORT_FANINT2    anint2_
#else
#  define   FORT_HGSRST     HGSRST
#  define   FORT_HGSCON     HGSCON
#  define   FORT_HGCEN      HGCEN
#  define   FORT_HGINTS     HGINTS
#  define   FORT_FACRST1    ACRST1
#  define   FORT_FANRST2    ANRST2
#  define   FORT_FANINT2    ANINT2
#endif

extern "C" 
{
    
#if (BL_SPACEDIM == 1)
#error not relevant
#elif (BL_SPACEDIM == 2 || BL_SPACEDIM == 3)
#  ifdef HG_TERRAIN
    void FORT_FACRST1(Real*, intS, intS, Real*, intS, intRS, const int&);
    void FORT_HGSRST(RealPS, intS, intS, RealPS, intS, intRS);
    void FORT_HGCEN(Real*, intS, Real*, intS, intS);
    void FORT_HGINTS(Real*, intS, intS, RealPS, intS, const Real*, intS, intS, intRS);
#  elif (! defined HG_CONSTANT)
    void FORT_HGSRST(RealPS, intS, intS, CRealPS, intS, intRS);
#    ifndef HG_SIGMA_NODE
    void FORT_HGCEN(Real*, intS, RealPS, intS, intS, RealRS,
	const int&, const int&);
    void FORT_HGINTS(Real*, intS, intS, CRealPS, intS, const Real*, intS, intS, intRS);
#    else
    void FORT_HGSCON(Real*, intS, RealPS, intS, intS, RealRS);
    void FORT_HGCEN(Real*, intS, Real*, intS, intS);
    void FORT_HGINTS(Real*, intS, intS, Real*, intS, const Real*, intS, intS, intRS);
#    endif
#  endif
    void FORT_FANRST2(Real*, intS, intS, const Real*, intS, intRS, const int&);
    void FORT_FANINT2(Real*, intS, intS, const Real*, intS, intS, intRS);
#endif
}

void 
holy_grail_amr_multigrid::alloc(PArray<MultiFab>& Dest,
				PArray<MultiFab>& Source,
				PArray<MultiFab>& Coarse_source,
				PArray<MultiFab>& Sigma,
				Real H[], int Lev_min, int Lev_max)
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
#ifdef HG_CONSTANT
	assert(H[i] == H[0]);
#endif
	for (int mglev = mglev_max - 1; mglev >= 0; mglev--) 
	{
	    int rat = mg_domain[mglev+1].length(i) / mg_domain[mglev].length(i);
	    h[mglev][i] = rat * h[mglev+1][i];
	}
    }
    
#ifdef HG_USE_CACHE
    
    corr_scache.resize(mglev_max + 1, (copy_cache*) 0);
    
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	int mglev = ml_index[lev];
	dest_bcache.set(lev, new copy_cache(dest[lev], lev_interface[mglev],
	    mg_boundary, 1));
    }
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	corr_bcache.set(mglev, new copy_cache(corr[mglev], lev_interface[mglev],
	    mg_boundary, 1));
	corr_scache.set(mglev, new copy_cache(corr[mglev], lev_interface[mglev],
	    mg_boundary));
    }
    for (int mglev = 1; mglev <= mglev_max; mglev++) 
    {
	work_bcache.set(mglev, new copy_cache(work[mglev], lev_interface[mglev],
	    mg_boundary, 1));
    }
#endif
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
    
#ifdef HG_USE_CACHE
    cgw1_bcache = new copy_cache(cgwork[1], lev_interface[0], mg_boundary, 1);
#endif
    
    assert(cgwork[3].nGrow() == ib &&
	cgwork[4].nGrow() == ib &&
	cgwork[5].nGrow() == ib);
    
#ifdef HG_CONSTANT
    cgw_ucache.resize(8);
    for (i = 0; i < 8; i++) 
    {
	cgw_ucache.set(i, new unroll_cache(cgwork[i]));
    }
#endif
    
    for (int igrid = 0; igrid < mg_mesh[0].length(); igrid++) 
    {
	FArrayBox& gtmp = cgwork[7][igrid];
	const Box& valid = cgwork[7].box(igrid);
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
    if (mg_boundary.singular()) 
    {
	long sng;
	for (int i = 0; i < mg_mesh[0].length(); i++) 
	{
	    sng += mg_mesh[0][i].numPts();
	}
	singular = (sng == mg_domain[0].numPts());
    }
    
#ifdef HG_TERRAIN
    integrate = 1;
#endif
}

#ifndef HG_CONSTANT

void 
holy_grail_sigma_restrictor_class::fill(FArrayBox& patch,
					const Box& region,
					const FArrayBox& fgr,
					const IntVect& rat) const
{
    assert(patch.box().cellCentered());
    assert(rat[0] == 2 && rat[1] == 2 ||
	rat[0] == 2 && rat[1] == 1 ||
	rat[0] == 1 && rat[1] == 2);
    
#ifdef HG_TERRAIN
    FORT_HGSRST(patch.dataPtr(0), patch.dataPtr(1),
#  if (BL_SPACEDIM == 3)
	patch.dataPtr(2),
#  endif
	DIMLIST(patch.box()),
	DIMLIST(region),
	fgr.dataPtr(0), fgr.dataPtr(1),
#  if (BL_SPACEDIM == 3)
	fgr.dataPtr(2),
#  endif
	DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]));
    FORT_FACRST1(patch.dataPtr(BL_SPACEDIM),
	DIMLIST(patch.box()),
	DIMLIST(region),
	fgr.dataPtr(BL_SPACEDIM),
	DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]), 0);
    
#  if (BL_SPACEDIM == 2)
    patch.mult((Real) rat[1] / rat[0], region, 0, 1);
    patch.mult((Real) rat[0] / rat[1], region, 1, 1);
    // component 2 remains unchanged
#  else
    FORT_FACRST1(patch.dataPtr(BL_SPACEDIM+1),
	DIMLIST(patch.box()),
	DIMLIST(region),
	fgr.dataPtr(BL_SPACEDIM+1),
	DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]), 0);
    patch.mult((Real) rat[1] * rat[2] / rat[0], region, 0, 1);
    patch.mult((Real) rat[0] * rat[2] / rat[1], region, 1, 1);
    patch.mult((Real) rat[0] * rat[1] / rat[2], region, 2, 1);
    patch.mult((Real) rat[1],                   region, 3, 1);
    patch.mult((Real) rat[0],                   region, 4, 1);
#  endif
    
#else
    
    if (fgr.nComp() == 1) 
    {
	FORT_HGSRST(patch.dataPtr(0), patch.dataPtr(1),
#  if (BL_SPACEDIM == 3)
	    patch.dataPtr(2),
#  endif
	    DIMLIST(patch.box()),
	    DIMLIST(region),
	    fgr.dataPtr(), fgr.dataPtr(),
#  if (BL_SPACEDIM == 3)
	    fgr.dataPtr(),
#  endif
	    DIMLIST(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]));
    }
    else 
    {
	FORT_HGSRST(patch.dataPtr(0), patch.dataPtr(1),
#  if (BL_SPACEDIM == 3)
	    patch.dataPtr(2),
#  endif
	    DIMLIST(patch.box()),
	    DIMLIST(region),
	    fgr.dataPtr(0), fgr.dataPtr(1),
#  if (BL_SPACEDIM == 3)
	    fgr.dataPtr(2),
#  endif
	    DIMLIST(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]));
    }
    
#endif
}

#endif

void 
holy_grail_amr_multigrid::build_sigma(PArray<MultiFab>& Sigma)
{
   
#ifdef HG_TERRAIN
    
    // For terrain stencils we have as many sigma arrays passed as
    // arguments and used at the lev_interface as we build for internal
    // multigrid purposes.  This simplifies handling as we do not
    // need to maintain separate arrays for different purposes.
    
    int ncomp = 2 * BL_SPACEDIM - 1;
    
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
	    for (int igrid = 0; igrid < target.length(); igrid++) 
	    {
		target[igrid].copy(s_in[igrid], s_in.box(igrid), 0,
		    target.box(igrid), 0, ncomp);
	    }
	}
    }
    
    for (int mglev = mglev_max; mglev > 0; mglev--) 
    {
	IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
	restrict_level(sigma[mglev-1], sigma[mglev], rat, 
#ifdef HG_USE_CACHE
	    0,
#endif
	    holy_grail_sigma_restrictor_class());
    }
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	fill_borders(sigma[mglev], 0, lev_interface[mglev], boundary.terrain_sigma());
    }
    
#elif (! defined HG_CONSTANT)
    
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
		for (int igrid = 0; igrid < target.length(); igrid++) 
		{
		    target[igrid].copy(s_comp[igrid], s_comp.box(igrid), 0,
			target.box(igrid), i, 1);
		}
	    }
	}
    }
    
    sigma[mglev_max].setVal(1.e20);
    for (int igrid = 0; igrid < mg_mesh[mglev_max].length(); igrid++) 
    {
	sigma[mglev_max][igrid].copy(Sigma[lev_max][igrid], mg_mesh[mglev_max][igrid], 0,
	    mg_mesh[mglev_max][igrid], 0, 1);
    }
    
    if (mglev_max > 0) 
    {
	IntVect rat = mg_domain[mglev_max].length() / mg_domain[mglev_max-1].length();
	restrict_level(sigma_split[mglev_max-1], sigma[mglev_max], rat, 
#ifdef HG_USE_CACHE
	    0,
#endif
	    holy_grail_sigma_restrictor_class());
    }
    fill_borders(sigma[mglev_max], 
#ifdef HG_USE_CACHE
	0, 
#endif
	lev_interface[mglev_max], boundary.scalar());
    for (int mglev = mglev_max - 1; mglev > 0; mglev--) 
    {
	IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
	restrict_level(sigma_split[mglev-1], sigma_split[mglev], rat, 
#ifdef HG_USE_CACHE
	    0,
#endif
	    holy_grail_sigma_restrictor_class());
    }
    for (int mglev = 0; mglev < mglev_max; mglev++) 
    {
	fill_borders(sigma_split[mglev], 
#ifdef HG_USE_CACHE
	    0, 
#endif
	    lev_interface[mglev], boundary.scalar());
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
	    for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) 
	    {
		d[igrid].copy(s[igrid], i, 0);
	    }
	}
	delete sigma_split.remove(mglev);
	sigma.set(mglev, &sigma_nd[0][mglev]);
    }
    
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	sigma_nd[i].set(mglev_max, &sigma[mglev_max]);
    }
    
#  ifdef HG_SIGMA_NODE
    
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
	const Real hx = h[mglev][0];
	const Real hy = h[mglev][1];
#    if (BL_SPACEDIM == 3)
	const Real hz = h[mglev][2];
#    endif
	for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) 
	{
	    const Box& scbox = sigma[mglev][igrid].box();
	    const Box& snbox = sigma_node[mglev][igrid].box();
	    const Box& reg = lev_interface[mglev].part_fine(igrid);
	    FORT_HGSCON(sigma_node[mglev][igrid].dataPtr(),
		DIMLIST(snbox),
		sigma_nd[0][mglev][igrid].dataPtr(),
		sigma_nd[1][mglev][igrid].dataPtr(),
#    if (BL_SPACEDIM == 3)
		sigma_nd[2][mglev][igrid].dataPtr(),
#    endif
		DIMLIST(scbox),
		DIMLIST(reg),
#    if (BL_SPACEDIM == 2)
		hx, hy
#    else
		hx, hy, hz
#    endif
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
#  endif  /* HG_SIGMA_NODE */
    
#endif
    
    cen.resize(mglev_max + 1);
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	cen.set(mglev, new MultiFab(corr[mglev].boxArray(), 1,
	    dest[lev_min].nGrow()));
	MultiFab& ctmp = cen[mglev];
	ctmp.setVal(0.0);
	
#ifdef HG_TERRAIN
	
	for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) 
	{
	    const Box& cenbox = ctmp[igrid].box();
	    const Box& reg = lev_interface[mglev].part_fine(igrid);
	    const Box& sigbox = sigma[mglev][igrid].box();
	    FORT_HGCEN(ctmp[igrid].dataPtr(), DIMLIST(cenbox),
		sigma[mglev][igrid].dataPtr(),
		DIMLIST(sigbox),
		DIMLIST(reg));
	}
	
#elif (defined HG_CONSTANT)
	
	for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) 
	{
	    ctmp[igrid].setVal(1.0, lev_interface[mglev].part_fine(igrid), 0);
	}
	
#else
	
#  ifndef HG_SIGMA_NODE
	const Real hx = h[mglev][0];
	const Real hy = h[mglev][1];
#    if (BL_SPACEDIM == 3)
	const Real hz = h[mglev][2];
#    endif
#  endif
	for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) 
	{
	    const Box& cenbox = cen[mglev][igrid].box();
	    const Box& reg = lev_interface[mglev].part_fine(igrid);
#  ifndef HG_SIGMA_NODE
	    const Box& sigbox = sigma[mglev][igrid].box();
	    FORT_HGCEN(cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
		sigma_nd[0][mglev][igrid].dataPtr(),
		sigma_nd[1][mglev][igrid].dataPtr(),
#    if (BL_SPACEDIM == 3)
		sigma_nd[2][mglev][igrid].dataPtr(),
#    endif
		DIMLIST(sigbox),
		DIMLIST(reg),
#    if (BL_SPACEDIM == 2)
		hx, hy,
		IsRZ(), mg_domain[mglev].bigEnd(0) + 1
#    else
		hx, hy, hz
#    endif
		);
#  else
	    const Box& sigbox = sigma_node[mglev][igrid].box();
	    FORT_HGCEN(cen[mglev][igrid].dataPtr(), DIMLIST(cenbox),
		sigma_node[mglev][igrid].dataPtr(),
		DIMLIST(sigbox),
		DIMLIST(reg));
#  endif
	}
	
#endif
	
	clear_part_interface(ctmp, lev_interface[mglev]);
    }
    
#ifdef HG_CONSTANT
    mask.resize(mglev_max + 1);
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	mask.set(mglev, &cen[mglev]);
    }
#endif
    
#ifdef HG_SIGMA_NODE
    
    mask.resize(mglev_max + 1);
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	mask.set(mglev, new MultiFab(corr[mglev].boxArray(), 1,
	    dest[lev_min].nGrow()));
	MultiFab& mtmp = mask[mglev];
	mtmp.setVal(0.0);
	for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) 
	{
	    mtmp[igrid].setVal(1.0, lev_interface[mglev].part_fine(igrid), 0);
	}
	clear_part_interface(mtmp, lev_interface[mglev]);
    }
    
#endif
}

void holy_grail_amr_multigrid::clear()
{
    line_order.clear();
    line_after.clear();
    
    delete_sync_caches();
    
#ifdef HG_USE_CACHE
    delete cgw1_bcache;
#endif
    
#ifdef HG_CONSTANT
    for (int i = 0; i < 8; i++) 
    {
	delete cgw_ucache[i];
	cgw_ucache[i] = 0;
    }
#endif
    
    delete cgwork.remove(0);
    delete cgwork.remove(1);
    delete cgwork.remove(2);
    cgwork.remove(3);
    cgwork.remove(4);
    cgwork.remove(5);
    delete cgwork.remove(6);
    delete cgwork.remove(7);
    
#ifndef HG_CONSTANT
#  ifdef HG_TERRAIN
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	delete sigma.remove(mglev);
    }
#  elif (defined HG_SIGMA_NODE)
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	delete sigma.remove(mglev);
	delete sigma_node.remove(mglev);
    }
#  else

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
#  endif
#endif
    
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	delete cen.remove(mglev);
#ifdef HG_CONSTANT
	mask.remove(mglev);
#endif
#ifdef HG_SIGMA_NODE
	delete mask.remove(mglev);
#endif
    }
    
    delete [] h;
    if (source_owned) 
    {
	for (int lev = lev_min; lev <= lev_max; lev++) 
	{
	    if (source.defined(lev)) delete source.remove(lev);
	}
    }
    
#ifdef HG_USE_CACHE
    for (int lev = lev_min; lev <= lev_max; lev++) 
    {
	delete dest_bcache[lev];
	dest_bcache[lev] = 0;
    }
    for (int mglev = 0; mglev <= mglev_max; mglev++) 
    {
	delete corr_bcache[mglev];
	corr_bcache[mglev] = 0;
	delete corr_scache[mglev];
	corr_scache[mglev] = 0;
    }
    for (int mglev = 1; mglev <= mglev_max; mglev++) 
    {
	delete work_bcache[mglev];
	work_bcache[mglev] = 0;
    }
#endif
    
    amr_multigrid::clear();
}

int holy_grail_amr_multigrid::can_coarsen(const BoxArray& mesh,
					  const Box& domain) const
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
    return retval;
}

void holy_grail_amr_multigrid::sync_interfaces()
{
    for (int lev = lev_min+1; lev <= lev_max; lev++) 
    {
	int mglev = ml_index[lev];
	int mgc = ml_index[lev-1];
	IntVect rat = mg_domain[mglev].length() / mg_domain[mgc].length();
	MultiFab& target = dest[lev];
	for (int iface = 0; iface < lev_interface[mglev].nfaces(); iface++) 
	{
	    // find a fine grid touching this face
	    int igrid = lev_interface[mglev].fgrid(iface, 0);
	    if (igrid < 0)
		igrid = lev_interface[mglev].fgrid(iface, 1);
	    unsigned geo = lev_interface[mglev].fgeo(iface);
	    // reject fine-fine interfaces and those without an interior fine grid
	    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].fflag(iface) == 1)
		continue;
	    interpolate_patch(target[igrid], lev_interface[mglev].node_face(iface),
		dest[lev-1], rat,
		bilinear_interpolator_class(), lev_interface[mgc]);
	}
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
	for (int iface = 0; iface < lev_interface[mglev].nfaces(); iface++) 
	{
	    // find a fine grid touching this face
	    int igrid = lev_interface[mglev].fgrid(iface, 0);
	    if (igrid < 0)
		igrid = lev_interface[mglev].fgrid(iface, 1);
	    unsigned geo = lev_interface[mglev].fgeo(iface);
	    // use only exterior coarse-fine faces with an interior fine grid
	    const Box& nbox = lev_interface[mglev].node_face(iface);
	    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].fflag(iface) == 1 ||
		idomain.intersects(nbox))
		continue;
	    interpolate_patch(target[igrid], nbox, dest[lev-1], rat,
		bilinear_interpolator_class(), lev_interface[mgc]);
	}
    }
}

void holy_grail_amr_multigrid::mg_restrict_level(int lto, int lfrom)
{
    IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
    if (get_amr_level(lto) >= 0) 
    {
	if (integrate == 0) 
	{
	    restrict_level(resid[lto], 0, work[lfrom], rat, 
#ifdef HG_USE_CACHE
		work_bcache[lfrom],
#endif
		bilinear_restrictor_coarse_class(0),
		lev_interface[lfrom], mg_boundary);
	}
	else 
	{
	    restrict_level(resid[lto], 0, work[lfrom], rat, 
#ifdef HG_USE_CACHE
		work_bcache[lfrom],
#endif
		bilinear_restrictor_coarse_class(1),
		lev_interface[lfrom], mg_boundary);
	}
    }
    else 
    {
	mg_restrict(lto, lfrom);
    }
}

void holy_grail_amr_multigrid::mg_restrict(int lto, int lfrom)
{
    int igrid;
    fill_borders(work[lfrom], 
#ifdef HG_USE_CACHE
	work_bcache[lfrom], 
#endif
	lev_interface[lfrom], mg_boundary);
    IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
    for (igrid = 0; igrid < resid[lto].length(); igrid++) 
    {
	const Box& fbox = work[lfrom][igrid].box();
	const Box& cbox = resid[lto][igrid].box();
	const Box& creg = lev_interface[lto].part_fine(igrid);
	FORT_FANRST2(resid[lto][igrid].dataPtr(), DIMLIST(cbox), DIMLIST(creg),
	    work[lfrom][igrid].dataPtr(), DIMLIST(fbox),
	    D_DECL(rat[0], rat[1], rat[2]), integrate);
    }
    
    clear_part_interface(resid[lto], lev_interface[lto]);
}

#ifndef HG_CONSTANT

void holy_grail_interpolator_class::fill(FArrayBox& patch,
					 const Box& region,
					 const FArrayBox& cgr,
					 const Box& cb,
					 const IntVect& rat) const
{
#ifndef HG_SIGMA_NODE
    
    FORT_HGINTS(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	sigptr[0], sigptr[1],
#if (BL_SPACEDIM == 3)
	sigptr[2],
#endif
	DIMLIST(sigbox),
	cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb),
	D_DECL(rat[0], rat[1], rat[2]));
    
#else
    
    FORT_HGINTS(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	sigptr, DIMLIST(sigbox),
	cgr.dataPtr(), DIMLIST(cgr.box()), DIMLIST(cb),
	D_DECL(rat[0], rat[1], rat[2]));
    
#endif
}

#endif

void holy_grail_amr_multigrid::mg_interpolate_level(int lto, int lfrom)
{
    if (get_amr_level(lfrom) >= 0) 
    {
#ifdef HG_CONSTANT
	// multilevel interpolation, use general form
	amr_multigrid::mg_interpolate_level(lto, lfrom);
#else
	// general version---attempt to use special stencils for multilevel
	int ltmp = lfrom + 1;
	MultiFab& target = work[ltmp];
	IntVect rat = mg_domain[ltmp].length() / mg_domain[lfrom].length();
	for (int igrid = 0; igrid < target.length(); igrid++) 
	{
#ifdef HG_TERRAIN
	    Real *sigptr[BL_SPACEDIM];
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
		sigptr[i] = sigma[ltmp][igrid].dataPtr(i);
	    }
	    const Box& sigbox = sigma[ltmp][igrid].box();
#elif (! defined HG_SIGMA_NODE)
	    Real *sigptr[BL_SPACEDIM];
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
		sigptr[i] = sigma_nd[i][ltmp][igrid].dataPtr();
	    }
	    const Box& sigbox = sigma[ltmp][igrid].box();
#else
	    Real *sigptr = sigma_node[ltmp][igrid].dataPtr();
	    const Box& sigbox = sigma_node[ltmp][igrid].box();
#endif
	    interpolate_patch(target[igrid], target.box(igrid),
		corr[lfrom], rat,
		holy_grail_interpolator_class(sigptr, sigbox),
		lev_interface[lfrom], error_boundary);
	}
	if (lto > ltmp) 
	{
	    corr[ltmp].copy(target);
	    mg_interpolate_level(lto, ltmp);
	}
#endif
    }
    else 
    {
    /*
    // this general version not currently used, but may need to be revived
    // if special stencils are ever needed for the multilevel iteration
    MultiFab& target = work[lto];
    for (int igrid = 0; igrid < target.length(); igrid++) {
    Real *sigptr[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
    sigptr[i] = sigma_nd[i][lto][igrid].dataPtr();
    }
    const Box& sigbox = sigma[lto][igrid].box();
    corr[lfrom].interpolate_patch(target[igrid], target.box(igrid),
    corr[lfrom], rat,
    error_boundary, holy_grail_interpolator_class(sigptr, sigbox));
    }
	*/
	// multigrid interpolation, grids known to match up
	// special stencil needed for multigrid convergence
	IntVect rat = mg_domain[lto].length() / mg_domain[lfrom].length();
	for (int igrid = 0; igrid < mg_mesh[lto].length(); igrid++) 
	{
	    const Box& fbox = work[lto][igrid].box();
	    const Box& freg = work[lto].box(igrid);
	    const Box& cbox = corr[lfrom][igrid].box();
	    const Box& creg = corr[lfrom].box(igrid);
#ifdef HG_CONSTANT
	    FANINT2(work[lto][igrid].dataPtr(), DIMLIST(fbox), DIMLIST(freg),
		corr[lfrom][igrid].dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		D_DECL(rat[0], rat[1], rat[2]));
#else
#ifndef HG_SIGMA_NODE
	    const Box& sigbox = sigma[lto][igrid].box();
#else
	    const Box& sigbox = sigma_node[lto][igrid].box();
#endif
	    FORT_HGINTS(work[lto][igrid].dataPtr(), DIMLIST(fbox), DIMLIST(freg),
#ifdef HG_TERRAIN
		sigma[lto][igrid].dataPtr(0),
		sigma[lto][igrid].dataPtr(1),
#  if (BL_SPACEDIM == 3)
		sigma[lto][igrid].dataPtr(2),
#  endif
#elif (! defined HG_SIGMA_NODE)
		sigma_nd[0][lto][igrid].dataPtr(),
		sigma_nd[1][lto][igrid].dataPtr(),
#  if (BL_SPACEDIM == 3)
		sigma_nd[2][lto][igrid].dataPtr(),
#  endif
#else
		sigma_node[lto][igrid].dataPtr(),
#endif
		DIMLIST(sigbox),
		corr[lfrom][igrid].dataPtr(), DIMLIST(cbox), DIMLIST(creg),
		D_DECL(rat[0], rat[1], rat[2]));
#endif
	}
    }
}
