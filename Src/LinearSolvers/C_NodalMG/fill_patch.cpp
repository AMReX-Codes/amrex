
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_FIPRODC   iprodc_
#  define FORT_FIPRODN   iprodn_
#  define FORT_FFCPY2    fcpy2_
#  define FORT_FFCPY     fcpy_
#else
#  define FORT_FIPRODC   IPRODC
#  define FORT_FIPRODN   IPRODN
#  define FORT_FFCPY2    FCPY2
#  define FORT_FFCPY     FCPY
#endif

extern "C"
{
    void FORT_FIPRODC(const Real*, intS, const Real*, intS, intS, Real*);
    void FORT_FIPRODN(const Real*, intS, const Real*, intS, intS, Real*);
    void FORT_FFCPY(Real*, intS, intS, const Real*, intS, const int&);
#if (BL_SPACEDIM == 2)
    void FORT_FFCPY2(Real*, intS, const Real*, intS, intS, const int*, const int&);
#else
    void FORT_FFCPY2(Real*, intS, const Real*, intS, intS, const int*, const int*, const int&);
#endif
}

void internal_copy(MultiFab& r, int destgrid, int srcgrid, const Box& b) 
{
    Real* dptr = r[destgrid].dataPtr();
    const Real* sptr = r[srcgrid].dataPtr();
    const Box& dbx = r[destgrid].box();
    const Box& sbx = r[srcgrid].box();
    FORT_FFCPY(dptr, DIMLIST(dbx), DIMLIST(b), sptr, DIMLIST(sbx), r.nComp());
}


Real inner_product(const MultiFab& r, const MultiFab& s)
{
    assert(r.ok() && s.ok());
    assert(r.nComp() == 1);
    assert(s.nComp() == 1);
    assert(type(r) == type(s));
    
    Real sum = 0.0;
    
    if (type(r) == IntVect::TheCellVector()) 
    {
	for (ConstMultiFabIterator rcmfi(r); rcmfi.isValid(); ++rcmfi) 
	{
	    ConstDependentMultiFabIterator scmfi(rcmfi, s);
	    const Box& rbox = rcmfi->box();
	    const Box& sbox = scmfi->box();
	    const Box& reg  = rcmfi.validbox();
	    FORT_FIPRODC(rcmfi->dataPtr(), DIMLIST(rbox), scmfi->dataPtr(), DIMLIST(sbox), DIMLIST(reg), &sum);
	}
    }
    else if (type(r) == IntVect::TheNodeVector()) 
    {
	for (ConstMultiFabIterator rcmfi(r); rcmfi.isValid(); ++rcmfi) 
	{
	    ConstDependentMultiFabIterator scmfi(rcmfi, s);
	    const Box& rbox = rcmfi->box();
	    const Box& sbox = scmfi->box();
	    const Box& reg  = rcmfi.validbox();
	    FORT_FIPRODN(rcmfi->dataPtr(), DIMLIST(rbox), scmfi->dataPtr(), DIMLIST(sbox), DIMLIST(reg), &sum);
	}
    }
    else 
    {
	BoxLib::Error("inner_product---only supported for CELL- or NODE-based data");
    }
    ParallelDescriptor::ReduceRealSum(sum);
    return sum;
}


int find_patch(const Box& region, const MultiFab& r)
{
    if (r.nGrow() == 0 ) 
    {
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    assert( r[igrid].box() == r.box(igrid) );
	    if (r[igrid].box().contains(region))
		return igrid;
	}
    }
    else 
    {
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    if (r.box(igrid).contains(region))
		return igrid;
	}
    }
    return -1;
}

static bool fill_patch_blindly(FArrayBox& patch, const Box& region, const MultiFab& r)
{
    if (r.nGrow() == 0) 
    {
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    if (r[igrid].box().contains(region)) 
	    {
		patch.copy(r[igrid], region, 0, region, 0, patch.nComp());
		return true;
	    }
	}
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    if (r[igrid].box().intersects(region)) 
	    {
		Box tb = region & r[igrid].box();
		patch.copy(r[igrid], tb, 0, tb, 0, patch.nComp());
	    }
	}
    }
    else 
    {
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    Box tb = grow(r[igrid].box(), -r.nGrow());
	    if (tb.contains(region)) 
	    {
		patch.copy(r[igrid], region, 0, region, 0, patch.nComp());
		return true;
	    }
	}
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    Box tb = grow(r[igrid].box(), -r.nGrow());
	    if (tb.intersects(region)) 
	    {
		tb &= region;
		patch.copy(r[igrid], tb, 0, tb, 0, patch.nComp());
	    }
	}
    }
    return false;
}

static bool fill_exterior_patch_blindly(FArrayBox& patch,
				const Box& region,
				const MultiFab& r,
				const level_interface& lev_interface,
				const amr_boundary_class* bdy)
{
    const BoxArray& em = lev_interface.exterior_mesh();
    for (int igrid = 0; igrid < em.length(); igrid++) 
    {
	int jgrid = lev_interface.direct_exterior_ref(igrid);
	if (jgrid >= 0) 
	{
	    Box tb;
	    tb = em[igrid];
	    tb.convert(type(r));
	    if (tb.contains(region)) 
	    {
		assert(bdy != 0);
		bdy->fill(patch, region, r[jgrid], jgrid, lev_interface.domain());
		return true;
	    }
	    if (tb.intersects(region)) 
	    {
		tb &= region;
		assert(bdy != 0);
		bdy->fill(patch, tb, r[jgrid], jgrid, lev_interface.domain());
	    }
	}
    }
    return false;
}

bool fill_patch(FArrayBox& patch, const Box& region,
	       const MultiFab& r,
	       const level_interface& lev_interface,
	       const amr_boundary_class* bdy,
	       int idim, int index)
{
    if (!region.ok())
	return true;
    
    assert(patch.nComp() == r.nComp());
    assert(type(patch) == type(r));
    assert(lev_interface.ok());
    assert( idim >= -1 && idim < BL_SPACEDIM );
    
    Box tdomain = lev_interface.domain();
    tdomain.convert(type(patch));
    Box idomain = grow(tdomain, IntVect::TheZeroVector() - type(r));
    
    if (idim == -1 ) 
    {
	if (idomain.contains(region) || bdy == 0) 
	{
	    return fill_patch_blindly(patch, region, r);
	}
	else if (!tdomain.intersects(region)) 
	{
	    return fill_exterior_patch_blindly(patch, region, r, lev_interface, bdy);
	}
	else if (idomain.intersects(region)) 
	{
	    if (fill_patch_blindly(patch, region, r) )
		return true;
	    else
		return fill_exterior_patch_blindly(patch, region, r, lev_interface, bdy);
	}
	else 
	{
	    if (fill_exterior_patch_blindly(patch, region, r, lev_interface, bdy) )
		return true;
	    else
		return fill_patch_blindly(patch, region, r);
	}
    }
    else
    {
	Array<int> gridnum(lev_interface.ngrids(idim)+1);
	gridnum[0] = -1;
	for (int i = 0; i < lev_interface.ngrids(idim); i++) 
	{
	    int igrid = lev_interface.grid(idim, index,i);
	    if (igrid != -1) 
	    {
		for (int j = 0; gridnum[j] != igrid; j++) 
		{
		    if (gridnum[j] == -1) 
		    {
			gridnum[j] = igrid;
			gridnum[j+1] = -1;
			if (igrid >= 0) 
			{
			    Box tb = r.box(igrid);
			    tb &= region;
			    const Box& rbox = r[igrid].box();
			    FORT_FFCPY(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(tb), r[igrid].dataPtr(), DIMLIST(rbox), patch.nComp());
			}
			else 
			{
			    igrid = -2 - igrid;
			    Box tb = lev_interface.exterior_mesh()[igrid];
			    tb.convert(type(r));
			    tb &= region;
			    bdy->fill(patch, tb, r[lev_interface.direct_exterior_ref(igrid)], lev_interface.direct_exterior_ref(igrid), lev_interface.domain());
			}
			break;
		    }
		}
	    }
	}
    }
    return true;
}

static void sync_internal_borders(MultiFab& r, const level_interface& lev_interface)
{
    assert(type(r) == IntVect::TheNodeVector());

    // PARALLEL
    for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
    {
	int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	int jgrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
	// only do interior faces with fine grid on both sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
	    break;
	internal_copy(r, jgrid, igrid, lev_interface.node_box(level_interface::FACEDIM, iface));
    }
#if (BL_SPACEDIM == 2)
    // PARALLEL
    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
    {
	int igrid = lev_interface.grid(0, icor, 0);
	int jgrid = lev_interface.grid(0, icor, 3);
	// only do interior corners with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.grid(0, icor, 1))
	    internal_copy(r, jgrid, igrid, lev_interface.box(0, icor));
    }
#else
    // PARALLEL
    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
    {
	int igrid = lev_interface.grid(1, iedge, 0);
	int jgrid = lev_interface.grid(1, iedge, 3);
	// only do interior edges with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(1, iedge) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.grid(1, iedge, 1))
	    internal_copy(r, jgrid, igrid, lev_interface.node_box(1, iedge));
    }
    // PARALLEL
    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
    {
	int igrid = lev_interface.grid(0, icor, 0);
	int jgrid = lev_interface.grid(0, icor, 7);
	// only do interior corners with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (lev_interface.grid(0, icor, 3) == lev_interface.grid(0, icor, 1)) 
	{
	    if (jgrid != lev_interface.grid(0, icor, 3)) 
	    {
		internal_copy(r, jgrid, igrid, lev_interface.box(0, icor));
		jgrid = lev_interface.grid(0, icor, 5);
		if (jgrid != lev_interface.grid(0, icor, 7))
		    internal_copy(r, jgrid, igrid, lev_interface.box(0, icor));
	    }
	}
	else if (lev_interface.grid(0, icor, 5) == lev_interface.grid(0, icor, 1)) 
	{
	    if (jgrid != lev_interface.grid(0, icor, 5)) 
	    {
		internal_copy(r, jgrid, igrid, lev_interface.box(0, icor));
		jgrid = lev_interface.grid(0, icor, 3);
		if (jgrid != lev_interface.grid(0, icor, 7)) 
		{
		    internal_copy(r, jgrid, igrid, lev_interface.box(0, icor));
		    if (jgrid == lev_interface.grid(0, icor, 2)) 
		    {
			jgrid = lev_interface.grid(0, icor, 6);
			if (jgrid != lev_interface.grid(0, icor, 7))
			    internal_copy(r, jgrid, igrid, lev_interface.box(0, icor));
		    }
		}
	    }
	}
    }
#endif
}

// local function used only by fill_internal_borders:

static inline void node_dirs(int dir[2], const IntVect& typ)
{
    if (typ[0] == IndexType::NODE) 
    {
	dir[0] = 0;
	if (typ[1] == IndexType::NODE)
	    dir[1] = 1;
	else
	    dir[1] = 2;
    }
    else 
    {
	dir[0] = 1;
	dir[1] = 2;
    }
}

// The sequencing used in fill_internal_borders, fcpy2 and set_border_cache
// (narrow x, medium y, wide z) is necessary to avoid overwrite problems
// like those seen in the sync routines.  Boundary copies are all wide
// regardless of direction and come after interior copies---overwrite
// difficulties are avoided since grids can't bridge a boundary.

// Modifications are necessary in 3D to deal with lack of diagonal
// communication across edges at the coarse-fine lev_interface.  These
// modifications take the form of narrowing certain copies to avoid
// overwriting good values with bad ones.

static void fill_internal_borders(MultiFab& r, const level_interface& lev_interface, int w, bool hg_terrain)
{
    assert(type(r) == IntVect::TheCellVector() || type(r) == IntVect::TheNodeVector() );
    w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;
    if ( type(r) == IntVect::TheNodeVector() ) 
    {
#if (BL_SPACEDIM == 3)
	if(hg_terrain)
	{
	    // attempt to deal with corner-coupling problem with 27-point stencils
	    //PARALLEL
	    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
	    {
		if (lev_interface.geo(1, iedge) == level_interface::ALL)
		    continue;
		int igrid = lev_interface.grid(1, iedge, 0);
		int jgrid = lev_interface.grid(1, iedge, 3);
		if (igrid >= 0 && jgrid >= 0) 
		{
		    int kgrid = lev_interface.grid(1, iedge, 1);
		    if (kgrid == -1)
			kgrid = lev_interface.grid(1, iedge, 2);
		    if (kgrid != -1 && kgrid != igrid && kgrid != jgrid)
		    {	
			Box b;
			int dir[2];
			node_dirs(dir, lev_interface.box(1, iedge).type());
			if (kgrid == lev_interface.grid(1, iedge, 1)) 
			{
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[0], -1));
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, igrid, jgrid, b.shift(dir[1],  1));
			}
			else 
			{
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[1], -1));
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, igrid, jgrid, b.shift(dir[0],  1));
			}
		    }
		}
		igrid = lev_interface.grid(1, iedge, 1);
		jgrid = lev_interface.grid(1, iedge, 2);
		if (igrid >= 0 && jgrid >= 0) 
		{
		    int kgrid = lev_interface.grid(1, iedge, 0);
		    if (kgrid == -1)
			kgrid = lev_interface.grid(1, iedge, 3);
		    if (kgrid != -1 && kgrid != igrid && kgrid != jgrid) 
		    {
			Box b;		    
			int dir[2];
			node_dirs(dir, lev_interface.box(1, iedge).type());
			if (kgrid == lev_interface.grid(1, iedge, 0)) 
			{
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[0],  1));
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, igrid, jgrid, b.shift(dir[1],  1));
			}
			else 
			{
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[1], -1));
			    b = lev_interface.node_box(1, iedge);
			    internal_copy(r, igrid, jgrid, b.shift(dir[0], -1));
			}
		    }
		}
	    }
	}
#endif
	// PARALLEL
	for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	{
	    const int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	    const int jgrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
	    if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
		break;
	    const Box& b = lev_interface.node_box(level_interface::FACEDIM, iface);
	    Real* ptra = r[igrid].dataPtr();
	    const Real* ptrb = r[jgrid].dataPtr();
	    const Box& boxa = r[igrid].box();
	    const Box& boxb = r[jgrid].box();
#  if (BL_SPACEDIM == 2)
	    FORT_FFCPY2(ptra, DIMLIST(boxa), ptrb, DIMLIST(boxb), DIMLIST(b), &w, r.nComp());
#  else
	    const int ibord = r.nGrow();
	    FORT_FFCPY2(ptra, DIMLIST(boxa), ptrb, DIMLIST(boxb), DIMLIST(b), &w, &ibord, r.nComp());
#  endif
	}
    }
    else if (type(r) == IntVect::TheCellVector()) 
    {
	// PARALLEL
	for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	{
	    const int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	    const int jgrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
	    if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
		break;
	    const int idim = lev_interface.fdim(iface);
#if (BL_SPACEDIM == 2)
	    Box b = lev_interface.box(level_interface::FACEDIM, iface);
	    if (idim == 1)
		b.grow(0, w);
	    b.growLo(idim, w).convert(IntVect::TheCellVector());
	    internal_copy(r, jgrid, igrid, b);
	    internal_copy(r, igrid, jgrid, b.shift(idim, w));
#else
	    Box bj = lev_interface.box(level_interface::FACEDIM, iface);
	    Box bi = lev_interface.box(level_interface::FACEDIM, iface);
	    for (int i = 0; i < idim; i++) 
	    {
		if (r.box(jgrid).smallEnd(i) == bj.smallEnd(i))
		    bj.growLo(i, w);
		if (r.box(jgrid).bigEnd(i) == bj.bigEnd(i))
		    bj.growHi(i, w);
		if (r.box(igrid).smallEnd(i) == bi.smallEnd(i))
		    bi.growLo(i, w);
		if (r.box(igrid).bigEnd(i) == bi.bigEnd(i))
		    bi.growHi(i, w);
	    }
	    bj.growLo(idim, w).convert(IntVect::TheCellVector());
	    bi.growHi(idim, w).convert(IntVect::TheCellVector());
	    internal_copy(r, jgrid, igrid, bj);
	    internal_copy(r, igrid, jgrid, bi);
#endif
	}
    }
}

void clear_part_interface(MultiFab& r, const level_interface& lev_interface)
{
    if (r.nComp() != 1)
	BoxLib::Error("clear_part_interface---only single components currently supported");
    
    assert(type(r) == IntVect::TheNodeVector());
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	// PARALLEL
	for (int ibox = 0; ibox < lev_interface.nboxes(i); ibox++) 
	{
	    // coarse-fine face contained in part_fine grid, or orphan edge/corner
	    int igrid = lev_interface.aux(i, ibox);
	    if ( igrid >= 0 )
		r[igrid].setVal(0.0, lev_interface.node_box(i, ibox), 0);
	}
    }
}

void interpolate_patch(FArrayBox& patch, const Box& region,
		       const MultiFab& r, const IntVect& rat,
		       const amr_interpolator_class& interp,
		       const level_interface& lev_interface)
{
    assert(region.sameType(patch.box()));
    const Box cb = interp.box(region, rat);
    const int igrid = find_patch(cb, r);
    if (igrid == -1) 
    {
	FArrayBox cgr(cb, r.nComp());
	fill_patch(cgr, cb, r, lev_interface, 0);
	interp.fill(patch, region, cgr, cb, rat);
    }
    else 
    {
	interp.fill(patch, region, r[igrid], cb, rat);
    }
}

void restrict_level(MultiFab& dest, 
		    MultiFab& r, const IntVect& rat,
		    const amr_restrictor_class& restric,
		    const level_interface& lev_interface,
		    const amr_boundary_class* bdy)
{
    assert(type(dest) == type(r));
    // PARALLEL
    for (int jgrid = 0; jgrid < dest.length(); jgrid++) 
    {
        // restrict_patch(dest[jgrid], dest.box(jgrid), r, rat, restric, lev_interface, bdy);
	const Box& region = dest.box(jgrid);
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    Box cbox = r.box(igrid);
	    cbox = restric.box(cbox, rat);
	    if (region.intersects(cbox)) 
	    {
		cbox &= region;
		restric.fill(dest[jgrid], cbox, r[igrid], rat);
	    }
	}
    }
    if ( lev_interface.ok() )
    {
	restric.fill_interface( dest, r, lev_interface, bdy, rat);
    }
}

void sync_borders(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy)
{
    sync_internal_borders(r, lev_interface);
    assert(bdy != 0);
    bdy->sync_borders(r, lev_interface);
}

void fill_borders(MultiFab& r,
		  const level_interface& lev_interface,
		  const amr_boundary_class* bdy, int w, bool hg_terrain)
{
    fill_internal_borders(r, lev_interface, w, hg_terrain);
    assert(bdy != 0);
    bdy->fill_borders(r, lev_interface, w);
}
