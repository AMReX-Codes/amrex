#include "cache.H"

copy_cache::copy_cache(int Nsets, Real *Dptr, Real *Sptr)
   : nsets(Nsets), dptr(Dptr), sptr(Sptr)
{
#if (BL_SPACEDIM == 2)
    dstart = new int[5 * nsets];
    sstart = dstart + nsets;
    dstrid = dstart + 2 * nsets;
    sstrid = dstart + 3 * nsets;
    nvals  = dstart + 4 * nsets;
    for (int i = 0; i < nsets; i++)
	nvals[i] = 0;
#else
    dstart = new int[8 * nsets];
    sstart = dstart + nsets;
    dstrid1 = dstart + 2 * nsets;
    dstrid2 = dstart + 3 * nsets;
    sstrid1 = dstart + 4 * nsets;
    sstrid2 = dstart + 5 * nsets;
    nvals1  = dstart + 6 * nsets;
    nvals2  = dstart + 7 * nsets;
    for (int i = 0; i < nsets; i++) 
    {
	nvals1[i] = 0;
	nvals2[i] = 0;
    }
#endif
}

// sync cache

copy_cache::copy_cache(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy)
{
    assert(r.length() > 0);
    assert(r.nComp() == 1);
    assert(type(r) == IntVect::TheNodeVector());
    
    nsets = 0;
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	for (int igrid = 0; igrid < lev_interface.nboxes(i); igrid++) 
	{
	    if (lev_interface.geo(i, igrid) != level_interface::ALL)
		break;
	    nsets++;
	}
    }
    
    // While it is possible for two copies to occur at a corner in 3D,
    // this implies that the edge from that corner in the +x direction
    // exists and does not copy, so the above counts are sufficient.
    
    Real *baseptr = r[0].dataPtr();
    dptr = sptr = baseptr;
    
#if (BL_SPACEDIM == 2)
    dstart = new int[5 * nsets];
    sstart = dstart + nsets;
    dstrid = dstart + 2 * nsets;
    sstrid = dstart + 3 * nsets;
    nvals  = dstart + 4 * nsets;
    for (int i = 0; i < nsets; i++)
	nvals[i] = 0;
#else
    dstart = new int[8 * nsets];
    sstart = dstart + nsets;
    dstrid1 = dstart + 2 * nsets;
    dstrid2 = dstart + 3 * nsets;
    sstrid1 = dstart + 4 * nsets;
    sstrid2 = dstart + 5 * nsets;
    nvals1  = dstart + 6 * nsets;
    nvals2  = dstart + 7 * nsets;
    for (int i = 0; i < nsets; i++) 
    {
	nvals1[i] = 0;
	nvals2[i] = 0;
    }
#endif
    
    int iset = 0;
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	const int igrid = lev_interface.fgrid(iface, 0);
	const int jgrid = lev_interface.fgrid(iface, 1);
	if (igrid < 0 || jgrid < 0 || lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	const Box& b = lev_interface.node_face(iface);
#if (BL_SPACEDIM == 2)
	int dstartj, sstarti, stridi, stridj, nvals;
	stridi = r[igrid].box().length(0);
	stridj = r[jgrid].box().length(0);
	sstarti = r[igrid].dataPtr() - baseptr +
	    b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	    stridi * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
	dstartj = r[jgrid].dataPtr() - baseptr +
	    b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	    stridj * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
	if (lev_interface.fdim(iface) == 0) 
	{
	    nvals = b.length(1);
	}
	else 
	{
	    nvals = b.length(0);
	    stridi = 1;
	    stridj = 1;
	}
	set(iset++, dstartj, sstarti, stridj, stridi, nvals);
#else
	int stridi1 = r[igrid].box().length(0);
	int stridi2 = stridi1 * r[igrid].box().length(1);
	int stridj1 = r[jgrid].box().length(0);
	int stridj2 = stridj1 * r[jgrid].box().length(1);
	int sstarti = r[igrid].dataPtr() - baseptr +
	    b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	    stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
	    stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	int dstartj = r[jgrid].dataPtr() - baseptr +
	    b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	    stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
	    stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	int nvals1, nvals2;
	if (lev_interface.fdim(iface) == 0) 
	{
	    nvals1 = b.length(1);
	    nvals2 = b.length(2);
	}
	else if (lev_interface.fdim(iface) == 1) 
	{
	    nvals1 = b.length(0);
	    nvals2 = b.length(2);
	    stridi1 = 1;
	    stridj1 = 1;
	}
	else 
	{
	    nvals1 = b.length(0);
	    nvals2 = b.length(1);
	    stridi2 = stridi1;
	    stridj2 = stridj1;
	    stridi1 = 1;
	    stridj1 = 1;
	}
	set(iset++, dstartj, sstarti, stridj1, stridj2, stridi1, stridi2, nvals1, nvals2);
#endif
    }
    
#if (BL_SPACEDIM == 2)
    for (int icor = 0; icor < lev_interface.ncorners(); icor++) 
    {
	const int igrid = lev_interface.cgrid(icor, 0);
	const int jgrid = lev_interface.cgrid(icor, 3);
	// only do interior corners with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.cgrid(icor, 1)) 
	{
	    const Box& b = lev_interface.corner(icor);
	    int dstartj, sstarti, stridi, stridj;
	    stridi = r[igrid].box().length(0);
	    stridj = r[jgrid].box().length(0);
	    sstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - r[igrid].box().smallEnd(0) +
		stridi * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
		stridj * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
	    set(iset++, dstartj, sstarti, 0, 0, 1);
	}
    }
#else
    for (int iedge = 0; iedge < lev_interface.nedges(); iedge++) 
    {
	const int igrid = lev_interface.egrid(iedge, 0);
	const int jgrid = lev_interface.egrid(iedge, 3);
	// only do interior edges with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(1, iedge) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.egrid(iedge, 1)) 
	{
	    const Box& b = lev_interface.node_edge(iedge);
	    int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2, nvals1;
	    stridi1 = r[igrid].box().length(0);
	    stridi2 = stridi1 * r[igrid].box().length(1);
	    stridj1 = r[jgrid].box().length(0);
	    stridj2 = stridj1 * r[jgrid].box().length(1);
	    sstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - r[igrid].box().smallEnd(0) +
		stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
		stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
		stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	    if ((nvals1 = b.length(0)) > 1) 
	    {
		stridi1 = 1;
		stridj1 = 1;
	    }
	    else if ((nvals1 = b.length(2)) > 1) 
	    {
		stridi1 = stridi2;
		stridj1 = stridj2;
	    }
	    else 
	    {
		nvals1 = b.length(1);
	    }
	    set(iset++, dstartj, sstarti, stridj1, 0, stridi1, 0, nvals1, 1);
	}
    }
    for (int icor = 0; icor < lev_interface.ncorners(); icor++) 
    {
	const int igrid = lev_interface.cgrid(icor, 0);
	int jgrid = lev_interface.cgrid(icor, 7);
	// only do interior corners with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (lev_interface.cgrid(icor, 3) == lev_interface.cgrid(icor, 1)) 
	{
	    if (jgrid != lev_interface.cgrid(icor, 3)) 
	    {
		const Box& b = lev_interface.corner(icor);
		int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2;
		stridi1 = r[igrid].box().length(0);
		stridi2 = stridi1 * r[igrid].box().length(1);
		stridj1 = r[jgrid].box().length(0);
		stridj2 = stridj1 * r[jgrid].box().length(1);
		sstarti = r[igrid].dataPtr() - baseptr +
		    b.smallEnd(0) - r[igrid].box().smallEnd(0) +
		    stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
		    stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		dstartj = r[jgrid].dataPtr() - baseptr +
		    b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
		    stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		    stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
		jgrid = lev_interface.cgrid(icor, 5);
		if (jgrid != lev_interface.cgrid(icor, 7)) 
		{
		    stridj1 = r[jgrid].box().length(0);
		    stridj2 = stridj1 * r[jgrid].box().length(1);
		    dstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
		}
	    }
	}
	else if (lev_interface.cgrid(icor, 5) == lev_interface.cgrid(icor, 1)) 
	{
	    if (jgrid != lev_interface.cgrid(icor, 5)) 
	    {
		const Box& b = lev_interface.corner(icor);
		int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2;
		stridi1 = r[igrid].box().length(0);
		stridi2 = stridi1 * r[igrid].box().length(1);
		stridj1 = r[jgrid].box().length(0);
		stridj2 = stridj1 * r[jgrid].box().length(1);
		sstarti = r[igrid].dataPtr() - baseptr +
		    b.smallEnd(0) - r[igrid].box().smallEnd(0) +
		    stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
		    stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		dstartj = r[jgrid].dataPtr() - baseptr +
		    b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
		    stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		    stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
		jgrid = lev_interface.cgrid(icor, 3);
		if (jgrid != lev_interface.cgrid(icor, 7)) 
		{
		    stridj1 = r[jgrid].box().length(0);
		    stridj2 = stridj1 * r[jgrid].box().length(1);
		    dstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
		    if (jgrid == lev_interface.cgrid(icor, 2)) 
		    {
			jgrid = lev_interface.cgrid(icor, 6);
			if (jgrid != lev_interface.cgrid(icor, 7)) 
			{
			    stridj1 = r[jgrid].box().length(0);
			    stridj2 = stridj1 * r[jgrid].box().length(1);
			    dstartj = r[jgrid].dataPtr() - baseptr +
				b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
				stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
				stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
			    set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
			}
		    }
		}
	    }
	}
    }
#endif
    assert(bdy != 0);
    bdy->set_sync_cache(this, nsets, iset, r, lev_interface);
    nsets = iset;
}

// border cache

// local function used by the copy_cache border constructor
// (identical to one used in fill_patch.C)
#ifdef HG_TERRAIN
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
#endif

copy_cache::copy_cache(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy, int w)
{
    assert(r.length() > 0);
    assert(r.nComp() == 1);
    assert(type(r) == IntVect::TheNodeVector());
    
    w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;
    
    assert(w == 1);
    
    nsets = 0;
#if ((BL_SPACEDIM == 3) && (defined HG_TERRAIN))
    // attempt to deal with corner-coupling problem with 27-point stencils
    int iedge, kgrid;
    for (iedge = 0; iedge < lev_interface.nedges(); iedge++) 
    {
	if (lev_interface.geo(1, iedge) == level_interface::ALL)
	    continue;
	int igrid = lev_interface.egrid(iedge, 0);
	int jgrid = lev_interface.egrid(iedge, 3);
	if (igrid >= 0 && jgrid >= 0) 
	{
	    kgrid = lev_interface.egrid(iedge, 1);
	    if (kgrid == -1)
		kgrid = lev_interface.egrid(iedge, 2);
	    if (kgrid != -1 && kgrid != igrid && kgrid != jgrid)
		nsets += 2;
	}
	igrid = lev_interface.egrid(iedge, 1);
	jgrid = lev_interface.egrid(iedge, 2);
	if (igrid >= 0 && jgrid >= 0) \
	{
	    kgrid = lev_interface.egrid(iedge, 0);
	    if (kgrid == -1)
		kgrid = lev_interface.egrid(iedge, 3);
	    if (kgrid != -1 && kgrid != igrid && kgrid != jgrid)
		nsets += 2;
	}
    }
#endif
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	if (lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	if (lev_interface.fgrid(iface, 0) >= 0)
	    nsets++;
	if (lev_interface.fgrid(iface, 1) >= 0)
	    nsets++;
    }
    
    Real *baseptr = r[0].dataPtr();
    dptr = sptr = baseptr;
    
#if (BL_SPACEDIM == 2)
    dstart = new int[5 * nsets];
    sstart = dstart + nsets;
    dstrid = dstart + 2 * nsets;
    sstrid = dstart + 3 * nsets;
    nvals  = dstart + 4 * nsets;
    for (int i = 0; i < nsets; i++)
	nvals[i] = 0;
#else
    dstart = new int[8 * nsets];
    sstart = dstart + nsets;
    dstrid1 = dstart + 2 * nsets;
    dstrid2 = dstart + 3 * nsets;
    sstrid1 = dstart + 4 * nsets;
    sstrid2 = dstart + 5 * nsets;
    nvals1  = dstart + 6 * nsets;
    nvals2  = dstart + 7 * nsets;
    for (int i = 0; i < nsets; i++) 
    {
	nvals1[i] = 0;
	nvals2[i] = 0;
    }
#endif
    
    int iset = 0;
#if ((BL_SPACEDIM == 3) && (defined HG_TERRAIN))
    // attempt to deal with corner-coupling problem with 27-point stencils
    int dir[2];
    Box b;
    int sstartj, dstartj, sstarti, dstarti;
    int stridi1, stridi2, stridj1, stridj2, nvals1;
    for (int iedge = 0; iedge < lev_interface.nedges(); iedge++) 
    {
	if (lev_interface.geo(1, iedge) == level_interface::ALL)
	    continue;
	int igrid = lev_interface.egrid(iedge, 0);
	int jgrid = lev_interface.egrid(iedge, 3);
	if (igrid >= 0 && jgrid >= 0) 
	{
	    int kgrid = lev_interface.egrid(iedge, 1);
	    if (kgrid == -1)
		kgrid = lev_interface.egrid(iedge, 2);
	    if (kgrid != -1 && kgrid != igrid && kgrid != jgrid) 
	    {
		node_dirs(dir, lev_interface.edge(iedge).type());
		stridi1 = r[igrid].box().length(0);
		stridi2 = stridi1 * r[igrid].box().length(1);
		stridj1 = r[jgrid].box().length(0);
		stridj2 = stridj1 * r[jgrid].box().length(1);
		if (kgrid == lev_interface.egrid(iedge, 1)) 
		{
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[0], -1);
		    sstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		    dstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[1],  1);
		    sstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    dstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		}
		else 
		{
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[1], -1);
		    sstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		    dstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[0],  1);
		    sstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    dstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		}
		if ((nvals1 = b.length(0)) > 1) 
		{
		    stridi1 = 1;
		    stridj1 = 1;
		}
		else if ((nvals1 = b.length(2)) > 1) 
		{
		    stridi1 = stridi2;
		    stridj1 = stridj2;
		}
		else 
		{
		    nvals1 = b.length(1);
		}
		set(iset++, dstartj, sstarti, stridj1, 0, stridi1, 0, nvals1, 1);
		set(iset++, dstarti, sstartj, stridi1, 0, stridj1, 0, nvals1, 1);
	    }
	}
	igrid = lev_interface.egrid(iedge, 1);
	jgrid = lev_interface.egrid(iedge, 2);
	if (igrid >= 0 && jgrid >= 0) 
	{
	    int kgrid = lev_interface.egrid(iedge, 0);
	    if (kgrid == -1)
		kgrid = lev_interface.egrid(iedge, 3);
	    if (kgrid != -1 && kgrid != igrid && kgrid != jgrid) 
	    {
		node_dirs(dir, lev_interface.edge(iedge).type());
		stridi1 = r[igrid].box().length(0);
		stridi2 = stridi1 * r[igrid].box().length(1);
		stridj1 = r[jgrid].box().length(0);
		stridj2 = stridj1 * r[jgrid].box().length(1);
		if (kgrid == lev_interface.egrid(iedge, 0)) 
		{
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[0],  1);
		    sstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		    dstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[1],  1);
		    sstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    dstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		}
		else 
		{
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[1], -1);
		    sstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		    dstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    b = lev_interface.node_edge(iedge);
		    b.shift(dir[0], -1);
		    sstartj = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    dstarti = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		}
		if ((nvals1 = b.length(0)) > 1) 
		{
		    stridi1 = 1;
		    stridj1 = 1;
		}
		else if ((nvals1 = b.length(2)) > 1) 
		{
		    stridi1 = stridi2;
		    stridj1 = stridj2;
		}
		else 
		{
		    nvals1 = b.length(1);
		}
		set(iset++, dstartj, sstarti, stridj1, 0, stridi1, 0, nvals1, 1);
		set(iset++, dstarti, sstartj, stridi1, 0, stridj1, 0, nvals1, 1);
	    }
	}
    }
#endif

    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	const int igrid = lev_interface.fgrid(iface, 0);
	const int jgrid = lev_interface.fgrid(iface, 1);
	if (igrid < 0 || jgrid < 0 || lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	const Box& b = lev_interface.node_face(iface);
#if (BL_SPACEDIM == 2)
	int dstarti, dstartj, sstarti, sstartj, stridi, stridj, nvals;
	if (lev_interface.fdim(iface) == 0) 
	{
	    nvals = b.length(1);
	    stridi = r[igrid].box().length(0);
	    stridj = r[jgrid].box().length(0);
	    dstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
		stridi * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
	    sstarti = dstarti - 2;
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
		stridj * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
	    sstartj = dstartj + 2;
	}
	else 
	{
	    nvals = b.length(0) + 2 * w;
	    stridi = r[igrid].box().length(0);
	    stridj = r[jgrid].box().length(0);
	    dstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
		stridi * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1));
	    sstarti = dstarti - 2 * stridi;
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
		stridj * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1));
	    sstartj = dstartj + 2 * stridj;
	    stridi = 1;
	    stridj = 1;
	}
	set(iset++, dstarti, sstartj, stridi, stridj, nvals);
	set(iset++, dstartj, sstarti, stridj, stridi, nvals);
#else
	int dstarti, dstartj, sstarti, sstartj;
	int stridi1, stridi2, stridj1, stridj2;
	if (lev_interface.fdim(iface) == 0) 
	{
	    int nvals1, nvals2;
	    nvals1 = b.length(1);
	    nvals2 = b.length(2);
	    stridi1 = r[igrid].box().length(0);
	    stridi2 = stridi1 * r[igrid].box().length(1);
	    stridj1 = r[jgrid].box().length(0);
	    stridj2 = stridj1 * r[jgrid].box().length(1);
	    dstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
		stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
		stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	    sstarti = dstarti - 2;
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
		stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	    sstartj = dstartj + 2;
	    set(iset++, dstarti, sstartj, stridi1, stridi2, stridj1, stridj2, nvals1, nvals2);
	    set(iset++, dstartj, sstarti, stridj1, stridj2, stridi1, stridi2, nvals1, nvals2);
	}
	else if (lev_interface.fdim(iface) == 1) 
	{
	    int il1 = 0, ih1 = 0, jl1 = 0, jh1 = 0;
	    if (r.box(jgrid).smallEnd(0) == b.smallEnd(0))
		jl1 = w;
	    if (r.box(jgrid).bigEnd(0) == b.bigEnd(0))
		jh1 = w;
	    if (r.box(igrid).smallEnd(0) == b.smallEnd(0))
		il1 = w;
	    if (r.box(igrid).bigEnd(0) == b.bigEnd(0))
		ih1 = w;
	    int nvals1i = b.length(0) + il1 + ih1;
	    int nvals1j = b.length(0) + jl1 + jh1;
	    int nvals2 = b.length(2);
	    stridi1 = r[igrid].box().length(0);
	    stridi2 = stridi1 * r[igrid].box().length(1);
	    stridj1 = r[jgrid].box().length(0);
	    stridj2 = stridj1 * r[jgrid].box().length(1);
	    dstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - il1 - r[igrid].box().smallEnd(0) +
		stridi1 * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1)) +
		stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	    sstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - jl1 - r[igrid].box().smallEnd(0) +
		stridi1 * (b.smallEnd(1) - 1 - r[igrid].box().smallEnd(1)) +
		stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - jl1 - r[jgrid].box().smallEnd(0) +
		stridj1 * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1)) +
		stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	    sstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - il1 - r[jgrid].box().smallEnd(0) +
		stridj1 * (b.smallEnd(1) + 1 - r[jgrid].box().smallEnd(1)) +
		stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	    stridi1 = 1;
	    stridj1 = 1;
	    set(iset++, dstarti, sstartj, stridi1, stridi2, stridj1, stridj2, nvals1i, nvals2);
	    set(iset++, dstartj, sstarti, stridj1, stridj2, stridi1, stridi2, nvals1j, nvals2);
	}
	else 
	{
	    int il1 = 0, ih1 = 0, jl1 = 0, jh1 = 0;
	    int il2 = 0, ih2 = 0, jl2 = 0, jh2 = 0;
	    if (r.box(jgrid).smallEnd(0) == b.smallEnd(0))
		jl1 = w;
	    if (r.box(jgrid).bigEnd(0) == b.bigEnd(0))
		jh1 = w;
	    if (r.box(igrid).smallEnd(0) == b.smallEnd(0))
		il1 = w;
	    if (r.box(igrid).bigEnd(0) == b.bigEnd(0))
		ih1 = w;
	    if (r.box(jgrid).smallEnd(1) == b.smallEnd(1))
		jl2 = w;
	    if (r.box(jgrid).bigEnd(1) == b.bigEnd(1))
		jh2 = w;
	    if (r.box(igrid).smallEnd(1) == b.smallEnd(1))
		il2 = w;
	    if (r.box(igrid).bigEnd(1) == b.bigEnd(1))
		ih2 = w;
	    int nvals1i = b.length(0) + il1 + ih1;
	    int nvals1j = b.length(0) + jl1 + jh1;
	    int nvals2i = b.length(1) + il2 + ih2;
	    int nvals2j = b.length(1) + jl2 + jh2;
	    stridi1 = r[igrid].box().length(0);
	    stridi2 = stridi1 * r[igrid].box().length(1);
	    stridj1 = r[jgrid].box().length(0);
	    stridj2 = stridj1 * r[jgrid].box().length(1);
	    dstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - il1 - r[igrid].box().smallEnd(0) +
		stridi1 * (b.smallEnd(1) - il2 - r[igrid].box().smallEnd(1)) +
		stridi2 * (b.smallEnd(2) + 1 - r[igrid].box().smallEnd(2));
	    sstarti = r[igrid].dataPtr() - baseptr +
		b.smallEnd(0) - jl1 - r[igrid].box().smallEnd(0) +
		stridi1 * (b.smallEnd(1) - jl2 - r[igrid].box().smallEnd(1)) +
		stridi2 * (b.smallEnd(2) - 1 - r[igrid].box().smallEnd(2));
	    dstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - jl1 - r[jgrid].box().smallEnd(0) +
		stridj1 * (b.smallEnd(1) - jl2 - r[jgrid].box().smallEnd(1)) +
		stridj2 * (b.smallEnd(2) - 1 - r[jgrid].box().smallEnd(2));
	    sstartj = r[jgrid].dataPtr() - baseptr +
		b.smallEnd(0) - il1 - r[jgrid].box().smallEnd(0) +
		stridj1 * (b.smallEnd(1) - il2 - r[jgrid].box().smallEnd(1)) +
		stridj2 * (b.smallEnd(2) + 1 - r[jgrid].box().smallEnd(2));
	    stridi2 = stridi1;
	    stridj2 = stridj1;
	    stridi1 = 1;
	    stridj1 = 1;
	    set(iset++, dstarti, sstartj, stridi1, stridi2, stridj1, stridj2, nvals1i, nvals2i);
	    set(iset++, dstartj, sstarti, stridj1, stridj2, stridi1, stridi2, nvals1j, nvals2j);
	}
#endif
  }
  assert( bdy != 0);
  bdy->set_border_cache(this, nsets, iset, r, lev_interface, w);
  nsets = iset;
}
