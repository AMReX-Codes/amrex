
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_FIPRODC   iprodc_
#define FORT_FIPRODN   iprodn_
#define FORT_FFCPY2    fcpy2_
#else
#define FORT_FIPRODC   IPRODC
#define FORT_FIPRODN   IPRODN
#define FORT_FFCPY2    FCPY2
#endif

extern "C"
{
    void FORT_FIPRODC(const Real*, intS, const Real*, intS, intS, Real*);
    void FORT_FIPRODN(const Real*, intS, const Real*, intS, intS, Real*);
#if (BL_SPACEDIM == 2)
    void FORT_FFCPY2(Real*, intS, const Real*, intS, intS, const int*, const int&);
#else
    void FORT_FFCPY2(Real*, intS, const Real*, intS, intS, const int*, const int*, const int&);
#endif
}

void internal_copy(MultiFab& r, int destgrid, int srcgrid, const Box& b) 
{
    r[destgrid].copy(r[srcgrid], b);
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
    for (int igrid = 0; igrid < r.length(); igrid++) 
    {
	if (r.box(igrid).contains(region))
	    return igrid;
    }
    return -1;
}


class task_fab : public task
{
public:
    virtual const FArrayBox& fab() const = 0;
};

class task_fab_get : public task_fab
{
public:
    task_fab_get(const MultiFab& r_, int grid_);
    virtual ~task_fab_get();
    virtual const FArrayBox& fab() const;
    virtual bool ready();
private:
    const MultiFab& r;
    const int grid;
};

task_fab_get::task_fab_get(const MultiFab& r_, int grid_) : r(r_), grid(grid_) {}
const FArrayBox& task_fab_get::fab() const
{
    return r[grid];
}
task_fab_get::~task_fab_get()
{
}

bool task_fab_get::ready()
{
    return true;
}

class task_fill_patch : public task_fab
{
public:
    task_fill_patch(FArrayBox& fab_, const Box& region_,
	const MultiFab& r_, const level_interface& lev_interface_, const amr_boundary_class* bdy_, int idim_, int index_);
    task_fill_patch(const Box& region_, int ncomp_, 
	const MultiFab& r_, const level_interface& lev_interface_, const amr_boundary_class* bdy_, int idim_, int index_);
    virtual ~task_fill_patch();
    virtual const FArrayBox& fab() const;
    virtual bool ready();
private:
    bool fill_patch_blindly();
    bool fill_exterior_patch_blindly();
    void fill_patch();
    bool newed;
    FArrayBox* target;
    const Box region;
    const MultiFab& r;
    const level_interface& lev_interface;
    const amr_boundary_class* bdy;
    const int idim;
    const int index;
    task_list tl;
};

bool
task_fill_patch::fill_patch_blindly()
{
    for (int igrid = 0; igrid < r.length(); igrid++) 
    {
	Box tb = grow(r[igrid].box(), -r.nGrow());
	if (tb.contains(region)) 
	{
	    tl.add_task(new task_copy(target, r, igrid, region));
	    // target->copy(r[igrid], region);
	    return true;
	}
    }
    for (int igrid = 0; igrid < r.length(); igrid++) 
    {
	Box tb = grow(r[igrid].box(), -r.nGrow());
	if (tb.intersects(region)) 
	{
	    tb &= region;
	    tl.add_task(new task_copy(target, r, igrid, tb));
	    // target->copy(r[igrid], tb);
	}
    }
    return false;
}

class task_bdy_fill : public task
{
public:
    task_bdy_fill(const amr_boundary_class* bdy_, FArrayBox& fab_, const Box& region_, const MultiFab& src_, int grid_, const Box& domain_);
    virtual bool ready();
private:
    const amr_boundary_class* bdy;
    FArrayBox& fab;
    const Box region;
    const MultiFab& src;
    const int grid;
    const Box& domain;
};

task_bdy_fill::task_bdy_fill(const amr_boundary_class* bdy_, FArrayBox& fab_, const Box& region_, const MultiFab& src_, int grid_, const Box& domain_)
: bdy(bdy_), fab(fab_), region(region_), src(src_), grid(grid_), domain(domain_)
{
}

bool task_bdy_fill::ready()
{
    bdy->fill(fab, region, src[grid], domain);
    return true;
}

bool task_fill_patch::fill_exterior_patch_blindly()
{
    const BoxArray& em = lev_interface.exterior_mesh();
    for (int igrid = 0; igrid < em.length(); igrid++) 
    {
	int jgrid = lev_interface.direct_exterior_ref(igrid);
	if (jgrid >= 0) 
	{
	    assert(bdy != 0);
	    Box tb = em[igrid];
	    tb.convert(type(r));
	    if (tb.contains(region)) 
	    {
		tl.add_task(new task_bdy_fill(bdy, *target, region, r, jgrid, lev_interface.domain()));
		// bdy->fill(*target, region, r[jgrid], lev_interface.domain());
		return true;
	    }
	    if (tb.intersects(region)) 
	    {
		tb &= region;
		tl.add_task(new task_bdy_fill(bdy, *target, tb, r, jgrid, lev_interface.domain()));
		//bdy->fill(*target, tb, r[jgrid], lev_interface.domain());
	    }
	}
    }
    return false;
}

void task_fill_patch::fill_patch()
{
    if ( !region.ok() ) return;

    assert(target->box() == region);
    assert(target->nComp() == r.nComp());
    assert(type(*target) == type(r));
    assert(lev_interface.ok());
    assert( idim >= -1 && idim < BL_SPACEDIM );
    
    Box tdomain = lev_interface.domain();
    tdomain.convert(type(*target));
    Box idomain = grow(tdomain, IntVect::TheZeroVector() - type(r));
    
    if (idim == -1 ) 
    {
	if (idomain.contains(region) || bdy == 0) 
	{
	    fill_patch_blindly();
	}
	else if (!tdomain.intersects(region)) 
	{
	    fill_exterior_patch_blindly();
	}
	else if (idomain.intersects(region)) 
	{
	    if ( !fill_patch_blindly() )
		fill_exterior_patch_blindly();
	}
	else 
	{
	    if ( !fill_exterior_patch_blindly() )
		fill_patch_blindly();
	}
    }
    else
    {
	Array<int> gridnum(lev_interface.ngrids(idim)+1);
	gridnum[0] = -1;
	for (int i = 0; i < lev_interface.ngrids(idim); i++) 
	{
	    int igrid = lev_interface.grid(idim, index, i);
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
			    target->copy(r[igrid], tb);
			}
			else 
			{
			    igrid = -2 - igrid;
			    Box tb = lev_interface.exterior_mesh()[igrid];
			    tb.convert(type(r));
			    tb &= region;
			    assert( bdy != 0 );
			    bdy->fill(*target, tb, r[lev_interface.direct_exterior_ref(igrid)], lev_interface.domain());
			}
			break;
		    }
		}
	    }
	}
    }
}

task_fill_patch::task_fill_patch(const Box& region_, int ncomp_,
				 const MultiFab& r_,
				 const level_interface& lev_interface_,
				 const amr_boundary_class* bdy_,
				 int idim_, int index_)
				 : target(0), region(region_),
				 r(r_), lev_interface(lev_interface_), bdy(bdy_), idim(idim_), index(index_), newed(true)
{
    target = new FArrayBox(region, ncomp_);
}

task_fill_patch::task_fill_patch(FArrayBox& fab_, const Box& region_,
				 const MultiFab& r_,
				 const level_interface& lev_interface_,
				 const amr_boundary_class* bdy_,
				 int idim_, int index_)
				 : target(0), region(region_),
				 r(r_), lev_interface(lev_interface_), bdy(bdy_), idim(idim_), index(index_), newed(false)
{
    target = &fab_;
}


task_fill_patch::~task_fill_patch()
{
    if ( newed) delete target;
}

bool task_fill_patch::ready()
{
    fill_patch();
    tl.execute();
    return true;
}

const FArrayBox& task_fill_patch::fab() const
{
    return *target;
}

void fill_patch(FArrayBox& patch, const Box& region,
	       const MultiFab& r,
	       const level_interface& lev_interface,
	       const amr_boundary_class* bdy,
	       int idim, int index)
{
    task_fill_patch* t = new task_fill_patch(patch, region, r, lev_interface, bdy, idim, index);
    assert(t->ready());
    delete t;
    return;
}

static void sync_internal_borders(MultiFab& r, const level_interface& lev_interface)
{
    assert(type(r) == IntVect::TheNodeVector());

    task_list tl;
    for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
    {
	int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	int jgrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
	// only do interior faces with fine grid on both sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
	    break;
	tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.node_box(level_interface::FACEDIM, iface)));
    }
#if (BL_SPACEDIM == 2)
    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
    {
	int igrid = lev_interface.grid(0, icor, 0);
	int jgrid = lev_interface.grid(0, icor, 3);
	// only do interior corners with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.grid(0, icor, 1))
	    tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.box(0, icor)));
    }
#else
    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
    {
	int igrid = lev_interface.grid(1, iedge, 0);
	int jgrid = lev_interface.grid(1, iedge, 3);
	// only do interior edges with fine grid on all sides
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(1, iedge) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.grid(1, iedge, 1))
	    tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.node_box(1, iedge)));
    }
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
		tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.box(0, icor)));
		jgrid = lev_interface.grid(0, icor, 5);
		if (jgrid != lev_interface.grid(0, icor, 7))
		    tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.box(0, icor)));
	    }
	}
	else if (lev_interface.grid(0, icor, 5) == lev_interface.grid(0, icor, 1)) 
	{
	    if (jgrid != lev_interface.grid(0, icor, 5)) 
	    {
		tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.box(0, icor)));
		jgrid = lev_interface.grid(0, icor, 3);
		if (jgrid != lev_interface.grid(0, icor, 7)) 
		{
		    tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.box(0, icor)));
		    if (jgrid == lev_interface.grid(0, icor, 2)) 
		    {
			jgrid = lev_interface.grid(0, icor, 6);
			if (jgrid != lev_interface.grid(0, icor, 7))
			    tl.add_task(new task_copy(&r, jgrid, r, igrid, lev_interface.box(0, icor)));
		    }
		}
	    }
	}
    }
#endif
    tl.execute();
}

void sync_borders(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy)
{
    sync_internal_borders(r, lev_interface);
    assert(bdy != 0);
    bdy->sync_borders(r, lev_interface);
}

#if BL_SPACEDIM == 3
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
#endif


class task_copy_2 : public task
{
public:
    task_copy_2(MultiFab* r1, int i1, const MultiFab& r2, int i2, const Box& bx, int w)
	: m_r1(r1), m_i1(i1), m_r2(r2), m_i2(i2), m_bx(bx), m_w(w)
    {
	    Real* ptra = (*m_r1)[m_i1].dataPtr();
	    const Real* ptrb = m_r2[m_i2].dataPtr();
	    const Box& boxa = (*m_r1)[m_i1].box();
	    const Box& boxb = m_r2[m_i2].box();
#if (BL_SPACEDIM == 2)
	    FORT_FFCPY2(ptra, DIMLIST(boxa), ptrb, DIMLIST(boxb), DIMLIST(m_bx), &m_w, m_r1->nComp());
#else
	    const int ibord = m_r1->nGrow();
	    FORT_FFCPY2(ptra, DIMLIST(boxa), ptrb, DIMLIST(boxb), DIMLIST(m_bx), &m_w, &ibord, m_r1->nComp());
#endif
    }
    virtual bool ready()
    {
	return true;
    }
private:
    MultiFab* m_r1;
    const int m_i1;
    const MultiFab& m_r2;
    const int m_i2;
    const Box m_bx;
    const int m_w;
};

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
    assert( w == 1 || w == 0 );

    task_list tl;
    if ( type(r) == IntVect::TheNodeVector() ) 
    {
#if (BL_SPACEDIM == 3)
	if(hg_terrain)
	{
	    // attempt to deal with corner-coupling problem with 27-point stencils
	    task_list tl;
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
			int dir[2];
			node_dirs(dir, lev_interface.box(1, iedge).type());
			if (kgrid == lev_interface.grid(1, iedge, 1)) 
			{
			    Box b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[0], -1));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(&r, igrid, r, jgrid, b.shift(dir[1],  1)));
			}
			else 
			{
			    Box b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[1], -1));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(&r, igrid, r, jgrid, b.shift(dir[0],  1)));
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
			int dir[2];
			node_dirs(dir, lev_interface.box(1, iedge).type());
			if (kgrid == lev_interface.grid(1, iedge, 0)) 
			{
			    Box b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[0],  1));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(&r, igrid, r, jgrid, b.shift(dir[1],  1)));
			}
			else 
			{
			    Box b = lev_interface.node_box(1, iedge);
			    internal_copy(r, jgrid, igrid, b.shift(dir[1], -1));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(&r, igrid, r, jgrid, b.shift(dir[0], -1)));
			}
		    }
		}
	    }
	}
#endif
	for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	{
	    const int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	    const int jgrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
	    if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
		break;
	    const Box& b = lev_interface.node_box(level_interface::FACEDIM, iface);
	    tl.add_task(new task_copy_2(&r, igrid, r, jgrid, b, w));
	}
    }
    else if (type(r) == IntVect::TheCellVector()) 
    {
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
	    tl.add_task(new task_copy(&r, jgrid, r, igrid, b));
	    tl.add_task(new task_copy(&r, igrid, r, jgrid, b.shift(idim, w)));
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
	    tl.add_task(new task_copy(&r, jgrid, r, igrid, bj));
	    tl.add_task(new task_copy(&r, igrid, r, jgrid, bi));
#endif
	}
    }
    tl.execute();
}

void fill_borders(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy, int w, bool hg_terrain)
{
    fill_internal_borders(r, lev_interface, w, hg_terrain);
    assert(bdy != 0);
    bdy->fill_borders(r, lev_interface, w);
}

void clear_part_interface(MultiFab& r, const level_interface& lev_interface)
{
    assert(r.nComp() == 1);
    assert(type(r) == IntVect::TheNodeVector());

    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	for (int ibox = 0; ibox < lev_interface.nboxes(i); ibox++) 
	{
	    // coarse-fine face contained in part_fine grid, or orphan edge/corner
	    int igrid = lev_interface.aux(i, ibox);
	    if ( igrid < 0  || is_remote(r, igrid) ) continue;
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
    const int jgrid = find_patch(cb, r);
    if (jgrid == -1) 
    {
	FArrayBox* cgr = new FArrayBox(cb, r.nComp());
	fill_patch(*cgr, cb, r, lev_interface, 0);
	interp.fill(patch, region, *cgr, cb, rat);
	delete cgr;
    }
    else 
    {
	interp.fill(patch, region, r[jgrid], cb, rat);
    }
}

class task_restric_fill : public task
{
public:
    task_restric_fill(const amr_restrictor_class& restric,
	MultiFab& dest, int dgrid, MultiFab& r, int rgrid, const Box& box, const IntVect& rat)
	: m_restric(restric), m_dest(dest), m_dgrid(dgrid), m_rgrid(rgrid), m_r(r), m_box(box), m_rat(rat)
    {
	m_restric.fill(m_dest[m_dgrid], m_box, m_r[m_rgrid], m_rat);
    }
    virtual bool ready();
private:
    const amr_restrictor_class& m_restric;
    MultiFab& m_dest;
    const int m_dgrid;
    MultiFab& m_r;
    const int m_rgrid;
    const Box m_box;
    const IntVect m_rat;
};

bool
task_restric_fill::ready()
{
    return true;
}

void restrict_level(MultiFab& dest, 
		    MultiFab& r, const IntVect& rat,
		    const amr_restrictor_class& restric,
		    const level_interface& lev_interface,
		    const amr_boundary_class* bdy)
{
    assert(type(dest) == type(r));
     task_list tl;
    for (int jgrid = 0; jgrid < dest.length(); jgrid++) 
    {
	const Box& region = dest.box(jgrid);
	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    Box cbox = r.box(igrid);
	    cbox = restric.box(cbox, rat);
	    if (region.intersects(cbox)) 
	    {
		cbox &= region;
		// restric.fill(dest[jgrid], cbox, r[igrid], rat);
		tl.add_task(new task_restric_fill(restric, dest, jgrid, r, igrid, cbox, rat));
	    }
	}
    }
    tl.execute();
    if ( lev_interface.ok() )
    {
	restric.fill_interface( dest, r, lev_interface, bdy, rat);
    }
}

