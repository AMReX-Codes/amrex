//BL_COPYRIGHT_NOTICE

#include <RunStats.H>

#include "fill_patch.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define FORT_FIPRODC   iprodc_
#define FORT_FIPRODN   iprodn_
#elif defined( BL_FORT_USE_UPPERCASE )
#define FORT_FIPRODC   IPRODC
#define FORT_FIPRODN   IPRODN
#elif defined( BL_FORT_USE_LOWERCASE )
#define FORT_FIPRODC   iprodc
#define FORT_FIPRODN   iprodn
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_FIPRODC(const Real*, intS, const Real*, intS, intS, Real*);
    void FORT_FIPRODN(const Real*, intS, const Real*, intS, intS, Real*);
}

Real
inner_product (const MultiFab& r,
               const MultiFab& s)
{
    BL_ASSERT(r.ok() && s.ok());
    BL_ASSERT(r.nComp() == 1);
    BL_ASSERT(s.nComp() == 1);
    BL_ASSERT(type(r) == type(s));
    
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
	BoxLib::Abort( "inner_product(): only supported for CELL- or NODE-based data" );
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

class task_bdy_fill : public task
{
public:

    task_bdy_fill (task_list&                tl_,
                   const amr_boundary_class* bdy_,
                   FArrayBox*                fab_,
                   int                       target_proc_id,
                   const Box&                region_,
                   const MultiFab&           smf_,
                   int                       grid_,
                   const Box&                domain_);

    virtual ~task_bdy_fill ();
    virtual bool ready ();
    virtual bool startup (long& sndcnt, long& rcvcnt);
    virtual bool work_to_do () const;
    virtual bool need_to_communicate (int& with) const;

private:
#ifdef BL_USE_MPI
    MPI_Request               m_request;
#endif
    const amr_boundary_class* m_bdy;
    FArrayBox*                m_fab;
    FArrayBox*                tmp;
    const Box                 m_region;
    const MultiFab&           m_smf;
    Box                       m_bx;
    const int                 m_sgrid;
    const Box&                m_domain;
    int                       m_target_proc_id;
    bool                      m_local;
};

task_bdy_fill::task_bdy_fill (task_list&                tl_,
                              const amr_boundary_class* bdy_,
                              FArrayBox*                fab_,
                              int                       target_proc_id,
                              const Box&                region_,
                              const MultiFab&           src_,
                              int                       grid_,
                              const Box&                domain_)
    :
    task(tl_),
    m_bdy(bdy_),
    m_fab(fab_),
    tmp(0),
    m_region(region_),
    m_smf(src_),
    m_sgrid(grid_),
    m_domain(domain_),
    m_target_proc_id(target_proc_id),
    m_local(false)
{
    BL_ASSERT(m_bdy != 0);

    m_bx = m_bdy->image(m_region,m_smf.boxArray()[m_sgrid],m_domain);

    if (m_fab != 0 && is_local(m_smf,m_sgrid))
    {
	m_local = true;

	m_bdy->fill(*m_fab, m_region, m_smf[m_sgrid], m_domain);
    }

    if (!m_smf.fabbox(m_sgrid).contains(m_bx))
        //
        // No work to do -- skip bad code in mixed_boundary_class::fill().
        //
        m_local = true;
}

task_bdy_fill::~task_bdy_fill ()
{
    delete tmp;
}

bool
task_bdy_fill::work_to_do () const
{
    return !m_local && (m_fab != 0 || (m_fab == 0 && is_local(m_smf,m_sgrid)));
}

bool
task_bdy_fill::need_to_communicate (int& with) const
{
    bool result = false;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (m_fab != 0)
        {
            with   = processor_number(m_smf,m_sgrid);
            result = true;
        }
        else if (m_fab == 0 && is_local(m_smf, m_sgrid)) 
        {
            with   = m_target_proc_id;
            result = true;
        }
    }
#endif

    return result;
}

bool
task_bdy_fill::startup (long& sndcnt, long& rcvcnt)
{
    m_started = true;

    bool result = true;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (m_fab != 0)
        {
            static RunStats irecv_stats("hg_irecv");
            tmp = new FArrayBox(m_bx,m_smf.nComp());
            irecv_stats.start();
            rcvcnt = tmp->box().numPts()*tmp->nComp();
            int res = MPI_Irecv(tmp->dataPtr(),
                                rcvcnt,
                                MPI_DOUBLE,
                                processor_number(m_smf,m_sgrid),
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            rcvcnt *= sizeof(double);
            irecv_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else if (m_fab == 0 && is_local(m_smf, m_sgrid)) 
        {
            static RunStats isend_stats("hg_isend");
            tmp = new FArrayBox(m_bx,m_smf.nComp());
            tmp->copy(m_smf[m_sgrid], m_bx);
            isend_stats.start();
            sndcnt = tmp->box().numPts()*tmp->nComp();
            int res = MPI_Isend(tmp->dataPtr(),
                                sndcnt,
                                MPI_DOUBLE,
                                m_target_proc_id,
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            sndcnt *= sizeof(double);
            isend_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else
        {
            result = false;
        }
    }
#endif

    return result;
}

bool
task_bdy_fill::ready ()
{
    BL_ASSERT(is_started());

    if (m_local) return true;

#ifdef BL_USE_MPI
    int flag, res;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);
    static RunStats test_stats("hg_test");
    test_stats.start();
    if ((res = MPI_Test(&m_request, &flag, &status)) != 0)
	ParallelDescriptor::Abort(res);
    test_stats.end();
    if (flag)
    {
	BL_ASSERT(m_request == MPI_REQUEST_NULL);
	if (m_fab)
	{
#ifndef NDEBUG
	    int count;
	    BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
	    BL_ASSERT(status.MPI_TAG    == m_sno);
	    if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
		ParallelDescriptor::Abort(res);
	    BL_ASSERT(count == tmp->box().numPts()*tmp->nComp());
#endif
	    m_bdy->fill(*m_fab, m_region, *tmp, m_domain);
	}
	return true;
    }
#endif

    return false;
}

task_fill_patch::task_fill_patch (task_list&                tl_,
                                  const MultiFab&           t_,
                                  int                       tt_,
                                  const Box&                region_,
                                  const MultiFab&           r_,
                                  const level_interface&    lev_interface_,
                                  const amr_boundary_class* bdy_,
                                  int                       idim_,
                                  int                       index_)
    :
    task_fab(tl_,t_,tt_,region_,r_.nComp()),
    r(r_),
    lev_interface(lev_interface_),
    bdy(bdy_),
    idim(idim_),
    index(index_)
{
    fill_patch();
}

task_fill_patch::~task_fill_patch () {}

bool
task_fill_patch::work_to_do () const
{
    return task_fab::work_to_do() || !dependencies.empty();
}

bool
task_fill_patch::fill_patch_blindly ()
{
    const BoxArray& r_ba = r.boxArray();

    for (int igrid = 0; igrid < r.length(); igrid++) 
    {
	if (is_local(r,igrid)) 
	{
	    BL_ASSERT(::grow(r[igrid].box(),-r.nGrow()) == r_ba[igrid]);
	}
	if (r_ba[igrid].contains(region)) 
	{
	    depend_on(m_task_list.add_task(new task_copy_local(m_task_list,
                                                               target,
                                                               target_proc_id(),
                                                               region,
                                                               r,
                                                               igrid)));
	    return true;
	}
    }
    for (int igrid = 0; igrid < r.length(); igrid++) 
    {
	if (is_local(r,igrid))
	{
	    BL_ASSERT(::grow(r[igrid].box(),-r.nGrow()) == r_ba[igrid]);
	}
	if (r_ba[igrid].intersects(region)) 
	{
            Box tb = r_ba[igrid] & region;
	    depend_on(m_task_list.add_task(new task_copy_local(m_task_list,
                                                               target,
                                                               target_proc_id(),
                                                               tb,
                                                               r,
                                                               igrid)));
	}
    }
    return false;
}

bool
task_fill_patch::fill_exterior_patch_blindly ()
{
    const BoxArray& em = lev_interface.exterior_mesh();

    for (int igrid = 0; igrid < em.length(); igrid++) 
    {
	const int jgrid = lev_interface.direct_exterior_ref(igrid);

	if (jgrid >= 0) 
	{
	    Box tb = em[igrid];
	    tb.convert(type(r));
	    if (tb.contains(region)) 
	    {
		depend_on(m_task_list.add_task(new task_bdy_fill(m_task_list,
                                                                 bdy,
                                                                 target,
                                                                 target_proc_id(),
                                                                 region,
                                                                 r,
                                                                 jgrid,
                                                                 lev_interface.domain())));
		return true;
	    }
	    if (tb.intersects(region)) 
	    {
		tb &= region;
		depend_on(m_task_list.add_task(new task_bdy_fill(m_task_list,
                                                                 bdy,
                                                                 target,
                                                                 target_proc_id(),
                                                                 tb,
                                                                 r,
                                                                 jgrid,
                                                                 lev_interface.domain())));
	    }
	}
    }

    return false;
}

void
task_fill_patch::fill_patch ()
{
    if (!region.ok()) return;

    if (target != 0)
    {
	BL_ASSERT(target->box() == region);
	BL_ASSERT(target->nComp() == r.nComp());
	BL_ASSERT(type(*target) == type(r));
    }
    BL_ASSERT(lev_interface.ok());
    BL_ASSERT(idim >= -1 && idim < BL_SPACEDIM);
    
    Box tdomain = lev_interface.domain();
    tdomain.convert(region.type());
    BL_ASSERT(target == 0 || type(*target) == region.type());
    Box idomain = ::grow(tdomain, IntVect::TheZeroVector() - type(r));

    if (idim == -1) 
    {
	if (idomain.contains(region) || bdy == 0) 
	{
	    fill_patch_blindly();
	}
	else if (!tdomain.intersects(region)) 
	{
	    fill_exterior_patch_blindly();
	}
	else if (idomain.intersects(region) && !fill_patch_blindly())
	{
            fill_exterior_patch_blindly();
	}
	else if (!fill_exterior_patch_blindly())
	{
            fill_patch_blindly();
	}
    }
    else
    {
        //
	// FIXME!!!
        //
        const BoxArray& r_ba = r.boxArray();
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
			    Box tb = r_ba[igrid] & region;
			    depend_on(m_task_list.add_task(new task_copy_local(m_task_list,
                                                                               target,
                                                                               target_proc_id(),
                                                                               tb,
                                                                               r,
                                                                               igrid)));
			}
			else 
			{
			    igrid = -2 - igrid;
			    Box tb = lev_interface.exterior_mesh()[igrid];
			    tb.convert(type(r));
			    tb &= region;
                            depend_on(m_task_list.add_task(new task_bdy_fill(m_task_list,
                                                                             bdy,
                                                                             target,
                                                                             target_proc_id(),
                                                                             tb,
                                                                             r,
                                                                             lev_interface.direct_exterior_ref(igrid),
                                                                             lev_interface.domain())));
			}
			break;
		    }
		}
	    }
	}
    }
}

static
void
sync_internal_borders (MultiFab&              r,
                       const level_interface& lev_interface)
{
    BL_ASSERT(type(r) == IntVect::TheNodeVector());

    task_list tl;
    for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
    {
	const int igrid = lev_interface.grid(level_interface::FACEDIM,iface,0);
	const int jgrid = lev_interface.grid(level_interface::FACEDIM,iface,1);
        //
	// Only do interior faces with fine grid on both sides.
        //
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
	    break;
	tl.add_task(new task_copy(tl,
                                  r,
                                  jgrid,
                                  r,
                                  igrid,
                                  lev_interface.node_box(level_interface::FACEDIM,iface)));
    }
#if (BL_SPACEDIM == 2)
    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
    {
	const int igrid = lev_interface.grid(0,icor,0);
	const int jgrid = lev_interface.grid(0,icor,3);
        //
	// Only do interior corners with fine grid on all sides.
        //
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.grid(0, icor, 1))
	    tl.add_task(new task_copy(tl,
                                      r,
                                      jgrid,
                                      r,
                                      igrid,
                                      lev_interface.box(0,icor)));
    }
#else
    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
    {
	const int igrid = lev_interface.grid(1,iedge,0);
	const int jgrid = lev_interface.grid(1,iedge,3);
        //
	// Only do interior edges with fine grid on all sides.
        //
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(1, iedge) != level_interface::ALL)
	    break;
	if (jgrid == lev_interface.grid(1, iedge, 1))
	    tl.add_task(new task_copy(tl,
                                      r,
                                      jgrid,
                                      r,
                                      igrid,
                                      lev_interface.node_box(1,iedge)));
    }
    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
    {
        int igrid = lev_interface.grid(0,icor,0);
        int jgrid = lev_interface.grid(0,icor,7);
        //
	// Only do interior corners with fine grid on all sides.
        //
	if (igrid < 0 || jgrid < 0 || lev_interface.geo(0, icor) != level_interface::ALL)
	    break;
	if (lev_interface.grid(0, icor, 3) == lev_interface.grid(0, icor, 1)) 
	{
	    if (jgrid != lev_interface.grid(0, icor, 3)) 
	    {
		tl.add_task(new task_copy(tl,
                                          r,
                                          jgrid,
                                          r,
                                          igrid,
                                          lev_interface.box(0,icor)));
		jgrid = lev_interface.grid(0, icor, 5);
		if (jgrid != lev_interface.grid(0, icor, 7))
		    tl.add_task(new task_copy(tl,
                                              r,
                                              jgrid,
                                              r,
                                              igrid,
                                              lev_interface.box(0,icor)));
	    }
	}
	else if (lev_interface.grid(0, icor, 5) == lev_interface.grid(0, icor, 1)) 
	{
	    if (jgrid != lev_interface.grid(0, icor, 5)) 
	    {
		tl.add_task(new task_copy(tl,
                                          r,
                                          jgrid,
                                          r,
                                          igrid,
                                          lev_interface.box(0,icor)));
		jgrid = lev_interface.grid(0, icor, 3);
		if (jgrid != lev_interface.grid(0, icor, 7)) 
		{
		    tl.add_task(new task_copy(tl,
                                              r,
                                              jgrid,
                                              r,
                                              igrid,
                                              lev_interface.box(0,icor)));
		    if (jgrid == lev_interface.grid(0, icor, 2)) 
		    {
			jgrid = lev_interface.grid(0, icor, 6);
			if (jgrid != lev_interface.grid(0, icor, 7))
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      jgrid,
                                                      r,
                                                      igrid,
                                                      lev_interface.box(0,icor)));
		    }
		}
	    }
	}
    }
#endif
    tl.execute("sync_internal_borders");
}

void
sync_borders (MultiFab&                 r,
              const level_interface&    lev_interface,
              const amr_boundary_class* bdy)
{
    sync_internal_borders(r, lev_interface);
    BL_ASSERT(bdy != 0);
    bdy->sync_borders(r, lev_interface);
}

#if BL_SPACEDIM == 3
//
// Local function used only by fill_internal_borders:
//
inline
void
node_dirs (int            dir[2],
           const IntVect& typ)
{
    if (typ[0] == IndexType::NODE) 
    {
	dir[0] = 0;
        dir[1] = (typ[1] == IndexType::NODE) ? 1 : 2;
    }
    else 
    {
	dir[0] = 1;
	dir[1] = 2;
    }
}
#endif

//
// The sequencing used in fill_internal_borders, fcpy2 and set_border_cache
// (narrow x, medium y, wide z) is necessary to avoid overwrite problems
// like those seen in the sync routines.  Boundary copies are all wide
// regardless of direction and come after interior copies---overwrite
// difficulties are avoided since grids can't bridge a boundary.

// Modifications are necessary in 3D to deal with lack of diagonal
// communication across edges at the coarse-fine lev_interface.  These
// modifications take the form of narrowing certain copies to avoid
// overwriting good values with bad ones.
//

static
void
fill_internal_borders (MultiFab&              r,
                       const level_interface& lev_interface,
                       int                    w,
                       bool                   hg_terrain)
{
    BL_ASSERT(type(r) == IntVect::TheCellVector() || type(r) == IntVect::TheNodeVector() );

    w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;

    BL_ASSERT(w == 1 || w == 0);

    task_list tl;
    if (type(r) == IntVect::TheNodeVector()) 
    {
#if (BL_SPACEDIM == 3)
	if(hg_terrain)
	{
            //
	    // Attempt to deal with corner-coupling problem with 27-point stencils
            //
	    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
	    {
		if (lev_interface.geo(1, iedge) == level_interface::ALL)
		    continue;
                int igrid = lev_interface.grid(1,iedge,0);
                int jgrid = lev_interface.grid(1,iedge,3);
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
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      jgrid,
                                                      r,
                                                      igrid,
                                                      b.shift(dir[0],-1)));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      igrid,
                                                      r,
                                                      jgrid,
                                                      b.shift(dir[1],1)));
			}
			else 
			{
			    Box b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      jgrid,
                                                      r,
                                                      igrid,
                                                      b.shift(dir[1],-1)));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      igrid,
                                                      r,
                                                      jgrid,
                                                      b.shift(dir[0],1)));
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
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      jgrid,
                                                      r,
                                                      igrid,
                                                      b.shift(dir[0],1)));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      igrid,
                                                      r,
                                                      jgrid,
                                                      b.shift(dir[1],1)));
			}
			else 
			{
			    Box b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      jgrid,
                                                      r,
                                                      igrid,
                                                      b.shift(dir[1],-1)));
			    b = lev_interface.node_box(1, iedge);
			    tl.add_task(new task_copy(tl,
                                                      r,
                                                      igrid,
                                                      r,
                                                      jgrid,
                                                      b.shift(dir[0],-1)));
			}
		    }
		}
	    }
	}
#endif
        const BoxArray& r_ba = r.boxArray();

	for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	{
	    const int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	    const int jgrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
	    if (igrid < 0 || jgrid < 0 || lev_interface.geo(level_interface::FACEDIM, iface) != level_interface::ALL)
		break;
	    const Box& b = lev_interface.node_box(level_interface::FACEDIM, iface);
	    // tl.add_task(new task_copy_2(r, igrid, r, jgrid, b, w));
            const int idim = lev_interface.fdim(iface);
            Box bj = lev_interface.node_box(level_interface::FACEDIM, iface);
            Box bi = lev_interface.node_box(level_interface::FACEDIM, iface);
            for (int i = 0; i < idim; i++) 
            {
                if (r_ba[jgrid].smallEnd(i) == bj.smallEnd(i)) bj.growLo(i,w);
                if (r_ba[jgrid].bigEnd(i)   == bj.bigEnd(i))   bj.growHi(i,w);
                if (r_ba[igrid].smallEnd(i) == bi.smallEnd(i)) bi.growLo(i,w);
                if (r_ba[igrid].bigEnd(i)   == bi.bigEnd(i))   bi.growHi(i,w);
            }
            bj.shift(idim, -1).growLo(idim, w-1);
            bi.shift(idim,  1).growHi(idim, w-1);
#if 1
	    tl.add_task(new task_copy(tl,r,jgrid,r,igrid,bj));
	    tl.add_task(new task_copy(tl,r,igrid,r,jgrid,bi));
#else
	    tl.add_task(new task_copy(tl,r,jgrid,r,igrid,w_shift(b,r_ba[jgrid],r.nGrow(), -w)));
	    tl.add_task(new task_copy(tl,r,igrid,r,jgrid,w_shift(b,r_ba[igrid],r.nGrow(),  w)));
#endif
	}
    }
    else if (type(r) == IntVect::TheCellVector()) 
    {
        const BoxArray& r_ba = r.boxArray();

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
	    tl.add_task(new task_copy(tl, r, jgrid, r, igrid, b));
	    tl.add_task(new task_copy(tl, r, igrid, r, jgrid, b.shift(idim, w)));
#else
	    Box bj = lev_interface.box(level_interface::FACEDIM, iface);
	    Box bi = lev_interface.box(level_interface::FACEDIM, iface);
	    for (int i = 0; i < idim; i++) 
	    {
		if (r_ba[jgrid].smallEnd(i) == bj.smallEnd(i)) bj.growLo(i,w);
		if (r_ba[jgrid].bigEnd(i)   == bj.bigEnd(i))   bj.growHi(i,w);
		if (r_ba[igrid].smallEnd(i) == bi.smallEnd(i)) bi.growLo(i,w);
		if (r_ba[igrid].bigEnd(i)   == bi.bigEnd(i))   bi.growHi(i,w);
	    }
	    bj.growLo(idim, w).convert(IntVect::TheCellVector());
	    bi.growHi(idim, w).convert(IntVect::TheCellVector());
	    tl.add_task(new task_copy(tl,r,jgrid,r,igrid,bj));
	    tl.add_task(new task_copy(tl,r,igrid,r,jgrid,bi));
#endif
	}
    }
    tl.execute("fill_internal_borders");
}

void
fill_borders (MultiFab&                 r,
              const level_interface&    lev_interface,
              const amr_boundary_class* bdy,
              int                       w,
              bool                      hg_terrain)
{
    HG_TEST_NORM(r, "fill_borders 0");
    fill_internal_borders(r, lev_interface, w, hg_terrain);
    HG_TEST_NORM(r, "fill_borders 1");
    BL_ASSERT(bdy != 0);
    bdy->fill_borders(r, lev_interface, w);
    HG_TEST_NORM(r, "fill_borders 2");
}

void
clear_part_interface (MultiFab&              r,
                      const level_interface& lev_interface)
{
    BL_ASSERT(r.nComp() == 1);
    BL_ASSERT(type(r) == IntVect::TheNodeVector());

    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	for (int ibox = 0; ibox < lev_interface.nboxes(i); ibox++) 
	{
            //
	    // coarse-fine face contained in part_fine grid, or orphan edge/corner
            //
	    const int igrid = lev_interface.aux(i, ibox);
	    if (igrid < 0  || is_remote(r, igrid))
                continue;
	    BL_ASSERT(is_local(r, igrid));
	    r[igrid].setVal(0.0, lev_interface.node_box(i, ibox), 0);
	}
    }
}

class task_restric_fill : public task
{
public:

    task_restric_fill (task_list&                  tl_,
                       const amr_restrictor_class& restric,
                       MultiFab&                   dest,
                       int                         dgrid,
                       const MultiFab&             r,
                       int                         rgrid,
                       const Box&                  box,
                       const IntVect&              rat);

    virtual ~task_restric_fill ();
    virtual bool ready ();
    virtual void hint () const;
    virtual bool startup (long& sndcnt, long& rcvcnt);
    virtual bool work_to_do () const;
    virtual bool need_to_communicate (int& with) const;

private:
    //
    // The data.
    //
#ifdef BL_USE_MPI
    MPI_Request                 m_request;
#endif
    const amr_restrictor_class& m_restric;
    FArrayBox*                  m_tmp;
    MultiFab&                   m_d;
    const MultiFab&             m_r;
    const int                   m_dgrid;
    const int                   m_rgrid;
    const Box                   m_box;
    const IntVect               m_rat;
    bool                        m_local;
};

task_restric_fill::task_restric_fill (task_list&                  tl_,
                                      const amr_restrictor_class& restric,
                                      MultiFab&                   dest,
                                      int                         dgrid,
                                      const MultiFab&             r,
                                      int                         rgrid,
                                      const Box&                  box,
                                      const IntVect&              rat)
    :
    task(tl_),
    m_restric(restric),
    m_d(dest),
    m_r(r),
    m_dgrid(dgrid),
    m_rgrid(rgrid),
    m_tmp(0),
    m_box(box),
    m_rat(rat),
    m_local(false)
{
    if (is_local(m_d, m_dgrid) && is_local(m_r, m_rgrid))
    {
	m_local = true;

	m_restric.fill(m_d[m_dgrid], m_box, m_r[m_rgrid], m_rat);
    }
}

task_restric_fill::~task_restric_fill ()
{
    delete m_tmp;
}

bool
task_restric_fill::work_to_do () const
{
    return !m_local && (is_local(m_d,m_dgrid) || is_local(m_r,m_rgrid));
}

bool
task_restric_fill::need_to_communicate (int& with) const
{
    bool result = false;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (is_local(m_d, m_dgrid))
        {
            with   = processor_number(m_r,m_rgrid);
            result = true;
        }
        else if (is_local(m_r, m_rgrid))
        {
            with   = processor_number(m_d,m_dgrid);
            result = true;
        }
    }
#endif

    return result;
}

bool
task_restric_fill::startup (long& sndcnt, long& rcvcnt)
{
    m_started = true;

    bool result = true;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (is_local(m_d, m_dgrid))
        {
            static RunStats irecv_stats("hg_irecv");
            m_tmp = new FArrayBox(m_box, m_d.nComp());
            irecv_stats.start();
            rcvcnt = m_tmp->box().numPts()*m_tmp->nComp();
            int res = MPI_Irecv(m_tmp->dataPtr(),
                                rcvcnt,
                                MPI_DOUBLE,
                                processor_number(m_r,m_rgrid),
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            rcvcnt *= sizeof(double);
            irecv_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_r, m_rgrid))
        {
            static RunStats isend_stats("hg_isend");
            m_tmp = new FArrayBox(m_box, m_d.nComp());
	    m_restric.fill(*m_tmp, m_box, m_r[m_rgrid], m_rat);
            isend_stats.start();
            sndcnt = m_tmp->box().numPts()*m_tmp->nComp();
            int res = MPI_Isend(m_tmp->dataPtr(),
                                sndcnt,
                                MPI_DOUBLE,
                                processor_number(m_d,m_dgrid),
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            sndcnt *= sizeof(double);
            isend_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else
        {
            result = false;
        }
    }
#endif

    return result;
}

bool
task_restric_fill::ready ()
{
    BL_ASSERT(is_started());

    if (m_local) return true;

#ifdef BL_USE_MPI
    int flag, res;
    MPI_Status status;
    static RunStats test_stats("hg_test");
    test_stats.start();
    if ((res = MPI_Test(&m_request, &flag, &status)) != 0)
	ParallelDescriptor::Abort(res);
    test_stats.end();
    if (flag)
    {
	if (is_local(m_d, m_dgrid))
	{
            m_d[m_dgrid].copy(*m_tmp);
	}
	return true;
    }
#endif

    return false;
}

void
task_restric_fill::hint () const
{
    task::_hint();
    if (is_local(m_r, m_rgrid) && is_local(m_d, m_dgrid))
    {
	HG_DEBUG_OUT( "L" );
    }
    else if (is_local(m_r, m_rgrid))
    {
	HG_DEBUG_OUT( "S" );
    }
    else if (is_local(m_d, m_dgrid))
    {
    	HG_DEBUG_OUT( "R" );
    }
    else
    {
	HG_DEBUG_OUT( "?" );
    }
    HG_DEBUG_OUT( ' ' << m_box  << ' ' << m_dgrid << ' '; );
    HG_DEBUG_OUT( ")" << endl );
}

void
restrict_level (MultiFab&                   dest,
                MultiFab&                   r,
                const IntVect&              rat,
                const amr_restrictor_class& restric,
                const level_interface&      lev_interface,
                const amr_boundary_class*   bdy)
{
    BL_ASSERT(type(dest) == type(r));

    HG_TEST_NORM( dest, "restrict_level a");

    const BoxArray& r_ba    = r.boxArray();
    const BoxArray& dest_ba = dest.boxArray();

    task_list tl;
    for (int jgrid = 0; jgrid < dest.length(); jgrid++) 
    {
	const Box& region = dest_ba[jgrid];

	for (int igrid = 0; igrid < r.length(); igrid++) 
	{
	    Box cbox = restric.box(r_ba[igrid], rat);

	    if (region.intersects(cbox)) 
	    {
		cbox &= region;
		tl.add_task(new task_restric_fill(tl,
                                                  restric,
                                                  dest,
                                                  jgrid,
                                                  r,
                                                  igrid,
                                                  cbox,
                                                  rat));
	    }
	}
    }
    tl.execute("restrict_level");

    if (lev_interface.ok())
    {
	restric.fill_interface( dest, r, lev_interface, bdy, rat);
    }

    HG_TEST_NORM( dest, "restrict_level a1");
}
