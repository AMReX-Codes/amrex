
#include "hg_multi.H"

#ifdef BL_FORT_USE_UNDERSCORE
#define   FORT_HGFRES		hgfres_
#define   FORT_HGFRES_TERRAIN   hgfres_terrain_
#define   FORT_HGFRES_FULL	hgfres_full_
#define   FORT_HGERES		hgeres_
#define   FORT_HGERES_TERRAIN   hgeres_terrain_
#define   FORT_HGCRES		hgcres_
#define   FORT_HGCRES_TERRAIN   hgcres_terrain_
#define   FORT_HGORES		hgores_full_
#define   FORT_HGIRES		hgires_full_
#define   FORT_HGDRES		hgdres_full_
#else
#define   FORT_HGFRES		HGFRES
#define   FORT_HGFRES_TERRAIN   HGFRES_TERRAIN
#define   FORT_HGFRES_FULL	HGFRES_FULL
#define   FORT_HGERES		HGERES
#define   FORT_HGERES_TERRAIN   HGERES_TERRAIN
#define   FORT_HGCRES		HGCRES
#define   FORT_HGCRES_TERRAIN   HGCRES_TERRAIN
#define   FORT_HGORES		HGORES_FULL
#define   FORT_HGIRES		HGIRES_FULL
#define   FORT_HGDRES		HGDRES_FULL
#endif

extern "C" 
{
    
#if (BL_SPACEDIM == 1)
#error not relevant
#endif
    void FORT_HGFRES_TERRAIN(Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGFRES        (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGCRES_TERRAIN(Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGCRES        (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGERES_TERRAIN(Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGERES        (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);
    void FORT_HGFRES_FULL   (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*, const int*, const int*);
#if (BL_SPACEDIM == 2)
    void FORT_HGIRES        (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*, const int*, const int*);
    void FORT_HGDRES        (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*, const int*, const int*);
    void FORT_HGORES        (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*, const int*, const int*);
#endif
}

typedef void (*OFDIRES)(Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*, const int*, const int*);
typedef void (*FCERES) (Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, const Real*, intS, intS, CRealPS, intRS, const int*, const int*);

struct task_ofdires : public task
{
private:
    OFDIRES func;
};

class task_fceres_2 : public task
{
public:
    task_fceres_2(FCERES f_, MultiFab& r_, const MultiFab& s_, const MultiFab& d_, const MultiFab& sg_, int igrid_,
	task_fab* c_, task_fab* sc_, const Box& creg_, const Real h_[BL_SPACEDIM], const IntVect& rat_, int idim_, int idir_)
	: func(f_), r(r_), s(s_), d(d_), igrid(igrid_), c(c_), sg(sg_), sc(sc_), creg(creg_), rat(rat_), idim(idim_), idir(idir_)
    {
	for(int i = 0; i < BL_SPACEDIM; ++i)
	{
	    h[i] = h_[i];
	}
    }
    virtual bool ready()
    {
	abort(); return false;
    };
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	abort(); return false;
    }
private:
    FCERES func;
    MultiFab& r;
    const MultiFab& s;
    const MultiFab& d;
    const MultiFab& sg;
    const int igrid;
    task_fab* c;
    task_fab* sc;
    const Box creg;
    Real h[BL_SPACEDIM];
    const IntVect rat;
    int idim;
    int idir;
};

class task_fceres_4 : public task
{
public:
    task_fceres_4(FCERES f_, MultiFab& r_, const MultiFab& s_, int igrid_,
	task_fab* c_, task_fab* sc_, task_fab* sigmaf_, task_fab* sigmac_, const Box& creg_, const Real h_[BL_SPACEDIM], const IntVect& rat_, const Array<int>& ga_, const IntVect& t_ = IntVect())
	: func(f_), r(r_), s(s_), igrid(igrid_), c(c_), sc(sc_), sigmaf(sigmaf_), sigmac(sigmac_), creg(creg_), rat(rat_), ga(ga_), t(BL_SPACEDIM)
    {
	for(int i = 0; i < BL_SPACEDIM; ++i)
	{
	    h[i] = h_[i];
	    t[i] = t_[i];
	}
    }
    virtual bool ready()
    {
	abort(); return false;
    };
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	abort(); return false;
    }
private:
    FCERES func;
    MultiFab& r;
    const MultiFab& s;
    const int igrid;
    task_fab* c;
    task_fab* sc;
    task_fab* sigmaf;
    task_fab* sigmac;
    const Box creg;
    Real h[BL_SPACEDIM];
    const IntVect rat;
    Array<int> ga;
    Array<int> t;
};

void holy_grail_amr_multigrid::alloc_sync_caches()
{
    if (lev_min < lev_max) 
    {
	fres_fbox  = new Box*[lev_max+1];
	fres_cbox  = new Box*[lev_max+1];
	fres_creg  = new Box*[lev_max+1];
	fres_sfbox = new Box*[lev_max+1];
	fres_scbox = new Box*[lev_max+1];
#if (BL_SPACEDIM == 3)
	eres_fbox  = new Box*[lev_max+1];
	eres_cbox  = new Box*[lev_max+1];
	eres_creg  = new Box*[lev_max+1];
	eres_sfbox = new Box*[lev_max+1];
	eres_scbox = new Box*[lev_max+1];
#endif
	cres_fbox  = new Box*[lev_max+1];
	cres_cbox  = new Box*[lev_max+1];
	cres_creg  = new Box*[lev_max+1];
	cres_sfbox = new Box*[lev_max+1];
	cres_scbox = new Box*[lev_max+1];
    }
    
    for (int lev = lev_min + 1; lev <= lev_max; lev++) 
    {
	int mglev = ml_index[lev];
	fres_fbox[lev]  = new Box[lev_interface[mglev].nboxes(level_interface::FACEDIM)];
	fres_cbox[lev]  = new Box[lev_interface[mglev].nboxes(level_interface::FACEDIM)];
	fres_creg[lev]  = new Box[lev_interface[mglev].nboxes(level_interface::FACEDIM)];
	fres_sfbox[lev] = new Box[lev_interface[mglev].nboxes(level_interface::FACEDIM)];
	fres_scbox[lev] = new Box[lev_interface[mglev].nboxes(level_interface::FACEDIM)];
#if (BL_SPACEDIM == 3)
	eres_fbox[lev]  = new Box[lev_interface[mglev].nboxes(1)];
	eres_cbox[lev]  = new Box[lev_interface[mglev].nboxes(1)];
	eres_creg[lev]  = new Box[lev_interface[mglev].nboxes(1)];
	eres_sfbox[lev] = new Box[lev_interface[mglev].nboxes(1)];
	eres_scbox[lev] = new Box[lev_interface[mglev].nboxes(1)];
#endif
	cres_fbox[lev]  = new Box[lev_interface[mglev].nboxes(0)];
	cres_cbox[lev]  = new Box[lev_interface[mglev].nboxes(0)];
	cres_creg[lev]  = new Box[lev_interface[mglev].nboxes(0)];
	cres_sfbox[lev] = new Box[lev_interface[mglev].nboxes(0)];
	cres_scbox[lev] = new Box[lev_interface[mglev].nboxes(0)];
	build_sync_cache(mglev, lev);
    }
}

void holy_grail_amr_multigrid::delete_sync_caches()
{
    for (int lev = lev_min + 1; lev <= lev_max; lev++) 
    {
	int mglev = ml_index[lev];
	delete [] fres_fbox[lev];
	delete [] fres_cbox[lev];
	delete [] fres_creg[lev];
	delete [] fres_sfbox[lev];
	delete [] fres_scbox[lev];
#if (BL_SPACEDIM == 3)
	delete [] eres_fbox[lev];
	delete [] eres_cbox[lev];
	delete [] eres_creg[lev];
	delete [] eres_sfbox[lev];
	delete [] eres_scbox[lev];
#endif
	delete [] cres_fbox[lev];
	delete [] cres_cbox[lev];
	delete [] cres_creg[lev];
	delete [] cres_sfbox[lev];
	delete [] cres_scbox[lev];
    }
    if (lev_min < lev_max) 
    {
	delete [] fres_fbox;
	delete [] fres_cbox;
	delete [] fres_creg;
	delete [] fres_sfbox;
	delete [] fres_scbox;
#if (BL_SPACEDIM == 3)
	delete [] eres_fbox;
	delete [] eres_cbox;
	delete [] eres_creg;
	delete [] eres_sfbox;
	delete [] eres_scbox;
#endif
	delete [] cres_fbox;
	delete [] cres_cbox;
	delete [] cres_creg;
	delete [] cres_sfbox;
	delete [] cres_scbox;
    }
}

void holy_grail_amr_multigrid::build_sync_cache(int mglev, int lev)
{
    const IntVect& rat = gen_ratio[lev-1];
    const int mglevc = ml_index[lev-1];
    
    const amr_boundary_class* bndry = (m_hg_terrain)? boundary.terrain_sigma() : boundary.scalar();
    
    for (int iface = 0; iface < lev_interface[mglev].nboxes(level_interface::FACEDIM); iface++) 
    {
	// find a fine grid touching this face
	int igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	if (igrid < 0)
	    igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	const unsigned int geo = lev_interface[mglev].geo(level_interface::FACEDIM, iface);
	// reject fine-fine interfaces and those without an interior fine grid
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) ) 
	{
	    continue;
	}
	// fine grid on just one side
	const int idim = lev_interface[mglev].fdim(iface);
	const int idir = (geo & level_interface::LOW) ? -1 : 1;
	Box& fbox = fres_fbox[lev][iface];
	Box& cbox = fres_cbox[lev][iface];
	Box& creg = fres_creg[lev][iface];
	fbox = dest[lev][igrid].box();
	cbox = lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	cbox.coarsen(rat);
	if (idir > 0)
	    cbox.growLo(idim, 1);
	else
	    cbox.growHi(idim, 1);
	Box& sigmafbox = fres_sfbox[lev][iface];
	Box& sigmacbox = fres_scbox[lev][iface];
	sigmafbox = sigma[mglev][igrid].box();
	sigmacbox = cbox;
	sigmacbox.convert(IntVect::TheCellVector());
	const IntVect t = lev_interface[mglev].box(level_interface::FACEDIM, iface).type();
	creg = lev_interface[mglev].node_box(level_interface::FACEDIM, iface);
	creg.coarsen(rat).grow(t - 1);
    }
    
#if (BL_SPACEDIM == 3)
    for (int iedge = 0; iedge < lev_interface[mglev].nboxes(1); iedge++) 
    {
	// find a fine grid touching this edge
	int igrid;
	for (int i = 0; i < lev_interface[mglev].ngrids(1); i++) 
	{
	    igrid = lev_interface[mglev].grid(1, iedge, i);
	    if (igrid >= 0)
		break;
	}
	const unsigned int geo = lev_interface[mglev].geo(1, iedge);
	// reject fine-fine interfaces and those without an interior fine grid
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(1, iedge) ) 
	{
	    continue;
	}
	Box& fbox = eres_fbox[lev][iedge];
	Box& cbox = eres_cbox[lev][iedge];
	Box& creg = eres_creg[lev][iedge];
	const IntVect t = lev_interface[mglev].box(1, iedge).type();
	cbox = lev_interface[mglev].node_box(1, iedge);
	cbox.coarsen(rat).grow(t);
	fbox = refine(cbox, rat);
	Box& sigmafbox = eres_sfbox[lev][iedge];
	Box& sigmacbox = eres_scbox[lev][iedge];
	sigmafbox = fbox;
	sigmafbox.convert(IntVect::TheCellVector());
	sigmacbox = cbox;
	sigmacbox.convert(IntVect::TheCellVector());
	creg = lev_interface[mglev].node_box(1, iedge);
	creg.coarsen(rat).grow(t - 1);
    }
    
#endif
    for (int icor = 0; icor < lev_interface[mglev].nboxes(0); icor++) 
    {
	// find a fine grid touching this corner
	int igrid;
	for (int i = 0; i < lev_interface[mglev].ngrids(0); i++) 
	{
	    igrid = lev_interface[mglev].grid(0, icor, i);
	    if (igrid >= 0)
		break;
	}
	const unsigned int geo = lev_interface[mglev].geo(0, icor);
	// reject fine-fine interfaces and those without an interior fine grid
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(0, icor) ) 
	{
	    continue;
	}
	Box& fbox = cres_fbox[lev][icor];
	Box& cbox = cres_cbox[lev][icor];
	Box& creg = cres_creg[lev][icor];
	cbox = lev_interface[mglev].box(0, icor);
	fbox = lev_interface[mglev].box(0, icor);
	cbox.coarsen(rat).grow(1);
	fbox.grow(rat);
	Box& sigmafbox = cres_sfbox[lev][icor];
	Box& sigmacbox = cres_scbox[lev][icor];
	sigmafbox = fbox;
	sigmafbox.convert(IntVect::TheCellVector());
	sigmacbox = cbox;
	sigmacbox.convert(IntVect::TheCellVector());
	creg = lev_interface[mglev].box(0, icor);
	creg.coarsen(rat);
    }
}


void holy_grail_amr_multigrid::interface_residual(int mglev, int lev)
{ 
    const amr_boundary_class* bndry = (m_hg_terrain)? boundary.terrain_sigma() : boundary.scalar();

    const IntVect& rat = gen_ratio[lev-1];
    const int mglevc = ml_index[lev-1];
    
    task_list tl;
    for (int iface = 0; iface < lev_interface[mglev].nboxes(level_interface::FACEDIM); iface++) 
    {
	// find a fine grid touching this face
	int igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 0);
	if (igrid < 0)
	    igrid = lev_interface[mglev].grid(level_interface::FACEDIM, iface, 1);
	const unsigned int geo = lev_interface[mglev].geo(level_interface::FACEDIM, iface);
	// reject fine-fine interfaces and those without an interior fine grid
	if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(level_interface::FACEDIM, iface) )
	    continue;
	// fine grid on just one side
	const int idim = lev_interface[mglev].fdim(iface);
	const int idir = (geo & level_interface::LOW) ? -1 : 1;
	const Box& cbox = fres_cbox[lev][iface];
	const Box& sigmacbox = fres_scbox[lev][iface];
	task_fab* sigmac = new task_fill_patch(sigmacbox, sigma[mglev], lev_interface[mglevc], bndry);
	const Box& creg = fres_creg[lev][iface];
	task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure());
	if (m_hg_terrain)
	{
	    tl.add_task(
		new task_fceres_2(&FORT_HGFRES_TERRAIN, resid[mglev], source[mglev], dest[mglev], sigma[mglev], igrid, cdst, sigmac, creg, h[mglev], rat, idim, idir)
		);
	}
	else if (m_hg_full_stencil)
	{
#if BL_SPACEDIM == 2
	    const int isRZ = IsRZ();
	    const int imax = mg_domain[mglevc].bigEnd(0) + 1;
	    tl.add_task(
		new task_fceres(&FORT_HGFRES_FULL, resid[mglev], source[mglev], dest[mglev], sigma[mglev], igrid, cdst, sigmac, creg, h[mglev], rat, idim, idir, isRZ, imax)
		);
#endif
	}
	else
	{
	    tl.add_task(
		new task_fceres_2(&FORT_HGFRES, resid[mglev], source[mglev], dest[mglev], sigma[mglev], igrid, cdst, sigmac, creg, h[mglev], rat, idim, idir)
		);
	}
    }
    tl.execute();

    if (m_hg_cross_stencil || m_hg_terrain)
    {
	
#if (BL_SPACEDIM == 3)
	
	task_list tl;
	for (int iedge = 0; iedge < lev_interface[mglev].nboxes(1); iedge++) 
	{
	    // find a fine grid touching this edge
	    int igrid;
	    for (int i = 0; i < lev_interface[mglev].ngrids(1); i++) 
	    {
		igrid = lev_interface[mglev].grid(1, iedge, i);
		if (igrid >= 0)
		    break;
	    }
	    const unsigned int geo = lev_interface[mglev].geo(1, iedge);
	    // reject fine-fine interfaces and those without an interior fine grid
	    if (geo != level_interface::ALL && igrid >= 0 && !lev_interface[mglev].flag(1, iedge) ) 
	    {
		const Box& fbox = eres_fbox[lev][iedge];
		const Box& cbox = eres_cbox[lev][iedge];
		const Box& sigmafbox = eres_sfbox[lev][iedge];
		const Box& sigmacbox = eres_scbox[lev][iedge];
		task_fab* sigmaf = new task_fill_patch(sigmafbox, sigma[mglev],  lev_interface[mglev],  bndry, 1, iedge);
		task_fab* sigmac = new task_fill_patch(sigmacbox, sigma[mglevc], lev_interface[mglevc], bndry);
		const Box& creg = eres_creg[lev][iedge];
		const IntVect t = lev_interface[mglev].box(1, iedge).type();
		task_fab* fdst = new task_fill_patch(fbox, dest[lev],   lev_interface[mglev],  boundary.pressure(), 1, iedge);
		task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure());
		Array<int> ga = lev_interface[mglev].geo_array(1, iedge);
		if (m_hg_terrain)
		{
		    tl.add_task(
			new task_fceres_4(&FORT_HGERES_TERRAIN, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, ga, t)
			);
		}
		else
		{
		    tl.add_task(
			new task_fceres_4(&FORT_HGERES, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, ga, t)
			);
		}
		// fill in the grids on the other sides, if any
		const Box& freg = lev_interface[mglev].node_box(1, iedge);
		for (int i = 1; i < lev_interface[mglev].ngrids(1); i++) 
		{
		    const int jgrid = lev_interface[mglev].grid(1, iedge, i);
		    if (jgrid >= 0 && jgrid != igrid)
			tl.add_task(new task_copy_link(resid[mglev], jgrid, igrid, freg, 0));
		}
	    }
	}
	tl.execute();
#endif
	for (int icor = 0; icor < lev_interface[mglev].nboxes(0); icor++) 
	{
	    // find a fine grid touching this corner
	    int igrid;
	    for (int i = 0; i < lev_interface[mglev].ngrids(0); i++) 
	    {
		igrid = lev_interface[mglev].grid(0, icor, i);
		if (igrid >= 0)
		    break;
	    }
	    const unsigned int geo = lev_interface[mglev].geo(0, icor);
	    // reject fine-fine interfaces and those without an interior fine grid
	    if (geo != level_interface::ALL && igrid >= 0 && !lev_interface[mglev].flag(0, icor) ) 
	    {
		const Box& fbox = cres_fbox[lev][icor];
		const Box& cbox = cres_cbox[lev][icor];
		const Box& sigmafbox = cres_sfbox[lev][icor];
		const Box& sigmacbox = cres_scbox[lev][icor];
		task_fab* sigmaf = new task_fill_patch(sigmafbox, sigma[mglev],  lev_interface[mglev],  bndry, 0, icor);
		task_fab* sigmac = new task_fill_patch(sigmacbox, sigma[mglevc], lev_interface[mglevc], bndry);
		const Box& creg = cres_creg[lev][icor];
		task_fab* fdst = new task_fill_patch(fbox, dest[lev],   lev_interface[mglev],  boundary.pressure(), 0, icor);
		task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure());
		Array<int> ga = lev_interface[mglev].geo_array(0, icor);
		if (m_hg_terrain)
		{
		    tl.add_task(
			new task_fceres_4(&FORT_HGCRES_TERRAIN, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, ga)
			);
		}
		else
		{
		    tl.add_task(
			new task_fceres_4(&FORT_HGCRES, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, ga)
			);
		}
		// fill in the grids on the other sides, if any
		const Box& freg = lev_interface[mglev].box(0, icor);
		for (int i = 1; i < lev_interface[mglev].ngrids(0); i++) 
		{
		    const int jgrid = lev_interface[mglev].grid(0, icor, i);
		    if (jgrid >= 0 && jgrid != igrid)
			tl.add_task(new task_copy_link(resid[mglev], jgrid, igrid, freg, 0));
		}
	    }
	}
	tl.execute();
    }
    else if (m_hg_full_stencil)
    {
#if BL_SPACEDIM == 2
	task_list tl;
	for (int icor = 0; icor < lev_interface[mglev].nboxes(0); icor++) 
	{
	    // find a fine grid touching this corner
	    int igrid;
	    for (int i = 0; i < lev_interface[mglev].ngrids(0); i++) 
	    {
		igrid = lev_interface[mglev].grid(0, icor, i);
		if (igrid >= 0)
		    break;
	    }
	    const unsigned int geo = lev_interface[mglev].geo(0, icor);
	    // reject fine-fine interfaces and those without an interior fine grid
	    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].flag(0, icor) )
		continue;
	    else if (geo == level_interface::XL || geo == level_interface::XH || geo == level_interface::YL || geo == level_interface::YH) 
	    {
		// fine grid on two adjacent sides
		const int idim = (geo == level_interface::XL || geo == level_interface::XH) ? 0 : 1;
		const int idir = (geo & level_interface::LL) ? -1 : 1;
		const Box& fbox = cres_fbox[lev][icor];
		const Box& cbox = cres_cbox[lev][icor];
		const Box& sigmafbox = cres_sfbox[lev][icor];
		const Box& sigmacbox = cres_scbox[lev][icor];
		task_fab* sigmaf = new task_fill_patch(sigmafbox, sigma[mglev],   lev_interface[mglev],  bndry, 0, icor);
		task_fab* sigmac = new task_fill_patch(sigmacbox, sigmac[mglevc], lev_interface[mglevc], bndry);
		const Box& creg = cres_creg[lev][icor];
		task_fab* fdst = new task_fill_patch(fbox, dest[lev],   lev_interface[mglev],  boundary.pressure(), 0, icor);
		task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure(), 0, -1);
		const int isRZ = IsRZ();
		const int imax = mg_domain[mglevc].bigEnd(0) + 1;
		tl.add_task(
		    new task_fceres_4(&FORT_HGFRES_FULL, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, idim, idir, isRZ, imax)
		    );
		// fill in the grids on the other sides, if any
		const Box& freg = lev_interface[mglev].box(0, icor);
		for (int i = 1; i < lev_interface[mglev].ngrids(0); i++) 
		{
		    const int jgrid = lev_interface[mglev].grid(0, icor, i);
		    if (jgrid >= 0 && jgrid != igrid)
			tl.add_task(new task_copy_link(resid[mglev], jgrid, igrid, freg, 0));
		}
	    }
	    else if (geo == level_interface::LL || geo == level_interface::HL || geo == level_interface::LH || geo == level_interface::HH) 
	    {
		// outside corner
		const int idir0 = (geo & level_interface::XL) ? -1 : 1;
		const int idir1 = (geo & level_interface::YL) ? -1 : 1;
		const Box& cbox = cres_cbox[lev][icor];
		const Box& sigmafbox = cres_sfbox[lev][icor];
		const Box& sigmacbox = cres_scbox[lev][icor];
		task_fab* sigmaf = new task_fill_patch(sigmafbox, sigma[mglev],  lev_interface[mglev],  bndry, 0, icor);
		task_fab* sigmac = new task_fill_patch(sigmacbox, sigma[mglevc], lev_interface[mglevc], bndry);
		const Box& creg = cres_creg[lev][icor];
		task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure());
		const Box& fbox = dest[lev][igrid].box();
		const int isRZ = IsRZ();
		//FIXME!!! might be a broken calling sequence.
		tl.add_task(new task_fceres_2(&FORT_HGORES, resid[mglev], source[lev], dest[lev], igrid, cdst, sigmaf, sigmac, creg, h[mglev], rat, idir0, idir1, isRZ, 0)
	    }
	    else if (geo == (level_interface::LL | level_interface::HH) || geo == (level_interface::LH | level_interface::HL)) 
	    {
		// diagonal corner
		const int jdir = (geo == (level_interface::LL | level_interface::HH)) ? 1 : -1;
		const Box& fbox = cres_fbox[lev][icor];
		const Box& cbox = cres_cbox[lev][icor];
		const Box& sigmafbox = cres_sfbox[lev][icor];
		const Box& sigmacbox = cres_scbox[lev][icor];
		task_fab* sigmaf = new task_fill_patch(sigmafbox, sigma[mglev],  lev_interface[mglev],  bndry, 0, icor);
		task_fab* sigmac = new task_fill_patch(sigmacbox, sigma[mglevc], lev_interface[mglevc], bndry);
		const Box& creg = cres_creg[lev][icor];
		task_fab* fdst = new task_fill_patch(fbox, dest[lev],   lev_interface[mglev],  boundary.pressure(), 0, icor);
		task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure(), 0, -1);
		const int isRZ = IsRZ();
		tl.add_task(
		    new task_fceres_4(&FORT_HGDRES, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, jdir,0, isRZ, 0)
		    );
		// fill in the grids on the other sides, if any
		const Box& freg = lev_interface[mglev].box(0, icor);
		for (int i = 1; i < lev_interface[mglev].ngrids(0); i++) 
		{
		    int jgrid = lev_interface[mglev].grid(0, icor, i);
		    if (jgrid >= 0 && jgrid != igrid)
			tl.add_task(new task_copy_link(resid[mglev], jgrid, igrid, freg, 0));
		}
	    }
	    else 
	    {
		// inside corner
		const int idir0 = ((geo & level_interface::XL) == level_interface::XL) ? -1 : 1;
		const int idir1 = ((geo & level_interface::YL) == level_interface::YL) ? -1 : 1;
		const Box& fbox = cres_fbox[lev][icor];
		const Box& cbox = cres_cbox[lev][icor];
		const Box& sigmafbox = cres_sfbox[lev][icor];
		const Box& sigmacbox = cres_scbox[lev][icor];
		task_fab* sigmaf = new task_fill_patch(sigmafbox, sigma[mglev],  lev_interface[mglev],  bndry, 0, icor);
		task_fab* sigmac = new task_fill_patch(sigmacbox, sigma[mglevc], lev_interface[mglevc], bndry);
		const Box& creg = cres_creg[lev][icor];
		task_fab* fdst = new task_fill_patch(fbox, dest[lev],   lev_interface[mglev],  boundary.pressure(), 0, icor);
		task_fab* cdst = new task_fill_patch(cbox, dest[lev-1], lev_interface[mglevc], boundary.pressure(), 0, -1);
		const int isRZ = IsRZ();
		tl.add_task(
		    new task_fceres_4(&FORT_HGIRES, resid[mglev], source[mglev], igrid, fdst, cdst, sigmaf, sigmac, creg, h[mglev], rat, idir0, idir1, isRZ, 0)
		    );
		// fill in the grids on the other sides, if any
		const Box& freg = lev_interface[mglev].box(0, icor);
		int kgrid = -1;
		for (int i = 1; i < lev_interface[mglev].ngrids(0); i++) 
		{
		    int jgrid = lev_interface[mglev].grid(0, icor, i);
		    if (jgrid >= 0 && jgrid != igrid && jgrid != kgrid) 
		    {
			tl.add_task(new task_copy_link(resid[mglev], jgrid, igrid, freg,0));
			kgrid = jgrid;
		    }
		}
	    }
        }
	tl.execute();
#endif
    }
}
