#ifndef _HG_MULTI_H_
#define _HG_MULTI_H_

#include "amr_multi.H"

#ifndef HG_CONSTANT

#ifndef HG_TERRAIN
#  if (BL_SPACEDIM == 3)
#    define HG_SIGMA_NODE
#  endif
#endif

#  ifndef HG_SIGMA_NODE

class holy_grail_interpolator_class : public bilinear_interpolator_class 
{
public:
    holy_grail_interpolator_class(Real **Sigptr, const Box& Sigbox)
	: sigptr(Sigptr), sigbox(Sigbox) { }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& cgr, const Box& cb, const IntVect& rat) const;
protected:
    Real **sigptr;
    const Box& sigbox;
};

#  else

class holy_grail_interpolator_class : public bilinear_interpolator_class 
{
public:
    holy_grail_interpolator_class(Real *Sigptr, const Box& Sigbox)
	: sigptr(Sigptr), sigbox(Sigbox) { }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& cgr, const Box& cb, const IntVect&) const;
protected:
    Real *sigptr;
    const Box& sigbox;
};

#  endif

class holy_grail_sigma_restrictor_class : public cell_average_restrictor_class 
{
public:
    holy_grail_sigma_restrictor_class() : cell_average_restrictor_class(0) { }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
};

//const holy_grail_sigma_restrictor_class holy_grail_sigma_restrictor_class;

#endif

class holy_grail_amr_multigrid : public amr_multigrid 
{
public:
    
    holy_grail_amr_multigrid(const Array<BoxArray>& Mesh, const Array<IntVect>& Gen_ratio,
	const Box& fdomain,
	int Lev_min_min, int Lev_min_max, int Lev_max_max,
	const amr_fluid_boundary_class& Boundary, int Pcode = 2)
	: amr_multigrid(Mesh, Gen_ratio, Lev_min_min, Lev_min_max, Lev_max_max,
	Boundary.pressure(), Pcode),
	boundary(Boundary), smoother_mode(2), line_solve_dim(-1)
    {
	build_mesh(fdomain);
    }
    
    void alloc(PArray<MultiFab>& Dest, PArray<MultiFab>& Source, PArray<MultiFab>& Coarse_source, PArray<MultiFab>& Sigma, 
	Real H[], int Lev_min, int Lev_max);
    void clear();
    
    int smoother_mode;
    int line_solve_dim;
protected:
    const amr_fluid_boundary_class& boundary;
    PArray<MultiFab> sigma, sigma_nd[BL_SPACEDIM], sigma_node, cen;
    Real (*h)[BL_SPACEDIM];
    
    PArray<MultiFab> mask;   // lives on mg_mesh
    PArray<MultiFab> cgwork; // all 8 components live on coarsest level
    bool singular;
    bool source_owned;
    
    // level caches which are not members of amr_multigrid:
#ifdef HG_USE_CACHE
    Array<copy_cache*> corr_scache;
    Array<unroll_cache*> cgw_ucache;
    copy_cache* cgw1_bcache;
#endif
    // sync caches:
    
    Box **fres_fbox, **fres_cbox, **fres_sfbox, **fres_scbox, **fres_creg;
    PArray<FArrayBox> *fres_sc, *fres_dc;
    int **fres_flag;
#if (BL_SPACEDIM == 3)
    Box **eres_fbox, **eres_cbox, **eres_sfbox, **eres_scbox, **eres_creg;
    PArray<FArrayBox> *eres_sf, *eres_sc, *eres_df, *eres_dc;
    int **eres_flag;
#endif
    Box **cres_fbox, **cres_cbox, **cres_sfbox, **cres_scbox, **cres_creg;
    PArray<FArrayBox> *cres_sf, *cres_sc, *cres_df, *cres_dc;
    int **cres_flag;
    
    // orderings for full-level line solves:
    
    Array< Array<int> > line_order;
    Array< Array< List<int> > > line_after;
    
    void build_line_order(int lsd);
    
    void build_sigma(PArray<MultiFab>& Sigma);
    
    void alloc_sync_caches();
    void delete_sync_caches();
    void build_sync_cache(int mglev, int lev);
    
    // all these inherited in some form from amr_multigrid
    virtual bool can_coarsen(const BoxArray& mesh, const Box& domain) const;
    virtual void sync_interfaces();
    
    // special case for periodic cleanup at the end:
    void sync_periodic_interfaces();
    
    virtual void mg_interpolate_level(int lto, int lfrom);
    virtual void mg_restrict_level(int lto, int lfrom);
    virtual void interface_residual(int mglev, int lev);
    virtual void level_residual(MultiFab& r,
		      MultiFab& s,
		      MultiFab& d,
#ifdef HG_USE_CACHE
		      copy_cache* dbc,
#endif
		      int mglev,
		      bool iclear);
    virtual void relax(int mglev, int i1, int is_zero);
    virtual void cgsolve(int mglev);
    
    virtual void mg_restrict(int lto, int lfrom);
    
};

#endif
