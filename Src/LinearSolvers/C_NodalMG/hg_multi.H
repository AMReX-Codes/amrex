#ifndef _HG_MULTI_H_
#define _HG_MULTI_H_

#include "amr_multi.H"

class holy_grail_interpolator_class_not_cross : public bilinear_interpolator_class 
{
public:
    holy_grail_interpolator_class_not_cross(Real **Sigptr, const Box& Sigbox)
	: sigptr(Sigptr), sigbox(Sigbox) { }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& cgr, const Box& cb, const IntVect& rat) const;
protected:
    Real **sigptr;
    const Box& sigbox;
};

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

class holy_grail_sigma_restrictor_class : public cell_average_restrictor_class 
{
public:
    holy_grail_sigma_restrictor_class() : cell_average_restrictor_class(0), m_hg_terrain(false) { }
    virtual void fill(FArrayBox& patch, const Box& region, const FArrayBox& fgr, const IntVect& rat) const;
private:
    bool m_hg_terrain;
};

//const holy_grail_sigma_restrictor_class holy_grail_sigma_restrictor_class;


class holy_grail_amr_multigrid : public amr_multigrid 
{
public:
    
    holy_grail_amr_multigrid(
	const Array<BoxArray>& Mesh, const Array<IntVect>& Gen_ratio,
	const Box& fdomain,
	int Lev_min_min, int Lev_min_max, int Lev_max_max,
	const amr_fluid_boundary_class& Boundary, 
	bool hg_terrain, bool hg_cross_stencil, bool hg_full_stencil,
	int Pcode = 2)
	: amr_multigrid(Mesh, Gen_ratio, Lev_min_min, Lev_min_max, Lev_max_max, Boundary.pressure(), Pcode),
	boundary(Boundary), smoother_mode(2), line_solve_dim(-1), 
	m_hg_terrain(hg_terrain),
	m_hg_cross_stencil(hg_cross_stencil),
	m_hg_full_stencil(hg_full_stencil)
    {
	assert( !(BL_SPACEDIM==3 && m_hg_full_stencil) );
	build_mesh(fdomain);
    }
    void set_line_solve_dimension(int dim)
    {
	if (dim != -1) BoxLib::Error("LineSolves not supported in parallel");
	line_solve_dim = dim;
    }
    void set_smoother_mode(int mode)
    {
	smoother_mode = mode;
    }
protected:
    void alloc(PArray<MultiFab>& Dest, PArray<MultiFab>& Source, PArray<MultiFab>& Coarse_source, PArray<MultiFab>& Sigma, 
	Real H[], int Lev_min, int Lev_max);
    void clear();
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
    virtual void level_residual(MultiFab& r, MultiFab& s, MultiFab& d, int mglev, bool iclear);
    virtual void relax(int mglev, int i1, bool is_zero);
    virtual void cgsolve(int mglev);
    
    void mg_restrict(int lto, int lfrom);
private:
    holy_grail_amr_multigrid(const holy_grail_amr_multigrid&);
    void operator=(const holy_grail_amr_multigrid&);
    
protected: 
    
    // these are mutually exclusive, should be a switch.
    // full stencil not implemented in 3D.
    bool m_hg_terrain;
    bool m_hg_cross_stencil;
    bool m_hg_full_stencil;
    int smoother_mode;
    int line_solve_dim;
    const amr_fluid_boundary_class& boundary;
    PArray<MultiFab> sigma;
    PArray<MultiFab> sigma_nd[BL_SPACEDIM];
    PArray<MultiFab> sigma_node;
    PArray<MultiFab> cen;
    Real (*h)[BL_SPACEDIM];
    
    PArray<MultiFab> mask;   // lives on mg_mesh
    PArray<MultiFab> cgwork; // all 8 components live on coarsest level
    bool singular;
    bool source_owned;
    
    // sync caches:
    Box** fres_fbox;
    Box** fres_cbox;
    Box** fres_sfbox;
    Box** fres_scbox;
    Box** fres_creg;
    PArray<FArrayBox>* fres_sc;
    PArray<FArrayBox>* fres_dc;
    bool** fres_flag;
#if (BL_SPACEDIM == 3)
    Box** eres_fbox;
    Box** eres_cbox;
    Box** eres_sfbox;
    Box** eres_scbox;
    Box** eres_creg;
    PArray<FArrayBox>* eres_sf;
    PArray<FArrayBox>* eres_sc;
    PArray<FArrayBox>* eres_df;
    PArray<FArrayBox>* eres_dc;
    bool** eres_flag;
#endif
    Box** cres_fbox;
    Box** cres_cbox;
    Box** cres_sfbox;
    Box** cres_scbox;
    Box** cres_creg;
    PArray<FArrayBox>* cres_sf;
    PArray<FArrayBox>* cres_sc;
    PArray<FArrayBox>* cres_df;
    PArray<FArrayBox>* cres_dc;
    bool** cres_flag;
    
    // orderings for full-level line solves:
    
    Array< Array<int> > line_order;
    Array< Array< List<int> > > line_after;
    
};

#endif
