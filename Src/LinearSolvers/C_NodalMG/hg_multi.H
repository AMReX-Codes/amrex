//BL_COPYRIGHT_NOTICE

#ifndef _HG_MULTI_H_
#define _HG_MULTI_H_

#include "amr_multi.H"

class holy_grail_amr_multigrid
    : public amr_multigrid 
{
public:
    enum stencil { terrain = 0, full = 1, cross = 2 };

    static bool is_dense(stencil sval)
    {
	return sval == terrain || sval == full;
    }

    holy_grail_amr_multigrid(const Array<BoxArray>& Mesh,
			     const Array<IntVect>& Gen_ratio,
			     const Box& fdomain,
			     int Lev_min_min,
			     int Lev_min_max,
			     int Lev_max_max,
			     const amr_fluid_boundary& Boundary,
			     stencil stencil_ = cross,
			     int Pcode = 2)
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
    void set_line_solve_dimension(int dim)
    {
	if (dim != -1) 
	{
	    BoxLib::Abort(
		"holy_grail_amr_multigrid::holy_grail_amr_multigrid():"
		"LineSolves not supported in parallel" );
	}
	line_solve_dim = dim;
    }
    void set_smoother_mode(int mode)
    {
	smoother_mode = mode;
    }

protected:
    void alloc(PArray<MultiFab>& Dest,
	       PArray<MultiFab>& Source,
	       PArray<MultiFab>& Coarse_source,
	       PArray<MultiFab>& Sigma,
	       Real H[],
	       int Lev_min,
	       int Lev_max,
	       int for_sync_reg);
    void sync_resid_clear();
    void clear();
    void build_line_order(int lsd);
    
    void build_sigma(PArray<MultiFab>& Sigma,
		     int for_sync_reg);
    
    void alloc_sync_caches();
    void delete_sync_caches();
    void build_sync_cache(int mglev,
			  int lev);
    
    // all these inherited in some form from amr_multigrid
    virtual bool can_coarsen(const BoxArray& mesh,
			     const Box& domain) const;
    virtual void sync_interfaces();
    
    // special case for periodic cleanup at the end:
    void sync_periodic_interfaces();
    
    virtual void mg_interpolate_level(int lto,
				      int lfrom);
    virtual void mg_restrict_level(int lto,
				   int lfrom);
    virtual void interface_residual(int mglev,
				    int lev);
    virtual void level_residual(MultiFab& r,
				MultiFab& s,
				MultiFab& d,
				int mglev,
				bool iclear,
				int for_sync_reg = 0);
    virtual void relax(int mglev,
		       int i1,
		       bool is_zero);
    virtual void cgsolve(int mglev);
    
    void mg_restrict(int lto,
		     int lfrom);
private:
    holy_grail_amr_multigrid(const holy_grail_amr_multigrid&);
    void operator=(const holy_grail_amr_multigrid&);
    
protected: 
    
    // these are mutually exclusive, should be a switch.
    // full stencil not implemented in 3D.
    stencil m_stencil;
    int smoother_mode;
    int line_solve_dim;
    const amr_fluid_boundary& boundary;
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
#if (BL_SPACEDIM == 3)
    Box** eres_fbox;
    Box** eres_cbox;
    Box** eres_sfbox;
    Box** eres_scbox;
    Box** eres_creg;
#endif
    Box** cres_fbox;
    Box** cres_cbox;
    Box** cres_sfbox;
    Box** cres_scbox;
    Box** cres_creg;
    
    // orderings for full-level line solves:
    
    Array< Array<int> > line_order;
    Array< Array< List<int> > > line_after;
    
};

#endif /*_HG_MULTI_H_*/
