#ifndef _AMR_MULTI_H_
#define _AMR_MULTI_H_

#include "fill_patch.H"

class amr_multigrid 
{
public:
    
    amr_multigrid(const Array<BoxArray>& Mesh, const Array<IntVect>& Gen_ratio,
	int Lev_min_min, int Lev_min_max, int Lev_max_max,
	const amr_boundary_class* Boundary, int Pcode) 
	: ml_mesh(Mesh), gen_ratio(Gen_ratio),
	lev_min_min(Lev_min_min), lev_min_max(Lev_min_max),
	lev_max_max(Lev_max_max), mg_boundary(Boundary), pcode(Pcode),
	integrate(0)
    {
    }
    
    virtual ~amr_multigrid();
    
#if BL_SPACEDIM==2
    static void SetRZ()		{ c_sys = 1; }
    static void SetCartesian()	{ c_sys = 0; }
    
    static int IsRZ()		{ return c_sys == 1; }
    static int IsCartesian()	{ return c_sys == 0; }
#endif
    
    static void mesh_read(Array<BoxArray>& m, Array<IntVect>& r, Array<Box>& d, istream& is);
    
    static void mesh_write(const Array<BoxArray>& m, const Array<Box>& d, ostream& os);
    static void mesh_write(const Array<BoxArray>& m, const Array<IntVect>& r, Box fd, ostream& os);
    
    void solve(Real reltol, Real abstol, int i1, int i2, int linesolvdim = -1);
    Real ml_cycle(int lev, int mglev, int i1, int i2, Real tol, Real res_norm_fine = 0.0);
    void mg_cycle(int mglev, int i1, int i2, bool is_zero);
    Real ml_residual(int mglev, int lev);

    virtual void mg_interpolate_level(int lto, int lfrom) = 0;	// defined
    virtual void mg_restrict_level(int lto, int lfrom) = 0;
    virtual void level_residual(MultiFab& r, MultiFab& s, MultiFab& d, int mglev, bool iclear) = 0;
    virtual void interface_residual(int mglev, int mgc) = 0;
    virtual void relax(int mglev, int i1, bool is_zero) = 0;
    virtual void cgsolve(int mglev) = 0;
protected:
    void build_index();
    
    int get_amr_level(int mglev) const 
    {
	for (int i = lev_min; i <= lev_max; i++) 
	{
	    if (ml_index[i] == mglev)
		return i;
	}
	return -1;
    }
    
    // All of these have built-in defaults except can_coarsen()
    void build_mesh(const Box& fdomain);
    int build_down(const BoxArray& l_mesh, const Box& l_domain, int flev, IntVect rat, int nlev);
    void make_coarser_level(BoxArray& mesh, Box& domain, int& flev, IntVect& rat);

    
    // Before using the multigrid, a derived class must call build_mesh()
    // followed by alloc(), or their equivalents.  The intended use is
    // that build_mesh() be called by the derived class constructor, while
    // alloc(), which actually creates large data objects, and clear(),
    // which deletes them, be called when the solver is actually being used.
    // These are made non-virtual so that a derived class can define a more
    // complicated lev_interface.
    void alloc(PArray<MultiFab>& Dest, PArray<MultiFab>& Source, PArray<MultiFab>& Coarse_source, int Lev_min, int Lev_max);
    void clear();
    
    virtual void sync_interfaces() = 0;
    virtual bool can_coarsen(const BoxArray& mesh, const Box& domain) const = 0;
    
    const int lev_min_min, lev_min_max, lev_max_max;
    int lev_min, lev_max, mglev_max;
    Array<IntVect> gen_ratio;            // assoc with ml_mesh
    Array<int> ml_index;                 // assoc with ml_mesh
    Array<Box> mg_domain;                // assoc with mg_mesh
    Array<BoxArray> ml_mesh, mg_mesh;
    Array< Array<Box> > mg_domain_array;
    Array< Array<BoxArray> > mg_mesh_array;
    const amr_boundary_class* mg_boundary;
    level_interface* lev_interface;
    level_interface** interface_array;
    PArray<MultiFab> source, dest;       // live on active levels of ml_mesh
    PArray<MultiFab> resid, corr, work;  // live on mg_mesh
    PArray<MultiFab> save;               // lives on middle levels of ml_mesh
    PArray<MultiFab> coarse_source;      // lives on ml_mesh, if it exists
    
#if BL_SPACEDIM==2
    static int c_sys;
#endif
    
    int pcode;
    int integrate;
};

#endif
