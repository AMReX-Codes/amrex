#ifndef _HG_MULTI_H_
#define _HG_MULTI_H_

#ifndef _AMR_MULTI_H_
#  include "amr_multi.H"
#endif

#include "hg_version.H"

#ifndef CONSTANT

#if (BL_SPACEDIM == 3)
#  define SIGMA_NODE
#endif

#  ifndef SIGMA_NODE

class holy_grail_interpolator_class : public bilinear_interpolator_class {
 protected:
  Real **sigptr;
  const Box& sigbox;
 public:
  holy_grail_interpolator_class(Real **Sigptr, const Box& Sigbox)
    : sigptr(Sigptr), sigbox(Sigbox) { }
  void fill(Fab& patch, const Box& region,
	    Fab& cgr, const Box& cb, const IntVect& rat) const;
};

#  else

class holy_grail_interpolator_class : public bilinear_interpolator_class {
 protected:
  Real *sigptr;
  const Box& sigbox;
 public:
  holy_grail_interpolator_class(Real *Sigptr, const Box& Sigbox)
    : sigptr(Sigptr), sigbox(Sigbox) { }
  void fill(Fab& patch, const Box& region,
	    Fab& cgr, const Box& cb, const IntVect&) const;
};

#  endif

class holy_grail_sigma_restrictor_class :
        public cell_average_restrictor_class {
 public:
  holy_grail_sigma_restrictor_class() { }
  void fill(Fab& patch, const Box& region,
	    Fab& fgr, const IntVect& rat) const;
};

const holy_grail_sigma_restrictor_class holy_grail_sigma_restrictor;

#endif

class holy_grail_amr_multigrid : public amr_multigrid {
 protected:
  amr_fluid_boundary boundary;
  PArray<MultiFab> sigma, sigma_nd[BL_SPACEDIM], sigma_node, cen;
  Real (*h)[BL_SPACEDIM];

  PArray<MultiFab> mask;   // lives on mg_mesh
  PArray<MultiFab> cgwork; // all 8 components live on coarsest level
  int singular, source_owned;

  // level caches which are not members of amr_multigrid:

  Array<copy_cache*> corr_scache;
  Array<unroll_cache*> cgw_ucache;
  copy_cache* cgw1_bcache;

  // sync caches:

  Box **fres_fbox, **fres_cbox, **fres_sfbox, **fres_scbox, **fres_creg;
  PArray<Fab> *fres_sc, *fres_dc;
  int **fres_flag;
#if (BL_SPACEDIM == 3)
  Box **eres_fbox, **eres_cbox, **eres_sfbox, **eres_scbox, **eres_creg;
  PArray<Fab> *eres_sf, *eres_sc, *eres_df, *eres_dc;
  int **eres_flag;
#endif
  Box **cres_fbox, **cres_cbox, **cres_sfbox, **cres_scbox, **cres_creg;
  PArray<Fab> *cres_sf, *cres_sc, *cres_df, *cres_dc;
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
  int can_coarsen(const BoxArray& mesh, const Box& domain);
  void sync_interfaces();

  // special case for periodic cleanup at the end:
  void sync_periodic_interfaces();

  void mg_interpolate_level(int lto, int lfrom);
  void mg_restrict_level(int lto, int lfrom);
  void interface_residual(int mglev, int lev);
  void level_residual(MultiFab& r,
		      MultiFab& s,
		      MultiFab& d,
		      copy_cache* dbc,
		      int mglev,
		      int iclear = 1);
  void relax(int mglev, int i1, int is_zero);
  void cgsolve(int mglev);

  void mg_restrict(int lto, int lfrom);

 public:

  holy_grail_amr_multigrid(Array<BoxArray>& Mesh, Array<IntVect>& Gen_ratio,
			   const Box& fdomain,
			   int Lev_min_min, int Lev_min_max, int Lev_max_max,
			   amr_fluid_boundary Boundary, int Pcode = 2)
    : amr_multigrid(Mesh, Gen_ratio, Lev_min_min, Lev_min_max, Lev_max_max,
		    Boundary.pressure(), Pcode),
      boundary(Boundary), smoother_mode(2), line_solve_dim(-1)
  {
    build_mesh(fdomain);
  }

  void alloc(PArray<MultiFab>& Dest, PArray<MultiFab>& Source,
	     PArray<MultiFab>& Coarse_source,
	     PArray<MultiFab>& Sigma, Real H[],
	     int Lev_min, int Lev_max);
  void clear();

  int smoother_mode;
  int line_solve_dim;
};

#endif
