
#include "hg_projector.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGDIV      hgdiv_
#  define   FORT_HGFDIV     hgfdiv_
#  define   FORT_HGEDIV     hgediv_
#  define   FORT_HGCDIV     hgcdiv_
#  define   FORT_HGODIV     hgodiv_
#  define   FORT_HGIDIV     hgidiv_
#  define   FORT_HGDDIV     hgddiv_
#  define   FORT_HGGRAD     hggrad_
#  define   FORT_HGGRAD     hggrad_
#  define   FORT_HGAVG      hgavg_
#  define   FORT_HGFAVG     hgfavg_
#  define   FORT_HGEAVG     hgeavg_
#  define   FORT_HGCAVG     hgcavg_
#else
#  define   FORT_HGDIV      HGDIV
#  define   FORT_HGFDIV     HGFDIV
#  define   FORT_HGEDIV     HGEDIV
#  define   FORT_HGCDIV     HGCDIV
#  define   FORT_HGODIV     HGODIV
#  define   FORT_HGIDIV     HGIDIV
#  define   FORT_HGDDIV     HGDDIV
#  define   FORT_HGGRAD     HGGRAD
#  define   FORT_HGGRAD     HGGRAD
#  define   FORT_HGAVG      HGAVG
#  define   FORT_HGFAVG     HGFAVG
#  define   FORT_HGEAVG     HGEAVG
#  define   FORT_HGCAVG     HGCAVG
#endif

extern "C" {

#if (BL_SPACEDIM == 1)
  ERROR, not relevant
#elif (BL_SPACEDIM == 2 || BL_SPACEDIM == 3)
#  if (BL_SPACEDIM == 2)
  void FORT_HGGRAD(RealPS, intS, Real*, intS, intS, RealRS, const int&);
  void FORT_HGDIV(Real*, intS, RealPS, intS, intS, RealRS,
		  const int&, const int&);
  void FORT_HGFDIV(Real*, intS, RealPS, intS, RealPS, intS,
		   intS, RealRS, intRS, int&, int&, const int&, const int&);
  void FORT_HGODIV(Real*, intS, Real*, Real*, intS, Real*, Real*, intS,
		   intS, Real&, Real&, intRS, int&, int&, const int&);
  void FORT_HGIDIV(Real*, intS, Real*, Real*, intS, Real*, Real*, intS,
		   intS, Real&, Real&, intRS, int&, int&, const int&);
  void FORT_HGDDIV(Real*, intS, Real*, Real*, intS, Real*, Real*, intS,
		   intS, Real&, Real&, intRS, int&, const int&);
  void FORT_HGAVG(Real*, intS, Real*, intS, intS,
                  const Real&, const int&, const int&);
  void FORT_HGFAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, intRS, int&, int&,
                   const Real&, const int&, const int&);
  void FORT_HGCAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, intRS, const int*,
                   const Real&, const int&, const int&);
#  elif (BL_SPACEDIM == 3)
  void FORT_HGGRAD(RealPS, intS, Real*, intS, intS, RealRS);
  void FORT_HGDIV(Real*, intS, RealPS, intS, intS, RealRS);
  void FORT_HGFDIV(Real*, intS, RealPS, intS, RealPS, intS,
		   intS, RealRS, intRS, int&, int&);
  void FORT_HGEDIV(Real*, intS, RealPS, intS, RealPS, intS,
		   intS, RealRS, intRS, const int*, const int*);
  void FORT_HGCDIV(Real*, intS, RealPS, intS, RealPS, intS,
		   intS, RealRS, intRS, const int*);
  void FORT_HGAVG(Real*, intS, Real*, intS, intS);
  void FORT_HGFAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, intRS, int&, int&);
  void FORT_HGEAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, intRS, const int*, const int*);
  void FORT_HGCAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, intRS, const int*);
#  endif
#endif
}

PArray<MultiFab> null_amr_real;

void holy_grail_amr_projector::project(PArray<MultiFab>* u,
				       PArray<MultiFab>& p,
				       PArray<MultiFab>& Coarse_source,
				       PArray<MultiFab>& Sigma,
				       Real H[], Real tol,
				       int Lev_min, int Lev_max,
				       Real scale)
{
  if (Lev_min < 0)
    Lev_min = lev_min_max;
  if (Lev_max < 0)
    Lev_max = Lev_min;

#ifndef CONSTANT
  assert(Sigma.length() > 0);
#endif

  assert(u[      0      ][Lev_min].nGrow() == 1);
  assert(u[      1      ][Lev_min].nGrow() == 1);
  assert(u[BL_SPACEDIM-1][Lev_min].nGrow() == 1);

  alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
  right_hand_side(u, null_amr_real);
  if (singular && Coarse_source.ready() && make_sparse_node_source_solvable) {
    sparse_node_source_adjustment(coarse_source);
  }

  solve(tol, scale, 2, 2);
  form_solution_vector(u, Sigma);
  clear();
}

void holy_grail_amr_projector::sync_project(PArray<MultiFab>* u,
					    PArray<MultiFab>& p,
					    PArray<MultiFab>& Coarse_source,
					    PArray<MultiFab>& Sigma,
					    Real H[], Real tol,
					    int Lev_min, int Lev_max,
					    Real scale)
{
  if (Lev_min < 0)
    Lev_min = lev_min_max;
  if (Lev_max < 0)
    Lev_max = Lev_min;

#ifndef CONSTANT
  assert(Sigma.length() > 0);
#endif

  assert(u[      0      ][Lev_min].nGrow() == 1);
  assert(u[      1      ][Lev_min].nGrow() == 1);
  assert(u[BL_SPACEDIM-1][Lev_min].nGrow() == 1);

  alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
  sync_right_hand_side(u);
  if (singular && Coarse_source.ready() && make_sparse_node_source_solvable) {
    sparse_node_source_adjustment(coarse_source);
  }

  solve(tol, scale, 2, 2);
  form_solution_vector(u, Sigma);
  clear();
}

void holy_grail_amr_projector::manual_project(PArray<MultiFab>* u,
					      PArray<MultiFab>& p,
					      PArray<MultiFab>& rhs,
					      PArray<MultiFab>& Coarse_source,
					      PArray<MultiFab>& Sigma,
					      int use_u,
					      Real H[], Real tol,
					      int Lev_min, int Lev_max,
					      Real scale)
{
  if (Lev_min < 0)
    Lev_min = lev_min_max;
  if (Lev_max < 0)
    Lev_max = Lev_min;

#ifndef CONSTANT
  assert(Sigma.length() > 0);
#endif

  if (use_u) {
    assert(u[      0      ][Lev_min].nGrow() == 1);
    assert(u[      1      ][Lev_min].nGrow() == 1);
    assert(u[BL_SPACEDIM-1][Lev_min].nGrow() == 1);

    alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
    if (rhs.length() > 0) {
      if (type(rhs[Lev_min]) == nodevect) {
	right_hand_side(u, null_amr_real);
	for (int lev = Lev_min; lev <= Lev_max; lev++)
	  source[lev].plus(rhs[lev], 0, 1, 0);
	if (singular && make_sparse_node_source_solvable) {
	  sparse_node_source_adjustment(rhs);
	}
      }
      else {
	assert(rhs[Lev_min].nGrow() == 1);
	right_hand_side(u, rhs);
      }
    }
    else {
      right_hand_side(u, null_amr_real);
    }
  }
  else {
    assert(rhs.length() > 0);
    assert(rhs[Lev_min].nGrow() == 1);

    if (type(rhs[Lev_min]) == nodevect) {
      alloc(p, rhs, Coarse_source, Sigma, H, Lev_min, Lev_max);
      if (singular && make_sparse_node_source_solvable) {
	sparse_node_source_adjustment(rhs);
      }
    }
    else {
      alloc(p, null_amr_real, Coarse_source, Sigma, H, Lev_min, Lev_max);
      // source is set to 0 at this point
      right_hand_side((PArray<MultiFab>*)0, rhs);
    }
  }
  if (singular && Coarse_source.ready() && make_sparse_node_source_solvable) {
    sparse_node_source_adjustment(coarse_source);
  }

  solve(tol, scale, 2, 2);
  form_solution_vector(u, Sigma);
  clear();
}

void holy_grail_amr_projector::sparse_node_source_adjustment(PArray<MultiFab>&
							     sparse_source)
{
  // This routine takes advantage of the sparse structure of
  // the sync source to avoid costly restriction operations.  It
  // is necessary to use the inner_product routine, which weights
  // boundary nodes, since the coarse-fine interface can touch
  // the boundary.  Otherwise a call to sum would suffice.

  // Note that the correction is applied to source, not to
  // sparse_source, since the sparse structure of the latter
  // may need to be preserved.

  assert(singular);
  assert(make_sparse_node_source_solvable);

  int lev, i;
  Real adjust = 0.0;
  for (lev = lev_max; lev >= lev_min; lev--) {
    if (sparse_source.defined(lev)) {
      int mglev = ml_index[lev];
      corr[mglev].setVal(1.0);
      adjust += inner_product(sparse_source[lev], corr[mglev]);
    }
    if (lev > lev_min && adjust != 0.0) {
      for (i = 0; i < BL_SPACEDIM; i++) {
	adjust /= gen_ratio[lev-1][i];
      }
    }
  }
  if (adjust != 0.0) {
    adjust /= mg_domain[ml_index[lev_min]].numPts();

    if (pcode >= 2)
      cout << "Sparse-source solvability adjustment: " << adjust << endl;

    for (lev = lev_min; lev <= lev_max; lev++) {
      source[lev].plus(-adjust, 0);
    }
  }
}

// This is a combination routine which combines sources from a divergence
// and from a cell-based right hand side S in the proper sequence.  The
// key feature is that both "grid" routines must be called before starting
// the interface calculation, since they trash some interface points.

void holy_grail_amr_projector::right_hand_side(PArray<MultiFab>* u,
					       PArray<MultiFab>& S)
{
  if (u) {
    grid_divergence(u);
  }
  if (S.length() > 0) {
    grid_average(S);
  }

  int lev;
  for (lev = lev_min; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];

    clear_part_interface(source[lev], interface[mglev]);

    if (lev > lev_min) {
      if (u) {
	interface_divergence(u, lev);
      }
      if (S.length() > 0) {
	interface_average(S, lev);
      }
    }
  }
}

// Averages cell-based S onto node-based source conservatively
// across the composite mesh.  S must be passed in with a one
// cell wide border.  At inflow boundaries border values should
// be set to correct inflow condition.  Other border values passed
// in may be meaningless, but should not be NaNs.

// This routine will modify the borders of S.  Also, if the problem
// being solved is singular, S will be adjusted so that it integrates
// to 0 to maximum precision.  (It is assumed that any additional
// contribution to the right hand side will also integrate to 0.)

// This is an incomplete routine---interface_average must also be called.

void holy_grail_amr_projector::grid_average(PArray<MultiFab>& S)
{
  assert(S[lev_min].nGrow() == 1);

  int lev, igrid;
  if (singular) {
    Real adjust = 0.0;
    for (lev = lev_max; lev > lev_min; lev--) {
      restrict_level(S[lev-1], 0, S[lev], gen_ratio[lev-1]);
    }
    for (igrid = 0; igrid < ml_mesh[lev_min].length(); igrid++) {
      adjust += S[lev_min][igrid].sum(S[lev_min].box(igrid), 0);
    }
    adjust /= mg_domain[ml_index[lev_min]].numPts();

    if (pcode >= 2)
      cout << "Cell-source solvability adjustment: " << adjust << endl;

    for (lev = lev_min; lev <= lev_max; lev++) {
      S[lev].plus(-adjust, 0);
    }
  }

  for (lev = lev_min; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    Real hx = h[mglev][0];

    fill_borders(S[lev], 0, interface[mglev], boundary.scalar());

    for (igrid = 0; igrid < ml_mesh[lev].length(); igrid++) {
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = S[lev][igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
      Real *const sptr = source[lev][igrid].dataPtr();
      Real *const csptr = S[lev][igrid].dataPtr();
#if (BL_SPACEDIM == 2)
      FORT_HGAVG(sptr, dimlist(sbox),
		 csptr, dimlist(fbox), dimlist(freg),
		 hx, IsRZ(), mg_domain[mglev].bigEnd(0) + 1);
#else
      FORT_HGAVG(sptr, dimlist(sbox),
		 csptr, dimlist(fbox), dimlist(freg));
#endif
    }
  }
}

// This is an incomplete routine---interface_divergence must also be called.

void holy_grail_amr_projector::grid_divergence(PArray<MultiFab>* u)
{ 
  int lev;
  for (lev = lev_min; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    Real hx = h[mglev][0];
    Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
    Real hz = h[mglev][2];
#endif
    int igrid, i;

    for (i = 0; i < BL_SPACEDIM; i++) {
      fill_borders(u[i][lev], 0, interface[mglev], boundary.velocity(i));
    }

    for (igrid = 0; igrid < ml_mesh[lev].length(); igrid++) {
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = u[0][lev][igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
      Real *const sptr = source[lev][igrid].dataPtr();
      Real *const u0ptr = u[0][lev][igrid].dataPtr();
      Real *const u1ptr = u[1][lev][igrid].dataPtr();
#if (BL_SPACEDIM == 2)
      FORT_HGDIV(sptr, dimlist(sbox),
		 u0ptr, u1ptr, dimlist(fbox), dimlist(freg), hx, hy,
		 IsRZ(), mg_domain[mglev].bigEnd(0) + 1);
#else
      Real *const u2ptr = u[2][lev][igrid].dataPtr();
      FORT_HGDIV(sptr, dimlist(sbox),
		 u0ptr, u1ptr, u2ptr, dimlist(fbox), dimlist(freg),
		 hx, hy, hz);
#endif
    }
  }
}

// Obsolete:
void holy_grail_amr_projector::sync_right_hand_side(PArray<MultiFab>* u)
{ 
  int lev;
  for (lev = lev_min; lev <= lev_max; lev++) {
    source[lev].setVal(0.0);
  }

  int mglev0 = ml_index[lev_min];
  interface_divergence(u, lev_min+1);

  if (singular) {
    int mglev1 = ml_index[lev_min+1];
    restrict_level(source[lev_min], 0, source[lev_min+1],
		   gen_ratio[lev_min], 0,
		   bilinear_restrictor_coarse,
		   interface[mglev1], mg_boundary);
    work[mglev0].setVal(1.0);
    Real adjustment = inner_product(source[lev_min], work[mglev0]) /
      mg_domain[ml_index[lev_min]].volume();
    if (pcode >= 2)
      cout << "Solvability adjustment is " << adjustment << endl;
    for (lev = lev_min; lev <= lev_max; lev++) {
      source[lev].plus(-adjustment, 0);
    }
  }
}

void holy_grail_amr_projector::interface_average(PArray<MultiFab>& S, int lev)
{ 
  int mglev = ml_index[lev];
  int mgc = ml_index[lev-1];
  Real hx = h[mglev][0];
  int igrid, jgrid, i;

  DECLARE_GEOMETRY_TYPES;

  const IntVect& rat = gen_ratio[lev-1];
  for (int iface = 0; iface < interface[mglev].nfaces(); iface++) {
    // find a fine grid touching this face
    igrid = interface[mglev].fgrid(iface, 0);
    if (igrid < 0)
      igrid = interface[mglev].fgrid(iface, 1);
    unsigned geo = interface[mglev].fgeo(iface);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].fflag(iface) == 1)
      continue;
    // fine grid on just one side
    int idim = interface[mglev].fdim(iface);
    int idir = (geo & LOW) ? -1 : 1;
    const Box& sbox = source[lev][igrid].box();
    const Box& fbox = S[lev][igrid].box();
    Box cbox = interface[mglev].face(iface);
    IntVect t = cbox.type();
    if (idir > 0)
      cbox.growLo(idim, rat[idim]);
    else
      cbox.growHi(idim, rat[idim]);
    cbox.convert(cellvect).coarsen(rat);
    Fab *Scp;
    jgrid = find_patch(cbox, S[lev-1]);
    if (jgrid < 0) {
      Scp = new Fab(cbox);
      fill_patch(*Scp, S[lev-1], interface[mgc], boundary.scalar());
    }
    else {
      Scp = &S[lev-1][jgrid];
      cbox = Scp->box();
    }
    Fab& Sc = *Scp;
    Box creg = interface[mglev].node_face(iface);
    creg.coarsen(rat).grow(t - unitvect);
    Real *const sptr = source[lev][igrid].dataPtr();
    Real *const Sfptr = S[lev][igrid].dataPtr();
#if (BL_SPACEDIM == 2)
    FORT_HGFAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sfptr, dimlist(fbox), dimlist(creg),
		rat[0], rat[1], idim, idir,
		hx, IsRZ(), mg_domain[mgc].bigEnd(0) + 1);
#else
    FORT_HGFAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sfptr, dimlist(fbox), dimlist(creg),
		rat[0], rat[1], rat[2], idim, idir);
#endif
    if (jgrid < 0) {
      delete Scp;
    }
  }

  int ga[N_CORNER_GRIDS];

#if (BL_SPACEDIM == 3)

  for (int iedge = 0; iedge < interface[mglev].nedges(); iedge++) {
    // find a fine grid touching this edge
    for (i = 0; i < N_EDGE_GRIDS; i++) {
      igrid = interface[mglev].egrid(iedge, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = interface[mglev].egeo(iedge);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].eflag(iedge) == 1)
      continue;
    // fine grid on just one side
    const Box& sbox = source[lev][igrid].box();
    Real *const sptr = source[lev][igrid].dataPtr();
    Box cbox = interface[mglev].edge(iedge);
    IntVect t = cbox.type();
    cbox.coarsen(rat).grow(t).convert(cellvect);
    Box fbox = cbox;
    fbox.refine(rat);
    Fab *Scp;
    jgrid = find_patch(cbox, S[lev-1]);
    if (jgrid < 0) {
      Scp = new Fab(cbox);
      fill_patch(*Scp, S[lev-1], interface[mgc], boundary.scalar());
    }
    else {
      Scp = &S[lev-1][jgrid];
      cbox = Scp->box();
    }
    Fab& Sc = *Scp;
    Fab Sf(fbox);
    fill_patch(Sf, S[lev], interface[mglev], boundary.scalar(),
	       0, 1, iedge);
    Box creg = interface[mglev].node_edge(iedge);
    creg.coarsen(rat).grow(t - unitvect);
    interface[mglev].geo_array(ga, 1, iedge);
    FORT_HGEAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sf.dataPtr(), dimlist(fbox),
		dimlist(creg), rat[0], rat[1], rat[2], t.getVect(), ga);
    if (jgrid < 0) {
      delete Scp;
    }
    // fill in the grids on the other sides, if any
    const Box& freg = interface[mglev].node_edge(iedge);
    for (i = 1; i < N_EDGE_GRIDS; i++) {
      jgrid = interface[mglev].egrid(iedge, i);
      if (jgrid >= 0 && jgrid != igrid)
	internal_copy(source[lev], jgrid, igrid, freg);
    }
  }

#endif

  for (int icor = 0; icor < interface[mglev].ncorners(); icor++) {
    // find a fine grid touching this corner
    for (i = 0; i < N_CORNER_GRIDS; i++) {
      igrid = interface[mglev].cgrid(icor, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = interface[mglev].cgeo(icor);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].cflag(icor) == 1)
      continue;
    // fine grid on just one side
    const Box& sbox = source[lev][igrid].box();
    Real *const sptr = source[lev][igrid].dataPtr();
    Box cbox = interface[mglev].corner(icor);
    cbox.coarsen(rat).grow(1).convert(cellvect);
    Box fbox = cbox;
    fbox.refine(rat);
    Fab *Scp;
    jgrid = find_patch(cbox, S[lev-1]);
    if (jgrid < 0) {
      Scp = new Fab(cbox);
      fill_patch(*Scp, S[lev-1], interface[mgc], boundary.scalar());
    }
    else {
      Scp = &S[lev-1][jgrid];
      cbox = Scp->box();
    }
    Fab& Sc = *Scp;
    Fab Sf(fbox);
    fill_patch(Sf, S[lev], interface[mglev], boundary.scalar(),
	       0, 0, icor);
    Box creg = interface[mglev].corner(icor);
    creg.coarsen(rat);
    interface[mglev].geo_array(ga, 0, icor);
#if (BL_SPACEDIM == 2)
    FORT_HGCAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sf.dataPtr(), dimlist(fbox),
		dimlist(creg), rat[0], rat[1], ga,
		hx, IsRZ(), mg_domain[mgc].bigEnd(0) + 1);
#else
    FORT_HGCAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sf.dataPtr(), dimlist(fbox),
		dimlist(creg), rat[0], rat[1], rat[2], ga);
#endif
    if (jgrid < 0) {
      delete Scp;
    }
    // fill in the grids on the other sides, if any
    const Box& freg = interface[mglev].corner(icor);
    for (i = 1; i < N_CORNER_GRIDS; i++) {
      jgrid = interface[mglev].cgrid(icor, i);
      if (jgrid >= 0 && jgrid != igrid)
	internal_copy(source[lev], jgrid, igrid, freg);
    }
  }
}

void holy_grail_amr_projector::interface_divergence(PArray<MultiFab>* u,
						    int lev)
{ 
  int mglev = ml_index[lev];
  int mgc = ml_index[lev-1];
  Real hx = h[mglev][0];
  Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
  Real hz = h[mglev][2];
#endif
  int igrid, jgrid, i;

  DECLARE_GEOMETRY_TYPES;

  const IntVect& rat = gen_ratio[lev-1];
  for (int iface = 0; iface < interface[mglev].nfaces(); iface++) {
    // find a fine grid touching this face
    igrid = interface[mglev].fgrid(iface, 0);
    if (igrid < 0)
      igrid = interface[mglev].fgrid(iface, 1);
    unsigned geo = interface[mglev].fgeo(iface);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].fflag(iface) == 1)
      continue;
    // fine grid on just one side
    int idim = interface[mglev].fdim(iface);
    int idir = (geo & LOW) ? -1 : 1;
    const Box& sbox = source[lev][igrid].box();
    const Box& fbox = u[0][lev][igrid].box();
    Box cbox = interface[mglev].face(iface);
    IntVect t = cbox.type();
    if (idir > 0)
      cbox.growLo(idim, rat[idim]);
    else
      cbox.growHi(idim, rat[idim]);
    cbox.convert(cellvect).coarsen(rat);
    Fab *ucp, *vcp;
    jgrid = find_patch(cbox, u[0][lev-1]);
    if (jgrid < 0) {
      ucp = new Fab(cbox);
      vcp = new Fab(cbox);
      fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
      fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
    }
    else {
      ucp = &u[0][lev-1][jgrid];
      vcp = &u[1][lev-1][jgrid];
      cbox = ucp->box();
    }
    Fab& uc = *ucp;
    Fab& vc = *vcp;
#if (BL_SPACEDIM == 3)
    Fab *wcp;
    if (jgrid < 0) {
      wcp = new Fab(cbox);
      fill_patch(*wcp, u[2][lev-1], interface[mgc], boundary.velocity(2));
    }
    else {
      wcp = &u[2][lev-1][jgrid];
    }
    Fab& wc = *wcp;
#endif
    Box creg = interface[mglev].node_face(iface);
    creg.coarsen(rat).grow(t - unitvect);
    Real *const sptr = source[lev][igrid].dataPtr();
    Real *const u0ptr = u[0][lev][igrid].dataPtr();
    Real *const u1ptr = u[1][lev][igrid].dataPtr();
#if (BL_SPACEDIM == 2)
    FORT_HGFDIV(sptr, dimlist(sbox),
		uc.dataPtr(), vc.dataPtr(), dimlist(cbox),
		u0ptr, u1ptr, dimlist(fbox), dimlist(creg),
		hx, hy, rat[0], rat[1], idim, idir,
		IsRZ(), mg_domain[mgc].bigEnd(0) + 1);
#else
    Real *const u2ptr = u[2][lev][igrid].dataPtr();
    FORT_HGFDIV(sptr, dimlist(sbox),
		uc.dataPtr(), vc.dataPtr(), wc.dataPtr(), dimlist(cbox),
		u0ptr, u1ptr, u2ptr, dimlist(fbox), dimlist(creg),
		hx, hy, hz, rat[0], rat[1], rat[2], idim, idir);
    if (jgrid < 0) {
      delete wcp;
    }
#endif
    if (jgrid < 0) {
      delete ucp;
      delete vcp;
    }
  }

#if (BL_SPACEDIM == 3)

  int ga[N_CORNER_GRIDS];

  for (int iedge = 0; iedge < interface[mglev].nedges(); iedge++) {
    // find a fine grid touching this edge
    for (i = 0; i < N_EDGE_GRIDS; i++) {
      igrid = interface[mglev].egrid(iedge, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = interface[mglev].egeo(iedge);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].eflag(iedge) == 1)
      continue;
    // fine grid on just one side
    const Box& sbox = source[lev][igrid].box();
    Real *const sptr = source[lev][igrid].dataPtr();
    Box cbox = interface[mglev].edge(iedge);
    IntVect t = cbox.type();
    cbox.coarsen(rat).grow(t).convert(cellvect);
    Box fbox = cbox;
    fbox.refine(rat);
    Fab *ucp, *vcp, *wcp;
    jgrid = find_patch(cbox, u[0][lev-1]);
    if (jgrid < 0) {
      ucp = new Fab(cbox);
      vcp = new Fab(cbox);
      wcp = new Fab(cbox);
      fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
      fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
      fill_patch(*wcp, u[2][lev-1], interface[mgc], boundary.velocity(2));
    }
    else {
      ucp = &u[0][lev-1][jgrid];
      vcp = &u[1][lev-1][jgrid];
      wcp = &u[2][lev-1][jgrid];
      cbox = ucp->box();
    }
    Fab& uc = *ucp;
    Fab& vc = *vcp;
    Fab& wc = *wcp;
    Fab uf(fbox), vf(fbox), wf(fbox);
    fill_patch(uf, u[0][lev], interface[mglev], boundary.velocity(0),
	       0, 1, iedge);
    fill_patch(vf, u[1][lev], interface[mglev], boundary.velocity(1),
	       0, 1, iedge);
    fill_patch(wf, u[2][lev], interface[mglev], boundary.velocity(2),
	       0, 1, iedge);
    Box creg = interface[mglev].node_edge(iedge);
    creg.coarsen(rat).grow(t - unitvect);
    interface[mglev].geo_array(ga, 1, iedge);
    FORT_HGEDIV(sptr, dimlist(sbox),
		uc.dataPtr(), vc.dataPtr(), wc.dataPtr(), dimlist(cbox),
		uf.dataPtr(), vf.dataPtr(), wf.dataPtr(), dimlist(fbox),
		dimlist(creg), hx, hy, hz, rat[0], rat[1], rat[2],
		t.getVect(), ga);
    if (jgrid < 0) {
      delete ucp;
      delete vcp;
      delete wcp;
    }
    // fill in the grids on the other sides, if any
    const Box& freg = interface[mglev].node_edge(iedge);
    for (i = 1; i < N_EDGE_GRIDS; i++) {
      jgrid = interface[mglev].egrid(iedge, i);
      if (jgrid >= 0 && jgrid != igrid)
	internal_copy(source[lev], jgrid, igrid, freg);
    }
  }

  for (int icor = 0; icor < interface[mglev].ncorners(); icor++) {
    // find a fine grid touching this corner
    for (i = 0; i < N_CORNER_GRIDS; i++) {
      igrid = interface[mglev].cgrid(icor, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = interface[mglev].cgeo(icor);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].cflag(icor) == 1)
      continue;
    // fine grid on just one side
    const Box& sbox = source[lev][igrid].box();
    Real *const sptr = source[lev][igrid].dataPtr();
    Box cbox = interface[mglev].corner(icor);
    cbox.coarsen(rat).grow(1).convert(cellvect);
    Box fbox = cbox;
    fbox.refine(rat);
    Fab *ucp, *vcp, *wcp;
    jgrid = find_patch(cbox, u[0][lev-1]);
    if (jgrid < 0) {
      ucp = new Fab(cbox);
      vcp = new Fab(cbox);
      wcp = new Fab(cbox);
      fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
      fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
      fill_patch(*wcp, u[2][lev-1], interface[mgc], boundary.velocity(2));
    }
    else {
      ucp = &u[0][lev-1][jgrid];
      vcp = &u[1][lev-1][jgrid];
      wcp = &u[2][lev-1][jgrid];
      cbox = ucp->box();
    }
    Fab& uc = *ucp;
    Fab& vc = *vcp;
    Fab& wc = *wcp;
    Fab uf(fbox), vf(fbox), wf(fbox);
    fill_patch(uf, u[0][lev], interface[mglev], boundary.velocity(0),
	       0, 0, icor);
    fill_patch(vf, u[1][lev], interface[mglev], boundary.velocity(1),
	       0, 0, icor);
    fill_patch(wf, u[2][lev], interface[mglev], boundary.velocity(2),
	       0, 0, icor);
    Box creg = interface[mglev].corner(icor);
    creg.coarsen(rat);
    interface[mglev].geo_array(ga, 0, icor);
    FORT_HGCDIV(sptr, dimlist(sbox),
		uc.dataPtr(), vc.dataPtr(), wc.dataPtr(), dimlist(cbox),
		uf.dataPtr(), vf.dataPtr(), wf.dataPtr(), dimlist(fbox),
		dimlist(creg), hx, hy, hz, rat[0], rat[1], rat[2], ga);
    if (jgrid < 0) {
      delete ucp;
      delete vcp;
      delete wcp;
    }
    // fill in the grids on the other sides, if any
    const Box& freg = interface[mglev].corner(icor);
    for (i = 1; i < N_CORNER_GRIDS; i++) {
      jgrid = interface[mglev].cgrid(icor, i);
      if (jgrid >= 0 && jgrid != igrid)
	internal_copy(source[lev], jgrid, igrid, freg);
    }
  }

#else

  for (int icor = 0; icor < interface[mglev].ncorners(); icor++) {
    // find a fine grid touching this corner
    for (i = 0; i < 4; i++) {
      igrid = interface[mglev].cgrid(icor, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = interface[mglev].cgeo(icor);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].cflag(icor) == 1)
      continue;
    else if (geo == XL || geo == XH || geo == YL || geo == YH) {
      // fine grid on two adjacent sides
      int idim = (geo == XL || geo == XH) ? 0 : 1;
      int idir = (geo & LL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      Box cbox = interface[mglev].corner(icor);
      Box fbox = interface[mglev].corner(icor);
      cbox.grow(1 - idim, rat[1-idim]);
      fbox.grow(1 - idim, rat[1-idim]);
      if (idir > 0) {
	cbox.growLo(idim, rat[idim]);
	fbox.growHi(idim, 1);
      }
      else {
	cbox.growHi(idim, rat[idim]);
	fbox.growLo(idim, 1);
      }
      cbox.convert(cellvect).coarsen(rat);
      Fab *ucp, *vcp;
      jgrid = find_patch(cbox, u[0][lev-1]);
      if (jgrid < 0) {
	ucp = new Fab(cbox);
	vcp = new Fab(cbox);
	fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
	fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
      }
      else {
	ucp = &u[0][lev-1][jgrid];
	vcp = &u[1][lev-1][jgrid];
	cbox = ucp->box();
      }
      Fab& uc = *ucp;
      Fab& vc = *vcp;
      fbox.convert(cellvect);
      Fab uf(fbox), vf(fbox);
      fill_patch(uf, u[0][lev], interface[mglev], boundary.velocity(0),
		 0, 0, icor);
      fill_patch(vf, u[1][lev], interface[mglev], boundary.velocity(1),
		 0, 0, icor);
      Box creg = interface[mglev].corner(icor);
      creg.coarsen(rat);
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGFDIV(sptr, dimlist(sbox),
		  uc.dataPtr(), vc.dataPtr(), dimlist(cbox),
		  uf.dataPtr(), vf.dataPtr(), dimlist(fbox),
		  dimlist(creg), hx, hy,
		  rat[0], rat[1], idim, idir,
		  IsRZ(), mg_domain[mgc].bigEnd(0) + 1);
      if (jgrid < 0) {
	delete ucp;
	delete vcp;
      }
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      for (i = 1; i < 4; i++) {
	jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(source[lev], jgrid, igrid, freg);
      }
    }
    else if (geo == LL || geo == HL || geo == LH || geo == HH) {
      // outside corner
      int idir0 = (geo & XL) ? -1 : 1;
      int idir1 = (geo & YL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      Box cbox = interface[mglev].corner(icor);
      cbox.grow(rat).convert(cellvect).coarsen(rat);
      Fab *ucp, *vcp;
      jgrid = find_patch(cbox, u[0][lev-1]);
      if (jgrid < 0) {
	ucp = new Fab(cbox);
	vcp = new Fab(cbox);
	fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
	fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
      }
      else {
	ucp = &u[0][lev-1][jgrid];
	vcp = &u[1][lev-1][jgrid];
	cbox = ucp->box();
      }
      Fab& uc = *ucp;
      Fab& vc = *vcp;
      const Box& fbox = u[0][lev][igrid].box();
      Box creg = interface[mglev].corner(icor);
      creg.coarsen(rat);
      Real *const sptr = source[lev][igrid].dataPtr();
      Real *const u0ptr = u[0][lev][igrid].dataPtr();
      Real *const u1ptr = u[1][lev][igrid].dataPtr();
      FORT_HGODIV(sptr, dimlist(sbox),
		  uc.dataPtr(), vc.dataPtr(), dimlist(cbox),
		  u0ptr, u1ptr,
		  dimlist(fbox), dimlist(creg),
		  hx, hy, rat[0], rat[1], idir0, idir1, IsRZ());
      if (jgrid < 0) {
	delete ucp;
	delete vcp;
      }
    }
    else if (geo == (LL | HH) || geo == (LH | HL)) {
      // diagonal corner
      const Box& sbox = source[lev][igrid].box();
      Box cbox = interface[mglev].corner(icor);
      cbox.grow(rat).convert(cellvect).coarsen(rat);
      Fab *ucp, *vcp;
      jgrid = find_patch(cbox, u[0][lev-1]);
      if (jgrid < 0) {
	ucp = new Fab(cbox);
	vcp = new Fab(cbox);
	fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
	fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
      }
      else {
	ucp = &u[0][lev-1][jgrid];
	vcp = &u[1][lev-1][jgrid];
	cbox = ucp->box();
      }
      Fab& uc = *ucp;
      Fab& vc = *vcp;
      Box fbox = interface[mglev].corner(icor);
      fbox.grow(rat).convert(cellvect);
      Fab uf(fbox), vf(fbox);
      fill_patch(uf, u[0][lev], interface[mglev], boundary.velocity(0),
		 0, 0, icor);
      fill_patch(vf, u[1][lev], interface[mglev], boundary.velocity(1),
		 0, 0, icor);
      int jdir = (geo == (LL | HH)) ? 1 : -1;
      Box creg = interface[mglev].corner(icor);
      creg.coarsen(rat);
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGDDIV(sptr, dimlist(sbox),
		  uc.dataPtr(), vc.dataPtr(), dimlist(cbox),
		  uf.dataPtr(), vf.dataPtr(), dimlist(fbox),
		  dimlist(creg), hx, hy, rat[0], rat[1], jdir, IsRZ());
      if (jgrid < 0) {
	delete ucp;
	delete vcp;
      }
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      for (i = 1; i < 4; i++) {
	jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(source[lev], jgrid, igrid, freg);
      }
    }
    else {
      // inside corner
      const Box& sbox = source[lev][igrid].box();
      Box cbox = interface[mglev].corner(icor);
      Box fbox = interface[mglev].corner(icor);
      cbox.grow(rat).convert(cellvect).coarsen(rat);
      fbox.grow(1).convert(cellvect);
      if ((geo & XL) == XL)
	fbox.growHi(0, rat[0]-1);
      else
	fbox.growLo(0, rat[0]-1);
      if ((geo & YL) == YL)
	fbox.growHi(1, rat[1]-1);
      else
	fbox.growLo(1, rat[1]-1);
      Fab uf(fbox), vf(fbox);
      fill_patch(uf, u[0][lev], interface[mglev], boundary.velocity(0),
		 0, 0, icor);
      fill_patch(vf, u[1][lev], interface[mglev], boundary.velocity(1),
		 0, 0, icor);
      int idir0, idir1;
      if ((geo & XL) == XL) {
	idir0 = -1;
	cbox.growLo(0, -1);
      }
      else {
	idir0 = 1;
	cbox.growHi(0, -1);
      }
      if ((geo & YL) == YL) {
	idir1 = -1;
	cbox.growLo(1, -1);
      }
      else {
	idir1 = 1;
	cbox.growHi(1, -1);
      }
      Fab *ucp, *vcp;
      jgrid = find_patch(cbox, u[0][lev-1]);
      if (jgrid < 0) {
	ucp = new Fab(cbox);
	vcp = new Fab(cbox);
	fill_patch(*ucp, u[0][lev-1], interface[mgc], boundary.velocity(0));
	fill_patch(*vcp, u[1][lev-1], interface[mgc], boundary.velocity(1));
      }
      else {
	ucp = &u[0][lev-1][jgrid];
	vcp = &u[1][lev-1][jgrid];
	cbox = ucp->box();
      }
      Fab& uc = *ucp;
      Fab& vc = *vcp;
      cbox = uc.box();
      Box creg = interface[mglev].corner(icor);
      creg.coarsen(rat);
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGIDIV(sptr, dimlist(sbox),
		  uc.dataPtr(), vc.dataPtr(), dimlist(cbox),
		  uf.dataPtr(), vf.dataPtr(), dimlist(fbox),
		  dimlist(creg), hx, hy, rat[0], rat[1], idir0, idir1, IsRZ());
      if (jgrid < 0) {
	delete ucp;
	delete vcp;
      }
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      for (i = 1; i < 4; i++) {
	jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(source[lev], jgrid, igrid, freg);
      }
    }
  }
#endif
}

void holy_grail_amr_projector::form_solution_vector(PArray<MultiFab>* u,
						    PArray<MultiFab>& sigma_in)
{
  int i, lev;

  if (u) {
    for (lev = lev_min; lev <= lev_max; lev++) {
      int mglev = ml_index[lev];
      Real hx = h[mglev][0];
      Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
      Real hz = h[mglev][2];
#endif
      for (int igrid = 0; igrid < ml_mesh[lev].length(); igrid++) {
	const Box& gbox = ml_mesh[lev][igrid];
	const Box& dbox = dest[lev][igrid].box();
	Fab gp[BL_SPACEDIM];
	for (i = 0; i < BL_SPACEDIM; i++) {
	  gp[i].resize(gbox);
	}
#if (BL_SPACEDIM == 2)
	FORT_HGGRAD(gp[0].dataPtr(), gp[1].dataPtr(), dimlist(gbox),
		    dest[lev][igrid].dataPtr(), dimlist(dbox),
		    dimlist(gbox), hx, hy, IsRZ());
#else
	FORT_HGGRAD(gp[0].dataPtr(), gp[1].dataPtr(), gp[2].dataPtr(),
		    dimlist(gbox),
		    dest[lev][igrid].dataPtr(), dimlist(dbox),
		    dimlist(gbox), hx, hy, hz);
#endif
	for (i = 0; i < BL_SPACEDIM; i++) {
#ifndef CONSTANT
	  gp[i].mult(sigma_in[lev][igrid]);
#endif
	  u[i][lev][igrid].minus(gp[i]);
	}
      }
    }

    sync_periodic_interfaces();
    for (lev = lev_max; lev > lev_min; lev--) {
      const IntVect& rat = gen_ratio[lev-1];
      restrict_level(dest[lev-1], 0, dest[lev], rat,
		     dest_bcache[lev], injection_restrictor);
      for (i = 0; i < BL_SPACEDIM; i++) {
	restrict_level(u[i][lev-1], 0, u[i][lev], rat);
      }
    }
  }
  else {
    sync_periodic_interfaces();
    for (lev = lev_max; lev > lev_min; lev--) {
      restrict_level(dest[lev-1], 0, dest[lev], gen_ratio[lev-1],
		     dest_bcache[lev], injection_restrictor);
    }
  }
}
