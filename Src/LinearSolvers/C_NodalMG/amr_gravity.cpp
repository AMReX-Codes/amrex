
#include "amr_gravity.H"

#ifdef CRAY
#  define   FORT_GVGRAD     GVGRAD
#  define   FORT_GVAVG      GVAVG
#  define   FORT_GVFAVG     GVFAVG
#  define   FORT_GVEAVG     GVEAVG
#  define   FORT_GVCAVG     GVCAVG
#else
#  define   FORT_GVGRAD     gvgrad_
#  define   FORT_GVAVG      gvavg_
#  define   FORT_GVFAVG     gvfavg_
#  define   FORT_GVEAVG     gveavg_
#  define   FORT_GVCAVG     gvcavg_
#endif

extern "C" {

#if (SPACEDIM == 1)
  ERROR, not relevant
#elif (SPACEDIM == 2 || SPACEDIM == 3)
  void FORT_GVGRAD(Real*, intS, const Real*, intS, intS, const Real*);
  void FORT_GVAVG(Real*, intS, Real*, intS, intS);
  void FORT_GVFAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, int&, int&, int&);
#  if (SPACEDIM == 3)
  void FORT_GVEAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, int&, const int*, const int*);
#  endif
  void FORT_GVCAVG(Real*, intS, Real*, intS, Real*, intS,
		   intS, int&, const int*);
#endif
};

DEFINE_ERROR_FUNC(amr_gravity_module)

void amr_gravity_module::poisson(const amr_real& Dest,
				 const amr_real& Source,
				 Real H[], Real reltol,
				 int Lev_min, int Lev_max,
				 Real abstol)
{
  if (Lev_min < 0)
    Lev_min = lev_min_min;
  if (Lev_max < 0)
    Lev_max = lev_max_max;

#ifndef CONSTANT
  error("poisson---only implemented for constant-coefficient version");
#endif

  alloc(Dest, null_amr_real, null_amr_real, H, Lev_min, Lev_max);
  right_hand_side(Source);
  solve(reltol, abstol, 2, 2);
  for (int lev = lev_max; lev > lev_min; lev--) {
    dest[lev].restrict_level(dest[lev-1]);
  }

  if ((lev = lev_min) == 0) {
    int mglev = mg_mesh.which_level(ml_mesh[lev].sig());
    corr[mglev].assign(1.0);
    Real adjust = dest[lev].inner_product(corr[mglev]) /
                  ml_mesh[lev].domain().numPts();
    dest.minus(adjust);
  }

  clear();
}

// Averages cell-based Source onto node-based source conservatively
// across the composite mesh.  Source must be passed in with a one
// cell wide border.  Border values passed in may be meaningless, but
// should not be NaNs.

void amr_gravity_module::right_hand_side(const amr_real& Source)
{ 
  if (Source.border() != 1)
    error("right_hand_side---width-1 border on source currently required");

  int lev, igrid;
  Real adjust = 0.0;
  if (singular) {
    for (lev = lev_max; lev > lev_min; lev--) {
      Source[lev].restrict_level(Source[lev-1]);
    }
    for (igrid = 0; igrid < ml_mesh[lev_min].ngrids(); igrid++) {
      const grid_real tmp(Source[lev_min][igrid]);
#ifdef GUTHAMR
      adjust +=	tmp.fab().sum(ml_mesh[lev_min][igrid], tmp.component());
#else
      adjust +=	tmp.fab().sum(tmp.component(), ml_mesh[lev_min][igrid]);
#endif
    }
    adjust /= ml_mesh[lev_min].domain().numPts();
    if (pcode >= 2)
      cout << "Cell solvability adjustment: " << adjust << endl;

    for (lev = lev_min; lev <= lev_max; lev++) {
      Source[lev].minus(adjust);
    }
  }

  for (lev = lev_min; lev <= lev_max; lev++) {
    int mglev = mg_mesh.which_level(ml_mesh[lev].sig());

    Source[lev].fill_borders(interface[mglev], boundary.scalar());

    for (igrid = 0; igrid < ml_mesh[lev].ngrids(); igrid++) {
      const Box& sbox = source.box(lev, igrid);
      const Box& fbox = Source.box(lev, igrid);
      const Box& freg = interface[mglev].part_fine(igrid);
      Real *const sptr = source.dataPtr(lev, igrid);
      Real *const csptr = Source.dataPtr(lev, igrid);
      FORT_GVAVG(sptr, dimlist(sbox),
		 csptr, dimlist(fbox), dimlist(freg));
    }

    source[lev].clear_part_interface(interface[mglev]);

    if (lev > lev_min)
      interface_average(Source, lev);
  }
}

void amr_gravity_module::interface_average(const amr_real& Source, int lev)
{ 
  int mglev = mg_mesh.which_level(ml_mesh[lev].sig());
  int mgc = mg_mesh.which_level(ml_mesh[lev-1].sig());
  int igrid, jgrid, i;

  DECLARE_GEOMETRY_TYPES;

  int rat = ml_mesh.ratio(lev);
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
    const Box& sbox = source.box(lev, igrid);
    const Box& fbox = Source.box(lev, igrid);
    Box cbox = interface[mglev].face(iface);
    Intvect t = cbox.type();
    if (idir > 0)
      cbox.growLo(idim, rat);
    else
      cbox.growHi(idim, rat);
    cbox.convert(cellvect).coarsen(rat);
    grid_real Sc = Source[lev-1].get_patch(cbox, interface[mgc],
					   boundary.scalar());
    cbox = Sc.box();
    Box creg = interface[mglev].node_face(iface);
    creg.coarsen(rat).grow(t - UNITV);
    Real *const sptr = source.dataPtr(lev, igrid);
    Real *const Sfptr = Source.dataPtr(lev, igrid);
    FORT_GVFAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sfptr, dimlist(fbox), dimlist(creg),
		rat, idim, idir);
  }

  int ga[N_CORNER_GRIDS];

#if (SPACEDIM == 3)

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
    const Box& sbox = source.box(lev, igrid);
    Real *const sptr = source.dataPtr(lev, igrid);
    Box cbox = interface[mglev].edge(iedge);
    Intvect t = cbox.type();
    cbox.coarsen(rat).grow(t).convert(cellvect);
    Box fbox = cbox;
    fbox.refine(rat);
    grid_real Sc = Source[lev-1].get_patch(cbox, interface[mgc],
					   boundary.scalar());
    cbox = Sc.box();
    grid_real Sf(fbox);
    Source[lev].fill_patch(Sf, interface[mglev], boundary.scalar(),
			   0, 1, iedge);
    Box creg = interface[mglev].node_edge(iedge);
    creg.coarsen(rat).grow(t - UNITV);
    interface[mglev].geo_array(ga, 1, iedge);
    FORT_GVEAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sf.dataPtr(), dimlist(fbox),
		dimlist(creg), rat, t.getVect(), ga);
    // fill in the grids on the other sides, if any
    const Box& freg = interface[mglev].node_edge(iedge);
    for (i = 1; i < N_EDGE_GRIDS; i++) {
      jgrid = interface[mglev].egrid(iedge, i);
      if (jgrid >= 0 && jgrid != igrid)
	source[lev].internal_copy(jgrid, igrid, freg);
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
    const Box& sbox = source.box(lev, igrid);
    Real *const sptr = source.dataPtr(lev, igrid);
    Box cbox = interface[mglev].corner(icor);
    cbox.coarsen(rat).grow(1).convert(cellvect);
    Box fbox = cbox;
    fbox.refine(rat);
    grid_real Sc = Source[lev-1].get_patch(cbox, interface[mgc],
					   boundary.scalar());
    cbox = Sc.box();
    grid_real Sf(fbox);
    Source[lev].fill_patch(Sf, interface[mglev], boundary.scalar(),
			   0, 0, icor);
    Box creg = interface[mglev].corner(icor);
    creg.coarsen(rat);
    interface[mglev].geo_array(ga, 0, icor);
    FORT_GVCAVG(sptr, dimlist(sbox),
		Sc.dataPtr(), dimlist(cbox),
		Sf.dataPtr(), dimlist(fbox),
		dimlist(creg), rat, ga);
    // fill in the grids on the other sides, if any
    const Box& freg = interface[mglev].corner(icor);
    for (i = 1; i < N_CORNER_GRIDS; i++) {
      jgrid = interface[mglev].cgrid(icor, i);
      if (jgrid >= 0 && jgrid != igrid)
	source[lev].internal_copy(jgrid, igrid, freg);
    }
  }
}

void amr_gravity_module::gradient(Real* grad_ptr,
                                  const Box& grad_box,
                                  const Real* potential_ptr,
                                  const Box& potential_box,
                                  const Box& region,
                                  const Real* dx)
{
#if (SPACEDIM == 2)
  error("gradient---not implemented for SPACEDIM == 2");
#else
  FORT_GVGRAD(grad_ptr, dimlist(grad_box),
              potential_ptr, dimlist(potential_box),
              dimlist(region), dx);
#endif
}
