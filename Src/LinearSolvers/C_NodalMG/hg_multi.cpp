
#include "hg_multi.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGSRST     hgsrst_
#  define   FORT_HGSCON     hgscon_
#  define   FORT_HGCEN      hgcen_
#  define   FORT_HGINTS     hgints_
#  define   FACRST1         acrst1_
#  define   FANRST2         anrst2_
#  define   FANINT2         anint2_
#else
#  define   FORT_HGSRST     HGSRST
#  define   FORT_HGSCON     HGSCON
#  define   FORT_HGCEN      HGCEN
#  define   FORT_HGINTS     HGINTS
#  define   FACRST1         ACRST1
#  define   FANRST2         ANRST2
#  define   FANINT2         ANINT2
#endif

extern "C" {

#if (BL_SPACEDIM == 1)
  ERROR, not relevant
#elif (BL_SPACEDIM == 2 || BL_SPACEDIM == 3)
#  ifdef TERRAIN
  void FACRST1(Real*, intS, intS, Real*, intS, intRS, const int&);
  void FORT_HGSRST(RealPS, intS, intS, RealPS, intS, intRS);
  void FORT_HGCEN(Real*, intS, Real*, intS, intS);
  void FORT_HGINTS(Real*, intS, intS, RealPS, intS, Real*, intS, intS, intRS);
#  elif (! defined CONSTANT)
  void FORT_HGSRST(RealPS, intS, intS, RealPS, intS, intRS);
#    ifndef SIGMA_NODE
  void FORT_HGCEN(Real*, intS, RealPS, intS, intS, RealRS,
		  const int&, const int&);
  void FORT_HGINTS(Real*, intS, intS, RealPS, intS, Real*, intS, intS, intRS);
#    else
  void FORT_HGSCON(Real*, intS, RealPS, intS, intS, RealRS);
  void FORT_HGCEN(Real*, intS, Real*, intS, intS);
  void FORT_HGINTS(Real*, intS, intS, Real*, intS, Real*, intS, intS, intRS);
#    endif
#  endif
  void FANRST2(Real*, intS, intS, Real*, intS, intRS, const int&);
  void FANINT2(Real*, intS, intS, Real*, intS, intS, intRS);
#endif
}

void holy_grail_amr_multigrid::alloc(PArray<MultiFab>& Dest,
				     PArray<MultiFab>& Source,
				     PArray<MultiFab>& Coarse_source,
				     PArray<MultiFab>& Sigma,
				     Real H[], int Lev_min, int Lev_max)
{
  int lev, mglev, i;

  assert(Dest.length() > Lev_max);
  assert(Dest[Lev_min].nGrow() == 1);

  if (Source.ready()) {
    source_owned = 0;
    amr_multigrid::alloc(Dest, Source, Coarse_source, Lev_min, Lev_max);
  }
  else {
    source_owned = 1;
    PArray<MultiFab> Src;
    Src.resize(Lev_max + 1);
    for (lev = Lev_min; lev <= Lev_max; lev++) {
      const BoxArray& mesh = Dest[lev].boxArray();
      Src.set(lev, new MultiFab(mesh, 1, Dest[Lev_min].nGrow()));
      Src[lev].setVal(0.0);
    }
    amr_multigrid::alloc(Dest, Src,    Coarse_source, Lev_min, Lev_max);
  }

  h = new Real[mglev_max + 1][BL_SPACEDIM];
  for (i = 0; i < BL_SPACEDIM; i++) {
    h[mglev_max][i] = H[i];
#ifdef CONSTANT
    assert(H[i] == H[0]);
#endif
    for (mglev = mglev_max - 1; mglev >= 0; mglev--) {
      int rat = mg_domain[mglev+1].length(i) / mg_domain[mglev].length(i);
      h[mglev][i] = rat * h[mglev+1][i];
    }
  }

  corr_scache.resize(mglev_max + 1, (copy_cache*) 0);

  for (lev = lev_min; lev <= lev_max; lev++) {
    mglev = ml_index[lev];
    dest_bcache.set(lev, new copy_cache(dest[lev], interface[mglev],
					mg_boundary, 1));
  }
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    corr_bcache.set(mglev, new copy_cache(corr[mglev], interface[mglev],
					  mg_boundary, 1));
    corr_scache.set(mglev, new copy_cache(corr[mglev], interface[mglev],
					  mg_boundary));
  }
  for (mglev = 1; mglev <= mglev_max; mglev++) {
    work_bcache.set(mglev, new copy_cache(work[mglev], interface[mglev],
					  mg_boundary, 1));
  }

  build_sigma(Sigma);

  alloc_sync_caches();

  int ib = dest[lev_min].nGrow();
  const BoxArray& mesh0 = corr[0].boxArray();

  cgwork.resize(8);
  cgwork.set(0, new MultiFab(mesh0, 1, ib));
  cgwork.set(1, new MultiFab(mesh0, 1, ib));
  cgwork.set(2, new MultiFab(mesh0, 1, ib));
  cgwork.set(3, &corr[0]);
  cgwork.set(4, &work[0]);
  cgwork.set(5, &cen[0]);
  cgwork.set(6, new MultiFab(mesh0, 1, ib));
  cgwork[6].setVal(0.0);
  cgwork.set(7, new MultiFab(mesh0, 1, ib));

  cgw1_bcache = new copy_cache(cgwork[1], interface[0], mg_boundary, 1);

  assert(cgwork[3].nGrow() == ib &&
	 cgwork[4].nGrow() == ib &&
	 cgwork[5].nGrow() == ib);

  cgw_ucache.resize(8);
  for (i = 0; i < 8; i++) {
    cgw_ucache.set(i, new unroll_cache(cgwork[i]));
  }

  for (int igrid = 0; igrid < mg_mesh[0].length(); igrid++) {
    Fab& gtmp = cgwork[7][igrid];
    const Box& valid = cgwork[7].box(igrid);
    gtmp.setVal(0.0);
    gtmp.setVal(1.0, valid, 0);
    Box b = bdryLo(valid, 1);
    gtmp.setVal(0.5, b, 0);
    b = bdryHi(valid, 1);
    gtmp.setVal(0.5, b, 0);
    b = bdryLo(valid, 0);
    gtmp.setVal(0.5, b, 0);
    gtmp.setVal(0.25, bdryLo(b, 1), 0);
    gtmp.setVal(0.25, bdryHi(b, 1), 0);
    b = bdryHi(valid, 0);
    gtmp.setVal(0.5, b, 0);
    gtmp.setVal(0.25, bdryLo(b, 1), 0);
    gtmp.setVal(0.25, bdryHi(b, 1), 0);
#if (BL_SPACEDIM == 3)
    b = bdryLo(valid, 2);
    gtmp.setVal(0.5, b, 0);
    gtmp.setVal(0.25, bdryLo(b, 0), 0);
    gtmp.setVal(0.25, bdryHi(b, 0), 0);
    Box bb = bdryLo(b, 1);
    gtmp.setVal(0.25, bb, 0);
    gtmp.setVal(0.125, bdryLo(bb, 0), 0);
    gtmp.setVal(0.125, bdryHi(bb, 0), 0);
    bb = bdryHi(b, 1);
    gtmp.setVal(0.25, bb, 0);
    gtmp.setVal(0.125, bdryLo(bb, 0), 0);
    gtmp.setVal(0.125, bdryHi(bb, 0), 0);
    b = bdryHi(valid, 2);
    gtmp.setVal(0.5, b, 0);
    gtmp.setVal(0.25, bdryLo(b, 0), 0);
    gtmp.setVal(0.25, bdryHi(b, 0), 0);
    bb = bdryLo(b, 1);
    gtmp.setVal(0.25, bb, 0);
    gtmp.setVal(0.125, bdryLo(bb, 0), 0);
    gtmp.setVal(0.125, bdryHi(bb, 0), 0);
    bb = bdryHi(b, 1);
    gtmp.setVal(0.25, bb, 0);
    gtmp.setVal(0.125, bdryLo(bb, 0), 0);
    gtmp.setVal(0.125, bdryHi(bb, 0), 0);
#endif
  }

  singular = 0;
  if (mg_boundary.singular()) {
    for (i = 0; i < mg_mesh[0].length(); i++) {
      singular += mg_mesh[0][i].numPts();
    }
    singular = (singular == mg_domain[0].numPts());
  }

#ifdef TERRAIN
  integrate = 1;
#endif
}

#ifndef CONSTANT

void holy_grail_sigma_restrictor_class::fill(Fab& patch,
					     const Box& region,
					     Fab& fgr,
					     const IntVect& rat) const
{
  assert(patch.box().cellCentered());
  assert(rat[0] == 2 && rat[1] == 2 ||
	 rat[0] == 2 && rat[1] == 1 ||
	 rat[0] == 1 && rat[1] == 2);

#ifdef TERRAIN
    FORT_HGSRST(patch.dataPtr(0), patch.dataPtr(1),
#  if (BL_SPACEDIM == 3)
		patch.dataPtr(2),
#  endif
		dimlist(patch.box()),
		dimlist(region),
		fgr.dataPtr(0), fgr.dataPtr(1),
#  if (BL_SPACEDIM == 3)
		fgr.dataPtr(2),
#  endif
		dimlist(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
    FACRST1(patch.dataPtr(BL_SPACEDIM),
		dimlist(patch.box()),
		dimlist(region),
		fgr.dataPtr(BL_SPACEDIM),
		dimlist(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]), 0);

#  if (BL_SPACEDIM == 2)
    patch.mult((Real) rat[1] / rat[0], region, 0, 1);
    patch.mult((Real) rat[0] / rat[1], region, 1, 1);
    // component 2 remains unchanged
#  else
    FACRST1(patch.dataPtr(BL_SPACEDIM+1),
		dimlist(patch.box()),
		dimlist(region),
		fgr.dataPtr(BL_SPACEDIM+1),
		dimlist(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]), 0);
    patch.mult((Real) rat[1] * rat[2] / rat[0], region, 0, 1);
    patch.mult((Real) rat[0] * rat[2] / rat[1], region, 1, 1);
    patch.mult((Real) rat[0] * rat[1] / rat[2], region, 2, 1);
    patch.mult((Real) rat[1],                   region, 3, 1);
    patch.mult((Real) rat[0],                   region, 4, 1);
#  endif

#else

  if (fgr.nVar() == 1) {
    FORT_HGSRST(patch.dataPtr(0), patch.dataPtr(1),
#  if (BL_SPACEDIM == 3)
		patch.dataPtr(2),
#  endif
		dimlist(patch.box()),
		dimlist(region),
		fgr.dataPtr(), fgr.dataPtr(),
#  if (BL_SPACEDIM == 3)
		fgr.dataPtr(),
#  endif
		dimlist(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
  }
  else {
    FORT_HGSRST(patch.dataPtr(0), patch.dataPtr(1),
#  if (BL_SPACEDIM == 3)
		patch.dataPtr(2),
#  endif
		dimlist(patch.box()),
		dimlist(region),
		fgr.dataPtr(0), fgr.dataPtr(1),
#  if (BL_SPACEDIM == 3)
		fgr.dataPtr(2),
#  endif
		dimlist(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
  }

#endif
}

#endif

void holy_grail_amr_multigrid::build_sigma(PArray<MultiFab>& Sigma)
{
  int mglev, igrid;

#ifdef TERRAIN

  // For terrain stencils we have as many sigma arrays passed as
  // arguments and used at the interface as we build for internal
  // multigrid purposes.  This simplifies handling as we do not
  // need to maintain separate arrays for different purposes.

  int lev;
  int ncomp = 2 * BL_SPACEDIM - 1;

  sigma.resize(mglev_max+1);

  for (mglev = 0; mglev <= mglev_max; mglev++) {
    sigma.set(mglev, new MultiFab(mg_mesh[mglev], ncomp, 1));
    MultiFab& target = sigma[mglev];
    target.setVal(1.e20);
    if ((lev = get_amr_level(mglev)) >= 0) {
      MultiFab& s_in = Sigma[lev];
      for (igrid = 0; igrid < target.length(); igrid++) {
	target[igrid].copy(s_in[igrid], s_in.box(igrid), 0,
			   target.box(igrid), 0, ncomp);
      }
    }
  }

  for (mglev = mglev_max; mglev > 0; mglev--) {
    IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
    restrict_level(sigma[mglev-1], sigma[mglev], rat, 0,
		   holy_grail_sigma_restrictor);
  }
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    fill_borders(sigma[mglev], 0, interface[mglev], boundary.terrain_sigma());
  }

#elif (! defined CONSTANT)

  // Intended functionality:  sigma_split exists only at coarser levels,
  // since only after coarsening is sigma different in different directions.
  // sigma exists at all levels, and is intended for use on fine grids
  // and at interface points, where all components are the same.  To
  // save storage it is aliased to the first component of sigma_split
  // on all but the finest level.

  // sigma_split replaced by sigma_nd in more recent version, used
  // only as a local variable here

  int lev, i;
  PArray<MultiFab> sigma_split;
  sigma_split.resize(mglev_max);
  for (mglev = 0; mglev < mglev_max; mglev++) {
    sigma_split.set(mglev, new MultiFab(mg_mesh[mglev], BL_SPACEDIM, 1));
  }
  sigma.resize(mglev_max + 1);
  mglev = mglev_max;
  sigma.set(mglev, new MultiFab(mg_mesh[mglev], 1, 1));

  // Level project:
  // Any value can fill values in the border cells that fill_borders
  // will not touch---those touching coarser grids.  The values in these
  // cells will only be seen by the interpolation, and the quantity
  // being interpolated will always be zero at these edges, but we
  // should insure that no NaN's or other garbage is there that could
  // cause a floating point fault.

  // Sync project:
  // Ghost values will be seen by multilevel interpolation, so put
  // a huge value in ghost cells so that coarse-fine interface
  // interpolation will linear, as finite element derivation requires.

  for (mglev = 0; mglev < mglev_max; mglev++) {
    MultiFab& target = sigma_split[mglev];
    target.setVal(1.e20);
    if ((lev = get_amr_level(mglev)) >= 0) {
      MultiFab& s_comp = Sigma[lev];
      for (i = 0; i < BL_SPACEDIM; i++) {
	for (igrid = 0; igrid < target.length(); igrid++) {
	  target[igrid].copy(s_comp[igrid], s_comp.box(igrid), 0,
			     target.box(igrid), i, 1);
	}
      }
    }
  }

  mglev = mglev_max;
  sigma[mglev].setVal(1.e20);
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    sigma[mglev][igrid].copy(Sigma[lev_max][igrid], mg_mesh[mglev][igrid], 0,
			     mg_mesh[mglev][igrid], 0, 1);
  }

  mglev = mglev_max;
  if (mglev_max > 0) {
    IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
    restrict_level(sigma_split[mglev-1], sigma[mglev], rat, 0,
		   holy_grail_sigma_restrictor);
  }
  fill_borders(sigma[mglev], 0, interface[mglev], boundary.scalar());
  for (mglev = mglev_max - 1; mglev > 0; mglev--) {
    IntVect rat = mg_domain[mglev].length() / mg_domain[mglev-1].length();
    restrict_level(sigma_split[mglev-1], sigma_split[mglev], rat, 0,
		   holy_grail_sigma_restrictor);
  }
  for (mglev = 0; mglev < mglev_max; mglev++) {
    fill_borders(sigma_split[mglev], 0, interface[mglev], boundary.scalar());
  }

  for (i = 0; i < BL_SPACEDIM; i++) {
    sigma_nd[i].resize(mglev_max + 1);
  }

  for (mglev = 0; mglev < mglev_max; mglev++) {
    MultiFab& s = sigma_split[mglev];
    for (i = 0; i < BL_SPACEDIM; i++) {
      sigma_nd[i].set(mglev, new MultiFab(mg_mesh[mglev], 1, 1));
      MultiFab& d = sigma_nd[i][mglev];
      for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
	d[igrid].copy(s[igrid], i, 0);
      }
    }
    delete sigma_split.remove(mglev);
    sigma.set(mglev, &sigma_nd[0][mglev]);
  }

  mglev = mglev_max;
  for (i = 0; i < BL_SPACEDIM; i++) {
    sigma_nd[i].set(mglev, &sigma[mglev]);
  }

#  ifdef SIGMA_NODE

  sigma_node.resize(mglev_max + 1);
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    BoxArray mesh = mg_mesh[mglev];
    mesh.convert(nodevect);
    sigma_node.set(mglev, new MultiFab(mesh, BL_SPACEDIM, 1));
    sigma_node[mglev].setVal(1.e20);
  }

  for (mglev = 0; mglev <= mglev_max; mglev++) {
    const Real hx = h[mglev][0];
    const Real hy = h[mglev][1];
#    if (BL_SPACEDIM == 3)
    const Real hz = h[mglev][2];
#    endif
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& scbox = sigma[mglev][igrid].box();
      const Box& snbox = sigma_node[mglev][igrid].box();
      const Box& reg = interface[mglev].part_fine(igrid);
      FORT_HGSCON(sigma_node[mglev][igrid].dataPtr(),
                  dimlist(snbox),
                  sigma_nd[0][mglev][igrid].dataPtr(),
                  sigma_nd[1][mglev][igrid].dataPtr(),
#    if (BL_SPACEDIM == 3)
                  sigma_nd[2][mglev][igrid].dataPtr(),
#    endif
                  dimlist(scbox),
                  dimlist(reg),
#    if (BL_SPACEDIM == 2)
                  hx, hy
#    else
                  hx, hy, hz
#    endif
		  );
    }

    if (mglev < mglev_max) {
      sigma_nd[0].remove(mglev);
      for (i = 1; i < BL_SPACEDIM; i++) {
	delete sigma_nd[i].remove(mglev);
      }
    }
    else {
      for (i = 0; i < BL_SPACEDIM; i++) {
	sigma_nd[i].remove(mglev);
      }
    }
  }
#  endif  // SIGMA_NODE

#endif

  DECLARE_GEOMETRY_TYPES;

  cen.resize(mglev_max + 1);
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    cen.set(mglev, new MultiFab(corr[mglev].boxArray(), 1,
				dest[lev_min].nGrow()));
    MultiFab& ctmp = cen[mglev];
    ctmp.setVal(0.0);

#ifdef TERRAIN

    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& cenbox = ctmp[igrid].box();
      const Box& reg = interface[mglev].part_fine(igrid);
      const Box& sigbox = sigma[mglev][igrid].box();
      FORT_HGCEN(ctmp[igrid].dataPtr(), dimlist(cenbox),
		 sigma[mglev][igrid].dataPtr(),
                 dimlist(sigbox),
		 dimlist(reg));
    }

#elif (defined CONSTANT)

    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      ctmp[igrid].setVal(1.0, interface[mglev].part_fine(igrid), 0);
    }

#else

#  ifndef SIGMA_NODE
    const Real hx = h[mglev][0];
    const Real hy = h[mglev][1];
#    if (BL_SPACEDIM == 3)
    const Real hz = h[mglev][2];
#    endif
#  endif
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& cenbox = cen[mglev][igrid].box();
      const Box& reg = interface[mglev].part_fine(igrid);
#  ifndef SIGMA_NODE
      const Box& sigbox = sigma[mglev][igrid].box();
      FORT_HGCEN(cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		 sigma_nd[0][mglev][igrid].dataPtr(),
		 sigma_nd[1][mglev][igrid].dataPtr(),
#    if (BL_SPACEDIM == 3)
		 sigma_nd[2][mglev][igrid].dataPtr(),
#    endif
                 dimlist(sigbox),
		 dimlist(reg),
#    if (BL_SPACEDIM == 2)
		 hx, hy,
		 IsRZ(), mg_domain[mglev].bigEnd(0) + 1
#    else
		 hx, hy, hz
#    endif
		 );
#  else
      const Box& sigbox = sigma_node[mglev][igrid].box();
      FORT_HGCEN(cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		 sigma_node[mglev][igrid].dataPtr(),
                 dimlist(sigbox),
		 dimlist(reg));
#  endif
    }

#endif

    clear_part_interface(ctmp, interface[mglev]);
  }

#ifdef CONSTANT
  mask.resize(mglev_max + 1);
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    mask.set(mglev, &cen[mglev]);
  }
#endif

#ifdef SIGMA_NODE

  mask.resize(mglev_max + 1);
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    mask.set(mglev, new MultiFab(corr[mglev].boxArray(), 1,
				 dest[lev_min].nGrow()));
    MultiFab& mtmp = mask[mglev];
    mtmp.setVal(0.0);
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      mtmp[igrid].setVal(1.0, interface[mglev].part_fine(igrid), 0);
    }
    clear_part_interface(mtmp, interface[mglev]);
  }

#endif
}

void holy_grail_amr_multigrid::clear()
{
  int i, lev, mglev;

  line_order.clear();
  line_after.clear();

  delete_sync_caches();

  delete cgw1_bcache;
  for (i = 0; i < 8; i++) {
    delete cgw_ucache[i];
    cgw_ucache[i] = 0;
  }
  delete cgwork.remove(0);
  delete cgwork.remove(1);
  delete cgwork.remove(2);
  cgwork.remove(3);
  cgwork.remove(4);
  cgwork.remove(5);
  delete cgwork.remove(6);
  delete cgwork.remove(7);

#ifndef CONSTANT
#  ifdef TERRAIN
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    delete sigma.remove(mglev);
  }
#  elif (defined SIGMA_NODE)
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    delete sigma.remove(mglev);
    delete sigma_node.remove(mglev);
  }
#  else
  mglev = mglev_max;
  delete sigma.remove(mglev);
  for (i = 0; i < BL_SPACEDIM; i++) {
    sigma_nd[i].remove(mglev);
  }
  for (mglev = 0; mglev < mglev_max; mglev++) {
    sigma.remove(mglev);
    for (i = 0; i < BL_SPACEDIM; i++) {
      delete sigma_nd[i].remove(mglev);
    }
  }
#  endif
#endif

  for (mglev = 0; mglev <= mglev_max; mglev++) {
    delete cen.remove(mglev);
#ifdef CONSTANT
    mask.remove(mglev);
#endif
#ifdef SIGMA_NODE
    delete mask.remove(mglev);
#endif
  }

  delete [] h;
  if (source_owned) {
    for (lev = lev_min; lev <= lev_max; lev++) {
      if (source.defined(lev)) delete source.remove(lev);
    }
  }

  for (lev = lev_min; lev <= lev_max; lev++) {
    delete dest_bcache[lev];
    dest_bcache[lev] = 0;
  }
  for (mglev = 0; mglev <= mglev_max; mglev++) {
    delete corr_bcache[mglev];
    corr_bcache[mglev] = 0;
    delete corr_scache[mglev];
    corr_scache[mglev] = 0;
  }
  for (mglev = 1; mglev <= mglev_max; mglev++) {
    delete work_bcache[mglev];
    work_bcache[mglev] = 0;
  }

  amr_multigrid::clear();
}

int holy_grail_amr_multigrid::can_coarsen(const BoxArray& mesh,
					  const Box& domain)
{
  int retval = 1;
  for (int i = 0; i < BL_SPACEDIM; i++) {
    retval &= ((domain.smallEnd(i)&1) == 0);
    retval &= ((domain.bigEnd(i)&1)   == 1);
    retval &= (domain.length(i) >= 4);
    for (int igrid = 0; igrid < mesh.length(); igrid++) {
      retval &= ((mesh[igrid].smallEnd(i)&1) == 0 &&
		 (mesh[igrid].bigEnd(i)&1)   == 1 &&
		 (mesh[igrid].length(i) >= 4));
    }
  }
  return retval;
}

void holy_grail_amr_multigrid::sync_interfaces()
{
  DECLARE_GEOMETRY_TYPES;

  for (int lev = lev_min+1; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    int mgc = ml_index[lev-1];
    IntVect rat = mg_domain[mglev].length() / mg_domain[mgc].length();
    MultiFab& target = dest[lev];
    for (int iface = 0; iface < interface[mglev].nfaces(); iface++) {
      // find a fine grid touching this face
      int igrid = interface[mglev].fgrid(iface, 0);
      if (igrid < 0)
	igrid = interface[mglev].fgrid(iface, 1);
      unsigned geo = interface[mglev].fgeo(iface);
      // reject fine-fine interfaces and those without an interior fine grid
      if (geo == ALL || igrid < 0 || interface[mglev].fflag(iface) == 1)
	continue;
      interpolate_patch(target[igrid], interface[mglev].node_face(iface),
			dest[lev-1], rat,
			bilinear_interpolator, interface[mgc]);
    }
  }
}

void holy_grail_amr_multigrid::sync_periodic_interfaces()
{
  DECLARE_GEOMETRY_TYPES;

  for (int lev = lev_min+1; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    int mgc = ml_index[lev-1];
    IntVect rat = mg_domain[mglev].length() / mg_domain[mgc].length();
    Box idomain = mg_domain[mglev];
    idomain.convert(type(dest[lev])).grow(-1);
    MultiFab& target = dest[lev];
    for (int iface = 0; iface < interface[mglev].nfaces(); iface++) {
      // find a fine grid touching this face
      int igrid = interface[mglev].fgrid(iface, 0);
      if (igrid < 0)
	igrid = interface[mglev].fgrid(iface, 1);
      unsigned geo = interface[mglev].fgeo(iface);
      // use only exterior coarse-fine faces with an interior fine grid
      const Box& nbox = interface[mglev].node_face(iface);
      if (geo == ALL || igrid < 0 || interface[mglev].fflag(iface) == 1 ||
	  idomain.intersects(nbox))
	continue;
      interpolate_patch(target[igrid], nbox, dest[lev-1], rat,
			bilinear_interpolator, interface[mgc]);
    }
  }
}

void holy_grail_amr_multigrid::mg_restrict_level(int lto, int lfrom)
{
  IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
  if (get_amr_level(lto) >= 0) {
    if (integrate == 0) {
      restrict_level(resid[lto], 0, work[lfrom], rat, work_bcache[lfrom],
		     bilinear_restrictor_coarse,
		     interface[lfrom], mg_boundary);
    }
    else {
      restrict_level(resid[lto], 0, work[lfrom], rat, work_bcache[lfrom],
		     bilinear_integrator_coarse,
		     interface[lfrom], mg_boundary);
    }
  }
  else {
    mg_restrict(lto, lfrom);
  }
}

void holy_grail_amr_multigrid::mg_restrict(int lto, int lfrom)
{
  int igrid;
  fill_borders(work[lfrom], work_bcache[lfrom], interface[lfrom], mg_boundary);
  IntVect rat = mg_domain[lfrom].length() / mg_domain[lto].length();
  for (igrid = 0; igrid < resid[lto].length(); igrid++) {
    const Box& fbox = work[lfrom][igrid].box();
    const Box& cbox = resid[lto][igrid].box();
    const Box& creg = interface[lto].part_fine(igrid);
    FANRST2(resid[lto][igrid].dataPtr(), dimlist(cbox), dimlist(creg),
	    work[lfrom][igrid].dataPtr(), dimlist(fbox),
	    D_DECL(rat[0], rat[1], rat[2]), integrate);
  }

  clear_part_interface(resid[lto], interface[lto]);
}

#ifndef CONSTANT

void holy_grail_interpolator_class::fill(Fab& patch,
					 const Box& region,
					 Fab& cgr,
					 const Box& cb,
					 const IntVect& rat) const
{
#ifndef SIGMA_NODE

  FORT_HGINTS(patch.dataPtr(), dimlist(patch.box()), dimlist(region),
	      sigptr[0], sigptr[1],
#if (BL_SPACEDIM == 3)
              sigptr[2],
#endif
              dimlist(sigbox),
	      cgr.dataPtr(), dimlist(cgr.box()), dimlist(cb),
	      D_DECL(rat[0], rat[1], rat[2]));

#else

  FORT_HGINTS(patch.dataPtr(), dimlist(patch.box()), dimlist(region),
	      sigptr, dimlist(sigbox),
	      cgr.dataPtr(), dimlist(cgr.box()), dimlist(cb),
	      D_DECL(rat[0], rat[1], rat[2]));

#endif
}

#endif

void holy_grail_amr_multigrid::mg_interpolate_level(int lto, int lfrom)
{
  if (get_amr_level(lfrom) >= 0) {
#ifdef CONSTANT
    // multilevel interpolation, use general form
    amr_multigrid::mg_interpolate_level(lto, lfrom);
#else
    // general version---attempt to use special stencils for multilevel
    int ltmp = lfrom + 1;
    MultiFab& target = work[ltmp];
    IntVect rat = mg_domain[ltmp].length() / mg_domain[lfrom].length();
    for (int igrid = 0; igrid < target.length(); igrid++) {
#ifdef TERRAIN
      Real *sigptr[BL_SPACEDIM];
      for (int i = 0; i < BL_SPACEDIM; i++) {
	sigptr[i] = sigma[ltmp][igrid].dataPtr(i);
      }
      const Box& sigbox = sigma[ltmp][igrid].box();
#elif (! defined SIGMA_NODE)
      Real *sigptr[BL_SPACEDIM];
      for (int i = 0; i < BL_SPACEDIM; i++) {
	sigptr[i] = sigma_nd[i][ltmp][igrid].dataPtr();
      }
      const Box& sigbox = sigma[ltmp][igrid].box();
#else
      Real *sigptr = sigma_node[ltmp][igrid].dataPtr();
      const Box& sigbox = sigma_node[ltmp][igrid].box();
#endif
      interpolate_patch(target[igrid], target.box(igrid),
			corr[lfrom], rat,
			holy_grail_interpolator_class(sigptr, sigbox),
			interface[lfrom], error_boundary);
    }
    if (lto > ltmp) {
      corr[ltmp].copy(target);
      mg_interpolate_level(lto, ltmp);
    }
#endif
  }
  else {
/*
    // this general version not currently used, but may need to be revived
    // if special stencils are ever needed for the multilevel iteration
    MultiFab& target = work[lto];
    for (int igrid = 0; igrid < target.length(); igrid++) {
      Real *sigptr[BL_SPACEDIM];
      for (int i = 0; i < BL_SPACEDIM; i++) {
	sigptr[i] = sigma_nd[i][lto][igrid].dataPtr();
      }
      const Box& sigbox = sigma[lto][igrid].box();
      corr[lfrom].interpolate_patch(target[igrid], target.box(igrid),
	corr[lfrom], rat,
	error_boundary, holy_grail_interpolator_class(sigptr, sigbox));
    }
*/
    // multigrid interpolation, grids known to match up
    // special stencil needed for multigrid convergence
    IntVect rat = mg_domain[lto].length() / mg_domain[lfrom].length();
    for (int igrid = 0; igrid < mg_mesh[lto].length(); igrid++) {
      const Box& fbox = work[lto][igrid].box();
      const Box& freg = work[lto].box(igrid);
      const Box& cbox = corr[lfrom][igrid].box();
      const Box& creg = corr[lfrom].box(igrid);
#ifdef CONSTANT
      FANINT2(work[lto][igrid].dataPtr(), dimlist(fbox), dimlist(freg),
	      corr[lfrom][igrid].dataPtr(), dimlist(cbox), dimlist(creg),
	      D_DECL(rat[0], rat[1], rat[2]));
#else
#ifndef SIGMA_NODE
      const Box& sigbox = sigma[lto][igrid].box();
#else
      const Box& sigbox = sigma_node[lto][igrid].box();
#endif
      FORT_HGINTS(work[lto][igrid].dataPtr(), dimlist(fbox), dimlist(freg),
#ifdef TERRAIN
		  sigma[lto][igrid].dataPtr(0),
		  sigma[lto][igrid].dataPtr(1),
#  if (BL_SPACEDIM == 3)
		  sigma[lto][igrid].dataPtr(2),
#  endif
#elif (! defined SIGMA_NODE)
		  sigma_nd[0][lto][igrid].dataPtr(),
		  sigma_nd[1][lto][igrid].dataPtr(),
#  if (BL_SPACEDIM == 3)
		  sigma_nd[2][lto][igrid].dataPtr(),
#  endif
#else
		  sigma_node[lto][igrid].dataPtr(),
#endif
                  dimlist(sigbox),
		  corr[lfrom][igrid].dataPtr(), dimlist(cbox), dimlist(creg),
		  D_DECL(rat[0], rat[1], rat[2]));
#endif
    }
  }
}

#include "hg_multi.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGFRES     hgfres_
#  define   FORT_HGERES     hgeres_
#  define   FORT_HGCRES     hgcres_
#  define   FORT_HGORES     hgores_
#  define   FORT_HGIRES     hgires_
#  define   FORT_HGDRES     hgdres_
#else
#  define   FORT_HGFRES     HGFRES
#  define   FORT_HGERES     HGERES
#  define   FORT_HGCRES     HGCRES
#  define   FORT_HGORES     HGORES
#  define   FORT_HGIRES     HGIRES
#  define   FORT_HGDRES     HGDRES
#endif

extern "C" {

#if (BL_SPACEDIM == 1)
  ERROR, not relevant
#elif (BL_SPACEDIM == 2)
#  ifdef CROSS_STENCIL
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   int&, int&);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   const int*);
#  elif (defined TERRAIN)
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   int&, int&);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   const int*);
#  else
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&,
		   const int&, int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, int&, const int&, const int&
#    endif
		   );
  void FORT_HGORES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&,
		   const int&, int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, int&, const int&
#    endif
		   );
  void FORT_HGIRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&,
		   const int&, int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, int&, const int&
#    endif
		   );
  void FORT_HGDRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&,
		   const int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, const int&
#    endif
		   );
#  endif
#elif (BL_SPACEDIM == 3)
#  if (defined TERRAIN)
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   int&, int&);
  void FORT_HGERES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   const int*, const int*);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   const int*);
#  else
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   int&, int&);
  void FORT_HGERES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   const int*, const int*);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   const int*);
#  endif
#endif
}

void holy_grail_amr_multigrid::alloc_sync_caches()
{
  if (lev_min < lev_max) {
    fres_fbox  = new Box*[lev_max+1];
    fres_cbox  = new Box*[lev_max+1];
    fres_creg  = new Box*[lev_max+1];
#ifndef CONSTANT
    fres_sfbox = new Box*[lev_max+1];
    fres_scbox = new Box*[lev_max+1];
    fres_sc = new PArray<Fab>[lev_max+1];
#endif
    fres_dc = new PArray<Fab>[lev_max+1];
    fres_flag = new int*[lev_max+1];
#if (BL_SPACEDIM == 3)
    eres_fbox  = new Box*[lev_max+1];
    eres_cbox  = new Box*[lev_max+1];
    eres_creg  = new Box*[lev_max+1];
#  ifndef CONSTANT
    eres_sfbox = new Box*[lev_max+1];
    eres_scbox = new Box*[lev_max+1];
    eres_sf = new PArray<Fab>[lev_max+1];
    eres_sc = new PArray<Fab>[lev_max+1];
#  endif
    eres_df = new PArray<Fab>[lev_max+1];
    eres_dc = new PArray<Fab>[lev_max+1];
    eres_flag = new int*[lev_max+1];
#endif
    cres_fbox  = new Box*[lev_max+1];
    cres_cbox  = new Box*[lev_max+1];
    cres_creg  = new Box*[lev_max+1];
#ifndef CONSTANT
    cres_sfbox = new Box*[lev_max+1];
    cres_scbox = new Box*[lev_max+1];
    cres_sf = new PArray<Fab>[lev_max+1];
    cres_sc = new PArray<Fab>[lev_max+1];
#endif
    cres_df = new PArray<Fab>[lev_max+1];
    cres_dc = new PArray<Fab>[lev_max+1];
    cres_flag = new int*[lev_max+1];
  }

  for (int lev = lev_min + 1; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    fres_fbox[lev]  = new Box[interface[mglev].nfaces()];
    fres_cbox[lev]  = new Box[interface[mglev].nfaces()];
    fres_creg[lev]  = new Box[interface[mglev].nfaces()];
#ifndef CONSTANT
    fres_sfbox[lev] = new Box[interface[mglev].nfaces()];
    fres_scbox[lev] = new Box[interface[mglev].nfaces()];
    fres_sc[lev].resize(interface[mglev].nfaces());
#endif
    fres_dc[lev].resize(interface[mglev].nfaces());
    fres_flag[lev] = new int[interface[mglev].nfaces()];
#if (BL_SPACEDIM == 3)
    eres_fbox[lev]  = new Box[interface[mglev].nedges()];
    eres_cbox[lev]  = new Box[interface[mglev].nedges()];
    eres_creg[lev]  = new Box[interface[mglev].nedges()];
#  ifndef CONSTANT
    eres_sfbox[lev] = new Box[interface[mglev].nedges()];
    eres_scbox[lev] = new Box[interface[mglev].nedges()];
    eres_sf[lev].resize(interface[mglev].nedges());
    eres_sc[lev].resize(interface[mglev].nedges());
#  endif
    eres_df[lev].resize(interface[mglev].nedges());
    eres_dc[lev].resize(interface[mglev].nedges());
    eres_flag[lev] = new int[interface[mglev].nedges()];
#endif
    cres_fbox[lev]  = new Box[interface[mglev].ncorners()];
    cres_cbox[lev]  = new Box[interface[mglev].ncorners()];
    cres_creg[lev]  = new Box[interface[mglev].ncorners()];
#ifndef CONSTANT
    cres_sfbox[lev] = new Box[interface[mglev].ncorners()];
    cres_scbox[lev] = new Box[interface[mglev].ncorners()];
    cres_sf[lev].resize(interface[mglev].ncorners());
    cres_sc[lev].resize(interface[mglev].ncorners());
#endif
    cres_df[lev].resize(interface[mglev].ncorners());
    cres_dc[lev].resize(interface[mglev].ncorners());
    cres_flag[lev] = new int[interface[mglev].ncorners()];
    build_sync_cache(mglev, lev);
  }
}

void holy_grail_amr_multigrid::delete_sync_caches()
{
  int i;
  for (int lev = lev_min + 1; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    delete [] fres_fbox[lev];
    delete [] fres_cbox[lev];
    delete [] fres_creg[lev];
#ifndef CONSTANT
    delete [] fres_sfbox[lev];
    delete [] fres_scbox[lev];
    for (i = 0; i < interface[mglev].nfaces(); i++) {
      if (fres_flag[lev][i] == 0) delete fres_sc[lev].remove(i);
    }
#endif
    for (i = 0; i < interface[mglev].nfaces(); i++) {
      if (fres_flag[lev][i] == 0) delete fres_dc[lev].remove(i);
    }
    delete [] fres_flag[lev];
#if (BL_SPACEDIM == 3)
    delete [] eres_fbox[lev];
    delete [] eres_cbox[lev];
    delete [] eres_creg[lev];
#  ifndef CONSTANT
    delete [] eres_sfbox[lev];
    delete [] eres_scbox[lev];
    for (i = 0; i < interface[mglev].nedges(); i++) {
      if (eres_sf[lev].defined(i)) delete eres_sf[lev].remove(i);
      if (eres_flag[lev][i] == 0)  delete eres_sc[lev].remove(i);
    }
#  endif
    for (i = 0; i < interface[mglev].nedges(); i++) {
      if (eres_df[lev].defined(i)) delete eres_df[lev].remove(i);
      if (eres_flag[lev][i] == 0)  delete eres_dc[lev].remove(i);
    }
    delete [] eres_flag[lev];
#endif
    delete [] cres_fbox[lev];
    delete [] cres_cbox[lev];
    delete [] cres_creg[lev];
#ifndef CONSTANT
    delete [] cres_sfbox[lev];
    delete [] cres_scbox[lev];
    for (i = 0; i < interface[mglev].ncorners(); i++) {
      if (cres_sf[lev].defined(i)) delete cres_sf[lev].remove(i);
      if (cres_flag[lev][i] == 0)  delete cres_sc[lev].remove(i);
    }
#endif
    for (i = 0; i < interface[mglev].ncorners(); i++) {
      if (cres_df[lev].defined(i)) delete cres_df[lev].remove(i);
      if (cres_flag[lev][i] == 0)  delete cres_dc[lev].remove(i);
    }
    delete [] cres_flag[lev];
  }
  if (lev_min < lev_max) {
    delete [] fres_fbox;
    delete [] fres_cbox;
    delete [] fres_creg;
#ifndef CONSTANT
    delete [] fres_sfbox;
    delete [] fres_scbox;
    delete [] fres_sc;
#endif
    delete [] fres_dc;
    delete [] fres_flag;
#if (BL_SPACEDIM == 3)
    delete [] eres_fbox;
    delete [] eres_cbox;
    delete [] eres_creg;
#  ifndef CONSTANT
    delete [] eres_sfbox;
    delete [] eres_scbox;
    delete [] eres_sf;
    delete [] eres_sc;
#  endif
    delete [] eres_df;
    delete [] eres_dc;
    delete [] eres_flag;
#endif
    delete [] cres_fbox;
    delete [] cres_cbox;
    delete [] cres_creg;
#ifndef CONSTANT
    delete [] cres_sfbox;
    delete [] cres_scbox;
    delete [] cres_sf;
    delete [] cres_sc;
#endif
    delete [] cres_df;
    delete [] cres_dc;
    delete [] cres_flag;
  }
}

void holy_grail_amr_multigrid::build_sync_cache(int mglev, int lev)
{
  int igrid, i;

  DECLARE_GEOMETRY_TYPES;

  const IntVect& rat = gen_ratio[lev-1];
  int mglevc = ml_index[lev-1];

#ifdef TERRAIN
  int ncomp = 2 * BL_SPACEDIM - 1;
  amr_boundary bndry = boundary.terrain_sigma();
#else
  int ncomp = 1;
  amr_boundary bndry = boundary.scalar();
#endif

  for (int iface = 0; iface < interface[mglev].nfaces(); iface++) {
    // find a fine grid touching this face
    igrid = interface[mglev].fgrid(iface, 0);
    if (igrid < 0)
      igrid = interface[mglev].fgrid(iface, 1);
    unsigned geo = interface[mglev].fgeo(iface);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == ALL || igrid < 0 || interface[mglev].fflag(iface) == 1) {
      fres_flag[lev][iface] = 1; // means "fres_dc[lev][iface] not allocated"
      continue;
    }
    // fine grid on just one side
    int idim = interface[mglev].fdim(iface);
    int idir = (geo & LOW) ? -1 : 1;
    Box& fbox = fres_fbox[lev][iface];
    Box& cbox = fres_cbox[lev][iface];
    Box& creg = fres_creg[lev][iface];
    fbox = dest[lev][igrid].box();
    cbox = interface[mglev].node_face(iface);
    cbox.coarsen(rat);
    if (idir > 0)
      cbox.growLo(idim, 1);
    else
      cbox.growHi(idim, 1);
    int cgrid = find_patch(cbox, dest[lev-1]);
#ifndef CONSTANT
    Box& sigmafbox = fres_sfbox[lev][iface];
    Box& sigmacbox = fres_scbox[lev][iface];
    sigmafbox = sigma[mglev][igrid].box();
    sigmacbox = cbox;
    sigmacbox.convert(cellvect);
    if (cgrid >= 0) {
      fres_sc[lev].set(iface, &sigma[mglevc][cgrid]);
      sigmacbox = fres_sc[lev][iface].box();
    }
    else {
      fres_sc[lev].set(iface, new Fab(sigmacbox, ncomp));
      fill_patch(fres_sc[lev][iface], sigma[mglevc],
		 interface[mglevc], bndry);
    }
#endif
    if (cgrid >= 0) {
      fres_dc[lev].set(iface, &dest[lev-1][cgrid]);
      cbox = fres_dc[lev][iface].box();
      fres_flag[lev][iface] = 1;
    }
    else {
      fres_dc[lev].set(iface, new Fab(cbox));
      fres_flag[lev][iface] = 0;
    }
    IntVect t = interface[mglev].face(iface).type();
    creg = interface[mglev].node_face(iface);
    creg.coarsen(rat).grow(t - unitvect);
  }

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
    if (geo == ALL || igrid < 0 || interface[mglev].eflag(iedge) == 1) {
      eres_flag[lev][iedge] = 1; // means "eres_dc[lev][iedge] not allocated"
      continue;
    }
    Box& fbox = eres_fbox[lev][iedge];
    Box& cbox = eres_cbox[lev][iedge];
    Box& creg = eres_creg[lev][iedge];
    IntVect t = interface[mglev].edge(iedge).type();
    cbox = interface[mglev].node_edge(iedge);
    cbox.coarsen(rat).grow(t);
    fbox = refine(cbox, rat);
    eres_df[lev].set(iedge, new Fab(fbox));
    int cgrid = find_patch(cbox, dest[lev-1]);
#  ifndef CONSTANT
    Box& sigmafbox = eres_sfbox[lev][iedge];
    Box& sigmacbox = eres_scbox[lev][iedge];
    sigmafbox = fbox;
    sigmafbox.convert(cellvect);
    eres_sf[lev].set(iedge, new Fab(sigmafbox, ncomp));
    fill_patch(eres_sf[lev][iedge], sigma[mglev], interface[mglev],
	       bndry, 0, 1, iedge);
    sigmacbox = cbox;
    sigmacbox.convert(cellvect);
    if (cgrid >= 0) {
      eres_sc[lev].set(iedge, &sigma[mglevc][cgrid]);
      sigmacbox = eres_sc[lev][iedge].box();
    }
    else {
      eres_sc[lev].set(iedge, new Fab(sigmacbox, ncomp));
      fill_patch(eres_sc[lev][iedge], sigma[mglevc],
		 interface[mglevc], bndry);
    }
#  endif
    if (cgrid >= 0) {
      eres_dc[lev].set(iedge, &dest[lev-1][cgrid]);
      cbox = eres_dc[lev][iedge].box();
      eres_flag[lev][iedge] = 1;
    }
    else {
      eres_dc[lev].set(iedge, new Fab(cbox));
      eres_flag[lev][iedge] = 0;
    }
    creg = interface[mglev].node_edge(iedge);
    creg.coarsen(rat).grow(t - unitvect);
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
    if (geo == ALL || igrid < 0 || interface[mglev].cflag(icor) == 1) {
      cres_flag[lev][icor] = 1; // means "cres_dc[lev][icor] not allocated"
      continue;
    }
    Box& fbox = cres_fbox[lev][icor];
    Box& cbox = cres_cbox[lev][icor];
    Box& creg = cres_creg[lev][icor];
    cbox = interface[mglev].corner(icor);
    fbox = interface[mglev].corner(icor);
    cbox.coarsen(rat).grow(1);
    fbox.grow(rat);
    cres_df[lev].set(icor, new Fab(fbox));
    int cgrid = find_patch(cbox, dest[lev-1]);
#ifndef CONSTANT
    Box& sigmafbox = cres_sfbox[lev][icor];
    Box& sigmacbox = cres_scbox[lev][icor];
    sigmafbox = fbox;
    sigmafbox.convert(cellvect);
    cres_sf[lev].set(icor, new Fab(sigmafbox, ncomp));
    fill_patch(cres_sf[lev][icor], sigma[mglev], interface[mglev],
	       bndry, 0, 0, icor);
    sigmacbox = cbox;
    sigmacbox.convert(cellvect);
    if (cgrid >= 0) {
      cres_sc[lev].set(icor, &sigma[mglevc][cgrid]);
      sigmacbox = cres_sc[lev][icor].box();
    }
    else {
      cres_sc[lev].set(icor, new Fab(sigmacbox, ncomp));
      fill_patch(cres_sc[lev][icor], sigma[mglevc],
		 interface[mglevc], bndry);
    }
#endif
    if (cgrid >= 0) {
      cres_dc[lev].set(icor, &dest[lev-1][cgrid]);
      cbox = cres_dc[lev][icor].box();
      cres_flag[lev][icor] = 1;
    }
    else {
      cres_dc[lev].set(icor, new Fab(cbox));
      cres_flag[lev][icor] = 0;
    }
    creg = interface[mglev].corner(icor);
    creg.coarsen(rat);
  }
}

void holy_grail_amr_multigrid::interface_residual(int mglev, int lev)
{ 
  Real hx = h[mglev][0];
  Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
  Real hz = h[mglev][2];
#endif
  int igrid, i;

  DECLARE_GEOMETRY_TYPES;

  const IntVect& rat = gen_ratio[lev-1];
  int mglevc = ml_index[lev-1];

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
    const Box& fbox = fres_fbox[lev][iface];
    const Box& cbox = fres_cbox[lev][iface];
#ifndef CONSTANT
    const Box& sigmafbox = fres_sfbox[lev][iface];
    const Box& sigmacbox = fres_scbox[lev][iface];
    Fab& sigmac = fres_sc[lev][iface];
    Real *const sigmafptr = sigma[mglev][igrid].dataPtr();
#endif
    const Box& creg = fres_creg[lev][iface];
    Fab& cdst = fres_dc[lev][iface];
    if (fres_flag[lev][iface] == 0) {
      fill_patch(cdst, dest[lev-1], interface[mglevc], boundary.pressure());
    }
    Real *const rptr = resid[mglev][igrid].dataPtr();
    Real *const sptr = source[lev][igrid].dataPtr();
    Real *const dptr = dest[lev][igrid].dataPtr();
    FORT_HGFRES(rptr, dimlist(sbox),
		sptr, dimlist(sbox),
		dptr, dimlist(fbox),
		cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		sigmafptr, dimlist(sigmafbox),
		sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		dimlist(creg),
#ifdef CONSTANT
		hx,
		rat[0], idim, idir
#elif (defined TERRAIN)
		D_DECL(rat[0], rat[1], rat[2]), idim, idir
#elif (BL_SPACEDIM == 2)
		hx, hy,
#  ifdef CROSS_STENCIL
		rat[0], rat[1], idim, idir
#  else
		rat[0], rat[1], idim, idir,
		IsRZ(), mg_domain[mglevc].bigEnd(0) + 1
#  endif
#else
		hx, hy, hz,
		rat[0], rat[1], rat[2], idim, idir
#endif
		);
  }

#if (defined CROSS_STENCIL) || (defined TERRAIN)

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
    if (geo != ALL && igrid >= 0 && interface[mglev].eflag(iedge) == 0) {
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = eres_fbox[lev][iedge];
      const Box& cbox = eres_cbox[lev][iedge];
#ifndef CONSTANT
      const Box& sigmafbox = eres_sfbox[lev][iedge];
      const Box& sigmacbox = eres_scbox[lev][iedge];
      Fab& sigmaf = eres_sf[lev][iedge];
      Fab& sigmac = eres_sc[lev][iedge];
#endif
      const Box& creg = eres_creg[lev][iedge];
      IntVect t = interface[mglev].edge(iedge).type();
      Fab& fdst = eres_df[lev][iedge];
      fill_patch(fdst, dest[lev], interface[mglev], boundary.pressure(),
                 0, 1, iedge);
      Fab& cdst = eres_dc[lev][iedge];
      if (eres_flag[lev][iedge] == 0) {
	fill_patch(cdst, dest[lev-1], interface[mglevc], boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      interface[mglev].geo_array(ga, 1, iedge);
      FORT_HGERES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef CONSTANT
                  hx, rat[0],
#elif (defined TERRAIN)
		  rat[0], rat[1], rat[2],
#else
		  hx, hy, hz, rat[0], rat[1], rat[2],
#endif
		  t.getVect(), ga);
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].node_edge(iedge);
      for (i = 1; i < N_EDGE_GRIDS; i++) {
	int jgrid = interface[mglev].egrid(iedge, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
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
    if (geo != ALL && igrid >= 0 && interface[mglev].cflag(icor) == 0) {
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      Fab& sigmaf = cres_sf[lev][icor];
      Fab& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      Fab& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], interface[mglev], boundary.pressure(),
		 0, 0, icor);
      Fab& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], interface[mglevc],
		   boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      interface[mglev].geo_array(ga, 0, icor);
      FORT_HGCRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef CONSTANT
                  hx, rat[0],
#elif (defined TERRAIN)
		  D_DECL(rat[0], rat[1], rat[2]),
#else
		  D_DECL(hx, hy, hz), D_DECL(rat[0], rat[1], rat[2]),
#endif
		  ga);
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      for (i = 1; i < N_CORNER_GRIDS; i++) {
	int jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
  }

#else // 2d 9-point stencils:

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
    else if (geo == XL || geo == XH || geo == YL || geo == YH) {
      // fine grid on two adjacent sides
      int idim = (geo == XL || geo == XH) ? 0 : 1;
      int idir = (geo & LL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      Fab& sigmaf = cres_sf[lev][icor];
      Fab& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      Fab& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], interface[mglev], boundary.pressure(),
		 0, 0, icor);
      Fab& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], interface[mglevc],
		   boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGFRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef CONSTANT
                  hx,
                  rat[0], idim, idir
#else
                  hx, hy,
                  rat[0], rat[1], idim, idir,
		  IsRZ(), mg_domain[mglevc].bigEnd(0) + 1
#endif
		  );
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      for (i = 1; i < N_CORNER_GRIDS; i++) {
	int jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
    else if (geo == LL || geo == HL || geo == LH || geo == HH) {
      // outside corner
      int idir0 = (geo & XL) ? -1 : 1;
      int idir1 = (geo & YL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      const Box& cbox = cres_cbox[lev][icor];
#ifndef CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      Fab& sigmaf = cres_sf[lev][icor];
      Fab& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      Fab& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], interface[mglevc], boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      const Box& fbox = dest[lev][igrid].box();
      Real *const dptr = dest[lev][igrid].dataPtr();
      FORT_HGORES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  dptr, dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef CONSTANT
                  hx,
                  rat[0], idir0, idir1
#else
                  hx, hy,
                  rat[0], rat[1], idir0, idir1, IsRZ()
#endif
		  );
    }
    else if (geo == (LL | HH) || geo == (LH | HL)) {
      // diagonal corner
      int jdir = (geo == (LL | HH)) ? 1 : -1;
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      Fab& sigmaf = cres_sf[lev][icor];
      Fab& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      Fab& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], interface[mglev], boundary.pressure(),
		 0, 0, icor);
      Fab& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], interface[mglevc],
		   boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGDRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef CONSTANT
                  hx,
                  rat[0], jdir
#else
                  hx, hy,
                  rat[0], rat[1], jdir, IsRZ()
#endif
		  );
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      for (i = 1; i < N_CORNER_GRIDS; i++) {
	int jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
    else {
      // inside corner
      int idir0 = ((geo & XL) == XL) ? -1 : 1;
      int idir1 = ((geo & YL) == YL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      Fab& sigmaf = cres_sf[lev][icor];
      Fab& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      Fab& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], interface[mglev], boundary.pressure(),
		 0, 0, icor);
      Fab& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], interface[mglevc], boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGIRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef CONSTANT
                  hx,
                  rat[0], idir0, idir1
#else
                  hx, hy,
                  rat[0], rat[1], idir0, idir1, IsRZ()
#endif
		  );
      // fill in the grids on the other sides, if any
      const Box& freg = interface[mglev].corner(icor);
      int kgrid = -1;
      for (i = 1; i < N_CORNER_GRIDS; i++) {
	int jgrid = interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid && jgrid != kgrid) {
	  internal_copy(resid[mglev], jgrid, igrid, freg);
	  kgrid = jgrid;
	}
      }
    }
  }
#endif
}

#include "hg_multi.H"

#ifdef CONSTANT
#  define CGOPT 2
#else
#  define CGOPT 1
#endif

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGRES      hgres_
#  define   FORT_HGRESU     hgresu_
#  define   FORT_HGRLX      hgrlx_
#  define   FORT_HGRLXU     hgrlxu_
#  define   FORT_HGRLXL     hgrlxl_
#  define   FORT_HGRLNF     hgrlnf_
#  define   FORT_HGRLNB     hgrlnb_
#  define   FORT_HGCG       hgcg_
#  define   FORT_HGCG1      hgcg1_
#  define   FORT_HGCG2      hgcg2_
#  define   FORT_HGIP       hgip_
#else
#  define   FORT_HGRES      HGRES
#  define   FORT_HGRESU     HGRESU
#  define   FORT_HGRLX      HGRLX
#  define   FORT_HGRLXU     HGRLXU
#  define   FORT_HGRLXL     HGRLXL
#  define   FORT_HGRLNF     HGRLNF
#  define   FORT_HGRLNB     HGRLNB
#  define   FORT_HGCG       HGCG
#  define   FORT_HGCG1      HGCG1
#  define   FORT_HGCG2      HGCG2
#  define   FORT_HGIP       HGIP
#endif

extern "C" {

#if (BL_SPACEDIM == 1)
  ERROR, not relevant
#else
#  ifdef CONSTANT
  void FORT_HGRES(Real*, intS, Real*, Real*, intS, Real&);
  void FORT_HGRESU(Real*, intS, Real*, Real*, intS, Real&);
#    ifdef CROSS_STENCIL
  void FORT_HGRLXU(Real*, Real*, intS, Real*, intS, Real&);
#    else
  void FORT_HGRLX(Real*, Real*, intS, Real*, intS, Real&);
#    endif
#  else
#    ifdef TERRAIN
  void FORT_HGRES(Real*, intS, Real*, intS, Real*, intS,
		  Real*, intS, Real*, intS, intS);
  void FORT_HGRLX(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
  void FORT_HGRLNF(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, intS, intS, const int&, const int&);
#    elif (defined SIGMA_NODE)
  void FORT_HGRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
  void FORT_HGRESU(Real*, intS, Real*, Real*, Real*, Real*, intS);
  void FORT_HGRLX(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
  void FORT_HGRLXU(Real*, Real*, Real*, Real*, intS, Real*, intS);
  void FORT_HGRLXL(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   intS, intS, const int&);
  void FORT_HGRLNF(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, intS, intS, const int&, const int&);
#    else
  void FORT_HGRES(Real*, intS, Real*, intS, Real*, intS, RealPS, intS, intS,

		  RealRS, const int&, const int&);
  void FORT_HGRLX(Real*, intS, Real*, intS, RealPS, intS, Real*, intS, intS,
		  RealRS, const int&, const int&);
  void FORT_HGRLXL(Real*, intS, Real*, intS, RealPS, intS, Real*, intS,
                   intS, intS, RealRS, const int&, const int&, const int&);
  void FORT_HGRLNF(Real*, intS, Real*, intS, Real*, intS,
		   RealPS, intS, Real*, intS, intS, intS, RealRS,
		   const int&, const int&, const int&, const int&);
#    endif
  void FORT_HGRLNB(Real*, intS, Real*, intS,
		   intS, const int&, const int&);
#  endif

#  if (CGOPT == 1)
  void FORT_HGCG1(Real*, Real*, Real*, Real*, Real*, Real*, Real*,
		  intS, const Real&, Real&);
  void FORT_HGCG2(Real*, Real*, intS, const Real&);
  void FORT_HGIP(Real*, Real*, Real*, intS, Real&);
#  elif (CGOPT == 2)
#    if (BL_SPACEDIM == 2)
  void FORT_HGCG(Real*, Real*, Real*, Real*, Real*, Real*, Real*,
		 const int&, int*, int*,
		 int*, int*, int*, int*, int*, int*, int*,
		 const int&, Real*,
		 int*, int*, int*, int*, int*,
		 const Real&, Real&, Real&, int&, const int&);
#    else
  void FORT_HGCG(Real*, Real*, Real*, Real*, Real*, Real*, Real*,
		 const int&, int*, int*, int*,
		 int*, int*, int*, int*, int*, int*, int*,
		 const int&, Real*,
		 int*, int*, int*, int*, int*, int*, int*, int*,
		 const Real&, Real&, Real&, int&, const int&);
#    endif
#  endif
#endif
}

void holy_grail_amr_multigrid::level_residual(MultiFab& r,
					      MultiFab& s,
					      MultiFab& d,
					      copy_cache* dbc,
					      int mglev,
					      int iclear)
{
  assert(r.boxArray() == s.boxArray());
  assert(r.boxArray() == d.boxArray());

  int igrid;

  fill_borders(d, dbc, interface[mglev], mg_boundary);

#ifdef TERRAIN

  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    const Box& rbox = r[igrid].box();
    const Box& sbox = s[igrid].box();
    const Box& dbox = d[igrid].box();
    const Box& cenbox = cen[mglev][igrid].box();
    const Box& sigbox = sigma[mglev][igrid].box();
    const Box& freg = interface[mglev].part_fine(igrid);
    FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
               s[igrid].dataPtr(), dimlist(sbox),
               d[igrid].dataPtr(), dimlist(dbox),
               sigma[mglev][igrid].dataPtr(), dimlist(sigbox),
               cen[mglev][igrid].dataPtr(), dimlist(cenbox),
               dimlist(freg));
  }

  if (iclear) {
    clear_part_interface(r, interface[mglev]);
  }

#elif (defined SIGMA_NODE)

  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    const Box& rbox = r[igrid].box();
    const Box& freg = interface[mglev].part_fine(igrid);
    FORT_HGRESU(r[igrid].dataPtr(), dimlist(rbox),
                s[igrid].dataPtr(),
                d[igrid].dataPtr(),
                sigma_node[mglev][igrid].dataPtr(),
                mask[mglev][igrid].dataPtr(),
                dimlist(freg));
  }

#else

  Real hx = h[mglev][0];
  Real hy = h[mglev][1];
#  if (BL_SPACEDIM == 3)
  Real hz = h[mglev][2];
#  endif

#  ifdef CONSTANT

  if (!iclear) {
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& rbox = r[igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
      FORT_HGRESU(r[igrid].dataPtr(), dimlist(rbox),
		  s[igrid].dataPtr(), d[igrid].dataPtr(),
		  dimlist(freg), hx);
    }
  }
  else {
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& rbox = r[igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
      FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
		 s[igrid].dataPtr(), d[igrid].dataPtr(),
		 dimlist(freg), hx);
    }
    clear_part_interface(r, interface[mglev]);
  }

#  else

  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    const Box& rbox = r[igrid].box();
    const Box& sbox = s[igrid].box();
    const Box& dbox = d[igrid].box();
    const Box& freg = interface[mglev].part_fine(igrid);
#    ifndef SIGMA_NODE
    // this branch is the only one that can be reached here
    const Box& sigbox = sigma[mglev][igrid].box();
    FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
	       s[igrid].dataPtr(), dimlist(sbox),
	       d[igrid].dataPtr(), dimlist(dbox),
	       sigma_nd[0][mglev][igrid].dataPtr(),
#      if (BL_SPACEDIM == 2)
	       sigma_nd[1][mglev][igrid].dataPtr(), dimlist(sigbox),
	       dimlist(freg), hx, hy,
	       IsRZ(), mg_domain[mglev].bigEnd(0) + 1
#      else
	       sigma_nd[1][mglev][igrid].dataPtr(),
	       sigma_nd[2][mglev][igrid].dataPtr(), dimlist(sigbox),
               dimlist(freg), hx, hy, hz
#      endif
	       );
#    else
    // this branch is unreachable
    const Box& sigbox = sigma_node[mglev][igrid].box();
    FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
               s[igrid].dataPtr(), dimlist(sbox),
               d[igrid].dataPtr(), dimlist(dbox),
               sigma_node[mglev][igrid].dataPtr(), dimlist(sigbox),
               dimlist(freg));
#    endif // SIGMA_NODE
  }

  if (iclear) {
    clear_part_interface(r, interface[mglev]);
  }

#  endif // CONSTANT
#endif // SIGMA_NODE
}

void holy_grail_amr_multigrid::relax(int mglev, int i1, int is_zero)
{
  Real hx = h[mglev][0];
  Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
  Real hz = h[mglev][2];
#endif

  DECLARE_GEOMETRY_TYPES;

  Box tdom = mg_domain[mglev];
  tdom.convert(nodevect);

  for (int icount = 0; icount < i1; icount++) {

    if (smoother_mode == 0 || smoother_mode == 1 || line_solve_dim == -1) {

      if (is_zero == 0)
	fill_borders(corr[mglev], corr_bcache[mglev],
		     interface[mglev], mg_boundary);
      else
	is_zero = 0;
      for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
	const Box& sbox = resid[mglev][igrid].box();
	const Box& freg = interface[mglev].part_fine(igrid);
	if (line_solve_dim == -1) {
	  // Gauss-Seidel section:
#ifdef TERRAIN
	  const Box& fbox = corr[mglev][igrid].box();
	  const Box& cenbox = cen[mglev][igrid].box();
	  const Box& sigbox = sigma[mglev][igrid].box();
	  FORT_HGRLX(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		     resid[mglev][igrid].dataPtr(), dimlist(sbox),
		     sigma[mglev][igrid].dataPtr(), dimlist(sigbox),
		     cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		     dimlist(freg));
#elif (defined CONSTANT)
#  ifdef CROSS_STENCIL
	  FORT_HGRLXU(corr[mglev][igrid].dataPtr(),
		      resid[mglev][igrid].dataPtr(), dimlist(sbox),
		      mask[mglev][igrid].dataPtr(),
		      dimlist(freg), hx);
#  else
	  FORT_HGRLX(corr[mglev][igrid].dataPtr(),
		     resid[mglev][igrid].dataPtr(), dimlist(sbox),
		     mask[mglev][igrid].dataPtr(),
		     dimlist(freg), hx);
#  endif
#else
#ifdef SIGMA_NODE
/*
	  const Box& fbox = corr[mglev][igrid].box();
	  const Box& cenbox = cen[mglev][igrid].box();
	  const Box& sigbox = sigma_node[mglev][igrid].box();
	  FORT_HGRLX(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		     resid[mglev][igrid].dataPtr(), dimlist(sbox),
		     sigma_node[mglev][igrid].dataPtr(), dimlist(sigbox),
		     cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		     dimlist(freg));
*/
	  FORT_HGRLXU(corr[mglev][igrid].dataPtr(),
		      resid[mglev][igrid].dataPtr(),
		      sigma_node[mglev][igrid].dataPtr(),
		      cen[mglev][igrid].dataPtr(), dimlist(sbox),
		      mask[mglev][igrid].dataPtr(),
		      dimlist(freg));
#else
	  const Box& fbox = corr[mglev][igrid].box();
	  const Box& cenbox = cen[mglev][igrid].box();
	  const Box& sigbox = sigma[mglev][igrid].box();
	  FORT_HGRLX(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		     resid[mglev][igrid].dataPtr(), dimlist(sbox),
		     sigma_nd[0][mglev][igrid].dataPtr(),
#  if (BL_SPACEDIM == 2)
		     sigma_nd[1][mglev][igrid].dataPtr(), dimlist(sigbox),
		     cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		     dimlist(freg), hx, hy,
		     IsRZ(), mg_domain[mglev].bigEnd(0) + 1
#  else
		     sigma_nd[1][mglev][igrid].dataPtr(),
		     sigma_nd[2][mglev][igrid].dataPtr(), dimlist(sigbox),
		     cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		     dimlist(freg), hx, hy, hz
#  endif
		     );
#endif
#endif
	}
	else {
	  // Grid-by-grid line solve section:
#ifdef TERRAIN
	  BoxLib::Error("Terrain line solves not implemented");
#elif (defined CONSTANT)
	  BoxLib::Error("Constant-coefficient line solves not implemented");
#else
	  const Box& fbox = corr[mglev][igrid].box();
	  const Box& cenbox = cen[mglev][igrid].box();
#  ifdef SIGMA_NODE
	  const Box& sigbox = sigma_node[mglev][igrid].box();
	  FORT_HGRLXL(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		      resid[mglev][igrid].dataPtr(), dimlist(sbox),
		      sigma_node[mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), line_solve_dim);
#  else
	  const Box& sigbox = sigma[mglev][igrid].box();
	  FORT_HGRLXL(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		      resid[mglev][igrid].dataPtr(), dimlist(sbox),
		      sigma_nd[0][mglev][igrid].dataPtr(),
#    if (BL_SPACEDIM == 2)
		      sigma_nd[1][mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), hx, hy,
		      IsRZ(), mg_domain[mglev].bigEnd(0) + 1, line_solve_dim
#    else
		      sigma_nd[1][mglev][igrid].dataPtr(),
		      sigma_nd[2][mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), hx, hy, hz
#    endif
		      );
#  endif
#endif
	}
      }
      sync_borders(corr[mglev], corr_scache[mglev],
		   interface[mglev], mg_boundary);
    }
    else {
      // Full-level line solve section:
      if (line_order.length() == 0) {
	build_line_order(line_solve_dim);
      }
      int lev = lev_min, i;
      while (ml_index[lev] < mglev)
	lev++;

      for (int ipass = 0; ipass <= 1; ipass++) {
	if (is_zero == 0)
	  fill_borders(corr[mglev], corr_bcache[mglev],
		       interface[mglev], mg_boundary);
	else
	  is_zero = 0;

	// Forward solve:
	for (i = 0; i < mg_mesh[mglev].length(); i++) {

	  // Do grids in order along line_solve_dim:
	  int igrid = line_order[lev][i];
	  const Box& fbox = corr[mglev][igrid].box();
	  const Box& sbox = resid[mglev][igrid].box();
	  const Box& wbox = work[mglev][igrid].box();
	  const Box& cenbox = cen[mglev][igrid].box();
	  const Box& freg = corr[mglev].box(igrid);
#ifdef TERRAIN
	  const Box& sigbox = sigma[mglev][igrid].box();
	  FORT_HGRLNF(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		      resid[mglev][igrid].dataPtr(), dimlist(sbox),
		      work[mglev][igrid].dataPtr(), dimlist(wbox),
		      sigma[mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), line_solve_dim, ipass);
#elif (defined CONSTANT)
	  BoxLib::Error("Constant-coefficient line solves not implemented");
#else
#  ifdef SIGMA_NODE
	  const Box& sigbox = sigma_node[mglev][igrid].box();
	  FORT_HGRLNF(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		      resid[mglev][igrid].dataPtr(), dimlist(sbox),
		      work[mglev][igrid].dataPtr(), dimlist(wbox),
		      sigma_node[mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), line_solve_dim, ipass);
#  else
	  const Box& sigbox = sigma[mglev][igrid].box();
	  FORT_HGRLNF(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		      resid[mglev][igrid].dataPtr(), dimlist(sbox),
		      work[mglev][igrid].dataPtr(), dimlist(wbox),
		      sigma_nd[0][mglev][igrid].dataPtr(),
#    if (BL_SPACEDIM == 2)
		      sigma_nd[1][mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), hx, hy,
		      IsRZ(), mg_domain[mglev].bigEnd(0) + 1,
		      line_solve_dim, ipass
#    else
		      sigma_nd[1][mglev][igrid].dataPtr(),
		      sigma_nd[2][mglev][igrid].dataPtr(), dimlist(sigbox),
		      cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		      dimlist(freg), dimlist(tdom), hx, hy, hz
#    endif
		      );
#  endif
#endif

	  // Copy work arrays to following grids:
	  ListIterator<int> j(line_after[lev][igrid]);
	  for ( ; j; j++) {
	    Box b = (freg & corr[mglev].box(j()));
	    internal_copy(corr[mglev], j(), igrid, b);
	    internal_copy(work[mglev], j(), igrid, b);
	  }
	}

	// Back substitution:
	for (i = mg_mesh[mglev].length() - 1; i >= 0; i--) {

	  // Do grids in reverse order along line_solve_dim:
	  int igrid = line_order[lev][i];
	  const Box& freg = corr[mglev].box(igrid);

	  // Copy solution array from following grids:
	  ListIterator<int> j(line_after[lev][igrid]);
	  for ( ; j; j++) {
	    Box b = (freg & corr[mglev].box(j()));
	    internal_copy(corr[mglev], igrid, j(), b);
	  }

	  const Box& fbox = corr[mglev][igrid].box();
	  const Box& wbox = work[mglev][igrid].box();
	  FORT_HGRLNB(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		      work[mglev][igrid].dataPtr(), dimlist(wbox),
		      dimlist(freg), line_solve_dim, ipass);
	}
      }
    }
  }
}

void holy_grail_amr_multigrid::build_line_order(int lsd)
{
  line_order.resize(lev_max + 1);
  line_after.resize(lev_max + 1);

  for (int lev = lev_min; lev <= lev_max; lev++) {
    int igrid, i, mglev = ml_index[lev], ngrids = mg_mesh[mglev].length();

    line_order[lev].resize(ngrids);
    line_after[lev].resize(ngrids);

    for (igrid = 0; igrid < ngrids; igrid++) {
      line_order[lev].set(igrid, igrid);

      // bubble sort, replace with something faster if necessary:
      for (i = igrid; i > 0; i--) {
	if (ml_mesh[lev][line_order[lev][i]].smallEnd(lsd) <
	    ml_mesh[lev][line_order[lev][i-1]].smallEnd(lsd)) {
	  int tmp              = line_order[lev][i-1];
	  line_order[lev][i-1] = line_order[lev][i];
	  line_order[lev][i]   = tmp;
	}
	else {
	  break;
	}
      }

      for (i = 0; i < ngrids; i++) {
	if (bdryLo(ml_mesh[lev][i], lsd).intersects
	      (bdryHi(ml_mesh[lev][igrid], lsd))) {
	  line_after[lev][igrid].append(i);
	}
      }
    }
/*
    for (igrid = 0; igrid < ngrids; igrid++) {
      cout << line_order[lev][igrid] << "    ";
      ListIterator<int> j(line_after[lev][igrid]);
      for ( ; j; j++) {
	cout << " " << j();
      }
      cout << endl;
    }
*/
  }
}

void holy_grail_amr_multigrid::cgsolve(int mglev)
{
  assert(mglev == 0);

  MultiFab& r = cgwork[0];
  MultiFab& p = cgwork[1];
  MultiFab& z = cgwork[2];
  MultiFab& x = cgwork[3];
  MultiFab& w = cgwork[4];
  MultiFab& c = cgwork[5];
  MultiFab& zero_array = cgwork[6];
  MultiFab& ipmask = cgwork[7];

  unroll_cache& ruc = *cgw_ucache[0];
  unroll_cache& puc = *cgw_ucache[1];
  unroll_cache& zuc = *cgw_ucache[2];
  unroll_cache& xuc = *cgw_ucache[3];
  unroll_cache& wuc = *cgw_ucache[4];
  unroll_cache& cuc = *cgw_ucache[5];
  unroll_cache& muc = *cgw_ucache[7];

  copy_cache& pbc = *cgw1_bcache;

  Real alpha, rho;
  int i = 0, igrid;

  // x (corr[0]) should be all 0.0 at this point
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    r[igrid].copy(resid[mglev][igrid]);
    r[igrid].negate();
  }
  //r.copy(resid[mglev]);
  //r.negate();

  if (singular) {
    // singular systems are very sensitive to solvability
    w.setVal(1.0);
    alpha = inner_product(r, w) / mg_domain[mglev].volume();
    r.plus(-alpha, 0);
  }

#if (CGOPT == 2)
  FORT_HGCG(ruc.ptr, puc.ptr,
	    zuc.ptr, xuc.ptr,
	    wuc.ptr, cuc.ptr,
	    muc.ptr, mg_mesh[0].length(),
#  if (BL_SPACEDIM == 2)
	    ruc.strid, ruc.nvals,
#  else
	    ruc.strid1, ruc.strid2,
	    ruc.nvals,
#  endif
	    ruc.start, puc.start,
	    zuc.start, xuc.start,
	    wuc.start, cuc.start,
	    muc.start,
	    pbc.nsets, pbc.dptr,
#  if (BL_SPACEDIM == 2)
	    pbc.nvals,
	    pbc.dstart, pbc.sstart,
	    pbc.dstrid, pbc.sstrid,
#  else
	    pbc.nvals1, pbc.nvals2,
	    pbc.dstart, pbc.sstart,
	    pbc.dstrid1, pbc.dstrid2,
	    pbc.sstrid1, pbc.sstrid2,
#  endif
	    h[0][0], alpha, rho, i, pcode);
#elif (CGOPT == 1)
  rho = 0.0;
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    z[igrid].copy(r[igrid]);
    z[igrid].mult(c[igrid]);
    const Box& reg = p[igrid].box();
    FORT_HGIP(z[igrid].dataPtr(), r[igrid].dataPtr(),
	      ipmask[igrid].dataPtr(),
	      dimlist(reg), rho);
    p[igrid].copy(z[igrid]);
  }
  Real tol = 1.e-3 * rho;

  while (tol > 0.0) {
    i++;
    if (i > 250 && pcode >= 2)
      cout << "Conjugate-gradient iteration failed to converge" << endl;
    Real rho_old = rho;
    // safe to set the clear flag to 0 here---bogus values make it
    // into r but are cleared from z by the mask in c
    level_residual(w, zero_array, p, &pbc, 0, 0);
    alpha = 0.0;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& reg = p[igrid].box();
      FORT_HGIP(p[igrid].dataPtr(), w[igrid].dataPtr(),
		ipmask[igrid].dataPtr(),
		dimlist(reg), alpha);
    }
    alpha = rho / alpha;
    rho = 0.0;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& reg = p[igrid].box();
      FORT_HGCG1(r[igrid].dataPtr(), p[igrid].dataPtr(),
		 z[igrid].dataPtr(), x[igrid].dataPtr(),
		 w[igrid].dataPtr(), c[igrid].dataPtr(),
		 ipmask[igrid].dataPtr(), dimlist(reg), alpha, rho);
    }
    if (pcode >= 3)
      cout << i << " " << rho << endl;
    if (rho <= tol || i > 250)
      break;
    alpha = rho / rho_old;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& reg = p[igrid].box();
      FORT_HGCG2(p[igrid].dataPtr(), z[igrid].dataPtr(),
		 dimlist(reg), alpha);
    }
  }
#else
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    z[igrid].copy(r[igrid]);
    z[igrid].mult(c[igrid]);
  }
  //z.assign(r).mult(c);
  rho = inner_product(z, r);
  Real tol = 1.e-3 * rho;
  //p.assign(0.0);
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    p[igrid].copy(z[igrid]);
  }

  while (tol > 0.0) {
    i++;
    if (i > 250 && pcode >= 2)
      cout << "Conjugate-gradient iteration failed to converge" << endl;
    Real rho_old = rho;
    // safe to set the clear flag to 0 here---bogus values make it
    // into r but are cleared from z by the mask in c
    level_residual(w, zero_array, p, &pbc, 0, 0);
    alpha = rho / inner_product(p, w);
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      w[igrid].mult(alpha);
      r[igrid].minus(w[igrid]);
      w[igrid].copy(p[igrid]);
      w[igrid].mult(alpha);
      x[igrid].plus(w[igrid]);
      z[igrid].copy(r[igrid]);
      z[igrid].mult(c[igrid]);
    }
    //r.minus(w.mult(alpha));
    //x.plus(w.assign(p).mult(alpha));
    //z.assign(r).mult(c);
    rho = inner_product(z, r);
    if (pcode >= 3)
      cout << i << " " << rho << endl;
    if (rho <= tol || i > 250)
      break;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      p[igrid].mult(rho / rho_old);
      p[igrid].plus(z[igrid]);
    }
    //p.mult(rho / rho_old).plus(z);
  }
#endif

  if (pcode >= 2)
    cout << i << " iterations required for conjugate-gradient" << endl;
}
