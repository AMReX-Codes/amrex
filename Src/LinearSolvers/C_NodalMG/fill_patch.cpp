
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_FIPRODC   iprodc_
#  define FORT_FIPRODN   iprodn_
#  define FORT_FFCPYU    fcpyu_
#  define FORT_FFCPY2    fcpy2_
#else
#  define FORT_FIPRODC   IPRODC
#  define FORT_FIPRODN   IPRODN
#  define FORT_FFCPYU    FCPYU
#  define FORT_FFCPY2    FCPY2
#endif

extern "C"
{
  void FORT_FIPRODC(const Real*, intS, const Real*, intS, intS, Real&);
  void FORT_FIPRODN(const Real*, intS, const Real*, intS, intS, Real&);
  void FORT_FFCPYU(Real*, Real*, intS, const int&);
#if (BL_SPACEDIM == 2)
  void FORT_FFCPY2(Real*, intS, Real*, intS, intS, const int&, const int&);
#else
  void FORT_FFCPY2(Real*, intS, Real*, intS, intS, const int&, const int&, const int&);
#endif
}

Real inner_product(const MultiFab& r, const MultiFab& s)
{
  assert(r.ok() && s.ok());
  assert(r.nComp() == 1);
  assert(s.nComp() == 1);
  assert(type(r) == type(s));

  Real sum = 0.0;

  if (type(r) == IntVect::TheCellVector()) {
    // for (igrid = 0; igrid < r.length(); igrid++) {
    for ( ConstMultiFabIterator rcmfi(r); rcmfi.isValid(); ++rcmfi) {
	ConstDependentMultiFabIterator scmfi(rcmfi, s);
      const Box& rbox = rcmfi->box();
      const Box& sbox = scmfi->box();
      const Box& reg  = rcmfi.validbox();
      FORT_FIPRODC(rcmfi->dataPtr(), DIMLIST(rbox),
	      scmfi->dataPtr(), DIMLIST(sbox),
	      DIMLIST(reg), sum);
    }
  }
  else if (type(r) == IntVect::TheNodeVector()) {
    // for (igrid = 0; igrid < r.length(); igrid++) {
    for ( ConstMultiFabIterator rcmfi(r); rcmfi.isValid(); ++rcmfi) {
	ConstDependentMultiFabIterator scmfi(rcmfi, s);
      const Box& rbox = rcmfi->box();
      const Box& sbox = scmfi->box();
      const Box& reg  = rcmfi.validbox();
      FORT_FIPRODN(rcmfi->dataPtr(), DIMLIST(rbox),
	      scmfi->dataPtr(), DIMLIST(sbox),
	      DIMLIST(reg), sum);
    }
  }
  else {
    BoxLib::Error("inner_product---only supported for CELL- or NODE-based data");
  }
    ParallelDescriptor::ReduceRealSum(sum);
  return sum;
}

/*
const MultiFab&
  initialize(Real (*f)(const Intvect&,const Intvect&,int,int))
{
  for (int i = 0; i < mesh().ngrids(); i++) {
    grid(i).initialize(f, mesh().sig());
  }
  return *this;
}
*/

#if 0
// Begin fillpatch stuff, still in unfinished state.
// Significant optimizations possible: avoid copying patches
// whenever an existing one will do.
// All Boxes and data objects must have same index type.

// Chops dest according to an intersection with source.  The returned
// box is the half completely outside of source, while dest is modified
// to contain the half that intersects source.  It is assumed that such a
// nontrivial chop is possible, or this routine would not have been called.

static Box box_chop(Box& dest, const Box& source)
{
  int i;
  for (i = BL_SPACEDIM - 1; i >= 0; i--) {
    if (dest.type(i) == IndexType::CELL) {
      if (dest.smallEnd(i) <= source.bigEnd(i) &&
	  source.bigEnd(i) < dest.bigEnd(i)) {
	return dest.chop(i,source.bigEnd(i)+1);
      }
      if (dest.smallEnd(i) < source.smallEnd(i) &&
	  source.smallEnd(i) <= dest.bigEnd(i)) {
	Box tmp(dest);
	dest = tmp.chop(i,source.smallEnd(i));
	return tmp;
      }
    }
    else {
      if (dest.smallEnd(i) < source.bigEnd(i) &&
	  source.bigEnd(i) < dest.bigEnd(i)) {
	return dest.chop(i,source.bigEnd(i)).growLo(i, -1);
      }
      if (dest.smallEnd(i) < source.smallEnd(i) &&
	  source.smallEnd(i) < dest.bigEnd(i)) {
	Box tmp(dest);
	dest = tmp.chop(i,source.smallEnd(i));
	return tmp.growHi(i, -1);
      }
    }
  }
  // this section only reached for node box outside but touching boundary
  for (i = BL_SPACEDIM - 1; i >= 0; i--) {
    if (dest.bigEnd(i) == source.smallEnd(i) &&
	dest.smallEnd(i) < source.smallEnd(i)) {
      Box tmp(dest);
      dest.setSmall(i, dest.bigEnd(i));
      return tmp.growHi(i, -1);
    }
    if (dest.smallEnd(i) == source.bigEnd(i) &&
	dest.bigEnd(i) > source.bigEnd(i)) {
      Box tmp(dest);
      dest.setBig(i, dest.smallEnd(i));
      return tmp.growLo(i, -1);
    }
  }
  assert(0);
  return Box();
}

static int best_match(MultiFab& r, const Box& region, int& igrid, int bord)
{
  int overlap = 0;
  if (bord == r.nGrow()) {
    for (int i = 0; i < r.length(); i++) {
      if (region.intersects(r[i].box())) {
	int overlap1 = (region & r[i].box()).numPts();
	if (overlap1 > overlap) {
	  igrid = i;
	  overlap = overlap1;
	}
      }
    }
  }
  else {
    for (int i = 0; i < r.length(); i++) {
      Box tb = grow(r[i].box(), bord - r.nGrow());
      if (region.intersects(tb)) {
	int overlap1 = (region & tb).numPts();
	if (overlap1 > overlap) {
	  igrid = i;
	  overlap = overlap1;
	}
      }
    }
  }
  return (overlap > 0) ? (overlap == region.numPts() ? 1 : 2) : 0;
}
#endif

/*
grid_real get_patch(const Box& region,
		    const level_interface& lev_interface,
		    const amr_boundary_class& bdy, int flags)
{
  int igrid;
  if (border() == 0 || (flags & 2)) {
    for (igrid = 0; igrid < mesh().ngrids(); igrid++) {
      if (box(igrid).contains(region)) {
	if (ncomp == fgrid(igrid).nComp())
	  return fgrid(igrid);
	else
	  return grid(igrid);
      }
    }
  }
  else {
    for (igrid = 0; igrid < mesh().ngrids(); igrid++) {
      Box tb = grow(box(igrid), -border());
      if (tb.contains(region)) {
	if (ncomp == fgrid(igrid).nComp())
	  return fgrid(igrid);
	else
	  return grid(igrid);
      }
    }
  }

  grid_real retgr(region, nComp());
  fill_patch(retgr, region, lev_interface, bdy, flags);
  return retgr;
}

int get_patch(FArrayBox& patch, const Box& region,
	      const level_interface& lev_interface,
	      const amr_boundary_class& bdy, int flags)
{
  if (flags & 8) {
    for (int iqq = 0; iqq < flags/16; iqq++)
      cout << "  ";
    cout << "Getting " << region << endl;
    flags += 16;
  }

  int igrid;
  if (border() == 0 || (flags & 2)) {
    for (igrid = 0; igrid < mesh().ngrids(); igrid++) {
      if (box(igrid).contains(region)) {
	if (ncomp == fgrid(igrid).nComp())
	  patch.alias(fgrid(igrid));
	else
	  patch.alias(grid(igrid));
	return 1;
      }
    }
  }
  else {
    for (igrid = 0; igrid < mesh().ngrids(); igrid++) {
      Box tb = grow(box(igrid), -border());
      if (tb.contains(region)) {
	if (ncomp == fgrid(igrid).nComp())
	  patch.alias(fgrid(igrid));
	else
	  patch.alias(grid(igrid));
	return 1;
      }
    }
  }

  patch.alloc(region, nComp());
  return fill_patch(patch, region, lev_interface, bdy, flags);
}
*/

int find_patch(const Box& region, const MultiFab& r, int flags)
{
  int igrid;
  if (r.nGrow() == 0 || (flags & 2)) {
    for (igrid = 0; igrid < r.length(); igrid++) {
      if (r[igrid].box().contains(region))
	return igrid;
    }
  }
  else {
    for (igrid = 0; igrid < r.length(); igrid++) {
      if (r.box(igrid).contains(region))
	return igrid;
    }
  }

  return -1;
}

int fill_patch_blindly(FArrayBox& patch,
		       const Box& region,
		       const MultiFab& r,
		       int flags)
{
  int igrid;
  if (r.nGrow() == 0 || (flags & 2)) {
    for (igrid = 0; igrid < r.length(); igrid++) {
      if (r[igrid].box().contains(region)) {
	patch.copy(r[igrid], region, 0, region, 0, patch.nComp());
	return 1;
      }
    }
    for (igrid = 0; igrid < r.length(); igrid++) {
      if (r[igrid].box().intersects(region)) {
	Box tb = region & r[igrid].box();
	patch.copy(r[igrid], tb, 0, tb, 0, patch.nComp());
      }
    }
  }
  else {
    for (igrid = 0; igrid < r.length(); igrid++) {
      Box tb = grow(r[igrid].box(), -r.nGrow());
      if (tb.contains(region)) {
	patch.copy(r[igrid], region, 0, region, 0, patch.nComp());
	return 1;
      }
    }
    for (igrid = 0; igrid < r.length(); igrid++) {
      Box tb = grow(r[igrid].box(), -r.nGrow());
      if (tb.intersects(region)) {
	tb &= region;
	patch.copy(r[igrid], tb, 0, tb, 0, patch.nComp());
      }
    }
  }
  return 0;
}

int fill_exterior_patch_blindly(FArrayBox& patch,
				const Box& region,
				const MultiFab& r,
				const level_interface& lev_interface,
				const amr_boundary_class& bdy,
				int flags)
{
  const BoxArray& em = lev_interface.exterior_mesh();
  int igrid;
  for (igrid = 0; igrid < em.length(); igrid++) {
    int jgrid = lev_interface.direct_exterior_ref(igrid);
    if (jgrid >= 0) {
      Box tb;
      tb = em[igrid];
      tb.convert(type(r));
      if (r.nGrow() > 0 && (flags & 2))
	tb.grow(r.nGrow());
      if (tb.contains(region)) {
	bdy.fill(patch, region, r, jgrid, lev_interface.domain());
	return 1;
      }
      if (tb.intersects(region)) {
	tb &= region;
	bdy.fill(patch, tb, r, jgrid, lev_interface.domain());
      }
    }
  }
  return 0;
}

int fill_patch(FArrayBox& patch, const Box& region,
	       const MultiFab& r,
	       const level_interface& lev_interface,
	       const amr_boundary_class& bdy, int flags,
	       int idim, int index)
{
  if (!region.ok())
    return 1;

  if (flags & 4)
    patch.setVal(0.0, region, 0, patch.nComp());

  assert(patch.nComp() == r.nComp());
  assert(type(patch) == type(r));
  assert(lev_interface.ok());

  Box tdomain = lev_interface.domain();
  tdomain.convert(type(patch));
  Box idomain = grow(tdomain, IntVect::TheZeroVector() - type(r));

  if ((flags & 1) == 0) {
    if (idim == -1 || (flags & 2)) {
      if (idomain.contains(region) || bdy.defined() == 0) {
	return fill_patch_blindly(patch, region, r, flags);
      }
      else if (!tdomain.intersects(region)) {
	return fill_exterior_patch_blindly(patch, region, r,
					   lev_interface, bdy, flags);
      }
      else if (idomain.intersects(region)) {
	if (fill_patch_blindly(patch, region, r, flags) == 1)
	  return 1;
	else
	  return fill_exterior_patch_blindly(patch, region, r,
					     lev_interface, bdy, flags);
      }
      else {
	if (fill_exterior_patch_blindly(patch, region, r,
					lev_interface, bdy, flags) == 1)
	  return 1;
	else
	  return fill_patch_blindly(patch, region, r, flags);
      }
    }
    else if (idim == 0) {
      int gridnum[level_interface::N_CORNER_GRIDS+1];
      gridnum[0] = -1;
      for (int i = 0; i < level_interface::N_CORNER_GRIDS; i++) {
	int igrid = lev_interface.cgrid(index,i);
	if (igrid != -1) {
	  for (int j = 0; gridnum[j] != igrid; j++) {
	    if (gridnum[j] == -1) {
	      gridnum[j] = igrid;
	      gridnum[j+1] = -1;
	      if (igrid >= 0) {
		Box tb = r.box(igrid);
		tb &= region;
		const Box& rbox = r[igrid].box();
                FORT_FFCPY(patch.dataPtr(), DIMLIST(patch.box()),
                      DIMLIST(tb),
                      r[igrid].dataPtr(), DIMLIST(rbox), patch.nComp());
	      }
	      else {
		igrid = -2 - igrid;
		Box tb = lev_interface.exterior_mesh()[igrid];
		tb.convert(type(r));
		tb &= region;
		bdy.fill(patch, tb, r,
			 lev_interface.direct_exterior_ref(igrid),
			 lev_interface.domain());
	      }
	      break;
	    }
	  }
	}
      }
    }
#if (BL_SPACEDIM == 3)
    else if (idim == 1) {
      int gridnum[level_interface::N_EDGE_GRIDS+1];
      gridnum[0] = -1;
      for (int i = 0; i < level_interface::N_EDGE_GRIDS; i++) {
	int igrid = lev_interface.egrid(index,i);
	if (igrid != -1) {
	  for (int j = 0; gridnum[j] != igrid; j++) {
	    if (gridnum[j] == -1) {
	      gridnum[j] = igrid;
	      gridnum[j+1] = -1;
	      if (igrid >= 0) {
		Box tb = r.box(igrid);
		tb &= region;
		const Box& rbox = r[igrid].box();
                FORT_FFCPY(patch.dataPtr(), DIMLIST(patch.box()),
                      DIMLIST(tb),
                      r[igrid].dataPtr(), DIMLIST(rbox), patch.nComp());
	      }
	      else {
		igrid = -2 - igrid;
		Box tb = lev_interface.exterior_mesh()[igrid];
		tb.convert(type(r));
		tb &= region;
		bdy.fill(patch, tb, r,
			 lev_interface.direct_exterior_ref(igrid),
			 lev_interface.domain());
	      }
	      break;
	    }
	  }
	}
      }
    }
#endif
    else if (idim == level_interface::FACEDIM) {
      int gridnum[level_interface::N_FACE_GRIDS+1];
      gridnum[0] = -1;
      for (int i = 0; i < level_interface::N_FACE_GRIDS; i++) {
	int igrid = lev_interface.fgrid(index,i);
	if (igrid != -1) {
	  for (int j = 0; gridnum[j] != igrid; j++) {
	    if (gridnum[j] == -1) {
	      gridnum[j] = igrid;
	      gridnum[j+1] = -1;
	      if (igrid >= 0) {
		Box tb = r.box(igrid);
		tb &= region;
		const Box& rbox = r[igrid].box();
                FORT_FFCPY(patch.dataPtr(), DIMLIST(patch.box()),
                      DIMLIST(tb),
                      r[igrid].dataPtr(), DIMLIST(rbox), patch.nComp());
	      }
	      else {
		igrid = -2 - igrid;
		Box tb = lev_interface.exterior_mesh()[igrid];
		tb.convert(type(r));
		tb &= region;
		bdy.fill(patch, tb, r,
			 lev_interface.direct_exterior_ref(igrid),
			 lev_interface.domain());
	      }
	      break;
	    }
	  }
	}
      }
    }
  }
  else {
    BoxLib::Error("fill_patch---lev_interface version only defined for blind mode");
  }
  return 1;
}

/*
int fill_patch(FArrayBox& patch,
	       const Box& region, const Box& active,
	       const amr_boundary_class& bdy, int flags,
	       int igrid)
{
  cout << "Warning:  using obsolete form of fill_patch" << endl;

  if (!region.ok())
    return 1;

  // This version is particularly inefficient for handling node-based
  // boundary conditions, and sometimes gives incorrect results for
  // node-based periodic boundaries.

  if (flags & 8) {
    for (int iqq = 0; iqq < flags/16; iqq++)
      cout << "  ";
    cout << "Filling " << region << endl;
    flags += 16;
  }

  if (flags & 4)
    patch.assign(0.0, region);

  chkcomp(patch.nComp());
  if (!region.sameType(patch.box()))
    BoxLib::Error("fill_patch---incompatible patch");

  Box idomain = grow(tdomain(), IntVect::TheZeroVector - type());

  if (igrid < 0 && (flags & 1) == 0 && idomain.contains(region)) {
    fill_patch_blindly(patch, region, flags);
    return 1;
  }

  if (tdomain().contains(region)) {
    // Interior patch;
    int retval;
    if (igrid >= 0) {
      if (flags & 2)
	retval = obj().gr[igrid].box().contains(region) ? 1 : 2;
      else if (region.cellCentered())
	retval = mesh()[igrid].contains(region) ? 1 : 2;
      else
	retval = mesh().box(igrid, type()).contains(region) ? 1 : 2;
    }
    else {
      retval = best_match(region, igrid, (flags & 2) ? border() : 0);
    }

    if (retval == 0) {
      // no intersections with grids on this level
      if (!idomain.intersects(region)) {
	// last-ditch chance, maybe boundary condition can do something
	cout << "Using fill_patch_special" << endl;
	return bdy.fill_patch_special(patch, region, *this, flags);
      }
//      else if (!idomain.contains(region)) {
//	Box side1 = region;
//	Box side2 = box_chop(side1, idomain);
//	cout << side1 << " " << side2 << endl;
//	fill_patch(patch, side1, active, bdy, flags);
//	fill_patch(patch, side2, active, bdy, flags);
//	return 0;
//      }
      else {
	return 0;
      }
    }
    else if (retval == 1) {
      patch.fab().copy(fgrid(igrid).fab(), region, comp, region,
		       patch.component(), ncomp);
      return 1;
    }
    else if (retval == 2) {
      Box side1 = region, side2;
      if (flags & 2)
	side2 = obj().gr[igrid].box();
      else
	side2 = mesh().box(igrid, type());
      side2 = box_chop(side1, side2);
      int ret0 = fill_patch(patch, side1, active, bdy, flags, igrid);
      int ret1 = fill_patch(patch, side2, active, bdy, flags);
      return (ret0 && ret1);
    }
  }
  else {
    Box side2 = region;
    int ret0 = 1, ret1 = 1;
    if (idomain.intersects(region)) {
      // patch crosses boundary
      Box side1 = region;
      side2 = box_chop(side1, tdomain());
      ret0 = fill_patch(patch, side1, active, bdy, (flags | 4), igrid);
    }
    // patch outside boundary
    int bdir = bdy.dir(side2, mesh().domain());
    Box bb = bdy.box(side2, mesh().domain(), bdir);
    // tdomain().contains(bb) can be false for certain inflow bc's.
    // If this is the case, bb must be fillable from ghost cells and
    // (flags & 2) must be set or an infinite recursion will occur.
    if (active.contains(bb) && tdomain().contains(bb)) {
      if (!patch.box().contains(active))
	BoxLib::Error("fill_patch---bogus active region");
      bdy.fill(patch, side2, patch, bb, mesh().domain(), bdir);
    }
    else {
      FArrayBox gb;
      ret1 = get_patch(gb, bb, level_interface(), bdy, (flags | 5));
      bdy.fill(patch, side2, gb, bb, mesh().domain(), bdir);
    }
    if (ret1 == 0 && tdomain().intersects(region)) {
      // bc failed and region just touches boundary---split
      Box side1 = region;
      side2 = box_chop(side1, tdomain());
      ret0 = fill_patch(patch, side1, active, bdy, flags, igrid);
      if (flags & 1) {
	// redo outside to get return value right
	bb = bdy.box(side2, mesh().domain(), bdir);
	FArrayBox gb;
	ret1 = get_patch(gb, bb, level_interface(), bdy, (flags | 5));
	bdy.fill(patch, side2, gb, bb, mesh().domain(), bdir);
      }
      else {
	ret1 = 1;
      }
    }
    return (ret0 && ret1);
  }
  BoxLib::Error("fill_patch---shouldn't get here");
  return 0;
}
*/

void sync_internal_borders(MultiFab& r, const level_interface& lev_interface)
{
  int igrid, jgrid;
  if (type(r) == IntVect::TheNodeVector()) {
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) {
      igrid = lev_interface.fgrid(iface, 0);
      jgrid = lev_interface.fgrid(iface, 1);
      // only do interior faces with fine grid on both sides
      if (igrid < 0 || jgrid < 0 || lev_interface.fgeo(iface) != level_interface::ALL)
	break;
      internal_copy(r, jgrid, igrid, lev_interface.node_face(iface));
    }
#if (BL_SPACEDIM == 2)
    for (int icor = 0; icor < lev_interface.ncorners(); icor++) {
      igrid = lev_interface.cgrid(icor, 0);
      jgrid = lev_interface.cgrid(icor, 3);
      // only do interior corners with fine grid on all sides
      if (igrid < 0 || jgrid < 0 || lev_interface.cgeo(icor) != level_interface::ALL)
	break;
      if (jgrid == lev_interface.cgrid(icor, 1))
	internal_copy(r, jgrid, igrid, lev_interface.corner(icor));
    }
#else
    for (int iedge = 0; iedge < lev_interface.nedges(); iedge++) {
      igrid = lev_interface.egrid(iedge, 0);
      jgrid = lev_interface.egrid(iedge, 3);
      // only do interior edges with fine grid on all sides
      if (igrid < 0 || jgrid < 0 || lev_interface.egeo(iedge) != level_interface::ALL)
	break;
      if (jgrid == lev_interface.egrid(iedge, 1))
	internal_copy(r, jgrid, igrid, lev_interface.node_edge(iedge));
    }
    for (int icor = 0; icor < lev_interface.ncorners(); icor++) {
      igrid = lev_interface.cgrid(icor, 0);
      jgrid = lev_interface.cgrid(icor, 7);
      // only do interior corners with fine grid on all sides
      if (igrid < 0 || jgrid < 0 || lev_interface.cgeo(icor) != level_interface::ALL)
	break;
      if (lev_interface.cgrid(icor, 3) == lev_interface.cgrid(icor, 1)) {
	if (jgrid != lev_interface.cgrid(icor, 3)) {
	  internal_copy(r, jgrid, igrid, lev_interface.corner(icor));
	  jgrid = lev_interface.cgrid(icor, 5);
	  if (jgrid != lev_interface.cgrid(icor, 7))
	    internal_copy(r, jgrid, igrid, lev_interface.corner(icor));
	}
      }
      else if (lev_interface.cgrid(icor, 5) == lev_interface.cgrid(icor, 1)) {
	if (jgrid != lev_interface.cgrid(icor, 5)) {
	  internal_copy(r, jgrid, igrid, lev_interface.corner(icor));
	  jgrid = lev_interface.cgrid(icor, 3);
	  if (jgrid != lev_interface.cgrid(icor, 7)) {
	    internal_copy(r, jgrid, igrid, lev_interface.corner(icor));
	    if (jgrid == lev_interface.cgrid(icor, 2)) {
	      jgrid = lev_interface.cgrid(icor, 6);
	      if (jgrid != lev_interface.cgrid(icor, 7))
		internal_copy(r, jgrid, igrid, lev_interface.corner(icor));
	    }
	  }
	}
      }
    }
#endif
  }
  else {
    BoxLib::Error("sync_internal_borders---only NODE-based sync defined");
  }
}

// The sequencing used in fill_internal_borders, fcpy2 and set_border_cache
// (narrow x, medium y, wide z) is necessary to avoid overwrite problems
// like those seen in the sync routines.  Boundary copies are all wide
// regardless of direction and come after interior copies---overwrite
// difficulties are avoided since grids can't bridge a boundary.

// Modifications are necessary in 3D to deal with lack of diagonal
// communication across edges at the coarse-fine lev_interface.  These
// modifications take the form of narrowing certain copies to avoid
// overwriting good values with bad ones.

void fill_internal_borders(MultiFab& r, const level_interface& lev_interface,
			   int w)
{
  w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;
  int igrid, jgrid;
  if (type(r) == IntVect::TheNodeVector()) {
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) {
      igrid = lev_interface.fgrid(iface, 0);
      jgrid = lev_interface.fgrid(iface, 1);
      if (igrid < 0 || jgrid < 0 || lev_interface.fgeo(iface) != level_interface::ALL)
	break;
#if 1
      const Box& b = lev_interface.node_face(iface);
      Real *const ptra = r[igrid].dataPtr();
      Real *const ptrb = r[jgrid].dataPtr();
      const Box& boxa = r[igrid].box();
      const Box& boxb = r[jgrid].box();
#  if (BL_SPACEDIM == 2)
      FORT_FFCPY2(ptra, DIMLIST(boxa), ptrb, DIMLIST(boxb),
	     DIMLIST(b), w, r.nComp());
#  else
      const int ibord = r.nGrow();
      FORT_FFCPY2(ptra, DIMLIST(boxa), ptrb, DIMLIST(boxb),
	     DIMLIST(b), w, ibord, r.nComp());
#  endif
#else
      const int idim = lev_interface.fdim(iface);
      Box bj = lev_interface.node_face(iface);
      Box bi = lev_interface.node_face(iface);
      for (int i = 0; i < idim; i++) {
	if (r.box(jgrid).smallEnd(i) == bj.smallEnd(i))
	  bj.growLo(i, w);
	if (r.box(jgrid).bigEnd(i) == bj.bigEnd(i))
	  bj.growHi(i, w);
	if (r.box(igrid).smallEnd(i) == bi.smallEnd(i))
	  bi.growLo(i, w);
	if (r.box(igrid).bigEnd(i) == bi.bigEnd(i))
	  bi.growHi(i, w);
      }
      bj.shift(idim, -1).growLo(idim, w-1);
      bi.shift(idim,  1).growHi(idim, w-1);
      internal_copy(r, jgrid, igrid, bj);
      internal_copy(r, igrid, jgrid, bi);
#endif
    }
  }
  else if (type(r) == IntVect::TheCellVector()) {
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) {
      igrid = lev_interface.fgrid(iface, 0);
      jgrid = lev_interface.fgrid(iface, 1);
      if (igrid < 0 || jgrid < 0 || lev_interface.fgeo(iface) != level_interface::ALL)
	break;
      const int idim = lev_interface.fdim(iface);
#if (BL_SPACEDIM == 2)
      Box b = lev_interface.face(iface);
      if (idim == 1)
        b.grow(0, w);
      b.growLo(idim, w).convert(IntVect::TheCellVector());
      internal_copy(r, jgrid, igrid, b);
      internal_copy(r, igrid, jgrid, b.shift(idim, w));
#else
      Box bj = lev_interface.face(iface);
      Box bi = lev_interface.face(iface);
      for (int i = 0; i < idim; i++) {
	if (r.box(jgrid).smallEnd(i) == bj.smallEnd(i))
	  bj.growLo(i, w);
	if (r.box(jgrid).bigEnd(i) == bj.bigEnd(i))
	  bj.growHi(i, w);
	if (r.box(igrid).smallEnd(i) == bi.smallEnd(i))
	  bi.growLo(i, w);
	if (r.box(igrid).bigEnd(i) == bi.bigEnd(i))
	  bi.growHi(i, w);
      }
      bj.growLo(idim, w).convert(IntVect::TheCellVector());
      bi.growHi(idim, w).convert(IntVect::TheCellVector());
      internal_copy(r, jgrid, igrid, bj);
      internal_copy(r, igrid, jgrid, bi);
#endif
    }
  }
  else {
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) {
      igrid = lev_interface.fgrid(iface, 0);
      jgrid = lev_interface.fgrid(iface, 1);
      if (igrid < 0 || jgrid < 0 || lev_interface.fgeo(iface) != level_interface::ALL)
	break;
      const int idim = lev_interface.fdim(iface);
      const int a = (type(r, idim) == IndexType::NODE);
#if (BL_SPACEDIM == 2)
      Box b = lev_interface.face(iface);
      if (idim == 1)
        b.grow(0, w);
      b.shift(idim, -a).growLo(idim, w-a).convert(type(r));
      internal_copy(r, jgrid, igrid, b);
      internal_copy(r, igrid, jgrid, b.shift(idim, w+a));
#else
      Box bj = lev_interface.face(iface);
      Box bi = lev_interface.face(iface);
#if 0
      for (int i = 0; i < idim; i++) {
	if (mesh().box(jgrid).smallEnd(i) == bj.smallEnd(i))
	  bj.growLo(i, w);
	if (mesh().box(jgrid).bigEnd(i) == bj.bigEnd(i))
	  bj.growHi(i, w);
	if (mesh().box(igrid).smallEnd(i) == bi.smallEnd(i))
	  bi.growLo(i, w);
	if (mesh().box(igrid).bigEnd(i) == bi.bigEnd(i))
	  bi.growHi(i, w);
      }
#else
      BoxLib::Error("fill_internal_borders---check index arithmetic for mixed types in 3D");
#endif
      bj.shift(idim, -a).growLo(idim, w-a).convert(type(r));
      bi.shift(idim,  a).growHi(idim, w-a).convert(type(r));
      internal_copy(r, jgrid, igrid, bj);
      internal_copy(r, igrid, jgrid, bi);
#endif
    }
  }
}

void clear_part_interface(MultiFab& r, const level_interface& lev_interface)
{
  if (r.nComp() != 1)
    BoxLib::Error("clear_part_interface---only single components currently supported");

  int igrid;
  if (type(r) == IntVect::TheNodeVector()) {
    for (int i = 0; i < BL_SPACEDIM; i++) {
      for (int ibox = 0; ibox < lev_interface.nboxes(i); ibox++) {
	// coarse-fine face contained in part_fine grid, or orphan edge/corner
	if ((igrid = lev_interface.aux(i, ibox)) >= 0)
	  r[igrid].setVal(0.0, lev_interface.node_box(i, ibox), 0);
      }
    }
  }
  else {
    BoxLib::Error("clear_part_interface---only NODE-based version defined");
  }
}

void interpolate_patch(FArrayBox& patch, const Box& region,
		       MultiFab& r, const IntVect& rat,
		       const amr_interpolator& interp,
		       const level_interface& lev_interface,
		       const amr_boundary_class& bdy)
{
  assert(region.sameType(patch.box()));

  Box cb = interp.box(region, rat);
  int igrid = find_patch(cb, r);
  if (igrid == -1) {
    FArrayBox cgr(cb, r.nComp());
    fill_patch(cgr, cb, r, lev_interface, bdy);
    interp.fill(patch, region, cgr, cb, rat);
  }
  else {
    interp.fill(patch, region, r[igrid], cb, rat);
  }
}

void restrict_patch(FArrayBox& patch, const Box& region,
		    MultiFab& r, const IntVect& rat,
#ifdef HG_USE_CACHE
		    const copy_cache* border_cache,
#endif
		    const amr_restrictor_class& restric,
		    const level_interface& lev_interface,
		    const amr_boundary_class& bdy)
{
  assert(region.sameType(patch.box()));
  assert(region.type() == type(r));

  for (int igrid = 0; igrid < r.length(); igrid++) {
    Box cbox = r.box(igrid);
    cbox = restric.box(cbox, rat);
    if (region.intersects(cbox)) {
      cbox &= region;
      restric.fill(patch, cbox, r[igrid], rat);
    }
  }

  // Interface restriction is sufficiently rare and specialized that
  // we will let the restrictor handle it---at least for now.

  if (!lev_interface.null()) {
    // This assertion difficult in BoxLib since r.mesh() is not cc:
    //assert(r.mesh() == lev_interface.interior_mesh());
    restric.lev_interface(patch, region, r,
#ifdef HG_USE_CACHE
	border_cache, 
#endif
	lev_interface, bdy, rat);
  }
}

void restrict_level(MultiFab& dest, int bflag,
		    MultiFab& r, const IntVect& rat,
#ifdef HG_USE_CACHE
		    const copy_cache* border_cache,
#endif
		    const amr_restrictor_class& restric,
		    const level_interface& lev_interface,
		    const amr_boundary_class& bdy)
{
  for (int igrid = 0; igrid < dest.length(); igrid++) {
    if (bflag) {
      restrict_patch(dest[igrid], r, rat,
#ifdef HG_USE_CACHE
	  border_cache,
#endif
		     restric, lev_interface, bdy);
    }
    else {
      restrict_patch(dest[igrid], dest.box(igrid), r, rat, 
#ifdef HG_USE_CACHE
	  border_cache,
#endif
		     restric, lev_interface, bdy);
    }
  }
}
