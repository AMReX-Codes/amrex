
#include "restrictor.H"
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FACRST1  acrst1_
#  define FANRST1  anrst1_
#  define FANRST2  anrst2_
#  define FANFR2   anfr2_
#  define FANER2   aner2_
#  define FANCR2   ancr2_
#  define FANOR2   anor2_
#  define FANIR2   anir2_
#  define FANDR2   andr2_
#else
#  define FACRST1  ACRST1
#  define FANRST1  ANRST1
#  define FANRST2  ANRST2
#  define FANFR2   ANFR2
#  define FANER2   ANER2
#  define FANCR2   ANCR2
#  define FANOR2   ANOR2
#  define FANIR2   ANIR2
#  define FANDR2   ANDR2
#endif

extern "C" {
  void FACRST1(Real*, intS, intS, Real*, intS, intRS, const int&);
  void FANRST1(Real*, intS, intS, Real*, intS, intRS);
  void FANRST2(Real*, intS, intS, Real*, intS, intRS, const int&);
  void FANFR2(Real*, intS, intS, Real*, intS,
	      intRS, const int&, const int&, const int&);
  void FANER2(Real*, intS, intS, Real*, intS,
	      intRS, const int*, const int*, const int&);
  void FANCR2(Real*, intS, intS, Real*, intS, intRS, const int*, const int&);
  void FANOR2(Real*, intS, intS, Real*, intS,
	      intRS, const int&, const int&, const int&);
  void FANIR2(Real*, intS, intS, Real*, intS,
	      intRS, const int&, const int&, const int&);
  void FANDR2(Real*, intS, intS, Real*, intS, intRS, const int&, const int&);
}

Box cell_average_restrictor_class::box(const Box& fb, const IntVect& rat) const
{
  Box retbox(fb);
  return retbox.coarsen(rat);
}

void cell_average_restrictor_class::fill(FArrayBox& patch,
					 const Box& region,
					 FArrayBox& fgr,
					 const IntVect& rat) const
{
  assert(patch.box().cellCentered());
  for (int i = 0; i < patch.nComp(); i++) {
    FACRST1(patch.dataPtr(i), dimlist(patch.box()), dimlist(region),
	    fgr.dataPtr(i), dimlist(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]), integrate);
  }
}

void terrain_velocity_restrictor_class::fill(FArrayBox& patch,
					     const Box& region,
					     FArrayBox& fgr,
					     const IntVect& rat) const
{
  assert(patch.box().cellCentered());
  assert(patch.nComp() == 1);
  FACRST1(patch.dataPtr(), dimlist(patch.box()), dimlist(region),
	  fgr.dataPtr(), dimlist(fgr.box()),
	  D_DECL(rat[0], rat[1], rat[2]), 1);
  Real fac = 1.0 / rat[integrate];
  patch.mult(fac, region);
}

Box injection_restrictor_class::box(const Box& fb, const IntVect& rat) const
{
  Box retbox(fb);
  return retbox.coarsen(rat);
}

void injection_restrictor_class::fill(FArrayBox& patch,
				      const Box& region,
				      FArrayBox& fgr,
				      const IntVect& rat) const
{
  if (patch.box().type() == IntVect::TheNodeVector()) {
    for (int i = 0; i < patch.nComp(); i++) {
      FANRST1(patch.dataPtr(i), dimlist(patch.box()), dimlist(region),
	      fgr.dataPtr(i), dimlist(fgr.box()),
	      D_DECL(rat[0], rat[1], rat[2]));
    }
  }
  else
    BoxLib::Error("injection_restrictor_class::fill---Injection only defined for NODE-based data");
}

Box default_restrictor_class::box(const Box& fb, const IntVect& rat) const
{
  Box retbox(fb);
  return retbox.coarsen(rat);
}

void default_restrictor_class::fill(FArrayBox& patch,
				    const Box& region,
				    FArrayBox& fgr,
				    const IntVect& rat) const
{
  if (patch.box().cellCentered())
    cell_average_restrictor.fill(patch, region, fgr, rat);
  else if (patch.box().type() == IntVect::TheNodeVector())
    injection_restrictor.fill(patch, region, fgr, rat);
  else
    BoxLib::Error("default_restrictor_class::fill---No default restriction defined for mixed data");
}

Box bilinear_restrictor_class::box(const Box& fb, const IntVect& rat) const
{
  Box retbox(fb);
  return retbox.coarsen(rat).grow(-1);
}

void bilinear_restrictor_class::fill(FArrayBox& patch,
				     const Box& region,
				     FArrayBox& fgr,
				     const IntVect& rat) const
{
  if (patch.box().type() == IntVect::TheNodeVector()) {
    for (int i = 0; i < patch.nComp(); i++) {
      FANRST2(patch.dataPtr(i), dimlist(patch.box()), dimlist(region),
	      fgr.dataPtr(i), dimlist(fgr.box()),
	      D_DECL(rat[0], rat[1], rat[2]), integrate);
    }
  }
  else
    BoxLib::Error("bilinear_restrictor_coarse_class::fill---Bilinear restriction only defined for NODE-based data");
}

void bilinear_restrictor_class::interface(FArrayBox& patch,
					  const Box& region,
					  MultiFab& fine,
					  const copy_cache* border_cache,
					  const level_interface& interface,
					  amr_boundary bdy,
					  const IntVect& rat) const
{
  if (patch.box().type() != IntVect::TheNodeVector())
    BoxLib::Error("bilinear_restrictor_coarse_class::interface---bilinear restriction only defined for NODE-based data");

  Box regplus = grow(region,1);
  const Box& pb = patch.box();

  int ratmax = rat[0];
  ratmax = (rat[1] > ratmax ? rat[1] : ratmax);
#if (BL_SPACEDIM == 3)
  ratmax = (rat[2] > ratmax ? rat[2] : ratmax);
#endif

  if (fine.nGrow() < ratmax - 1) {
    for (int iface = 0; iface < interface.nfaces(); iface++) {
      if (interface.fgeo(iface) == level_interface::ALL && interface.fflag(iface) == 0) {
	// fine grid on both sides
	Box cbox = interface.node_face(iface);
	IntVect t = interface.face(iface).type();
	cbox.coarsen(rat);
	if (region.intersects(cbox)) {
	  // This extends fine face by one coarse cell past coarse face:
	  cbox &= regplus;
	  // Note:  Uses numerical values of index types:
	  cbox.grow(t - IntVect::TheUnitVector());
	  FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), patch.nComp());
	  fill_patch(fgr, fine, interface, bdy, 0, level_interface::FACEDIM, iface);
	  const Box& fb = fgr.box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fgr.dataPtr(i), dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
    }
#if (BL_SPACEDIM == 3)
    for (int iedge = 0; iedge < interface.nedges(); iedge++) {
      if (interface.egeo(iedge) == level_interface::ALL && interface.eflag(iedge) == 0) {
	// fine grid on all sides
	Box cbox = interface.node_edge(iedge);
	IntVect t = interface.edge(iedge).type();
	cbox.coarsen(rat);
	if (region.intersects(cbox)) {
	  // This extends fine edge by one coarse cell past coarse edge:
	  cbox &= regplus;
	  // Note:  Uses numerical values of index types:
	  cbox.grow(t - IntVect::TheUnitVector());
	  FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), patch.nComp());
	  fill_patch(fgr, fine, interface, bdy, 0, 1, iedge);
	  const Box& fb = fgr.box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fgr.dataPtr(i), dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
    }
#endif
    for (int icor = 0; icor < interface.ncorners(); icor++) {
      if (interface.cgeo(icor) == level_interface::ALL && interface.cflag(icor) == 0) {
	// fine grid on all sides
	Box cbox = interface.corner(icor);
	cbox.coarsen(rat);
	if (region.intersects(cbox)) {
	  FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), patch.nComp());
	  fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	  const Box& fb = fgr.box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fgr.dataPtr(i), dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
    }
  }
  else {
    fill_borders(fine, border_cache, interface, bdy, ratmax - 1);
    for (int iface = 0; iface < interface.nfaces(); iface++) {
      if (interface.fgeo(iface) == level_interface::ALL && interface.fflag(iface) == 0) {
	// fine grid on both sides
	Box cbox = interface.node_face(iface);
	IntVect t = interface.face(iface).type();
	cbox.coarsen(rat);
	if (region.intersects(cbox)) {
	  // This extends fine face by one coarse cell past coarse face:
	  cbox &= regplus;
	  // Note:  Uses numerical values of index types:
	  cbox.grow(t - IntVect::TheUnitVector());
	  int igrid = interface.fgrid(iface, 0);
	  if (igrid < 0)
	    igrid = interface.fgrid(iface, 1);
	  const Box& fb = fine[igrid].box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    Real *const fptr = fine[igrid].dataPtr(i);
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fptr, dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
    }
#if (BL_SPACEDIM == 3)
    for (int iedge = 0; iedge < interface.nedges(); iedge++) {
      if (interface.egeo(iedge) == level_interface::ALL && interface.eflag(iedge) == 0) {
	// fine grid on both sides
	Box cbox = interface.node_edge(iedge);
	IntVect t = interface.edge(iedge).type();
	cbox.coarsen(rat);
	if (region.intersects(cbox)) {
	  // This extends fine edge by one coarse cell past coarse edge:
	  cbox &= regplus;
	  // Note:  Uses numerical values of index types:
	  cbox.grow(t - IntVect::TheUnitVector());
	  int igrid = interface.egrid(iedge, 0);
	  for (int itmp = 1; igrid < 0; itmp++)
	    igrid = interface.egrid(iedge, itmp);
	  const Box& fb = fine[igrid].box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    Real *const fptr = fine[igrid].dataPtr(i);
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fptr, dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
    }
#endif
    for (int icor = 0; icor < interface.ncorners(); icor++) {
      if (interface.cgeo(icor) == level_interface::ALL && interface.cflag(icor) == 0) {
	// fine grid on all sides
	Box cbox = interface.corner(icor);
	cbox.coarsen(rat);
	if (region.intersects(cbox)) {
	  int igrid = interface.cgrid(icor, 0);
	  for (int itmp = 1; igrid < 0; itmp++)
	    igrid = interface.cgrid(icor, itmp);
	  const Box& fb = fine[igrid].box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    Real *const fptr = fine[igrid].dataPtr(i);
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fptr, dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
    }
  }
}

void bilinear_restrictor_coarse_class::interface(FArrayBox& patch,
						 const Box& region,
						 MultiFab& fine,
						 const copy_cache*
						   border_cache,
						 const level_interface&
						   interface,
						 amr_boundary bdy,
						 const IntVect& rat) const
{
  if (patch.box().type() != IntVect::TheNodeVector())
    BoxLib::Error("bilinear_restrictor_coarse_class::interface---bilinear restriction only defined for NODE-based data");

  Box regplus = grow(region,1);
  const Box& pb = patch.box();

  int ratmax = rat[0];
  ratmax = (rat[1] > ratmax ? rat[1] : ratmax);
#if (BL_SPACEDIM == 3)
  ratmax = (rat[2] > ratmax ? rat[2] : ratmax);
#endif

  if (fine.nGrow() >= ratmax - 1)
    fill_borders(fine, border_cache, interface, bdy, ratmax - 1);

  for (int iface = 0; iface < interface.nfaces(); iface++) {
    if (interface.fflag(iface) == 1)
      continue;
    Box cbox = interface.node_face(iface);
    IntVect t = interface.face(iface).type();
    unsigned geo = interface.fgeo(iface);
    cbox.coarsen(rat);
    if (region.intersects(cbox)) {
      // This extends fine face by one coarse cell past coarse face:
      cbox &= regplus;
      int idim = interface.fdim(iface);
      cbox.grow(t - IntVect::TheUnitVector());
      if (geo == level_interface::ALL) { // fine grid on both sides
	if (fine.nGrow() >= ratmax - 1) {
	  int igrid = interface.fgrid(iface, 0);
	  if (igrid < 0)
	    igrid = interface.fgrid(iface, 1);
	  const Box& fb = fine[igrid].box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    Real *const fptr = fine[igrid].dataPtr(i);
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fptr, dimlist(fb),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
	else {
	  Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
	  FArrayBox fgr(fbox, patch.nComp());
	  fill_patch(fgr, fine, interface, bdy, 0, level_interface::FACEDIM, iface);
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fgr.dataPtr(i), dimlist(fbox),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
      }
      else { // fine grid on just one side
	int idir = (geo & level_interface::LOW) ? -1 : 1;
	int igrid = (idir < 0) ? interface.fgrid(iface, 0) :
	                         interface.fgrid(iface, 1) ;
	if (igrid >= 0) {
	  // Usual case, a fine grid extends all along the face.
	  const Box& fb = fine[igrid].box();
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANFR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		   fine[igrid].dataPtr(i), dimlist(fb),
		   D_DECL(rat[0], rat[1], rat[2]), idim, idir, integrate);
	  }
	}
	else {
	  // A virtual fine grid is on the other side of the boundary.
	  Box fbox = refine(cbox, rat).grow(rat - IntVect::TheUnitVector());
	  if (geo & level_interface::LOW)
	    fbox.growHi(idim, 1 - rat[idim]);
	  else
	    fbox.growLo(idim, 1 - rat[idim]);
	  FArrayBox fgr(fbox, patch.nComp());
	  fill_patch(fgr, fine, interface, bdy, 0, level_interface::FACEDIM, iface);
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANFR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		   fgr.dataPtr(i), dimlist(fbox),
		   D_DECL(rat[0], rat[1], rat[2]), idim, idir, integrate);
	  }
	}
      }
    }
  }

#if (BL_SPACEDIM == 2)

  for (int icor = 0; icor < interface.ncorners(); icor++) {
    if (interface.cflag(icor) == 1)
      continue;
    Box cbox = interface.corner(icor);
    cbox.coarsen(rat);
    if (region.intersects(cbox)) {
      unsigned geo = interface.cgeo(icor);
      if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) { // fine grid on all sides
	int igrid = interface.cgrid(icor, 0);
	for (int itmp = 1; igrid < 0; itmp++)
	  igrid = interface.cgrid(icor, itmp);
	const Box& fb = fine[igrid].box();
	for (int i = 0; i < patch.nComp(); i++) {
	  Real *const fptr = fine[igrid].dataPtr(i);
	  FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		  fptr, dimlist(fb),
		  D_DECL(rat[0], rat[1], rat[2]), integrate);
	}
      }
      else if (geo == level_interface::ALL) { // fine grid on all sides
	FArrayBox fgr(refine(grow(cbox, 1), rat), patch.nComp());
	fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	const Box& fb = fgr.box();
	for (int i = 0; i < patch.nComp(); i++) {
	  FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		  fgr.dataPtr(i), dimlist(fb),
		  D_DECL(rat[0], rat[1], rat[2]), integrate);
	}
      }
      else if (geo == XL || geo == XH || geo == YL || geo == YH) {
	// fine grid on two adjacent sides
	int idim = (geo == XL || geo == XH) ? 0 : 1;
	int idir = (geo & LL) ? -1 : 1;
	Box fbox = refine(cbox, rat).grow(1 - idim, rat[1-idim]);
	if (geo & LL)
	  fbox.growLo(idim, rat[idim]);
	else
	  fbox.growHi(idim, rat[idim]);
	FArrayBox fgr(fbox, patch.nComp());
	fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	for (int i = 0; i < patch.nComp(); i++) {
	  FANFR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		 fgr.dataPtr(i), dimlist(fbox),
		 D_DECL(rat[0], rat[1], rat[2]), idim, idir, integrate);
	}
      }
      else if (geo == LL || geo == HL || geo == LH || geo == HH) {
	// outside corner
	Box fbox = refine(cbox, rat);
	int idir0, idir1;
	if (geo & XL) {
	  fbox.growLo(0, rat[0]);
	  idir0 = -1;
	}
	else {
	  fbox.growHi(0, rat[0]);
	  idir0 = 1;
	}
	if (geo & YL) {
	  fbox.growLo(1, rat[1]);
	  idir1 = -1;
	}
	else {
	  fbox.growHi(1, rat[1]);
	  idir1 = 1;
	}
	FArrayBox fgr(fbox, patch.nComp());
	fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	for (int i = 0; i < patch.nComp(); i++) {
	  FANOR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		 fgr.dataPtr(i), dimlist(fbox),
		 D_DECL(rat[0], rat[1], rat[2]), idir0, idir1, integrate);
	}
      }
      else if (geo == (LL | HH) || geo == (LH | HL)) {
	// diagonal corner
	Box fbox = refine(cbox, rat).grow(rat);
	FArrayBox fgr(fbox, patch.nComp());
	int idir1 = (geo == (LL | HH)) ? 1 : -1;
	fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	for (int i = 0; i < patch.nComp(); i++) {
	  FANDR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		 fgr.dataPtr(i), dimlist(fbox),
		 D_DECL(rat[0], rat[1], rat[2]), idir1, integrate);
	}
      }
      else {
	// inside corner
	Box fbox = refine(cbox, rat).grow(rat);
	FArrayBox fgr(fbox, patch.nComp());
	int idir0 = ((geo & XL) == XL) ? -1 : 1;
	int idir1 = ((geo & YL) == YL) ? -1 : 1;
	fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	for (int i = 0; i < patch.nComp(); i++) {
	  FANIR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		 fgr.dataPtr(i), dimlist(fbox),
		 D_DECL(rat[0], rat[1], rat[2]), idir0, idir1, integrate);
	}
      }
    }
  }

#elif (BL_SPACEDIM == 3)

  int ga[level_interface::N_CORNER_GRIDS];

  for (int iedge = 0; iedge < interface.nedges(); iedge++) {
    if (interface.eflag(iedge) == 1)
      continue;
    Box cbox = interface.node_edge(iedge);
    IntVect t = interface.edge(iedge).type();
    cbox.coarsen(rat);
    if (region.intersects(cbox)) {
      // This extends fine edge by one coarse cell past coarse face:
      cbox &= regplus;
      cbox.grow(t - IntVect::TheUnitVector());
      unsigned geo = interface.egeo(iedge);
      if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) {
	int igrid = interface.egrid(iedge, 0);
	for (int itmp = 1; igrid < 0; itmp++)
	  igrid = interface.egrid(iedge, itmp);
	const Box& fb = fine[igrid].box();
	for (int i = 0; i < patch.nComp(); i++) {
	  Real *const fptr = fine[igrid].dataPtr(i);
	  FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		  fptr, dimlist(fb),
		  D_DECL(rat[0], rat[1], rat[2]), integrate);
	}
      }
      else {
	Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
	FArrayBox fgr(fbox, patch.nComp());
	fill_patch(fgr, fine, interface, bdy, 0, 1, iedge);
	if (geo == level_interface::ALL) { // fine grid on all sides
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fgr.dataPtr(i), dimlist(fbox),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
	else {
	  interface.geo_array(ga, 1, iedge);
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANER2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		   fgr.dataPtr(i), dimlist(fbox),
		   D_DECL(rat[0], rat[1], rat[2]), t.getVect(), ga, integrate);
	  }
	}
      }
    }
  }

  for (int icor = 0; icor < interface.ncorners(); icor++) {
    if (interface.cflag(icor) == 1)
      continue;
    Box cbox = interface.corner(icor);
    cbox.coarsen(rat);
    if (region.intersects(cbox)) {
      unsigned geo = interface.cgeo(icor);
      if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) {
	int igrid = interface.cgrid(icor, 0);
	for (int itmp = 1; igrid < 0; itmp++)
	  igrid = interface.cgrid(icor, itmp);
	const Box& fb = fine[igrid].box();
	for (int i = 0; i < patch.nComp(); i++) {
	  Real *const fptr = fine[igrid].dataPtr(i);
	  FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		  fptr, dimlist(fb),
		  D_DECL(rat[0], rat[1], rat[2]), integrate);
	}
      }
      else {
	Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
	FArrayBox fgr(fbox, patch.nComp());
	fill_patch(fgr, fine, interface, bdy, 0, 0, icor);
	if (geo == level_interface::ALL) { // fine grid on all sides
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANRST2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		    fgr.dataPtr(i), dimlist(fbox),
		    D_DECL(rat[0], rat[1], rat[2]), integrate);
	  }
	}
	else {
	  interface.geo_array(ga, 0, icor);
	  for (int i = 0; i < patch.nComp(); i++) {
	    FANCR2(patch.dataPtr(i), dimlist(pb), dimlist(cbox),
		   fgr.dataPtr(i), dimlist(fbox),
		   D_DECL(rat[0], rat[1], rat[2]), ga, integrate);
	  }
	}
      }
    }
  }
#endif
}
