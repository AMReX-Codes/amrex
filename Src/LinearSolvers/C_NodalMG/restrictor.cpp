
#include "restrictor.H"
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_FACRST1  acrst1_
#  define FORT_FANRST1  anrst1_
#  define FORT_FANRST2  anrst2_
#  define FORT_FANFR2   anfr2_
#  define FORT_FANER2   aner2_
#  define FORT_FANCR2   ancr2_
#  define FORT_FANOR2   anor2_
#  define FORT_FANIR2   anir2_
#  define FORT_FANDR2   andr2_
#else
#  define FORT_FACRST1  ACRST1
#  define FORT_FANRST1  ANRST1
#  define FORT_FANRST2  ANRST2
#  define FORT_FANFR2   ANFR2
#  define FORT_FANER2   ANER2
#  define FORT_FANCR2   ANCR2
#  define FORT_FANOR2   ANOR2
#  define FORT_FANIR2   ANIR2
#  define FORT_FANDR2   ANDR2
#endif

extern "C" 
{
    void FORT_FACRST1(Real*, intS, intS, const Real*, intS, intRS, const int&);
    void FORT_FANRST1(Real*, intS, intS, const Real*, intS, intRS);
    void FORT_FANRST2(Real*, intS, intS, const Real*, intS, intRS, const int&);
    void FORT_FANFR2(Real*, intS, intS, const Real*, intS,
	intRS, const int&, const int&, const int&);
    void FORT_FANER2(Real*, intS, intS, const Real*, intS,
	intRS, const int*, const int*, const int&);
    void FORT_FANCR2(Real*, intS, intS, const Real*, intS, intRS, const int*, const int&);
    void FORT_FANOR2(Real*, intS, intS, const Real*, intS,
	intRS, const int&, const int&, const int&);
    void FORT_FANIR2(Real*, intS, intS, const Real*, intS,
	intRS, const int&, const int&, const int&);
    void FORT_FANDR2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int&);
}

Box cell_average_restrictor_class::box(const Box& fb, const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat);
}

void cell_average_restrictor_class::fill(FArrayBox& patch,
					 const Box& region,
					 const FArrayBox& fgr,
					 const IntVect& rat) const
{
    assert(patch.box().cellCentered());
    for (int i = 0; i < patch.nComp(); i++) 
    {
	FORT_FACRST1(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
	    fgr.dataPtr(i), DIMLIST(fgr.box()),
	    D_DECL(rat[0], rat[1], rat[2]), integrate);
    }
}

void terrain_velocity_restrictor_class::fill(FArrayBox& patch,
					     const Box& region,
					     const FArrayBox& fgr,
					     const IntVect& rat) const
{
    assert(patch.box().cellCentered());
    assert(patch.nComp() == 1);
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
	fgr.dataPtr(), DIMLIST(fgr.box()),
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
				      const FArrayBox& fgr,
				      const IntVect& rat) const
{
    if (patch.box().type() == IntVect::TheNodeVector()) 
    {
	for (int i = 0; i < patch.nComp(); i++) 
	{
	    FORT_FANRST1(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		fgr.dataPtr(i), DIMLIST(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]));
	}
    }
    else
	BoxLib::Error("injection_restrictor_class::fill---Injection only defined for NODE-based data");
}

Box default_restrictor::box(const Box& fb, const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat);
}

void default_restrictor::fill(FArrayBox& patch,
			      const Box& region,
			      const FArrayBox& fgr,
			      const IntVect& rat) const
{
    if (patch.box().cellCentered())
	cell_average_restrictor_class(0).fill(patch, region, fgr, rat);
    else if (patch.box().type() == IntVect::TheNodeVector())
	injection_restrictor_class().fill(patch, region, fgr, rat);
    else
	BoxLib::Error("default_restrictor::fill---No default restriction defined for mixed data");
}

Box bilinear_restrictor_class::box(const Box& fb, const IntVect& rat) const
{
    Box retbox(fb);
    return retbox.coarsen(rat).grow(-1);
}

void bilinear_restrictor_class::fill(FArrayBox& patch,
				     const Box& region,
				     const FArrayBox& fgr,
				     const IntVect& rat) const
{
    if (patch.box().type() == IntVect::TheNodeVector()) 
    {
	for (int i = 0; i < patch.nComp(); i++) 
	{
	    FORT_FANRST2(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		fgr.dataPtr(i), DIMLIST(fgr.box()),
		D_DECL(rat[0], rat[1], rat[2]), integrate);
	}
    }
    else
	BoxLib::Error("bilinear_restrictor_coarse_class::fill---Bilinear restriction only defined for NODE-based data");
}

void bilinear_restrictor_class::lev_interface(FArrayBox& patch,
					      const Box& region,
					      MultiFab& fine,
#ifdef HG_USE_CACHE
					      const copy_cache* border_cache,
#endif
					      const level_interface& lev_interface,
					      const amr_boundary_class& bdy,
					      const IntVect& rat) const
{
    if (patch.box().type() != IntVect::TheNodeVector())
	BoxLib::Error("bilinear_restrictor_coarse_class::lev_interface---bilinear restriction only defined for NODE-based data");
    
    Box regplus = grow(region,1);
    const Box& pb = patch.box();
    
    int ratmax = rat[0];
    ratmax = (rat[1] > ratmax ? rat[1] : ratmax);
#if (BL_SPACEDIM == 3)
    ratmax = (rat[2] > ratmax ? rat[2] : ratmax);
#endif
    
    if (fine.nGrow() < ratmax - 1) 
    {
	for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
	{
	    if (lev_interface.fgeo(iface) == level_interface::ALL && lev_interface.fflag(iface) == 0) 
	    {
		// fine grid on both sides
		Box cbox = lev_interface.node_face(iface);
		IntVect t = lev_interface.face(iface).type();
		cbox.coarsen(rat);
		if (region.intersects(cbox)) 
		{
		    // This extends fine face by one coarse cell past coarse face:
		    cbox &= regplus;
		    // Note:  Uses numerical values of index types:
		    cbox.grow(t - IntVect::TheUnitVector());
		    FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), patch.nComp());
		    fill_patch(fgr, fine, lev_interface, bdy, 0, level_interface::FACEDIM, iface);
		    const Box& fb = fgr.box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fgr.dataPtr(i), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	}
#if (BL_SPACEDIM == 3)
	for (int iedge = 0; iedge < lev_interface.nedges(); iedge++) 
	{
	    if (lev_interface.egeo(iedge) == level_interface::ALL && lev_interface.eflag(iedge) == 0) 
	    {
		// fine grid on all sides
		Box cbox = lev_interface.node_edge(iedge);
		IntVect t = lev_interface.edge(iedge).type();
		cbox.coarsen(rat);
		if (region.intersects(cbox)) 
		{
		    // This extends fine edge by one coarse cell past coarse edge:
		    cbox &= regplus;
		    // Note:  Uses numerical values of index types:
		    cbox.grow(t - IntVect::TheUnitVector());
		    FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), patch.nComp());
		    fill_patch(fgr, fine, lev_interface, bdy, 0, 1, iedge);
		    const Box& fb = fgr.box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fgr.dataPtr(i), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	}
#endif
	for (int icor = 0; icor < lev_interface.ncorners(); icor++) 
	{
	    if (lev_interface.cgeo(icor) == level_interface::ALL && lev_interface.cflag(icor) == 0) 
	    {
		// fine grid on all sides
		Box cbox = lev_interface.corner(icor);
		cbox.coarsen(rat);
		if (region.intersects(cbox)) 
		{
		    FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), patch.nComp());
		    fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
		    const Box& fb = fgr.box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fgr.dataPtr(i), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	}
    }
    else 
    {
	fill_borders(fine, 
#ifdef HG_USE_CACHE
	    border_cache, 
#endif
	    lev_interface, bdy, ratmax - 1);
	for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
	{
	    if (lev_interface.fgeo(iface) == level_interface::ALL && lev_interface.fflag(iface) == 0) 
	    {
		// fine grid on both sides
		Box cbox = lev_interface.node_face(iface);
		IntVect t = lev_interface.face(iface).type();
		cbox.coarsen(rat);
		if (region.intersects(cbox)) 
		{
		    // This extends fine face by one coarse cell past coarse face:
		    cbox &= regplus;
		    // Note:  Uses numerical values of index types:
		    cbox.grow(t - IntVect::TheUnitVector());
		    int igrid = lev_interface.fgrid(iface, 0);
		    if (igrid < 0)
			igrid = lev_interface.fgrid(iface, 1);
		    const Box& fb = fine[igrid].box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			Real *const fptr = fine[igrid].dataPtr(i);
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fptr, DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	}
#if (BL_SPACEDIM == 3)
	for (int iedge = 0; iedge < lev_interface.nedges(); iedge++) 
	{
	    if (lev_interface.egeo(iedge) == level_interface::ALL && lev_interface.eflag(iedge) == 0) 
	    {
		// fine grid on both sides
		Box cbox = lev_interface.node_edge(iedge);
		IntVect t = lev_interface.edge(iedge).type();
		cbox.coarsen(rat);
		if (region.intersects(cbox)) 
		{
		    // This extends fine edge by one coarse cell past coarse edge:
		    cbox &= regplus;
		    // Note:  Uses numerical values of index types:
		    cbox.grow(t - IntVect::TheUnitVector());
		    int igrid = lev_interface.egrid(iedge, 0);
		    for (int itmp = 1; igrid < 0; itmp++)
			igrid = lev_interface.egrid(iedge, itmp);
		    const Box& fb = fine[igrid].box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			Real *const fptr = fine[igrid].dataPtr(i);
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fptr, DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	}
#endif
	for (int icor = 0; icor < lev_interface.ncorners(); icor++) 
	{
	    if (lev_interface.cgeo(icor) == level_interface::ALL && lev_interface.cflag(icor) == 0) 
	    {
		// fine grid on all sides
		Box cbox = lev_interface.corner(icor);
		cbox.coarsen(rat);
		if (region.intersects(cbox)) 
		{
		    int igrid = lev_interface.cgrid(icor, 0);
		    for (int itmp = 1; igrid < 0; itmp++)
			igrid = lev_interface.cgrid(icor, itmp);
		    const Box& fb = fine[igrid].box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			Real *const fptr = fine[igrid].dataPtr(i);
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fptr, DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	}
    }
}

void bilinear_restrictor_coarse_class::lev_interface(FArrayBox& patch,
						     const Box& region,
						     MultiFab& fine,
#ifdef HG_USE_CACHE
						     const copy_cache* border_cache,
#endif
						     const level_interface& lev_interface,
						     const amr_boundary_class& bdy,
						     const IntVect& rat) const
{
    if (patch.box().type() != IntVect::TheNodeVector())
	BoxLib::Error("bilinear_restrictor_coarse_class::lev_interface---bilinear restriction only defined for NODE-based data");
    
    Box regplus = grow(region,1);
    const Box& pb = patch.box();
    
    int ratmax = rat[0];
    ratmax = (rat[1] > ratmax ? rat[1] : ratmax);
#if (BL_SPACEDIM == 3)
    ratmax = (rat[2] > ratmax ? rat[2] : ratmax);
#endif
    
    if (fine.nGrow() >= ratmax - 1)
	fill_borders(fine, 
#ifdef HG_USE_CACHE
	border_cache, 
#endif
	lev_interface, bdy, ratmax - 1);
    
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	if (lev_interface.fflag(iface) == 1)
	    continue;
	Box cbox = lev_interface.node_face(iface);
	IntVect t = lev_interface.face(iface).type();
	unsigned geo = lev_interface.fgeo(iface);
	cbox.coarsen(rat);
	if (region.intersects(cbox)) 
	{
	    // This extends fine face by one coarse cell past coarse face:
	    cbox &= regplus;
	    int idim = lev_interface.fdim(iface);
	    cbox.grow(t - IntVect::TheUnitVector());
	    if (geo == level_interface::ALL) 
	    { 
		// fine grid on both sides
		if (fine.nGrow() >= ratmax - 1) 
		{
		    int igrid = lev_interface.fgrid(iface, 0);
		    if (igrid < 0)
			igrid = lev_interface.fgrid(iface, 1);
		    const Box& fb = fine[igrid].box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			Real *const fptr = fine[igrid].dataPtr(i);
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fptr, DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
		else 
		{
		    Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
		    FArrayBox fgr(fbox, patch.nComp());
		    fill_patch(fgr, fine, lev_interface, bdy, 0, level_interface::FACEDIM, iface);
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fgr.dataPtr(i), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), integrate);
		    }
		}
	    }
	    else 
	    { // fine grid on just one side
		int idir = (geo & level_interface::LOW) ? -1 : 1;
		int igrid = (idir < 0) ? lev_interface.fgrid(iface, 0) :
		lev_interface.fgrid(iface, 1) ;
		if (igrid >= 0) 
		{
		    // Usual case, a fine grid extends all along the face.
		    const Box& fb = fine[igrid].box();
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			FORT_FANFR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fine[igrid].dataPtr(i), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), idim, idir, integrate);
		    }
		}
		else 
		{
		    // A virtual fine grid is on the other side of the boundary.
		    Box fbox = refine(cbox, rat).grow(rat - IntVect::TheUnitVector());
		    if (geo & level_interface::LOW)
			fbox.growHi(idim, 1 - rat[idim]);
		    else
			fbox.growLo(idim, 1 - rat[idim]);
		    FArrayBox fgr(fbox, patch.nComp());
		    fill_patch(fgr, fine, lev_interface, bdy, 0, level_interface::FACEDIM, iface);
		    for (int i = 0; i < patch.nComp(); i++) 
		    {
			FORT_FANFR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			    fgr.dataPtr(i), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), idim, idir, integrate);
		    }
		}
	    }
	}
    }
    
#if (BL_SPACEDIM == 2)
    
    for (int icor = 0; icor < lev_interface.ncorners(); icor++) 
    {
	if (lev_interface.cflag(icor) == 1)
	    continue;
	Box cbox = lev_interface.corner(icor);
	cbox.coarsen(rat);
	if (region.intersects(cbox)) 
	{
	    unsigned geo = lev_interface.cgeo(icor);
	    if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
	    { 
		// fine grid on all sides
		int igrid = lev_interface.cgrid(icor, 0);
		for (int itmp = 1; igrid < 0; itmp++)
		    igrid = lev_interface.cgrid(icor, itmp);
		const Box& fb = fine[igrid].box();
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    Real *const fptr = fine[igrid].dataPtr(i);
		    FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			fptr, DIMLIST(fb),
			D_DECL(rat[0], rat[1], rat[2]), integrate);
		}
	    }
	    else if (geo == level_interface::ALL) 
	    {
		// fine grid on all sides
		FArrayBox fgr(refine(grow(cbox, 1), rat), patch.nComp());
		fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
		const Box& fb = fgr.box();
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			fgr.dataPtr(i), DIMLIST(fb),
			D_DECL(rat[0], rat[1], rat[2]), integrate);
		}
	    }
	    else if (geo == level_interface::XL || geo == level_interface::XH || geo == level_interface::YL || geo == level_interface::YH) 
	    {
		// fine grid on two adjacent sides
		int idim = (geo == level_interface::XL || geo == level_interface::XH) ? 0 : 1;
		int idir = (geo & level_interface::LL) ? -1 : 1;
		Box fbox = refine(cbox, rat).grow(1 - idim, rat[1-idim]);
		if (geo & level_interface::LL)
		    fbox.growLo(idim, rat[idim]);
		else
		    fbox.growHi(idim, rat[idim]);
		FArrayBox fgr(fbox, patch.nComp());
		fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FANFR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			fgr.dataPtr(i), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idim, idir, integrate);
		}
	    }
	    else if (geo == level_interface::LL || geo == level_interface::HL || geo == level_interface::LH || geo == level_interface::HH) 
	    {
		// outside corner
		Box fbox = refine(cbox, rat);
		int idir0, idir1;
		if (geo & level_interface::XL) 
		{
		    fbox.growLo(0, rat[0]);
		    idir0 = -1;
		}
		else 
		{
		    fbox.growHi(0, rat[0]);
		    idir0 = 1;
		}
		if (geo & level_interface::YL) 
		{
		    fbox.growLo(1, rat[1]);
		    idir1 = -1;
		}
		else 
		{
		    fbox.growHi(1, rat[1]);
		    idir1 = 1;
		}
		FArrayBox fgr(fbox, patch.nComp());
		fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FANOR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			fgr.dataPtr(i), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idir0, idir1, integrate);
		}
	    }
	    else if (geo == (level_interface::LL | level_interface::HH) || geo == (level_interface::LH | level_interface::HL)) 
	    {
		// diagonal corner
		Box fbox = refine(cbox, rat).grow(rat);
		FArrayBox fgr(fbox, patch.nComp());
		int idir1 = (geo == (level_interface::LL | level_interface::HH)) ? 1 : -1;
		fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FANDR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			fgr.dataPtr(i), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idir1, integrate);
		}
	    }
	    else 
	    {
		// inside corner
		Box fbox = refine(cbox, rat).grow(rat);
		FArrayBox fgr(fbox, patch.nComp());
		int idir0 = ((geo & level_interface::XL) == level_interface::XL) ? -1 : 1;
		int idir1 = ((geo & level_interface::YL) == level_interface::YL) ? -1 : 1;
		fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FANIR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			fgr.dataPtr(i), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idir0, idir1, integrate);
		}
	    }
    }
  }
  
#elif (BL_SPACEDIM == 3)
  
  int ga[level_interface::N_CORNER_GRIDS];
  
  for (int iedge = 0; iedge < lev_interface.nedges(); iedge++) 
  {
      if (lev_interface.eflag(iedge) == 1)
	  continue;
      Box cbox = lev_interface.node_edge(iedge);
      IntVect t = lev_interface.edge(iedge).type();
      cbox.coarsen(rat);
      if (region.intersects(cbox)) 
      {
	  // This extends fine edge by one coarse cell past coarse face:
	  cbox &= regplus;
	  cbox.grow(t - IntVect::TheUnitVector());
	  unsigned geo = lev_interface.egeo(iedge);
	  if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
	  {
	      int igrid = lev_interface.egrid(iedge, 0);
	      for (int itmp = 1; igrid < 0; itmp++)
		  igrid = lev_interface.egrid(iedge, itmp);
	      const Box& fb = fine[igrid].box();
	      for (int i = 0; i < patch.nComp(); i++) 
	      {
		  Real *const fptr = fine[igrid].dataPtr(i);
		  FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
		      fptr, DIMLIST(fb),
		      D_DECL(rat[0], rat[1], rat[2]), integrate);
	      }
	  }
	  else 
	  {
	      Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
	      FArrayBox fgr(fbox, patch.nComp());
	      fill_patch(fgr, fine, lev_interface, bdy, 0, 1, iedge);
	      if (geo == level_interface::ALL) 
	      { 
		  // fine grid on all sides
		  for (int i = 0; i < patch.nComp(); i++) 
		  {
		      FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			  fgr.dataPtr(i), DIMLIST(fbox),
			  D_DECL(rat[0], rat[1], rat[2]), integrate);
		  }
	      }
	      else 
	      {
		  lev_interface.geo_array(ga, 1, iedge);
		  for (int i = 0; i < patch.nComp(); i++) 
		  {
		      FORT_FANER2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			  fgr.dataPtr(i), DIMLIST(fbox),
			  D_DECL(rat[0], rat[1], rat[2]), t.getVect(), ga, integrate);
		  }
	      }
	  }
      }
  }
  
  for (int icor = 0; icor < lev_interface.ncorners(); icor++) 
  {
      if (lev_interface.cflag(icor) == 1)
	  continue;
      Box cbox = lev_interface.corner(icor);
      cbox.coarsen(rat);
      if (region.intersects(cbox)) 
      {
	  unsigned geo = lev_interface.cgeo(icor);
	  if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
	  {
	      int igrid = lev_interface.cgrid(icor, 0);
	      for (int itmp = 1; igrid < 0; itmp++)
		  igrid = lev_interface.cgrid(icor, itmp);
	      const Box& fb = fine[igrid].box();
	      for (int i = 0; i < patch.nComp(); i++) 
	      {
		  Real *const fptr = fine[igrid].dataPtr(i);
		  FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
		      fptr, DIMLIST(fb),
		      D_DECL(rat[0], rat[1], rat[2]), integrate);
	      }
	  }
	  else 
	  {
	      Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
	      FArrayBox fgr(fbox, patch.nComp());
	      fill_patch(fgr, fine, lev_interface, bdy, 0, 0, icor);
	      if (geo == level_interface::ALL) 
	      { 
		  // fine grid on all sides
		  for (int i = 0; i < patch.nComp(); i++) 
		  {
		      FORT_FANRST2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			  fgr.dataPtr(i), DIMLIST(fbox),
			  D_DECL(rat[0], rat[1], rat[2]), integrate);
		  }
	      }
	      else 
	      {
		  lev_interface.geo_array(ga, 0, icor);
		  for (int i = 0; i < patch.nComp(); i++) 
		  {
		      FORT_FANCR2(patch.dataPtr(i), DIMLIST(pb), DIMLIST(cbox),
			  fgr.dataPtr(i), DIMLIST(fbox),
			  D_DECL(rat[0], rat[1], rat[2]), ga, integrate);
		  }
	      }
	  }
      }
  }
#endif
}
