
#include "restrictor.H"
#include "fill_patch.H"

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_FACRST1  acrst1_
#define FORT_FANRST1  anrst1_
#define FORT_FANRST2  anrst2_
#define FORT_FANFR2   anfr2_
#define FORT_FANER2   aner2_
#define FORT_FANCR2   ancr2_
#define FORT_FANOR2   anor2_
#define FORT_FANIR2   anir2_
#define FORT_FANDR2   andr2_
#else
#define FORT_FACRST1  ACRST1
#define FORT_FANRST1  ANRST1
#define FORT_FANRST2  ANRST2
#define FORT_FANFR2   ANFR2
#define FORT_FANER2   ANER2
#define FORT_FANCR2   ANCR2
#define FORT_FANOR2   ANOR2
#define FORT_FANIR2   ANIR2
#define FORT_FANDR2   ANDR2
#endif

extern "C" 
{
    void FORT_FACRST1(Real*, intS, intS, const Real*, intS, intRS, const int&, const int*);
    void FORT_FANRST1(Real*, intS, intS, const Real*, intS, intRS, const int&);
    void FORT_FANRST2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int*);
    void FORT_FANFR2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int&, const int&, const int*);
    void FORT_FANER2(Real*, intS, intS, const Real*, intS, intRS, const int*, const int*, const int&, const int*);
    void FORT_FANCR2(Real*, intS, intS, const Real*, intS, intRS, const int*, const int&, const int*);
#if BL_SPACEDIM == 2
    void FORT_FANOR2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int&, const int&, const int *);
    void FORT_FANIR2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int&, const int&, const int *);
    void FORT_FANDR2(Real*, intS, intS, const Real*, intS, intRS, const int&, const int&, const int*);
#endif
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
    assert(patch.nComp() == 1);
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]), patch.nComp(), &integrate);
}

void terrain_velocity_restrictor_class::fill(FArrayBox& patch,
					     const Box& region,
					     const FArrayBox& fgr,
					     const IntVect& rat) const
{
    assert(patch.box().cellCentered());
    assert(patch.nComp() == 1);
    const int integ = 1;
    FORT_FACRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]), 1, &integ);
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
    assert(patch.box().type() == IntVect::TheNodeVector());
    FORT_FANRST1(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]), patch.nComp());
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
    assert(patch.box().cellCentered() || patch.box().type() == IntVect::TheNodeVector());
    if (patch.box().cellCentered())
    {
	cell_average_restrictor_class(0).fill(patch, region, fgr, rat);
    }
    else if (patch.box().type() == IntVect::TheNodeVector())
    {
	injection_restrictor_class().fill(patch, region, fgr, rat);
    }
}

bilinear_restrictor_class::bilinear_restrictor_class(int i, bool hg_terrain) : integrate(i), m_hg_terrain(hg_terrain)
{
    assert(i == 0 || i == 1);
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
    assert(patch.box().type() == IntVect::TheNodeVector());
    FORT_FANRST2(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region), fgr.dataPtr(), DIMLIST(fgr.box()),
	D_DECL(rat[0], rat[1], rat[2]), patch.nComp(), &integrate);
}

void bilinear_restrictor_class::fill_interface(MultiFab& dest,
					       MultiFab& fine,
					       const level_interface& lev_interface,
					       const amr_boundary_class* bdy,
					       const IntVect& rat) const
{
    assert(type(dest) == IntVect::TheNodeVector());
    for (int jgrid = 0; jgrid < dest.length(); jgrid++)
    {
	const Box& region = dest.box(jgrid);
	// Interface restriction is sufficiently rare and specialized that
	// we will let the restrictor handle it---at least for now.
	
	// This assertion difficult in BoxLib since r.mesh() is not cc:
	//assert(r.mesh() == lev_interface.interior_mesh());
	
	Box regplus = grow(region,1);
	const Box& pb = dest[jgrid].box();
	
	int ratmax = rat[0];
	for ( int i = 1; i < BL_SPACEDIM; ++i )
	    ratmax = (rat[i] > ratmax) ? rat[i] : ratmax;
	
	if (fine.nGrow() < ratmax - 1) 
	{
	    // PARALLEL
	    for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	    {
		if (lev_interface.geo(level_interface::FACEDIM, iface) == level_interface::ALL && !lev_interface.flag(level_interface::FACEDIM, iface) ) 
		{
		    // fine grid on both sides
		    Box cbox = lev_interface.node_box(level_interface::FACEDIM, iface);
		    IntVect t = lev_interface.box(level_interface::FACEDIM, iface).type();
		    cbox.coarsen(rat);
		    if (region.intersects(cbox)) 
		    {
			// This extends fine face by one coarse cell past coarse face:
			cbox &= regplus;
			// Note:  Uses numerical values of index types:
			cbox.grow(t - IntVect::TheUnitVector());
			FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), dest[jgrid].nComp());
			fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, level_interface::FACEDIM, iface);
			const Box& fb = fgr.box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
	    }
#if (BL_SPACEDIM == 3)
	    // PARALLEL
	    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
	    {
		if (lev_interface.geo(1, iedge) == level_interface::ALL && !lev_interface.flag(1, iedge) ) 
		{
		    // fine grid on all sides
		    Box cbox = lev_interface.node_box(1, iedge);
		    IntVect t = lev_interface.box(1, iedge).type();
		    cbox.coarsen(rat);
		    if (region.intersects(cbox)) 
		    {
			// This extends fine edge by one coarse cell past coarse edge:
			cbox &= regplus;
			// Note:  Uses numerical values of index types:
			cbox.grow(t - IntVect::TheUnitVector());
			FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), dest[jgrid].nComp());
			fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 1, iedge);
			const Box& fb = fgr.box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
	    }
#endif
	    // PARALLEL
	    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
	    {
		if (lev_interface.geo(0, icor) == level_interface::ALL && !lev_interface.flag(0, icor) ) 
		{
		    // fine grid on all sides
		    Box cbox = lev_interface.box(0, icor);
		    cbox.coarsen(rat);
		    if (region.intersects(cbox)) 
		    {
			FArrayBox fgr(grow(refine(cbox, rat), rat - IntVect::TheUnitVector()), dest[jgrid].nComp());
			fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
			const Box& fb = fgr.box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
	    }
	}
	else 
	{
	    fill_borders(fine, lev_interface, bdy, ratmax - 1, m_hg_terrain);
	    // PARALLEL
	    for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	    {
		if (lev_interface.geo(level_interface::FACEDIM, iface) == level_interface::ALL && !lev_interface.flag(level_interface::FACEDIM, iface) ) 
		{
		    // fine grid on both sides
		    Box cbox = lev_interface.node_box(level_interface::FACEDIM, iface);
		    IntVect t = lev_interface.box(level_interface::FACEDIM, iface).type();
		    cbox.coarsen(rat);
		    if (region.intersects(cbox)) 
		    {
			// This extends fine face by one coarse cell past coarse face:
			cbox &= regplus;
			// Note:  Uses numerical values of index types:
			cbox.grow(t - IntVect::TheUnitVector());
			int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
			if (igrid < 0)
			    igrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
			const Box& fb = fine[igrid].box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
	    }
#if (BL_SPACEDIM == 3)
	    // PARALLEL
	    for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
	    {
		if (lev_interface.geo(1, iedge) == level_interface::ALL && !lev_interface.flag(1, iedge) ) 
		{
		    // fine grid on both sides
		    Box cbox = lev_interface.node_box(1, iedge);
		    IntVect t = lev_interface.box(1, iedge).type();
		    cbox.coarsen(rat);
		    if (region.intersects(cbox)) 
		    {
			// This extends fine edge by one coarse cell past coarse edge:
			cbox &= regplus;
			// Note:  Uses numerical values of index types:
			cbox.grow(t - IntVect::TheUnitVector());
			int igrid = lev_interface.grid(1, iedge, 0);
			for (int itmp = 1; igrid < 0; itmp++)
			    igrid = lev_interface.grid(1, iedge, itmp);
			const Box& fb = fine[igrid].box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
	    }
#endif
	    // PARALLEL
	    for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
	    {
		if (lev_interface.geo(0, icor) == level_interface::ALL && !lev_interface.flag(0, icor) ) 
		{
		    // fine grid on all sides
		    Box cbox = lev_interface.box(0, icor);
		    cbox.coarsen(rat);
		    if (region.intersects(cbox)) 
		    {
			int igrid = lev_interface.grid(0, icor, 0);
			for (int itmp = 1; igrid < 0; itmp++)
			    igrid = lev_interface.grid(0, icor, itmp);
			const Box& fb = fine[igrid].box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
	    }
	}
    }
}

void bilinear_restrictor_coarse_class::fill_interface(MultiFab& dest,
						      MultiFab& fine,
						      const level_interface& lev_interface,
						      const amr_boundary_class* bdy,
						      const IntVect& rat) const
{
    assert(type(dest) == IntVect::TheNodeVector());
    for (int jgrid = 0; jgrid < dest.length(); jgrid++)
    {
	const Box& region = dest.box(jgrid);
	// Interface restriction is sufficiently rare and specialized that
	// we will let the restrictor handle it---at least for now.
	
	// This assertion difficult in BoxLib since r.mesh() is not cc:
	//assert(r.mesh() == lev_interface.interior_mesh());
	
	const Box regplus = grow(region,1);
	const Box& pb = dest[jgrid].box();
	
	int ratmax = rat[0];
	for (int i = 1; i < BL_SPACEDIM; ++i)
	    ratmax = (rat[i] > ratmax) ? rat[i] : ratmax;
	
	if (fine.nGrow() >= ratmax - 1)
	    fill_borders(fine, lev_interface, bdy, ratmax - 1, m_hg_terrain);
	
	// PARALLEL
	for (int iface = 0; iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
	{
	    if ( lev_interface.flag(level_interface::FACEDIM, iface) )
		continue;
	    Box cbox = lev_interface.node_box(level_interface::FACEDIM, iface);
	    IntVect t = lev_interface.box(level_interface::FACEDIM, iface).type();
	    const unsigned int geo = lev_interface.geo(level_interface::FACEDIM, iface);
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
			int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
			if (igrid < 0)
			    igrid = lev_interface.grid(level_interface::FACEDIM, iface, 1);
			const Box& fb = fine[igrid].box();
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb), 
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		    else 
		    {
			Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
			FArrayBox fgr(fbox, dest[jgrid].nComp());
			fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, level_interface::FACEDIM, iface);
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox), 
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		}
		else 
		{ // fine grid on just one side
		    const int idir = (geo & level_interface::LOW) ? -1 : 1;
		    const int igrid = (idir < 0) ? lev_interface.grid(level_interface::FACEDIM, iface, 0) :
		    lev_interface.grid(level_interface::FACEDIM, iface, 1) ;
		    if (igrid >= 0) 
		    {
			// Usual case, a fine grid extends all along the face.
			const Box& fb = fine[igrid].box();
			FORT_FANFR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb), 
			    D_DECL(rat[0], rat[1], rat[2]), idim, idir, dest.nComp(), &integrate);
		    }
		    else 
		    {
			// A virtual fine grid is on the other side of the boundary.
			Box fbox = refine(cbox, rat).grow(rat - IntVect::TheUnitVector());
			if (geo & level_interface::LOW)
			    fbox.growHi(idim, 1 - rat[idim]);
			else
			    fbox.growLo(idim, 1 - rat[idim]);
			FArrayBox fgr(fbox, dest[jgrid].nComp());
			fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, level_interface::FACEDIM, iface);
			FORT_FANFR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), idim, idir, dest.nComp(), &integrate);
		    }
		}
	    }
	}
	
#if (BL_SPACEDIM == 2)
	
	// PARALLEL
	for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
	{
	    if ( lev_interface.flag(0, icor) )
		continue;
	    Box cbox = lev_interface.box(0, icor);
	    cbox.coarsen(rat);
	    if (region.intersects(cbox)) 
	    {
		const unsigned int geo = lev_interface.geo(0, icor);
		if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
		{ 
		    // fine grid on all sides
		    int igrid = lev_interface.grid(0, icor, 0);
		    for (int itmp = 1; igrid < 0; itmp++)
			igrid = lev_interface.grid(0, icor, itmp);
		    const Box& fb = fine[igrid].box();
		    FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb),
			D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		}
		else if (geo == level_interface::ALL) 
		{
		    // fine grid on all sides
		    FArrayBox fgr(refine(grow(cbox, 1), rat), dest[jgrid].nComp());
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
		    const Box& fb = fgr.box();
		    FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fb),
			D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		}
		else if (geo == level_interface::XL || geo == level_interface::XH || geo == level_interface::YL || geo == level_interface::YH) 
		{
		    // fine grid on two adjacent sides
		    const int idim = (geo == level_interface::XL || geo == level_interface::XH) ? 0 : 1;
		    const int idir = (geo & level_interface::LL) ? -1 : 1;
		    Box fbox = refine(cbox, rat).grow(1 - idim, rat[1-idim]);
		    if (geo & level_interface::LL)
			fbox.growLo(idim, rat[idim]);
		    else
			fbox.growHi(idim, rat[idim]);
		    FArrayBox fgr(fbox, dest[jgrid].nComp());
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
		    FORT_FANFR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idim, idir, dest.nComp(), &integrate);
		}
		else if (geo == level_interface::LL || geo == level_interface::HL || geo == level_interface::LH || geo == level_interface::HH) 
		{
		    // outside corner
		    Box fbox = refine(cbox, rat);
		    int idir0;
		    int idir1;
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
		    FArrayBox fgr(fbox, dest[jgrid].nComp());
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
		    FORT_FANOR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idir0, idir1, dest.nComp(), &integrate);
		}
		else if (geo == (level_interface::LL | level_interface::HH) || geo == (level_interface::LH | level_interface::HL)) 
		{
		    // diagonal corner
		    Box fbox = refine(cbox, rat).grow(rat);
		    FArrayBox fgr(fbox, dest[jgrid].nComp());
		    const int idir1 = (geo == (level_interface::LL | level_interface::HH)) ? 1 : -1;
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
		    FORT_FANDR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idir1, dest.nComp(), &integrate);
		}
		else 
		{
		    // inside corner
		    Box fbox = refine(cbox, rat).grow(rat);
		    FArrayBox fgr(fbox, dest[jgrid].nComp());
		    const int idir0 = ((geo & level_interface::XL) == level_interface::XL) ? -1 : 1;
		    const int idir1 = ((geo & level_interface::YL) == level_interface::YL) ? -1 : 1;
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
		    FORT_FANIR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			D_DECL(rat[0], rat[1], rat[2]), idir0, idir1, dest.nComp(), &integrate);
		}
	    }
	}
	
#elif (BL_SPACEDIM == 3)
	
	// PARALLEL  
	for (int iedge = 0; iedge < lev_interface.nboxes(1); iedge++) 
	{
	    if ( lev_interface.flag(1, iedge) )
		continue;
	    Box cbox = lev_interface.node_box(1, iedge);
	    IntVect t = lev_interface.box(1, iedge).type();
	    cbox.coarsen(rat);
	    if (region.intersects(cbox)) 
	    {
		// This extends fine edge by one coarse cell past coarse face:
		cbox &= regplus;
		cbox.grow(t - IntVect::TheUnitVector());
		const unsigned int geo = lev_interface.geo(1, iedge);
		if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
		{
		    int igrid = lev_interface.grid(1, iedge, 0);
		    for (int itmp = 1; igrid < 0; itmp++)
			igrid = lev_interface.grid(1, iedge, itmp);
		    const Box& fb = fine[igrid].box();
		    FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb),
			D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		}
		else 
		{
		    Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
		    FArrayBox fgr(fbox, dest[jgrid].nComp());
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 1, iedge);
		    if (geo == level_interface::ALL) 
		    { 
			// fine grid on all sides
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		    else 
		    {
			Array<int> ga = lev_interface.geo_array(1, iedge);
			FORT_FANER2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), t.getVect(), ga.dataPtr(), dest.nComp(), &integrate);
		    }
		}
	    }
	}
	// PARALLEL  
	for (int icor = 0; icor < lev_interface.nboxes(0); icor++) 
	{
	    if ( lev_interface.flag(0, icor) )
		continue;
	    Box cbox = lev_interface.box(0, icor);
	    cbox.coarsen(rat);
	    if (region.intersects(cbox)) 
	    {
		const unsigned int geo = lev_interface.geo(0, icor);
		if (geo == level_interface::ALL && fine.nGrow() >= ratmax - 1) 
		{
		    int igrid = lev_interface.grid(0, icor, 0);
		    for (int itmp = 1; igrid < 0; itmp++)
			igrid = lev_interface.grid(0, icor, itmp);
		    const Box& fb = fine[igrid].box();
		    FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fine[igrid].dataPtr(), DIMLIST(fb),
			D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		}
		else 
		{
		    Box fbox = grow(refine(cbox, rat), rat - IntVect::TheUnitVector());
		    FArrayBox fgr(fbox, dest[jgrid].nComp());
		    fill_patch(fgr, fgr.box(), fine, lev_interface, bdy, 0, icor);
		    if (geo == level_interface::ALL) 
		    { 
			// fine grid on all sides
			FORT_FANRST2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), dest.nComp(), &integrate);
		    }
		    else 
		    {
			Array<int> ga = lev_interface.geo_array(0, icor);
			FORT_FANCR2(dest[jgrid].dataPtr(), DIMLIST(pb), DIMLIST(cbox), fgr.dataPtr(), DIMLIST(fbox),
			    D_DECL(rat[0], rat[1], rat[2]), ga.dataPtr(), dest.nComp(), &integrate);
		    }
		}
	    }
	}
#endif
    }
}
