
#include "boundary.H"
#include "cache.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_FBREF    bref_
#  define FORT_FBREFM   brefm_
#  define FORT_FBNEG    bneg_
#  define FORT_FBNEGM   bnegm_
#  define FORT_FBINFLO  binflo_
#  define FORT_FBINFIL  binfil_
#else
#  define FORT_FBREF    BREF
#  define FORT_FBREFM   BREFM
#  define FORT_FBNEG    BNEG
#  define FORT_FBNEGM   BNEGM
#  define FORT_FBINFLO  BINFLO
#  define FORT_FBINFIL  BINFIL
#endif

extern "C" 
{
    void FORT_FBREF(Real*, intS, intS, const Real*, intS, intS, const int&);
    void FORT_FBREFM(Real*, intS, intS, const Real*, intS, intS, const int*);
    void FORT_FBNEG(Real*, intS, intS, const Real*, intS, intS, const int&);
    void FORT_FBNEGM(Real*, intS, intS, const Real*, intS, intS, const int*);
    void FORT_FBINFLO(Real*, intS, intS, const Real*, intS, intS, const int&);
    void FORT_FBINFIL(Real*, intS, intS, const Real*, intS, intS, const int&);
}

Box amr_boundary_class::box(const Box& region, const Box& domain, int idir) const
{
    const int idim = abs(idir) - 1;
    Box retbox(region);
    if (idir < 0) 
    {
	if (region.type(idim) == IndexType::CELL)
	    retbox.shift(idim, 2 * domain.smallEnd(idim) - 1 - region.bigEnd(idim) - region.smallEnd(idim));
	else
	    retbox.shift(idim, 2 * domain.smallEnd(idim) - region.bigEnd(idim) - region.smallEnd(idim));
    }
    else if (idir > 0) 
    {
	if (region.type(idim) == IndexType::CELL)
	    retbox.shift(idim, 2 * domain.bigEnd(idim) + 1 - region.bigEnd(idim) - region.smallEnd(idim));
	else
	    retbox.shift(idim, 2 * domain.bigEnd(idim) + 2 - region.bigEnd(idim) - region.smallEnd(idim));
    }
    else 
    {
	BoxLib::Error("amr_boundary_class::box---undefined boundary direction");
    }
    return retbox;
}

Box mixed_boundary_class::box(const Box& region, const Box& domain, int idir) const
{
    const int idim = abs(idir) - 1;
    RegType t = ptr->bc[idim][idir > 0];
    Box retbox(region);
    
    if (t == refWall || t == outflow || (t == inflow && idim != flowdim)) 
    {
	// all these cases use a reflected box
	if (idir < 0) 
	{
	    if (region.type(idim) == IndexType::CELL)
		retbox.shift(idim, 2 * domain.smallEnd(idim) - 1 - region.bigEnd(idim) - region.smallEnd(idim));
	    else
		retbox.shift(idim, 2 * domain.smallEnd(idim) - region.bigEnd(idim) - region.smallEnd(idim));
	}
	else if (idir > 0) 
	{
	    if (region.type(idim) == IndexType::CELL)
		retbox.shift(idim, 2 * domain.bigEnd(idim) + 1 - region.bigEnd(idim) - region.smallEnd(idim));
	    else
		retbox.shift(idim, 2 * domain.bigEnd(idim) + 2 - region.bigEnd(idim) - region.smallEnd(idim));
	}
	else 
	{
	    BoxLib::Error("mixed_boundary_class::box---undefined boundary direction");
	}
    }
    else if (t == inflow) 
    {
	// This case needs the boundary node or the first cell outside
	// the domain to extrapolate.
	// For cell-based, the fill_patch call must set the use-ghost-cell
	// option (flags & 2) or an infinite recursion will result.
	// It's a kludge, hopefully temporary, so fill_patch does not
	// test to see if this is the case.
	if (idir < 0) 
	{
	    if (region.type(idim) == IndexType::CELL) 
	    {
		retbox.shift(idim, 2 * domain.smallEnd(idim) - 1 - region.bigEnd(idim) - region.smallEnd(idim));
		retbox.setSmall(idim, domain.smallEnd(idim) - 1);
	    }
	    else 
	    {
		retbox.shift(idim, 2 * domain.smallEnd(idim) - region.bigEnd(idim) - region.smallEnd(idim));
		retbox.setSmall(idim, domain.smallEnd(idim));
	    }
	}
	else if (idir > 0) 
	{
	    if (region.type(idim) == IndexType::CELL) 
	    {
		retbox.shift(idim, 2 * domain.bigEnd(idim) + 1 - region.bigEnd(idim) - region.smallEnd(idim));
		retbox.setBig(idim, domain.bigEnd(idim) + 1);
	    }
	    else 
	    {
		retbox.shift(idim, 2 * domain.bigEnd(idim) + 2 - region.bigEnd(idim) - region.smallEnd(idim));
		retbox.setBig(idim, domain.bigEnd(idim) + 1);
	    }
	}
	else 
	{
	    BoxLib::Error("mixed_boundary_class::box---undefined boundary direction");
	}
    }
    else if (t == periodic) 
    {
	if (idir < 0)
	    retbox.shift(idim, domain.length(idim));
	else if (idir > 0)
	    retbox.shift(idim, -domain.length(idim));
	else
	    BoxLib::Error("mixed_boundary_class::box---undefined boundary direction");
    }
    else 
    {
	BoxLib::Error("mixed_boundary_class::box---boundary type not supported");
    }
    return retbox;
}

// Reflects on all outflow cases (which aren't called anyway).
// On velocity inflow, uses box function which extends interior
// box just past edge of domain.

//void mixed_boundary_class::fill(FArrayBox& patch,
//				const Box& region,
//				const MultiFab& src,
//				int igrid, const Box& domain) const
void mixed_boundary_class::fill(FArrayBox& patch,
				const Box& region,
				const FArrayBox& src, int igrid,
				const Box& domain) const
{
    Box tdomain = domain;
    tdomain.convert(type(src));
    Box idomain = grow(tdomain, IntVect::TheZeroVector() - type(src));
    Box image = region;
    int refarray[BL_SPACEDIM], negflag = 1;
    int idir = 0;
    
    int negarray[BL_SPACEDIM-1];
    for (int i = 0; i < BL_SPACEDIM - 1; i++) 
    {
	negarray[i] = 1;
    }
    
    for (int idim = 0; idim < BL_SPACEDIM; idim++) 
    {
	refarray[idim] = 0;
	if (region.bigEnd(idim) < idomain.smallEnd(idim)) 
	{
	    RegType t = ptr->bc[idim][0];
	    if (t == inflow && idim == flowdim) 
	    {
		refarray[idim] = 0;
		idir = -1 - idim;
	    }
	    else if (t == refWall || t == inflow || t == outflow) 
	    {
		refarray[idim] = 1;
		image.shift(idim, tdomain.smallEnd(idim) + idomain.smallEnd(idim) - 1 - region.bigEnd(idim) - region.smallEnd(idim));
		if (flowdim == -3 || t == refWall && idim == flowdim)
		    negflag = -negflag;
		if (flowdim == -4) 
		{
		    if (idim < BL_SPACEDIM - 1) 
		    {
			negarray[idim] = -negarray[idim];
		    }
		    else 
		    {
			for (int i = 0; i < BL_SPACEDIM - 1; i++) 
			{
			    negarray[i] = -negarray[i];
			}
		    }
		}
	    }
	    else if (t == periodic) 
	    {
		refarray[idim] = 0;
		image.shift(idim, domain.length(idim));
	    }
	}
	else if (region.smallEnd(idim) > idomain.bigEnd(idim)) 
	{
	    RegType t = ptr->bc[idim][1];
	    if (t == inflow && idim == flowdim) 
	    {
		refarray[idim] = 0;
		idir = 1 + idim;
	    }
	    if (t == refWall || t == inflow || t == outflow) 
	    {
		refarray[idim] = 1;
		image.shift(idim, tdomain.bigEnd(idim) + idomain.bigEnd(idim) + 1 - region.bigEnd(idim) - region.smallEnd(idim));
		if (flowdim == -3 || t == refWall && idim == flowdim)
		    negflag = -negflag;
		if (flowdim == -4) 
		{
		    if (idim < BL_SPACEDIM - 1) 
		    {
			negarray[idim] = -negarray[idim];
		    }
		    else 
		    {
			for (int i = 0; i < BL_SPACEDIM - 1; i++) 
			{
			    negarray[i] = -negarray[i];
			}
		    }
		}
	    }
	    else if (t == periodic) 
	    {
		refarray[idim] = 0;
		image.shift(idim, -domain.length(idim));
	    }
	}
    }
    
    if (idir != 0) 
    {
	// normal-component inflow section, assume patch.nComp() == 1
	if (igrid >= 0) 
	{
	    Box bb = box(image, domain, idir);
	    if (image == region) 
	    {
		// only bdy involved, can fill directly from interior
		fill(patch, region, src, bb, domain, idir);
	    }
	    else 
	    {
		// multiple bdys, fill intermediate patch
		FArrayBox gb(image);
		fill(gb, image, src, bb, domain, idir);
		if (negflag == 1) 
		{
		    FORT_FBREFM(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
			gb.dataPtr(), DIMLIST(image),
			DIMLIST(image), refarray);
		}
		else if (negflag == -1) 
		{
		    FORT_FBNEGM(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
			gb.dataPtr(), DIMLIST(image),
			DIMLIST(image), refarray);
		}
	    }
	}
	else 
	{
	    BoxLib::Error("mixed_boundary_class::fill---exterior ref undefined");
	}
    }
    else if (flowdim == -4) 
    {
	assert(igrid >= 0);
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    FORT_FBREFM(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		src.dataPtr(i), DIMLIST(src.box()),
		DIMLIST(image), refarray);
	}
	for (int idim = 0; idim < BL_SPACEDIM - 1; idim++) 
	{
	    int i = idim + BL_SPACEDIM;
	    if (negarray[idim] == 1) 
	    {
		FORT_FBREFM(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    src.dataPtr(i), DIMLIST(src.box()),
		    DIMLIST(image), refarray);
	    }
	    else if (negarray[idim] == -1) 
	    {
		FORT_FBNEGM(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    src.dataPtr(i), DIMLIST(src.box()),
		    DIMLIST(image), refarray);
	    }
	}
    }
    else 
    {
	// all cases other than normal-component inflow
	if (negflag == 1) 
	{
	    if (igrid >= 0) 
	    {
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FBREFM(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
			src.dataPtr(i), DIMLIST(src.box()),
			DIMLIST(image), refarray);
		}
	    }
	    else 
	    {
		BoxLib::Error("mixed_boundary_class::fill---exterior ref undefined");
	    }
	}
	else if (negflag == -1) 
	{
	    if (igrid >= 0) 
	    {
		for (int i = 0; i < patch.nComp(); i++) 
		{
		    FORT_FBNEGM(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
			src.dataPtr(i), DIMLIST(src.box()),
			DIMLIST(image), refarray);
		}
	    }
	    else 
	    {
		BoxLib::Error("mixed_boundary_class::fill---exterior ref undefined");
	    }
	}
    }
}

void mixed_boundary_class::fill(FArrayBox& patch,
				const Box& region,
				const FArrayBox& bgr,
				const Box& bb,
				const Box& domain,
				int idir) const
{
    const int idim = abs(idir) - 1;
    RegType t = ptr->bc[idim][idir > 0];
    
    if (flowdim == -4 && (t == refWall || t == inflow)) 
    {
	// terrain sigma
	BoxLib::Error("mixed_boundary_class::fill---terrain undefined");
    }
    else if (t == refWall) 
    {
	if (idim == flowdim || flowdim == -3) 
	{
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBNEG(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
	else 
	{
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBREF(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
    }
    else if (t == periodic) 
    {
	patch.copy(bgr, bb, 0, region, 0, patch.nComp());
    }
    else if (t == inflow) 
    {
	if (flowdim == -2) 
	{
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBREF(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
	else if (flowdim == -1) 
	{
	    //BoxLib::Error("mixed_boundary_class::Don't know how to do inflow density");
	    // Inflow density---just reflect interior for now
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBREF(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
	else if (flowdim == -3) 
	{
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBNEG(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
	else if (idim == flowdim) 
	{
	    // For this to work, fill_borders must already have been called
	    // to initialize values in the first ghost cell outside the domain.
	    FORT_FBINFIL(patch.dataPtr(), DIMLIST(patch.box()), DIMLIST(region),
		bgr.dataPtr(), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	}
	else if (flowdim >= 0) 
	{
	    // transverse velocity components
	    //patch.assign(0.0, region);
	    // we now believe this looks like a refWall to transverse components
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBREF(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
    }
    else if (t == outflow) 
    {
	// Do nothing if NODE-based, reflect if CELL-based 
	if (type(patch,idim) == IndexType::CELL) 
	{
	    for (int i = 0; i < patch.nComp(); i++) 
	    {
		FORT_FBREF(patch.dataPtr(i), DIMLIST(patch.box()), DIMLIST(region),
		    bgr.dataPtr(i), DIMLIST(bgr.box()), DIMLIST(bb), idim);
	    }
	}
    }
    else 
    {
	BoxLib::Error("mixed_boundary_class::fill---boundary type not supported");
    }
}

void mixed_boundary_class::sync_borders(MultiFab& r, const level_interface& lev_interface) const
{
    assert(type(r) == IntVect::TheNodeVector());
    
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	if (lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	int igrid = lev_interface.fgrid(iface, 0);
	if (igrid < 0) 
	{
	    int idim = lev_interface.fdim(iface);
	    if (ptr->bc[idim][0] == periodic) 
	    {
		int jgrid = lev_interface.fgrid(iface, 1);
		igrid = lev_interface.exterior_ref(igrid);
		const Box& b = lev_interface.node_face(iface);
		Box bb = b;
		bb.shift(idim, lev_interface.domain().length(idim));
		r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());
	    }
	}
    }
}

void mixed_boundary_class::fill_borders(MultiFab& r, const level_interface& lev_interface, int w) const
{
    w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;
    assert( w == 1 || w == 0 );
    const Box& domain = lev_interface.domain();
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	if (lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	const int igrid = lev_interface.fgrid(iface, 0);
	const int jgrid = lev_interface.fgrid(iface, 1);
	if (igrid < 0 || jgrid < 0) 
	{
	    Box b = lev_interface.face(iface);
	    const int idim = lev_interface.fdim(iface);
	    const int a = (type(r,idim) == IndexType::NODE);
	    // need to do on x borders too in case y border is an interior face
#if (BL_SPACEDIM == 2)
	    //b.grow(1 - idim, w); // replaced by grow stmts below?
#else
	    //could cause overwrite problem for periodic boundaries:
	    //b.grow((idim + 1) % 3, w).grow((idim + 2) % 3, w);
#endif
	    if (igrid < 0) 
	    {
		for (int i = 0; i < BL_SPACEDIM; i++) 
		{
		    if (i != idim) 
		    {
			if (lev_interface.interior_mesh()[jgrid].smallEnd(i) == b.smallEnd(i))
			    b.growLo(i, w);
			if (lev_interface.interior_mesh()[jgrid].bigEnd(i) == b.bigEnd(i))
			    b.growHi(i, w);
		    }
		}
		b.shift(idim, -a).growLo(idim, w-a).convert(type(r));
		RegType t = ptr->bc[idim][0];
		Box bb = b;
		if (flowdim == -4 && (t == refWall || t == inflow)) 
		{
		    // terrain sigma
		    bb.shift(idim, 2 * domain.smallEnd(idim) - 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
		    const Box& rbox = r[jgrid].box();
		    for (int i = 0; i < r.nComp(); i++) 
		    {
			Real *const rptr = r[jgrid].dataPtr(i);
			if ((i == idim + BL_SPACEDIM) ||
			    (i >= BL_SPACEDIM && idim == BL_SPACEDIM - 1)) 
			{
			    FORT_FBNEG(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
			}
			else 
			{
			    FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
			}
		    }
		}
		else if (t == refWall) 
		{
		    bb.shift(idim, 2 * domain.smallEnd(idim) - 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
		    const Box& rbox = r[jgrid].box();
		    if (idim == flowdim || flowdim == -3) 
		    {
			for (int i = 0; i < r.nComp(); i++) 
			{
			    Real *const rptr = r[jgrid].dataPtr(i);
			    FORT_FBNEG(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
			}
		    }
		    else 
		    {
			for (int i = 0; i < r.nComp(); i++) 
			{
			    Real *const rptr = r[jgrid].dataPtr(i);
			    FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
			}
		    }
		}
		else if (t == periodic) 
		{
		    int isrc = lev_interface.exterior_ref(igrid);
		    bb.shift(idim, domain.length(idim));
		    r[jgrid].copy(r[isrc], bb, 0, b, 0, r.nComp());
		}
		else if (t == inflow) 
		{
		    bb.shift(idim, 2 * domain.smallEnd(idim) - 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
		    const Box& rbox = r[jgrid].box();
		    if (flowdim == -2) 
		    {
			Real *const rptr = r[jgrid].dataPtr();
			FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		    }
		    else if (flowdim == -1) 
		    {
			//BoxLib::Error("mixed_boundary_class::Don't know how to do inflow density");
			// Inflow density---just reflect interior for now
			for (int i = 0; i < r.nComp(); i++) 
			{
			    Real *const rptr = r[jgrid].dataPtr(i);
			    FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
			}
		    }
		    else if (flowdim == -3) 
		    {
			Real *const rptr = r[jgrid].dataPtr();
			FORT_FBNEG(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		    }
		    else if (idim == flowdim) 
		    {
			Real *const rptr = r[jgrid].dataPtr();
			// For this to work, fill_borders must be called exactly
			// once for each level of this variable.
			FORT_FBINFLO(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		    }
		    else if (flowdim >= 0) 
		    {
			// transverse velocity components
			//r[jgrid].assign(0.0, b);
			// we now believe this looks like a refWall to transverse comps
			Real *const rptr = r[jgrid].dataPtr();
			FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		    }
		}
		else if (t == outflow) 
		{
		    // Do nothing if NODE-based, reflect if CELL-based 
		    if (type(r,idim) == IndexType::CELL) 
		    {
			bb.shift(idim, 2 * domain.smallEnd(idim) - 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
			const Box& rbox = r[jgrid].box();
			for (int i = 0; i < r.nComp(); i++) 
			{
			    Real *const rptr = r[jgrid].dataPtr(i);
			    FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
			}
		    }
		}
      }
      else if (jgrid < 0) 
      {
	  for (int i = 0; i < BL_SPACEDIM; i++) 
	  {
	      if (i != idim) 
	      {
		  if (lev_interface.interior_mesh()[igrid].smallEnd(i) == b.smallEnd(i))
		      b.growLo(i, w);
		  if (lev_interface.interior_mesh()[igrid].bigEnd(i) == b.bigEnd(i))
		      b.growHi(i, w);
	      }
	  }
	  b.shift(idim, a).growHi(idim, w-a).convert(type(r));
	  RegType t = ptr->bc[idim][1];
	  Box bb = b;
	  if (flowdim == -4 && (t == refWall || t == inflow)) 
	  {
	      // terrain sigma
	      bb.shift(idim, 2 * domain.bigEnd(idim) + 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
	      const Box& rbox = r[igrid].box();
	      for (int i = 0; i < r.nComp(); i++) 
	      {
		  Real *const rptr = r[igrid].dataPtr(i);
		  if ((i == idim + BL_SPACEDIM) ||
		      (i >= BL_SPACEDIM && idim == BL_SPACEDIM - 1)) 
		  {
		      FORT_FBNEG(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		  }
		  else 
		  {
		      FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		  }
	      }
	  }
	  else if (t == refWall) 
	  {
	      bb.shift(idim, 2 * domain.bigEnd(idim) + 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
	      const Box& rbox = r[igrid].box();
	      if (idim == flowdim || flowdim == -3) 
	      {
		  for (int i = 0; i < r.nComp(); i++) 
		  {
		      Real *const rptr = r[igrid].dataPtr(i);
		      FORT_FBNEG(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		  }
	      }
	      else 
	      {
		  for (int i = 0; i < r.nComp(); i++) 
		  {
		      Real *const rptr = r[igrid].dataPtr(i);
		      FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		  }
	      }
	  }
	  else if (t == periodic) 
	  {
	      int isrc = lev_interface.exterior_ref(jgrid);
	      bb.shift(idim, -domain.length(idim));
	      r[igrid].copy(r[isrc], bb, 0, b, 0, r.nComp());
	  }
	  else if (t == inflow) 
	  {
	      bb.shift(idim, 2 * domain.bigEnd(idim) + 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
	      const Box& rbox = r[igrid].box();
	      if (flowdim == -2) 
	      {
		  Real *const rptr = r[igrid].dataPtr();
		  FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
	      }
	      else if (flowdim == -1) 
	      {
		  //BoxLib::Error("mixed_boundary_class::Don't know how to do inflow density");
		  // Inflow density---just reflect interior for now
		  for (int i = 0; i < r.nComp(); i++) 
		  {
		      Real *const rptr = r[igrid].dataPtr(i);
		      FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		  }
	      }
	      else if (flowdim == -3) 
	      {
		  Real *const rptr = r[igrid].dataPtr();
		  FORT_FBNEG(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
	      }
	      else if (idim == flowdim) 
	      {
		  Real *const rptr = r[igrid].dataPtr();
		  // For this to work, fill_borders must be called exactly
		  // once for each level of this variable.
		  FORT_FBINFLO(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
	      }
	      else if (flowdim >= 0) 
	      {
		  // transverse velocity components
		  //r[igrid].assign(0.0, rbox);
		  // we now believe this looks like a refWall to transverse comps
		  Real *const rptr = r[igrid].dataPtr();
		  FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
	      }
	  }
	  else if (t == outflow) 
	  {
	      // Do nothing if NODE-based, reflect if CELL-based 
	      if (type(r,idim) == IndexType::CELL) 
	      {
		  bb.shift(idim, 2 * domain.bigEnd(idim) + 1 + a - b.bigEnd(idim) - b.smallEnd(idim));
		  const Box& rbox = r[igrid].box();
		  for (int i = 0; i < r.nComp(); i++) 
		  {
		      Real *const rptr = r[igrid].dataPtr(i);
		      FORT_FBREF(rptr, DIMLIST(rbox), DIMLIST(b), rptr, DIMLIST(rbox), DIMLIST(bb), idim);
		  }
	      }
	  }
      }
    }
  }
}

void mixed_boundary_class::set_sync_cache(copy_cache* cache,
					  int nsets, int& iset,
					  const MultiFab& r,
					  const level_interface&
					  lev_interface) const
{
    if (r.nComp() != 1)
	BoxLib::Error("mixed_boundary_class::set_sync_cache---only single components currently supported");
    if (type(r) != IntVect::TheNodeVector())
	BoxLib::Error("mixed_boundary_class::set_sync_cache---only NODE-based sync defined");
    
    assert( cache != 0);
    Real *const baseptr = cache->dataPtr();
    
    const Box& domain = lev_interface.domain();
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	if (lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	int igrid = lev_interface.fgrid(iface, 0);
	if (igrid < 0) 
	{
	    int idim = lev_interface.fdim(iface);
	    if (ptr->bc[idim][0] == periodic) 
	    {
		int jgrid = lev_interface.fgrid(iface, 1);
		igrid = lev_interface.exterior_ref(igrid);
		const Box& b = lev_interface.node_face(iface);
#if (BL_SPACEDIM == 2)
		int dstart, sstart, dstrid, sstrid, nvals;
		dstrid = r[jgrid].box().length(0);
		sstrid = r[igrid].box().length(0);
		if (idim == 0) 
		{
		    nvals = b.length(1);
		    dstart = r[jgrid].dataPtr() - baseptr + 
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			dstrid * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
		    sstart = r[igrid].dataPtr() - baseptr +
			(b.smallEnd(0) + domain.length(0) -
			r[igrid].box().smallEnd(0)) +
			sstrid * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
		}
		else 
		{
		    nvals = b.length(0);
		    dstart = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			dstrid * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
		    sstart = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			sstrid * (b.smallEnd(1) + domain.length(1) -
			r[igrid].box().smallEnd(1));
		    dstrid = 1;
		    sstrid = 1;
		}
		cache->set(iset++, dstart, sstart, dstrid, sstrid, nvals);
#else
		int dstart, sstart, dstrid1, dstrid2, sstrid1, sstrid2, nvals1, nvals2;
		dstrid1 = r[jgrid].box().length(0);
		dstrid2 = dstrid1 * r[jgrid].box().length(1);
		sstrid1 = r[igrid].box().length(0);
		sstrid2 = sstrid1 * r[igrid].box().length(1);
		if (idim == 0) 
		{
		    nvals1 = b.length(1);
		    nvals2 = b.length(2);
		    dstart = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			dstrid1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			dstrid2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    sstart = r[igrid].dataPtr() - baseptr +
			(b.smallEnd(0) + domain.length(0) -
			r[igrid].box().smallEnd(0)) +
			sstrid1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			sstrid2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		}
		else if (idim == 1) 
		{
		    nvals1 = b.length(0);
		    nvals2 = b.length(2);
		    dstart = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			dstrid1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			dstrid2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    sstart = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			sstrid1 * (b.smallEnd(1) + domain.length(1) -
			r[igrid].box().smallEnd(1)) +
			sstrid2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
		    dstrid1 = 1;
		    sstrid1 = 1;
		}
		else 
		{
		    nvals1 = b.length(0);
		    nvals2 = b.length(1);
		    dstart = r[jgrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			dstrid1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
			dstrid2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
		    sstart = r[igrid].dataPtr() - baseptr +
			b.smallEnd(0) - r[igrid].box().smallEnd(0) +
			sstrid1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
			sstrid2 * (b.smallEnd(2) + domain.length(2) -
			r[igrid].box().smallEnd(2));
		    dstrid2 = dstrid1;
		    sstrid2 = sstrid1;
		    dstrid1 = 1;
		    sstrid1 = 1;
		}
		cache->set(iset++, dstart, sstart, dstrid1, dstrid2, sstrid1, sstrid2, nvals1, nvals2);
#endif
	    }
	}
    }
    if (iset > nsets)
	BoxLib::Error("mixed_boundary_class::set_sync_cache---too many boundary edges to cache");
}

void mixed_boundary_class::set_border_cache(copy_cache* cache,
					    int nsets, int& iset,
					    const MultiFab& r,
					    const level_interface& lev_interface,
					    int w) const
{
    assert( w == 0 || w == 1 );
    if (r.nComp() != 1)
	BoxLib::Error("mixed_boundary_class::set_border_cache---only single components currently supported");
    if (w != 1 || type(r) != IntVect::TheNodeVector())
	BoxLib::Error("mixed_boundary_class::set_border_cache---only width 1 borders currently supported");
    if (flowdim >= 0 || flowdim == -3)
	BoxLib::Error("mixed_boundary_class::set_border_cache---negation borders not currently supported");
  
    assert(cache != 0);
    Real *const baseptr = cache->dataPtr();
    
    const Box& domain = lev_interface.domain();
    for (int iface = 0; iface < lev_interface.nfaces(); iface++) 
    {
	if (lev_interface.fgeo(iface) != level_interface::ALL)
	    break;
	int igrid = lev_interface.fgrid(iface, 0);
	int jgrid = lev_interface.fgrid(iface, 1);
	if (igrid < 0 || jgrid < 0) 
	{
	    const Box& b = lev_interface.node_face(iface);
	    int idim = lev_interface.fdim(iface);
#if (BL_SPACEDIM == 2)
	    int dstart, sstart, dstrid, sstrid, nvals;
#else
	    int dstart, sstart, dstrid1, dstrid2, sstrid1, sstrid2, nvals1, nvals2;
#endif
	    if (igrid < 0) 
	    {
		RegType t = ptr->bc[idim][0];
		if (t == refWall || t == inflow) 
		{
#if (BL_SPACEDIM == 2)
		    dstrid = r[jgrid].box().length(0);
		    if (idim == 0) 
		    {
			nvals = b.length(1) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
			    dstrid * (b.smallEnd(1) - w - r[jgrid].box().smallEnd(1));
			sstart = dstart + 2;
		    }
		    else 
		    {
			nvals = b.length(0) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
			    dstrid * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1));
			sstart = dstart + 2 * dstrid;
			dstrid = 1;
		    }
		    cache->set(iset++, dstart, sstart, dstrid, dstrid, nvals);
#else
		    dstrid1 = r[jgrid].box().length(0);
		    dstrid2 = dstrid1 * r[jgrid].box().length(1);
		    if (idim == 0) 
		    {
			nvals1 = b.length(1) + 2 * w;
			nvals2 = b.length(2) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
			    dstrid1 * (b.smallEnd(1) - w - r[jgrid].box().smallEnd(1)) +
			    dstrid2 * (b.smallEnd(2) - w - r[jgrid].box().smallEnd(2));
			sstart = dstart + 2;
		    }
		    else if (idim == 1) 
		    {
			nvals1 = b.length(0) + 2 * w;
			nvals2 = b.length(2) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
			    dstrid1 * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1)) +
			    dstrid2 * (b.smallEnd(2) - w - r[jgrid].box().smallEnd(2));
			sstart = dstart + 2 * dstrid1;
			dstrid1 = 1;
		    }
		    else 
		    {
			nvals1 = b.length(0) + 2 * w;
			nvals2 = b.length(1) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
			    dstrid1 * (b.smallEnd(1) - w - r[jgrid].box().smallEnd(1)) +
			    dstrid2 * (b.smallEnd(2) - 1 - r[jgrid].box().smallEnd(2));
			sstart = dstart + 2 * dstrid2;
			dstrid2 = dstrid1;
			dstrid1 = 1;
		    }
		    cache->set(iset++, dstart, sstart, dstrid1, dstrid2, dstrid1, dstrid2, nvals1, nvals2);
#endif
		}
		else if (t == periodic) 
		{
		    igrid = lev_interface.exterior_ref(igrid);
#if (BL_SPACEDIM == 2)
		    dstrid = r[jgrid].box().length(0);
		    sstrid = r[igrid].box().length(0);
		    if (idim == 0) 
		    {
			nvals = b.length(1) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
			    dstrid * (b.smallEnd(1) - w - r[jgrid].box().smallEnd(1));
			sstart = r[igrid].dataPtr() - baseptr +
			    (b.smallEnd(0) - 1 + domain.length(0) -
			    r[igrid].box().smallEnd(0)) +
			    sstrid * (b.smallEnd(1) - w - r[igrid].box().smallEnd(1));
		    }
		    else 
		    {
			nvals = b.length(0) + 2 * w;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
			    dstrid * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1));
			sstart = r[igrid].dataPtr() - baseptr +
			    b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
			    sstrid * (b.smallEnd(1) - 1 + domain.length(1) -
			    r[igrid].box().smallEnd(1));
			dstrid = 1;
			sstrid = 1;
		    }
		    cache->set(iset++, dstart, sstart, dstrid, sstrid, nvals);
#else
		    dstrid1 = r[jgrid].box().length(0);
		    dstrid2 = dstrid1 * r[jgrid].box().length(1);
		    sstrid1 = r[igrid].box().length(0);
		    sstrid2 = sstrid1 * r[igrid].box().length(1);
		    if (idim == 0) 
		    {
			int jl1 = 0, jh1 = 0, jl2 = 0, jh2 = 0;
			if (r.box(jgrid).smallEnd(1) == b.smallEnd(1))
			    jl1 = w;
			if (r.box(jgrid).bigEnd(1) == b.bigEnd(1))
			    jh1 = w;
			if (r.box(jgrid).smallEnd(2) == b.smallEnd(2))
			    jl2 = w;
			if (r.box(jgrid).bigEnd(2) == b.bigEnd(2))
			    jh2 = w;
			nvals1 = b.length(1) + jl1 + jh1;
			nvals2 = b.length(2) + jl2 + jh2;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
			    dstrid1 * (b.smallEnd(1) - jl1 - r[jgrid].box().smallEnd(1)) +
			    dstrid2 * (b.smallEnd(2) - jl2 - r[jgrid].box().smallEnd(2));
			sstart = r[igrid].dataPtr() - baseptr +
			    (b.smallEnd(0) - 1 + domain.length(0) -
			    r[igrid].box().smallEnd(0)) +
			    sstrid1 * (b.smallEnd(1) - jl1 - r[igrid].box().smallEnd(1)) +
			    sstrid2 * (b.smallEnd(2) - jl2 - r[igrid].box().smallEnd(2));
		    }
		    else if (idim == 1) 
		    {
			int jl1 = 0, jh1 = 0, jl2 = 0, jh2 = 0;
			if (r.box(jgrid).smallEnd(0) == b.smallEnd(0))
			    jl1 = w;
			if (r.box(jgrid).bigEnd(0) == b.bigEnd(0))
			    jh1 = w;
			if (r.box(jgrid).smallEnd(2) == b.smallEnd(2))
			    jl2 = w;
			if (r.box(jgrid).bigEnd(2) == b.bigEnd(2))
			    jh2 = w;
			nvals1 = b.length(0) + jl1 + jh1;
			nvals2 = b.length(2) + jl2 + jh2;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - jl1 - r[jgrid].box().smallEnd(0) +
			    dstrid1 * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1)) +
			    dstrid2 * (b.smallEnd(2) - jl2 - r[jgrid].box().smallEnd(2));
			sstart = r[igrid].dataPtr() - baseptr +
			    b.smallEnd(0) - jl1 - r[igrid].box().smallEnd(0) +
			    sstrid1 * (b.smallEnd(1) - 1 + domain.length(1) -
			    r[igrid].box().smallEnd(1)) +
			    sstrid2 * (b.smallEnd(2) - jl2 - r[igrid].box().smallEnd(2));
			dstrid1 = 1;
			sstrid1 = 1;
		    }
		    else 
		    {
			int jl1 = 0, jh1 = 0, jl2 = 0, jh2 = 0;
			if (r.box(jgrid).smallEnd(0) == b.smallEnd(0))
			    jl1 = w;
			if (r.box(jgrid).bigEnd(0) == b.bigEnd(0))
			    jh1 = w;
			if (r.box(jgrid).smallEnd(1) == b.smallEnd(1))
			    jl2 = w;
			if (r.box(jgrid).bigEnd(1) == b.bigEnd(1))
			    jh2 = w;
			nvals1 = b.length(0) + jl1 + jh1;
			nvals2 = b.length(1) + jl2 + jh2;
			dstart = r[jgrid].dataPtr() - baseptr +
			    b.smallEnd(0) - jl1 - r[jgrid].box().smallEnd(0) +
			    dstrid1 * (b.smallEnd(1) - jl2 - r[jgrid].box().smallEnd(1)) +
			    dstrid2 * (b.smallEnd(2) - 1 - r[jgrid].box().smallEnd(2));
			sstart = r[igrid].dataPtr() - baseptr +
			    b.smallEnd(0) - jl1 - r[igrid].box().smallEnd(0) +
			    sstrid1 * (b.smallEnd(1) - jl2 - r[igrid].box().smallEnd(1)) +
			    sstrid2 * (b.smallEnd(2) - 1 + domain.length(2) -
			    r[igrid].box().smallEnd(2));
			dstrid2 = dstrid1;
			sstrid2 = sstrid1;
			dstrid1 = 1;
			sstrid1 = 1;
		    }
		    cache->set(iset++, dstart, sstart, dstrid1, dstrid2, sstrid1, sstrid2, nvals1, nvals2);
#endif
	}
	else if (t == outflow) 
	{
	    // do nothing for node-based outflow
	}
	else 
	{
	    BoxLib::Error("mixed_boundary_class::set_boundary_cache---boundary type not supported");
	}
      }
      else if (jgrid < 0) 
      {
	  RegType t = ptr->bc[idim][1];
	  if (t == refWall || t == inflow) 
	  {
#if (BL_SPACEDIM == 2)
	      dstrid = r[igrid].box().length(0);
	      if (lev_interface.fdim(iface) == 0) 
	      {
		  nvals = b.length(1) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
		      dstrid * (b.smallEnd(1) - w - r[igrid].box().smallEnd(1));
		  sstart = dstart - 2;
	      }
	      else 
	      {
		  nvals = b.length(0) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
		      dstrid * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1));
		  sstart = dstart - 2 * dstrid;
		  dstrid = 1;
	      }
	      cache->set(iset++, dstart, sstart, dstrid, dstrid, nvals);
#else
	      dstrid1 = r[igrid].box().length(0);
	      dstrid2 = dstrid1 * r[igrid].box().length(1);
	      if (idim == 0) 
	      {
		  nvals1 = b.length(1) + 2 * w;
		  nvals2 = b.length(2) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
		      dstrid1 * (b.smallEnd(1) - w - r[igrid].box().smallEnd(1)) +
		      dstrid2 * (b.smallEnd(2) - w - r[igrid].box().smallEnd(2));
		  sstart = dstart - 2;
	      }
	      else if (idim == 1) 
	      {
		  nvals1 = b.length(0) + 2 * w;
		  nvals2 = b.length(2) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
		      dstrid1 * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1)) +
		      dstrid2 * (b.smallEnd(2) - w - r[igrid].box().smallEnd(2));
		  sstart = dstart - 2 * dstrid1;
		  dstrid1 = 1;
	      }
	      else 
	      {
		  nvals1 = b.length(0) + 2 * w;
		  nvals2 = b.length(1) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
		      dstrid1 * (b.smallEnd(1) - w - r[igrid].box().smallEnd(1)) +
		      dstrid2 * (b.smallEnd(2) + 1 - r[igrid].box().smallEnd(2));
		  sstart = dstart - 2 * dstrid2;
		  dstrid2 = dstrid1;
		  dstrid1 = 1;
	      }
	      cache->set(iset++, dstart, sstart,
		  dstrid1, dstrid2, dstrid1, dstrid2, nvals1, nvals2);
#endif
	  }
	  else if (t == periodic) 
	  {
	      jgrid = lev_interface.exterior_ref(jgrid);
#if (BL_SPACEDIM == 2)
	      dstrid = r[igrid].box().length(0);
	      sstrid = r[jgrid].box().length(0);
	      if (lev_interface.fdim(iface) == 0) 
	      {
		  nvals = b.length(1) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
		      dstrid * (b.smallEnd(1) - w - r[igrid].box().smallEnd(1));
		  sstart = r[jgrid].dataPtr() - baseptr +
		      (b.smallEnd(0) + 1 - domain.length(0) -
		      r[jgrid].box().smallEnd(0)) +
		      sstrid * (b.smallEnd(1) - w - r[jgrid].box().smallEnd(1));
	      }
	      else 
	      {
		  nvals = b.length(0) + 2 * w;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
		      dstrid * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1));
		  sstart = r[jgrid].dataPtr() - baseptr +
		      b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
		      sstrid * (b.smallEnd(1) + 1 - domain.length(1) -
		      r[jgrid].box().smallEnd(1));
		  dstrid = 1;
		  sstrid = 1;
	      }
	      cache->set(iset++, dstart, sstart, dstrid, sstrid, nvals);
#else
	      dstrid1 = r[igrid].box().length(0);
	      dstrid2 = dstrid1 * r[igrid].box().length(1);
	      sstrid1 = r[jgrid].box().length(0);
	      sstrid2 = sstrid1 * r[jgrid].box().length(1);
	      if (idim == 0) 
	      {
		  int il1 = 0, ih1 = 0, il2 = 0, ih2 = 0;
		  if (r.box(igrid).smallEnd(1) == b.smallEnd(1))
		      il1 = w;
		  if (r.box(igrid).bigEnd(1) == b.bigEnd(1))
		      ih1 = w;
		  if (r.box(igrid).smallEnd(2) == b.smallEnd(2))
		      il2 = w;
		  if (r.box(igrid).bigEnd(2) == b.bigEnd(2))
		      ih2 = w;
		  nvals1 = b.length(1) + il1 + ih1;
		  nvals2 = b.length(2) + il2 + ih2;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
		      dstrid1 * (b.smallEnd(1) - il1 - r[igrid].box().smallEnd(1)) +
		      dstrid2 * (b.smallEnd(2) - il2 - r[igrid].box().smallEnd(2));
		  sstart = r[jgrid].dataPtr() - baseptr +
		      (b.smallEnd(0) + 1 - domain.length(0) -
		      r[jgrid].box().smallEnd(0)) +
		      sstrid1 * (b.smallEnd(1) - il1 - r[jgrid].box().smallEnd(1)) +
		      sstrid2 * (b.smallEnd(2) - il2 - r[jgrid].box().smallEnd(2));
	      }
	      else if (idim == 1) 
	      {
		  int il1 = 0, ih1 = 0, il2 = 0, ih2 = 0;
		  if (r.box(igrid).smallEnd(0) == b.smallEnd(0))
		      il1 = w;
		  if (r.box(igrid).bigEnd(0) == b.bigEnd(0))
		      ih1 = w;
		  if (r.box(igrid).smallEnd(2) == b.smallEnd(2))
		      il2 = w;
		  if (r.box(igrid).bigEnd(2) == b.bigEnd(2))
		      ih2 = w;
		  nvals1 = b.length(0) + il1 + ih1;
		  nvals2 = b.length(2) + il2 + ih2;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) - il1 - r[igrid].box().smallEnd(0) +
		      dstrid1 * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1)) +
		      dstrid2 * (b.smallEnd(2) - il2 - r[igrid].box().smallEnd(2));
		  sstart = r[jgrid].dataPtr() - baseptr +
		      b.smallEnd(0) - il1 - r[jgrid].box().smallEnd(0) +
		      sstrid1 * (b.smallEnd(1) + 1 - domain.length(1) -
		      r[jgrid].box().smallEnd(1)) +
		      sstrid2 * (b.smallEnd(2) - il2 - r[jgrid].box().smallEnd(2));
		  dstrid1 = 1;
		  sstrid1 = 1;
	      }
	      else 
	      {
		  int il1 = 0, ih1 = 0, il2 = 0, ih2 = 0;
		  if (r.box(igrid).smallEnd(0) == b.smallEnd(0))
		      il1 = w;
		  if (r.box(igrid).bigEnd(0) == b.bigEnd(0))
		      ih1 = w;
		  if (r.box(igrid).smallEnd(1) == b.smallEnd(1))
		      il2 = w;
		  if (r.box(igrid).bigEnd(1) == b.bigEnd(1))
		      ih2 = w;
		  nvals1 = b.length(0) + il1 + ih1;
		  nvals2 = b.length(1) + il2 + ih2;
		  dstart = r[igrid].dataPtr() - baseptr +
		      b.smallEnd(0) - il1 - r[igrid].box().smallEnd(0) +
		      dstrid1 * (b.smallEnd(1) - il2 - r[igrid].box().smallEnd(1)) +
		      dstrid2 * (b.smallEnd(2) + 1 - r[igrid].box().smallEnd(2));
		  sstart = r[jgrid].dataPtr() - baseptr +
		      b.smallEnd(0) - il1 - r[jgrid].box().smallEnd(0) +
		      sstrid1 * (b.smallEnd(1) - il2 - r[jgrid].box().smallEnd(1)) +
		      sstrid2 * (b.smallEnd(2) + 1 - domain.length(2) -
		      r[jgrid].box().smallEnd(2));
		  dstrid2 = dstrid1;
		  sstrid2 = sstrid1;
		  dstrid1 = 1;
		  sstrid1 = 1;
	      }
	      cache->set(iset++, dstart, sstart, dstrid1, dstrid2, sstrid1, sstrid2, nvals1, nvals2);
#endif
	}
	else if (t == outflow) 
	{
	    // do nothing for node-based outflow
	}
	else 
	{
	    BoxLib::Error("mixed_boundary_class::set_boundary_cache---boundary type not supported");
	}
      }
    }
  }
  if (iset > nsets)
      BoxLib::Error("mixed_boundary_class::set_boundary_cache---too many boundary edges to cache");
}

void amr_boundary_class::boundary_mesh(BoxArray& exterior_mesh,
				       int *&grid_ref,
				       const BoxArray& interior_mesh,
				       const Box& domain) const
{
    BoxList bl;
    List<int> il;
    const Box &d = domain;
    for (int igrid = 0; igrid < interior_mesh.length(); igrid++) 
    {
	check_against_boundary(bl, il, interior_mesh[igrid], igrid, d, 0);
    }
    exterior_mesh.define(bl);
    bl.clear();
    
    if (il.length() != exterior_mesh.length())
	BoxLib::Error("amr_boundary_class::boundary_mesh---List<int> wrong length");
    
    grid_ref = new int[exterior_mesh.length()];
    ListIterator<int> in(il);
    for (int igrid = 0; in; in++, igrid++) 
    {
	grid_ref[igrid] = in();
    }
    il.clear();
}

void amr_boundary_class::check_against_boundary(BoxList& bl, List<int>& il,
						const Box& b, int ib,
						const Box& d, int dim1) const
{
    for (int i = dim1; i < BL_SPACEDIM; i++) 
    {
	if (b.smallEnd(i) == d.smallEnd(i)) 
	{
	    Box bn = adjCellLo(b,i);
	    bl.append(bn);
	    il.append(ib);
	    check_against_boundary(bl, il, bn, ib, d, i+1);
	}
	if (b.bigEnd(i) == d.bigEnd(i)) 
	{
	    Box bn = adjCellHi(b,i);
	    bl.append(bn);
	    il.append(ib);
	    check_against_boundary(bl, il, bn, ib, d, i+1);
	}
    }
}

void mixed_boundary_class::check_against_boundary(BoxList& bl, List<int>& il,
						  const Box& b, int ib,
						  const Box& d, int dim1) const
{
    for (int i = dim1; i < BL_SPACEDIM; i++) 
    {
	if (b.smallEnd(i) == d.smallEnd(i)) 
	{
	    if (ptr->bc[i][0] == refWall || ptr->bc[i][0] == inflow) 
	    {
		Box bn = b;
		bn.shift(i,-b.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][0] == periodic) 
	    {
		Box bn = b;
		bn.shift(i,d.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][0] == outflow) 
	    {
		Box bn = b;
		bn.shift(i,-b.length(i));
		bl.append(bn);
		il.append(-2);
		check_against_boundary(bl, il, bn, -1, d, i+1);
	    }
	    else 
	    {
		BoxLib::Error("mixed_boundary_class::check_against_boundary()  Boundary type not supported");
	    }
	}
	if (b.bigEnd(i) == d.bigEnd(i)) 
	{
	    if (ptr->bc[i][1] == refWall ||
		ptr->bc[i][1] == inflow) 
	    {
		Box bn = b;
		bn.shift(i,b.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][1] == periodic) 
	    {
		Box bn = b;
		bn.shift(i,-d.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][1] == outflow) 
	    {
		Box bn = b;
		bn.shift(i,b.length(i));
		bl.append(bn);
		il.append(-2);
		check_against_boundary(bl, il, bn, -1, d, i+1);
	    }
	    else 
	    {
		BoxLib::Error("mixed_boundary_class::check_against_boundary()  Boundary type not supported");
	    }
	}
    }
}

void mixed_boundary_class::duplicate(List<Box>& bl, const Box& domain) const
{
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	if (ptr->bc[i][0] == periodic) 
	{
	    for ( ListIterator<Box> bn(bl.last()); bn; bn--) 
	    {
		if (bn().type(i) == IndexType::NODE) 
		{
		    if (bn().smallEnd(i) == domain.smallEnd(i)) 
		    {
			Box btmp = bn();
			btmp.shift(i, domain.length(i));
			if(!bl.includes(btmp))
			    bl.append(btmp);
		    }
		    else if (bn().bigEnd(i) - 1 == domain.bigEnd(i)) 
		    {
			Box btmp = bn();
			btmp.shift(i, -domain.length(i));
			if(!bl.includes(btmp))
			    bl.append(btmp);
		    }
		}
	    }
	}
    }
}

bool mixed_boundary_class::singular() const
{
    bool retval = true;
    if (flowdim == -2) 
    {
	for (int idim = 0; idim < BL_SPACEDIM; idim++) 
	{
	    if (ptr->bc[idim][0] == outflow || ptr->bc[idim][1] == outflow)
		retval = false;
	}
    }
    else 
    {
	BoxLib::Error("mixed_boundary_class::singular---only defined for pressure boundary");
    }
    return retval;
}

amr_fluid_boundary_class::amr_fluid_boundary_class()
{
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	v[i] = 0;
    }
    s = 0;
    p = 0;
    ts = 0;
}

inviscid_fluid_boundary_class::inviscid_fluid_boundary_class(RegType Bc[BL_SPACEDIM][2])
{
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	bc[i][0] = Bc[i][0];
	bc[i][1] = Bc[i][1];
	if ((bc[i][0] == periodic || bc[i][1] == periodic) && bc[i][1] != bc[i][0])
	    BoxLib::Error("inviscid_fluid_boundary_class::ctor---periodic bc's don't match");
	v[i] = new mixed_boundary_class(this, i);
    }
    s = new mixed_boundary_class(this, -1);
    p = new mixed_boundary_class(this, -2);
    ts = new mixed_boundary_class(this, -4);
}

inviscid_fluid_boundary_class::~inviscid_fluid_boundary_class()
{
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {
	delete v[i];
    }
    delete s;
    delete p;
    delete ts;
}
