//BL_COPYRIGHT_NOTICE

#include "boundary.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define FORT_FBREF    bref_
#define FORT_FBREFM   brefm_
#define FORT_FBNEG    bneg_
#define FORT_FBNEGM   bnegm_
#define FORT_FBINFLO  binflo_
#define FORT_FBINFIL  binfil_
#elif defined( BL_FORT_USE_UPPERCASE )
#define FORT_FBREF    BREF
#define FORT_FBREFM   BREFM
#define FORT_FBNEG    BNEG
#define FORT_FBNEGM   BNEGM
#define FORT_FBINFLO  BINFLO
#define FORT_FBINFIL  BINFIL
#elif defined( BL_FORT_USE_LOWERCASE )
#define FORT_FBREF    bref
#define FORT_FBREFM   brefm
#define FORT_FBNEG    bneg
#define FORT_FBNEGM   bnegm
#define FORT_FBINFLO  binflo
#define FORT_FBINFIL  binfil
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_FBREF  (Real*, intS, intS, const Real*, intS, intS,
		      const int*, const int&);
    void FORT_FBREFM (Real*, intS, intS, const Real*, intS, intS,
		      const int*, const int&);
    void FORT_FBNEG  (Real*, intS, intS, const Real*, intS, intS,
		      const int*, const int&);
    void FORT_FBNEGM (Real*, intS, intS, const Real*, intS, intS,
		      const int*, const int&);
    void FORT_FBINFLO(Real*, intS, intS, const Real*, intS, intS,
		      const int*, const int&);
    void FORT_FBINFIL(Real*, intS, intS, const Real*, intS, intS,
		      const int*, const int&);
}

amr_boundary::~amr_boundary () {}

int
amr_boundary::dir (const Box& region,
		   const Box& domain) const
{
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
	if (region.bigEnd(i) < domain.smallEnd(i))
	    return -(i+1);
	if (region.smallEnd(i) > domain.bigEnd(i))
	    return +(i+1);
    }
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
	if (region.bigEnd(i) == domain.smallEnd(i))
	    return -(i+1);
	if (region.smallEnd(i) == domain.bigEnd(i))
	    return +(i+1);
    }
    BL_ASSERT("amr_boundary::dir---boundary box not outside domain." == 0);
    return 0;
}

void
amr_boundary::boundary_mesh (BoxArray&       exterior_mesh,
			     int*&           grid_ref,
			     const BoxArray& interior_mesh,
			     const Box&      domain) const
{
    BoxList bl;
    List<int> il;
    const Box& d = domain;
    for (int igrid = 0; igrid < interior_mesh.length(); igrid++)
    {
	check_against_boundary_(bl, il, interior_mesh[igrid], igrid, d, 0);
    }
    exterior_mesh.define(bl);
    bl.clear();

    BL_ASSERT(il.length() == exterior_mesh.length());

    grid_ref = new int[exterior_mesh.length()];
    ListIterator<int> in(il);
    for (int igrid = 0; in; in++, igrid++)
    {
	grid_ref[igrid] = in();
    }
    il.clear();
}


//
mixed_boundary::mixed_boundary (inviscid_fluid_boundary* Ptr,
				int                      idim)
    :
    ptr(Ptr),
    flowdim(idim)
{
}

Box
mixed_boundary::box_ (const Box& image,
		      const Box& domain,
		      int        idir) const
{
    HG_DEBUG_OUT("BOX_:"
		 << "image(" << image << ") "
		 << "domain(" << domain << ") "
		 << "idir(" << idir << ") ");
    const int idim = abs(idir) - 1;
    const RegType t = ptr->bc[idim][idir > 0];
    Box retbox(image);

    BL_ASSERT(idir != 0);
    if (t == refWall || t == outflow || t == inflow)
    {
        //
	// All these cases use a reflected box.
        //
	if (idir < 0)
	{
	    if (image.type(idim) == IndexType::CELL)
	    {
		retbox.shift(
		    idim,
		    (2 * domain.smallEnd(idim) - 1
		     - image.bigEnd(idim) - image.smallEnd(idim)));
		if ( t == inflow && idim == flowdim )
		{
		    retbox.setSmall(idim, domain.smallEnd(idim) - 1);
		}
	    }
	    else
	    {
		retbox.shift(
		    idim,
		    (2 * domain.smallEnd(idim)
		     - image.bigEnd(idim) - image.smallEnd(idim)));
		if ( t == inflow && idim == flowdim )
		{
		    retbox.setSmall(idim, domain.smallEnd(idim));
		}
	    }
	}
	else if (idir > 0)
	{
	    if (image.type(idim) == IndexType::CELL)
	    {
		retbox.shift(
		    idim,
		    (2 * domain.bigEnd(idim) + 1
		     - image.bigEnd(idim) - image.smallEnd(idim)));
	    }
	    else
	    {
		retbox.shift(
		    idim,
		    (2 * domain.bigEnd(idim) + 2
		     - image.bigEnd(idim) - image.smallEnd(idim)));
	    }
	    if ( t == inflow && idim == flowdim )
	    {
		retbox.setBig(idim, domain.bigEnd(idim) + 1);
	    }
	}
    }
    else if (t == periodic)
    {
	if (idir < 0)
	{
	    retbox.shift(idim, domain.length(idim));
	}
	else if (idir > 0)
	{
	    retbox.shift(idim, -domain.length(idim));
	}
    }
    else
    {
	BoxLib::Error( "mixed_boundary::box---boundary type not supported");
    }
    HG_DEBUG_OUT("==>retbox(" << retbox << ")"
		 << endl);
    return retbox;
}

Box
mixed_boundary::anImage (const Box& region,
			 const Box& srcbox,
			 const Box& domain) const
{
    HG_DEBUG_OUT("anImage:"
		 << "region(" << region << ") "
		 << "srcbox(" << srcbox << ") "
		 << "domain(" << domain << ") ");
    Box tdomain = domain;
    tdomain.convert(srcbox.type());
    Box idomain = ::grow(tdomain, IntVect::TheZeroVector()-srcbox.type());
    Box image   = region;

    for (int idim = 0; idim < BL_SPACEDIM; idim++)
    {
	if (region.bigEnd(idim) < idomain.smallEnd(idim))
	{
	    const RegType t = ptr->bc[idim][0];

	    if ( t == inflow  && idim == flowdim )
	    {
	    }
	    else if (t == refWall || t == inflow || t == outflow)
	    {
		image.shift(
		    idim,
		    (tdomain.smallEnd(idim) + idomain.smallEnd(idim) - 1
		     - region.bigEnd(idim) - region.smallEnd(idim)));
	    }
	    else if (t == periodic)
	    {
		image.shift(idim, domain.length(idim));
	    }
	}
	else if (region.smallEnd(idim) > idomain.bigEnd(idim))
	{
	    const RegType t = ptr->bc[idim][1];

	    if ( t == refWall && idim == flowdim )
	    {
	    }
	    if (t == refWall || t == inflow || t == outflow)
	    {
		image.shift(
		    idim,
		    (tdomain.bigEnd(idim) + idomain.bigEnd(idim) + 1
		     - region.bigEnd(idim) - region.smallEnd(idim)));
	    }
	    else if (t == periodic)
	    {
		image.shift(idim, -domain.length(idim));
	    }
	}
    }

    BL_ASSERT(image.type() == srcbox.type());

    HG_DEBUG_OUT("==>image(" << image << ")"
		 << endl);
    return image;
}

//
// Reflects on all outflow cases (which aren't called anyway).
// On velocity inflow, uses box function which extends interior
// box just past edge of domain.
//

void
mixed_boundary::fill (FArrayBox&       patch,
		      const Box&       region,
		      const FArrayBox& src,
		      const Box&       domain) const
{
    // BL_ASSERT(domain.type() == IntVect::TheZeroVector());
    HG_DEBUG_OUT("FILL: "
		 << "patch.box(" << patch.box() << ") "
		 << "region(" << region << ") "
		 << "src.box(" << src.box() << ") "
		 << "domain(" << domain << ") "
		 <<  endl);
    Box tdomain = domain;
    tdomain.convert(type(src));
    Box idomain = ::grow(tdomain, IntVect::TheZeroVector()-type(src));
    Box img     = anImage(region, src.box(), domain);
    int idir = 0;
    int refarray[BL_SPACEDIM] = {0};
    bool negflag = true;
    bool negarray[BL_SPACEDIM-1];

    for (int i = 0; i < BL_SPACEDIM - 1; i++)
    {
	negarray[i] = true;
    }

    for (int idim = 0; idim < BL_SPACEDIM; idim++)
    {
	if (region.bigEnd(idim) < idomain.smallEnd(idim))
	{
	    const RegType t = ptr->bc[idim][0];

	    if (t == inflow && idim == flowdim)
	    {
		idir = -1 - idim;
	    }
	    else if (t == refWall || t == inflow || t == outflow)
	    {
		refarray[idim] = 1;

		if (flowdim == -3 || t == refWall && idim == flowdim)
		{
		    negflag = !negflag;
		}
		if (flowdim == -4)
		{
		    if (idim < BL_SPACEDIM - 1)
		    {
			negarray[idim] = !negarray[idim];
		    }
		    else
		    {
			for (int i = 0; i < BL_SPACEDIM - 1; i++)
			{
			    negarray[i] = !negarray[i];
			}
		    }
		}
	    }
	}
	else if (region.smallEnd(idim) > idomain.bigEnd(idim))
	{
	    const RegType t = ptr->bc[idim][1];

	    if (t == inflow && idim == flowdim)
	    {
		idir = 1 + idim;
	    }
	    if (t == refWall || t == inflow || t == outflow)
	    {
		refarray[idim] = 1;

		if (flowdim == -3 || t == refWall && idim == flowdim)
		{
		    negflag = !negflag;
		}
		if (flowdim == -4)
		{
		    if (idim < BL_SPACEDIM - 1)
		    {
			negarray[idim] = !negarray[idim];
		    }
		    else
		    {
			for (int i = 0; i < BL_SPACEDIM - 1; i++)
			{
			    negarray[i] = !negarray[i];
			}
		    }
		}
	    }
	}
    }

    HG_DEBUG_OUT("idir = " << idir
		 << " flowdim = " << flowdim
		 << " reg = " << region
		 << " img = " << img
		 << " src.box = " << src.box()
		 << endl );

    if (idir != 0)
    {
        //
	// Normal-component inflow section, assume patch.nComp() == 1
        //
	Box bb = box_(img, domain, idir);

	if (img == region)
        {
            //
	    // Only bdy involved, can fill directly from interior
            //
	    HG_DEBUG_OUT("IMG==REGION" << img << region << endl);
	    fill_(patch, region, src, bb, domain, idir);
        }
	else
	{
            //
	    // Multiple bdys, fill intermediate patch.
            //
	    FArrayBox gb(img);
	    fill_(gb, img, src, bb, domain, idir);
	    HG_DEBUG_OUT("IMG!=REGION negflag" << negflag << img << region << endl);
	    if (negflag)
	    {
		FORT_FBREFM(patch.dataPtr(), DIMLIST(patch.box()),
                            DIMLIST(region),
                            gb.dataPtr(), DIMLIST(img), DIMLIST(img),
                            refarray, 1);
	    }
	    else
	    {
		FORT_FBNEGM(patch.dataPtr(), DIMLIST(patch.box()),
                            DIMLIST(region),
                            gb.dataPtr(), DIMLIST(img), DIMLIST(img),
                            refarray, 1);
	    }
	}
    }
    else if (flowdim == -4)
    {
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    FORT_FBREFM(patch.dataPtr(i), DIMLIST(patch.box()),
                        DIMLIST(region),
                        src.dataPtr(i), DIMLIST(src.box()), DIMLIST(img),
                        refarray, 1);
	}
	for (int idim = 0; idim < BL_SPACEDIM - 1; idim++)
	{
	    const int i = idim + BL_SPACEDIM;

	    if (negarray[idim])
	    {
		FORT_FBREFM(patch.dataPtr(i), DIMLIST(patch.box()),
                            DIMLIST(region),
                            src.dataPtr(i), DIMLIST(src.box()), DIMLIST(img),
                            refarray, 1);
	    }
	    else
	    {
		FORT_FBNEGM(patch.dataPtr(i), DIMLIST(patch.box()),
                            DIMLIST(region),
                            src.dataPtr(i), DIMLIST(src.box()), DIMLIST(img),
                            refarray, 1);
	    }
	}
    }
    else
    {
        //
	// All cases other than normal-component inflow.
        //
        if (src.box().contains(img))  // FIXME!!!
        {
            if (negflag)
            {
		HG_DEBUG_OUT("GOT HERE negflag\n");
                FORT_FBREFM(patch.dataPtr(), DIMLIST(patch.box()),
                            DIMLIST(region),
                            src.dataPtr(), DIMLIST(src.box()), DIMLIST(img),
                            refarray, patch.nComp());
            }
            else
            {
		HG_DEBUG_OUT("GOT HERE not negflag\n");
                FORT_FBNEGM(patch.dataPtr(), DIMLIST(patch.box()),
                            DIMLIST(region),
                            src.dataPtr(), DIMLIST(src.box()), DIMLIST(img),
                            refarray, patch.nComp());
            }
        }
    }
}

void
mixed_boundary::fill_ (FArrayBox&       patch,
		       const Box&       region,
		       const FArrayBox& bgr,
		       const Box&       bb,
		       const Box&       domain,
		       int              idir) const
{
    const int idim = abs(idir) - 1;
    const RegType t = ptr->bc[idim][idir > 0];

    if (flowdim == -4 && (t == refWall || t == inflow))
    {
        //
	// Terrain sigma.
        //
	BoxLib::Abort( "mixed_boundary::fill(): terrain undefined" );
    }
    else if (t == refWall)
    {
	if (idim == flowdim || flowdim == -3)
	{
	    FORT_FBNEG(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
	}
	else
	{
	    FORT_FBREF(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
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
	    FORT_FBREF(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
	}
	else if (flowdim == -1)
	{
            //
	    //BoxLib::Error("mixed_boundary::Don't know how to do inflow density");
	    // Inflow density---just reflect interior for now
            //
	    FORT_FBREF(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
	}
	else if (flowdim == -3)
	{
	    FORT_FBNEG(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
	}
	else if (idim == flowdim)
	{
            //
	    // For this to work, fill_borders must already have been called
	    // to initialize values in the first ghost cell outside the domain.
            //
	    HG_DEBUG_OUT("fill:"
			 << "patch.box(" << patch.box() << ") "
			 << "region(" << region << ") "
			 << "bgr.box(" << bgr.box() << ") "
			 << "bb(" << bb << ") "
			 << endl);
	    FORT_FBINFIL(patch.dataPtr(), DIMLIST(patch.box()),
			 DIMLIST(region),
			 bgr.dataPtr(), DIMLIST(bgr.box()),
			 DIMLIST(bb), &idim, 1);
	}
	else if (flowdim >= 0)
	{
            //
	    // transverse velocity components
	    // patch.assign(0.0, region);
	    // we now believe this looks like a refWall to transverse components
            //
	    FORT_FBREF(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
	}
    }
    else if (t == outflow)
    {
        //
	// Do nothing if NODE-based, reflect if CELL-based.
        //
	if (type(patch, idim) == IndexType::CELL)
	{
	    FORT_FBREF(patch.dataPtr(), DIMLIST(patch.box()),
		       DIMLIST(region),
		       bgr.dataPtr(), DIMLIST(bgr.box()),
		       DIMLIST(bb), &idim, patch.nComp());
	}
    }
    else
    {
	BoxLib::Abort( "mixed_boundary::fill(): boundary type not supported" );
    }
}

void
mixed_boundary::sync_borders (MultiFab&              r,
			      const level_interface& lev_interface) const
{
    BL_ASSERT(type(r) == IntVect::TheNodeVector());

    task_list tl;
    // we are looping over only the fine-fine faces
    for (int iface = 0;
	 iface < lev_interface.nboxes(level_interface::FACEDIM); iface++)
    {
	if (lev_interface.geo(level_interface::FACEDIM, iface)
	    != level_interface::ALL)
	    break;
        int igrid = lev_interface.grid(level_interface::FACEDIM, iface, 0);
	if (igrid < 0)
	{
	    const int idim = lev_interface.fdim(iface);
	    if (ptr->bc[idim][0] == periodic)
	    {
		const int jgrid =
		    lev_interface.grid(level_interface::FACEDIM, iface, 1);
		igrid = lev_interface.exterior_ref(igrid);
		const Box& b =
		    lev_interface.node_box(level_interface::FACEDIM, iface);
		Box bb = b;
		bb.shift(idim, lev_interface.domain().length(idim));
		tl.add_task(new task_copy(tl, r, jgrid, b, r, igrid, bb));
	    }
	}
    }
    tl.execute("mixed_boundary::sync_borders");
}

void
mixed_boundary::fill_borders (MultiFab&              r,
			      const level_interface& lev_interface,
			      int                    w) const
{
    w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;

    BL_ASSERT(w == 1 || w == 0);

    const Box& domain = lev_interface.domain();

    task_list tl;
    // we are looping over only the fine-fine faces
    for (int iface = 0;
	 iface < lev_interface.nboxes(level_interface::FACEDIM); iface++)
    {
	if (lev_interface.geo(level_interface::FACEDIM, iface)
	    != level_interface::ALL)
	    break;

	const int igrid =
	    lev_interface.grid(level_interface::FACEDIM, iface, 0);
	const int jgrid =
	    lev_interface.grid(level_interface::FACEDIM, iface, 1);

	if (igrid < 0 || jgrid < 0)
	{
	    Box b = lev_interface.box(level_interface::FACEDIM, iface);
	    const int idim = lev_interface.fdim(iface);
	    const int a = (type(r, idim) == IndexType::NODE);
            //
	    // Need to do on x borders too in case y border is an interior face
            //
	    if (igrid < 0)
	    {
		for (int i = 0; i < BL_SPACEDIM; i++)
		{
		    if (i != idim)
		    {
			if (lev_interface.interior_mesh()[jgrid].smallEnd(i)
			    == b.smallEnd(i))
			{
			    b.growLo(i, w);
			}
			if (lev_interface.interior_mesh()[jgrid].bigEnd(i)
			    == b.bigEnd(i))
			{
			    b.growHi(i, w);
			}
		    }
		}
		b.shift(idim, -a).growLo(idim, w-a).convert(type(r));
		const RegType t = ptr->bc[idim][0];
		Box bb = b;
		if (flowdim == -4 && (t == refWall || t == inflow))
		{
                    //
		    // Terrain sigma
                    //
		    if (is_remote(r, jgrid)) continue;
		    bb.shift(
			idim,
			2 * domain.smallEnd(idim) - 1 + a
			- b.bigEnd(idim) - b.smallEnd(idim));
		    const Box& rbox = r[jgrid].box();
		    for (int i = 0; i < r.nComp(); i++)
		    {
			Real* rptr = r[jgrid].dataPtr(i);
			if ((i == idim + BL_SPACEDIM)
			    || (i >= BL_SPACEDIM && idim == BL_SPACEDIM - 1))
			{
			    FORT_FBNEG(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
			else
			{
			    FORT_FBREF(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
		    }
		}
		else if (t == refWall)
		{
		    if ( is_remote(r, jgrid) ) continue;
		    bb.shift(
			idim,
			(2 * domain.smallEnd(idim) - 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[jgrid].box();
		    if (idim == flowdim || flowdim == -3)
		    {
			FORT_FBNEG(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else
		    {
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
		else if (t == periodic)
		{
		    int isrc = lev_interface.exterior_ref(igrid);
		    bb.shift(idim, domain.length(idim));
		    //r[jgrid].copy(r[isrc], bb, 0, b, 0, r.nComp());
		    tl.add_task(new task_copy(tl, r, jgrid, b, r, isrc, bb));
		}
		else if (t == inflow)
		{
		    if ( is_remote(r, jgrid) ) continue;
		    bb.shift(
			idim,
			(2 * domain.smallEnd(idim) - 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[jgrid].box();
		    if (flowdim == -2)
		    {
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (flowdim == -1)
		    {
                        //
			//BoxLib::Error("mixed_boundary::Don't know how to do inflow density");
			// Inflow density---just reflect interior for now
                        //
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (flowdim == -3)
		    {
			FORT_FBNEG(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, 1);
		    }
		    else if (idim == flowdim)
		    {
                        //
			// For this to work, fill_borders must be called exactly
			// once for each level of this variable.
                        //
			FORT_FBINFLO(r[jgrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(b),
				     r[jgrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(bb), &idim, 1);
		    }
		    else if (flowdim >= 0)
		    {
                        //
			// transverse velocity components
			//r[jgrid].assign(0.0, b);
			// we now believe this looks like a refWall to transverse comps
                        //
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, 1);
		    }
		}
		else if (t == outflow)
		{
                    //
		    // Do nothing if NODE-based, reflect if CELL-based.
                    //
		    if (is_remote(r, jgrid))
                        continue;
		    if (type(r, idim) == IndexType::CELL)
		    {
			bb.shift(
			    idim,
			    (2 * domain.smallEnd(idim) - 1 + a
			     - b.bigEnd(idim) - b.smallEnd(idim)));
			const Box& rbox = r[jgrid].box();
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
	    }
	    else if (jgrid < 0)
	    {
		for (int i = 0; i < BL_SPACEDIM; i++)
		{
		    if (i != idim)
		    {
			if (lev_interface.interior_mesh()[igrid].smallEnd(i)
			    == b.smallEnd(i))
			{
			    b.growLo(i, w);
			}
			if (lev_interface.interior_mesh()[igrid].bigEnd(i) == b.bigEnd(i))
			{
			    b.growHi(i, w);
			}
		    }
		}
		b.shift(idim, a).growHi(idim, w-a).convert(type(r));
		const RegType t = ptr->bc[idim][1];
		Box bb = b;
		if (flowdim == -4 && (t == refWall || t == inflow))
		{
		    if (is_remote(r, igrid))
                        continue;
                    //
		    // Terrain sigma
                    //
		    bb.shift(
			idim,
			(2 * domain.bigEnd(idim) + 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[igrid].box();
		    for (int i = 0; i < r.nComp(); i++)
		    {
			Real* rptr = r[igrid].dataPtr(i);
			if ((i == idim + BL_SPACEDIM)
			    || (i >= BL_SPACEDIM && idim == BL_SPACEDIM - 1))
			{
			    FORT_FBNEG(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
			else
			{
			    FORT_FBREF(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
		    }
		}
		else if (t == refWall)
		{
		    if (is_remote(r, igrid))
                        continue;
		    bb.shift(
			idim,
			(2 * domain.bigEnd(idim) + 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[igrid].box();
		    if (idim == flowdim || flowdim == -3)
		    {
			FORT_FBNEG(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else
		    {
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
		else if (t == periodic)
		{
		    int isrc = lev_interface.exterior_ref(jgrid);
		    bb.shift(idim, -domain.length(idim));
		    //r[igrid].copy(r[isrc], bb, 0, b, 0, r.nComp());
		    tl.add_task(new task_copy(tl, r, igrid, b, r, isrc, bb));
		}
		else if (t == inflow)
		{
		    if (is_remote(r, igrid)) continue;
		    bb.shift(
			idim,
			(2 * domain.bigEnd(idim) + 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[igrid].box();
		    if (flowdim == -2)
		    {
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(),
				   DIMLIST(rbox), DIMLIST(bb), &idim, 1);
		    }
		    else if (flowdim == -1)
		    {
                        //
			//BoxLib::Error("mixed_boundary::Don't know how to do inflow density");
			// Inflow density---just reflect interior for now
                        //
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (flowdim == -3)
		    {
			FORT_FBNEG(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (idim == flowdim)
		    {
                        //
			// For this to work, fill_borders must be called exactly
			// once for each level of this variable.
                        //
			FORT_FBINFLO(r[igrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(b),
				     r[igrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(bb), &idim, 1);
		    }
		    else if (flowdim >= 0)
		    {
                        //
			// transverse velocity components
			//r[igrid].assign(0.0, rbox);
			// we now believe this looks like a refWall to transverse comps
                        //
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, 1);
		    }
		}
		else if (t == outflow)
		{
		    if (is_remote(r, igrid))
                        continue;
                    //
		    // Do nothing if NODE-based, reflect if CELL-based.
                    //
		    if (type(r, idim) == IndexType::CELL)
		    {
			bb.shift(
			    idim,
			    (2 * domain.bigEnd(idim) + 1 + a
			     - b.bigEnd(idim) - b.smallEnd(idim)));
			const Box& rbox = r[igrid].box();
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
	    }
        }
    }
    tl.execute("mixed_boundary::fill_borders");
}

void
mixed_boundary::fill_sync_reg_borders (MultiFab&              r,
			               const level_interface& lev_interface,
          			       int                    w) const
{
//  This is the same as the fill_borders routine except that it
//    doesn't fill outside periodic boundaries

    w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;

    BL_ASSERT(w == 1 || w == 0);

    const Box& domain = lev_interface.domain();

    task_list tl;
    // we are looping over only the fine-fine faces
    for (int iface = 0;
	 iface < lev_interface.nboxes(level_interface::FACEDIM); iface++) 
    {
	if (lev_interface.geo(level_interface::FACEDIM, iface)
	    != level_interface::ALL)
	    break;

        Box c = lev_interface.box(level_interface::FACEDIM, iface);

	const int igrid =
	    lev_interface.grid(level_interface::FACEDIM, iface, 0);
	const int jgrid =
	    lev_interface.grid(level_interface::FACEDIM, iface, 1);

	if (igrid < 0 || jgrid < 0) 
	{
	    Box b = lev_interface.box(level_interface::FACEDIM, iface);
	    const int idim = lev_interface.fdim(iface);
	    const int a = (type(r,idim) == IndexType::NODE);
            //
	    // Need to do on x borders too in case y border is an interior face
            //
	    for (int i = 0; i < BL_SPACEDIM; i++) 
	    {
	        if (i != idim) 
	        {
		  if (domain.smallEnd(i) == b.smallEnd(i)) b.growLo(i, w);
		  if (domain.bigEnd(i)   == b.bigEnd(i)  ) b.growHi(i, w);

	        }
	    }
	    if (igrid < 0) 
	    {
		b.shift(idim, -a).growLo(idim, w-a).convert(type(r));
		const RegType t = ptr->bc[idim][0];
		Box bb = b;
		if (flowdim == -4 && (t == refWall || t == inflow)) 
		{
                    //
		    // Terrain sigma
                    //
		    if (is_remote(r, jgrid)) continue;
		    bb.shift(
			idim,
			2 * domain.smallEnd(idim) - 1 + a
			- b.bigEnd(idim) - b.smallEnd(idim));
		    const Box& rbox = r[jgrid].box();
		    for (int i = 0; i < r.nComp(); i++) 
		    {
			Real* rptr = r[jgrid].dataPtr(i);
			if ((i == idim + BL_SPACEDIM)
			    || (i >= BL_SPACEDIM && idim == BL_SPACEDIM - 1)) 
			{
			    FORT_FBNEG(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
			else 
			{
			    FORT_FBREF(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
		    }
		}
		else if (t == refWall) 
		{
		    if ( is_remote(r, jgrid) ) continue;
		    bb.shift(
			idim,
			(2 * domain.smallEnd(idim) - 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[jgrid].box();
		    if (idim == flowdim || flowdim == -3) 
		    {
			FORT_FBNEG(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else 
		    {
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
		else if (t == periodic) 
		{
#if 0
		    int isrc = lev_interface.exterior_ref(igrid);
		    bb.shift(idim, domain.length(idim));
		    //r[jgrid].copy(r[isrc], bb, 0, b, 0, r.nComp());
		    tl.add_task(new task_copy(tl, r, jgrid, b, r, isrc, bb));
#endif
		}
		else if (t == inflow) 
		{
		    if ( is_remote(r, jgrid) ) continue;
		    bb.shift(
			idim,
			(2 * domain.smallEnd(idim) - 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[jgrid].box();
		    if (flowdim == -2) 
		    {
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (flowdim == -1) 
		    {
                        //
			//BoxLib::Error("mixed_boundary::Don't know how to do inflow density");
			// Inflow density---just reflect interior for now
                        //
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (flowdim == -3) 
		    {
			FORT_FBNEG(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, 1);
		    }
		    else if (idim == flowdim) 
		    {
                        //
			// For this to work, fill_borders must be called exactly
			// once for each level of this variable.
                        //
			FORT_FBINFLO(r[jgrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(b),
				     r[jgrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(bb), &idim, 1);
		    }
		    else if (flowdim >= 0) 
		    {
                        //
			// transverse velocity components
			//r[jgrid].assign(0.0, b);
			// we now believe this looks like a refWall to transverse comps
                        //
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, 1);
		    }
		}
		else if (t == outflow) 
		{
                    //
		    // Do nothing if NODE-based, reflect if CELL-based.
                    //
		    if (is_remote(r, jgrid))
                        continue;
		    if (type(r,idim) == IndexType::CELL) 
		    {
			bb.shift(
			    idim,
			    (2 * domain.smallEnd(idim) - 1 + a
			     - b.bigEnd(idim) - b.smallEnd(idim)));
			const Box& rbox = r[jgrid].box();
			FORT_FBREF(r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[jgrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
	    }
	    else if (jgrid < 0) 
	    {
		b.shift(idim, a).growHi(idim, w-a).convert(type(r));
		const RegType t = ptr->bc[idim][1];
		Box bb = b;
		if (flowdim == -4 && (t == refWall || t == inflow)) 
		{
		    if (is_remote(r, igrid))
                        continue;
                    //
		    // Terrain sigma
                    //
		    bb.shift(
			idim,
			(2 * domain.bigEnd(idim) + 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[igrid].box();
		    for (int i = 0; i < r.nComp(); i++) 
		    {
			Real* rptr = r[igrid].dataPtr(i);
			if ((i == idim + BL_SPACEDIM)
			    || (i >= BL_SPACEDIM && idim == BL_SPACEDIM - 1)) 
			{
			    FORT_FBNEG(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
			else 
			{
			    FORT_FBREF(rptr, DIMLIST(rbox),
				       DIMLIST(b),
				       rptr, DIMLIST(rbox),
				       DIMLIST(bb), &idim, 1);
			}
		    }
		}
		else if (t == refWall) 
		{
		    if (is_remote(r, igrid))
                        continue;
		    bb.shift(
			idim,
			(2 * domain.bigEnd(idim) + 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[igrid].box();
		    if (idim == flowdim || flowdim == -3) 
		    {
			FORT_FBNEG(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else 
		    {
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
		else if (t == periodic) 
		{
#if 0
		    int isrc = lev_interface.exterior_ref(jgrid);
		    bb.shift(idim, -domain.length(idim));
		    //r[igrid].copy(r[isrc], bb, 0, b, 0, r.nComp());
		    tl.add_task(new task_copy(tl, r, igrid, b, r, isrc, bb));
#endif
		}
		else if (t == inflow) 
		{
		    if (is_remote(r, igrid)) continue;
		    bb.shift(
			idim,
			(2 * domain.bigEnd(idim) + 1 + a
			 - b.bigEnd(idim) - b.smallEnd(idim)));
		    const Box& rbox = r[igrid].box();
		    if (flowdim == -2) 
		    {
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(),
				   DIMLIST(rbox), DIMLIST(bb), &idim, 1);
		    }
		    else if (flowdim == -1) 
		    {
                        //
			//BoxLib::Error("mixed_boundary::Don't know how to do inflow density");
			// Inflow density---just reflect interior for now
                        //
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (flowdim == -3) 
		    {
			FORT_FBNEG(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		    else if (idim == flowdim) 
		    {
                        //
			// For this to work, fill_borders must be called exactly
			// once for each level of this variable.
                        //
			FORT_FBINFLO(r[igrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(b),
				     r[igrid].dataPtr(), DIMLIST(rbox),
				     DIMLIST(bb), &idim, 1);
		    }
		    else if (flowdim >= 0) 
		    {
                        //
			// transverse velocity components
			//r[igrid].assign(0.0, rbox);
			// we now believe this looks like a refWall to transverse comps
                        //
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, 1);
		    }
		}
		else if (t == outflow) 
		{
		    if (is_remote(r, igrid))
                        continue;
                    //
		    // Do nothing if NODE-based, reflect if CELL-based.
                    //
		    if (type(r,idim) == IndexType::CELL) 
		    {
			bb.shift(
			    idim,
			    (2 * domain.bigEnd(idim) + 1 + a
			     - b.bigEnd(idim) - b.smallEnd(idim)));
			const Box& rbox = r[igrid].box();
			FORT_FBREF(r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(b),
				   r[igrid].dataPtr(), DIMLIST(rbox),
				   DIMLIST(bb), &idim, r.nComp());
		    }
		}
	    }
        }
    }
    tl.execute("mixed_boundary::fill_sync_reg_borders");
}


void
mixed_boundary::check_against_boundary_ (BoxList&   bl,
					 List<int>& il,
					 const Box& b,
					 int        ib,
					 const Box& d,
					 int        dim1) const
{
    for (int i = dim1; i < BL_SPACEDIM; i++)
    {
	if (b.smallEnd(i) == d.smallEnd(i))
	{
	    if (ptr->bc[i][0] == refWall || ptr->bc[i][0] == inflow)
	    {
		Box bn = b;
		bn.shift(i, -b.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary_(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][0] == periodic)
	    {
		Box bn = b;
		bn.shift(i, d.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary_(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][0] == outflow)
	    {
		Box bn = b;
		bn.shift(i, -b.length(i));
		bl.append(bn);
		il.append(-2);
		check_against_boundary_(bl, il, bn, -1, d, i+1);
	    }
	    else
	    {
		BoxLib::Abort( "mixed_boundary::check_against_boundary():"
			       "Boundary type not supported" );
	    }
	}
	if (b.bigEnd(i) == d.bigEnd(i))
	{
	    if (ptr->bc[i][1] == refWall || ptr->bc[i][1] == inflow)
	    {
		Box bn = b;
		bn.shift(i, b.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary_(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][1] == periodic)
	    {
		Box bn = b;
		bn.shift(i, -d.length(i));
		bl.append(bn);
		il.append(ib);
		check_against_boundary_(bl, il, bn, ib, d, i+1);
	    }
	    else if (ptr->bc[i][1] == outflow)
	    {
		Box bn = b;
		bn.shift(i, b.length(i));
		bl.append(bn);
		il.append(-2);
		check_against_boundary_(bl, il, bn, -1, d, i+1);
	    }
	    else
	    {
		BoxLib::Abort( "mixed_boundary::check_against_boundary():"
			       "Boundary type not supported" );
	    }
	}
    }
}

void
mixed_boundary::duplicate (List<Box>& bl,
			   const Box& domain) const
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
			if (!bl.includes(btmp))
			{
			    bl.append(btmp);
			}
		    }
		    else if (bn().bigEnd(i) - 1 == domain.bigEnd(i))
		    {
			Box btmp = bn();
			btmp.shift(i, -domain.length(i));
			if (!bl.includes(btmp))
			{
			    bl.append(btmp);
			}
		    }
		}
	    }
	}
    }
}

bool
mixed_boundary::singular () const
{
    BL_ASSERT(flowdim == -2); // pressure boundary only

    for (int idim = 0; idim < BL_SPACEDIM; idim++)
    {
	if (ptr->bc[idim][0] == outflow || ptr->bc[idim][1] == outflow)
	{
            return false;
	}
    }

    return true;
}

amr_fluid_boundary::amr_fluid_boundary ()
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	v[i] = 0;
    }

    s  = 0;
    p  = 0;
    ts = 0;
}

amr_fluid_boundary::~amr_fluid_boundary () {}

const amr_boundary*
amr_fluid_boundary::velocity (int i) const
{
    return v[i];
}

const amr_boundary*
amr_fluid_boundary::scalar () const
{
    return s;
}

const amr_boundary*
amr_fluid_boundary::pressure () const
{
    return p;
}

const amr_boundary*
amr_fluid_boundary::terrain_sigma () const
{
    return ts;
}

inviscid_fluid_boundary::inviscid_fluid_boundary (const RegType Bc[BL_SPACEDIM][2])
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	bc[i][0] = Bc[i][0];
	bc[i][1] = Bc[i][1];
	if ((bc[i][0] == periodic || bc[i][1] == periodic)
	    && bc[i][1] != bc[i][0])
	{
	    BoxLib::Abort( "inviscid_fluid_boundary::inviscid_fluid_boundary():"
			   "periodic bc's don't match" );
	}
	v[i] = new mixed_boundary(this, i);
    }
    s  = new mixed_boundary(this, -1);
    p  = new mixed_boundary(this, -2);
    ts = new mixed_boundary(this, -4);
}

inviscid_fluid_boundary::~inviscid_fluid_boundary ()
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	delete v[i];
    }
    delete s;
    delete p;
    delete ts;
}

RegType
inviscid_fluid_boundary::getLoBC(int idim) const
{
    return bc[idim][0];
}

RegType
inviscid_fluid_boundary::getHiBC(int idim) const
{
    return bc[idim][1];
}
