#include <winstd.H>

#include <algorithm>

#include <Profiler.H>

#include "interface.H"
#include "boundary.H"

const level_interface default_level_interface;

static
void
ins (level_interface::BoxMSet& bmset, const Box& b)
{
    level_interface::BoxMSetConstIterPair er_it = bmset.equal_range(b);

    for (level_interface::BoxMSetConstIter it = er_it.first; it != er_it.second; ++it)
        if (*it == b)
            return;

    bmset.insert(b);
}

typedef std::pair<int,Box> IBpair;

struct IBpairComp
{
    bool operator () (const IBpair& lhs, const IBpair& rhs) const
    {
        return lhs.first < rhs.first;
    }
};

void
add (level_interface::BoxMSet& bmset,
     const BoxArray& bim, 
     Box      b,
     int      startgrid)
{
    BL_PROFILE("interface::add()");

    const IntVect t = b.type();

    Box tb(b.smallEnd(), b.bigEnd());
    
    std::vector< std::pair<int,Box> > prs = bim.intersections(BoxLib::grow(tb,1));

    std::sort(prs.begin(), prs.end(), IBpairComp());

    for ( int j = 0; j < prs.size(); ++j ) 
      {
	int igrid = prs[j].first;
	if ( igrid < startgrid ) continue;
        Box ibox = bim[igrid];
	ibox.convert(t);
	if (ibox.intersects(b) && !ibox.contains(b))
	{
	    for (int i = 0; i < BL_SPACEDIM; i++)
	    {
		if (t[i] == IndexType::CELL)
		{
		    if (ibox.smallEnd(i) > b.smallEnd(i))
		    {
			Box c = b.chop(i, ibox.smallEnd(i));
			add(bmset, bim, b, igrid + 1);
			add(bmset, bim, c, igrid);
			return;
		    }
		    if (ibox.bigEnd(i) < b.bigEnd(i))
		    {
			Box c = b.chop(i, ibox.bigEnd(i) + 1);
			add(bmset, bim, b, igrid);
			add(bmset, bim, c, igrid + 1);
			return;
		    }
		}
	    }
	    BoxLib::Abort("level_interface::add(): Can't happen.");
	}
    }
    ins(bmset, b);
}

level_interface::~level_interface()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::~level_interface()");

    if (!ok())
	return;

    if (status & 1)
    {
        //
	// Owns boxes.
        //
	delete [] pf;
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    delete [] bx[i];
	    if (i > 0)
                delete [] nodebx[i];
	}
    }

    if (status & 2)
    {
        //
	// Owns flag arrays.
        //
	delete [] grid_ref;
	delete [] fdm;
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    delete [] ge[i];
	    delete [] ax[i];
	    delete [] flg[i];
	}
	delete [] fgr;
#if (BL_SPACEDIM == 3)
	delete [] egr;
#endif
	delete [] cgr;
    }
}

void
level_interface::copy (const level_interface& src)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::copy()");

    if (ok()) BoxLib::Error( "level_interface::copy: this object already allocated" );

    m_fill_internal_borders_fn = src.m_fill_internal_borders_fn;
    m_fill_internal_borders_fc = src.m_fill_internal_borders_fc;

    status = 0;

    dom = src.dom;
    im  = src.im;
    em  = src.em;
    grid_ref = src.grid_ref;
    fdm = src.fdm;
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	nbx[i] = src.nbx[i];
	ge[i]  = src.ge[i];
	ax[i]  = src.ax[i];
	flg[i] = src.flg[i];
    }
    fgr = src.fgr;
#if BL_SPACEDIM==3
    egr = src.egr;
#endif
    cgr = src.cgr;

    pf = src.pf;
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	bx[i] = src.bx[i];
	nodebx[i] = src.nodebx[i];
    }
}

void
level_interface::alloc_coarsened (const BoxArray&           Im,
                                  const amr_boundary* /*bdy*/,
                                  const level_interface&    src,
                                  const IntVect&            rat)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::coarsened()");

    if (ok()) BoxLib::Error( "level_interface::alloc_coarsened: this object already allocated" );

    m_fill_internal_borders_fn.resize(src.nboxes(level_interface::FACEDIM));
    m_fill_internal_borders_fc.resize(src.nboxes(level_interface::FACEDIM));

    for (int i = 0; i < src.nboxes(level_interface::FACEDIM); i++)
    {
        m_fill_internal_borders_fn[i] = 1;
        m_fill_internal_borders_fc[i] = 1;
    }

    status = 1;

    dom = BoxLib::coarsen(src.dom, rat);
    im  = Im;

    grid_ref = src.grid_ref;
    fdm = src.fdm;
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	nbx[i] = src.nbx[i];
	ge[i]  = src.ge[i];
	ax[i]  = src.ax[i];
	flg[i] = src.flg[i];
    }
    fgr = src.fgr;
#if BL_SPACEDIM == 3
    egr = src.egr;
#endif
    cgr = src.cgr;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	bx[i] = new Box[nbx[i]];
	for (int igrid = 0; igrid < nbx[i]; igrid++)
	{
	    bx[i][igrid] = BoxLib::coarsen(src.bx[i][igrid], rat);
	}
    }

    nodebx[0] = bx[0];
    for (int i = 1; i < BL_SPACEDIM; i++)
    {
	nodebx[i] = new Box[nbx[i]];
	for (int igrid = 0; igrid < nbx[i]; igrid++)
	{
	    nodebx[i][igrid] = BoxLib::coarsen(src.nodebx[i][igrid], rat);
	}
    }

    pf = new Box[im.size()];

    for (int igrid = 0; igrid < im.size(); igrid++)
    {
	pf[igrid] = im[igrid];
	pf[igrid].convert(IntVect::TheNodeVector()).grow(-1);
    }

    int idim = FACEDIM;
    for (int iface = 0; iface < nbx[idim]; iface++)
    {
	if ( ge[idim][iface] == ALL && !flg[idim][iface] )
	{
	    int igrid;
	    if ((igrid = fgr[iface][0]) >= 0)
	    {
		if (!pf[igrid].intersects(nodebx[idim][iface]))
		    pf[igrid].growHi(fdm[iface], 1);
	    }
	    if ((igrid = fgr[iface][1]) >= 0)
	    {
		if (!pf[igrid].intersects(nodebx[idim][iface]))
		    pf[igrid].growLo(fdm[iface], 1);
	    }
	}
    }
}

void
level_interface::alloc (const BoxArray&     Im,
                        const Box&          Domain,
                        const amr_boundary* bdy)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::alloc()");

    if (ok()) BoxLib::Error( "level_interface::alloc: this object already allocated" );

    status = 3;

    dom = Domain;
    im  = Im;
    BL_ASSERT( bdy != 0 );
    bdy->boundary_mesh(em, grid_ref, im, dom);

    BoxArray bim;
    {
      std::vector<Box> tbim;
      for ( int i = 0; i < im.size(); ++i )
	{
	  tbim.push_back(im[i]);
	}
      for ( int i = 0; i < em.size(); ++i )
	{
	  tbim.push_back(em[i]);
	}
      bim = BoxArray(&tbim[0], tbim.size());
    }
    //
    // Add edges in 2D or faces in 3D:
    //
    BoxMSet bmset;
    for (int igrid = 0; igrid < im.size(); igrid++)
    {
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    IntVect t = IntVect::TheCellVector();
	    t[i] = IndexType::NODE;
	    Box lo = BoxLib::bdryLo(im[igrid],i).convert(t);
	    Box hi = BoxLib::bdryHi(im[igrid],i).convert(t);
	    add(bmset,bim,lo,0);
	    add(bmset,bim,hi,0);
	}
    }

    bdy->duplicate(bmset, dom);
    xfer(bmset, FACEDIM);
    bmset.clear();

#if (BL_SPACEDIM == 3)
    //
    // Add edges in 3D:
    //
    for (int iface = 0; iface < nbx[2]; iface++)
    {
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    IntVect t = bx[2][iface].type();
	    if (t[i] == IndexType::NODE)
		continue;
	    else
		t[i] = IndexType::NODE;
	    Box lo = BoxLib::bdryLo(bx[2][iface], i).convert(t);
	    Box hi = BoxLib::bdryHi(bx[2][iface], i).convert(t);
	    add(bmset, bim,lo, 0);
	    add(bmset, bim,hi, 0);
	}
    }

    bdy->duplicate(bmset, dom);
    xfer(bmset, 1);
    bmset.clear();
#endif
    //
    // Add corners:
    //
    for (int iedge = 0; iedge < nbx[1]; iedge++)
    {
	const IntVect t = bx[1][iedge].type();
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    if (t[i] == IndexType::NODE)
		continue;
	    ins(bmset, BoxLib::bdryLo(bx[1][iedge], i));
	    ins(bmset, BoxLib::bdryHi(bx[1][iedge], i));
	}
    }

    bdy->duplicate(bmset, dom);
    xfer(bmset, 0);
    bmset.clear();
    //
    // Initialize face direction array.
    //
    fdm = new int[nbx[FACEDIM]];
    for (int iface = 0; iface < nbx[FACEDIM]; iface++)
    {
	const IntVect t = bx[FACEDIM][iface].type();
	fdm[iface] = -1;
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    if (t[i] == IndexType::NODE)
	    {
                //
		// One and only one face will be designated as a the direction.
                //
		BL_ASSERT(fdm[iface] == -1);
		fdm[iface] = i;
	    }
	}
	BL_ASSERT(fdm[iface] >= 0 && fdm[iface] < BL_SPACEDIM);
    }
    //
    // Initialize face grid array.
    //
    int idim = FACEDIM;
    fgr = new int[nbx[idim]][N_FACE_GRIDS];
    for (int iface = 0; iface < nbx[idim]; iface++)
    {
	Box b = bx[idim][iface];
	int id = fdm[iface];
	b.growLo(id, 1).convert(IntVect::TheCellVector());
	int imask = 1;
	flg[idim][iface] = false;
	for (int i = 0; i < N_FACE_GRIDS; i++)
	{
	    fgr[iface][i] = -1;
	    if ( ge[idim][iface] & imask )
	    {
		if (dom.contains(b))
		{
                    std::vector< std::pair<int,Box> > isects = im.intersections(b);

                    for (int j = 0; j < isects.size(); j++)
                    {
                        const int igrid = isects[j].first;

                        if (im[igrid].contains(b))
                        {
			    fgr[iface][i] = igrid;
			    break;
                        }
                    }
		}
		else
		{
                    std::vector< std::pair<int,Box> > isects = em.intersections(b);

                    for (int j = 0; j < isects.size(); j++)
                    {
                        const int igrid = isects[j].first;

			if (em[igrid].contains(b))
			{
			    fgr[iface][i] = -2 - igrid;
			    if (grid_ref[igrid] == -2)
				flg[idim][iface] = true;
			    break;
			}
                    }
		}
	    }
	    b.shift(id, 1);
	    imask <<= (1 << id);
	}
    }

#if (BL_SPACEDIM == 2)
    // egr = fgr;
#else
    //
    // Initialize edge grid array.
    //
    idim = 1;
    egr = new int[nbx[idim]][N_EDGE_GRIDS];
    for (int iedge = 0; iedge < nbx[idim]; iedge++)
    {
	Box b = bx[idim][iedge];
	const IntVect t = b.type();
	int id = 0;
	if (t[id] == IndexType::CELL)
	    id++;
	int jd = id + 1;
	if (t[jd] == IndexType::CELL)
	    jd++;
	b.growLo(id, 1).growLo(jd, 1).convert(IntVect::TheCellVector());
	int imask = 1;
	flg[idim][iedge] = false;
	for (int i = 0; i < N_EDGE_GRIDS; i++)
	{
	    egr[iedge][i] = -1;
	    if ( ge[idim][iedge] & imask )
	    {
		if (dom.contains(b))
		{
                    std::vector< std::pair<int,Box> > isects = im.intersections(b);

                    for (int j = 0; j < isects.size(); j++)
                    {
                        const int igrid = isects[j].first;

			if (im[igrid].contains(b))
			{
			    egr[iedge][i] = igrid;
			    break;
			}
		    }
		}
		else
		{
                    std::vector< std::pair<int,Box> > isects = em.intersections(b);

                    for (int j = 0; j < isects.size(); j++)
                    {
                        const int igrid = isects[j].first;

			if (em[igrid].contains(b))
			{
			    egr[iedge][i] = -2 - igrid;
			    if (grid_ref[igrid] == -2)
				flg[idim][iedge] = true;
			    break;
			}
		    }
		}
	    }
	    if ((i & 1) == 1)
	    {
		b.shift(id, -1).shift(jd, 1);
		imask <<= ((1 << jd) - (1 << id));
	    }
	    else
	    {
		b.shift(id, 1);
		imask <<= (1 << id);
	    }
	}
    }
#endif
    //
    // Initialize corner grid array.
    //
    idim = 0;
    cgr = new int[nbx[idim]][N_CORNER_GRIDS];
    for (int icor = 0; icor < nbx[idim]; icor++)
    {
	Box b = bx[idim][icor];
#if (BL_SPACEDIM == 3)
	b.growLo(2, 1);
#endif
	b.growLo(0, 1).growLo(1, 1).convert(IntVect::TheCellVector());
	int imask = 1;
	flg[idim][icor] = false;
	for (int i = 0; i < N_CORNER_GRIDS; i++)
	{
	    cgr[icor][i] = -1;
	    if ( ge[idim][icor] & imask )
	    {
		if (dom.contains(b))
		{
                    std::vector< std::pair<int,Box> > isects = im.intersections(b);

                    for (int j = 0; j < isects.size(); j++)
                    {
                        const int igrid = isects[j].first;

			if (im[igrid].contains(b))
			{
			    cgr[icor][i] = igrid;
			    break;
			}
		    }
		}
		else
		{
                    std::vector< std::pair<int,Box> > isects = em.intersections(b);

                    for (int j = 0; j < isects.size(); j++)
                    {
                        const int igrid = isects[j].first;

			if (em[igrid].contains(b))
			{
			    cgr[icor][i] = -2 - igrid;
			    if (grid_ref[igrid] == -2)
				flg[idim][icor] = true;
			    break;
			}
		    }
		}
	    }
#if (BL_SPACEDIM == 3)
	    if ((i & 3) == 3)
		b.shift(0, -1).shift(1, -1).shift(2, 1);
	    else
#endif
		if ((i & 1) == 1)
		    b.shift(0, -1).shift(1, 1);
		else
		    b.shift(0, 1);
		imask <<= 1;
	}
    }

    nodebx[0] = bx[0];
    for (int i = 1; i < BL_SPACEDIM; i++)
    {
	nodebx[i] = new Box[nbx[i]];
	for (int iface = 0; iface < nbx[i]; iface++)
	{
	    nodebx[i][iface] = bx[i][iface];
	    nodebx[i][iface].convert(IntVect::TheNodeVector());
	}
    }

    pf = new Box[im.size()];

    for (int igrid = 0; igrid < im.size(); igrid++)
    {
	pf[igrid] = im[igrid];
	pf[igrid].convert(IntVect::TheNodeVector()).grow(-1);
    }

    idim = FACEDIM;
    for (int iface = 0; iface < nbx[idim]; iface++)
    {
	if ( ge[idim][iface] == ALL && !flg[idim][iface] )
	{
	    int igrid;
	    if ((igrid = fgr[iface][0]) >= 0)
	    {
		if (!pf[igrid].intersects(nodebx[idim][iface]))
		    pf[igrid].growHi(fdm[iface], 1);
	    }
	    if ((igrid = fgr[iface][1]) >= 0)
	    {
		if (!pf[igrid].intersects(nodebx[idim][iface]))
		    pf[igrid].growLo(fdm[iface], 1);
	    }
	}
    }

    ax[idim] = new int[nbx[idim]];
    for (int iface = 0; iface < nbx[idim]; iface++)
    {
	ax[idim][iface] = -1;
	if ( ge[idim][iface] != ALL )
	{
	    for (int i = 0; i < 2; i++)
	    {
		int igrid;
		if ((igrid = fgr[iface][i]) >= 0)
		{
		    if (pf[igrid].intersects(nodebx[idim][iface]))
			ax[idim][iface] = igrid;
		}
	    }
	}
    }

#if (BL_SPACEDIM == 3)
    ax[1] = new int[nbx[1]];

    for (int iedge = 0; iedge < nbx[1]; iedge++)
    {
	ax[1][iedge] = -1;
	if ( ge[1][iedge] != ALL )
	{
	    for (int i = 0; i < N_EDGE_GRIDS && ax[1][iedge] == -1; i++)
	    {
		int igrid;
		if ((igrid = egr[iedge][i]) >= 0)
		{
		    if (pf[igrid].intersects(nodebx[1][iedge]))
		    {
			ax[1][iedge] = igrid;

                        IntVect iv(1,1,1);
                        IndexType ityp(IndexType(iv));

                        int first = -1;
			for (int iface = 0; iface < nbx[2] && first==-1; iface++)
                        {
                            if (ax[2][iface] == igrid) first = iface;
                        }
                        if (first != -1)
                        {
                            BoxList face_list(nodebx[2][first]);
                            for (int iface = 0; iface < nbx[2]; iface++)
                                if (ax[2][iface] == igrid && iface != first) 
                                    face_list.push_back(nodebx[2][iface]);
                            BoxArray face_array(face_list);

                            if (face_array.contains(nodebx[1][iedge]))
                                ax[1][iedge] = -1;
                        }
		    }
		}
	    }
	}
    }

#endif

    ax[0] = new int[nbx[0]];
    for (int icor = 0; icor < nbx[0]; icor++)
    {
	ax[0][icor] = -1;
	if ( ge[0][icor] != ALL )
	{
	    for (int i = 0; i < N_CORNER_GRIDS && ax[0][icor] == -1; i++)
	    {
		int igrid;
		if ((igrid = cgr[icor][i]) >= 0)
		{
		    if (pf[igrid].intersects(nodebx[0][icor]))
		    {
			ax[0][icor] = igrid;
#if (BL_SPACEDIM == 3)
			for (int iface = 0; iface < nbx[2]; iface++)
			{
			    if (ax[2][iface] == igrid)
			    {
				if (nodebx[2][iface].contains(nodebx[0][icor]))
				    ax[0][icor] = -1;
			    }
			}
#endif
			for (int iedge = 0; iedge < nbx[1]; iedge++)
			{
			    if (ax[1][iedge] == igrid)
			    {
				if (nodebx[1][iedge].contains(nodebx[0][icor]))
				    ax[0][icor] = -1;
			    }
			}
		    }
		}
	    }
	}
    }

    m_fill_internal_borders_fn.resize(nboxes(level_interface::FACEDIM));
    m_fill_internal_borders_fc.resize(nboxes(level_interface::FACEDIM));

    for (int i = 0; i < nboxes(level_interface::FACEDIM); i++)
    {
        m_fill_internal_borders_fn[i] = 1;
        m_fill_internal_borders_fc[i] = 1;
    }
}


void
level_interface::xfer (const BoxMSet& bmset,
                       int            idim)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::xfer()");
    nbx[idim] = bmset.size();
    bx[idim]  = new Box[nbx[idim]];
    ge[idim]  = new unsigned int[nbx[idim]];
    flg[idim] = new bool[nbx[idim]];

    BoxMSetConstIter bn = bmset.begin();

    for (int i = 0; bn != bmset.end(); ++bn, ++i)
    {
	bx[idim][i] = *bn;
	const Box btmp = BoxLib::grow(*bn,bn->type()).convert(IntVect::TheCellVector());
	IntVect tmp = btmp.smallEnd();
	if (dom.contains(btmp))
	{
	    ge[idim][i]  = im.contains(tmp);
#if (BL_SPACEDIM == 2)
	    tmp += IntVect(1, 0);
	    ge[idim][i] |= im.contains(tmp) << 1;
	    tmp += IntVect(-1, 1);
	    ge[idim][i] |= im.contains(tmp) << 2;
	    tmp += IntVect(1, 0);
	    ge[idim][i] |= im.contains(tmp) << 3;
#else
	    tmp += IntVect(1, 0, 0);
	    ge[idim][i] |= im.contains(tmp) << 1;
	    tmp += IntVect(-1, 1, 0);
	    ge[idim][i] |= im.contains(tmp) << 2;
	    tmp += IntVect(1, 0, 0);
	    ge[idim][i] |= im.contains(tmp) << 3;
	    tmp += IntVect(-1, -1, 1);
	    ge[idim][i] |= im.contains(tmp) << 4;
	    tmp += IntVect(1, 0, 0);
	    ge[idim][i] |= im.contains(tmp) << 5;
	    tmp += IntVect(-1, 1, 0);
	    ge[idim][i] |= im.contains(tmp) << 6;
	    tmp += IntVect(1, 0, 0);
	    ge[idim][i] |= im.contains(tmp) << 7;
#endif
	}
	else
	{
	    bool is_in = dom.contains(tmp);
	    ge[idim][i]  = ( is_in && im.contains(tmp) || !is_in && em.contains(tmp));
#if (BL_SPACEDIM == 2)
	    tmp += IntVect(1, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 1;
	    tmp += IntVect(-1, 1);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 2;
	    tmp += IntVect(1, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 3;
#else
	    tmp += IntVect(1, 0, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 1;
	    tmp += IntVect(-1, 1, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 2;
	    tmp += IntVect(1, 0, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 3;
	    tmp += IntVect(-1, -1, 1);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 4;
	    tmp += IntVect(1, 0, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 5;
	    tmp += IntVect(-1, 1, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 6;
	    tmp += IntVect(1, 0, 0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 7;
#endif
	}
    }
    //
    // Sort fine-fine boxes to beginning of list.
    //
    int j = nbx[idim];
    for (int i = 0; i < j; i++)
    {
	if (ge[idim][i] != ALL)
	{
	    while (--j > i)
	    {
		if (ge[idim][j] == ALL)
		{
		    Box btmp = bx[idim][j];
		    bx[idim][j] = bx[idim][i];
		    bx[idim][i] = btmp;
		    ge[idim][j] = ge[idim][i];
		    ge[idim][i] = ALL;
		    break;
		}
	    }
	}
    }

    Box idomain = dom;
    idomain.convert(IntVect::TheNodeVector()).grow(-1);

    j = -1;
    while (++j < nbx[idim] && ge[idim][j] == ALL)
	/* nothing*/;
    const int nff = j;
    //
    // Sort interior fine-fine boxes to beginning of list.
    //
    if (idim == 0)
    {
	for (int i = 0; i < j; i++)
	{
	    if (!bx[idim][i].intersects(idomain))
	    {
		while (--j > i)
		{
		    Box btmp = bx[idim][j];
		    if (btmp.intersects(idomain))
		    {
			bx[idim][j] = bx[idim][i];
			bx[idim][i] = btmp;
			break;
		    }
		}
	    }
	}
    }
    else
    {
	for (int i = 0; i < j; i++)
	{
	    Box btmp = bx[idim][i];
	    btmp.convert(IntVect::TheNodeVector());
	    if (!btmp.intersects(idomain))
	    {
		while (--j > i)
		{
		    Box btmp = bx[idim][j];
		    btmp.convert(IntVect::TheNodeVector());
		    if (btmp.intersects(idomain))
		    {
			Box btmp = bx[idim][j];
			bx[idim][j] = bx[idim][i];
			bx[idim][i] = btmp;
			break;
		    }
		}
	    }
	}
    }

    if (idim == FACEDIM)
    {
	j = -1;
	while (++j < nff)
	{
	    Box btmp = bx[idim][j];
	    btmp.convert(IntVect::TheNodeVector());
	    if (!btmp.intersects(idomain))
		break;
	}
	const int nin = j;
	//
	// Sort interior faces according to orientation, x first.
        //
	for (int i = 0; i < j; i++)
	{
	    if (bx[idim][i].type(0) == IndexType::CELL)
	    {
		while (--j > i)
		{
		    if (bx[idim][j].type(0) == IndexType::NODE)
		    {
			Box btmp = bx[idim][j];
			bx[idim][j] = bx[idim][i];
			bx[idim][i] = btmp;
			break;
		    }
		}
	    }
	}
#if (BL_SPACEDIM == 3)
	j = nin;
	for (int i = 0; i < j; i++)
	{
	    if (bx[idim][i].type(0) == IndexType::CELL)
	    {
		if (bx[idim][i].type(1) == IndexType::CELL)
		{
		    while (--j > i)
		    {
			if (bx[idim][j].type(1) == IndexType::NODE)
			{
			    Box btmp = bx[idim][j];
			    bx[idim][j] = bx[idim][i];
			    bx[idim][i] = btmp;
			    break;
			}
		    }
		}
	    }
	}
#endif
	//
	// Sort exterior faces according to orientation, x first.
        //
	j = nff;
	for (int i = nin; i < j; i++)
	{
	    if (bx[idim][i].type(0) == IndexType::CELL)
	    {
		while (--j > i)
		{
		    if (bx[idim][j].type(0) == IndexType::NODE)
		    {
			Box btmp = bx[idim][j];
			bx[idim][j] = bx[idim][i];
			bx[idim][i] = btmp;
			break;
		    }
		}
	    }
	}
#if (BL_SPACEDIM == 3)
	j = nff;
	for (int i = nin; i < j; i++)
	{
	    if (bx[idim][i].type(0) == IndexType::CELL)
	    {
		if (bx[idim][i].type(1) == IndexType::CELL)
		{
		    while (--j > i)
		    {
			if (bx[idim][j].type(1) == IndexType::NODE)
			{
			    Box btmp = bx[idim][j];
			    bx[idim][j] = bx[idim][i];
			    bx[idim][i] = btmp;
			    break;
			}
		    }
		}
	    }
	}
#endif
    }
}

Array<int>
level_interface::geo_array (int idim,
                            int i) const
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::geo_array()");
    Array<int> ga(N_CORNER_GRIDS);
    unsigned int gtmp = geo(idim, i);
    for (int k = 0; k < N_CORNER_GRIDS; k++)
    {
	ga[k] = (gtmp & 1);
	gtmp >>= 1;
    }
    return ga;
}
