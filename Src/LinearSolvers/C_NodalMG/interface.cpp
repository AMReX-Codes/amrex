
#include "interface.H"
#include "boundary.H"

static inline void ins(List<Box>& bl, const Box& b) 
{
    if (!bl.includes(b))
	bl.append(b);
}

level_interface::~level_interface()
{
    if (!ok())
	return;
    
    if (status & 1) 
    {
	// owns boxes
	delete [] pf;
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    delete [] bx[i];
	    if (i > 0) delete [] nodebx[i];
	}
    }
    
    if (status & 2) 
    {
	// owns flag arrays
	delete [] processorMap;
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

void level_interface::copy(const level_interface& src)
{
    if ( ok() )
	BoxLib::Error("level_interface::copy---this object already allocated");
    
    status = 0;

    dom = src.dom;
    im  = src.im;
    em  = src.em;
    grid_ref = src.grid_ref;
    fdm = src.fdm;
    processorMap = src.processorMap;
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

void level_interface::alloc_coarsened(const BoxArray& Im,
				      const amr_boundary_class* /*bdy*/,
				      const level_interface& src,
				      const IntVect& rat)
{
    if (ok())
	BoxLib::Error("level_interface::alloc_coarsened---this object already allocated");
    
    status = 1;
    
    dom = coarsen(src.dom, rat);
    im  = Im;
    
    grid_ref = src.grid_ref;
    fdm = src.fdm;
    processorMap = src.processorMap;
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
	    bx[i][igrid] = coarsen(src.bx[i][igrid], rat);
	}
    }
    
    nodebx[0] = bx[0];
    for (int i = 1; i < BL_SPACEDIM; i++) 
    {
	nodebx[i] = new Box[nbx[i]];
	for (int igrid = 0; igrid < nbx[i]; igrid++) 
	{
	    nodebx[i][igrid] = coarsen(src.nodebx[i][igrid], rat);
	}
    }
    
    pf = new Box[im.length()];
    
    for (int igrid = 0; igrid < im.length(); igrid++) 
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

void level_interface::alloc(const BoxArray& Im, const Box& Domain, const amr_boundary_class* bdy, const Array<int>& pd)
{
    if (ok())
	BoxLib::Error("level_interface::alloc---this object already allocated");
    
    status = 3;
    
    dom = Domain;
    im  = Im;
    assert( bdy != 0 );
    bdy->boundary_mesh(em, grid_ref, im, dom);

    processorMap = new int[im.length()];
    
    List<Box> bl;
    
    // Add edges in 2D or faces in 3D:
    
    for (int igrid = 0; igrid < im.length(); igrid++) 
    {
	processorMap[igrid] = pd[igrid];
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    IntVect t = IntVect::TheCellVector();
	    t.setVal(i, IndexType::NODE);
	    add(bl, bdryLo(im[igrid], i).convert(t), 0);
	    add(bl, bdryHi(im[igrid], i).convert(t), 0);
	}
    }
    
    assert(bdy != 0);
    bdy->duplicate(bl, dom);
    xfer(bl, FACEDIM);
    bl.clear();
    
#if (BL_SPACEDIM == 3)
    
    // Add edges in 3D:
    
    for (int iface = 0; iface < nbx[2]; iface++) 
    {
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    IntVect t = bx[2][iface].type();
	    if (t[i] == IndexType::NODE)
		continue;
	    else
		t.setVal(i, IndexType::NODE);
	    add(bl, bdryLo(bx[2][iface], i).convert(t), 0);
	    add(bl, bdryHi(bx[2][iface], i).convert(t), 0);
	}
    }
    
    assert(bdy != 0);
    bdy->duplicate(bl, dom);
    xfer(bl, 1);
    bl.clear();
#endif
    
    // Add corners:
    
    for (int iedge = 0; iedge < nbx[1]; iedge++) 
    {
	IntVect t = bx[1][iedge].type();
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    if (t[i] == IndexType::NODE)
		continue;
	    ins(bl, bdryLo(bx[1][iedge], i));
	    ins(bl, bdryHi(bx[1][iedge], i));
	}
    }
    
    bdy->duplicate(bl, dom);
    assert(bdy != 0);
    xfer(bl, 0);
    bl.clear();
    
    // initialize face direction array
    fdm = new int[nbx[FACEDIM]];
    for (int iface = 0; iface < nbx[FACEDIM]; iface++) 
    {
	IntVect t = bx[FACEDIM][iface].type();
	fdm[iface] = -1;
	for (int i = 0; i < BL_SPACEDIM; i++) 
	{
	    if (t[i] == IndexType::NODE)
	    {
		// one and only one face will be designated as a the direction
		assert(fdm[iface] == -1);
		fdm[iface] = i;
	    }
	}
	assert(fdm[iface] >= 0 && fdm[iface] < BL_SPACEDIM);
    }
    
    // initialize face grid array
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
		    for (int igrid = 0; igrid < im.length(); igrid++) 
		    {
			if (im[igrid].contains(b)) 
			{
			    fgr[iface][i] = igrid;
			    break;
			}
		    }
		}
		else 
		{
		    for (int igrid = 0; igrid < em.length(); igrid++) 
		    {
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
    // initialize edge grid array
    idim = 1;
    egr = new int[nbx[idim]][N_EDGE_GRIDS];
    for (int iedge = 0; iedge < nbx[idim]; iedge++) 
    {
	Box b = bx[idim][iedge];
	IntVect t = b.type();
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
		    for (int igrid = 0; igrid < im.length(); igrid++) 
		    {
			if (im[igrid].contains(b)) 
			{
			    egr[iedge][i] = igrid;
			    break;
			}
		    }
		}
		else 
		{
		    for (int igrid = 0; igrid < em.length(); igrid++) 
		    {
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
    
    // initialize corner grid array
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
		    for (int igrid = 0; igrid < im.length(); igrid++) 
		    {
			if (im[igrid].contains(b)) 
			{
			    cgr[icor][i] = igrid;
			    break;
			}
		    }
		}
		else 
		{
		    for (int igrid = 0; igrid < em.length(); igrid++) 
		    {
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
    
    pf = new Box[im.length()];
    
    for (int igrid = 0; igrid < im.length(); igrid++) 
    {
	//pf[igrid] = im.boxn(igrid);
	//pf[igrid].grow(-1);
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
			for (int iface = 0; iface < nbx[2]; iface++) 
			{
			    if (ax[2][iface] == igrid) 
			    {
				if (nodebx[2][iface].contains(nodebx[1][iedge]))
				    ax[1][iedge] = -1;
			    }
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
}

void level_interface::add(List<Box>& bl, Box b, int startgrid)
{
    const IntVect t = b.type();
    for (int igrid = startgrid; igrid < im.length() + em.length(); igrid++) 
    {
	Box ibox;
	if (igrid < im.length())
	    ibox = im[igrid];
	else
	    ibox = em[igrid-im.length()];
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
			add(bl, b, igrid + 1);
			add(bl, c, igrid);
			return;
		    }
		    if (ibox.bigEnd(i) < b.bigEnd(i)) 
		    {
			Box c = b.chop(i, ibox.bigEnd(i) + 1);
			add(bl, b, igrid);
			add(bl, c, igrid + 1);
			return;
		    }
		}
	    }
	    BoxLib::Error("level_interface::add()---shouldn't ever get here.");
	}
    }
    ins(bl, b);
}

void level_interface::xfer(const List<Box>& bl, int idim)
{
    nbx[idim] = bl.length();
    bx[idim]  = new Box[nbx[idim]];
    ge[idim]  = new unsigned int[nbx[idim]];
    flg[idim] = new bool[nbx[idim]];
    
    ListIterator<Box> bn(bl);
    for ( int i = 0; bn; bn++, i++) 
    {
	bx[idim][i] = bn();
	const Box btmp = grow(bn(), bn().type()).convert(IntVect::TheCellVector());
	IntVect tmp = btmp.smallEnd();
	if (dom.contains(btmp)) 
	{
	    ge[idim][i]  = im.contains(tmp);
#if (BL_SPACEDIM == 2)
	    tmp += IntVect(1,0);
	    ge[idim][i] |= im.contains(tmp) << 1;
	    tmp += IntVect(-1,1);
	    ge[idim][i] |= im.contains(tmp) << 2;
	    tmp += IntVect(1,0);
	    ge[idim][i] |= im.contains(tmp) << 3;
#else
	    tmp += IntVect(1,0,0);
	    ge[idim][i] |= im.contains(tmp) << 1;
	    tmp += IntVect(-1,1,0);
	    ge[idim][i] |= im.contains(tmp) << 2;
	    tmp += IntVect(1,0,0);
	    ge[idim][i] |= im.contains(tmp) << 3;
	    tmp += IntVect(-1,-1,1);
	    ge[idim][i] |= im.contains(tmp) << 4;
	    tmp += IntVect(1,0,0);
	    ge[idim][i] |= im.contains(tmp) << 5;
	    tmp += IntVect(-1,1,0);
	    ge[idim][i] |= im.contains(tmp) << 6;
	    tmp += IntVect(1,0,0);
	    ge[idim][i] |= im.contains(tmp) << 7;
#endif
	}
	else 
	{
	    bool is_in = dom.contains(tmp);
	    ge[idim][i]  = ( is_in && im.contains(tmp) || !is_in && em.contains(tmp));
#if (BL_SPACEDIM == 2)
	    tmp += IntVect(1,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 1;
	    tmp += IntVect(-1,1);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 2;
	    tmp += IntVect(1,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 3;
#else
	    tmp += IntVect(1,0,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 1;
	    tmp += IntVect(-1,1,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 2;
	    tmp += IntVect(1,0,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 3;
	    tmp += IntVect(-1,-1,1);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 4;
	    tmp += IntVect(1,0,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 5;
	    tmp += IntVect(-1,1,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 6;
	    tmp += IntVect(1,0,0);
	    is_in = dom.contains(tmp);
	    ge[idim][i] |= ( is_in && im.contains(tmp) || !is_in && em.contains(tmp)) << 7;
#endif
	}
    }
    
    // Sort fine-fine boxes to beginning of list
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
    int nff = j;
    
    // Sort interior fine-fine boxes to beginning of list
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
	int nin = j;
	
	// Sort interior faces according to orientation, x first
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
	
	// Sort exterior faces according to orientation, x first
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

Array<int> level_interface::geo_array(int idim, int i) const
{
    Array<int> ga(N_CORNER_GRIDS);
    unsigned int gtmp = geo(idim, i);
    for (int k = 0; k < N_CORNER_GRIDS; k++) 
    {
	ga[k] = (gtmp & 1);
	gtmp >>= 1;
    }
    return ga;
}
