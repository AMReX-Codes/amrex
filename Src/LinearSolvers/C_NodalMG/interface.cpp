
#include "interface.H"
#include "boundary.H"

const level_interface null_level_interface;

level_interface::~level_interface()
{
  if (null())
    return;

  if (status & 1) {
    // owns boxes
    delete [] pf;
    for (int i = 0; i < BL_SPACEDIM; i++) {
      delete [] bx[i];
      if (i > 0) delete [] nodebx[i];
    }
  }

  if (status & 2) {
    // owns flag arrays
    delete [] grid_ref;
    delete [] fdm;
    for (int i = 0; i < BL_SPACEDIM; i++) {
      delete [] ge[i]; delete [] ax[i]; delete [] flg[i];
    }
#if (BL_SPACEDIM == 3)
    delete [] fgr;
#endif
    delete [] egr;
    delete [] cgr;
  }
}

void level_interface::copy(const level_interface& src)
{
  if (ok())
    BoxLib::Error("level_interface::copy---this object already allocated");

  status = 0;

  int i;

  dom = src.dom;
  im = src.im;
  em = src.em;
  grid_ref = src.grid_ref;
  fdm = src.fdm;
  for (i = 0; i < BL_SPACEDIM; i++) {
    nbx[i] = src.nbx[i];
    ge[i]  = src.ge[i];
    ax[i]  = src.ax[i];
    flg[i] = src.flg[i];
  }
  fgr = src.fgr;
  egr = src.egr;
  cgr = src.cgr;

  pf = src.pf;
  for (i = 0; i < BL_SPACEDIM; i++) {
    bx[i] = src.bx[i];
    nodebx[i] = src.nodebx[i];
  }
}

void level_interface::alloc_coarsened(const BoxArray& Im,
				      const amr_boundary_class& /*bdy*/,
				      const level_interface& src,
				      const IntVect& rat)
{
  if (ok())
    BoxLib::Error("level_interface::alloc_coarsened---this object already allocated");

  status = 1;

  int igrid, iface, idim, i;

  dom = coarsen(src.dom, rat);
  im = Im;

  grid_ref = src.grid_ref;
  fdm = src.fdm;
  for (i = 0; i < BL_SPACEDIM; i++) {
    nbx[i] = src.nbx[i];
    ge[i]  = src.ge[i];
    ax[i]  = src.ax[i];
    flg[i] = src.flg[i];
  }
  fgr = src.fgr;
  egr = src.egr;
  cgr = src.cgr;

  for (i = 0; i < BL_SPACEDIM; i++) {
    bx[i] = new Box[nbx[i]];
    for (igrid = 0; igrid < nbx[i]; igrid++) {
      bx[i][igrid] = coarsen(src.bx[i][igrid], rat);
    }
  }

  nodebx[0] = bx[0];
  for (i = 1; i < BL_SPACEDIM; i++) {
    nodebx[i] = new Box[nbx[i]];
    for (igrid = 0; igrid < nbx[i]; igrid++) {
      nodebx[i][igrid] = coarsen(src.nodebx[i][igrid], rat);
    }
  }

  pf = new Box[im.length()];

  for (igrid = 0; igrid < im.length(); igrid++) {
    pf[igrid] = im[igrid];
    pf[igrid].convert(IntVect::TheNodeVector()).grow(-1);
  }

  idim = FACEDIM;
  for (iface = 0; iface < nbx[idim]; iface++) {
    if (ge[idim][iface] == ALL && flg[idim][iface] == 0) {
      if ((igrid = fgr[iface][0]) >= 0) {
	if (!pf[igrid].intersects(nodebx[idim][iface]))
	  pf[igrid].growHi(fdm[iface], 1);
      }
      if ((igrid = fgr[iface][1]) >= 0) {
	if (!pf[igrid].intersects(nodebx[idim][iface]))
	  pf[igrid].growLo(fdm[iface], 1);
      }
    }
  }
}

void level_interface::alloc(const BoxArray& Im, const Box& Domain,
			    const amr_boundary_class& bdy)
{
  if (ok())
    BoxLib::Error("level_interface::alloc---this object already allocated");

  status = 3;

  int igrid, iface, iedge, icor, idim, i;

  dom = Domain;
  im = Im;
  bdy.boundary_mesh(em, grid_ref, im, dom);

  List<Box> bl;

  // Add edges in 2D or faces in 3D:

  for (igrid = 0; igrid < im.length(); igrid++) {
    for (i = 0; i < BL_SPACEDIM; i++) {
      IntVect t = IntVect::TheCellVector();
      t.setVal(i, IndexType::NODE);
      add(bl, bdryLo(im[igrid], i).convert(t));
      add(bl, bdryHi(im[igrid], i).convert(t));
    }
  }

  bdy.duplicate(bl, dom);

  xfer(bl, FACEDIM);

#if (BL_SPACEDIM == 3)

  // Add edges in 3D:

  for (iface = 0; iface < nbx[2]; iface++) {
    for (i = 0; i < BL_SPACEDIM; i++) {
      IntVect t = bx[2][iface].type();
      if (t[i] == IndexType::NODE)
	continue;
      else
	t.setVal(i, IndexType::NODE);
      add(bl, bdryLo(bx[2][iface], i).convert(t));
      add(bl, bdryHi(bx[2][iface], i).convert(t));
    }
  }

  bdy.duplicate(bl, dom);

  xfer(bl, 1);
#endif

  // Add corners:

  for (iedge = 0; iedge < nbx[1]; iedge++) {
    IntVect t = bx[1][iedge].type();
    for (i = 0; i < BL_SPACEDIM; i++) {
      if (t[i] == IndexType::NODE)
	continue;
      ins(bl, bdryLo(bx[1][iedge], i));
      ins(bl, bdryHi(bx[1][iedge], i));
    }
  }

  bdy.duplicate(bl, dom);

  xfer(bl, 0);

  // initialize face direction array
  fdm = new int[nbx[FACEDIM]];
  for (iface = 0; iface < nbx[FACEDIM]; iface++) {
    IntVect t = bx[FACEDIM][iface].type();
    for (i = 0; i < BL_SPACEDIM; i++) {
      if (t[i] == IndexType::NODE)
	fdm[iface] = i;
    }
  }

  // initialize face grid array
  idim = FACEDIM;
  fgr = new int[nbx[idim]][N_FACE_GRIDS];
  for (iface = 0; iface < nbx[idim]; iface++) {
    Box b = bx[idim][iface];
    int id = fdm[iface];
    b.growLo(id, 1).convert(IntVect::TheCellVector());
    int imask = 1;
    flg[idim][iface] = 0;
    for (i = 0; i < N_FACE_GRIDS; i++) {
      fgr[iface][i] = -1;
      if (ge[idim][iface] & imask) {
	if (dom.contains(b)) {
	  for (igrid = 0; igrid < im.length(); igrid++) {
	    if (im[igrid].contains(b)) {
	      fgr[iface][i] = igrid;
	      break;
	    }
	  }
	}
	else {
	  for (igrid = 0; igrid < em.length(); igrid++) {
	    if (em[igrid].contains(b)) {
	      fgr[iface][i] = -2 - igrid;
	      if (grid_ref[igrid] == -2)
		flg[idim][iface] = 1;
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
  egr = fgr;
#else
  // initialize edge grid array
  idim = 1;
  egr = new int[nbx[idim]][N_EDGE_GRIDS];
  for (iedge = 0; iedge < nbx[idim]; iedge++) {
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
    flg[idim][iedge] = 0;
    for (i = 0; i < N_EDGE_GRIDS; i++) {
      egr[iedge][i] = -1;
      if (ge[idim][iedge] & imask) {
	if (dom.contains(b)) {
	  for (igrid = 0; igrid < im.length(); igrid++) {
	    if (im[igrid].contains(b)) {
	      egr[iedge][i] = igrid;
	      break;
	    }
	  }
	}
	else {
	  for (igrid = 0; igrid < em.length(); igrid++) {
	    if (em[igrid].contains(b)) {
	      egr[iedge][i] = -2 - igrid;
	      if (grid_ref[igrid] == -2)
		flg[idim][iedge] = 1;
	      break;
	    }
	  }
	}
      }
      if ((i & 1) == 1) {
	b.shift(id, -1).shift(jd, 1);
	imask <<= ((1 << jd) - (1 << id));
      }
      else {
	b.shift(id, 1);
	imask <<= (1 << id);
      }
    }
  }
#endif

  // initialize corner grid array
  idim = 0;
  cgr = new int[nbx[idim]][N_CORNER_GRIDS];
  for (icor = 0; icor < nbx[idim]; icor++) {
    Box b = bx[idim][icor];
#if (BL_SPACEDIM == 3)
    b.growLo(2, 1);
#endif
    b.growLo(0, 1).growLo(1, 1).convert(IntVect::TheCellVector());
    int imask = 1;
    flg[idim][icor] = 0;
    for (i = 0; i < N_CORNER_GRIDS; i++) {
      cgr[icor][i] = -1;
      if (ge[idim][icor] & imask) {
	if (dom.contains(b)) {
	  for (igrid = 0; igrid < im.length(); igrid++) {
	    if (im[igrid].contains(b)) {
	      cgr[icor][i] = igrid;
	      break;
	    }
	  }
	}
	else {
	  for (igrid = 0; igrid < em.length(); igrid++) {
	    if (em[igrid].contains(b)) {
	      cgr[icor][i] = -2 - igrid;
	      if (grid_ref[igrid] == -2)
		flg[idim][icor] = 1;
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
  for (i = 1; i < BL_SPACEDIM; i++) {
    nodebx[i] = new Box[nbx[i]];
    for (iface = 0; iface < nbx[i]; iface++) {
      nodebx[i][iface] = bx[i][iface];
      nodebx[i][iface].convert(IntVect::TheNodeVector());
    }
  }

  pf = new Box[im.length()];

  for (igrid = 0; igrid < im.length(); igrid++) {
    //pf[igrid] = im.boxn(igrid);
    //pf[igrid].grow(-1);
    pf[igrid] = im[igrid];
    pf[igrid].convert(IntVect::TheNodeVector()).grow(-1);
  }

  idim = FACEDIM;
  for (iface = 0; iface < nbx[idim]; iface++) {
    if (ge[idim][iface] == ALL && flg[idim][iface] == 0) {
      if ((igrid = fgr[iface][0]) >= 0) {
	if (!pf[igrid].intersects(nodebx[idim][iface]))
	  pf[igrid].growHi(fdm[iface], 1);
      }
      if ((igrid = fgr[iface][1]) >= 0) {
	if (!pf[igrid].intersects(nodebx[idim][iface]))
	  pf[igrid].growLo(fdm[iface], 1);
      }
    }
  }

  ax[idim] = new int[nbx[idim]];
  for (iface = 0; iface < nbx[idim]; iface++) {
    ax[idim][iface] = -1;
    if (ge[idim][iface] != ALL) {
      for (i = 0; i < 2; i++) {
	if ((igrid = fgr[iface][i]) >= 0) {
	  if (pf[igrid].intersects(nodebx[idim][iface]))
	    ax[idim][iface] = igrid;
	}
      }
    }
  }

#if (BL_SPACEDIM == 3)
  ax[1] = new int[nbx[1]];
  for (iedge = 0; iedge < nbx[1]; iedge++) {
    ax[1][iedge] = -1;
    if (ge[1][iedge] != ALL) {
      for (i = 0; i < N_EDGE_GRIDS && ax[1][iedge] == -1; i++) {
	if ((igrid = egr[iedge][i]) >= 0) {
	  if (pf[igrid].intersects(nodebx[1][iedge])) {
	    ax[1][iedge] = igrid;
	    for (iface = 0; iface < nbx[2]; iface++) {
	      if (ax[2][iface] == igrid) {
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
  for (icor = 0; icor < nbx[0]; icor++) {
    ax[0][icor] = -1;
    if (ge[0][icor] != ALL) {
      for (i = 0; i < N_CORNER_GRIDS && ax[0][icor] == -1; i++) {
	if ((igrid = cgr[icor][i]) >= 0) {
	  if (pf[igrid].intersects(nodebx[0][icor])) {
	    ax[0][icor] = igrid;
#if (BL_SPACEDIM == 3)
	    for (iface = 0; iface < nbx[2]; iface++) {
	      if (ax[2][iface] == igrid) {
		if (nodebx[2][iface].contains(nodebx[0][icor]))
		  ax[0][icor] = -1;
	      }
	    }
#endif
	    for (iedge = 0; iedge < nbx[1]; iedge++) {
	      if (ax[1][iedge] == igrid) {
		if (nodebx[1][iedge].contains(nodebx[0][icor]))
		  ax[0][icor] = -1;
	      }
	    }
	  }
	}
      }
    }
  }
/*
  cout << im.length() << " interior grids, "
       << em.length() << " exterior grids, ";
#if (BL_SPACEDIM == 3)
  cout << nbx[2] << " faces, ";
#endif
  cout << nbx[1] << " edges, " << nbx[0] << " corners" << endl;

  cout << "Exterior grids are:" << endl;
  for (igrid = 0; igrid < em.length(); igrid++) {
    cout << em[igrid] << " " << grid_ref[igrid] << endl;
  }
  cout << "Part-fine grids are:" << endl;
  for (igrid = 0; igrid < im.length(); igrid++) {
    cout << pf[igrid] << endl;
  }
  cout << "box    geo   flag   aux   grids" << endl;
  cout << "Faces are:" << endl;
  for (iface = 0; iface < nbx[FACEDIM]; iface++) {
    cout << bx[FACEDIM][iface] << " " << ge[FACEDIM][iface] << " "
	 << flg[FACEDIM][iface] << " " << ax[FACEDIM][iface];
    for (i = 0; i < N_FACE_GRIDS; i++)
      cout << " " << fgr[iface][i];
    cout << endl;
  }
#if (BL_SPACEDIM == 3)
  cout << "Edges are:" << endl;
  for (iedge = 0; iedge < nbx[1]; iedge++) {
    cout << bx[1][iedge] << " " << ge[1][iedge] << " "
	 << flg[1][iedge] << " " << ax[1][iedge];
    for (i = 0; i < N_EDGE_GRIDS; i++)
      cout << " " << egr[iedge][i];
    cout << endl;
  }
#endif
  cout << "Corners are:" << endl;
  for (icor = 0; icor < nbx[0]; icor++) {
    cout << bx[0][icor] << " " << ge[0][icor] << " "
	 << flg[0][icor] << " " << ax[0][icor];
    for (i = 0; i < N_CORNER_GRIDS; i++)
      cout << " " << cgr[icor][i];
    cout << endl;
  }
*/
}

void level_interface::add(List<Box>& bl, Box b, int startgrid)
{
  Box ibox;
  IntVect t = b.type();
  for (int igrid = startgrid; igrid < im.length() + em.length(); igrid++) {
    if (igrid < im.length())
      ibox = im[igrid];
    else
      ibox = em[igrid-im.length()];
    ibox.convert(t);
    if (ibox.intersects(b) && !ibox.contains(b)) {
      for (int i = 0; i < BL_SPACEDIM; i++) {
	if (t[i] == IndexType::CELL) {
	  if (ibox.smallEnd(i) > b.smallEnd(i)) {
	    Box c = b.chop(i, ibox.smallEnd(i));
	    add(bl, b, igrid + 1);
	    add(bl, c, igrid);
	    return;
	  }
	  if (ibox.bigEnd(i) < b.bigEnd(i)) {
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

void level_interface::xfer(List<Box>& bl, int idim)
{
  nbx[idim] = bl.length();
  bx[idim]  = new Box[nbx[idim]];
  ge[idim]  = new unsigned[nbx[idim]];
  flg[idim] = new unsigned[nbx[idim]];

  Box btmp;
  int i = 0;
  //Boxnode *bn;
  //for (bn = bl.first(); bn; bn = bl.next(bn), i++) {
  ListIterator<Box> bn(bl);
  for ( ; bn; bn++, i++) {
    bx[idim][i] = bn();
    btmp = grow(bn(), bn().type()).convert(IntVect::TheCellVector());
    IntVect tmp = btmp.smallEnd();
    if (dom.contains(btmp)) {
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
    else {
      int is_in = dom.contains(tmp);
      ge[idim][i]  = ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp));
#if (BL_SPACEDIM == 2)
      tmp += IntVect(1,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 1;
      tmp += IntVect(-1,1);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 2;
      tmp += IntVect(1,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 3;
#else
      tmp += IntVect(1,0,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 1;
      tmp += IntVect(-1,1,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 2;
      tmp += IntVect(1,0,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 3;
      tmp += IntVect(-1,-1,1);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 4;
      tmp += IntVect(1,0,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 5;
      tmp += IntVect(-1,1,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 6;
      tmp += IntVect(1,0,0);
      is_in = dom.contains(tmp);
      ge[idim][i] |= ( is_in && im.contains(tmp) ||
		      !is_in && em.contains(tmp)) << 7;
#endif
    }
  }

  bl.clear();

  // Sort fine-fine boxes to beginning of list
  int j = nbx[idim];
  for (i = 0; i < j; i++) {
    if (ge[idim][i] != ALL) {
      while (--j > i) {
	if (ge[idim][j] == ALL) {
	  btmp = bx[idim][j];
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
  while (++j < nbx[idim] && ge[idim][j] == ALL);
  int nff = j;

  // Sort interior fine-fine boxes to beginning of list
  if (idim == 0) {
    for (i = 0; i < j; i++) {
      if (!bx[idim][i].intersects(idomain)) {
	while (--j > i) {
	  btmp = bx[idim][j];
	  if (btmp.intersects(idomain)) {
	    bx[idim][j] = bx[idim][i];
	    bx[idim][i] = btmp;
	    break;
	  }
	}
      }
    }
  }
  else {
    for (i = 0; i < j; i++) {
      btmp = bx[idim][i];
      btmp.convert(IntVect::TheNodeVector());
      if (!btmp.intersects(idomain)) {
	while (--j > i) {
	  btmp = bx[idim][j];
	  btmp.convert(IntVect::TheNodeVector());
	  if (btmp.intersects(idomain)) {
	    btmp = bx[idim][j];
	    bx[idim][j] = bx[idim][i];
	    bx[idim][i] = btmp;
	    break;
	  }
	}
      }
    }
  }

  if (idim == FACEDIM) {
    j = -1;
    while (++j < nff) {
      btmp = bx[idim][j];
      btmp.convert(IntVect::TheNodeVector());
      if (!btmp.intersects(idomain))
	break;
    }
    int nin = j;

    // Sort interior faces according to orientation, x first
    for (i = 0; i < j; i++) {
      if (bx[idim][i].type(0) == IndexType::CELL) {
	while (--j > i) {
	  if (bx[idim][j].type(0) == IndexType::NODE) {
	    btmp = bx[idim][j];
	    bx[idim][j] = bx[idim][i];
	    bx[idim][i] = btmp;
	    break;
	  }
	}
      }
    }
#if (BL_SPACEDIM == 3)
    j = nin;
    for (i = 0; i < j; i++) {
      if (bx[idim][i].type(0) == IndexType::CELL) {
	if (bx[idim][i].type(1) == IndexType::CELL) {
	  while (--j > i) {
	    if (bx[idim][j].type(1) == IndexType::NODE) {
	      btmp = bx[idim][j];
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
    for (i = nin; i < j; i++) {
      if (bx[idim][i].type(0) == IndexType::CELL) {
	while (--j > i) {
	  if (bx[idim][j].type(0) == IndexType::NODE) {
	    btmp = bx[idim][j];
	    bx[idim][j] = bx[idim][i];
	    bx[idim][i] = btmp;
	    break;
	  }
	}
      }
    }
#if (BL_SPACEDIM == 3)
    j = nff;
    for (i = nin; i < j; i++) {
      if (bx[idim][i].type(0) == IndexType::CELL) {
	if (bx[idim][i].type(1) == IndexType::CELL) {
	  while (--j > i) {
	    if (bx[idim][j].type(1) == IndexType::NODE) {
	      btmp = bx[idim][j];
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

void level_interface::geo_array(int ga[], int idim, int i) const
{
  unsigned gtmp = geo(idim, i);
  for (int k = 0; k < N_CORNER_GRIDS; k++) {
    ga[k] = (gtmp & 1);
    gtmp >>= 1;
  }
}
