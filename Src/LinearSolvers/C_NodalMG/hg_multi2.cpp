
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
#  ifdef HG_CROSS_STENCIL
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   int&, int&);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   const int*);
#  elif (defined HG_TERRAIN)
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   int&, int&);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
		   Real*, intS, Real*, intS, intS, intRS,
		   const int*);
#  else
  void FORT_HGFRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&,
		   const int&, int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, int&, const int&, const int&
#    endif
		   );
  void FORT_HGORES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&,
		   const int&, int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, int&, const int&
#    endif
		   );
  void FORT_HGIRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&,
		   const int&, int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, int&, const int&
#    endif
		   );
  void FORT_HGDRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&,
		   const int&, int&
#    else
		   Real*, intS, Real*, intS, intS, RealRS,
		   intRS, int&, const int&
#    endif
		   );
#  endif
#elif (BL_SPACEDIM == 3)
#  if (defined HG_TERRAIN)
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
#    ifdef HG_CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   int&, int&);
  void FORT_HGERES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
		   intS, Real&, const int&,
#    else
		   Real*, intS, Real*, intS, intS, RealRS, intRS,
#    endif
		   const int*, const int*);
  void FORT_HGCRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS,
#    ifdef HG_CONSTANT
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
#ifndef HG_CONSTANT
    fres_sfbox = new Box*[lev_max+1];
    fres_scbox = new Box*[lev_max+1];
    fres_sc = new PArray<FArrayBox>[lev_max+1];
#endif
    fres_dc = new PArray<FArrayBox>[lev_max+1];
    fres_flag = new int*[lev_max+1];
#if (BL_SPACEDIM == 3)
    eres_fbox  = new Box*[lev_max+1];
    eres_cbox  = new Box*[lev_max+1];
    eres_creg  = new Box*[lev_max+1];
#  ifndef HG_CONSTANT
    eres_sfbox = new Box*[lev_max+1];
    eres_scbox = new Box*[lev_max+1];
    eres_sf = new PArray<FArrayBox>[lev_max+1];
    eres_sc = new PArray<FArrayBox>[lev_max+1];
#  endif
    eres_df = new PArray<FArrayBox>[lev_max+1];
    eres_dc = new PArray<FArrayBox>[lev_max+1];
    eres_flag = new int*[lev_max+1];
#endif
    cres_fbox  = new Box*[lev_max+1];
    cres_cbox  = new Box*[lev_max+1];
    cres_creg  = new Box*[lev_max+1];
#ifndef HG_CONSTANT
    cres_sfbox = new Box*[lev_max+1];
    cres_scbox = new Box*[lev_max+1];
    cres_sf = new PArray<FArrayBox>[lev_max+1];
    cres_sc = new PArray<FArrayBox>[lev_max+1];
#endif
    cres_df = new PArray<FArrayBox>[lev_max+1];
    cres_dc = new PArray<FArrayBox>[lev_max+1];
    cres_flag = new int*[lev_max+1];
  }

  for (int lev = lev_min + 1; lev <= lev_max; lev++) {
    int mglev = ml_index[lev];
    fres_fbox[lev]  = new Box[lev_interface[mglev].nfaces()];
    fres_cbox[lev]  = new Box[lev_interface[mglev].nfaces()];
    fres_creg[lev]  = new Box[lev_interface[mglev].nfaces()];
#ifndef HG_CONSTANT
    fres_sfbox[lev] = new Box[lev_interface[mglev].nfaces()];
    fres_scbox[lev] = new Box[lev_interface[mglev].nfaces()];
    fres_sc[lev].resize(lev_interface[mglev].nfaces());
#endif
    fres_dc[lev].resize(lev_interface[mglev].nfaces());
    fres_flag[lev] = new int[lev_interface[mglev].nfaces()];
#if (BL_SPACEDIM == 3)
    eres_fbox[lev]  = new Box[lev_interface[mglev].nedges()];
    eres_cbox[lev]  = new Box[lev_interface[mglev].nedges()];
    eres_creg[lev]  = new Box[lev_interface[mglev].nedges()];
#  ifndef HG_CONSTANT
    eres_sfbox[lev] = new Box[lev_interface[mglev].nedges()];
    eres_scbox[lev] = new Box[lev_interface[mglev].nedges()];
    eres_sf[lev].resize(lev_interface[mglev].nedges());
    eres_sc[lev].resize(lev_interface[mglev].nedges());
#  endif
    eres_df[lev].resize(lev_interface[mglev].nedges());
    eres_dc[lev].resize(lev_interface[mglev].nedges());
    eres_flag[lev] = new int[lev_interface[mglev].nedges()];
#endif
    cres_fbox[lev]  = new Box[lev_interface[mglev].ncorners()];
    cres_cbox[lev]  = new Box[lev_interface[mglev].ncorners()];
    cres_creg[lev]  = new Box[lev_interface[mglev].ncorners()];
#ifndef HG_CONSTANT
    cres_sfbox[lev] = new Box[lev_interface[mglev].ncorners()];
    cres_scbox[lev] = new Box[lev_interface[mglev].ncorners()];
    cres_sf[lev].resize(lev_interface[mglev].ncorners());
    cres_sc[lev].resize(lev_interface[mglev].ncorners());
#endif
    cres_df[lev].resize(lev_interface[mglev].ncorners());
    cres_dc[lev].resize(lev_interface[mglev].ncorners());
    cres_flag[lev] = new int[lev_interface[mglev].ncorners()];
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
#ifndef HG_CONSTANT
    delete [] fres_sfbox[lev];
    delete [] fres_scbox[lev];
    for (i = 0; i < lev_interface[mglev].nfaces(); i++) {
      if (fres_flag[lev][i] == 0) delete fres_sc[lev].remove(i);
    }
#endif
    for (i = 0; i < lev_interface[mglev].nfaces(); i++) {
      if (fres_flag[lev][i] == 0) delete fres_dc[lev].remove(i);
    }
    delete [] fres_flag[lev];
#if (BL_SPACEDIM == 3)
    delete [] eres_fbox[lev];
    delete [] eres_cbox[lev];
    delete [] eres_creg[lev];
#  ifndef HG_CONSTANT
    delete [] eres_sfbox[lev];
    delete [] eres_scbox[lev];
    for (i = 0; i < lev_interface[mglev].nedges(); i++) {
      if (eres_sf[lev].defined(i)) delete eres_sf[lev].remove(i);
      if (eres_flag[lev][i] == 0)  delete eres_sc[lev].remove(i);
    }
#  endif
    for (i = 0; i < lev_interface[mglev].nedges(); i++) {
      if (eres_df[lev].defined(i)) delete eres_df[lev].remove(i);
      if (eres_flag[lev][i] == 0)  delete eres_dc[lev].remove(i);
    }
    delete [] eres_flag[lev];
#endif
    delete [] cres_fbox[lev];
    delete [] cres_cbox[lev];
    delete [] cres_creg[lev];
#ifndef HG_CONSTANT
    delete [] cres_sfbox[lev];
    delete [] cres_scbox[lev];
    for (i = 0; i < lev_interface[mglev].ncorners(); i++) {
      if (cres_sf[lev].defined(i)) delete cres_sf[lev].remove(i);
      if (cres_flag[lev][i] == 0)  delete cres_sc[lev].remove(i);
    }
#endif
    for (i = 0; i < lev_interface[mglev].ncorners(); i++) {
      if (cres_df[lev].defined(i)) delete cres_df[lev].remove(i);
      if (cres_flag[lev][i] == 0)  delete cres_dc[lev].remove(i);
    }
    delete [] cres_flag[lev];
  }
  if (lev_min < lev_max) {
    delete [] fres_fbox;
    delete [] fres_cbox;
    delete [] fres_creg;
#ifndef HG_CONSTANT
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
#  ifndef HG_CONSTANT
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
#ifndef HG_CONSTANT
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

  const IntVect& rat = gen_ratio[lev-1];
  int mglevc = ml_index[lev-1];

#ifdef HG_TERRAIN
  int ncomp = 2 * BL_SPACEDIM - 1;
  const amr_boundary_class& bndry = boundary.terrain_sigma();
#else
  int ncomp = 1;
  const amr_boundary_class& bndry = boundary.scalar();
#endif

  for (int iface = 0; iface < lev_interface[mglev].nfaces(); iface++) {
    // find a fine grid touching this face
    igrid = lev_interface[mglev].fgrid(iface, 0);
    if (igrid < 0)
      igrid = lev_interface[mglev].fgrid(iface, 1);
    unsigned geo = lev_interface[mglev].fgeo(iface);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].fflag(iface) == 1) {
      fres_flag[lev][iface] = 1; // means "fres_dc[lev][iface] not allocated"
      continue;
    }
    // fine grid on just one side
    int idim = lev_interface[mglev].fdim(iface);
    int idir = (geo & level_interface::LOW) ? -1 : 1;
    Box& fbox = fres_fbox[lev][iface];
    Box& cbox = fres_cbox[lev][iface];
    Box& creg = fres_creg[lev][iface];
    fbox = dest[lev][igrid].box();
    cbox = lev_interface[mglev].node_face(iface);
    cbox.coarsen(rat);
    if (idir > 0)
      cbox.growLo(idim, 1);
    else
      cbox.growHi(idim, 1);
    int cgrid = find_patch(cbox, dest[lev-1]);
#ifndef HG_CONSTANT
    Box& sigmafbox = fres_sfbox[lev][iface];
    Box& sigmacbox = fres_scbox[lev][iface];
    sigmafbox = sigma[mglev][igrid].box();
    sigmacbox = cbox;
    sigmacbox.convert(IntVect::TheCellVector());
    if (cgrid >= 0) {
      fres_sc[lev].set(iface, &sigma[mglevc][cgrid]);
      sigmacbox = fres_sc[lev][iface].box();
    }
    else {
      fres_sc[lev].set(iface, new FArrayBox(sigmacbox, ncomp));
      fill_patch(fres_sc[lev][iface], sigma[mglevc],
		 lev_interface[mglevc], bndry);
    }
#endif
    if (cgrid >= 0) {
      fres_dc[lev].set(iface, &dest[lev-1][cgrid]);
      cbox = fres_dc[lev][iface].box();
      fres_flag[lev][iface] = 1;
    }
    else {
      fres_dc[lev].set(iface, new FArrayBox(cbox));
      fres_flag[lev][iface] = 0;
    }
    IntVect t = lev_interface[mglev].face(iface).type();
    creg = lev_interface[mglev].node_face(iface);
    creg.coarsen(rat).grow(t - IntVect::TheUnitVector());
  }

#if (BL_SPACEDIM == 3)

  for (int iedge = 0; iedge < lev_interface[mglev].nedges(); iedge++) {
    // find a fine grid touching this edge
      for (i = 0; i < level_interface::N_EDGE_GRIDS; i++) {
      igrid = lev_interface[mglev].egrid(iedge, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = lev_interface[mglev].egeo(iedge);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].eflag(iedge) == 1) {
      eres_flag[lev][iedge] = 1; // means "eres_dc[lev][iedge] not allocated"
      continue;
    }
    Box& fbox = eres_fbox[lev][iedge];
    Box& cbox = eres_cbox[lev][iedge];
    Box& creg = eres_creg[lev][iedge];
    IntVect t = lev_interface[mglev].edge(iedge).type();
    cbox = lev_interface[mglev].node_edge(iedge);
    cbox.coarsen(rat).grow(t);
    fbox = refine(cbox, rat);
    eres_df[lev].set(iedge, new FArrayBox(fbox));
    int cgrid = find_patch(cbox, dest[lev-1]);
#  ifndef HG_CONSTANT
    Box& sigmafbox = eres_sfbox[lev][iedge];
    Box& sigmacbox = eres_scbox[lev][iedge];
    sigmafbox = fbox;
    sigmafbox.convert(IntVect::TheCellVector());
    eres_sf[lev].set(iedge, new FArrayBox(sigmafbox, ncomp));
    fill_patch(eres_sf[lev][iedge], sigma[mglev], lev_interface[mglev],
	       bndry, 0, 1, iedge);
    sigmacbox = cbox;
    sigmacbox.convert(IntVect::TheCellVector());
    if (cgrid >= 0) {
      eres_sc[lev].set(iedge, &sigma[mglevc][cgrid]);
      sigmacbox = eres_sc[lev][iedge].box();
    }
    else {
      eres_sc[lev].set(iedge, new FArrayBox(sigmacbox, ncomp));
      fill_patch(eres_sc[lev][iedge], sigma[mglevc],
		 lev_interface[mglevc], bndry);
    }
#  endif
    if (cgrid >= 0) {
      eres_dc[lev].set(iedge, &dest[lev-1][cgrid]);
      cbox = eres_dc[lev][iedge].box();
      eres_flag[lev][iedge] = 1;
    }
    else {
      eres_dc[lev].set(iedge, new FArrayBox(cbox));
      eres_flag[lev][iedge] = 0;
    }
    creg = lev_interface[mglev].node_edge(iedge);
    creg.coarsen(rat).grow(t - IntVect::TheUnitVector());
  }

#endif

  for (int icor = 0; icor < lev_interface[mglev].ncorners(); icor++) {
    // find a fine grid touching this corner
    for (i = 0; i < level_interface::N_CORNER_GRIDS; i++) {
      igrid = lev_interface[mglev].cgrid(icor, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = lev_interface[mglev].cgeo(icor);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].cflag(icor) == 1) {
      cres_flag[lev][icor] = 1; // means "cres_dc[lev][icor] not allocated"
      continue;
    }
    Box& fbox = cres_fbox[lev][icor];
    Box& cbox = cres_cbox[lev][icor];
    Box& creg = cres_creg[lev][icor];
    cbox = lev_interface[mglev].corner(icor);
    fbox = lev_interface[mglev].corner(icor);
    cbox.coarsen(rat).grow(1);
    fbox.grow(rat);
    cres_df[lev].set(icor, new FArrayBox(fbox));
    int cgrid = find_patch(cbox, dest[lev-1]);
#ifndef HG_CONSTANT
    Box& sigmafbox = cres_sfbox[lev][icor];
    Box& sigmacbox = cres_scbox[lev][icor];
    sigmafbox = fbox;
    sigmafbox.convert(IntVect::TheCellVector());
    cres_sf[lev].set(icor, new FArrayBox(sigmafbox, ncomp));
    fill_patch(cres_sf[lev][icor], sigma[mglev], lev_interface[mglev],
	       bndry, 0, 0, icor);
    sigmacbox = cbox;
    sigmacbox.convert(IntVect::TheCellVector());
    if (cgrid >= 0) {
      cres_sc[lev].set(icor, &sigma[mglevc][cgrid]);
      sigmacbox = cres_sc[lev][icor].box();
    }
    else {
      cres_sc[lev].set(icor, new FArrayBox(sigmacbox, ncomp));
      fill_patch(cres_sc[lev][icor], sigma[mglevc],
		 lev_interface[mglevc], bndry);
    }
#endif
    if (cgrid >= 0) {
      cres_dc[lev].set(icor, &dest[lev-1][cgrid]);
      cbox = cres_dc[lev][icor].box();
      cres_flag[lev][icor] = 1;
    }
    else {
      cres_dc[lev].set(icor, new FArrayBox(cbox));
      cres_flag[lev][icor] = 0;
    }
    creg = lev_interface[mglev].corner(icor);
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

  const IntVect& rat = gen_ratio[lev-1];
  int mglevc = ml_index[lev-1];

  for (int iface = 0; iface < lev_interface[mglev].nfaces(); iface++) {
    // find a fine grid touching this face
    igrid = lev_interface[mglev].fgrid(iface, 0);
    if (igrid < 0)
      igrid = lev_interface[mglev].fgrid(iface, 1);
    unsigned geo = lev_interface[mglev].fgeo(iface);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].fflag(iface) == 1)
      continue;
    // fine grid on just one side
    int idim = lev_interface[mglev].fdim(iface);
    int idir = (geo & level_interface::LOW) ? -1 : 1;
    const Box& sbox = source[lev][igrid].box();
    const Box& fbox = fres_fbox[lev][iface];
    const Box& cbox = fres_cbox[lev][iface];
#ifndef HG_CONSTANT
    const Box& sigmafbox = fres_sfbox[lev][iface];
    const Box& sigmacbox = fres_scbox[lev][iface];
    FArrayBox& sigmac = fres_sc[lev][iface];
    Real *const sigmafptr = sigma[mglev][igrid].dataPtr();
#endif
    const Box& creg = fres_creg[lev][iface];
    FArrayBox& cdst = fres_dc[lev][iface];
    if (fres_flag[lev][iface] == 0) {
      fill_patch(cdst, dest[lev-1], lev_interface[mglevc], boundary.pressure());
    }
    Real *const rptr = resid[mglev][igrid].dataPtr();
    Real *const sptr = source[lev][igrid].dataPtr();
    Real *const dptr = dest[lev][igrid].dataPtr();
    FORT_HGFRES(rptr, dimlist(sbox),
		sptr, dimlist(sbox),
		dptr, dimlist(fbox),
		cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		sigmafptr, dimlist(sigmafbox),
		sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		dimlist(creg),
#ifdef HG_CONSTANT
		hx,
		rat[0], idim, idir
#elif (defined HG_TERRAIN)
		D_DECL(rat[0], rat[1], rat[2]), idim, idir
#elif (BL_SPACEDIM == 2)
		hx, hy,
#  ifdef HG_CROSS_STENCIL
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

#if (defined HG_CROSS_STENCIL) || (defined HG_TERRAIN)

  int ga[level_interface::N_CORNER_GRIDS];

#if (BL_SPACEDIM == 3)

  for (int iedge = 0; iedge < lev_interface[mglev].nedges(); iedge++) {
    // find a fine grid touching this edge
    for (i = 0; i < level_interface::N_EDGE_GRIDS; i++) {
      igrid = lev_interface[mglev].egrid(iedge, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = lev_interface[mglev].egeo(iedge);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo != level_interface::ALL && igrid >= 0 && lev_interface[mglev].eflag(iedge) == 0) {
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = eres_fbox[lev][iedge];
      const Box& cbox = eres_cbox[lev][iedge];
#ifndef HG_CONSTANT
      const Box& sigmafbox = eres_sfbox[lev][iedge];
      const Box& sigmacbox = eres_scbox[lev][iedge];
      FArrayBox& sigmaf = eres_sf[lev][iedge];
      FArrayBox& sigmac = eres_sc[lev][iedge];
#endif
      const Box& creg = eres_creg[lev][iedge];
      IntVect t = lev_interface[mglev].edge(iedge).type();
      FArrayBox& fdst = eres_df[lev][iedge];
      fill_patch(fdst, dest[lev], lev_interface[mglev], boundary.pressure(),
                 0, 1, iedge);
      FArrayBox& cdst = eres_dc[lev][iedge];
      if (eres_flag[lev][iedge] == 0) {
	fill_patch(cdst, dest[lev-1], lev_interface[mglevc], boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      lev_interface[mglev].geo_array(ga, 1, iedge);
      FORT_HGERES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef HG_CONSTANT
                  hx, rat[0],
#elif (defined HG_TERRAIN)
		  rat[0], rat[1], rat[2],
#else
		  hx, hy, hz, rat[0], rat[1], rat[2],
#endif
		  t.getVect(), ga);
      // fill in the grids on the other sides, if any
      const Box& freg = lev_interface[mglev].node_edge(iedge);
      for (i = 1; i < level_interface::N_EDGE_GRIDS; i++) {
	int jgrid = lev_interface[mglev].egrid(iedge, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
  }

#endif

  for (int icor = 0; icor < lev_interface[mglev].ncorners(); icor++) {
    // find a fine grid touching this corner
    for (i = 0; i < level_interface::N_CORNER_GRIDS; i++) {
      igrid = lev_interface[mglev].cgrid(icor, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = lev_interface[mglev].cgeo(icor);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo != level_interface::ALL && igrid >= 0 && lev_interface[mglev].cflag(icor) == 0) {
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef HG_CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      FArrayBox& sigmaf = cres_sf[lev][icor];
      FArrayBox& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      FArrayBox& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], lev_interface[mglev], boundary.pressure(),
		 0, 0, icor);
      FArrayBox& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], lev_interface[mglevc],
		   boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      lev_interface[mglev].geo_array(ga, 0, icor);
      FORT_HGCRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef HG_CONSTANT
                  hx, rat[0],
#elif (defined HG_TERRAIN)
		  D_DECL(rat[0], rat[1], rat[2]),
#else
		  D_DECL(hx, hy, hz), D_DECL(rat[0], rat[1], rat[2]),
#endif
		  ga);
      // fill in the grids on the other sides, if any
      const Box& freg = lev_interface[mglev].corner(icor);
      for (i = 1; i < level_interface::N_CORNER_GRIDS; i++) {
	int jgrid = lev_interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
  }

#else // 2d 9-point stencils:

  for (int icor = 0; icor < lev_interface[mglev].ncorners(); icor++) {
    // find a fine grid touching this corner
    for (i = 0; i < level_interface::N_CORNER_GRIDS; i++) {
      igrid = lev_interface[mglev].cgrid(icor, i);
      if (igrid >= 0)
	break;
    }
    unsigned geo = lev_interface[mglev].cgeo(icor);
    // reject fine-fine interfaces and those without an interior fine grid
    if (geo == level_interface::ALL || igrid < 0 || lev_interface[mglev].cflag(icor) == 1)
      continue;
    else if (geo == level_interface::XL || geo == level_interface::XH || geo == level_interface::YL || geo == level_interface::YH) {
      // fine grid on two adjacent sides
      int idim = (geo == level_interface::XL || geo == level_interface::XH) ? 0 : 1;
      int idir = (geo & level_interface::LL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef HG_CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      FArrayBox& sigmaf = cres_sf[lev][icor];
      FArrayBox& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      FArrayBox& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], lev_interface[mglev], boundary.pressure(),
		 0, 0, icor);
      FArrayBox& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], lev_interface[mglevc],
		   boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGFRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef HG_CONSTANT
                  hx,
                  rat[0], idim, idir
#else
                  hx, hy,
                  rat[0], rat[1], idim, idir,
		  IsRZ(), mg_domain[mglevc].bigEnd(0) + 1
#endif
		  );
      // fill in the grids on the other sides, if any
      const Box& freg = lev_interface[mglev].corner(icor);
      for (i = 1; i < level_interface::N_CORNER_GRIDS; i++) {
	int jgrid = lev_interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
    else if (geo == level_interface::LL || geo == level_interface::HL || geo == level_interface::LH || geo == level_interface::HH) {
      // outside corner
      int idir0 = (geo & level_interface::XL) ? -1 : 1;
      int idir1 = (geo & level_interface::YL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      const Box& cbox = cres_cbox[lev][icor];
#ifndef HG_CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      FArrayBox& sigmaf = cres_sf[lev][icor];
      FArrayBox& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      FArrayBox& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], lev_interface[mglevc], boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      const Box& fbox = dest[lev][igrid].box();
      Real *const dptr = dest[lev][igrid].dataPtr();
      FORT_HGORES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  dptr, dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef HG_CONSTANT
                  hx,
                  rat[0], idir0, idir1
#else
                  hx, hy,
                  rat[0], rat[1], idir0, idir1, IsRZ()
#endif
		  );
    }
    else if (geo == (level_interface::LL | level_interface::HH) || geo == (level_interface::LH | level_interface::HL)) {
      // diagonal corner
      int jdir = (geo == (level_interface::LL | level_interface::HH)) ? 1 : -1;
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef HG_CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      FArrayBox& sigmaf = cres_sf[lev][icor];
      FArrayBox& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      FArrayBox& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], lev_interface[mglev], boundary.pressure(),
		 0, 0, icor);
      FArrayBox& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], lev_interface[mglevc],
		   boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGDRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef HG_CONSTANT
                  hx,
                  rat[0], jdir
#else
                  hx, hy,
                  rat[0], rat[1], jdir, IsRZ()
#endif
		  );
      // fill in the grids on the other sides, if any
      const Box& freg = lev_interface[mglev].corner(icor);
      for (i = 1; i < level_interface::N_CORNER_GRIDS; i++) {
	int jgrid = lev_interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid)
	  internal_copy(resid[mglev], jgrid, igrid, freg);
      }
    }
    else {
      // inside corner
      int idir0 = ((geo & level_interface::XL) == level_interface::XL) ? -1 : 1;
      int idir1 = ((geo & level_interface::YL) == level_interface::YL) ? -1 : 1;
      const Box& sbox = source[lev][igrid].box();
      const Box& fbox = cres_fbox[lev][icor];
      const Box& cbox = cres_cbox[lev][icor];
#ifndef HG_CONSTANT
      const Box& sigmafbox = cres_sfbox[lev][icor];
      const Box& sigmacbox = cres_scbox[lev][icor];
      FArrayBox& sigmaf = cres_sf[lev][icor];
      FArrayBox& sigmac = cres_sc[lev][icor];
#endif
      const Box& creg = cres_creg[lev][icor];
      FArrayBox& fdst = cres_df[lev][icor];
      fill_patch(fdst, dest[lev], lev_interface[mglev], boundary.pressure(),
		 0, 0, icor);
      FArrayBox& cdst = cres_dc[lev][icor];
      if (cres_flag[lev][icor] == 0) {
	fill_patch(cdst, dest[lev-1], lev_interface[mglevc], boundary.pressure());
      }
      Real *const rptr = resid[mglev][igrid].dataPtr();
      Real *const sptr = source[lev][igrid].dataPtr();
      FORT_HGIRES(rptr, dimlist(sbox),
		  sptr, dimlist(sbox),
		  fdst.dataPtr(), dimlist(fbox),
		  cdst.dataPtr(), dimlist(cbox),
#ifndef HG_CONSTANT
		  sigmaf.dataPtr(), dimlist(sigmafbox),
		  sigmac.dataPtr(), dimlist(sigmacbox),
#endif
		  dimlist(creg),
#ifdef HG_CONSTANT
                  hx,
                  rat[0], idir0, idir1
#else
                  hx, hy,
                  rat[0], rat[1], idir0, idir1, IsRZ()
#endif
		  );
      // fill in the grids on the other sides, if any
      const Box& freg = lev_interface[mglev].corner(icor);
      int kgrid = -1;
      for (i = 1; i < level_interface::N_CORNER_GRIDS; i++) {
	int jgrid = lev_interface[mglev].cgrid(icor, i);
	if (jgrid >= 0 && jgrid != igrid && jgrid != kgrid) {
	  internal_copy(resid[mglev], jgrid, igrid, freg);
	  kgrid = jgrid;
	}
      }
    }
  }
#endif
}
