
#include "hg_multi.H"

#ifdef CONSTANT
#  define CGOPT 2
#else
#  define CGOPT 1
#endif

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HGRES      hgres_
#  define   FORT_HGRESU     hgresu_
#  define   FORT_HGRLX      hgrlx_
#  define   FORT_HGRLXU     hgrlxu_
#  define   FORT_HGCG       hgcg_
#  define   FORT_HGCG1      hgcg1_
#  define   FORT_HGCG2      hgcg2_
#  define   FORT_HGIP       hgip_
#else
#  define   FORT_HGRES      HGRES
#  define   FORT_HGRESU     HGRESU
#  define   FORT_HGRLX      HGRLX
#  define   FORT_HGRLXU     HGRLXU
#  define   FORT_HGCG       HGCG
#  define   FORT_HGCG1      HGCG1
#  define   FORT_HGCG2      HGCG2
#  define   FORT_HGIP       HGIP
#endif

extern "C" {

#if (BL_SPACEDIM == 1)
  ERROR, not relevant
#else
#  ifdef CONSTANT
  void FORT_HGRES(Real*, intS, Real*, Real*, intS, Real&);
  void FORT_HGRESU(Real*, intS, Real*, Real*, intS, Real&);
#    ifdef CROSS_STENCIL
  void FORT_HGRLXU(Real*, Real*, intS, Real*, intS, Real&);
#    else
  void FORT_HGRLX(Real*, Real*, intS, Real*, intS, Real&);
#    endif
#  else
#    ifdef SIGMA_NODE
  void FORT_HGRES(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
  void FORT_HGRESU(Real*, intS, Real*, Real*, Real*, Real*, intS);
  void FORT_HGRLX(Real*, intS, Real*, intS, Real*, intS, Real*, intS, intS);
  void FORT_HGRLXU(Real*, Real*, Real*, Real*, intS, Real*, intS);
#    else
  void FORT_HGRES(Real*, intS, Real*, intS, Real*, intS, RealPS, intS, intS,
		  RealRS, const int&, const int&);
  void FORT_HGRLX(Real*, intS, Real*, intS, RealPS, intS, Real*, intS, intS,
		  RealRS, const int&, const int&);
#    endif
#  endif

#  if (CGOPT == 1)
  void FORT_HGCG1(Real*, Real*, Real*, Real*, Real*, Real*, Real*,
		  intS, const Real&, Real&);
  void FORT_HGCG2(Real*, Real*, intS, const Real&);
  void FORT_HGIP(Real*, Real*, Real*, intS, Real&);
#  elif (CGOPT == 2)
#    if (BL_SPACEDIM == 2)
  void FORT_HGCG(Real*, Real*, Real*, Real*, Real*, Real*, Real*,
		 const int&, int*, int*,
		 int*, int*, int*, int*, int*, int*, int*,
		 const int&, Real*,
		 int*, int*, int*, int*, int*,
		 const Real&, Real&, Real&, int&, const int&);
#    else
  void FORT_HGCG(Real*, Real*, Real*, Real*, Real*, Real*, Real*,
		 const int&, int*, int*, int*,
		 int*, int*, int*, int*, int*, int*, int*,
		 const int&, Real*,
		 int*, int*, int*, int*, int*, int*, int*, int*,
		 const Real&, Real&, Real&, int&, const int&);
#    endif
#  endif
#endif
}

void holy_grail_amr_multigrid::level_residual(MultiFab& r,
					      MultiFab& s,
					      MultiFab& d,
					      copy_cache* dbc,
					      int mglev,
					      int iclear)
{
  assert(r.boxArray() == s.boxArray());
  assert(r.boxArray() == d.boxArray());

  int igrid;

  fill_borders(d, dbc, interface[mglev], mg_boundary);

#ifdef SIGMA_NODE

  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    const Box& rbox = r[igrid].box();
    const Box& freg = interface[mglev].part_fine(igrid);
    FORT_HGRESU(r[igrid].dataPtr(), dimlist(rbox),
                s[igrid].dataPtr(),
                d[igrid].dataPtr(),
                sigma_node[mglev][igrid].dataPtr(),
                mask[mglev][igrid].dataPtr(),
                dimlist(freg));
  }

#else

  Real hx = h[mglev][0];
  Real hy = h[mglev][1];
#  if (BL_SPACEDIM == 3)
  Real hz = h[mglev][2];
#  endif

#  ifdef CONSTANT

  if (!iclear) {
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& rbox = r[igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
      FORT_HGRESU(r[igrid].dataPtr(), dimlist(rbox),
		  s[igrid].dataPtr(), d[igrid].dataPtr(),
		  dimlist(freg), hx);
    }
  }
  else {
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& rbox = r[igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
      FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
		 s[igrid].dataPtr(), d[igrid].dataPtr(),
		 dimlist(freg), hx);
    }
    clear_part_interface(r, interface[mglev]);
  }

#  else

  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    const Box& rbox = r[igrid].box();
    const Box& sbox = s[igrid].box();
    const Box& dbox = d[igrid].box();
    const Box& freg = interface[mglev].part_fine(igrid);
#    ifndef SIGMA_NODE
    // this branch is the only one that can be reached here
    const Box& sigbox = sigma[mglev][igrid].box();
    FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
	       s[igrid].dataPtr(), dimlist(sbox),
	       d[igrid].dataPtr(), dimlist(dbox),
	       sigma_nd[0][mglev][igrid].dataPtr(),
#      if (BL_SPACEDIM == 2)
	       sigma_nd[1][mglev][igrid].dataPtr(), dimlist(sigbox),
	       dimlist(freg), hx, hy,
	       IsRZ(), mg_mesh[mglev].domain().bigEnd(0) + 1
#      else
	       sigma_nd[1][mglev][igrid].dataPtr(),
	       sigma_nd[2][mglev][igrid].dataPtr(), dimlist(sigbox),
               dimlist(freg), hx, hy, hz
#      endif
	       );
#    else
    // this branch is unreachable
    const Box& sigbox = sigma_node[mglev][igrid].box();
    FORT_HGRES(r[igrid].dataPtr(), dimlist(rbox),
               s[igrid].dataPtr(), dimlist(sbox),
               d[igrid].dataPtr(), dimlist(dbox),
               sigma_node[mglev][igrid].dataPtr(), dimlist(sigbox),
               dimlist(freg));
#    endif // SIGMA_NODE
  }

  if (iclear) {
    clear_part_interface(r, interface[mglev]);
  }

#  endif // CONSTANT
#endif // SIGMA_NODE
}

void holy_grail_amr_multigrid::relax(int mglev, int i1, int is_zero)
{
  Real hx = h[mglev][0];
  Real hy = h[mglev][1];
#if (BL_SPACEDIM == 3)
  Real hz = h[mglev][2];
#endif

  DECLARE_GEOMETRY_TYPES;

  for (int icount = 0; icount < i1; icount++) {
    if (is_zero == 0)
      fill_borders(corr[mglev], corr_bcache[mglev],
		   interface[mglev], mg_boundary);
    else
      is_zero = 0;
    for (int igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& sbox = resid[mglev][igrid].box();
      const Box& freg = interface[mglev].part_fine(igrid);
#ifdef CONSTANT
#  ifdef CROSS_STENCIL
      FORT_HGRLXU(corr[mglev][igrid].dataPtr(),
		  resid[mglev][igrid].dataPtr(), dimlist(sbox),
		  mask[mglev][igrid].dataPtr(),
		  dimlist(freg), hx);
#  else
      FORT_HGRLX(corr[mglev][igrid].dataPtr(),
		 resid[mglev][igrid].dataPtr(), dimlist(sbox),
		 mask[mglev][igrid].dataPtr(),
		 dimlist(freg), hx);
#  endif
#else
#ifdef SIGMA_NODE
/*
      const Box& fbox = corr[mglev][igrid].box();
      const Box& cenbox = cen[mglev][igrid].box();
      const Box& sigbox = sigma_node[mglev][igrid].box();
      FORT_HGRLX(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		 resid[mglev][igrid].dataPtr(), dimlist(sbox),
		 sigma_node[mglev][igrid].dataPtr(), dimlist(sigbox),
		 cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		 dimlist(freg));
*/
      FORT_HGRLXU(corr[mglev][igrid].dataPtr(),
		  resid[mglev][igrid].dataPtr(),
		  sigma_node[mglev][igrid].dataPtr(),
		  cen[mglev][igrid].dataPtr(), dimlist(sbox),
                  mask[mglev][igrid].dataPtr(),
                  dimlist(freg));
#else
      const Box& fbox = corr[mglev][igrid].box();
      const Box& cenbox = cen[mglev][igrid].box();
      const Box& sigbox = sigma[mglev][igrid].box();
      FORT_HGRLX(corr[mglev][igrid].dataPtr(), dimlist(fbox),
		 resid[mglev][igrid].dataPtr(), dimlist(sbox),
		 sigma_nd[0][mglev][igrid].dataPtr(),
#  if (BL_SPACEDIM == 2)
		 sigma_nd[1][mglev][igrid].dataPtr(), dimlist(sigbox),
		 cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		 dimlist(freg), hx, hy,
		 IsRZ(), mg_mesh[mglev].domain().bigEnd(0) + 1
#  else
		 sigma_nd[1][mglev][igrid].dataPtr(),
		 sigma_nd[2][mglev][igrid].dataPtr(), dimlist(sigbox),
		 cen[mglev][igrid].dataPtr(), dimlist(cenbox),
		 dimlist(freg), hx, hy, hz
#  endif
		 );
#endif
#endif
    }
    sync_borders(corr[mglev], corr_scache[mglev],
		 interface[mglev], mg_boundary);
  }
}

void holy_grail_amr_multigrid::cgsolve(int mglev)
{
  assert(mglev == 0);

  MultiFab& r = cgwork[0];
  MultiFab& p = cgwork[1];
  MultiFab& z = cgwork[2];
  MultiFab& x = cgwork[3];
  MultiFab& w = cgwork[4];
  MultiFab& c = cgwork[5];
  MultiFab& zero_array = cgwork[6];
  MultiFab& ipmask = cgwork[7];

  unroll_cache& ruc = *cgw_ucache[0];
  unroll_cache& puc = *cgw_ucache[1];
  unroll_cache& zuc = *cgw_ucache[2];
  unroll_cache& xuc = *cgw_ucache[3];
  unroll_cache& wuc = *cgw_ucache[4];
  unroll_cache& cuc = *cgw_ucache[5];
  unroll_cache& muc = *cgw_ucache[7];

  copy_cache& pbc = *cgw1_bcache;

  Real alpha, rho;
  int i = 0, igrid;

  // x (corr[0]) should be all 0.0 at this point
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    r[igrid].copy(resid[mglev][igrid]);
    r[igrid].negate();
  }
  //r.copy(resid[mglev]);
  //r.negate();

  if (singular) {
    // singular systems are very sensitive to solvability
    w.setVal(1.0);
    alpha = inner_product(r, w) / mg_domain[mglev].volume();
    r.plus(-alpha, 0);
  }

#if (CGOPT == 2)
  FORT_HGCG(ruc.ptr, puc.ptr,
	    zuc.ptr, xuc.ptr,
	    wuc.ptr, cuc.ptr,
	    muc.ptr, mg_mesh[0].length(),
#  if (BL_SPACEDIM == 2)
	    ruc.strid, ruc.nvals,
#  else
	    ruc.strid1, ruc.strid2,
	    ruc.nvals,
#  endif
	    ruc.start, puc.start,
	    zuc.start, xuc.start,
	    wuc.start, cuc.start,
	    muc.start,
	    pbc.nsets, pbc.dptr,
#  if (BL_SPACEDIM == 2)
	    pbc.nvals,
	    pbc.dstart, pbc.sstart,
	    pbc.dstrid, pbc.sstrid,
#  else
	    pbc.nvals1, pbc.nvals2,
	    pbc.dstart, pbc.sstart,
	    pbc.dstrid1, pbc.dstrid2,
	    pbc.sstrid1, pbc.sstrid2,
#  endif
	    h[0][0], alpha, rho, i, pcode);
#elif (CGOPT == 1)
  rho = 0.0;
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    z[igrid].copy(r[igrid]);
    z[igrid].mult(c[igrid]);
    const Box& reg = p[igrid].box();
    FORT_HGIP(z[igrid].dataPtr(), r[igrid].dataPtr(),
	      ipmask[igrid].dataPtr(),
	      dimlist(reg), rho);
    p[igrid].copy(z[igrid]);
  }
  Real tol = 1.e-3 * rho;

  while (tol > 0.0) {
    i++;
    if (i > 250)
      BoxLib::Error("holy_grail_amr_multigrid::cgsolve---conjugate-gradient iteration failed");
    Real rho_old = rho;
    // safe to set the clear flag to 0 here---bogus values make it
    // into r but are cleared from z by the mask in c
    level_residual(w, zero_array, p, &pbc, 0, 0);
    alpha = 0.0;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& reg = p[igrid].box();
      FORT_HGIP(p[igrid].dataPtr(), w[igrid].dataPtr(),
		ipmask[igrid].dataPtr(),
		dimlist(reg), alpha);
    }
    alpha = rho / alpha;
    rho = 0.0;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& reg = p[igrid].box();
      FORT_HGCG1(r[igrid].dataPtr(), p[igrid].dataPtr(),
		 z[igrid].dataPtr(), x[igrid].dataPtr(),
		 w[igrid].dataPtr(), c[igrid].dataPtr(),
		 ipmask[igrid].dataPtr(), dimlist(reg), alpha, rho);
    }
    if (pcode >= 3)
      cout << i << " " << rho << endl;
    if (rho <= tol)
      break;
    alpha = rho / rho_old;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      const Box& reg = p[igrid].box();
      FORT_HGCG2(p[igrid].dataPtr(), z[igrid].dataPtr(),
		 dimlist(reg), alpha);
    }
  }
#else
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    z[igrid].copy(r[igrid]);
    z[igrid].mult(c[igrid]);
  }
  //z.assign(r).mult(c);
  rho = inner_product(z, r);
  Real tol = 1.e-3 * rho;
  //p.assign(0.0);
  for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
    p[igrid].copy(z[igrid]);
  }

  while (tol > 0.0) {
    i++;
    if (i > 250)
      BoxLib::Error("holy_grail_amr_multigrid::cgsolve---conjugate-gradient iteration failed");
    Real rho_old = rho;
    // safe to set the clear flag to 0 here---bogus values make it
    // into r but are cleared from z by the mask in c
    level_residual(w, zero_array, p, &pbc, 0, 0);
    alpha = rho / inner_product(p, w);
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      w[igrid].mult(alpha);
      r[igrid].minus(w[igrid]);
      w[igrid].copy(p[igrid]);
      w[igrid].mult(alpha);
      x[igrid].plus(w[igrid]);
      z[igrid].copy(r[igrid]);
      z[igrid].mult(c[igrid]);
    }
    //r.minus(w.mult(alpha));
    //x.plus(w.assign(p).mult(alpha));
    //z.assign(r).mult(c);
    rho = inner_product(z, r);
    if (pcode >= 3)
      cout << i << " " << rho << endl;
    if (rho <= tol)
      break;
    for (igrid = 0; igrid < mg_mesh[mglev].length(); igrid++) {
      p[igrid].mult(rho / rho_old);
      p[igrid].plus(z[igrid]);
    }
    //p.mult(rho / rho_old).plus(z);
  }
#endif

  if (pcode >= 2)
    cout << i << " iterations required for conjugate-gradient" << endl;
}
