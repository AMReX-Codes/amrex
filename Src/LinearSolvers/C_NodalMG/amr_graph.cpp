
#include "amr_graph.H"

const Box      null_Box;
const BoxArray null_level_mesh;

#ifdef BL_FORT_USE_UNDERSCORE
#  define FCONT2D cont2d_
#  define FCONT3D cont3d_
#  define FGETSL0 getsl0_
#  define FGETSL1 getsl1_
#  define FGETSL2 getsl2_
#else
#  define FCONT2D CONT2D
#  define FCONT3D CONT3D
#  define FGETSL0 GETSL0
#  define FGETSL1 GETSL1
#  define FGETSL2 GETSL2
#endif

extern "C" void FCONT2D(const Real*,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			Real&, Real&, Real&, Real&,
			const int&, Real&, const int&);

#if (BL_SPACEDIM == 3)

extern "C" void FCONT3D(const Real*,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			const int&,
			const int&, const int&, const int&,
			const int&, const int&, const int&,
			Real&, Real&, Real&, Real&, Real&, Real&, const int&);

extern "C" void FGETSL0(const Real*,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			Real*, const int&);
extern "C" void FGETSL1(const Real*,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			Real*, const int&);
extern "C" void FGETSL2(const Real*,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			const int&, const int&, const int&, const int&,
			Real*, const int&);

static Real theta, phi, ct, st, cp, sp, ctsp, stsp, xo, yo, zo;

void set_graphics_knobs(Real Theta, Real Phi)
{
  theta = Theta;
  phi = Phi;

  Real fac = 3.14159265358979323846 / 180.0;
  ct = cos(theta*fac);
  st = sin(theta*fac);
  cp = cos(phi*fac);
  sp = sin(phi*fac);
  ctsp = ct * sp;
  stsp = st * sp;
}

#endif

#if (BL_SPACEDIM == 1)

void real_coords(const IntVect& iv, const IntVect& type, const IntVect& ratio,
		 Real& x)
{
  if (type[0] == BOX_NODE) {
    x = iv[0];
  }
  else if (type[0] == BOX_CELL) {
    x = iv[0] + 0.5;
  }
  else {
    cerr << "real_coords():  invalid type" << endl;
    exit(1);
  }

  x /= ratio[0];
}

void real_coords(const Box& b, const IntVect& ratio, Real& x0, Real& x1)
{
  assert(b.ok());

  real_coords(b.smallEnd(), b.type(), ratio, x0);
  real_coords(b.bigEnd(),   b.type(), ratio, x1);
}

#elif (BL_SPACEDIM == 2)

void real_coords(const IntVect& iv, const IntVect& type, const IntVect& ratio,
		 Real& x, Real& y)
{
  if (type[0] == BOX_NODE) {
    x = iv[0];
  }
  else if (type[0] == BOX_CELL) {
    x = iv[0] + 0.5;
  }
  else {
    cerr << "real_coords():  invalid type" << endl;
    exit(1);
  }

  if (type[1] == BOX_NODE) {
    y = iv[1];
  }
  else if (type[1] == BOX_CELL) {
    y = iv[1] + 0.5;
  }
  else {
    cerr << "real_coords():  invalid type" << endl;
    exit(1);
  }
  x /= ratio[0];
  y /= ratio[1];
}

void real_coords(const Box& b, const IntVect& ratio,
		 Real& x0, Real& x1, Real& y0, Real& y1)
{
  assert(b.ok());

  real_coords(b.smallEnd(), b.type(), ratio, x0, y0);
  real_coords(b.bigEnd(),   b.type(), ratio, x1, y1);
}

#elif (BL_SPACEDIM == 3)

void real_coords(const IntVect& iv, const IntVect& type, const IntVect& ratio,
		 Real& x, Real& y, Real& z)
{
  if (type[0] == BOX_NODE) {
    x = iv[0];
  }
  else if (type[0] == BOX_CELL) {
    x = iv[0] + 0.5;
  }
  else {
    cerr << "real_coords():  invalid type" << endl;
    exit(1);
  }

  if (type[1] == BOX_NODE) {
    y = iv[1];
  }
  else if (type[1] == BOX_CELL) {
    y = iv[1] + 0.5;
  }
  else {
    cerr << "real_coords():  invalid type" << endl;
    exit(1);
  }

  if (type[2] == BOX_NODE) {
    z = iv[2];
  }
  else if (type[2] == BOX_CELL) {
    z = iv[2] + 0.5;
  }
  else {
    cerr << "real_coords():  invalid type" << endl;
    exit(1);
  }
  x /= ratio[0];
  y /= ratio[1];
  z /= ratio[2];
}

void real_coords(const Box& b, const IntVect& ratio,
		 Real& x0, Real& x1, Real& y0, Real& y1, Real& z0, Real& z1)
{
  assert(b.ok());

  real_coords(b.smallEnd(), b.type(), ratio, x0, y0, z0);
  real_coords(b.bigEnd(),   b.type(), ratio, x1, y1, z1);
}

#endif

void fit(Box b, const IntVect& ratio, Real margin)
{
  b.convert(nodevect);
#if (BL_SPACEDIM == 2)
  Real x0, x1, y0, y1;
  real_coords(b, ratio, x0, x1, y0, y1);
  fitbox(x0, x1, y0, y1, margin);
#else
  Real x0, x1, y0, y1, z0, z1;
  real_coords(b, ratio, x0, x1, y0, y1, z0, z1);
  Real diag = 0.5 * sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
  fitbox(-diag, diag, -diag, diag, margin);
  xo = 0.5 * (x0 + x1);
  yo = 0.5 * (y0 + y1);
  zo = 0.5 * (z0 + z1);
#endif
}

void plot(Box b, const IntVect& ratio)
{
  b.convert(nodevect);
#if (BL_SPACEDIM == 2)
  Real x0, x1, y0, y1;
  real_coords(b, ratio, x0, x1, y0, y1);
  mkbox(x0, x1, y0, y1);
#else
  Real x0, x1, y0, y1, z0, z1;
  real_coords(b, ratio, x0, x1, y0, y1, z0, z1);
  Real x, y;
  x = (x0 - xo) * ct   - (y0 - yo) * st;
  y = (x0 - xo) * stsp + (y0 - yo) * ctsp + (z0 - zo) * cp;
  pm(x,y);
  x = (x0 - xo) * ct   - (y0 - yo) * st;
  y = (x0 - xo) * stsp + (y0 - yo) * ctsp + (z1 - zo) * cp;
  pc(x,y);
  x = (x0 - xo) * ct   - (y1 - yo) * st;
  y = (x0 - xo) * stsp + (y1 - yo) * ctsp + (z1 - zo) * cp;
  pc(x,y);
  x = (x0 - xo) * ct   - (y1 - yo) * st;
  y = (x0 - xo) * stsp + (y1 - yo) * ctsp + (z0 - zo) * cp;
  pc(x,y);
  x = (x0 - xo) * ct   - (y0 - yo) * st;
  y = (x0 - xo) * stsp + (y0 - yo) * ctsp + (z0 - zo) * cp;
  pc(x,y);

  x = (x1 - xo) * ct   - (y0 - yo) * st;
  y = (x1 - xo) * stsp + (y0 - yo) * ctsp + (z0 - zo) * cp;
  pc(x,y);
  x = (x1 - xo) * ct   - (y0 - yo) * st;
  y = (x1 - xo) * stsp + (y0 - yo) * ctsp + (z1 - zo) * cp;
  pc(x,y);
  x = (x1 - xo) * ct   - (y1 - yo) * st;
  y = (x1 - xo) * stsp + (y1 - yo) * ctsp + (z1 - zo) * cp;
  pc(x,y);
  x = (x1 - xo) * ct   - (y1 - yo) * st;
  y = (x1 - xo) * stsp + (y1 - yo) * ctsp + (z0 - zo) * cp;
  pc(x,y);
  x = (x1 - xo) * ct   - (y0 - yo) * st;
  y = (x1 - xo) * stsp + (y0 - yo) * ctsp + (z0 - zo) * cp;
  pc(x,y);

  x = (x0 - xo) * ct   - (y0 - yo) * st;
  y = (x0 - xo) * stsp + (y0 - yo) * ctsp + (z1 - zo) * cp;
  pm(x,y);
  x = (x1 - xo) * ct   - (y0 - yo) * st;
  y = (x1 - xo) * stsp + (y0 - yo) * ctsp + (z1 - zo) * cp;
  pc(x,y);

  x = (x0 - xo) * ct   - (y1 - yo) * st;
  y = (x0 - xo) * stsp + (y1 - yo) * ctsp + (z1 - zo) * cp;
  pm(x,y);
  x = (x1 - xo) * ct   - (y1 - yo) * st;
  y = (x1 - xo) * stsp + (y1 - yo) * ctsp + (z1 - zo) * cp;
  pc(x,y);

  x = (x0 - xo) * ct   - (y1 - yo) * st;
  y = (x0 - xo) * stsp + (y1 - yo) * ctsp + (z0 - zo) * cp;
  pm(x,y);
  x = (x1 - xo) * ct   - (y1 - yo) * st;
  y = (x1 - xo) * stsp + (y1 - yo) * ctsp + (z0 - zo) * cp;
  pc(x,y);
#endif
}

void plot(const BoxArray& mesh, const IntVect& ratio)
{
  gflush(0);
  for (int i = 0; i < mesh.length(); i++)
    plot(mesh[i], ratio);
  gflush(1);
}

void plot(const Array<BoxArray>& mesh, const Array<IntVect>& gen_ratio)
{
  IntVect ratio = unitvect;
  plot(mesh[0], ratio);
  for (int i = 1; i < mesh.length(); i++) {
    ratio *= gen_ratio[i-1];
    plot(mesh[i], ratio);
  }
}

// Splits Box in side1 into two pieces, one which is completely
// outside side2 and one which is at least partly covered by side2.
// If the cut dimension is cell-based, the piece covered by side2 is
// expanded by 1 cell into the other piece---this overlap prevents
// gaps in contour lines crossing between the two pieces.
// On call, side1 and side2 should have the same type.
// On return, side1 and side2 contain the two fragments.  If
// no cut point was found side1 and side2 are unchanged.
// Returns 2 if a cut was found, 1 otherwise.
static int mask_box_chop(Box& side1, Box& side2)
{
  int retval = 1;
  for (int i = BL_SPACEDIM - 1; i >= 0; i--) {
    if (side1.type(i) == BOX_CELL) {
      if (side1.smallEnd(i) <= side2.bigEnd(i) &&
	  side2.bigEnd(i) < side1.bigEnd(i) - 1) {
	side2 = side1.chop(i,side2.bigEnd(i) + 1);
	side1.growHi(i,1);
	retval = 2;
	break;
      }
      if (side1.smallEnd(i) + 1 < side2.smallEnd(i) &&
	  side2.smallEnd(i) <= side1.bigEnd(i)) {
	side2 = side1.chop(i,side2.smallEnd(i));
	side2.growLo(i,1);
	retval = 2;
	break;
      }
    }
    else {
      if (side1.smallEnd(i) < side2.bigEnd(i) &&
	  side2.bigEnd(i) < side1.bigEnd(i)) {
	side2 = side1.chop(i,side2.bigEnd(i));
	retval = 2;
	break;
      }
      if (side1.smallEnd(i) < side2.smallEnd(i) &&
	  side2.smallEnd(i) < side1.bigEnd(i)) {
	side2 = side1.chop(i,side2.smallEnd(i));
	retval = 2;
	break;
      }
    }
  }
  return retval;
}

static int best_match(const BoxArray& a, const Box& region, int& igrid)
{
  int overlap = 0;
  for (int i = 0; i < a.length(); i++) {
    Box b = a[i];
    b.convert(region.type());
    if (region.intersects(b)) {
      int overlap1 = (region & b).numPts();
      if (overlap1 > overlap) {
	igrid = i;
	overlap = overlap1;
      }
    }
  }
  return (overlap > 0) ? (overlap == region.numPts() ? 1 : 2) : 0;
}

// Checks the Box in side1 against the mask for (cell-based) overlap.
// There are three possible cases:
// side1 is obscured-----function returns 1, side1 unchanged, side2 trashed.
// side1 is unobscured---function returns 0, side1 unchanged, side2 trashed.
// side1 is partially obscured---function returns 2, side1 and side2
// contain smaller boxes appropriate for recursive descent.
static int mask_box(Box& side1, Box& side2,
		    const BoxArray& mask, const IntVect& mrat)
{
  int igrid;
  side2 = side1;
  side2.refine(mrat).convert(cellvect);
  int retval = best_match(mask, side2, igrid);
  if (retval == 2) {
    side2 = mask[igrid];
    side2.coarsen(mrat).convert(side1.type());
    retval = mask_box_chop(side1, side2);
  }
  return retval;
}

void contour(const Fab& r, Box subbox, const IntVect& ratio,
	     Real value, int cflag, const Box& domain,
	     void (*set_contour_color)(int,int),
	     const BoxArray& mask, const IntVect& mrat)
{
  if (!r.box().sameType(subbox)) {
    subbox.convert(r.box().type());
  }

  if (cflag) {
    if (domain.ok()) {
      fit(domain, ratio);
    }
    newpage();
    plot(subbox, ratio);
  }

  if (mask.ready()) {
    Box side2;
    int retval = mask_box(subbox, side2, mask, mrat);
    if (retval == 1)
      return;
    if (retval == 2) {
      contour(r, subbox, ratio, value, 0, null_Box,
	      set_contour_color, mask, mrat);
      contour(r, side2,  ratio, value, 0, null_Box,
	      set_contour_color, mask, mrat);
      return;
    }
  }

  int ifail;
  int red, green, blue;
  getcolor(&red, &green, &blue);
  set_contour_color(1, 1);
  gflush(0);
#if (BL_SPACEDIM == 2)
  Real x0, x1, y0, y1;
  real_coords(subbox, ratio, x0, x1, y0, y1);
  FCONT2D(r.dataPtr(),
	  r.box().smallEnd(0), r.box().bigEnd(0),
	  r.box().smallEnd(1), r.box().bigEnd(1),
	  subbox.smallEnd(0), subbox.bigEnd(0),
	  subbox.smallEnd(1), subbox.bigEnd(1),
	  x0, x1, y0, y1,
	  1, value, ifail);
#else
  Real xoff = xo - 0.5 * (subbox.type(0) == BOX_CELL) / ratio[0];
  Real yoff = yo - 0.5 * (subbox.type(1) == BOX_CELL) / ratio[1];
  Real zoff = zo - 0.5 * (subbox.type(2) == BOX_CELL) / ratio[2];
  FCONT3D(r.dataPtr(),
	  r.box().smallEnd(0), r.box().bigEnd(0),
	  r.box().smallEnd(1), r.box().bigEnd(1),
	  r.box().smallEnd(2), r.box().bigEnd(2),
	  subbox.smallEnd(0), subbox.bigEnd(0),
	  subbox.smallEnd(1), subbox.bigEnd(1),
	  subbox.smallEnd(2), subbox.bigEnd(2),
	  ratio[0], ratio[1], ratio[2], xoff, yoff, zoff,
	  theta, phi, value, ifail);
#endif
  color(red, green, blue);
  gflush(1);
}

void contour(MultiFab& r, const IntVect& ratio, Real value,
	     int cflag, const Box& domain,
	     void (*set_contour_color)(int,int),
	     const BoxArray& mask, const IntVect& mrat)
{
  if (cflag) {
    if (domain.ok()) {
      fit(domain, ratio);
    }
    newpage();
    plot(r.boxArray(), ratio);
  }

  for (int i = 0; i < r.length(); i++) {
    contour(r[i], r.box(i), ratio, value, 0, null_Box,
	    set_contour_color, mask, mrat);
  }
}

void contour(PArray<MultiFab>& r, const Array<IntVect>& gen_ratio, Real value,
	     int cflag, const Box& domain0,
	     void (*set_contour_color)(int,int))
{
  int i;
  if (cflag) {
    if (domain0.ok()) {
      fit(domain0, unitvect);
    }
    newpage();
    for (i = 0; i < r.length(); i++)
      plot(r[i].boxArray());
  }

  IntVect ratio = unitvect;
  for (i = 0; i < r.length() - 1; i++) {
    contour(r[i], ratio, value, 0, null_Box,
	    set_contour_color, r[i+1].boxArray(), gen_ratio[i]);
    ratio *= gen_ratio[i];
  }
  i = r.length() - 1;
  contour(r[i], ratio, value, 0, null_Box, set_contour_color);
}

void contour(const Fab& r, Box subbox, const IntVect& ratio,
	     int nval, Real value[], int cflag, const Box& domain,
	     void (*set_contour_color)(int,int),
	     const BoxArray& mask, const IntVect& mrat)
{
  if (!r.box().sameType(subbox)) {
    subbox.convert(r.box().type());
  }

  if (cflag) {
    if (domain.ok()) {
      fit(domain, ratio);
    }
    newpage();
    plot(subbox, ratio);
  }

  if (mask.ready()) {
    Box side2;
    int retval = mask_box(subbox, side2, mask, mrat);
    if (retval == 1)
      return;
    if (retval == 2) {
      contour(r, subbox, ratio, nval, value, 0, null_Box,
	      set_contour_color, mask, mrat);
      contour(r, side2,  ratio, nval, value, 0, null_Box,
	      set_contour_color, mask, mrat);
      return;
    }
  }

  int red, green, blue;
  getcolor(&red, &green, &blue);
  for (int i = 0; i < nval; i++) {
    set_contour_color(i, nval);
    contour(r, subbox, ratio, value[i], 0);
  }
  color(red, green, blue);
}

void contour(MultiFab& r, const IntVect& ratio, int nval, Real value[],
	     int cflag, const Box& domain,
	     void (*set_contour_color)(int,int),
	     const BoxArray& mask, const IntVect& mrat)
{
  if (cflag) {
    if (domain.ok()) {
      fit(domain, ratio);
    }
    newpage();
    plot(r.boxArray(), ratio);
  }

  for (int i = 0; i < r.length(); i++) {
    contour(r[i], r.box(i), ratio, nval, value, 0, null_Box,
	    set_contour_color, mask, mrat);
  }
}

void contour(PArray<MultiFab>& r, const Array<IntVect>& gen_ratio,
	     int nval, Real value[],
	     int cflag, const Box& domain0,
	     void (*set_contour_color)(int,int))
{
  int i;
  IntVect ratio;
  if (cflag) {
    if (domain0.ok()) {
      fit(domain0, unitvect);
    }
    newpage();
    ratio = unitvect;
    for (i = 0; i < r.length(); i++) {
      if (i > 0)
	ratio *= gen_ratio[i-1];
      plot(r[i].boxArray(), ratio);
    }
  }

  ratio = unitvect;
  for (i = 0; i < r.length() - 1; i++) {
    contour(r[i], ratio, nval, value, 0, null_Box,
	    set_contour_color, r[i+1].boxArray(), gen_ratio[i]);
    ratio *= gen_ratio[i];
  }
  i = r.length() - 1;
  contour(r[i], ratio, nval, value, 0, null_Box, set_contour_color);
}

void contour(const Fab& r, Box subbox, const IntVect& ratio,
	     int nval, int flag,
	     Real a, Real b, Real p,
	     int cflag, const Box& domain,
	     void (*set_contour_color)(int,int),
	     const BoxArray& mask, const IntVect& mrat)
{
  if (!r.box().sameType(subbox)) {
    subbox.convert(r.box().type());
  }

  Real *value = new Real[nval];
  if (flag == 0 || flag == 1) {
    if (a == 0.0)
      a = r.norm(subbox, 0) * (nval-1) / nval;
    Real mid = 0.5 * (nval-1);
    for (int i = 0; i < nval; i++)
      value[i] = a * (i - mid) / mid;
    if (flag == 1 && (nval & 1))
      value[nval/2] = value[nval/2+1];
  }
  else {
    if (a == b) {
	a = r.min(subbox);
	b = r.max(subbox);
	for (int i = 0; i < nval; i++)
	  value[i] = a + (b - a) * pow((i+0.5)/nval, p);
    }
    else {
      for (int i = 0; i < nval; i++)
	value[i] = a + (b - a) * pow((Real)i/(nval-1), p);
    }
  }

  contour(r, subbox, ratio, nval, value, cflag, domain,
	  set_contour_color, mask, mrat);

  delete [] value;
}

void contour(MultiFab& r, const IntVect& ratio, int nval, int flag,
	     Real a, Real b, Real p,
	     int cflag, const Box& domain,
	     void (*set_contour_color)(int,int),
	     const BoxArray& mask, const IntVect& mrat)
{
  Real *value = new Real[nval];
  if (flag == 0 || flag == 1) {
    if (a == 0.0)
      a = mfnorm(r) * (nval-1) / nval;
    Real mid = 0.5 * (nval-1);
    for (int i = 0; i < nval; i++)
      value[i] = a * (i - mid) / mid;
    if (flag == 1 && (nval & 1))
      value[nval/2] = value[nval/2+1];
  }
  else {
    if (a == b) {
	a = r.min(0);
	b = r.max(0);
	for (int i = 0; i < nval; i++)
	  value[i] = a + (b - a) * pow((i+0.5)/nval, p);
    }
    else {
      for (int i = 0; i < nval; i++)
	value[i] = a + (b - a) * pow((Real)i/(nval-1), p);
    }
  }

  contour(r, ratio, nval, value, cflag, domain,
	  set_contour_color, mask, mrat);

  delete [] value;
}

void contour(PArray<MultiFab>& r, const Array<IntVect>& gen_ratio,
	     int nval, int flag,
	     Real a, Real b, Real p,
	     int cflag, const Box& domain0,
	     void (*set_contour_color)(int,int))
{
  Real *value = new Real[nval];
  if (flag == 0 || flag == 1) {
    if (a == 0.0)
      a = pmfnorm(r) * (nval-1) / nval;
    Real mid = 0.5 * (nval-1);
    for (int i = 0; i < nval; i++)
      value[i] = a * (i - mid) / mid;
    if (flag == 1 && (nval & 1))
      value[nval/2] = value[nval/2+1];
  }
  else {
    if (a == b) {
	a = pmfmin(r);
	b = pmfmax(r);
	for (int i = 0; i < nval; i++)
	  value[i] = a + (b - a) * pow((i+0.5)/nval, p);
    }
    else {
      for (int i = 0; i < nval; i++)
	value[i] = a + (b - a) * pow((Real)i/(nval-1), p);
    }
  }

  contour(r, gen_ratio, nval, value, cflag, domain0, set_contour_color);

  delete [] value;
}

#if (BL_SPACEDIM == 3)

void fit_slice(Box b, int idim, const IntVect& ratio, Real margin)
{
  b.convert(nodevect);
  Real x0, x1, y0, y1, z0, z1;
  real_coords(b, ratio, x0, x1, y0, y1, z0, z1);
  if (idim == 0)
    fitbox(y0, y1, z0, z1, margin);
  else if (idim == 1)
    fitbox(x0, x1, z0, z1, margin);
  else
    fitbox(x0, x1, y0, y1, margin);
}

void plot_slice(Box b, int idim, const IntVect& ratio)
{
  b.convert(nodevect);
  Real x0, x1, y0, y1, z0, z1;
  real_coords(b, ratio, x0, x1, y0, y1, z0, z1);
  if (idim == 0)
    mkbox(y0, y1, z0, z1);
  else if (idim == 1)
    mkbox(x0, x1, z0, z1);
  else
    mkbox(x0, x1, y0, y1);
}

void slice(const Fab& r, Box subbox, const IntVect& ratio,
	   int idim, int islice,
	   Real value, int cflag, const Box& domain,
	   void (*set_contour_color)(int,int))
{
  if (!r.box().sameType(subbox)) {
    subbox.convert(r.box().type());
  }

  if (cflag) {
    if (domain.ok()) {
      fit_slice(domain, idim, ratio);
    }
    newpage();
    plot_slice(subbox, idim, ratio);
  }

  if (islice < subbox.smallEnd(idim) || islice > subbox.bigEnd(idim))
    return;

  int ifail;
  int red, green, blue;
  getcolor(&red, &green, &blue);
  set_contour_color(1, 1);
  gflush(0);

  Real x0, x1, y0, y1, z0, z1;
  real_coords(subbox, ratio, x0, x1, y0, y1, z0, z1);
  if (idim == 0) {
    Real* sl = new Real[subbox.length(1) * subbox.length(2)];
    FGETSL0(r.dataPtr(),
            r.box().smallEnd(0), r.box().bigEnd(0),
            r.box().smallEnd(1), r.box().bigEnd(1),
            r.box().smallEnd(2), r.box().bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            sl, islice);
    FCONT2D(sl,
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            y0, y1, z0, z1,
            1, value, ifail);
    delete sl;
  }
  else if (idim == 1) {
    Real* sl = new Real[subbox.length(0) * subbox.length(2)];
    FGETSL1(r.dataPtr(),
            r.box().smallEnd(0), r.box().bigEnd(0),
            r.box().smallEnd(1), r.box().bigEnd(1),
            r.box().smallEnd(2), r.box().bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            sl, islice);
    FCONT2D(sl,
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(2), subbox.bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(2), subbox.bigEnd(2),
            x0, x1, z0, z1,
            1, value, ifail);
    delete sl;
  }
  else {
    Real* sl = new Real[subbox.length(0) * subbox.length(1)];
    FGETSL2(r.dataPtr(),
            r.box().smallEnd(0), r.box().bigEnd(0),
            r.box().smallEnd(1), r.box().bigEnd(1),
            r.box().smallEnd(2), r.box().bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            sl, islice);
    FCONT2D(sl,
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            x0, x1, y0, y1,
            1, value, ifail);
    delete sl;
  }

  color(red, green, blue);
  gflush(1);
}

void slice(const Fab& r, Box subbox, const IntVect& ratio,
           int idim, int islice,
	   int nval, Real value[], int cflag, const Box& domain,
	   void (*set_contour_color)(int,int))
{
  if (!r.box().sameType(subbox)) {
    subbox.convert(r.box().type());
  }

  if (cflag) {
    if (domain.ok()) {
      fit_slice(domain, idim, ratio);
    }
    newpage();
    plot_slice(subbox, idim, ratio);
  }

  if (islice < subbox.smallEnd(idim) || islice > subbox.bigEnd(idim))
    return;

  int ifail;
  int red, green, blue;
  getcolor(&red, &green, &blue);
  gflush(0);

  Real x0, x1, y0, y1, z0, z1;
  real_coords(subbox, ratio, x0, x1, y0, y1, z0, z1);
  if (idim == 0) {
    Real* sl = new Real[subbox.length(1) * subbox.length(2)];
    FGETSL0(r.dataPtr(),
            r.box().smallEnd(0), r.box().bigEnd(0),
            r.box().smallEnd(1), r.box().bigEnd(1),
            r.box().smallEnd(2), r.box().bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            sl, islice);
    for (int i = 0; i < nval; i++) {
      set_contour_color(i, nval);
      FCONT2D(sl,
              subbox.smallEnd(1), subbox.bigEnd(1),
              subbox.smallEnd(2), subbox.bigEnd(2),
              subbox.smallEnd(1), subbox.bigEnd(1),
              subbox.smallEnd(2), subbox.bigEnd(2),
              y0, y1, z0, z1,
              1, value[i], ifail);
    }
    delete sl;
  }
  else if (idim == 1) {
    Real* sl = new Real[subbox.length(0) * subbox.length(2)];
    FGETSL1(r.dataPtr(),
            r.box().smallEnd(0), r.box().bigEnd(0),
            r.box().smallEnd(1), r.box().bigEnd(1),
            r.box().smallEnd(2), r.box().bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            sl, islice);
    for (int i = 0; i < nval; i++) {
      set_contour_color(i, nval);
      FCONT2D(sl,
              subbox.smallEnd(0), subbox.bigEnd(0),
              subbox.smallEnd(2), subbox.bigEnd(2),
              subbox.smallEnd(0), subbox.bigEnd(0),
              subbox.smallEnd(2), subbox.bigEnd(2),
              x0, x1, z0, z1,
              1, value[i], ifail);
    }
    delete sl;
  }
  else {
    Real* sl = new Real[subbox.length(0) * subbox.length(1)];
    FGETSL2(r.dataPtr(),
            r.box().smallEnd(0), r.box().bigEnd(0),
            r.box().smallEnd(1), r.box().bigEnd(1),
            r.box().smallEnd(2), r.box().bigEnd(2),
            subbox.smallEnd(0), subbox.bigEnd(0),
            subbox.smallEnd(1), subbox.bigEnd(1),
            subbox.smallEnd(2), subbox.bigEnd(2),
            sl, islice);
    for (int i = 0; i < nval; i++) {
      set_contour_color(i, nval);
      FCONT2D(sl,
              subbox.smallEnd(0), subbox.bigEnd(0),
              subbox.smallEnd(1), subbox.bigEnd(1),
              subbox.smallEnd(0), subbox.bigEnd(0),
              subbox.smallEnd(1), subbox.bigEnd(1),
              x0, x1, y0, y1,
              1, value[i], ifail);
    }
    delete sl;
  }

  color(red, green, blue);
  gflush(1);
}

void slice(const Fab& r, Box subbox, const IntVect& ratio,
           int idim, int islice,
	   int nval, int flag,
	   Real a, Real b, Real p,
	   int cflag, const Box& domain,
	   void (*set_contour_color)(int,int))
{
  if (!r.box().sameType(subbox)) {
    subbox.convert(r.box().type());
  }

  Real *value = new Real[nval];
  if (flag == 0 || flag == 1) {
    if (a == 0.0)
      a = r.norm(subbox, 0) * (nval-1) / nval;
    Real mid = 0.5 * (nval-1);
    for (int i = 0; i < nval; i++)
      value[i] = a * (i - mid) / mid;
    if (flag == 1 && (nval & 1))
      value[nval/2] = value[nval/2+1];
  }
  else {
    if (a == b) {
	a = r.min();
	b = r.max();
	for (int i = 0; i < nval; i++)
	  value[i] = a + (b - a) * pow((i+0.5)/nval, p);
    }
    else {
      for (int i = 0; i < nval; i++)
	value[i] = a + (b - a) * pow((Real)i/(nval-1), p);
    }
  }

  slice(r, subbox, ratio, idim, islice, nval, value,
	cflag, domain, set_contour_color);

  delete [] value;
}

#endif
