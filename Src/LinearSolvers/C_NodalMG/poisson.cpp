
#include "amr_graph.H"
#include "hg_elliptic.H"
#include "amr_gravity.H"

void malloc_info();

#ifndef UNICOS
#  define USE_GRAPHICS
#endif

void poisson_test(const amr_mesh& m);

main(int argc, char **argv)
{
#ifdef USE_GRAPHICS
  gopen(5);
  black();
#endif

  amr_mesh m;

  fstream grid;
  if (argc > 1)
    grid.open(argv[1], ios::in);
  else {
    cout << "usage:  poisson <gridfile>" << endl;
    exit(1);
  }
  grid >> m;
  grid.close();

  //cout << m << endl;
  //cin.get();

  //malloc_info();
  poisson_test(m);
  //malloc_info();

#ifdef USE_GRAPHICS
  gclose();
#endif
}

void node_init(const amr_real& src, const amr_mesh& m, amr_boundary bdy,
               const amr_real& work)
{
  src.assign(0.0);

  Real radius = 0.04;
  Real sig2 = radius * radius / 4.0;
  Real pi = 3.14159265358979323846;
  Real fac2d = 0.5 / (pi * sig2);
  Real fac3d = pow(fac2d, 1.5);
  int ilev;
  for (ilev = 0 ; ilev < m.nlevels() ; ilev++) {
    Real h = 1.0 / m[ilev].domain().length(0);
    for (int igrid = 0; igrid < m[ilev].ngrids() ; igrid++) {
      for (int i = m[ilev][igrid].smallEnd(0);
	   i <= m[ilev][igrid].bigEnd(0) + 1; i++) {
	Real x = i * h;
	for (int j = m[ilev][igrid].smallEnd(1);
	     j <= m[ilev][igrid].bigEnd(1) + 1; j++) {
          Real y = j * h;
#if (SPACEDIM == 2)
          src[ilev][igrid](Iv(i,j)) = fac2d *
	    (exp(-((x-0.4)*(x-0.4) + (y-0.7)*(y-0.7)) / (2.0 * sig2)) -
	     exp(-((x-0.6)*(x-0.6) + (y-0.3)*(y-0.3)) / (2.0 * sig2)));
#else
	  for (int k = m[ilev][igrid].smallEnd(2);
	       k <= m[ilev][igrid].bigEnd(2) + 1; k++) {
	    Real z = k * h;
	    src[ilev][igrid](Iv(i,j,k)) = fac3d *
	      (exp(-((x-0.4)*(x-0.4) + (y-0.7)*(y-0.7) + (z-0.6)*(z-0.6)) /
		   (2.0 * sig2)) -
	       exp(-((x-0.6)*(x-0.6) + (y-0.3)*(y-0.3) + (z-0.4)*(z-0.4)) /
		   (2.0 * sig2)));
	  }
#endif
        }
      }
    }
  }

  for (ilev = m.nlevels() - 1; ilev > 0; ilev--) {
    level_interface interface(m[ilev], bdy);
    src[ilev].restrict_level(src[ilev-1], 0, bilinear_restrictor_coarse,
			     interface, bdy);
  }
  work[0].assign(1.0);
  Real adjust = src[0].inner_product(work[0]) / m[0].domain().volume();
  cout << "Solvability adjustment: " << adjust << endl;
  src.minus(adjust);
}

void cell_init(const amr_real& src, const amr_mesh& m)
{
  src.assign(0.0);

  Real radius = 0.04;
  Real sig2 = radius * radius / 4.0;
  Real pi = 3.14159265358979323846;
  Real fac2d = 0.5 / (pi * sig2);
  Real fac3d = pow(fac2d, 1.5);
  int ilev;
  for (ilev = 0 ; ilev < m.nlevels() ; ilev++) {
    Real h = 1.0 / m[ilev].domain().length(0);
    for (int igrid = 0; igrid < m[ilev].ngrids() ; igrid++) {
      for (int i = m[ilev][igrid].smallEnd(0);
	   i <= m[ilev][igrid].bigEnd(0); i++) {
	Real x = (i + 0.5) * h;
	for (int j = m[ilev][igrid].smallEnd(1);
	     j <= m[ilev][igrid].bigEnd(1); j++) {
          Real y = (j + 0.5) * h;
#if (SPACEDIM == 2)
          src[ilev][igrid](Iv(i,j)) = fac2d *
	    (exp(-((x-0.4)*(x-0.4) + (y-0.7)*(y-0.7)) / (2.0 * sig2)) -
	     exp(-((x-0.6)*(x-0.6) + (y-0.3)*(y-0.3)) / (2.0 * sig2)));
#else
	  for (int k = m[ilev][igrid].smallEnd(2);
	       k <= m[ilev][igrid].bigEnd(2); k++) {
	    Real z = (k + 0.5) * h;
#  if (SPACEDIM == 4)
	    src[ilev][igrid](Iv(i,j,k)) = y;
#  else
/*
	    src[ilev][igrid](Iv(i,j,k)) = fac3d *
	      (exp(-((x-0.4)*(x-0.4) + (y-0.7)*(y-0.7) + (z-0.6)*(z-0.6)) /
		   (2.0 * sig2)) -
	       exp(-((x-0.6)*(x-0.6) + (y-0.3)*(y-0.3) + (z-0.4)*(z-0.4)) /
		   (2.0 * sig2)));
*/
	    src[ilev][igrid](Iv(i,j,k)) = fac3d *
	      (exp(-((x-0.4)*(x-0.4) + (y-0.7)*(y-0.7) + (z-0.6)*(z-0.6)) /
		   (2.0 * sig2)) -
	 2.0 * exp(-((x-0.6)*(x-0.6) + (y-0.01)*(y-0.01) + (z-0.4)*(z-0.4)) /
		   (2.0 * sig2)));
#  endif
	  }
#endif
        }
      }
    }
  }
}

void poisson_test(const amr_mesh& m)
{
  int ilev, i;

  Real h[SPACEDIM];
  for (i = 0; i < SPACEDIM; i++)
    h[i] = 1;

  RegType bc[SPACEDIM][2];

#if (SPACEDIM == 2)
  int ncont = 11;
  bc[0][0] = periodic;
  bc[0][1] = periodic;
  //bc[0][0] = refWall;
  //bc[0][1] = refWall;
  bc[1][0] = periodic;
  bc[1][1] = periodic;
  //bc[1][0] = refWall;
  //bc[1][1] = refWall;
#else
  set_graphics_knobs(30.0, 10.0);
  int ncont = 10;
  for (i = 0; i < SPACEDIM; i++) {
    bc[i][0] = periodic;
    bc[i][1] = periodic;
  }
#endif

  amr_real dest(m, nodevect, 1);
  amr_real cell_src(m, cellvect, 1);

  cell_init(cell_src, m);
  dest.assign(0.0);

#ifdef USE_GRAPHICS
  fit(m.domain());
  plot(m);
  //contour(cell_src, ncont);
  //cin.get();
#endif

  int sum = 0;
  cout << "Cells by level: ";
  for (ilev = 0; ilev < m.nlevels(); ilev++) {
    cout << " " << m[ilev].numPts();
    sum += m[ilev].numPts();
  }
  cout << "\nTotal cells:  " << sum << endl;

  clock_t t0, t1, t2;

#ifdef UNICOS
  int pcode = 1;
  Real tol = 1.e-6;
#else
  int pcode = 2;
  //Real tol = 1.e-6;
  //Real tol = 2.e-12;
  Real tol = 2.e-10;
  //Real tol = 5.e-9;
#endif
  inviscid_fluid_boundary afb(bc);
  t0 = clock();
  amr_gravity_module cell_solver(m, 0, m.nlevels() - 1, m.nlevels() - 1,
                                 afb, pcode);
  t1 = clock();
  cell_solver.poisson(dest, cell_src, h, tol, 0, m.nlevels() - 1);
  t2 = clock();
  cell_init(cell_src, m);
  cell_solver.poisson(dest, cell_src, h, tol, 0, m.nlevels() - 1);
  cell_solver.poisson(dest, cell_src, h, tol, 0, m.nlevels() - 1);
  cout << "Init time was     " << t1 - t0 << endl;
  cout << "Solution time was " << t2 - t1 << endl;
  cout << "Speed was " << double(t2 - t1) / sum << endl;
#ifdef USE_GRAPHICS
  contour(dest, ncont);
  cin.get();
  slice(dest[1][0], m[1].boxn(0), m[1].sig(), 0, 6, 15,
        0, 0.0, 0.0, 1.0, 1, m[1].domain());
  cin.get();
  slice(dest[1][0], m[1].boxn(0), m[1].sig(), 1, 9, 15,
        0, 0.0, 0.0, 1.0, 1, m[1].domain());
  cin.get();
  slice(dest[1][0], m[1].boxn(0), m[1].sig(), 2, 18, 15,
        0, 0.0, 0.0, 1.0, 1, m[1].domain());
  cin.get();
#endif

/*
  amr_real node_src(m, nodevect, 1);
  node_init(node_src, m, afb.scalar(), dest);
  dest.assign(0.0);

#ifdef USE_GRAPHICS
  //plot(m);
  contour(node_src, ncont);
  cin.get();
#endif

  t0 = clock();
  holy_grail_amr_elliptic node_solver(m, 0, m.nlevels() - 1, m.nlevels() - 1,
                                      afb, pcode);
  t1 = clock();
  node_solver.poisson(dest, node_src, h, tol, 0, m.nlevels() - 1);
  t2 = clock();
  cout << "Init time was     " << t1 - t0 << endl;
  cout << "Solution time was " << t2 - t1 << endl;
  cout << "Speed was " << double(t2 - t1) / sum << endl;
#ifdef USE_GRAPHICS
  contour(dest, ncont);
  cin.get();
#endif
*/
/*
  // Single grid solver
  init(src, m);
  dest.assign(0.0);

  t0 = clock();
  holy_grail_amr_elliptic grid_solver(m[0].domain(), afb, pcode);

  t1 = clock();
  grid_solver.poisson(dest[0][0].fab(), src[0][0].fab(), h, tol);
  t2 = clock();
  cout << "Init time was     " << t1 - t0 << endl;
  cout << "Solution time was " << t2 - t1 << endl;
  cout << "Speed was " << double(t2 - t1) / sum << endl;
#ifdef USE_GRAPHICS
  contour(dest[0][0], m[0].boxn(0), m[0].sig(), ncont);
  cin.get();
#endif
*/
}
