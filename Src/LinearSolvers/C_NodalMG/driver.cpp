
#include "amr_real.H"
#include "amr_graph.H"

void malloc_info();

Real fct(const Intvect& it, const Intvect& t, int sig, int)
{
#if (SPACEDIM == 2)
  Real x, y;
  real_coords(it, t, sig, x, y);
  //return x * log(1+y);
  //return x;
  return 1 + 0.000001 * x;
#else
  Real x, y, z;
  real_coords(it, t, sig, x, y, z);
  return x * log(1+y) * log(1+z);
  //return 1 + 0.000001 * x;
  //return 1;
#endif
}

void graphtest(const amr_mesh& m);

main(int argc, char **argv)
{
  gopen(1);
  black();

  amr_mesh m;

  fstream grid;
  if (argc > 1)
    grid.open(argv[1], ios::in);
  else {
    cout << "usage:  driver <gridfile>" << endl;
    exit(1);
  }
  grid >> m;
  grid.close();

  graphtest(m);

  gclose();
}

void graphtest(const amr_mesh& m)
{
#if (SPACEDIM == 2)
  //amr_real r(m, Intvect(BOX_CELL, BOX_CELL));
  //amr_real r(m, Intvect(BOX_NODE, BOX_NODE), 3);
  amr_real r(m, Intvect(BOX_NODE, BOX_NODE));
  int ncont = 11;
#else
  //amr_real r(m, Intvect(BOX_CELL, BOX_CELL, BOX_CELL));
  amr_real r(m, Intvect(BOX_NODE, BOX_NODE, BOX_NODE));
  set_graphics_knobs(210.0, 10.0);
  int ncont = 4;
#endif

  {
    cout << "Initializing..." << flush;
    r.initialize(fct);
    cout << "done." << endl;

    fit(m[1].domain());
    //contour(r[1][0], m[1][0], 1, ncont, 2, 0.0, 0.0, 1.0, 0);
    //cin.get();
    //contour(r[1][1], m[1][1], 1, ncont, 2, 0.0, 0.0, 1.0, 1, m[1].domain());
    //fit(m[0].domain(), m[0].sig());
    //plot(m);
    //cin.get();
    {
      level_interface li(m[1], reflection_boundary);
      //r[0][0].assign(0.0);
      r[1].restrict_patch(r[0][0], r[0].mesh().sig(),
			  //bilinear_restrictor,
			  bilinear_restrictor_coarse,
			  li, reflection_boundary);
      //r[1].restrict_patch(r[0][0], r[0].mesh().sig(), injection_restrictor);
      //r[1].restrict_patch(r[0][0], r[0].mesh().sig(), bilinear_restrictor);
      //r[1].restrict_patch(r[0][0], r[0].mesh().sig());
      contour(r[0], ncont, 2);
      plot(m[1]);
      cin.get();
    }
    contour(r[1], ncont, 2, 0.0, 0.0, 1.0, 1);
    cin.get();
    contour(r, ncont, 2, 0.0, 0.0, 1.0);
    cin.get();

    //r.interpolate_patch(r[1][1], r[1][1].box(), m[1].sig(), 0);
    r[0].interpolate_patch(r[1][1], r[1][1].box(), m[1].sig());
    contour(r, ncont, 2, 0.0, 0.0, 1.0);
    cin.get();

    r.mult(r);
    contour(r, ncont, 2, 0.0, 0.0, 1.0);
    cin.get();
    //malloc_info();
    r.detach();
    //malloc_info();
  }
  //malloc_info();
}
