
#include "hg_elliptic.H"

DEFINE_ERROR_FUNC(holy_grail_amr_elliptic)

void holy_grail_amr_elliptic::general(const amr_real& Dest,
				      const amr_real& Source,
				      const amr_real& Sigma,
				      Real H[], Real reltol,
				      int Lev_min, int Lev_max,
				      Real abstol)
{
  if (Lev_min < 0)
    Lev_min = lev_min_max;
  if (Lev_max < 0)
    Lev_max = Lev_min;

#ifndef CONSTANT
  if (Sigma.null())
    error("general---null Sigma input to variable-coefficient version");
#endif

  alloc(Dest, Source, Sigma, H, Lev_min, Lev_max);
  solve(reltol, abstol, 2, 2);
  clear();
}

void holy_grail_amr_elliptic::poisson(Fab& FDest,
				      Fab& FSource,
				      Real H[], Real reltol,
				      Real abstol)
{
  if (FDest.box().type() != nodevect)
    error("FDest not NODE-based");
  if (FSource.box().type() != nodevect)
    error("FSource not NODE-based");

  Box tmp = grow(ml_mesh[0].boxn(0), 1);
  if (tmp != FDest.box())
    error("FDest box not compatible with problem domain + 1 ghost cell");
  if (tmp != FSource.box())
    error("FSource box not compatible with problem domain + 1 ghost cell");

  level_real LDest(ml_mesh[0], nodevect, 1, 0);
  LDest.set_grid(0, FDest);
  amr_real Dest(ml_mesh, nodevect, 1, 0);
  Dest.set_level(0, LDest);

  level_real LSource(ml_mesh[0], nodevect, 1, 0);
  LSource.set_grid(0, FSource);
  amr_real Source(ml_mesh, nodevect, 1, 0);
  Source.set_level(0, LSource);

  alloc(Dest, Source, null_amr_real, H, 0, 0);
  solve(reltol, abstol, 2, 2);
  clear();
}

