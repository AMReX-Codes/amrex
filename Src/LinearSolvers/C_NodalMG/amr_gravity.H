#ifndef _AMR_GRAVITY_H_
#define _AMR_GRAVITY_H_

#ifndef _HG_MULTI_H_
#  include "hg_multi.H"
#endif

class amr_gravity_module : public holy_grail_amr_multigrid {
  DECLARE_ERROR_FUNC(amr_gravity_module)
 protected:
  void right_hand_side(const amr_real& Source);
  void interface_average(const amr_real& Source, int lev);

 public:

  // No, I'm not pulling your leg.  Lev_min and Lev_max, supplied with
  // project, tell the range of levels over which to project.  The
  // arguments to the constructor give the permissible range for Lev_min
  // and Lev_max.  Lev_min can be anything from Lev_min_min to Lev_min_max;
  // Lev_max can then be anything from Lev_min to Lev_max_max.

  amr_gravity_module(const amr_mesh& Mesh,
		     int Lev_min_min, int Lev_min_max, int Lev_max_max,
		     amr_fluid_boundary Boundary, int Pcode = 2)
    : holy_grail_amr_multigrid(Mesh, Lev_min_min, Lev_min_max, Lev_max_max,
			       Boundary, Pcode) { }

  const amr_mesh& mesh() {
    return ml_mesh;
  }

  // Let start_norm = max[norm(source), norm(first residual)].
  // Solver will attempt to reduce the Holy Grail residual to
  // less than max[reltol * start_norm, abstol].  At least one multilevel
  // cycle will be executed down to lev_min.  The last cycle will omit
  // any actions on or below lev_min.

  // H is an array of SPACEDIM Reals containing the fine mesh spacing
  // in each dimension.  "fine" means lev_max.  If a positive value is
  // not supplied for lev_min, lev_min_min will be used.  If a
  // positive value is not supplied for lev_max, lev_max_max will be used.
  // Note: this is different than defaults in the projection code!

  void poisson(const amr_real& Dest, const amr_real& Src,
	       Real H[], Real reltol,
	       int Lev_min = -1, int Lev_max = -1, Real abstol = 0.0);

  void gradient(Real* grad_ptr, const Box& grad_box,
		const Real* potential_ptr, const Box& potential_box,
                const Box& region, const Real* dx);
};

#endif
