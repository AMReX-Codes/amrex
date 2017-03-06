#ifndef _dirichletTest_F_H_
#define _dirichletTest_F_H_

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_IntVect.H>

using namespace amrex;

struct Stencil
{
  Real data[9];
  Real& operator[](int i) {return data[i];}
};

extern "C"
{
  void get_bndry_grad_stencil(Stencil* stencil, Real* stencilB, IntVect* base,
			      const RealVect* normals, const RealVect* centroids, const IntVect* ivs,
			      const int* num, const Real* dx);

  void apply_bndry_grad_stencil(Real* g, Real* b,
                                const Real* u, ARLIM_P(ulo), ARLIM_P(uhi),
                                const Stencil* stencil, const Real* stencilB, const IntVect *baseiv,
                                const int * num);
  
  void init_phi(Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
		const int* reglo, const int* reghi,
		const Real* plo, const Real* dx);
}

#endif // _dirichletTest_F_H_
