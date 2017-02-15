#ifndef _stencilTest_F_H_
#define _stencilTest_F_H_

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

extern "C"
{
    void lapleb_MSD(Real* lphi,      ARLIM_P(lphlo), ARLIM_P(lphhi),
		    const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
		    const Real*  cd, ARLIM_P(cdlo),  ARLIM_P(cdhi),
		    const Real* fd0, ARLIM_P(f0lo),  ARLIM_P(f0hi),
#if BL_SPACEDIM>1
		    const Real* fd1, ARLIM_P(f1lo),  ARLIM_P(f1hi),
#if BL_SPACEDIM>2
		    const Real* fd2, ARLIM_P(f2lo),  ARLIM_P(f2hi),
#endif
#endif
		    const int* reglo, const int* reghi, const Real* dx);

    void lapl_MSD(Real* lphi,      ARLIM_P(lphlo), ARLIM_P(lphhi),
		  const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
		  const int* reglo, const int* reghi, const Real* dx);

    void init_phi(Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
		  const int* reglo, const int* reghi,
                  const Real* plo, const Real* dx);
}

#endif // _stencilTest_F_H_
