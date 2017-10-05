
#include <AMReX_MultiFab.H>

#include <fort.H>

using namespace amrex;

void init_prob (const Array<Geometry>& geom, Array<MultiFab>& alpha, Array<MultiFab>& beta,
                Array<MultiFab>& rhs, Array<MultiFab>& exact)
{
    const Real a = 1.e-3;
    const Real b = 1.0;
    const Real sigma = 10.0;
    const Real w     = 0.05;

    const int nlevels = geom.size();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        const Real* problo = geom[ilev].ProbLo();
        const Real* probhi = geom[ilev].ProbHi();
        const Real* dx     = geom[ilev].CellSize();

        for (MFIter mfi(alpha[ilev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            fort_set_coef(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_ANYD(exact[ilev][mfi]),
                          BL_TO_FORTRAN_ANYD(alpha[ilev][mfi]),
                          BL_TO_FORTRAN_ANYD(beta[ilev][mfi]),
                          BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                          dx, problo, probhi, &a, &b, &sigma, &w);
        }
    }
}
