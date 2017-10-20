
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <fort.H>
#include <prob_par.H>

using namespace amrex;

namespace prob {

    amrex::Real a     = 1.e-3;
    amrex::Real b     = 1.0;
    amrex::Real sigma = 10.0;
    amrex::Real w     = 0.05;

}

using namespace prob;

void init_prob (const Vector<Geometry>& geom, Vector<MultiFab>& alpha, Vector<MultiFab>& beta,
                Vector<MultiFab>& rhs, Vector<MultiFab>& exact)
{
    {
        ParmParse pp("prob");
        pp.query("a"    , a);
        pp.query("b"    , b);
        pp.query("sigma", sigma);
        pp.query("w"    , w);
    }


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
                          dx, problo, probhi,
                          &a, &b, &sigma, &w);
        }
    }
}
