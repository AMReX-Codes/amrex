
#include "MyTest.H"
#include "MyTest_F.H"

using namespace amrex;

void
MyTest::initProbPoisson ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
#ifdef _OPENMP
#oragma omp parallel
#endif
        for (MFIter mfi(rhs[ilev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            actual_init_poisson(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                                BL_TO_FORTRAN_ANYD(exact_solution[ilev][mfi]),
                                geom[ilev].ProbLo(), geom[ilev].ProbHi(),
                                geom[ilev].CellSize());
        }

        solution[ilev].setVal(0.0);
    }
}

void
MyTest::initProbABecLaplacian ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
#ifdef _OPENMP
#oragma omp parallel
#endif
        for (MFIter mfi(rhs[ilev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = mfi.growntilebox(1);
            actual_init_abeclap(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_BOX(gbx),
                                BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                                BL_TO_FORTRAN_ANYD(exact_solution[ilev][mfi]),
                                BL_TO_FORTRAN_ANYD(acoef[ilev][mfi]),
                                BL_TO_FORTRAN_ANYD(bcoef[ilev][mfi]),
                                ascalar, bscalar,
                                geom[ilev].ProbLo(), geom[ilev].ProbHi(),
                                geom[ilev].CellSize());
        }

        solution[ilev].setVal(0.0);
    }
}
