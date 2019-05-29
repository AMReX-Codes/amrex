
#include <WarpX.H>
#include <AMReX_BoxIterator.H>
#include <array>
#include <algorithm>

using namespace amrex;

void
WarpX::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const Real* problo = Geom(lev).ProbLo();
    const Real* dx = Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        auto& fab = tags[mfi];
        const Box& bx = fab.box();
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const IntVect& cell = bi();
            RealVect pos {AMREX_D_DECL((cell[0]+0.5)*dx[0]+problo[0],
                                       (cell[1]+0.5)*dx[1]+problo[1],
                                       (cell[2]+0.5)*dx[2]+problo[2])};
            if (pos > fine_tag_lo && pos < fine_tag_hi) {
                fab(cell) = TagBox::SET;
            }
        }
    }
}

