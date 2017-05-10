
#include <WarpX.H>
#include <AMReX_BoxIterator.H>
#include <array>
#include <algorithm>

using namespace amrex;

void
WarpX::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const auto& gm = Geom(lev);
    const Box& domain = gm.Domain();
    const IntVect& sz = domain.size();
    IntVect ctr = sz / 2;

    // for testing, let's tag a sphere with a radius of
    const Real R = 0.25*std::min({D_DECL(sz[0],sz[1],sz[2])});

    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        auto& fab = tags[mfi];
        const Box& bx = fab.box();
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const IntVect& cell = bi();
            Real rsq = 0.0;
            for (int idim=0; idim<BL_SPACEDIM; ++idim) {
                Real d = cell[idim]-ctr[idim];
                rsq += d*d;
            }
            if (rsq < R*R) {
                fab(cell) = TagBox::SET;
            }
        }
    }
}

