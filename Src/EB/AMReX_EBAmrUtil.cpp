
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBAmrUtil_F.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

void
TagCutCells (TagBoxArray& tags, const MultiFab& state)
{
    const char   tagval = TagBox::SET;
    const char clearval = TagBox::CLEAR;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto& sfab = dynamic_cast<EBFArrayBox const&>(state[mfi]);
        const auto& flag = sfab.getEBCellFlagFab();

        const FabType typ = flag.getType(bx);
        if (typ != FabType::regular && typ != FabType::covered)
        {
            TagBox& tagfab = tags[mfi];
            amrex_tag_cutcells(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(tagfab),
                               BL_TO_FORTRAN_ANYD(flag),
                               tagval, clearval);
        }
    }
}

}
