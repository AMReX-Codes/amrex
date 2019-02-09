
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

    auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
    auto const& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto& flag = flags[mfi];

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


void
TagVolfrac (TagBoxArray& tags, const MultiFab& volfrac, Real tol)
{
    BL_PROFILE("amrex::TagVolfrac()");

    const char clearval = TagBox::CLEAR;
    const char   tagval = TagBox::SET;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(volfrac, true); mfi.isValid(); ++mfi) {
        const Box &    tilebox = mfi.tilebox();
              TagBox & tagfab  = tags[mfi];

        //_______________________________________________________________________
        // Tag cells for refinement
        amrex_tag_volfrac ( BL_TO_FORTRAN_BOX(tilebox),
                            BL_TO_FORTRAN_ANYD(tagfab),
                            BL_TO_FORTRAN_ANYD(volfrac[mfi]),
                            tagval, clearval, tol );
    }
}

}
