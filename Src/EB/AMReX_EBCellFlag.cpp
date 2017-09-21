#include <AMReX_EBCellFlag.H>
#include <AMReX_EBCellFlag_F.H>

namespace amrex {

constexpr std::array<std::array<std::array<int,3>,3>,3> EBCellFlag::pos_ngbr;

FabType
EBCellFlagFab::getType (const Box& bx_in) const
{
    FabType thistype = getType();

    if (thistype == FabType::regular)
    {
        return FabType::regular;
    }
    else if (thistype == FabType::covered)
    {
        return FabType::covered;
    }
    else
    {
        const Box& bx = amrex::enclosedCells(bx_in);
        BL_ASSERT(this->box().contains(bx));

        int nregular, nsingle, nmulti, ncovered;
        amrex_ebcellflag_count(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(*this),
                               &nregular, &nsingle, &nmulti, &ncovered);

        int ncells = bx.numPts();

        FabType t = FabType::undefined;
        if (nregular == ncells) {
            t = FabType::regular;
        } else if (ncovered == ncells) {
            t = FabType::covered;
        } else if (nmulti > 0) {
            t = FabType::multivalued;
        } else {
            t = FabType::singlevalued;
        }

        return t;
    }
}

}
