#include <AMReX_EBCellFlag.H>

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

        FabType t = FabType::undefined;
        int nregular=0, nsingle=0, nmulti=0, ncovered=0;
        int ncells = bx.numPts();
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const auto& cellflag = (*this)(bi());
            if (cellflag.isRegular()) {
                ++nregular;
            } else if (cellflag.isSingleValued()) {
                ++nsingle;
            } else if (cellflag.isMultiValued()) {
                ++nmulti;
            } else {
                ++ncovered;
            }
        }

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
