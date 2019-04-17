#include <AMReX_EBCellFlag.H>
#include <AMReX_EBCellFlag_F.H>

namespace amrex {

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

std::ostream&
operator<< (std::ostream& os, const EBCellFlag& flag)
{
    std::ios_base::fmtflags old_fmt = os.flags();
    os << std::hex << flag.getValue() << ":" << std::dec;

    if (flag.isRegular()) {
        os << "R";
    } else if (flag.isSingleValued()) {
        os << "S";
    } else if (flag.isCovered()) {
        os << "C";
    } else {
        os << "M";
    }

#if (AMREX_SPACEDIM == 3)
    for (int k = -1; k <= 1; ++k) {
#endif
        for (int j = -1; j <= 1; ++j) {
            for (int i = -1; i <= 1; ++i) {
                os << static_cast<int>(flag.isConnected({AMREX_D_DECL(i,j,k)}));
            }
        }
#if (AMREX_SPACEDIM == 3)
    }
#endif

    os.flags(old_fmt);

    return os;
}

}
