#include "Poisson.H"
#include "Poisson_F.H"

using namespace amrex;

void InitData (MultiFab& State) {

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(State.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    int ng = State.nGrow();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        const auto& sfab = State[mfi];
        const auto& flag = flags[mfi];

        if (flag.getType(bx) != FabType::covered) {
            init_data(BL_TO_FORTRAN_BOX(bx),
                      BL_TO_FORTRAN_ANYD(State[mfi]));
        }
    }
}
