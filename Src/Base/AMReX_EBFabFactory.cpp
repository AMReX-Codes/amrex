
#include <AMReX_FabFactory.H>
#include <AMReX_FArrayBox.H>

namespace amrex
{

template <>
FArrayBox*
DefaultFabFactory<FArrayBox>::create (const Box& box, int ncomps,
                                      const FabInfo& info, int box_index) const
{
    return new FArrayBox(box, ncomps, info.alloc, info.shared);
}

}
