#include <AMReX_GpuLaunch.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_LayoutData.H>

namespace amrex {
namespace Gpu {

#ifdef AMREX_USE_GPU
void getGridSize (FabArrayBase const& fa, int ngrow, LayoutData<GridSize>& gs, int& ntotblocks)
{
    gs = LayoutData<GridSize>(fa.boxArray(),fa.DistributionMap());
    ntotblocks = 0;
    for (MFIter mfi(gs); mfi.isValid(); ++mfi) {
        const auto& bx = amrex::grow(mfi.validbox(),ngrow);
        auto ec = ExecutionConfig(bx);
        gs[mfi].numBlocks = ec.numBlocks.x;
        gs[mfi].numThreads = ec.numThreads.x;
        gs[mfi].globalBlockId = ntotblocks;
        ntotblocks += ec.numBlocks.x;
    }
}
#endif

}}
