#include <AMReX_HypreNodeLap.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <cmath>
#include <numeric>
#include <limits>
#include <type_traits>

namespace amrex {

HypreNodeLap::HypreNodeLap (const BoxArray& grids, const DistributionMapping& dmap,
                            const Geometry& geom_, MPI_Comm comm_)
{
}

HypreNodeLap::~HypreNodeLap ()
{
}

std::unique_ptr<HypreNodeLap>
makeHypreNodeLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom, MPI_Comm comm_)
{
    return std::unique_ptr<HypreNodeLap>(new HypreNodeLap(grids,dmap,geom,comm_));
}

}

