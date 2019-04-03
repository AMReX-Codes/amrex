#ifndef AMREX_HYPRE_NODE_LAP_H_
#define AMREX_HYPRE_NODE_LAP_H_

#include <memory>
#include <algorithm>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BndryData.H>

#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"

namespace amrex {

class HypreNodeLap
{
public:
    HypreNodeLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);
    ~HypreNodeLap ();

private:

    HYPRE_IJMatrix A = NULL;
    HYPRE_IJVector b = NULL;
    HYPRE_IJVector x = NULL;
    HYPRE_Solver solver = NULL;

    LayoutData<HYPRE_Int> ncells_grid;
    LayoutData<Vector<HYPRE_Int> > cell_id_vec;
    FabArray<BaseFab<HYPRE_Int> > cell_id;
};

std::unique_ptr<HypreNodeLap> makeHypreNodeLap (const BoxArray& grids, const DistributionMapping& dmap,
                                                const Geometry& geom, MPI_Comm comm_);

}

#endif
