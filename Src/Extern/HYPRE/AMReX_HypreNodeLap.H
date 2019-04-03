#ifndef AMREX_HYPRE_NODE_LAP_H_
#define AMREX_HYPRE_NODE_LAP_H_

#include <memory>
#include <algorithm>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BndryData.H>
#include <AMReX_MultiFabUtil.H>

#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"

namespace amrex {

class HypreNodeLap
{
public:
    HypreNodeLap (const BoxArray& grids_, const DistributionMapping& dmap_,
                  const Geometry& geom_, const FabFactory<FArrayBox>& factory_,
                  const iMultiFab& owner_mask_, const iMultiFab& dirichlet_mask_,
                  MPI_Comm comm_);
    ~HypreNodeLap ();

private:

    BoxArray grids;
    DistributionMapping dmap;
    Geometry geom;
    FabFactory<FArrayBox> const* factory = nullptr;
    MPI_Comm comm = MPI_COMM_NULL;
    iMultiFab const* owner_mask = nullptr;
    iMultiFab const* dirichlet_mask = nullptr;

    HYPRE_IJMatrix A = NULL;
    HYPRE_IJVector b = NULL;
    HYPRE_IJVector x = NULL;
    HYPRE_Solver solver = NULL;

    LayoutData<HYPRE_Int> nnodes_grid;
    LayoutData<Vector<HYPRE_Int> > node_id_vec;
    FabArray<BaseFab<HYPRE_Int> > node_id;
};

}

#endif
