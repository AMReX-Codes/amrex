#ifndef AMREX_HYPREABECLAP3_H_
#define AMREX_HYPREABECLAP3_H_

#include <AMReX_Hypre.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>

#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"

#include <algorithm>

namespace amrex
{

class HypreABecLap3
    : public Hypre
{
public:

    HypreABecLap3 (const BoxArray& grids, const DistributionMapping& dmap,
                   const Geometry& geom_, MPI_Comm comm_);
    virtual ~HypreABecLap3 ();

    virtual void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                        int max_iter, const BndryData& bndry, int max_bndry_order) final;

#ifdef AMREX_USE_EB
    void setEBDirichlet (MultiFab const* beb) { m_eb_b_coeffs = beb; }
#endif

private :

    HYPRE_IJMatrix A = NULL;
    HYPRE_IJVector b = NULL;
    HYPRE_IJVector x = NULL;
    HYPRE_Solver solver = NULL;

    LayoutData<HYPRE_Int> ncells_grid;
    LayoutData<Gpu::ManagedDeviceVector<HYPRE_Int> > cell_id_vec;
    FabArray<BaseFab<HYPRE_Int> > cell_id;

    MultiFab const* m_eb_b_coeffs = nullptr;
    
    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

}

#endif
