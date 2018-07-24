#ifndef AMREX_HYPREABECLAP3_H_
#define AMREX_HYPREABECLAP3_H_


#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LayoutData.H>

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include <algorithm>

namespace amrex
{

class HypreABecLap3
{
public:

    HypreABecLap3 (const BoxArray& grids,
                   const DistributionMapping& dmap,
                   const Geometry& geom_,
                   MPI_Comm comm_);
    ~HypreABecLap3 ();

    void setScalars (Real sa, Real sb);
    void setACoeffs (const MultiFab& alpha);
    void setBCoeffs (const Array<const MultiFab*,AMREX_SPACEDIM>& beta);
    void setVerbose (int _verbose);
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                int max_iter, const BndryData& bndry, int max_bndry_order);

private :

    static constexpr HYPRE_Int regular_stencil_size = 2*AMREX_SPACEDIM + 1;
    static constexpr HYPRE_Int eb_stencil_size = AMREX_D_TERM(3, *3, *3);
    
    MPI_Comm comm = MPI_COMM_NULL;
    Geometry geom;

    int verbose = 0;

    MultiFab acoefs;
    Array<MultiFab,AMREX_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;

    HYPRE_IJMatrix A = NULL;
    HYPRE_IJVector b = NULL;
    HYPRE_IJVector x = NULL;
    HYPRE_Solver solver = NULL;

    FabFactory<FArrayBox> const* m_factory = nullptr;
    BndryData const* m_bndry = nullptr;
    int m_maxorder = -1;
    
    LayoutData<HYPRE_Int> ncells_grid;
    LayoutData<Vector<HYPRE_Int> > cell_id_vec;
    FabArray<BaseFab<HYPRE_Int> > cell_id;
    
    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

}

#endif
