#ifndef AMREX_HYPREABECLAP3_H_
#define AMREX_HYPREABECLAP3_H_
#include <AMReX_Config.H>

#include <AMReX_Hypre.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_HypreIJIface.H>

#include <algorithm>

namespace amrex
{

class HypreABecLap3 final
    : public Hypre
{
public:

    HypreABecLap3 (const BoxArray& grids, const DistributionMapping& dmap,
                   const Geometry& geom_, MPI_Comm comm_,
                   const iMultiFab* overset_mask = nullptr);
    ~HypreABecLap3 () override;

    HypreABecLap3 (HypreABecLap3 const&) = delete;
    HypreABecLap3 (HypreABecLap3 &&) = delete;
    HypreABecLap3& operator= (HypreABecLap3 const&) = delete;
    HypreABecLap3& operator= (HypreABecLap3 &&) = delete;

    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                int max_iter, const BndryData& bndry, int max_bndry_order) override;

#ifdef AMREX_USE_EB
    void setEBDirichlet (MultiFab const* beb) { m_eb_b_coeffs = beb; }
#endif

private :
    std::unique_ptr<HypreIJIface> hypre_ij;

    // Non-owning references to hypre matrix, rhs, and solution data
    HYPRE_IJMatrix A = nullptr;
    HYPRE_IJVector b = nullptr;
    HYPRE_IJVector x = nullptr;

    LayoutData<HYPRE_Int> ncells_grid;
    FabArray<BaseFab<HYPRE_Int> > cell_id;
    FabArray<BaseFab<HYPRE_Int> > cell_id_vec;

    MultiFab const* m_eb_b_coeffs = nullptr;

    iMultiFab const* m_overset_mask = nullptr;

public: // for CUDA
    void prepareSolver ();
    void getSolution (MultiFab& soln);
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
};

}

#endif
