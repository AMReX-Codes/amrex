#ifndef AMREX_PETSC_H_
#define AMREX_PETSC_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BndryData.H>

// Do NOT include any other PETSc headers here.  They should be included in
// .cpp files before AMReX_PETSc.H is included.
#include <petscsystypes.h>

namespace amrex {

class amrex_KSP;
class amrex_Mat;
class amrex_Vec;

class PETScABecLap
{
public:

    PETScABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);

    ~PETScABecLap ();

    PETScABecLap (PETScABecLap const&) = delete;
    PETScABecLap (PETScABecLap &&) noexcept = delete;
    PETScABecLap& operator= (PETScABecLap const&) = delete;
    PETScABecLap& operator= (PETScABecLap &&) noexcept = delete;

    void setScalars (Real sa, Real sb);
    void setACoeffs (const MultiFab& alpha);
    void setBCoeffs (const Array<const MultiFab*,BL_SPACEDIM>& beta);
    void setVerbose (int _verbose);
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                int max_iter, const BndryData& bndry, int max_bndry_order);

#ifdef AMREX_USE_EB
    void setEBDirichlet (MultiFab const* beb) { m_eb_b_coeffs = beb; }
#endif

private:

    static constexpr PetscInt regular_stencil_size = 2*AMREX_SPACEDIM + 1;
    static constexpr PetscInt eb_stencil_size = AMREX_D_TERM(3, *3, *3);

    MPI_Comm comm = MPI_COMM_NULL;
    Geometry geom;

    int verbose = 0;

    MultiFab acoefs;
    Array<MultiFab,AMREX_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;

    MultiFab diaginv;

    FabFactory<FArrayBox> const* m_factory = nullptr;
    BndryData const* m_bndry = nullptr;
    int m_maxorder = -1;

    std::unique_ptr<amrex_KSP> solver;
    std::unique_ptr<amrex_Mat> A;
    std::unique_ptr<amrex_Vec> b;
    std::unique_ptr<amrex_Vec> x;

    LayoutData<PetscInt> ncells_grid;
    FabArray<BaseFab<PetscInt> > cell_id;
    FabArray<BaseFab<PetscInt> > cell_id_vec;

    MultiFab const* m_eb_b_coeffs = nullptr;

public: // for cuda
    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

[[nodiscard]] std::unique_ptr<PETScABecLap>
makePetsc (const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, MPI_Comm comm_);

}

#endif
