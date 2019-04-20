#ifndef AMREX_PETSC_H_
#define AMREX_PETSC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BndryData.H>

// Do NOT include any PETSc headers here.  They should be included in
// .cpp files before AMReX_PETSc.H is included.

namespace amrex {

class PETScABecLap       
{
public:

    
    PETScABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);

    ~PETScABecLap ();

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

    KSP solver = nullptr;
    Mat A = nullptr;
    Vec b = nullptr;
    Vec x = nullptr;

    LayoutData<PetscInt> ncells_grid;
    LayoutData<Vector<PetscInt> > cell_id_vec;
    FabArray<BaseFab<PetscInt> > cell_id;

    MultiFab const* m_eb_b_coeffs = nullptr;

    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

std::unique_ptr<PETScABecLap> makePetsc (const BoxArray& grids, const DistributionMapping& damp,
                                  const Geometry& geom, MPI_Comm comm_);

}

#endif
