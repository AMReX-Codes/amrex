
#include <AMReX_HypreABecLap.H>
#include <AMReX_Habec_K.H>

#include <_hypre_struct_mv.h>

#include <string>
#include <algorithm>

namespace amrex {

HypreABecLap::HypreABecLap(const BoxArray& grids, const DistributionMapping& dmap,
                           const Geometry& geom_, MPI_Comm comm_)
    : Hypre(grids, dmap, geom_, comm_)
{
}

HypreABecLap::~HypreABecLap ()
{
    HYPRE_StructPFMGDestroy(solver);
    solver = NULL;
    HYPRE_StructMatrixDestroy(A);
    A = NULL;
//    HYPRE_StructVectorDestroy(b);
//    b = NULL;
//    HYPRE_StructVectorDestroy(x);
//    x = NULL;
    HYPRE_StructGridDestroy(grid);
    grid = NULL;
}

void
HypreABecLap::solve(MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
                    int maxiter, const BndryData& bndry, int max_bndry_order)
{
    if (solver == NULL || m_bndry != &bndry || m_maxorder != max_bndry_order)
    {
        m_bndry = &bndry;
        m_maxorder = max_bndry_order;
        m_factory = &(rhs.Factory());
        prepareSolver();
    }
    else
    {
        m_factory = &(rhs.Factory());
    }

    // do this repeatedly to avoid memory leak
    HYPRE_StructVectorCreate(comm, grid, &b);
    HYPRE_StructVectorCreate(comm, grid, &x);
    //
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);
    //
    loadVectors(soln, rhs);
    //
    HYPRE_StructVectorAssemble(x);
    HYPRE_StructVectorAssemble(b);

    HYPRE_StructPFMGSetMaxIter(solver, maxiter);
    HYPRE_StructPFMGSetTol(solver, reltol);
    if (abstol > 0.0)
    {
        Real bnorm;
        bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
                                      (hypre_StructVector *) b);
        bnorm = std::sqrt(bnorm);

        const BoxArray& grids = rhs.boxArray();
        Real volume = grids.d_numPts();
        Real reltol_new = (bnorm > 0.0) ? (abstol / bnorm * std::sqrt(volume)) : reltol;

        if (reltol_new > reltol) {
            HYPRE_StructPFMGSetTol(solver, reltol_new);
        }
    }

    HYPRE_StructPFMGSolve(solver, A, b, x);

    if (verbose >= 2)
    {
        HYPRE_Int num_iterations;
        Real res;
        HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
        HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &res);

        amrex::Print() << "\n" << num_iterations
                       << " Hypre PFMG Iterations, Relative Residual "
                       << res << std::endl;
    }

    getSolution(soln);

    // do this repeatedly to avoid memory leak
    HYPRE_StructVectorDestroy(b);
    b = NULL;
    HYPRE_StructVectorDestroy(x);
    x = NULL;
}

void
HypreABecLap::getSolution (MultiFab& a_soln)
{
    MultiFab* soln = &a_soln;
    MultiFab tmp;
    if (a_soln.nGrowVect() != 0) {
        tmp.define(a_soln.boxArray(), a_soln.DistributionMap(), 1, 0);
        soln = &tmp;
    }

    for (MFIter mfi(*soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();
        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_StructVectorGetBoxValues(x, reglo.data(), reghi.data(), (*soln)[mfi].dataPtr());
    }

    if (a_soln.nGrowVect() != 0) {
        MultiFab::Copy(a_soln, tmp, 0, 0, 1, 0);
    }
}

void
HypreABecLap::prepareSolver ()
{
    BL_PROFILE("HypreABecLap::prepareSolver()");

    HYPRE_StructGridCreate(comm, AMREX_SPACEDIM, &grid);

    Array<HYPRE_Int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (geom.isPeriodic(i)) {
            is_periodic[i] = geom.period(i);
            AMREX_ASSERT(Hypre::ispow2(is_periodic[i]));
            AMREX_ASSERT(geom.Domain().smallEnd(i) == 0);
        }
    }
    if (geom.isAnyPeriodic()) {
        HYPRE_StructGridSetPeriodic(grid, is_periodic.data());
    }

    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto lo = Hypre::loV(bx);
        auto hi = Hypre::hiV(bx);
        HYPRE_StructGridSetExtents(grid, lo.data(), hi.data());
    }

    HYPRE_StructGridAssemble(grid);

#if (AMREX_SPACEDIM == 2)
    HYPRE_Int offsets[regular_stencil_size][2] = {{ 0,  0},    // 0
                                                  {-1,  0},    // 1
                                                  { 1,  0},    // 2
                                                  { 0, -1},    // 3
                                                  { 0,  1}};   // 4
#elif (AMREX_SPACEDIM == 3)
    HYPRE_Int offsets[regular_stencil_size][3] = {{ 0,  0,  0},   // 0
                                                  {-1,  0,  0},   // 1
                                                  { 1,  0,  0},   // 2
                                                  { 0, -1,  0},   // 3
                                                  { 0,  1,  0},   // 4
                                                  { 0,  0, -1},   // 5
                                                  { 0,  0,  1}};  // 6
#endif

    HYPRE_StructStencil stencil;
    HYPRE_StructStencilCreate(AMREX_SPACEDIM, regular_stencil_size, &stencil);

    for (int i = 0; i < regular_stencil_size; ++i) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

    HYPRE_StructMatrixCreate(comm, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    HYPRE_StructStencilDestroy(stencil);

    HYPRE_StructVectorCreate(comm, grid, &b);
    HYPRE_StructVectorCreate(comm, grid, &x);

    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);

    // A.SetValues() & A.assemble()
    Array<HYPRE_Int,regular_stencil_size> stencil_indices;
    std::iota(stencil_indices.begin(), stencil_indices.end(), 0);
    const auto dx = geom.CellSizeArray();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    BaseFab<GpuArray<Real,regular_stencil_size> > rfab;
    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();
        rfab.resize(reg);

        Array4<Real const> const& afab = acoefs.const_array(mfi);
        GpuArray<Array4<Real const>, AMREX_SPACEDIM> bfabs {
            AMREX_D_DECL(bcoefs[0].const_array(mfi),
                         bcoefs[1].const_array(mfi),
                         bcoefs[2].const_array(mfi))};
        Array4<Real> const& diaginvfab = diaginv.array(mfi);
        GpuArray<int,AMREX_SPACEDIM*2> bctype;
        GpuArray<Real,AMREX_SPACEDIM*2> bcl;
        GpuArray<Array4<int const>, AMREX_SPACEDIM*2> msk;
        for (OrientationIter oit; oit; oit++)
        {
            Orientation ori = oit();
            int cdir(ori);
            bctype[cdir] = m_bndry->bndryConds(mfi)[cdir][0];
            bcl[cdir] = m_bndry->bndryLocs(mfi)[cdir];
            msk[cdir] = m_bndry->bndryMasks(ori)[mfi].const_array();
        }

        Real sa = scalar_a;
        Real sb = scalar_b;
        const auto boxlo = amrex::lbound(reg);
        const auto boxhi = amrex::ubound(reg);

        amrex::fill(rfab,
        [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,regular_stencil_size>& sten,
                                   int i, int j, int k)
        {
            habec_mat(sten, i, j, k, boxlo, boxhi, sa, afab, sb, dx, bfabs,
                      bctype, bcl, bho, msk, diaginvfab);
        });

        Real* mat = (Real*) rfab.dataPtr();
        Gpu::streamSynchronize();

        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_StructMatrixSetBoxValues(A, reglo.data(), reghi.data(),
                                       regular_stencil_size, stencil_indices.data(),
                                       mat);
        Gpu::synchronize();
    }
    HYPRE_StructMatrixAssemble(A);

    // create solver
    HYPRE_StructPFMGCreate(comm, &solver);
    int logging = (verbose >= 2) ? 1 : 0;
    HYPRE_StructPFMGSetLogging(solver, logging);
    HYPRE_StructPFMGSetup(solver, A, b, x);

    HYPRE_StructVectorDestroy(b);
    b = NULL;
    HYPRE_StructVectorDestroy(x);
    x = NULL;
}


void
HypreABecLap::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("HypreABecLap::loadVectors()");

    soln.setVal(0.0);

    MultiFab rhs_diag(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
#ifdef AMREX_USE_OMP
#pragma omp paralle if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs_diag,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& reg = mfi.validbox();
        Array4<Real> const& rhs_diag_a = rhs_diag.array(mfi);
        Array4<Real const> const& rhs_a = rhs.const_array(mfi);
        Array4<Real const> const& diaginv_a = diaginv.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
        {
            rhs_diag_a(i,j,k) = rhs_a(i,j,k) * diaginv_a(i,j,k);
        });
    }

    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();
        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_StructVectorSetBoxValues(x, reglo.data(), reghi.data(), soln[mfi].dataPtr());
        HYPRE_StructVectorSetBoxValues(b, reglo.data(), reghi.data(), rhs_diag[mfi].dataPtr());
    }
}

}
