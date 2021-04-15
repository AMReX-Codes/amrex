#include <AMReX_HypreABecLap3.H>
#include <AMReX_VisMF.H>

#include <AMReX_Habec_K.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <type_traits>

namespace amrex {

HypreABecLap3::HypreABecLap3 (const BoxArray& grids, const DistributionMapping& dmap,
                              const Geometry& geom_, MPI_Comm comm_,
                              const iMultiFab* overset_mask)
    : Hypre(grids, dmap, geom_, comm_),
      m_overset_mask(overset_mask)
{}

HypreABecLap3::~HypreABecLap3 ()
{}

void
HypreABecLap3::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                      int max_iter, const BndryData& bndry, int max_bndry_order)
{
    BL_PROFILE("HypreABecLap3::solve()");

    if (!hypre_ij || m_bndry != &bndry || m_maxorder != max_bndry_order)
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

    HYPRE_IJVectorInitialize(b);
    HYPRE_IJVectorInitialize(x);
    //
    loadVectors(soln, rhs);
    //
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorAssemble(b);

    hypre_ij->solve(rel_tol, abs_tol, max_iter);

    getSolution(soln);
}

void
HypreABecLap3::getSolution (MultiFab& a_soln)
{
    bool use_tmp_mf = (a_soln.nGrowVect() != 0);
    MultiFab* l_soln = &a_soln;
    MultiFab tmp;
    if (use_tmp_mf) {
        tmp.define(a_soln.boxArray(), a_soln.DistributionMap(), 1, 0);
        l_soln = &tmp;
    }

    for (MFIter mfi(*l_soln); mfi.isValid(); ++mfi)
    {
        const HYPRE_Int nrows = ncells_grid[mfi];
        if (nrows > 0) {
            HYPRE_IJVectorGetValues(x, nrows, cell_id_vec[mfi].dataPtr(), (*l_soln)[mfi].dataPtr());
        } else {
            (*l_soln)[mfi].setVal<RunOn::Device>(0.0);
        }
    }

    if (use_tmp_mf) {
        MultiFab::Copy(a_soln, tmp, 0, 0, 1, 0);
    }
}

void
HypreABecLap3::prepareSolver ()
{
    BL_PROFILE("HypreABecLap3::prepareSolver()");

    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myid);

    const BoxArray& ba = acoefs.boxArray();
    const DistributionMapping& dm = acoefs.DistributionMap();

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(HYPRE_Int) < sizeof(Long)) {
        Long ncells_grids = ba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncells_grids < static_cast<Long>(std::numeric_limits<HYPRE_Int>::max()),
                                         "You might need to configure Hypre with --enable-bigint");
    }
#endif

    static_assert(std::is_signed<HYPRE_Int>::value, "HYPRE_Int is assumed to be signed");

    // how many non-covered cells do we have?
    ncells_grid.define(ba,dm);
    cell_id.define(ba,dm,1,1);
    cell_id_vec.define(ba,dm,1,0);

#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_overset_mask == nullptr || ebfactory == nullptr,
                                     "Cannot have both EB and overset");
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (ebfactory) ? &(ebfactory->getVolFrac()) : nullptr;
    auto area = (ebfactory) ? ebfactory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (ebfactory) ? ebfactory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto barea = (ebfactory) ? &(ebfactory->getBndryArea()) : nullptr;
    auto bcent = (ebfactory) ? &(ebfactory->getBndryCent()) : nullptr;
#endif

    HYPRE_Int ncells_proc = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+:ncells_proc)
#endif
    for (MFIter mfi(cell_id); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const Box& gbx = amrex::grow(bx,1);
        Array4<HYPRE_Int> const& cid_arr = cell_id.array(mfi);
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered)
        {
            ncells_grid[mfi] = 0;
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(gbx, i, j, k,
            {
                cid_arr(i,j,k) = std::numeric_limits<HYPRE_Int>::lowest();
            });
        }
        else
#endif
        {
            Long npts = bx.numPts();
            ncells_grid[mfi] = npts;
            ncells_proc += npts;

            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(gbx, i, j, k,
            {
                if (bx.contains(i,j,k)) {
                    cid_arr(i,j,k) = bx.index(IntVect{AMREX_D_DECL(i,j,k)});
                } else {
                    cid_arr(i,j,k) = std::numeric_limits<HYPRE_Int>::lowest();
                }
            });
        }
    }

    Vector<HYPRE_Int> ncells_allprocs(num_procs);
    MPI_Allgather(&ncells_proc, sizeof(HYPRE_Int), MPI_CHAR,
                  ncells_allprocs.data(), sizeof(HYPRE_Int), MPI_CHAR,
                  comm);
    HYPRE_Int proc_begin = 0;
    for (int i = 0; i < myid; ++i) {
        proc_begin += ncells_allprocs[i];
    }

    LayoutData<HYPRE_Int> offset(ba,dm);
    HYPRE_Int proc_end = proc_begin;
    for (MFIter mfi(ncells_grid); mfi.isValid(); ++mfi)
    {
        offset[mfi] = proc_end;
        proc_end += ncells_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end == proc_begin+ncells_proc,
                                     "HypreABecLap3::prepareSolver: how did this happen?");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cell_id,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto os = offset[mfi];
        Array4<HYPRE_Int> const& cid_arr = cell_id.array(mfi);
        Array4<HYPRE_Int> const& cid_vec = cell_id_vec.array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(bx, i, j, k,
        {
            cid_arr(i,j,k) += os;
            cid_vec(i,j,k) = cid_arr(i,j,k);
        });
    }

    cell_id.FillBoundary(geom.periodicity());

    // Create and initialize A, b & x
    HYPRE_Int ilower = proc_begin;
    HYPRE_Int iupper = proc_end-1;

    hypre_ij = std::make_unique<HypreIJIface>(comm, ilower, iupper, verbose);
    hypre_ij->parse_inputs(options_namespace);

    // Obtain non-owning references to the matrix, rhs, and solution data
    A = hypre_ij->A();
    b = hypre_ij->b();
    x = hypre_ij->x();

    const auto dx = geom.CellSizeArray();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    BaseFab<HYPRE_Int> ncols_fab;

    BaseFab<Real> mat_aos_fab, mat_vec_fab;
    BaseFab<HYPRE_Int> cols_aos_fab, cols_vec_fab;

    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

#ifdef AMREX_USE_EB
        const auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
#else
        const auto fabtyp = FabType::regular;
#endif
        if (fabtyp != FabType::covered)
        {
            ncols_fab.resize(bx);

            int max_stencil_size;
            if  (fabtyp == FabType::regular) {
                max_stencil_size = 2*AMREX_SPACEDIM+1;
            } else {
                max_stencil_size = AMREX_D_TERM(3,*3,*3);
            }

            mat_aos_fab.resize(bx,max_stencil_size);
            cols_aos_fab.resize(bx,max_stencil_size);

            Array4<HYPRE_Int> const& ncols_a = ncols_fab.array();
            Array4<HYPRE_Int const> const& cid_a = cell_id.const_array(mfi);

            Array4<Real const> const& afab = acoefs.const_array(mfi);
            GpuArray<Array4<Real const>, AMREX_SPACEDIM> bfabs {
                AMREX_D_DECL(bcoefs[0].const_array(mfi),
                             bcoefs[1].const_array(mfi),
                             bcoefs[2].const_array(mfi))};
            Array4<Real> const& diaginvfab = diaginv.array(mfi);
            GpuArray<int,AMREX_SPACEDIM*2> bctype;
            GpuArray<Real,AMREX_SPACEDIM*2> bcl;
            for (OrientationIter oit; oit; oit++)
            {
                int cdir = oit();
                bctype[cdir] = m_bndry->bndryConds(mfi)[cdir][0];
                bcl[cdir] = m_bndry->bndryLocs(mfi)[cdir];
            }

            Real sa = scalar_a;
            Real sb = scalar_b;

            if (fabtyp == FabType::regular)
            {
                auto osmsk = (m_overset_mask) ? m_overset_mask->const_array(mfi)
                                              : Array4<int const>();
                constexpr int stencil_size = 2*AMREX_SPACEDIM+1;
                BaseFab<GpuArray<Real,stencil_size> > tmpmatfab
                    (bx, 1, (GpuArray<Real,stencil_size>*)mat_aos_fab.dataPtr());

                amrex::fill(tmpmatfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_ijmat(sten, ncols_a, diaginvfab, i, j, k, cid_a,
                                sa, afab, sb, dx, bfabs, bctype, bcl, bho, osmsk);
                });

                BaseFab<GpuArray<HYPRE_Int,stencil_size> > tmpcolfab
                    (bx, 1, (GpuArray<HYPRE_Int,stencil_size>*)cols_aos_fab.dataPtr());

                amrex::fill(tmpcolfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<HYPRE_Int,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_cols(sten, i, j, k, cid_a);
                });
            }
#ifdef AMREX_USE_EB
            else
            {
                auto const& flag_a = flags->const_array(mfi);
                auto const& vfrac_a = vfrac->const_array(mfi);
                AMREX_D_TERM(auto const& apx = area[0]->const_array(mfi);,
                             auto const& apy = area[1]->const_array(mfi);,
                             auto const& apz = area[2]->const_array(mfi);)
                AMREX_D_TERM(auto const& fcx = fcent[0]->const_array(mfi);,
                             auto const& fcy = fcent[1]->const_array(mfi);,
                             auto const& fcz = fcent[2]->const_array(mfi);)
                auto const& barea_a = barea->const_array(mfi);
                auto const& bcent_a = bcent->const_array(mfi);
                Array4<Real const> beb = (m_eb_b_coeffs) ? m_eb_b_coeffs->const_array(mfi)
                                                         : Array4<Real const>();

                constexpr int stencil_size = AMREX_D_TERM(3,*3,*3);
                BaseFab<GpuArray<Real,stencil_size> > tmpmatfab
                    (bx, 1, (GpuArray<Real,stencil_size>*)mat_aos_fab.dataPtr());

                amrex::fill(tmpmatfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_ijmat_eb(sten, ncols_a, diaginvfab, i, j, k, cid_a,
                                   sa, afab, sb, dx, bfabs, bctype, bcl, bho,
                                   flag_a, vfrac_a, AMREX_D_DECL(apx,apy,apz),
                                   AMREX_D_DECL(fcx,fcy,fcz),barea_a,bcent_a,beb);
                });

                BaseFab<GpuArray<HYPRE_Int,stencil_size> > tmpcolfab
                    (bx, 1, (GpuArray<HYPRE_Int,stencil_size>*)cols_aos_fab.dataPtr());

                amrex::fill(tmpcolfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<HYPRE_Int,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_cols_eb(sten, i, j, k, cid_a);
                });
            }
#endif

            const HYPRE_Int nrows = ncells_grid[mfi];
            AMREX_ASSERT(nrows == static_cast<HYPRE_Int>(bx.numPts()));
            HYPRE_Int* rows = cell_id_vec[mfi].dataPtr();
            HYPRE_Int* ncols = ncols_fab.dataPtr();

            // Remove invalid elements
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
            if (sizeof(HYPRE_Int) < sizeof(Long)) {
                Long ntot = static_cast<Long>(nrows)*max_stencil_size;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ntot <  static_cast<Long>(std::numeric_limits<HYPRE_Int>::max()),
                                                 "Integer overflow: please configure Hypre with --enable-bigint");
            }
#endif
            HYPRE_Int nelems = nrows * max_stencil_size;
            HYPRE_Int const* cols_in = cols_aos_fab.dataPtr();
            Real const* mat_in = mat_aos_fab.dataPtr();
            HYPRE_Int* cols;
            Real* mat;
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
                cols_vec_fab.resize(bx,max_stencil_size);
                mat_vec_fab.resize(bx,max_stencil_size);
                cols = cols_vec_fab.dataPtr();
                mat = mat_vec_fab.dataPtr();
                Scan::PrefixSum<HYPRE_Int>(nelems,
                    [=] AMREX_GPU_DEVICE (HYPRE_Int i) -> HYPRE_Int { return (mat_in[i] != Real(0.0)); },
                    [=] AMREX_GPU_DEVICE (HYPRE_Int i, HYPRE_Int const& x) {
                        if (mat_in[i] != Real(0.0)) {
                            cols[x] = cols_in[i];
                            mat[x] = mat_in[i];
                        }
                    }, Scan::Type::exclusive);
            } else
#endif
            {
                cols = const_cast<HYPRE_Int*>(cols_in);
                mat = const_cast<Real*>(mat_in);
                HYPRE_Int m = 0;
                for (HYPRE_Int ielem = 0; ielem < nelems; ++ielem) {
                    if (mat_in[ielem] != Real(0.0)) {
                        if (m != ielem) {
                            cols[m] = cols_in[ielem];
                            mat[m] = mat_in[ielem];
                        }
                        ++m;
                    }
                }
            }

            // For singular matrices set reference solution on one row
            if (hypre_ij->adjustSingularMatrix() && is_matrix_singular) {
                AMREX_ASSERT_WITH_MESSAGE(fabtyp == FabType::regular,
                                          "adjustSingularMatrix not supported for EB");
                AMREX_HOST_DEVICE_FOR_1D(1, m,
                {
                    amrex::ignore_unused(m);
                    if (rows[0] == 0) {
                        const int num_cols = ncols[0];
                        for (int ic = 0; ic < num_cols; ++ic) {
                            mat[ic] = (cols[ic] == rows[0]) ? mat[ic] : 0.0;
                        }
                    }
                });
            }

            Gpu::synchronize();
            HYPRE_IJMatrixSetValues(A,nrows,ncols,rows,cols,mat);
            Gpu::synchronize();
        }
    }
    HYPRE_IJMatrixAssemble(A);
}

void
HypreABecLap3::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("HypreABecLap3::loadVectors()");

#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

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
        auto osm = (m_overset_mask) ? m_overset_mask->const_array(mfi)
                                    : Array4<int const>();
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(reg) : FabType::regular;
        if (fabtyp == FabType::singlevalued) {
            auto const& flag = flags->const_array(mfi);
            if (osm) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
                {
                    rhs_diag_a(i,j,k) = (osm(i,j,k) == 0 || flag(i,j,k).isCovered()) ?
                        Real(0.0) : rhs_a(i,j,k) * diaginv_a(i,j,k);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
                {
                    rhs_diag_a(i,j,k) = (flag(i,j,k).isCovered()) ?
                        Real(0.0) : rhs_a(i,j,k) * diaginv_a(i,j,k);
                });
            }
        } else if (fabtyp == FabType::regular)
#endif
        {
            if (osm) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
                {
                    rhs_diag_a(i,j,k) = (osm(i,j,k) == 0) ?
                        Real(0.0) : rhs_a(i,j,k) * diaginv_a(i,j,k);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
                {
                    rhs_diag_a(i,j,k) = rhs_a(i,j,k) * diaginv_a(i,j,k);
                });
            }
        }
    }

    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const HYPRE_Int nrows = ncells_grid[mfi];
        if (nrows > 0)
        {
            // soln has been set to zero.
            HYPRE_IJVectorSetValues(x, nrows, cell_id_vec[mfi].dataPtr(), soln[mfi].dataPtr());

            if (hypre_ij->adjustSingularMatrix() && is_matrix_singular) {
                HYPRE_Int const* rows = cell_id_vec[mfi].dataPtr();
                Real* bp = rhs_diag[mfi].dataPtr();
                AMREX_HOST_DEVICE_FOR_1D(1, m,
                {
                    amrex::ignore_unused(m);
                    if (rows[0] == 0) {
                        bp[0] = Real(0.0);
                    }
                });
                Gpu::synchronize();
            }

            HYPRE_IJVectorSetValues(b, nrows, cell_id_vec[mfi].dataPtr(), rhs_diag[mfi].dataPtr());
        }
    }
}

}  // namespace amrex
