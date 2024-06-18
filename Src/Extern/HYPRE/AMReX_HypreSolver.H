/**
 * \file AMReX_HypreSolver.H
 *
 * \addtogroup HYPRE
 *  @{
 */

#ifndef AMREX_HYPRE_SOLVER_H_
#define AMREX_HYPRE_SOLVER_H_

#include <AMReX_Geometry.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_HypreIJIface.H>
#include <AMReX_BLProfiler.H>

#include "HYPRE.h"
#include "_hypre_utilities.h"

#include <string>

namespace amrex
{

/**
 * @brief Solve Ax = b using HYPRE's generic IJ matrix format where A is a sparse matrix
 * specified using the compressed sparse row (CSR) format.
 *
 * An example of using HypreSolver to solve Poisson's equation is located in
 * Tests/LinearSolvers/Hypre
 *
 */
template <int MSS>
class HypreSolver
{
public:

/**
 * Constructor
 *
 * \param [in]      a_index_type        nodality of the data
 * \param [in]      a_nghost            number of ghosts cells
 * \param [in]      a_geom              Geometry of problem domain
 * \param [in]      a_grids             BoxArray of computational grids
 * \param [in]      a_dmap              DistributionMapping
 * \param [in]      a_marker            functor that returns whether the variable n at (i,j,k) in Box boxno (local index) is valid
 *                                      (i.e., not exactly on Dirichlet boundary). function signature:
 *                                      \code
 *                                          (int boxno, int i, int j, int k, int n) -> bool
 *                                      \endcode
 * \param [in]      a_filler            functor that fills the row in the matrix A for variable n at (i,j,k) in Box boxno (local index)
 *                                      using the CSR format. function signature:
 *                                      \code
 *                                          // [in ] gid    gid[n] is the id for variable n at (i,j,k)
 *                                          // [out] ncols  # of non-zero columns in this row.
 *                                          // [out] cols   array of indices of columns with a non-zero matrix element in this row.
 *                                          // [out] mat    array of (non-zero) matrix elements in this row.
 *                                          (int boxno, int i, int j, int k, int n, Array4<HYPRE_Int const> const* gid, HYPRE_Int& ncols, HYPRE_Int* cols, HYPRE_Real* mat)
 *                                      \endcode
 * \param [in]      a_verbose           HYPRE verbosity (default 0)
 * \param [in]      a_options_namespace namespace to use when parsing runtime input parameters (default "hypre")
 *
 */
    template <class Marker, class Filler>
    HypreSolver (Vector<IndexType>   const& a_index_type,
                 IntVect             const& a_nghost,
                 Geometry            const& a_geom,
                 BoxArray            const& a_grids,
                 DistributionMapping const& a_dmap, // NOLINT(modernize-pass-by-value)
                 Marker                  && a_marker,
                 Filler                  && a_filler,
                 int                        a_verbose = 0,
                 std::string                a_options_namespace = "hypre");

    template <class MF, std::enable_if_t<IsFabArray<MF>::value &&
                                         std::is_same_v<typename MF::value_type,
                                                        HYPRE_Real>, int> = 0>
/**
 * Solve Ax=b
 *
 * \param [out]      a_soln        the solution (i.e. x)
 * \param [in]       a_rhs         the right-hand side (i.e. b)
 * \param [in]       rel_tol       relative convergence tolerance
 * \param [in]       abs_tol       absolute convergence tolerance
 * \param [in]       max_iter      maximum number of iterations
 *
 */
    void solve (Vector<MF      *> const& a_soln,
                Vector<MF const*> const& a_rhs,
                HYPRE_Real rel_tol, HYPRE_Real abs_tol, int max_iter);

    int getNumIters () const { return m_hypre_ij->getNumIters(); }

    HYPRE_Real getFinalResidualNorm () const {
        return m_hypre_ij->getFinalResidualNorm();
    }

    HYPRE_IJMatrix getA () const { return m_hypre_ij->A(); }
    HYPRE_IJVector getb () const { return m_hypre_ij->b(); }
    HYPRE_IJVector getx () const { return m_hypre_ij->x(); }

// public: // for cuda

    template <class Marker>
#ifdef AMREX_USE_CUDA
    std::enable_if_t<IsCallable<Marker,int,int,int,int,int>::value>
#else
    std::enable_if_t<IsCallableR<bool,Marker,int,int,int,int,int>::value>
#endif
    fill_local_id (Marker const& marker);

    void fill_global_id ();

    template <class Filler,
              std::enable_if_t<IsCallable<Filler,int,int,int,int,int,
                                          Array4<HYPRE_Int const> const*,
                                          HYPRE_Int&, HYPRE_Int*,
                                          HYPRE_Real*>::value,int> FOO = 0>
    void fill_matrix (Filler const& filler);

    template <class MF, std::enable_if_t<IsFabArray<MF>::value &&
                                         std::is_same_v<typename MF::value_type,
                                                        HYPRE_Real>, int> = 0>
    void load_vectors (Vector<MF      *> const& a_soln,
                       Vector<MF const*> const& a_rhs);

    template <class MF, std::enable_if_t<IsFabArray<MF>::value &&
                                         std::is_same_v<typename MF::value_type,
                                                        HYPRE_Real>, int> = 0>
    void get_solution (Vector<MF*> const& a_soln);

private:

    int                 m_nvars;
    Vector<IndexType>   m_index_type;
    IntVect             m_nghost;
    Geometry            m_geom;
    Vector<BoxArray>    m_grids;
    DistributionMapping m_dmap;

    int         m_verbose;
    std::string m_options_namespace;

    MPI_Comm m_comm = MPI_COMM_NULL;

    Vector<std::unique_ptr<iMultiFab>>        m_owner_mask;
    Vector<iMultiFab>                         m_local_id;
    Vector<FabArray<BaseFab<HYPRE_Int>>>      m_global_id;
    LayoutData<Gpu::DeviceVector<HYPRE_Int>>  m_global_id_vec;

#ifdef AMREX_USE_GPU
    LayoutData<Gpu::DeviceVector<int>> m_cell_offset;
#endif

    Vector<LayoutData<HYPRE_Int>> m_nrows_grid;
    Vector<LayoutData<HYPRE_Int>> m_id_offset;
    LayoutData<HYPRE_Int>         m_nrows;
    HYPRE_Int                     m_nrows_proc;

    std::unique_ptr<HypreIJIface> m_hypre_ij;

    // Non-owning references to Hypre matrix, rhs, and solution data
    HYPRE_IJMatrix m_A = nullptr;
    HYPRE_IJVector m_b = nullptr;
    HYPRE_IJVector m_x = nullptr;
};

template <int MSS>
template <class Marker, class Filler>
HypreSolver<MSS>::HypreSolver (Vector<IndexType>   const& a_index_type,
                               IntVect             const& a_nghost,
                               Geometry            const& a_geom,
                               BoxArray            const& a_grids,
                               DistributionMapping const& a_dmap, // NOLINT(modernize-pass-by-value)
                               Marker                  && a_marker,
                               Filler                  && a_filler,
                               int                        a_verbose,
                               std::string                a_options_namespace)
    : m_nvars            (int(a_index_type.size())),
      m_index_type       (a_index_type),
      m_nghost           (a_nghost),
      m_geom             (a_geom),
      m_dmap             (a_dmap),
      m_verbose          (a_verbose),
      m_options_namespace(std::move(a_options_namespace))
{
    BL_PROFILE("HypreSolver()");

#ifdef AMREX_USE_MPI
    m_comm = ParallelContext::CommunicatorSub(); // NOLINT(cppcoreguidelines-prefer-member-initializer)
#endif

    m_grids.resize(m_nvars);
    m_local_id.resize(m_nvars);
    m_global_id.resize(m_nvars);
    m_nrows_grid.resize(m_nvars);
    m_id_offset.resize(m_nvars);
    Long nrows_max = 0;
    for (int ivar = 0; ivar < m_nvars; ++ivar) {
        m_grids     [ivar] = amrex::convert(a_grids,m_index_type[ivar]);
        m_local_id  [ivar].define(m_grids[ivar], m_dmap, 1, 0);
        m_global_id [ivar].define(m_grids[ivar], m_dmap, 1, m_nghost);
        m_nrows_grid[ivar].define(m_grids[0], m_dmap);
        m_id_offset [ivar].define(m_grids[0], m_dmap);
        nrows_max += m_grids[ivar].numPts();
    }
    m_global_id_vec.define(m_grids[0], m_dmap);
    m_nrows.define        (m_grids[0], m_dmap);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nrows_max < static_cast<Long>(std::numeric_limits<HYPRE_Int>::max()-1),
                                     "Need to configure Hypre with --enable-bigint");

    m_owner_mask.resize(m_nvars);
    for (int ivar = 0; ivar < m_nvars; ++ivar) {
        m_owner_mask[ivar] = amrex::OwnerMask(m_local_id[ivar], m_geom.periodicity());
    }

#ifdef AMREX_USE_GPU
    m_cell_offset.define(m_grids[0], m_dmap);
#endif

    fill_local_id(std::forward<Marker>(a_marker));

    // At this point, m_local_id stores the ids local to each box.
    // m_nrows_grid stores the number of unique points in each box.
    // m_nrows_proc is the number of rowss for all variables on this MPI
    // process.  If a point is invalid, its id is invalid (i.e., a very
    // negative number).  Note that the data type of local_node_id is int,
    // not HYPRE_Int for performance on GPU.

    const int nprocs = ParallelContext::NProcsSub();
    const int myproc = ParallelContext::MyProcSub();

    Vector<HYPRE_Int> nrows_allprocs(nprocs);
#ifdef AMREX_USE_MPI
    if (nrows_allprocs.size() > 1) {
        MPI_Allgather(&m_nrows_proc, sizeof(HYPRE_Int), MPI_CHAR,
                      nrows_allprocs.data(), sizeof(HYPRE_Int), MPI_CHAR, m_comm);
    } else
#endif
    {
        nrows_allprocs[0] = m_nrows_proc;
    }

    HYPRE_Int proc_begin = 0;
    for (int i = 0; i < myproc; ++i) {
        proc_begin += nrows_allprocs[i];
    }

    HYPRE_Int proc_end = proc_begin;
    for (MFIter mfi(m_nrows_grid[0]); mfi.isValid(); ++mfi) {
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            m_id_offset[ivar][mfi] = proc_end;
            proc_end += m_nrows_grid[ivar][mfi];
        }
    }
    AMREX_ASSERT(proc_end == proc_begin + m_nrows_proc);

    fill_global_id();

    // Create and initialize A, b & x
    HYPRE_Int ilower = proc_begin;
    HYPRE_Int iupper = proc_end-1;
    m_hypre_ij = std::make_unique<HypreIJIface>(m_comm, ilower, iupper, m_verbose);
    m_hypre_ij->parse_inputs(m_options_namespace);

    // Obtain non-owning references to the matrix, rhs, and solution data
    m_A = m_hypre_ij->A();
    m_b = m_hypre_ij->b();
    m_x = m_hypre_ij->x();

    fill_matrix(std::forward<Filler>(a_filler));
}

template <int MSS>
template <class Marker>
#ifdef AMREX_USE_CUDA
    std::enable_if_t<IsCallable<Marker,int,int,int,int,int>::value>
#else
    std::enable_if_t<IsCallableR<bool,Marker,int,int,int,int,int>::value>
#endif
HypreSolver<MSS>::fill_local_id (Marker const& marker)
{
    BL_PROFILE("HypreSolver::fill_local_id()");

#ifdef AMREX_USE_GPU

    for (MFIter mfi(m_local_id[0]); mfi.isValid(); ++mfi) {
        int boxno = mfi.LocalIndex();
        Long npts_tot = 0;
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            Box const& bx = amrex::convert(mfi.validbox(),m_index_type[ivar]);
            npts_tot += bx.numPts();
        }
        m_cell_offset[mfi].resize(npts_tot);
        npts_tot = 0;
        int* p_cell_offset = m_cell_offset[mfi].data();
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            Box const& bx = amrex::convert(mfi.validbox(),m_index_type[ivar]);
            auto const& lid = m_local_id[ivar].array(mfi);
            auto const& owner = m_owner_mask[ivar]->const_array(mfi);
            AMREX_ASSERT(bx.numPts() < static_cast<Long>(std::numeric_limits<int>::max()));
            const auto npts = static_cast<int>(bx.numPts());
            int npts_box = amrex::Scan::PrefixSum<int>(npts,
                [=] AMREX_GPU_DEVICE (int offset) noexcept -> int
                {
                    const Dim3 cell = bx.atOffset(offset).dim3();
                    int id = (owner (      cell.x,cell.y,cell.z     ) &&
                              marker(boxno,cell.x,cell.y,cell.z,ivar)) ? 1 : 0;
                    lid(cell.x,cell.y,cell.z) = id;
                    return id;
                },
                [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
                {
                    const Dim3 cell = bx.atOffset(offset).dim3();
                    if (lid(cell.x,cell.y,cell.z)) {
                        lid(cell.x,cell.y,cell.z) = ps;
                        p_cell_offset[ps] = offset;
                    } else {
                        lid(cell.x,cell.y,cell.z) = std::numeric_limits<int>::lowest();
                    }
                },
                amrex::Scan::Type::exclusive);
            m_nrows_grid[ivar][mfi] = npts_box;
            npts_tot += npts_box;
            p_cell_offset += npts_box;
        }
        m_cell_offset[mfi].resize(npts_tot);
    }

#else

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(m_local_id[0]); mfi.isValid(); ++mfi) {
        int boxno = mfi.LocalIndex();
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            Box const& bx = amrex::convert(mfi.validbox(),m_index_type[ivar]);
            auto const& lid = m_local_id[ivar].array(mfi);
            auto const& owner = m_owner_mask[ivar]->const_array(mfi);
            int id = 0;
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (owner(i,j,k) && marker(boxno,i,j,k,ivar)) {
                    lid(i,j,k) = id++;
                } else {
                    lid(i,j,k) = std::numeric_limits<int>::lowest();
                }
            }}}
            m_nrows_grid[ivar][mfi] = id;
        }
    }
#endif

    m_nrows_proc = 0;
    for (MFIter mfi(m_nrows); mfi.isValid(); ++mfi) {
        int nrows = 0;
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            nrows += m_nrows_grid[ivar][mfi];
        }
        m_nrows[mfi] = nrows;
        m_nrows_proc += nrows;
    }
}

template <int MSS>
void
HypreSolver<MSS>::fill_global_id ()
{
    BL_PROFILE("HypreSolver::fill_global_id()");

    // To generate global ids for Hypre, we need to remove duplicates on
    // nodes shared by multiple Boxes with OverrideSync.  So we need to use
    // a type that supports atomicAdd.  HYPRE_Int is either int or long
    // long.  The latter does not have native atomicAdd support in CUDA/HIP.
    using AtomicInt = std::conditional_t<sizeof(HYPRE_Int) == 4,
                                         HYPRE_Int, unsigned long long>;

    Vector<FabArray<BaseFab<AtomicInt>>> global_id_raii;
    Vector<FabArray<BaseFab<AtomicInt>>*> p_global_id;
    if constexpr (std::is_same<HYPRE_Int, AtomicInt>()) {
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            p_global_id.push_back(&(m_global_id[ivar]));
        }
    } else {
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            global_id_raii[ivar].define(m_global_id[ivar].boxArray(),
                                        m_global_id[ivar].DistributionMap(),
                                        1, m_global_id[ivar].nGrowVect());
            p_global_id.push_back(&(global_id_raii[ivar]));
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_global_id[0]); mfi.isValid(); ++mfi) {
        auto& rows_vec = m_global_id_vec[mfi];
        rows_vec.resize(m_nrows[mfi]);

        HYPRE_Int nrows = 0;
        for (int ivar = 0; ivar < m_nvars; ++ivar) {
            HYPRE_Int const os = m_id_offset[ivar][mfi];
            Box bx = mfi.validbox();
            bx.convert(m_index_type[ivar]).grow(m_nghost);
            Array4<AtomicInt> const& gid = p_global_id[ivar]->array(mfi);
            auto const& lid = m_local_id[ivar].const_array(mfi);
            HYPRE_Int* rows = rows_vec.data() + nrows;
            nrows += m_nrows_grid[ivar][mfi];
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (lid.contains(i,j,k) && lid(i,j,k) >= 0) {
                    const auto id = lid(i,j,k) + os;
                    rows[lid(i,j,k)] = id;
                    gid(i,j,k) = static_cast<AtomicInt>(id);
                } else {
                    gid(i,j,k) = static_cast<AtomicInt>
                        (std::numeric_limits<HYPRE_Int>::max());
                }
            });
        }
    }

    for (int ivar = 0; ivar < m_nvars; ++ivar) {
        amrex::OverrideSync(*p_global_id[ivar], *m_owner_mask[ivar],
                            m_geom.periodicity());
        p_global_id[ivar]->FillBoundary(m_geom.periodicity());

        if constexpr (!std::is_same<HYPRE_Int, AtomicInt>()) {
            auto const& dst = m_global_id[ivar].arrays();
            auto const& src = p_global_id[ivar]->const_arrays();
            amrex::ParallelFor(m_global_id[ivar], m_global_id[ivar].nGrowVect(),
                               [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                dst[b](i,j,k) = static_cast<HYPRE_Int>(src[b](i,j,k));
            });
        }
    }
}

#ifdef AMREX_USE_GPU
namespace detail {
template <typename T>
void pack_matrix_gpu (Gpu::DeviceVector<HYPRE_Int>& cols_tmp,
                      Gpu::DeviceVector<HYPRE_Real> mat_tmp,
                      Gpu::DeviceVector<HYPRE_Int>& cols,
                      Gpu::DeviceVector<HYPRE_Real>& mat)
{
    auto* p_cols_tmp = cols_tmp.data();
    auto* p_mat_tmp  =  mat_tmp.data();
    auto const* p_cols = cols.data();
    auto const* p_mat  =  mat.data();
    const auto N = Long(cols.size());
    Scan::PrefixSum<T>(N,
        [=] AMREX_GPU_DEVICE (Long i) -> T
        {
            return static_cast<T>(p_cols[i] >= 0);
        },
        [=] AMREX_GPU_DEVICE (Long i, T s)
        {
            if (p_cols[i] >= 0) {
                p_cols_tmp[s] = p_cols[i];
                p_mat_tmp[s] = p_mat[i];
            }
        },
        Scan::Type::exclusive, Scan::noRetSum);
    std::swap(cols_tmp, cols);
    std::swap(mat_tmp, mat);
}
}
#endif

template <int MSS>
template <class Filler,
          std::enable_if_t<IsCallable<Filler,int,int,int,int,int,
                                      Array4<HYPRE_Int const> const*,
                                      HYPRE_Int&, HYPRE_Int*,
                                      HYPRE_Real*>::value,int> FOO>
void
HypreSolver<MSS>::fill_matrix (Filler const& filler)
{
    BL_PROFILE("HypreSolver::fill_matrix()");

    Gpu::DeviceVector<HYPRE_Int> ncols_vec;
    Gpu::DeviceVector<HYPRE_Int> cols_vec;
    Gpu::DeviceVector<HYPRE_Real> mat_vec;

    MFItInfo mfitinfo;
    mfitinfo.UseDefaultStream().DisableDeviceSync();
    for (MFIter mfi(m_local_id[0],mfitinfo); mfi.isValid(); ++mfi)
    {
        int boxno = mfi.LocalIndex();
        const HYPRE_Int nrows = m_nrows[mfi];
        if (nrows > 0)
        {
            ncols_vec.clear();
            ncols_vec.resize(nrows);
            HYPRE_Int* ncols = ncols_vec.data();

            cols_vec.clear();
            cols_vec.resize(Long(nrows)*MSS, -1);
            HYPRE_Int* cols = cols_vec.data();

            mat_vec.clear();
            mat_vec.resize(Long(nrows)*MSS);
            HYPRE_Real* mat = mat_vec.data();

            Vector<Array4<HYPRE_Int const>> gid_v(m_nvars);
            for (int ivar = 0; ivar < m_nvars; ++ivar) {
                gid_v[ivar] = m_global_id[ivar].const_array(mfi);
            }

#ifdef AMREX_USE_GPU
            const Gpu::Buffer<Array4<HYPRE_Int const>> gid_buf
                (gid_v.data(), gid_v.size());
            auto const* pgid = gid_buf.data();
            auto const* p_cell_offset = m_cell_offset[mfi].data();
            Long ntot = 0;
            for (int ivar = 0; ivar < m_nvars; ++ivar) {
                const HYPRE_Int nrows_var = m_nrows_grid[ivar][mfi];
                if (nrows_var > 0) {
                    Box const& bx = amrex::convert(mfi.validbox(),m_index_type[ivar]);
                    ntot += Reduce::Sum<Long>(nrows_var,
                                              [=] AMREX_GPU_DEVICE (HYPRE_Int offset)
                    {
                        const Dim3 cell = bx.atOffset(p_cell_offset[offset]).dim3();
                        filler(boxno, cell.x, cell.y, cell.z, ivar, pgid,
                               ncols[offset], cols+Long(offset)*MSS,
                               mat+Long(offset)*MSS);
                        return ncols[offset];
                    });
                    p_cell_offset += nrows_var;
                    ncols += nrows_var;
                    cols += Long(nrows_var)*MSS;
                    mat  += Long(nrows_var)*MSS;
                }
            }
            Gpu::DeviceVector<HYPRE_Int> cols_tmp(ntot);
            Gpu::DeviceVector<HYPRE_Real> mat_tmp(ntot);
            if (ntot >= Long(std::numeric_limits<int>::max())) {
                detail::pack_matrix_gpu<Long>(cols_tmp, mat_tmp, cols_vec, mat_vec);
            } else {
                detail::pack_matrix_gpu<int>(cols_tmp, mat_tmp, cols_vec, mat_vec);
            }
#else
            auto* pgid = gid_v.data();
            for (int ivar = 0; ivar < m_nvars; ++ivar) {
                if (m_nrows_grid[ivar][mfi] > 0) {
                    auto const& lid = m_local_id[ivar].const_array(mfi);
                    amrex::Loop(amrex::convert(mfi.validbox(),m_index_type[ivar]),
                                [=,&ncols,&cols,&mat] (int i, int j, int k)
                    {
                        if (lid(i,j,k) >= 0) {
                            filler(boxno, i, j, k, ivar, pgid, *ncols, cols, mat);
                            cols += (*ncols);
                            mat  += (*ncols);
                            ++ncols;
                        }
                    });
                }
            }
#endif

            const auto& rows_vec = m_global_id_vec[mfi];
            HYPRE_Int const* rows = rows_vec.data();

            Gpu::streamSynchronize();
            HYPRE_IJMatrixSetValues(m_A, nrows, ncols_vec.data(), rows,
                                    cols_vec.data(), mat_vec.data());
            Gpu::synchronize();
        }
    }
    HYPRE_IJMatrixAssemble(m_A);
}

template <int MSS>
template <class MF, std::enable_if_t<IsFabArray<MF>::value &&
                                     std::is_same_v<typename MF::value_type,
                                                    HYPRE_Real>, int> FOO>
void
HypreSolver<MSS>::solve (Vector<MF      *> const& a_soln,
                         Vector<MF const*> const& a_rhs,
                         HYPRE_Real rel_tol, HYPRE_Real abs_tol, int max_iter)
{
    BL_PROFILE("HypreSolver::solve()");

    AMREX_ASSERT(a_soln.size() == m_nvars && a_rhs.size() == m_nvars);

    HYPRE_IJVectorInitialize(m_b);
    HYPRE_IJVectorInitialize(m_x);

    load_vectors(a_soln, a_rhs);

    HYPRE_IJVectorAssemble(m_x);
    HYPRE_IJVectorAssemble(m_b);

    m_hypre_ij->solve(rel_tol, abs_tol, max_iter);

    get_solution(a_soln);
}

template <int MSS>
template <class MF, std::enable_if_t<IsFabArray<MF>::value &&
                                     std::is_same_v<typename MF::value_type,
                                                    HYPRE_Real>, int> FOO>
void
HypreSolver<MSS>::load_vectors (Vector<MF      *> const& a_soln,
                                Vector<MF const*> const& a_rhs)
{
    BL_PROFILE("HypreSolver::load_vectors()");

    MFItInfo mfitinfo;
    mfitinfo.UseDefaultStream().DisableDeviceSync();

    Gpu::DeviceVector<Real> xvec;
    Gpu::DeviceVector<Real> bvec;
    for (MFIter mfi(*a_soln[0],mfitinfo); mfi.isValid(); ++mfi)
    {
        const HYPRE_Int nrows = m_nrows[mfi];
        if (nrows > 0)
        {
            xvec.clear();
            xvec.resize(nrows);
            bvec.clear();
            bvec.resize(nrows);
            auto* xp = xvec.data();
            auto* bp = bvec.data();

            HYPRE_Int const* rows = m_global_id_vec[mfi].data();

            HYPRE_Int offset = 0;
            for (int ivar = 0; ivar < m_nvars; ++ivar) {
                if (m_nrows_grid[ivar][mfi] > 0) {
                    auto const& xfab = a_soln[ivar]->const_array(mfi);
                    auto const& bfab = a_rhs [ivar]->const_array(mfi);
                    auto const& lid = m_local_id[ivar].const_array(mfi);
                    HYPRE_Real* x = xp + offset;
                    HYPRE_Real* b = bp + offset;
                    Box box = amrex::convert(mfi.validbox(),m_index_type[ivar]);
                    amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (lid(i,j,k) >= 0) {
                            x[lid(i,j,k)] = xfab(i,j,k);
                            b[lid(i,j,k)] = bfab(i,j,k);
                        }
                    });
                    offset += m_nrows_grid[ivar][mfi];
                }
            }

            Gpu::streamSynchronize();
            HYPRE_IJVectorSetValues(m_x, nrows, rows, xp);
            HYPRE_IJVectorSetValues(m_b, nrows, rows, bp);
            Gpu::synchronize();
        }
    }
}

template <int MSS>
template <class MF, std::enable_if_t<IsFabArray<MF>::value &&
                                     std::is_same_v<typename MF::value_type,
                                                    HYPRE_Real>, int> FOO>
void
HypreSolver<MSS>::get_solution (Vector<MF*> const& a_soln)
{
    BL_PROFILE("HypreSolver::get_solution()");

    MFItInfo mfitinfo;
    mfitinfo.UseDefaultStream().DisableDeviceSync();

    Gpu::DeviceVector<Real> xvec;
    for (MFIter mfi(*a_soln[0],mfitinfo); mfi.isValid(); ++mfi)
    {
        const HYPRE_Int nrows = m_nrows[mfi];
        if (nrows > 0)
        {
            xvec.clear();
            xvec.resize(nrows);
            auto* xp = xvec.data();

            HYPRE_Int const* rows = m_global_id_vec[mfi].data();

            HYPRE_IJVectorGetValues(m_x, nrows, rows, xp);
            Gpu::synchronize();

            HYPRE_Int offset = 0;
            for (int ivar = 0; ivar < m_nvars; ++ivar) {
                if (m_nrows_grid[ivar][mfi] > 0) {
                    auto const& xfab = a_soln[ivar]->array(mfi);
                    auto const& lid = m_local_id[ivar].const_array(mfi);
                    HYPRE_Real* x = xp + offset;
                    Box box = amrex::convert(mfi.validbox(),m_index_type[ivar]);
                    amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (lid(i,j,k) >= 0) {
                            xfab(i,j,k) = x[lid(i,j,k)];
                        }
                    });
                    offset += m_nrows_grid[ivar][mfi];
                }
            }
            Gpu::streamSynchronize();
        }
    }

    for (int ivar = 0; ivar < m_nvars; ++ivar) {
        amrex::OverrideSync(*a_soln[ivar], *m_owner_mask[ivar],
                            m_geom.periodicity());
    }
}

}

#endif
/** @} */
