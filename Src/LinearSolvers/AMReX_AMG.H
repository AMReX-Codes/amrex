#ifndef AMREX_AMG_H_
#define AMREX_AMG_H_
#include <AMReX_Config.H>

#include <AMReX_Algebra.H>
#include <AMReX_BiCGStab_MV.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_CArena.H>
#include <AMReX_GMRES_MV.H>
#include <AMReX_PCG_MV.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Random.H>
#include <AMReX_Scan.H>
#include <AMReX_Smoother_MV.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace amrex {

namespace detail {

//! True if row `row` of `csr` has an undecided neighbor with a larger PMIS
//! weight than `wi`. Ties go to the larger global index. `gmap` maps the
//! column index to the global index, or, if null, `goff` is added.
template <typename TC, typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool pmis_blocked (CsrView<TC const,I> const& csr, Long row, T wi, Long ig,
                   T const* AMREX_RESTRICT w, int const* AMREX_RESTRICT state,
                   Long const* AMREX_RESTRICT gmap, Long goff, int undecided)
{
    for (Long idx = csr.row_offset[row]; idx < csr.row_offset[row+1]; ++idx) {
        Long j = csr.col_index[idx];
        Long jg = gmap ? gmap[j] : j + goff;
        if (jg != ig && state[j] == undecided &&
            (w[j] > wi || (w[j] == wi && jg > ig))) {
            return true;
        }
    }
    return false;
}

//! True if row `row` of `csr` has a neighbor with state `c`.
template <typename T, typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool csr_row_has_state (CsrView<T const,I> const& csr, Long row,
                        int const* AMREX_RESTRICT state, int c)
{
    for (Long idx = csr.row_offset[row]; idx < csr.row_offset[row+1]; ++idx) {
        if (state[csr.col_index[idx]] == c) { return true; }
    }
    return false;
}

//! Bytes in use in The_Arena, or 0 if it does not report usage.
inline std::size_t arena_bytes_in_use ()
{
    if (auto* carena = dynamic_cast<CArena*>(The_Arena())) {
        return carena->heap_space_actually_used();
    }
    return 0;
}

//! True if local row `i` has a remote part in `csr1`.
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool has_remote_row (ParCsr<T const> const& a, Long i)
{
    return a.csr1.nnz > 0 && a.row_map[i] >= 0;
}

//! Walks a row of A together with the matching row of S (both sorted
//! by column) to find the strong C-point neighbors.
template <typename T, typename TS>
struct DirectInterpRow
{
    ParCsr<T const> a;
    ParCsr<TS const> s;
    int const* AMREX_RESTRICT state;    // local rows
    Long const* AMREX_RESTRICT cidx;    // local coarse index, local rows
    int const* AMREX_RESTRICT state_r;  // A's remote columns
    int nlocal_c;                        // number of local coarse columns
    int coarse;

    //! Calls f(compact coarse column, a_ij, is strong C) for each
    //! off-diagonal entry of row i and sets `diag` to a_ii. The compact
    //! column is the local coarse index, or nlocal_c plus A's remote index.
    template <typename F>
    AMREX_GPU_HOST_DEVICE void walk (Long i, T& diag, F const& f) const
    {
        diag = T(0);
        {
            auto const& ac = a.csr0;
            auto const& sc = s.csr0;
            Long sp = sc.row_offset[i];
            Long const se = sc.row_offset[i+1];
            for (Long idx = ac.row_offset[i]; idx < ac.row_offset[i+1]; ++idx) {
                Long const j = ac.col_index[idx];
                T const v = ac.mat[idx];
                if (j == i) { diag = v; continue; }
                while (sp < se && sc.col_index[sp] < j) { ++sp; }
                bool const strong_c = (sp < se && sc.col_index[sp] == j
                                       && state[j] == coarse);
                f(int(cidx[j]), v, strong_c);
            }
        }
        if (has_remote_row(a, i)) {
            auto const& ac = a.csr1;
            auto const& sc = s.csr1;
            Long const ii = a.row_map[i];
            Long sp = 0, se = 0;
            if (has_remote_row(s, i)) {
                Long const jj = s.row_map[i];
                sp = sc.row_offset[jj];
                se = sc.row_offset[jj+1];
            }
            for (Long idx = ac.row_offset[ii]; idx < ac.row_offset[ii+1]; ++idx) {
                Long const j = ac.col_index[idx];
                Long const jg = a.col_map[j];
                T const v = ac.mat[idx];
                while (sp < se && s.col_map[sc.col_index[sp]] < jg) { ++sp; }
                bool const strong_c = (sp < se && s.col_map[sc.col_index[sp]] == jg
                                       && state_r[j] == coarse);
                f(int(nlocal_c + j), v, strong_c);
            }
        }
    }
};

//! Walks row i of the strong matrix S together with row i of S^T (both
//! sorted by column). For each strong neighbor j it provides a_ij, a_ji
//! (0 if j does not strongly depend on i), and j's state, D_beta and
//! global coarse index.
template <typename T, typename TS>
struct MMInterpRow
{
    ParCsr<TS const> s;
    ParCsr<TS const> st;
    int const* AMREX_RESTRICT state;
    T const* AMREX_RESTRICT dbeta;
    Long const* AMREX_RESTRICT cidx;    // local coarse index, local rows
    int const* AMREX_RESTRICT state_r;  // S's remote columns
    T const* AMREX_RESTRICT dbeta_r;
    int nlocal_f;                        // local fine columns
    int nlocal_c;                        // local coarse columns

    //! f(compact fine column, compact coarse column, a_ij, a_ji, state_j,
    //! dbeta_j). Compact columns are the local index, or nlocal plus S's
    //! remote index.
    template <typename F>
    AMREX_GPU_HOST_DEVICE void walk (Long i, F const& f) const
    {
        {
            auto const& sc = s.csr0;
            auto const& tc = st.csr0;
            Long tp = tc.row_offset[i];
            Long const te = tc.row_offset[i+1];
            for (Long idx = sc.row_offset[i]; idx < sc.row_offset[i+1]; ++idx) {
                Long const j = sc.col_index[idx];
                while (tp < te && tc.col_index[tp] < j) { ++tp; }
                T const aji = (tp < te && tc.col_index[tp] == j) ? T(tc.mat[tp]) : T(0);
                f(int(j), int(cidx[j]), T(sc.mat[idx]), aji, state[j], dbeta[j]);
            }
        }
        if (has_remote_row(s, i)) {
            auto const& sc = s.csr1;
            auto const& tc = st.csr1;
            Long const ii = s.row_map[i];
            Long tp = 0, te = 0;
            if (has_remote_row(st, i)) {
                Long const jj = st.row_map[i];
                tp = tc.row_offset[jj];
                te = tc.row_offset[jj+1];
            }
            for (Long idx = sc.row_offset[ii]; idx < sc.row_offset[ii+1]; ++idx) {
                Long const j = sc.col_index[idx];
                Long const jg = s.col_map[j];
                while (tp < te && st.col_map[tc.col_index[tp]] < jg) { ++tp; }
                T const aji = (tp < te && st.col_map[tc.col_index[tp]] == jg)
                    ? T(tc.mat[tp]) : T(0);
                f(int(nlocal_f + j), int(nlocal_c + j), T(sc.mat[idx]), aji,
                  state_r[j], dbeta_r[j]);
            }
        }
    }
};

}

/**
 * \brief Algebraic multigrid solver for SpMatrix / AlgVector systems.
 *
 * Setup uses the classical strength of connection [1], PMIS coarsening
 * [2], distance-two interpolation [3] in the matrix-matrix form of [4]
 * (or classical direct interpolation [1]), R = P^T and A_c = R A P,
 * repeated until the coarsest level has at most `max_coarse_size` rows.
 * On the first `aggressive_num_levels` levels, aggressive coarsening [5]
 * applies PMIS a second time on the distance-two graph of the C-points,
 * with a two-stage interpolation P = P1 P2 built with the same matrix-matrix
 * formulas. The solve runs V-cycles with Chebyshev [7] (default),
 * l1-Jacobi [6], weighted Jacobi or, on CPUs, l1 hybrid Gauss-Seidel [6]
 * smoothing. The coarsest level is solved with smoother sweeps (default),
 * BiCGStab or GMRES. With `setKrylovSolver`, the V-cycle is used as a
 * preconditioner for BiCGStab, GMRES or conjugate gradients (PCG)
 * instead of being iterated on its own, which roughly halved the number
 * of V-cycles in the tests in Tests/Algebra/AMG.
 *
 * Singular matrices with the constant null vector (`setSingular`) are
 * handled by removing the mean of the right-hand side and of the
 * coarsest-level residuals; interpolation preserves constants, so the
 * null vector is the same on every level.
 *
 * Interpolation rows are truncated to at most `p_max_elmts` entries (4 by
 * default, as in hypre) and entries below `trunc_factor` times the row
 * maximum are dropped; kept weights are rescaled to preserve the row
 * sums [3].
 *
 * Parameters are set per object with the setters below.
 *
 * The matrix must be square with rows sorted by column index (see
 * SpMatrix), and it must stay alive while the solver is in use.
 *
 * References:
 * [1] J. W. Ruge and K. Stüben, Algebraic multigrid, in Multigrid
 *     Methods, S. F. McCormick, ed., SIAM, 1987, pp. 73-130.
 *     doi:10.1137/1.9781611971057.ch4
 * [2] H. De Sterck, U. M. Yang and J. J. Heys, Reducing complexity in
 *     parallel algebraic multigrid preconditioners, SIAM J. Matrix Anal.
 *     Appl. 27 (2006), pp. 1019-1039. doi:10.1137/040615729
 * [3] H. De Sterck, R. D. Falgout, J. W. Nolting and U. M. Yang,
 *     Distance-two interpolation for parallel algebraic multigrid, Numer.
 *     Linear Algebra Appl. 15 (2008), pp. 115-139. doi:10.1002/nla.559
 * [4] R. Li, B. Sjögreen and U. M. Yang, A new class of AMG interpolation
 *     methods based on matrix-matrix multiplications, SIAM J. Sci. Comput.
 *     43 (2021), pp. S540-S564. doi:10.1137/20M134931X
 * [5] U. M. Yang, On long-range interpolation operators for aggressive
 *     coarsening, Numer. Linear Algebra Appl. 17 (2010), pp. 453-472.
 *     doi:10.1002/nla.689
 * [6] A. H. Baker, R. D. Falgout, T. V. Kolev and U. M. Yang, Multigrid
 *     smoothers for ultraparallel computing, SIAM J. Sci. Comput. 33
 *     (2011), pp. 2864-2887. doi:10.1137/100798806
 * [7] M. Adams, M. Brezina, J. Hu and R. Tuminaro, Parallel multigrid
 *     smoothing: polynomial versus Gauss-Seidel, J. Comput. Phys. 188
 *     (2003), pp. 593-610. doi:10.1016/S0021-9991(03)00194-3
 */
template <typename T>
class AMG
{
public:
    //! Value type of the strength matrices: float when T is double. Their
    //! entries are the matrix entries divided by `strength_scale`.
    using TS = std::conditional_t<std::is_same_v<T,double>, float, T>;

    enum class BottomSolver { Jacobi, BiCGStab, GMRES };
    //! None: iterate V-cycles. Otherwise use one V-cycle as the
    //! preconditioner of the Krylov solver. PCG needs a symmetric positive
    //! definite matrix and a symmetric cycle (same pre- and post-smoother,
    //! smoother sweeps at the bottom).
    enum class KrylovSolver { None, BiCGStab, GMRES, PCG };
    //! Direct: classical direct interpolation. MMExt and MMExtI:
    //! extended and extended+i interpolation in matrix-matrix form.
    enum class InterpType { Direct, MMExt, MMExtI };
    //! Jacobi: weighted Jacobi, weight 2/3. L1Jacobi: l1-Jacobi [6], weight
    //! 4/3 by default: the l1 scaling bounds the eigenvalues of the scaled
    //! operator by one, so 4/3 is the optimal damping for the upper half of
    //! the spectrum and the sweep converges for any SPD matrix.
    //! Chebyshev (default): polynomial of degree `cheby_degree` (2) on the
    //! l1-scaled operator, damping [lambda_max / cheby_ratio, lambda_max]
    //! with ratio 6; one sweep costs `cheby_degree` matrix-vector products.
    //! One sweep per pre- and post-smoothing step by default, versus two
    //! for the Jacobi smoothers, so the cost per cycle is the same.
    //! L1GaussSeidel: l1 hybrid Gauss-Seidel, forward before and backward
    //! after the coarse correction; CPU builds only.
    enum class Smoother { Jacobi, L1Jacobi, Chebyshev, L1GaussSeidel };

    AMG () = default;
    explicit AMG (SpMatrix<T>& a_mat);

    //! Define the solver with a square matrix. Must be called before solve.
    void define (SpMatrix<T>& a_mat);

    //! Solve A x = b. `a_sol` holds the initial guess on entry.
    void solve (AlgVector<T>& a_sol, AlgVector<T> const& a_rhs);

    void setVerbose (int v) { m_verbose = v; }
    void setMaxIter (int n) { m_maxiter = n; }
    //! If positive, run exactly this many V-cycles regardless of tolerance.
    void setFixedIter (int n) { m_fixediter = n; }
    void setRelTol (T t) { m_reltol = t; }
    //! Takes effect at the next setup.
    void setSmoother (Smoother sm) { m_smoother_type = sm; m_prepared = false; }
    //! Relaxation weight of the smoother; negative restores the default. Next setup.
    void setRelaxWeight (T w) { m_relax_weight = w; m_prepared = false; }
    //! Chebyshev polynomial degree. Next setup.
    void setChebyshevDegree (int d) { m_cheby_degree = d; m_prepared = false; }
    //! Chebyshev eigenvalue ratio lambda_max / lambda_min. Next setup.
    void setChebyshevRatio (T r) { m_cheby_ratio = r; m_prepared = false; }
    //! Number of pre-smoothing sweeps; negative restores the smoother's default.
    void setPreSmooth (int nu) { m_nu1 = nu; }
    //! Number of post-smoothing sweeps; negative restores the smoother's default.
    void setPostSmooth (int nu) { m_nu2 = nu; }
    //! Number of smoother sweeps on the coarsest level with the Jacobi
    //! bottom solver; negative restores the default.
    void setBottomSmooth (int nu) { m_nu_bottom = nu; }
    //! Takes effect at the next setup, i.e., the next solve after define.
    void setStrongThreshold (T theta) { m_strong_threshold = theta; m_prepared = false; }
    //! Maximum number of levels. Two means a two-level method. Takes
    //! effect at the next setup, i.e., the next solve after define.
    void setMaxLevels (int n) { m_max_levels = n; m_prepared = false; }
    //! Stop coarsening when a level has at most this many rows. Next setup.
    void setMaxCoarseSize (Long n) { m_max_coarse_size = n; m_prepared = false; }
    //! Number of levels with aggressive coarsening, starting at the finest. Next setup.
    void setAggressiveNumLevels (int n) { m_aggressive_num_levels = n; m_prepared = false; }
    //! Use classical direct interpolation for the first stage of the
    //! two-stage aggressive interpolation (cheaper setup, more cycles).
    void setAggressiveDirectInterp (bool b) { m_aggressive_direct_interp = b; m_prepared = false; }
    //! Singular matrix with the constant vector in its null space (e.g.
    //! Poisson with periodic or Neumann boundaries). The mean is removed
    //! from the right-hand side and the coarsest residuals, and the
    //! solution is returned with zero mean.
    void setSingular (bool b) { m_singular = b; }
    //! Takes effect at the next setup.
    void setInterpType (InterpType it) { m_interp_type = it; m_prepared = false; }
    //! Maximum number of entries per row of P (0: no limit). Next setup.
    void setPMaxElmts (int n) { m_p_max_elmts = n; m_prepared = false; }
    //! Drop interpolation weights below this fraction of the row maximum. Next setup.
    void setTruncFactor (T f) { m_trunc_factor = f; m_prepared = false; }
    //! Krylov solver with the V-cycle as preconditioner (see KrylovSolver). Default None.
    void setKrylovSolver (KrylovSolver a) {
        m_krylov = a;
        m_krylov_bicgstab.reset();
        m_krylov_gmres.reset();
        m_krylov_pcg.reset();
    }
    //! Apply one V-cycle to `b` with zero initial guess: x = M^{-1} b.
    //! This is the preconditioner used by the Krylov solver; usable
    //! stand-alone after the first solve or setup.
    void precond (AlgVector<T>& x, AlgVector<T> const& b);
    void setBottomSolver (BottomSolver bs) {
        m_bottom_solver = bs;
        m_gmres.reset();
        m_bicgstab.reset();
    }
    void setBottomTol (T t) { m_bottom_reltol = t; }
    void setBottomMaxIter (int n) { m_bottom_maxiter = n; }
    void setBottomVerbose (int v) { m_bottom_verbose = v; }
    void setThrowException (bool b) { m_throw_exception = b; }

    [[nodiscard]] int getNumIters () const { return m_niters; }
    [[nodiscard]] T getResidualNorm () const { return m_rnorm; }
    [[nodiscard]] int numLevels () const { return int(m_a.size()); }
    //! Peak bytes in use in The_Arena during setup (0 if not available).
    [[nodiscard]] std::size_t getSetupPeakBytes () const { return m_setup_peak_bytes; }
    //! Global number of rows on a level. Valid after setup.
    [[nodiscard]] Long numGlobalRows (int lev) const { return m_a[lev]->numGlobalRows(); }
    [[nodiscard]] SpMatrix<T> const& getMatrix (int lev) const { return *m_a[lev]; }
    [[nodiscard]] SpMatrix<T> const& getInterp (int lev) const { return m_P[lev]; }

private:

    enum PointType : int { undecided = 0, fine = 1, coarse = 2 };

    bool m_prepared = false;
    bool m_throw_exception = false;
    int m_verbose = 0;
    int m_maxiter = 200;
    int m_fixediter = 0;
    int m_niters = 0;

    Smoother m_smoother_type = Smoother::Chebyshev;
    T m_relax_weight = T(-1);
    int m_cheby_degree = 2;
    T m_cheby_ratio = T(6);
    int m_nu1 = -1;
    int m_nu2 = -1;
    int m_nu_bottom = -1;

    [[nodiscard]] int numPreSmooth () const {
        return (m_nu1 >= 0) ? m_nu1 : ((m_smoother_type == Smoother::Chebyshev) ? 1 : 2);
    }
    [[nodiscard]] int numPostSmooth () const {
        return (m_nu2 >= 0) ? m_nu2 : ((m_smoother_type == Smoother::Chebyshev) ? 1 : 2);
    }
    [[nodiscard]] int numBottomSmooth () const {
        return (m_nu_bottom >= 0) ? m_nu_bottom : ((m_smoother_type == Smoother::Chebyshev) ? 4 : 8);
    }

    int m_max_levels = 25;
    Long m_max_coarse_size = 9;
    int m_aggressive_num_levels = 0;
    bool m_aggressive_direct_interp = false;
    bool m_singular = false;
    std::unique_ptr<AlgVector<T>> m_bproj; // rhs with the mean removed
    T m_strong_threshold = T(0.25);

    InterpType m_interp_type = InterpType::MMExtI;
    int m_p_max_elmts = 4;
    T m_trunc_factor = T(0);
    BottomSolver m_bottom_solver = BottomSolver::Jacobi;
    T m_bottom_reltol = T(1.e-4);
    int m_bottom_maxiter = 200;
    int m_bottom_verbose = 0;

    T m_reltol = std::is_same_v<T,float> ? T(1.e-4) : T(1.e-8);
    T m_init_rnorm = std::numeric_limits<T>::max();
    T m_rnorm = std::numeric_limits<T>::max();
    T m_bnorm = std::numeric_limits<T>::max();

    std::string m_print_ident;
    std::size_t m_setup_peak_bytes = 0;

    //! Records the arena usage and, at verbose >= 2, prints it with a label.
    void note_memory (int lev, char const* stage)
    {
        auto const bytes = detail::arena_bytes_in_use();
        m_setup_peak_bytes = std::max(m_setup_peak_bytes, bytes);
        if (m_verbose >= 2 && bytes > 0) {
            std::size_t reserved = 0;
            if (auto* carena = dynamic_cast<CArena*>(The_Arena())) { reserved = carena->heap_space_used(); }
            amrex::Print() << m_print_ident << "AMG: level " << lev << ": after " << stage
                           << ": arena in use " << double(bytes)/double(1<<20)
                           << " MB, reserved " << double(reserved)/double(1<<20) << " MB"
#ifdef AMREX_USE_GPU
                           << ", device free " << double(Gpu::Device::freeMemAvailable())/double(1<<20) << " MB"
#endif
                           << "\n";
        }
    }

    SpMatrix<T>* m_mat = nullptr;
    Vector<SpMatrix<T>*> m_a;                        // matrix on each level
    Vector<std::unique_ptr<SpMatrix<T>>> m_crsemat;  // owns levels 1..n-1
    Vector<SpMatrix<T>> m_P;   // m_P[lev] interpolates lev+1 -> lev
    Vector<SpMatrix<T>> m_R;   // m_R[lev] restricts lev -> lev+1

    Vector<AlgVector<T>> m_res;    // rhs - A(sol)
    Vector<AlgVector<T>> m_cor;    // cor = A^-1 (res)
    Vector<AlgVector<T>> m_rescor; // res - A(cor)

    Vector<std::unique_ptr<JacobiSmoother<T>>> m_smoother;
    Vector<std::unique_ptr<ChebyshevSmoother<T>>> m_cheby;
    Vector<std::unique_ptr<L1GaussSeidelSmoother<T>>> m_l1gs;
    std::unique_ptr<GMRES_MV<T>> m_gmres;
    std::unique_ptr<BiCGStab_MV<T>> m_bicgstab;

    KrylovSolver m_krylov = KrylovSolver::None;
    std::unique_ptr<BiCGStab_MV<T>> m_krylov_bicgstab;
    std::unique_ptr<GMRES_MV<T>> m_krylov_gmres;
    std::unique_ptr<PCG_MV<T>> m_krylov_pcg;

    void vcycle ();
    void solve_krylov (AlgVector<T>& a_sol, AlgVector<T> const& a_rhs, T res_target);
    void smooth (int lev, AlgVector<T>& xvec, AlgVector<T> const& bvec, int nu,
                 bool with_initial_guess, bool backward = false);
    void bottom_solve ();

    static std::pair<T,T> compute_rbnorms (SpMatrix<T> const& mat,
                                           AlgVector<T> const& xvec,
                                           AlgVector<T> const& bvec,
                                           AlgVector<T>& rvec);

    static T compute_rnorm (SpMatrix<T> const& mat, AlgVector<T> const& xvec,
                            AlgVector<T> const& bvec, AlgVector<T>& rvec);

//Public for CUDA
public: // NOLINT

    void setup ();

    //! Strong part of A: S_ij = a_ij / strength_scale if i strongly depends
    //! on j, else absent.
    SpMatrix<TS> create_soc (SpMatrix<T> const& A) const;

    //! Largest |a_ii|, the unit of the strength matrix entries.
    static T strength_scale (SpMatrix<T> const& A);

    //! PMIS C/F selection. `state` gets PointType values. Points without
    //! any strong connection become F, or C if `isolated_as_coarse`.
    void pmis (SpMatrix<TS>& S, SpMatrix<TS>& ST, Gpu::DeviceVector<int>& state,
               bool isolated_as_coarse = false) const;

    //! Local index among the rows of type `which` and their partition.
    static AlgPartition coarse_numbering (Gpu::DeviceVector<int> const& state,
                                         Gpu::DeviceVector<Long>& cidx,
                                         int which = PointType::coarse);

    //! Singular systems: remove the mean of `x` (the constant null vector
    //! is the same on every level since interpolation preserves constants).
    void project_out (AlgVector<T>& x)
    {
        if (!m_singular) { return; }
        Long const n = x.numLocalRows();
        auto* AMREX_RESTRICT px = x.data();
        T sum = Reduce::Sum<T>(n, [=] AMREX_GPU_DEVICE (Long i) -> T { return px[i]; });
        ParallelAllReduce::Sum(sum, ParallelContext::CommunicatorSub());
        T const mean = sum / T(x.partition().numGlobalRows());
        ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { px[i] -= mean; });
        Gpu::streamSynchronize();
    }

    //! Classical direct interpolation.
    static SpMatrix<T> build_interp (SpMatrix<T>& A, SpMatrix<TS> const& S,
                                     Gpu::DeviceVector<int> const& state,
                                     Gpu::DeviceVector<Long> const& cidx,
                                     AlgPartition const& cpart);

    //! Keep the `max_elmts` largest weights of each row of P (0: all) and
    //! drop weights below `trunc_factor` times the row maximum. Kept
    //! weights are rescaled to preserve the row's positive and negative sums.
    static SpMatrix<T> truncate_interp (SpMatrix<T>& P, AlgPartition const& cpart,
                                        int max_elmts, T trunc_factor);

    //! Strength graph between C-points connected through a strong path of
    //! length two, in C-point numbering (rows partitioned by cpart).
    static SpMatrix<TS> second_pass_strength (SpMatrix<TS>& S,
                                             Gpu::DeviceVector<int> const& state,
                                             Gpu::DeviceVector<Long> const& cidx,
                                             AlgPartition const& cpart);

    //! Rows of P at the C-points, in C-point numbering.
    static SpMatrix<T> select_coarse_rows (SpMatrix<T> const& P,
                                           Gpu::DeviceVector<int> const& state,
                                           Gpu::DeviceVector<Long> const& cidx,
                                           AlgPartition const& cpart);

    //! Interpolation of the selected type from the C-points in `state`.
    SpMatrix<T> interpolation (SpMatrix<T>& A, SpMatrix<TS>& S, SpMatrix<TS> const& ST,
                               Gpu::DeviceVector<int> const& state,
                               Gpu::DeviceVector<Long> const& cidx,
                               AlgPartition const& cpart) const;

    //! Extended (plus_i = false) or extended+i interpolation in
    //! matrix-matrix form. ST is the transpose of S.
    static SpMatrix<T> build_interp_mm (SpMatrix<T>& A, SpMatrix<TS>& S,
                                        SpMatrix<TS> const& ST,
                                        Gpu::DeviceVector<int> const& state,
                                        Gpu::DeviceVector<Long> const& cidx,
                                        AlgPartition const& cpart, bool plus_i);
};

template <typename T>
AMG<T>::AMG (SpMatrix<T>& a_mat)
{
    define(a_mat);
}

template <typename T>
void AMG<T>::define (SpMatrix<T>& a_mat)
{
    m_mat = &a_mat;
    m_prepared = false;
    m_bproj.reset();
    m_krylov_bicgstab.reset();
    m_krylov_gmres.reset();
    m_krylov_pcg.reset();
}

template <typename T>
void AMG<T>::solve (AlgVector<T>& a_sol, AlgVector<T> const& a_rhs)
{
    BL_PROFILE("AMG::solve");

    if (!m_prepared) { setup(); }

    // Singular: iterate on the rhs with its mean removed.
    AlgVector<T> const* prhs = &a_rhs;
    if (m_singular) {
        if (!m_bproj) { m_bproj = std::make_unique<AlgVector<T>>(m_mat->partition()); }
        m_bproj->copy(a_rhs);
        project_out(*m_bproj);
        prhs = m_bproj.get();
    }
    AlgVector<T> const& rhs = *prhs;

    auto norms = compute_rbnorms(*m_mat, a_sol, rhs, m_res[0]);
    m_bnorm = norms.first;
    m_init_rnorm = m_rnorm = norms.second;
    m_niters = 0;

    if (m_bnorm == 0) {
        a_sol.setVal(0); // also the zero-mean solution of a singular system
        m_rnorm = 0;
        return;
    }

    if (m_verbose >= 1) {
        amrex::Print() << m_print_ident << "AMG: Initial rhs      2-norm = " << m_bnorm << "\n"
                       << m_print_ident << "AMG: Initial residual 2-norm = " << m_rnorm << "\n";
    }

    T res_target = m_bnorm * m_reltol;
    if (m_fixediter <= 0 && m_rnorm <= res_target) {
        if (m_verbose >= 1) {
            amrex::Print() << m_print_ident << "AMG: No iterations needed\n";
        }
        project_out(a_sol);
        return;
    }

    if (m_krylov != KrylovSolver::None) {
        solve_krylov(a_sol, rhs, res_target);
        project_out(a_sol);
        return;
    }

    int const niters = (m_fixediter > 0) ? m_fixediter : m_maxiter;
    for (int iter = 0; iter < niters; ++iter)
    {
        vcycle();
        a_sol.plusAsync(m_cor[0]);

        m_rnorm = compute_rnorm(*m_mat, a_sol, rhs, m_res[0]);
        m_niters = iter+1;

        if (m_verbose >= 2) {
            amrex::Print() << m_print_ident << "AMG: Iteration " << std::setw(3)
                           << iter+1 << ": rnorm, rnorm/bnorm = "
                           << m_rnorm << ", " << m_rnorm/m_bnorm << "\n";
        }

        bool converged = (m_fixediter <= 0) && (m_rnorm <= res_target);
        if (converged) {
            if (m_verbose >= 1) {
                amrex::Print() << m_print_ident << "AMG: Final Iter. " << iter+1
                               << ": rnorm, rnorm/bnorm = "
                               << m_rnorm << ", " << m_rnorm/m_bnorm << "\n";
            }
            break;
        } else if (m_rnorm > T(1.e20) * m_init_rnorm || std::isnan(m_rnorm)) {
            if (m_verbose > 0) {
                amrex::Print() << m_print_ident << "AMG: Failing to converge after " << iter+1 << " iterations."
                               << " rnorm, rnorm/bnorm = "
                               << m_rnorm << ", " << m_rnorm/m_bnorm << "\n";
            }
            if (m_throw_exception) {
                throw std::runtime_error("AMG blew up.");
            } else {
                amrex::Abort("AMG failing so lets stop here");
            }
        }
    }

    if (m_verbose >= 1 && m_fixediter <= 0 && m_rnorm > res_target) {
        amrex::Print() << m_print_ident << "AMG: Not converged after " << m_niters
                       << " iterations. rnorm/bnorm = " << m_rnorm/m_bnorm << "\n";
    }
    project_out(a_sol);
}

template <typename T>
void AMG<T>::precond (AlgVector<T>& x, AlgVector<T> const& b)
{
    if (!m_prepared) { setup(); }
    m_res[0].copyAsync(b);
    vcycle();
    x.copyAsync(m_cor[0]);
    project_out(x);
}

template <typename T>
void AMG<T>::solve_krylov (AlgVector<T>& a_sol, AlgVector<T> const& a_rhs, T res_target)
{
    BL_PROFILE("AMG::solve_krylov");

    // Bound on every solve so that a moved AMG object stays usable.
    auto pc = [this] (AlgVector<T>& x, AlgVector<T> const& b) { this->precond(x, b); };
    int const niters = (m_fixediter > 0) ? m_fixediter : m_maxiter;
    T const atol = (m_fixediter > 0) ? T(0) : res_target;
    int const verbose = std::max(0, m_verbose-1);
    T krylov_rnorm = T(0);

    // Relative tolerance 0: convergence is decided by the absolute target
    // reltol*|b|, the same criterion as the V-cycle iteration. The Krylov
    // residual is a recurrence, which can converge before the true one when
    // the V-cycle is not a fixed linear operator (Krylov bottom solver). In
    // that case, restart from the current solution.
    auto run = [&] (auto& s) {
        s.setVerbose(verbose);
        s.setInitialGuessNonzero(true);
        m_niters = 0;
        while (true) {
            s.setMaxIters(niters - m_niters);
            s.solve(a_sol, a_rhs, T(0), atol);
            m_niters += s.getNumIters();
            krylov_rnorm = s.getResidualNorm();
            m_rnorm = compute_rnorm(*m_mat, a_sol, a_rhs, m_res[0]);
            if (m_fixediter > 0 || m_rnorm <= res_target || krylov_rnorm > res_target ||
                m_niters >= niters || s.getNumIters() == 0 || std::isnan(m_rnorm)) {
                break;
            }
            if (m_verbose >= 2) {
                amrex::Print() << m_print_ident << "AMG: restarting Krylov solver at iteration "
                               << m_niters << ", rnorm/bnorm = " << m_rnorm/m_bnorm << "\n";
            }
        }
    };

    char const* name = "";
    if (m_krylov == KrylovSolver::BiCGStab) {
        name = "BiCGStab";
        if (!m_krylov_bicgstab) {
            m_krylov_bicgstab = std::make_unique<BiCGStab_MV<T>>(m_mat);
        }
        m_krylov_bicgstab->setPrecond(pc);
        run(m_krylov_bicgstab->getSolver());
    } else if (m_krylov == KrylovSolver::GMRES) {
        name = "GMRES";
        if (!m_krylov_gmres) {
            m_krylov_gmres = std::make_unique<GMRES_MV<T>>(m_mat);
        }
        m_krylov_gmres->setPrecond(pc);
        run(m_krylov_gmres->getSolver());
    } else {
        name = "PCG";
        if (!m_krylov_pcg) {
            m_krylov_pcg = std::make_unique<PCG_MV<T>>(m_mat);
        }
        m_krylov_pcg->setPrecond(pc);
        run(m_krylov_pcg->getSolver());
    }

    if (std::isnan(m_rnorm) || m_rnorm > T(1.e20) * m_init_rnorm) {
        if (m_throw_exception) {
            throw std::runtime_error("AMG blew up.");
        } else {
            amrex::Abort("AMG failing so lets stop here");
        }
    }

    if (m_verbose >= 1) {
        amrex::Print() << m_print_ident << "AMG(" << name << "): Final Iter. " << m_niters
                       << ": rnorm, rnorm/bnorm = " << m_rnorm << ", " << m_rnorm/m_bnorm;
        if (m_verbose >= 2) {
            amrex::Print() << " (Krylov recurrence " << krylov_rnorm << ")";
        }
        amrex::Print() << "\n";
        if (m_fixediter <= 0 && m_rnorm > res_target) {
            amrex::Print() << m_print_ident << "AMG(" << name << "): Not converged after "
                           << m_niters << " iterations.\n";
        }
    }
}

template <typename T>
std::pair<T,T> AMG<T>::compute_rbnorms (SpMatrix<T> const& mat,
                                        AlgVector<T> const& xvec,
                                        AlgVector<T> const& bvec,
                                        AlgVector<T>& rvec)
{
    auto bnorm = bvec.norm2(true);

    amrex::computeResidual(rvec, mat, xvec, bvec);
    auto rnorm = rvec.norm2(true);

    if (ParallelContext::NProcsSub() > 1) {
        bnorm = Math::powi<2>(bnorm);
        rnorm = Math::powi<2>(rnorm);
        ParallelAllReduce::Sum<T>({bnorm,rnorm}, ParallelContext::CommunicatorSub());
        bnorm = std::sqrt(bnorm);
        rnorm = std::sqrt(rnorm);
    }

    return {bnorm, rnorm};
}

template <typename T>
T AMG<T>::compute_rnorm (SpMatrix<T> const& mat, AlgVector<T> const& xvec,
                         AlgVector<T> const& bvec, AlgVector<T>& rvec)
{
    amrex::computeResidual(rvec, mat, xvec, bvec);
    return rvec.norm2();
}

template <typename T>
void AMG<T>::vcycle () // Solve A(m_cor[0]) = m_res[0] with a V-cycle
{
    BL_PROFILE("AMG::vcycle");

    int const nlev = numLevels();
    int const bottom = nlev-1;

    if (nlev == 1) {
        bottom_solve();
        return;
    }

    smooth(0, m_cor[0], m_res[0], numPreSmooth(), false); // cor = A^-1 (res) w/o initial guess

    for (int lev = 0; lev < bottom; ++lev) {
        computeResidual(m_rescor[lev], *m_a[lev], m_cor[lev], m_res[lev]);
        SpMV(m_res[lev+1], m_R[lev], m_rescor[lev]);
        if (lev+1 == bottom) {
            bottom_solve();
        } else {
            smooth(lev+1, m_cor[lev+1], m_res[lev+1], numPreSmooth(), false);
        }
    }

    for (int lev = bottom-1; lev >= 0; --lev) {
        SpMV(m_rescor[lev], m_P[lev], m_cor[lev+1]);
        m_cor[lev].plusAsync(m_rescor[lev]);
        smooth(lev, m_cor[lev], m_res[lev], numPostSmooth(), true, true);
    }
}

template <typename T>
void AMG<T>::smooth (int lev, AlgVector<T>& xvec, AlgVector<T> const& bvec, int nu,
                     bool with_initial_guess, bool backward)
{
    if (m_smoother_type == Smoother::L1GaussSeidel) {
        auto& smoother = *m_l1gs[lev];
        smoother.setNumIters(nu);
        smoother(xvec, bvec, with_initial_guess, backward);
    } else if (m_smoother_type == Smoother::Chebyshev) {
        auto& smoother = *m_cheby[lev];
        smoother.setNumIters(nu);
        smoother(xvec, bvec, with_initial_guess);
    } else {
        auto& smoother = *m_smoother[lev];
        smoother.setNumIters(nu);
        smoother(xvec, bvec, with_initial_guess);
    }
}

template <typename T>
void AMG<T>::bottom_solve ()
{
    BL_PROFILE("AMG::bottom_solve");

    int const lev = numLevels()-1;
    project_out(m_res[lev]); // keep the coarsest solve consistent
    switch (m_bottom_solver) {
    case BottomSolver::Jacobi:
        smooth(lev, m_cor[lev], m_res[lev], numBottomSmooth(), false);
        break;
    case BottomSolver::GMRES:
        if (!m_gmres) {
            m_gmres = std::make_unique<GMRES_MV<T>>(m_a[lev]);
            m_gmres->setPrecond(JacobiSmoother<T>(m_a[lev], m_smoother_type != Smoother::Jacobi));
        }
        m_gmres->setVerbose(m_bottom_verbose);
        m_gmres->getSolver().setMaxIters(m_bottom_maxiter);
        m_cor[lev].setVal(0);
        m_gmres->solve(m_cor[lev], m_res[lev], m_bottom_reltol, T(0));
        break;
    case BottomSolver::BiCGStab:
        if (!m_bicgstab) {
            m_bicgstab = std::make_unique<BiCGStab_MV<T>>(m_a[lev]);
            m_bicgstab->setPrecond(JacobiSmoother<T>(m_a[lev], m_smoother_type != Smoother::Jacobi));
        }
        m_bicgstab->setVerbose(m_bottom_verbose);
        m_bicgstab->getSolver().setMaxIters(m_bottom_maxiter);
        m_cor[lev].setVal(0);
        m_bicgstab->solve(m_cor[lev], m_res[lev], m_bottom_reltol, T(0));
        break;
    }
    project_out(m_cor[lev]);
}

template <typename T>
void AMG<T>::setup ()
{
    BL_PROFILE("AMG::setup");

    m_a.clear();
    m_crsemat.clear();
    m_P.clear();
    m_R.clear();
    m_res.clear();
    m_cor.clear();
    m_rescor.clear();
    m_smoother.clear();
    m_cheby.clear();
    m_l1gs.clear();
    m_gmres.reset();
    m_bicgstab.reset();

    m_a.push_back(m_mat);
    m_mat->setColumnPartition(m_mat->partition()); // square: local column indexing
    m_setup_peak_bytes = 0;
    note_memory(0, "start");

    // Reference: Reducing Complexity in Parallel Algebraic Multigrid
    // Preconditioners, De Sterck, Yang & Heys, SIAM J. Matrix Anal. Appl.,
    // 27(4), 2006. https://doi.org/10.1137/040615729

    for (int lev = 0; lev+1 < m_max_levels; ++lev)
    {
        auto& A = *m_a[lev];
        if (A.numGlobalRows() <= m_max_coarse_size) { break; }

        // The strength matrices are released before the triple product to
        // keep the peak memory down.
        SpMatrix<T> P;
        AlgPartition cpart;
        bool aggressive = false;
        Long ncoarse = 0;
        {
        auto S = create_soc(A);
        S.setColumnPartition(S.partition());
        note_memory(lev, "S");
        auto ST = amrex::transpose(S, S.partition());
        ST.setColumnPartition(ST.partition());

        note_memory(lev, "strength");
        Gpu::DeviceVector<int> state;
        pmis(S, ST, state);
        note_memory(lev, "PMIS");

        Gpu::DeviceVector<Long> cidx;
        cpart = coarse_numbering(state, cidx);
        ncoarse = cpart.numGlobalRows();
        if (ncoarse == 0 || ncoarse == A.numGlobalRows()) { break; }

        if (lev < m_aggressive_num_levels) {
            // Second PMIS pass on the C-points connected at distance two.
            auto S2 = second_pass_strength(S, state, cidx, cpart);
            S2.setColumnPartition(cpart);
            auto S2T = amrex::transpose(S2, cpart);
            S2T.setColumnPartition(cpart);
            note_memory(lev, "second-pass strength");
            Gpu::DeviceVector<int> state2;
            pmis(S2, S2T, state2, true); // isolated C-points stay coarse
            Gpu::DeviceVector<Long> cidx2;
            auto cpart2 = coarse_numbering(state2, cidx2);
            Long const nc2 = cpart2.numGlobalRows();
            if (nc2 > 0 && nc2 < ncoarse) {
                // C/F split of the fine points with respect to C2.
                Long const nrows = A.numLocalRows();
                Gpu::DeviceVector<int> statef(nrows);
                Gpu::DeviceVector<Long> cidxf(nrows);
                auto const* p1 = state.data();
                auto const* c1 = cidx.data();
                auto const* p2 = state2.data();
                auto const* c2 = cidx2.data();
                auto* pf = statef.data();
                auto* cf = cidxf.data();
                ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
                {
                    bool const c = (p1[i] == PointType::coarse) &&
                                   (p2[c1[i]] == PointType::coarse);
                    pf[i] = c ? PointType::coarse : PointType::fine;
                    cf[i] = (p1[i] == PointType::coarse) ? c2[c1[i]] : Long(0);
                });
                Gpu::streamSynchronize();

                bool const trunc = (m_p_max_elmts > 0 || m_trunc_factor > T(0));
                auto P1 = m_aggressive_direct_interp ? build_interp(A, S, state, cidx, cpart)
                                                     : interpolation(A, S, ST, state, cidx, cpart);
                if (trunc) { P1 = truncate_interp(P1, cpart, m_p_max_elmts, m_trunc_factor); }
                note_memory(lev, "P1");
                auto Pf = build_interp_mm(A, S, ST, statef, cidxf, cpart2,
                                          m_interp_type != InterpType::MMExt);
                if (trunc) { Pf = truncate_interp(Pf, cpart2, m_p_max_elmts, m_trunc_factor); }
                auto P2 = select_coarse_rows(Pf, state, cidx, cpart);
                note_memory(lev, "P2");
                P = amrex::SpGEMM(P1, P2, cpart2);
                cpart = cpart2;
                aggressive = true;
            }
        }
        if (!aggressive) {
            P = interpolation(A, S, ST, state, cidx, cpart);
        }
        note_memory(lev, "interpolation");
        }
        // The strength matrices were large; hand their chunks back so the
        // triple product can get contiguous memory.
        The_Arena()->freeUnused();
        note_memory(lev, "strength release");
        if (m_p_max_elmts > 0 || m_trunc_factor > T(0)) {
            Long nnz0 = P.numLocalNonZeros();
            P = truncate_interp(P, cpart, m_p_max_elmts, m_trunc_factor);
            note_memory(lev, "truncation");
            if (m_verbose >= 2) {
                Long nnz[2] = {nnz0, P.numLocalNonZeros()};
                ParallelAllReduce::Sum(nnz, 2, ParallelContext::CommunicatorSub());
                amrex::Print() << m_print_ident << "AMG: level " << lev
                               << ": P truncated from " << nnz[0] << " to "
                               << nnz[1] << " nonzeros\n";
            }
        }
        if (m_verbose >= 2 && aggressive) {
            amrex::Print() << m_print_ident << "AMG: level " << lev
                           << ": aggressive coarsening, " << ncoarse << " -> "
                           << cpart.numGlobalRows() << " C-points\n";
        }
        if (m_verbose >= 3) {
            // P should reproduce constants: row sums close to one.
            AlgVector<T> ones(cpart), rs(A.partition());
            ones.setVal(T(1));
            SpMV(rs, P, ones);
            Long nzero = Reduce::Sum<Long>(rs.numLocalRows(),
                [q = rs.data()] AMREX_GPU_DEVICE (Long i) -> Long { return q[i] == T(0); });
            T rmin = Reduce::Min<T>(rs.numLocalRows(),
                [q = rs.data()] AMREX_GPU_DEVICE (Long i) -> T { return q[i]; });
            T rmax = Reduce::Max<T>(rs.numLocalRows(),
                [q = rs.data()] AMREX_GPU_DEVICE (Long i) -> T { return q[i]; });
            ParallelAllReduce::Sum(nzero, ParallelContext::CommunicatorSub());
            ParallelAllReduce::Min(rmin, ParallelContext::CommunicatorSub());
            ParallelAllReduce::Max(rmax, ParallelContext::CommunicatorSub());
            amrex::Print() << m_print_ident << "AMG: level " << lev << ": P row sums in ["
                           << rmin << ", " << rmax << "], " << nzero << " empty rows\n";
        }
        auto R  = amrex::transpose(P, cpart);
        note_memory(lev, "R");
        // R (A P): measured smaller and faster than (R A) P, whose
        // intermediate gathers several fine rows per coarse row.
        SpMatrix<T> Ac;
        {
            auto AP = amrex::SpGEMM(A, P, cpart);
            note_memory(lev, "A*P");
            Ac = amrex::SpGEMM(R, AP, cpart);
        }
        note_memory(lev, "coarse operator");

        m_P.push_back(std::move(P));
        m_R.push_back(std::move(R));
        m_crsemat.push_back(std::make_unique<SpMatrix<T>>(std::move(Ac)));
        m_a.push_back(m_crsemat.back().get());
        The_Arena()->freeUnused();
    }

    int const nlev = numLevels();
    for (int lev = 0; lev < nlev; ++lev) {
        m_res.emplace_back(m_a[lev]->partition());
        m_cor.emplace_back(m_a[lev]->partition());
        if (lev < nlev-1) {
            m_rescor.emplace_back(m_a[lev]->partition());
        }
        if (m_smoother_type == Smoother::L1GaussSeidel) {
#ifdef AMREX_USE_GPU
            amrex::Abort("AMG: L1GaussSeidel smoother is only available in CPU builds");
#endif
            m_l1gs.push_back(std::make_unique<L1GaussSeidelSmoother<T>>(m_a[lev]));
        } else if (m_smoother_type == Smoother::Chebyshev) {
            m_cheby.push_back(std::make_unique<ChebyshevSmoother<T>>(m_a[lev], true));
            m_cheby.back()->setDegree(m_cheby_degree);
            m_cheby.back()->setEigRatio(m_cheby_ratio);
            if (m_verbose >= 2) {
                amrex::Print() << m_print_ident << "AMG: level " << lev
                               << ": Chebyshev lambda_max = " << m_cheby.back()->lambdaMax() << "\n";
            }
        } else {
            bool const l1 = (m_smoother_type == Smoother::L1Jacobi);
            m_smoother.push_back(std::make_unique<JacobiSmoother<T>>(m_a[lev], l1));
            T const w = (m_relax_weight >= T(0)) ? m_relax_weight : (l1 ? T(4./3.) : T(2./3.));
            m_smoother.back()->setWeight(w);
        }
    }

    if (m_verbose >= 1) {
        Long nnz0 = 0, nnz_total = 0;
        for (int lev = 0; lev < nlev; ++lev) {
            Long nnz = m_a[lev]->numLocalNonZeros();
            ParallelAllReduce::Sum(nnz, ParallelContext::CommunicatorSub());
            if (lev == 0) { nnz0 = nnz; }
            nnz_total += nnz;
            amrex::Print() << m_print_ident << "AMG: level " << lev << ": "
                           << m_a[lev]->numGlobalRows() << " rows, "
                           << nnz << " nonzeros\n";
        }
        amrex::Print() << m_print_ident << "AMG: operator complexity = "
                       << double(nnz_total)/double(std::max(nnz0,Long(1))) << "\n";
        if (m_setup_peak_bytes > 0) {
            auto peak = m_setup_peak_bytes;
            ParallelAllReduce::Max(peak, ParallelContext::CommunicatorSub());
            amrex::Print() << m_print_ident << "AMG: setup peak arena use = "
                           << double(peak)/double(1<<20) << " MB ("
                           << double(peak)/double(m_mat->numGlobalRows())
                           << " B per row)\n";
        }
    }

    m_prepared = true;
}

template <typename T>
T AMG<T>::strength_scale (SpMatrix<T> const& A)
{
    auto const& Ad = A.diagonalVector();
    auto const* AMREX_RESTRICT p = Ad.data();
    T dmax = Reduce::Max<T>(Ad.numLocalRows(),
        [=] AMREX_GPU_DEVICE (Long i) -> T { return std::abs(p[i]); });
    ParallelAllReduce::Max(dmax, ParallelContext::CommunicatorSub());
    return (dmax > T(0)) ? dmax : T(1);
}

template <typename T>
SpMatrix<typename AMG<T>::TS> AMG<T>::create_soc (SpMatrix<T> const& A) const
{
    BL_PROFILE("AMG::create_soc");

    auto const& a = A.const_parcsr();
    Long const nrows = a.csr0.nrows;
    int const nlocal = int(A.numLocalRows());

    auto const& Ad = A.diagonalVector();
    auto const* AMREX_RESTRICT A_ii = Ad.data();
    auto const threshold = m_strong_threshold;
    T const scale = strength_scale(A);

    // Strength test for one entry of row i.
    auto row_amin = [=] AMREX_GPU_DEVICE (Long i) -> T {
        auto const& c0 = a.csr0;
        T const sgn = std::copysign(T(1), A_ii[i]);
        T amin = std::numeric_limits<T>::max();
        for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
            if (c0.col_index[idx] != i) {
                amin = std::min(amin, sgn*c0.mat[idx]);
            }
        }
        if (detail::has_remote_row(a, i)) {
            auto const& c1 = a.csr1;
            Long const ii = a.row_map[i];
            for (Long idx = c1.row_offset[ii]; idx < c1.row_offset[ii+1]; ++idx) {
                amin = std::min(amin, sgn*c1.mat[idx]);
            }
        }
        return amin * threshold;
    };

    using local_csr_type = typename SpMatrix<TS>::local_csr_type;
    local_csr_type csr;
    csr.row_offset.resize(nrows+1);
    auto* AMREX_RESTRICT psrow = csr.row_offset.data();
    Long const nnz = Scan::PrefixSum<Long>
        (nrows+1,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (i >= nrows) { return 0; }
             T const amin = row_amin(i);
             T const sgn = std::copysign(T(1), A_ii[i]);
             Long n = 0;
             auto const& c0 = a.csr0;
             for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
                 if (c0.col_index[idx] != i && sgn*c0.mat[idx] < amin) { ++n; }
             }
             if (detail::has_remote_row(a, i)) {
                 auto const& c1 = a.csr1;
                 Long const ii = a.row_map[i];
                 for (Long idx = c1.row_offset[ii]; idx < c1.row_offset[ii+1]; ++idx) {
                     if (sgn*c1.mat[idx] < amin) { ++n; }
                 }
             }
             return n;
         },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { psrow[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    csr.mat.resize(nnz);
    csr.col_index.resize(nnz);
    csr.nnz = nnz;
    auto* AMREX_RESTRICT ps = csr.mat.data();
    auto* AMREX_RESTRICT pc = csr.col_index.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        T const amin = row_amin(i);
        T const sgn = std::copysign(T(1), A_ii[i]);
        Long p = psrow[i];
        auto const& c0 = a.csr0;
        for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
            if (c0.col_index[idx] != i && sgn*c0.mat[idx] < amin) {
                ps[p] = TS(c0.mat[idx]/scale);
                pc[p] = c0.col_index[idx];
                ++p;
            }
        }
        if (detail::has_remote_row(a, i)) {
            auto const& c1 = a.csr1;
            Long const ii = a.row_map[i];
            for (Long idx = c1.row_offset[ii]; idx < c1.row_offset[ii+1]; ++idx) {
                if (sgn*c1.mat[idx] < amin) {
                    ps[p] = TS(c1.mat[idx]/scale);
                    pc[p] = nlocal + c1.col_index[idx];
                    ++p;
                }
            }
        }
    });
    Gpu::streamSynchronize();

    SpMatrix<TS> S;
    S.define_split(A.partition(), A.partition(), std::move(csr), nlocal, a.col_map, 0);
    return S;
}

template <typename T>
void AMG<T>::pmis (SpMatrix<TS>& S, SpMatrix<TS>& ST, Gpu::DeviceVector<int>& state,
                   bool isolated_as_coarse) const
{
    BL_PROFILE("AMG::pmis");

    Long const nrows = S.numLocalRows();
    int const U = PointType::undecided;
    int const C = PointType::coarse;

    // Measure: number of points that strongly depend on i. A negative
    // measure marks a point with no strong connection at all.
    auto const& s = S.const_parcsr();
    auto const& st = ST.const_parcsr();
    Gpu::DeviceVector<Long> measure(nrows);
    auto* AMREX_RESTRICT pm = measure.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        Long n = st.csr0.row_offset[i+1] - st.csr0.row_offset[i];
        if (detail::has_remote_row(st, i)) {
            Long const ii = st.row_map[i];
            n += st.csr1.row_offset[ii+1] - st.csr1.row_offset[ii];
        }
        Long nout = s.csr0.row_offset[i+1] - s.csr0.row_offset[i];
        if (detail::has_remote_row(s, i)) {
            Long const ii = s.row_map[i];
            nout += s.csr1.row_offset[ii+1] - s.csr1.row_offset[ii];
        }
        pm[i] = (n == 0 && nout == 0) ? Long(-1) : n;
    });

    Gpu::DeviceVector<T> weight(nrows);
    state.resize(nrows);
    Gpu::DeviceVector<int> state2(nrows);
    auto* AMREX_RESTRICT pw = weight.data();
    auto* AMREX_RESTRICT ps = state.data();
    auto* AMREX_RESTRICT ps2 = state2.data();

    ParallelForRNG(nrows, [=] AMREX_GPU_DEVICE (Long i, RandomEngine const& eng) noexcept
    {
        if (pm[i] <= 0) {
            pw[i] = T(0);
            ps[i] = (pm[i] < 0 && isolated_as_coarse) ? PointType::coarse
                                                       : PointType::fine;
        } else {
            pw[i] = T(pm[i]) + T(Random(eng));
            ps[i] = PointType::undecided;
        }
    });

    Long const row_begin = s.row_begin;

    // Ghost weights are fixed; ghost states are gathered every sweep.
    auto const wS  = S.gatherRemote(pw);
    auto const wST = ST.gatherRemote(pw);
    auto const* AMREX_RESTRICT pwS  = wS.data();
    auto const* AMREX_RESTRICT pwST = wST.data();

    Gpu::DeviceVector<int> sS, sST, s2S;

    while (true) {
        sS  = S.gatherRemote(ps);
        sST = ST.gatherRemote(ps);
        auto const* AMREX_RESTRICT psS  = sS.data();
        auto const* AMREX_RESTRICT psST = sST.data();

        // Independent set in S union S^T: keep only own index writes.
        ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
        {
            int st_i = ps[i];
            if (st_i == U) {
                T const wi = pw[i];
                Long const ig = i + row_begin;
                bool blocked =
                    detail::pmis_blocked(s.csr0,  i, wi, ig, pw, ps, nullptr, s.col_begin,  U) ||
                    detail::pmis_blocked(st.csr0, i, wi, ig, pw, ps, nullptr, st.col_begin, U);
                if (!blocked && detail::has_remote_row(s, i)) {
                    blocked = detail::pmis_blocked(s.csr1, s.row_map[i], wi, ig,
                                                   pwS, psS, s.col_map, 0, U);
                }
                if (!blocked && detail::has_remote_row(st, i)) {
                    blocked = detail::pmis_blocked(st.csr1, st.row_map[i], wi, ig,
                                                   pwST, psST, st.col_map, 0, U);
                }
                if (!blocked) { st_i = C; }
            }
            ps2[i] = st_i;
        });

        s2S = S.gatherRemote(ps2);
        auto const* AMREX_RESTRICT ps2S = s2S.data();

        // Undecided points that strongly depend on a C-point become F.
        ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
        {
            int st_i = ps2[i];
            if (st_i == U) {
                if (detail::csr_row_has_state(s.csr0, i, ps2, C) ||
                    (detail::has_remote_row(s, i) &&
                     detail::csr_row_has_state(s.csr1, s.row_map[i], ps2S, C))) {
                    st_i = PointType::fine;
                }
            }
            ps[i] = st_i;
        });

        Long nund = Reduce::Sum<Long>(nrows, [=] AMREX_GPU_DEVICE (Long i) -> Long {
            return (ps[i] == U) ? 1 : 0;
        });
        ParallelAllReduce::Sum(nund, ParallelContext::CommunicatorSub());
        if (nund == 0) { break; }
    }

    if (m_verbose >= 2) {
        sS = S.gatherRemote(ps);
        auto const* AMREX_RESTRICT psS = sS.data();
        Long nc = Reduce::Sum<Long>(nrows, [=] AMREX_GPU_DEVICE (Long i) -> Long {
            return (ps[i] == C) ? 1 : 0;
        });
        // Two C-points must not be strongly connected.
        Long ncc = Reduce::Sum<Long>(nrows, [=] AMREX_GPU_DEVICE (Long i) -> Long {
            if (ps[i] != C) { return 0; }
            return (detail::csr_row_has_state(s.csr0, i, ps, C) ||
                    (detail::has_remote_row(s, i) &&
                     detail::csr_row_has_state(s.csr1, s.row_map[i], psS, C))) ? 1 : 0;
        });
        // F-points without a strong C neighbor get no interpolation.
        Long nf0 = Reduce::Sum<Long>(nrows, [=] AMREX_GPU_DEVICE (Long i) -> Long {
            if (ps[i] != PointType::fine) { return 0; }
            return (detail::csr_row_has_state(s.csr0, i, ps, C) ||
                    (detail::has_remote_row(s, i) &&
                     detail::csr_row_has_state(s.csr1, s.row_map[i], psS, C))) ? 0 : 1;
        });
        Long counts[3] = {nc, ncc, nf0};
        ParallelAllReduce::Sum(counts, 3, ParallelContext::CommunicatorSub());
        amrex::Print() << m_print_ident << "AMG: PMIS: " << counts[0] << " C-points of "
                       << S.numGlobalRows() << ", C-C strong pairs " << counts[1]
                       << ", F-points w/o C neighbor " << counts[2] << "\n";
    }
}

template <typename T>
AlgPartition AMG<T>::coarse_numbering (Gpu::DeviceVector<int> const& state,
                                       Gpu::DeviceVector<Long>& cidx, int which)
{
    Long const nrows = Long(state.size());
    auto const* AMREX_RESTRICT ps = state.data();
    cidx.resize(nrows);
    auto* AMREX_RESTRICT pc = cidx.data();
    Long const nc = Scan::PrefixSum<Long>
        (nrows,
         [=] AMREX_GPU_DEVICE (Long i) -> Long { return (ps[i] == which) ? 1 : 0; },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { pc[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    int const nprocs = ParallelContext::NProcsSub();
    Vector<Long> counts(nprocs, nc);
#ifdef AMREX_USE_MPI
    // One value per rank: an all-gather, not a reduction of a P-vector.
    ParallelAllGather::AllGather(nc, counts.data(), ParallelContext::CommunicatorSub());
#endif
    Vector<Long> rows(nprocs+1);
    rows[0] = 0;
    for (int i = 0; i < nprocs; ++i) { rows[i+1] = rows[i] + counts[i]; }
    return AlgPartition(std::move(rows));
}

template <typename T>
SpMatrix<T> AMG<T>::build_interp (SpMatrix<T>& A, SpMatrix<TS> const& S,
                                  Gpu::DeviceVector<int> const& state,
                                  Gpu::DeviceVector<Long> const& cidx,
                                  AlgPartition const& cpart)
{
    BL_PROFILE("AMG::build_interp");

    // Reference: hypre's hypre_BoomerAMGBuildDirectInterp.

    auto const& a = A.const_parcsr();
    Long const nrows = a.csr0.nrows;
    auto const* AMREX_RESTRICT ps = state.data();
    auto const* AMREX_RESTRICT pc = cidx.data();
    Long const cbegin = cpart.globalRowBegin();
    int const nlocal_c = int(cpart.numLocalRows());

    // Global coarse index of local rows, and ghost copies for A's remote
    // columns (the remote column list of P).
    Gpu::DeviceVector<Long> gcidx(nrows);
    auto* AMREX_RESTRICT pgc = gcidx.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept { pgc[i] = pc[i] + cbegin; });
    auto const state_r = A.gatherRemote(ps);
    auto const gcidx_r = A.gatherRemote(pgc);

    detail::DirectInterpRow<T,TS> const w{a, S.const_parcsr(), ps, pc,
                                       state_r.data(), nlocal_c,
                                       int(PointType::coarse)};

    using local_csr_type = typename SpMatrix<T>::local_csr_type;
    local_csr_type csr;
    csr.row_offset.resize(nrows+1);
    auto* AMREX_RESTRICT prow = csr.row_offset.data();

    Long const nnz = Scan::PrefixSum<Long>
        (nrows+1,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (i >= nrows) { return 0; }
             if (ps[i] == PointType::coarse) { return 1; }
             Long n = 0;
             T diag;
             w.walk(i, diag, [&] (int, T, bool strong_c) { if (strong_c) { ++n; } });
             return n;
         },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { prow[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    csr.mat.resize(nnz);
    csr.col_index.resize(nnz);
    csr.nnz = nnz;
    auto* AMREX_RESTRICT pmat = csr.mat.data();
    auto* AMREX_RESTRICT pcol = csr.col_index.data();

    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        Long p = prow[i];
        if (ps[i] == PointType::coarse) {
            pcol[p] = int(pc[i]);
            pmat[p] = T(1);
            return;
        }

        T sum_n_neg = T(0), sum_n_pos = T(0), sum_p_neg = T(0), sum_p_pos = T(0);
        T diag;
        w.walk(i, diag, [&] (int, T v, bool strong_c) {
            if (v < T(0)) {
                sum_n_neg += v;
                if (strong_c) { sum_p_neg += v; }
            } else {
                sum_n_pos += v;
                if (strong_c) { sum_p_pos += v; }
            }
        });

        T alfa = T(0), beta = T(0);
        if (sum_p_neg != T(0)) { alfa = sum_n_neg / sum_p_neg; } else { diag += sum_n_neg; }
        if (sum_p_pos != T(0)) { beta = sum_n_pos / sum_p_pos; } else { diag += sum_n_pos; }
        if (diag == T(0)) { diag = T(1); } // no interpolation weights anyway

        T unused;
        w.walk(i, unused, [&] (int cj, T v, bool strong_c) {
            if (strong_c) {
                pcol[p] = cj;
                pmat[p] = ((v < T(0)) ? -alfa*v : -beta*v) / diag;
                ++p;
            }
        });
    });
    Gpu::streamSynchronize();

    SpMatrix<T> P;
    P.define_split(A.partition(), cpart, std::move(csr), nlocal_c, gcidx_r.data(), 0);
    return P;
}

namespace detail {

//! Selects the entries of a row of P that survive truncation. Entries come
//! from the local and remote blocks; columns are compact: the local index
//! or nlocal plus the remote index.
template <typename T, int KMAX>
struct TruncRow
{
    ParCsr<T const> p;
    int max_elmts;
    T trunc_factor;
    int nlocal_cols;

    //! Calls f(compact column, value) for each entry of row i.
    template <typename F>
    AMREX_GPU_HOST_DEVICE void walk (Long i, F const& f) const
    {
        auto const& c0 = p.csr0;
        for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
            f(c0.col_index[idx], c0.mat[idx]);
        }
        if (has_remote_row(p, i)) {
            auto const& c1 = p.csr1;
            Long const ii = p.row_map[i];
            for (Long idx = c1.row_offset[ii]; idx < c1.row_offset[ii+1]; ++idx) {
                f(nlocal_cols + c1.col_index[idx], c1.mat[idx]);
            }
        }
    }

    //! Drop threshold of row i, with the positive and negative row sums.
    AMREX_GPU_HOST_DEVICE T threshold (Long i, T& sum_pos, T& sum_neg) const
    {
        T amax = T(0);
        sum_pos = T(0);
        sum_neg = T(0);
        walk(i, [&] (int, T v) {
            amax = std::max(amax, std::abs(v));
            if (v > T(0)) { sum_pos += v; } else { sum_neg += v; }
        });
        return trunc_factor * amax;
    }

    //! Number of entries of row i at or above the threshold.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE int count (Long i, T thresh) const
    {
        int n = 0;
        walk(i, [&] (int, T v) { if (std::abs(v) >= thresh) { ++n; } });
        return n;
    }

    //! For max_elmts > 0: returns the number of kept entries and fills
    //! cols/vals by decreasing magnitude.
    AMREX_GPU_HOST_DEVICE
    int select (Long i, int* cols, T* vals, T& sum_pos, T& sum_neg) const
    {
        T const thresh = threshold(i, sum_pos, sum_neg);
        int const kmax = std::min(max_elmts, KMAX);
        int n = 0;
        walk(i, [&] (int col, T v) {
            T const av = std::abs(v);
            if (av < thresh) { return; }
            if (n == kmax && av <= std::abs(vals[n-1])) { return; }
            int pos = (n < kmax) ? n : n-1;
            while (pos > 0 && std::abs(vals[pos-1]) < av) {
                cols[pos] = cols[pos-1];
                vals[pos] = vals[pos-1];
                --pos;
            }
            cols[pos] = col;
            vals[pos] = v;
            if (n < kmax) { ++n; }
        });
        return n;
    }
};

}

namespace detail {

//! CSR made of the rows i with state[i] == c, ordered by cidx (the
//! exclusive count of such rows). `colmap(j, remote)` maps a column index
//! of the local or remote block to the new compact column, or returns -1
//! to drop the entry.
template <typename T, typename F>
typename SpMatrix<T>::local_csr_type
compact_rows (ParCsr<T const> const& p, Long nc, int const* AMREX_RESTRICT state,
              Long const* AMREX_RESTRICT cidx, int c, F const& colmap)
{
    Long const nrows = p.csr0.nrows;

    // Offsets by fine row, then gathered at the selected rows.
    Gpu::DeviceVector<Long> off(nrows+1);
    auto* AMREX_RESTRICT poff = off.data();
    Long const nnz = Scan::PrefixSum<Long>
        (nrows+1,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (i >= nrows || state[i] != c) { return 0; }
             Long n = 0;
             auto const& c0 = p.csr0;
             for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
                 if (colmap(c0.col_index[idx], false) >= 0) { ++n; }
             }
             if (has_remote_row(p, i)) {
                 auto const& c1 = p.csr1;
                 Long const ii = p.row_map[i];
                 for (Long idx = c1.row_offset[ii]; idx < c1.row_offset[ii+1]; ++idx) {
                     if (colmap(c1.col_index[idx], true) >= 0) { ++n; }
                 }
             }
             return n;
         },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { poff[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    typename SpMatrix<T>::local_csr_type csr;
    csr.resize(nc, nnz);
    auto* AMREX_RESTRICT prow = csr.row_offset.data();
    auto* AMREX_RESTRICT pcol = csr.col_index.data();
    auto* AMREX_RESTRICT pmat = csr.mat.data();
    ParallelForOMP(nrows+1, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        if (i == nrows) { prow[nc] = nnz; return; }
        if (state[i] != c) { return; }
        Long q = poff[i];
        prow[cidx[i]] = q;
        auto const& c0 = p.csr0;
        for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
            int const g = colmap(c0.col_index[idx], false);
            if (g >= 0) { pcol[q] = g; pmat[q] = c0.mat[idx]; ++q; }
        }
        if (has_remote_row(p, i)) {
            auto const& c1 = p.csr1;
            Long const ii = p.row_map[i];
            for (Long idx = c1.row_offset[ii]; idx < c1.row_offset[ii+1]; ++idx) {
                int const g = colmap(c1.col_index[idx], true);
                if (g >= 0) { pcol[q] = g; pmat[q] = c1.mat[idx]; ++q; }
            }
        }
    });
    Gpu::streamSynchronize();
    return csr;
}

}

template <typename T>
SpMatrix<T> AMG<T>::interpolation (SpMatrix<T>& A, SpMatrix<TS>& S, SpMatrix<TS> const& ST,
                                   Gpu::DeviceVector<int> const& state,
                                   Gpu::DeviceVector<Long> const& cidx,
                                   AlgPartition const& cpart) const
{
    if (m_interp_type == InterpType::Direct) {
        return build_interp(A, S, state, cidx, cpart);
    } else {
        return build_interp_mm(A, S, ST, state, cidx, cpart,
                               m_interp_type == InterpType::MMExtI);
    }
}

template <typename T>
SpMatrix<typename AMG<T>::TS> AMG<T>::second_pass_strength (SpMatrix<TS>& S,
                                          Gpu::DeviceVector<int> const& state,
                                          Gpu::DeviceVector<Long> const& cidx,
                                          AlgPartition const& cpart)
{
    BL_PROFILE("AMG::second_pass_strength");

    // The C-C block of S^2 is S_CF S_FC after PMIS (S_CC is empty). Both
    // factors are built in compact C and F numberings, so the product never
    // forms the full S^2.
    auto const& s = S.const_parcsr();
    Long const nrows = s.csr0.nrows;
    Long const nc = cpart.numLocalRows();
    int const C = PointType::coarse;
    int const F = PointType::fine;

    Gpu::DeviceVector<Long> fidx;
    auto const fpart = coarse_numbering(state, fidx, PointType::fine);
    Long const nf = fpart.numLocalRows();
    int const nlocal_c = int(nc);
    int const nlocal_f = int(nf);

    auto const* AMREX_RESTRICT ps = state.data();
    auto const* AMREX_RESTRICT pc = cidx.data();
    auto const* AMREX_RESTRICT pf = fidx.data();
    Long const cbegin = cpart.globalRowBegin();
    Long const fbegin = fpart.globalRowBegin();
    Gpu::DeviceVector<Long> gcidx(nrows);
    Gpu::DeviceVector<Long> gfidx(nrows);
    auto* AMREX_RESTRICT pgc = gcidx.data();
    auto* AMREX_RESTRICT pgf = gfidx.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        pgc[i] = pc[i] + cbegin;
        pgf[i] = pf[i] + fbegin;
    });
    auto const state_r = S.gatherRemote(ps);
    auto const gcidx_r = S.gatherRemote(pgc);
    auto const gfidx_r = S.gatherRemote(pgf);
    auto const* AMREX_RESTRICT ps_r = state_r.data();

    // S_CF: C rows, F columns in F numbering (remote: nlocal_f + remote index).
    auto csr_cf = detail::compact_rows<TS>
        (s, nc, ps, pc, C,
         [=] AMREX_GPU_DEVICE (int j, bool remote) -> int {
             if (remote) {
                 return (ps_r[j] == F) ? nlocal_f + j : -1;
             } else {
                 return (ps[j] == F) ? int(pf[j]) : -1;
             }
         });
    SpMatrix<TS> SCF;
    SCF.define_split(cpart, fpart, std::move(csr_cf), nlocal_f, gfidx_r.data(), 0);

    // S_FC: F rows, C columns in C numbering.
    auto csr_fc = detail::compact_rows<TS>
        (s, nf, ps, pf, F,
         [=] AMREX_GPU_DEVICE (int j, bool remote) -> int {
             if (remote) {
                 return (ps_r[j] == C) ? nlocal_c + j : -1;
             } else {
                 return (ps[j] == C) ? int(pc[j]) : -1;
             }
         });
    SpMatrix<TS> SFC;
    SFC.define_split(fpart, cpart, std::move(csr_fc), nlocal_c, gcidx_r.data(), 0);

    auto S2 = amrex::SpGEMM(SCF, SFC, cpart);
    auto& csr = S2.m_csr_local; // AMG is a friend of SpMatrix
    Long const nnz_before = csr.nnz;

    // Drop the diagonal (local block only).
    {
        auto* AMREX_RESTRICT pcol = csr.col_index.data();
        auto* AMREX_RESTRICT prow = csr.row_offset.data();
        auto* AMREX_RESTRICT pmat = csr.mat.data();
        ParallelForOMP(nc, [=] AMREX_GPU_DEVICE (Long r) noexcept
        {
            Long q = prow[r];
            for (Long idx = prow[r]; idx < prow[r+1]; ++idx) {
                if (pcol[idx] != int(r)) {
                    pcol[q] = pcol[idx];
                    pmat[q] = pmat[idx];
                    ++q;
                }
            }
            for (; q < prow[r+1]; ++q) { pcol[q] = -1; }
        });
        Gpu::streamSynchronize();
    }
    // Remove the -1 holes.
    {
        auto const* AMREX_RESTRICT pcol = csr.col_index.data();
        auto const* AMREX_RESTRICT pmat = csr.mat.data();
        auto const* AMREX_RESTRICT prow = csr.row_offset.data();
        Long const nnz0 = csr.nnz;
        typename SpMatrix<TS>::local_csr_type out;
        out.row_offset.resize(nc+1);
        auto* AMREX_RESTRICT orow = out.row_offset.data();
        Long const nnz = Scan::PrefixSum<Long>
            (nc+1,
             [=] AMREX_GPU_DEVICE (Long r) -> Long {
                 if (r >= nc) { return 0; }
                 Long n = 0;
                 for (Long idx = prow[r]; idx < prow[r+1]; ++idx) { n += (pcol[idx] >= 0); }
                 return n;
             },
             [=] AMREX_GPU_DEVICE (Long r, Long const& x) { orow[r] = x; },
             Scan::Type::exclusive, Scan::retSum);
        amrex::ignore_unused(nnz0);
        out.mat.resize(nnz);
        out.col_index.resize(nnz);
        out.nnz = nnz;
        auto* AMREX_RESTRICT ocol = out.col_index.data();
        auto* AMREX_RESTRICT omat = out.mat.data();
        ParallelForOMP(nc, [=] AMREX_GPU_DEVICE (Long r) noexcept
        {
            Long q = orow[r];
            for (Long idx = prow[r]; idx < prow[r+1]; ++idx) {
                if (pcol[idx] >= 0) { ocol[q] = pcol[idx]; omat[q] = pmat[idx]; ++q; }
            }
        });
        Gpu::streamSynchronize();
        csr = std::move(out);
        S2.m_nnz -= nnz_before - csr.nnz;
    }
    return S2;
}

template <typename T>
SpMatrix<T> AMG<T>::select_coarse_rows (SpMatrix<T> const& P,
                                        Gpu::DeviceVector<int> const& state,
                                        Gpu::DeviceVector<Long> const& cidx,
                                        AlgPartition const& cpart)
{
    auto const& p = P.const_parcsr();
    int const nlocal = int(P.columnPartition().numLocalRows());
    auto csr = detail::compact_rows<T>
        (p, cpart.numLocalRows(), state.data(), cidx.data(), int(PointType::coarse),
         [=] AMREX_GPU_DEVICE (int j, bool remote) -> int {
             return remote ? nlocal + j : j;
         });
    SpMatrix<T> Q;
    Q.define_split(cpart, P.columnPartition(), std::move(csr), nlocal, p.col_map, 0);
    return Q;
}

template <typename T>
SpMatrix<T> AMG<T>::truncate_interp (SpMatrix<T>& P, AlgPartition const& cpart,
                                     int max_elmts, T trunc_factor)
{
    BL_PROFILE("AMG::truncate_interp");

    constexpr int KMAX = 32; // rows of P are short
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_elmts <= KMAX,
                                     "AMG: p_max_elmts must be at most 32");

    P.setColumnPartition(cpart);
    auto const& p = P.const_parcsr();
    Long const nrows = p.csr0.nrows;
    int const nlocal = int(cpart.numLocalRows());
    detail::TruncRow<T,KMAX> const tr{p, max_elmts, trunc_factor, nlocal};

    using local_csr_type = typename SpMatrix<T>::local_csr_type;
    local_csr_type csr;
    csr.row_offset.resize(nrows+1);
    auto* AMREX_RESTRICT prow = csr.row_offset.data();
    Long const nnz = Scan::PrefixSum<Long>
        (nrows+1,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (i >= nrows) { return 0; }
             T sp, sn;
             if (max_elmts > 0) {
                 int cols[KMAX];
                 T vals[KMAX];
                 return tr.select(i, cols, vals, sp, sn);
             } else {
                 return tr.count(i, tr.threshold(i, sp, sn));
             }
         },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { prow[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    csr.mat.resize(nnz);
    csr.col_index.resize(nnz);
    csr.nnz = nnz;
    auto* AMREX_RESTRICT pmat = csr.mat.data();
    auto* AMREX_RESTRICT pcol = csr.col_index.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        T sum_pos, sum_neg;
        Long q = prow[i];
        if (max_elmts > 0) {
            int cols[KMAX];
            T vals[KMAX];
            int const n = tr.select(i, cols, vals, sum_pos, sum_neg);
            T kept_pos = T(0), kept_neg = T(0);
            for (int k = 0; k < n; ++k) {
                if (vals[k] > T(0)) { kept_pos += vals[k]; } else { kept_neg += vals[k]; }
            }
            T const fpos = (kept_pos != T(0)) ? sum_pos/kept_pos : T(1);
            T const fneg = (kept_neg != T(0)) ? sum_neg/kept_neg : T(1);
            for (int k = 0; k < n; ++k) {
                pcol[q] = cols[k];
                pmat[q] = vals[k] * ((vals[k] > T(0)) ? fpos : fneg);
                ++q;
            }
        } else { // threshold only, no limit on the count
            T const thresh = tr.threshold(i, sum_pos, sum_neg);
            T kept_pos = T(0), kept_neg = T(0);
            tr.walk(i, [&] (int, T v) {
                if (std::abs(v) >= thresh) {
                    if (v > T(0)) { kept_pos += v; } else { kept_neg += v; }
                }
            });
            T const fpos = (kept_pos != T(0)) ? sum_pos/kept_pos : T(1);
            T const fneg = (kept_neg != T(0)) ? sum_neg/kept_neg : T(1);
            tr.walk(i, [&] (int col, T v) {
                if (std::abs(v) >= thresh) {
                    pcol[q] = col;
                    pmat[q] = v * ((v > T(0)) ? fpos : fneg);
                    ++q;
                }
            });
        }
    });
    Gpu::streamSynchronize();

    // Compact columns sort like the global ones within each block.
    csr.sort();
    SpMatrix<T> Pt;
    Pt.define_split(P.partition(), cpart, std::move(csr), nlocal, p.col_map, 0);
    return Pt;
}

template <typename T>
SpMatrix<T> AMG<T>::build_interp_mm (SpMatrix<T>& A, SpMatrix<TS>& S,
                                     SpMatrix<TS> const& ST,
                                     Gpu::DeviceVector<int> const& state,
                                     Gpu::DeviceVector<Long> const& cidx,
                                     AlgPartition const& cpart, bool plus_i)
{
    BL_PROFILE("AMG::build_interp_mm");

    // Reference: Li, Sjögreen & Yang (2021), Eqs. (3.4)-(3.12). With
    // P = (W; I) in F/C ordering, W = -M_F B_F, where for MM-ext
    //   M_F = (D_FF + D_gamma)^-1 (A^s_FF + D_beta),  B_F = D_beta^-1 A^s_FC,
    // and for MM-ext+i
    //   M_F = (D_FF + D_gamma + D_theta)^-1 (Ahat^s_FF + I),  B_F = A^s_FC.
    // Here both factors keep the natural row numbering; C rows are unit
    // rows, so P = SpGEMM(M, B) in one product.

    auto const& s = S.const_parcsr();
    Long const nrows = s.csr0.nrows;
    auto const* AMREX_RESTRICT ps = state.data();
    auto const* AMREX_RESTRICT pc = cidx.data();
    Long const cbegin = cpart.globalRowBegin();
    int const nlocal_f = int(A.numLocalRows());
    int const nlocal_c = int(cpart.numLocalRows());
    int const C = PointType::coarse;
    int const F = PointType::fine;

    // alpha_i = a_ii + weak off-diagonal row sum = rowsum(A) - rowsum(S),
    // in the units of S. The weights do not depend on the scaling.
    T const scale = strength_scale(A);
    auto const rsA = A.rowSum();
    auto const rsS = S.rowSum();
    auto const* AMREX_RESTRICT prsA = rsA.data();
    auto const* AMREX_RESTRICT prsS = rsS.data();

    Gpu::DeviceVector<Long> gcidx(nrows);
    auto* AMREX_RESTRICT pgc = gcidx.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept { pgc[i] = pc[i] + cbegin; });
    auto const state_r = S.gatherRemote(ps);
    auto const gcidx_r = S.gatherRemote(pgc);
    auto const* AMREX_RESTRICT pstate_r = state_r.data();

    // D_beta: row sum of strong connections to C-points.
    Gpu::DeviceVector<T> dbeta(nrows);
    auto* AMREX_RESTRICT pdb = dbeta.data();
    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        T d = T(0);
        auto const& sc = s.csr0;
        for (Long idx = sc.row_offset[i]; idx < sc.row_offset[i+1]; ++idx) {
            if (ps[sc.col_index[idx]] == C) { d += T(sc.mat[idx]); }
        }
        if (detail::has_remote_row(s, i)) {
            auto const& sc1 = s.csr1;
            Long const ii = s.row_map[i];
            for (Long idx = sc1.row_offset[ii]; idx < sc1.row_offset[ii+1]; ++idx) {
                if (pstate_r[sc1.col_index[idx]] == C) { d += T(sc1.mat[idx]); }
            }
        }
        pdb[i] = d;
    });
    Gpu::streamSynchronize();
    auto const dbeta_r = S.gatherRemote(pdb);

    detail::MMInterpRow<T,TS> const w{s, ST.const_parcsr(), ps, pdb, pc,
                                   pstate_r, dbeta_r.data(), nlocal_f, nlocal_c};

    // Counts: M gets a diagonal plus the usable strong F neighbors; B gets
    // the strong C neighbors. C rows get one unit entry each.
    using local_csr_type = typename SpMatrix<T>::local_csr_type;
    local_csr_type mcsr, bcsr;
    mcsr.row_offset.resize(nrows+1);
    bcsr.row_offset.resize(nrows+1);
    auto* AMREX_RESTRICT pmrow = mcsr.row_offset.data();
    auto* AMREX_RESTRICT pbrow = bcsr.row_offset.data();

    auto usable_f = [=] AMREX_GPU_DEVICE (T aji, int st_j, T db_j) -> bool {
        return st_j == F && (db_j + (plus_i ? aji : T(0))) != T(0);
    };

    Long const mnnz = Scan::PrefixSum<Long>
        (nrows+1,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (i >= nrows) { return 0; }
             Long n = 1;
             if (ps[i] == F) {
                 w.walk(i, [&] (int, int, T, T aji, int st_j, T db_j) {
                     if (usable_f(aji, st_j, db_j)) { ++n; }
                 });
             }
             return n;
         },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { pmrow[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    Long const bnnz = Scan::PrefixSum<Long>
        (nrows+1,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (i >= nrows) { return 0; }
             if (ps[i] == C) { return 1; }
             Long n = 0;
             w.walk(i, [&] (int, int, T, T, int st_j, T) {
                 if (st_j == C) { ++n; }
             });
             return n;
         },
         [=] AMREX_GPU_DEVICE (Long i, Long const& x) { pbrow[i] = x; },
         Scan::Type::exclusive, Scan::retSum);

    mcsr.mat.resize(mnnz);
    mcsr.col_index.resize(mnnz);
    mcsr.nnz = mnnz;
    bcsr.mat.resize(bnnz);
    bcsr.col_index.resize(bnnz);
    bcsr.nnz = bnnz;
    auto* AMREX_RESTRICT pmmat = mcsr.mat.data();
    auto* AMREX_RESTRICT pmcol = mcsr.col_index.data();
    auto* AMREX_RESTRICT pbmat = bcsr.mat.data();
    auto* AMREX_RESTRICT pbcol = bcsr.col_index.data();

    ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        Long mp = pmrow[i];
        Long bp = pbrow[i];
        if (ps[i] != F) {
            pmcol[mp] = int(i);
            pmmat[mp] = T(1);
            pbcol[bp] = int(pc[i]);
            pbmat[bp] = T(1);
            return;
        }

        // Strong F neighbors without a usable denominator are lumped into
        // the diagonal as weak connections (remark after Eq. (3.8)).
        T alpha = prsA[i]/scale - T(prsS[i]);
        T theta = T(0);
        w.walk(i, [&] (int, int, T aij, T aji, int st_j, T db_j) {
            if (st_j == F) {
                T const den = db_j + (plus_i ? aji : T(0));
                if (den != T(0)) {
                    if (plus_i) { theta += aij * aji / den; }
                } else {
                    alpha += aij;
                }
            }
        });
        T const denom = alpha + theta;
        T const scale = (denom != T(0)) ? T(-1)/denom : T(0);
        T const db_i = pdb[i];

        pmcol[mp] = int(i);
        pmmat[mp] = plus_i ? scale : scale*db_i;
        ++mp;
        T const bscale = plus_i ? T(1) : ((db_i != T(0)) ? T(1)/db_i : T(0));
        w.walk(i, [&] (int jf, int jc, T aij, T aji, int st_j, T db_j) {
            if (st_j == C) {
                pbcol[bp] = jc;
                pbmat[bp] = bscale*aij;
                ++bp;
            } else if (usable_f(aji, st_j, db_j)) {
                pmcol[mp] = jf;
                pmmat[mp] = plus_i ? scale*aij/(db_j + aji) : scale*aij;
                ++mp;
            }
        });
    });
    Gpu::streamSynchronize();

    // The diagonal of M is emitted first and must be sorted into place.
    mcsr.sort();
    SpMatrix<T> M, B;
    M.define_split(A.partition(), A.partition(), std::move(mcsr), nlocal_f, s.col_map, 0);
    B.define_split(A.partition(), cpart, std::move(bcsr), nlocal_c, gcidx_r.data(), 0);
    return amrex::SpGEMM(M, B, cpart);
}

}

#endif
