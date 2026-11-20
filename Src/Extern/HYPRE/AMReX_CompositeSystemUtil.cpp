#include <AMReX_CompositeSystemUtil.H>
#include <AMReX_Hypre.H>
#include <AMReX_Reduce.H>

#include <_hypre_parcsr_mv.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace amrex::composite {

bool ij_entry_slots (HYPRE_IJMatrix A, Gpu::DeviceVector<HYPRE_BigInt> const& rows,
                     Gpu::DeviceVector<HYPRE_BigInt> const& cols, MPI_Comm comm,
                     EntryGather& g)
{
    BL_PROFILE("composite::ij_entry_slots");
    g.clear();
    HYPRE_ParCSRMatrix parA = nullptr;
    HYPRE_IJMatrixGetObject(A, (void**) &parA);
    auto* pA = (hypre_ParCSRMatrix*) parA;
    hypre_CSRMatrix* diag = hypre_ParCSRMatrixDiag(pA);
    hypre_CSRMatrix* offd = hypre_ParCSRMatrixOffd(pA);
    const auto nr = Long(hypre_CSRMatrixNumRows(diag));
    auto host_copy = [] (hypre_CSRMatrix* m, HYPRE_Int const* p, Long n) {
        Vector<HYPRE_Int> h(n);
        if (n > 0) {
            hypre_TMemcpy(h.data(), p, HYPRE_Int, n, HYPRE_MEMORY_HOST,
                          hypre_CSRMatrixMemoryLocation(m));
        }
        return h;
    };
    auto const di = host_copy(diag, hypre_CSRMatrixI(diag), nr+1);
    auto const dj = host_copy(diag, hypre_CSRMatrixJ(diag), di[nr]);
    auto const oi = host_copy(offd, hypre_CSRMatrixI(offd), nr+1);
    auto const oj = host_copy(offd, hypre_CSRMatrixJ(offd), oi[nr]);
    HYPRE_BigInt const* colmap = hypre_ParCSRMatrixColMapOffd(pA); // host
    const auto ncolmap = Long(hypre_CSRMatrixNumCols(offd));
    const HYPRE_BigInt row0 = hypre_ParCSRMatrixFirstRowIndex(pA);
    const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(pA);
    const HYPRE_BigInt col1 = hypre_ParCSRMatrixLastColDiag(pA);

    const auto ntot = Long(rows.size());
    Vector<HYPRE_BigInt> hr(ntot), hc(ntot);
    Gpu::copyAsync(Gpu::deviceToHost, rows.begin(), rows.end(), hr.begin());
    Gpu::copyAsync(Gpu::deviceToHost, cols.begin(), cols.end(), hc.begin());
    Gpu::streamSynchronize();
    Vector<Long> entry(ntot);
    bool all_found = true;
    for (Long n = 0; n < ntot && all_found; ++n) {
        const Long r = Long(hr[n] - row0);
        Long s = -1;
        bool found = false;
        if (r < 0 || r >= nr) {
            // Only rows owned by this rank are assembled here
            all_found = false;
            break;
        }
        if (hc[n] >= col0 && hc[n] <= col1) {
            const auto lc = HYPRE_Int(hc[n] - col0);
            for (HYPRE_Int k = di[r]; k < di[r+1]; ++k) {
                if (dj[k] == lc) { s = k; found = true; break; }
            }
        } else {
            auto const* it = std::lower_bound(colmap, colmap+ncolmap, hc[n]);
            if (it != colmap+ncolmap && *it == hc[n]) {
                const auto lc = HYPRE_Int(it - colmap);
                for (HYPRE_Int k = oi[r]; k < oi[r+1]; ++k) {
                    if (oj[k] == lc) { s = Long(di[nr]) + k; found = true; break; }
                }
            }
        }
        all_found = found;
        entry[n] = s;
    }
    // hypre may have dropped an entry: then all ranks keep the general path,
    // whose assembly is collective
    ParallelAllReduce::And(all_found, comm);
    if (!all_found) { return false; }
    define_gather(g, entry, Long(di[nr]) + Long(oi[nr]));
    return true;
}

void ij_update_values (HYPRE_IJMatrix A, EntryGather const& g, Gpu::DeviceVector<Real> const& vals)
{
    BL_PROFILE("composite::ij_update_values");
    HYPRE_ParCSRMatrix parA = nullptr;
    HYPRE_IJMatrixGetObject(A, (void**) &parA);
    auto* pA = (hypre_ParCSRMatrix*) parA;
    hypre_CSRMatrix* diag = hypre_ParCSRMatrixDiag(pA);
    hypre_CSRMatrix* offd = hypre_ParCSRMatrixOffd(pA);
    Real* dd = hypre_CSRMatrixData(diag);
    Real* od = hypre_CSRMatrixData(offd);
    const auto nd = Long(hypre_CSRMatrixNumNonzeros(diag));
    const auto no = Long(hypre_CSRMatrixNumNonzeros(offd));
    AMREX_ALWAYS_ASSERT(Long(g.offset.size()) == nd+no+1);
    // hypre may still be using the matrix on its own stream
    Gpu::hypreSynchronize();
    auto const* go = g.offset.data();
    auto const* gt = g.trip.data();
    auto const* v = vals.data();
    amrex::ParallelFor(nd+no, [=] AMREX_GPU_DEVICE (Long e) noexcept
    {
        Real sum = 0.0;
        for (Long k = go[e]; k < go[e+1]; ++k) { sum += v[gt[k]]; }
        if (e < nd) { dd[e] = sum; } else { od[e-nd] = sum; }
    });
    Gpu::streamSynchronize(); // before hypre uses the matrix
}

ReusedSolve ij_solve_reused (HypreIJIface& ij, bool reused, bool throw_exception,
                             Real reltol, Real abstol, int maxiter, MPI_Comm comm,
                             int verbose, std::string const& prefix,
                             std::function<bool(Real)> const& converged,
                             std::function<void()> const& restore_x)
{
    ReusedSolve r;
    // With a reused setup, a hypre error or a non-finite residual is caught
    // and retried below like a solve that does not converge.  hypre's error
    // flag is per rank, so the decision is made collective: the solve itself
    // ran on every rank, and the iteration count comes from a rank that
    // reports one (HypreIJIface::solve counts nothing after an error)
    ij.setThrowException(throw_exception || reused);
    int blew_up = 0;
    try {
        ij.solve(reltol, abstol, maxiter);
    } catch (std::runtime_error const& e) {
        if (!reused) { // the caller's policy: it throws after its bookkeeping
            ij.setThrowException(throw_exception);
            r.error = e.what();
            return r;
        }
        blew_up = 1;
    }
    ij.setThrowException(throw_exception);
    r.num_iters = static_cast<int>(ij.getNumIters());
    r.final_res = ij.getFinalResidualNorm();
    if (reused) {
        ParallelAllReduce::Max(blew_up, comm);
        ParallelAllReduce::Max(r.num_iters, comm);
    }
    if (reused && (blew_up || !converged(r.final_res))) {
        // The reused setup may have become too poor a preconditioner for the
        // current matrix: redo it and continue from the current x (from the
        // initial guess after a blow-up)
        if (verbose) {
            amrex::Print() << prefix << ": no convergence with the reused setup ("
                           << r.num_iters << " iterations); redoing the setup\n";
        }
        if (blew_up || !std::isfinite(r.final_res)) { restore_x(); }
        ij.setNeedSetup();
        r.redone = true; // also if the solve throws: the setup is done first
        try {
            ij.solve(reltol, abstol, maxiter);
        } catch (std::runtime_error const& e) {
            r.error = e.what();
        }
        r.num_iters += static_cast<int>(ij.getNumIters());
        r.final_res = ij.getFinalResidualNorm();
    }
    return r;
}

ReusedSolve algmg_solve_kept (AlgMG<Real>& mg, SpMatrix<Real>& mat, bool kept,
                              bool throw_exception, AlgVector<Real>& x,
                              AlgVector<Real> const& b, Real bnorm, MPI_Comm comm,
                              int verbose, std::string const& prefix,
                              std::function<bool(Real)> const& converged,
                              std::function<void()> const& restore_x)
{
    ReusedSolve r;
    auto residual = [&] () {
        return (bnorm > Real(0.0)) ? mg.getResidualNorm() / bnorm : Real(0.0);
    };
    // With a kept hierarchy, a blow-up is caught and retried below like a
    // solve that does not converge (AlgMG decides it on global norms; the
    // reduction keeps every rank on the same path regardless)
    mg.setThrowException(throw_exception || kept);
    int blew_up = 0;
    try {
        mg.solve(x, b);
    } catch (std::runtime_error const& e) {
        if (!kept) { // the caller's policy: it throws after its bookkeeping
            mg.setThrowException(throw_exception);
            r.error = e.what();
            return r;
        }
        blew_up = 1;
    }
    mg.setThrowException(throw_exception);
    r.num_iters = mg.getNumIters();
    r.final_res = residual();
    if (kept) { ParallelAllReduce::Max(blew_up, comm); }
    if (kept && (blew_up || !converged(r.final_res))) {
        // The kept hierarchy may be too poor for the current matrix: build a
        // new one and continue from the current x (from the initial guess
        // after a blow-up)
        if (verbose) {
            amrex::Print() << prefix << ": no convergence with the kept AlgMG hierarchy ("
                           << r.num_iters << " iterations); rebuilding it\n";
        }
        if (blew_up) { restore_x(); }
        mg.define(mat);
        r.redone = true; // also if the solve throws: the hierarchy is built first
        try {
            mg.solve(x, b);
        } catch (std::runtime_error const& e) {
            r.error = e.what();
        }
        r.num_iters += mg.getNumIters();
        r.final_res = residual();
    }
    return r;
}

}
