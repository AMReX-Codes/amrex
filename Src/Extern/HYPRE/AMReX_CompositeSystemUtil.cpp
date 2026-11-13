#include <AMReX_CompositeSystemUtil.H>
#include <AMReX_Hypre.H>

#include <_hypre_parcsr_mv.h>

#include <algorithm>

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

}
