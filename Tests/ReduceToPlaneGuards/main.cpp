#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <cmath>
#include <vector>

using namespace amrex;

// Coverage for the guard-cell summing option of ReduceToPlaneMF2:
//
//   ReduceToPlaneMF2<Op>(dir, domain, mf, f, nghost, period)
//
// reduces `mf` to a plane normal to `dir`, summing `nghost` guard cells in the
// reduction direction into the plane and folding the transverse guard cells
// into the owning valid nodes (via ParallelAdd with `period`). On *raw*
// (un-synchronized) data this conserves a summed quantity (e.g. a deposited
// charge) regardless of the domain decomposition -- in particular it recovers
// the charge that higher-order deposition leaves in the guard cells at internal
// box seams and outside the single valid cell along `dir`.
//
// The test "deposits" a unit source straddling an internal box seam and the
// single valid cell in `dir`, so part of it lands in guard cells, and checks:
//   (1) ReduceToPlaneMF2 with nghost recovers the full source (conservation),
//   (2) the result is independent of the number of boxes (MPI-safe),
//   (3) the default (nghost = 0) drops the guard-resident part.

#if (AMREX_SPACEDIM >= 2)

namespace {

constexpr int n = 8;            // cells per transverse direction
constexpr int dir = AMREX_SPACEDIM - 1; // reduce over the last dimension

// the single source sits on the dim-0 box seam (at n/2), interior in the other
// transverse directions, and in the only valid cell along `dir`
IntVect source_cell ()
{
    IntVect s(0);
    s[0] = n/2;
    for (int d = 1; d < AMREX_SPACEDIM - 1; ++d) { s[d] = n/4; }
    return s; // s[dir] == 0
}

// build a cell-centered domain with a single cell along `dir`
Box make_domain ()
{
    IntVect hi(n-1);
    hi[dir] = 0;
    return Box(IntVect(0), hi);
}

// deposit a unit source as a 3^SPACEDIM stencil of equal weights; each source
// is added to the *owning* box only (the box whose valid region holds the
// central cell), so the deposition is unique even when the stencil spills into
// guard cells -- exactly like particle charge deposition.
void deposit_unit_source (MultiFab& mf, IntVect const& s)
{
    Real const w = Real(1.0) / Real(std::pow(3, AMREX_SPACEDIM));
    Box const stencil(s - IntVect(1), s + IntVect(1));
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(s)) {
            Box const dep = stencil & mfi.fabbox();
            auto const& a = mf.array(mfi);
            amrex::ParallelFor(dep, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i,j,k) += w;
            });
        }
    }
}

// total of the reduced-plane result, summing each unique cell once across ranks
Real reduced_total (MultiFab const& mf, Box const& domain, IntVect const& nghost)
{
    auto const& ma = mf.const_arrays();
    auto [plane_local, plane_unique] = ReduceToPlaneMF2<ReduceOpSum>(
        dir, domain, mf,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> Real
        {
            return ma[b](i,j,k);
        },
        nghost);
    amrex::ignore_unused(plane_local);
    return plane_unique.sum_unique(0, false);
}

// run the deposit + reduction for a given maximum grid size and return both the
// guard-summing total (nghost = nGrow) and the default total (nghost = 0)
std::pair<Real,Real> run (IntVect const& max_grid_size)
{
    Box const domain = make_domain();
    IntVect const ng(1);

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    MultiFab mf(ba, dm, 1, ng);
    mf.setVal(0.0);
    deposit_unit_source(mf, source_cell());

    Real const with_guards = reduced_total(mf, domain, ng);
    Real const without_guards = reduced_total(mf, domain, IntVect(0));
    return {with_guards, without_guards};
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        Real const q_in = 1.0; // a single source, fully inside the domain

        // multi-box: split dim 0 to create an internal seam under the source
        IntVect mgs_multi(n);
        mgs_multi[0] = n/2;
        mgs_multi[dir] = 1;
        auto const [multi_guards, multi_default] = run(mgs_multi);

        // single box: no internal seams
        IntVect mgs_single(1024);
        mgs_single[dir] = 1;
        auto const [single_guards, single_default] = run(mgs_single);

        // (1) conservation: the guard-summing reduction recovers the full source
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(multi_guards, q_in));
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(single_guards, q_in));

        // (2) decomposition independence (MPI-safe): same answer with seams
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(multi_guards, single_guards));

        // (3) the default (no guards) drops the guard-resident charge: at least
        //     the cells along `dir` outside the single valid cell are lost
        AMREX_ALWAYS_ASSERT(multi_default < q_in);
        AMREX_ALWAYS_ASSERT(single_default < q_in);

        amrex::Print() << "ReduceToPlaneGuards: PASSED\n"
                       << "  conserved (multi/single) = " << multi_guards
                       << " / " << single_guards << " (expected " << q_in << ")\n"
                       << "  default  (multi/single)  = " << multi_default
                       << " / " << single_default << " (guard charge dropped)\n";
    }
    amrex::Finalize();
    return 0;
}

#else // AMREX_SPACEDIM == 1

int main (int argc, char* argv[])
{
    // The seam/guard scenario needs at least one transverse direction.
    amrex::Initialize(argc, argv);
    amrex::Print() << "ReduceToPlaneGuards: skipped in 1D\n";
    amrex::Finalize();
    return 0;
}

#endif
