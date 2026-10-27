#include <AMReX.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

// Integer values stored as Real, so that sums are exact in any order.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real value (int i, int j, int k, int n)
{
    return Real((((i*7 + j*13 + k*29 + n*31) % 101) + 101) % 101 - 50);
}

// Compare ParReducePerComp and ParReducePerBoxComp against ParReduce
// restricted to one component (and one box).  Return the number of mismatches.
int check (BoxArray const& ba, int ncomp, int ng)
{
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, ncomp, ng);
    iMultiFab imf(ba, dm, ncomp, ng);
    auto const& ma = mf.arrays();
    auto const& ima = imf.arrays();
    ParallelFor(mf, mf.nGrowVect(), ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
    {
        ma[b](i,j,k,n) = value(i,j,k,n);
        ima[b](i,j,k,n) = (value(i,j,k,n) > Real(45)) ? 1 : 0;
    });
    Gpu::streamSynchronize();

    auto const& cma = mf.const_arrays();
    auto const& cima = imf.const_arrays();
    IntVect const nghost(ng);
    int nbad = 0;

    auto pc = ParReducePerComp(TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum>{},
                               TypeList<Real,Real,Real>{}, mf, nghost, ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
        -> GpuTuple<Real,Real,Real>
    {
        Real v = cma[b](i,j,k,n);
        return {v, v, v};
    });
    auto pcsum = ParReducePerComp(TypeList<ReduceOpSum>{}, TypeList<Long>{}, imf, nghost, ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept -> GpuTuple<Long>
    {
        return { Long(cima[b](i,j,k,n)) };
    });
    auto pclog = ParReducePerComp(TypeList<ReduceOpLogicalOr,ReduceOpLogicalAnd>{},
                                  TypeList<int,int>{}, imf, nghost, ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept -> GpuTuple<int,int>
    {
        int v = cima[b](i,j,k,n);
        return {v, v};
    });
    if (int(pc.size()) != ncomp || int(pcsum.size()) != ncomp || int(pclog.size()) != ncomp) {
        return 1;
    }

    for (int n = 0; n < ncomp; ++n) {
        auto ref = ParReduce(TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum>{},
                             TypeList<Real,Real,Real>{}, mf, nghost,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real,Real,Real>
        {
            Real v = cma[b](i,j,k,n);
            return {v, v, v};
        });
        Long refsum = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Long>{}, imf, nghost,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Long>
        {
            return { Long(cima[b](i,j,k,n)) };
        });
        auto reflog = ParReduce(TypeList<ReduceOpLogicalOr,ReduceOpLogicalAnd>{},
                                TypeList<int,int>{}, imf, nghost,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<int,int>
        {
            int v = cima[b](i,j,k,n);
            return {v, v};
        });
        if (amrex::get<0>(pc[n]) != amrex::get<0>(ref) ||
            amrex::get<1>(pc[n]) != amrex::get<1>(ref) ||
            amrex::get<2>(pc[n]) != amrex::get<2>(ref)) { ++nbad; }
        if (pcsum[n] != refsum) { ++nbad; }
        if (amrex::get<0>(pclog[n]) != amrex::get<0>(reflog) ||
            amrex::get<1>(pclog[n]) != amrex::get<1>(reflog)) { ++nbad; }
    }

    auto pb = ParReducePerBoxComp(TypeList<ReduceOpMin,ReduceOpMax,ReduceOpSum>{},
                                  TypeList<Real,Real,Real>{}, mf, nghost, ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
        -> GpuTuple<Real,Real,Real>
    {
        Real v = cma[b](i,j,k,n);
        return {v, v, v};
    });
    auto pbsum = ParReducePerBoxComp(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mf, nghost, ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept -> GpuTuple<Real>
    {
        return { cma[b](i,j,k,n) };
    });
    const int nlocal = mf.local_size();
    if (int(pb.size()) != nlocal*ncomp || int(pbsum.size()) != nlocal*ncomp) { return 1; }

    for (int L = 0; L < nlocal; ++L) {
        for (int n = 0; n < ncomp; ++n) {
            ReduceOps<ReduceOpMin,ReduceOpMax,ReduceOpSum> rop;
            ReduceData<Real,Real,Real> rdata(rop);
            rop.eval(amrex::grow(mf.box(mf.IndexArray()[L]), nghost), rdata,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real,Real,Real>
            {
                Real v = cma[L](i,j,k,n);
                return {v, v, v};
            });
            auto ref = rdata.value(rop);
            auto const& r = pb[L*ncomp+n];
            if (amrex::get<0>(r) != amrex::get<0>(ref) ||
                amrex::get<1>(r) != amrex::get<1>(ref) ||
                amrex::get<2>(r) != amrex::get<2>(ref) ||
                pbsum[L*ncomp+n] != amrex::get<2>(ref)) { ++nbad; }
        }
    }

    // MultiFab norms over a list of components, in reverse order.
    Vector<int> comps;
    for (int n = ncomp-1; n >= 0; --n) { comps.push_back(n); }
    auto nm0 = mf.norm0(comps, ng);
    for (int m = 0; m < ncomp; ++m) {
        if (nm0[m] != mf.norm0(comps[m], ng)) { ++nbad; }
    }
    if (mf.is_cell_centered()) {
        auto nm1 = mf.norm1(comps, ng);
        for (int m = 0; m < ncomp; ++m) {
            if (nm1[m] != mf.norm1(comps[m], ng)) { ++nbad; }
        }
        if (ng == 0) {
            // Sums of squares can exceed the exact integer range of float.
            auto nm2 = mf.norm2(comps);
            for (int m = 0; m < ncomp; ++m) {
                Real ref = mf.norm2(comps[m]);
                if (std::abs(nm2[m] - ref) > Real(100)*std::numeric_limits<Real>::epsilon()*ref) {
                    ++nbad;
                }
            }
        }
    }
    return nbad;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int nfail = 0, ntest = 0;
        // Uniform boxes, many small boxes, and mixed box sizes; cell and nodal.
        for (int nodal = 0; nodal <= 1; ++nodal) {
        for (int n : {16, 64}) {
        for (IntVect const& mgs : {IntVect(32), IntVect(8), IntVect(AMREX_D_DECL(24,17,13))}) {
            BoxArray ba(Box(IntVect(-3), IntVect(n-4)));
            ba.maxSize(mgs);
            if (nodal) { ba.surroundingNodes(); }
            for (int ncomp : {1, 3}) {
            for (int ng : {0, 1}) {
                ++ntest;
                if (int nbad = check(ba, ncomp, ng); nbad > 0) {
                    ++nfail;
                    amrex::Print() << "FAIL: domain " << ba.minimalBox() << " nboxes "
                                   << ba.size() << " ncomp " << ncomp << " ng " << ng
                                   << " nbad " << nbad << "\n";
                }
            }}
        }}}
        amrex::Print() << ntest << " cases, " << nfail << " failures\n";
        AMREX_ALWAYS_ASSERT(nfail == 0);
    }
    amrex::Finalize();
}
