#include <AMReX.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

// Return the number of points visited a wrong number of times.
int check (BoxArray const& ba, IntVect const& stride, IntVect const& offset)
{
    DistributionMapping dm(ba);
    iMultiFab cnt(ba, dm, 1, 1);
    cnt.setVal(0);
    auto const& ma = cnt.arrays();

    ParallelForStrided(cnt, stride, offset,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        Gpu::Atomic::AddNoRet(&ma[b](i,j,k), 1);
    });

    Vector<Box> hvb;
    for (MFIter mfi(cnt); mfi.isValid(); ++mfi) { hvb.push_back(mfi.validbox()); }
    Gpu::DeviceVector<Box> dvb(hvb.size());
    Gpu::copyAsync(Gpu::hostToDevice, hvb.begin(), hvb.end(), dvb.begin());
    Box const* pvb = dvb.data();

    int nbad = ParReduce(TypeList<ReduceOpSum>{}, TypeList<int>{}, cnt, IntVect(1),
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<int>
    {
        IntVect iv(AMREX_D_DECL(i,j,k));
        bool match = pvb[b].contains(iv);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            match = match && ((iv[d]-offset[d]) % stride[d] + stride[d]) % stride[d] == 0;
        }
        return { int(ma[b](i,j,k) != (match ? 1 : 0)) };
    });
    ParallelDescriptor::ReduceIntSum(nbad);
    return nbad;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int nfail = 0, ntest = 0;
        for (int nodal = 0; nodal <= 1; ++nodal) {
        for (int lo : {-7, -4, 0, 3}) {
        for (int n : {1, 5, 16, 33}) {
        for (int mgs : {4, 7, 16}) {
            BoxArray ba(Box(IntVect(lo), IntVect(lo+n-1)));
            ba.maxSize(mgs);
            if (nodal) { ba.surroundingNodes(); }
            for (int s : {1, 2, 3}) {
            for (int o : {-1, 0, 1, 2}) {
                IntVect stride(s), offset(o);
#if (AMREX_SPACEDIM >= 2)
                stride[1] = (s == 3) ? 2 : s;
                offset[1] = o + 1;
#endif
                ++ntest;
                if (int nbad = check(ba, stride, offset); nbad > 0) {
                    ++nfail;
                    amrex::Print() << "FAIL: domain " << ba.minimalBox() << " nboxes "
                                   << ba.size() << " stride " << stride << " offset "
                                   << offset << " nbad " << nbad << "\n";
                }
            }}
        }}}}
        amrex::Print() << ntest << " cases, " << nfail << " failures\n";
        AMREX_ALWAYS_ASSERT(nfail == 0);
    }
    amrex::Finalize();
}
