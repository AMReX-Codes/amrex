#include <AMReX.H>
#include <AMReX_GpuParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

//
// Unit test for the Gpu::DeviceVector overloads of the MPI collectives in
// AMReX_GpuParallelReduce.H:
//   - ParallelAllReduce::Sum(Gpu::DeviceVector&, comm)
//   - ParallelReduce::Sum   (Gpu::DeviceVector&, root, comm)
//   - ParallelDescriptor::Bcast(Gpu::DeviceVector&, root, comm)
//
// Each rank seeds a device vector with rank-dependent values. The test checks
// the reduced/broadcast results against the analytic expectation. It is
// meaningful with >1 MPI rank but also passes serially, and exercises the
// device<->host staging path when built for GPU without GPU-aware MPI.
//

using namespace amrex;

namespace {

Gpu::DeviceVector<Real> to_device (Vector<Real> const& h)
{
    Gpu::DeviceVector<Real> d(h.size());
    Gpu::copy(Gpu::hostToDevice, h.begin(), h.end(), d.begin());
    return d;
}

Vector<Real> to_host (Gpu::DeviceVector<Real> const& d)
{
    Vector<Real> h(d.size());
    Gpu::copy(Gpu::deviceToHost, d.begin(), d.end(), h.begin());
    return h;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        const int nprocs = ParallelDescriptor::NProcs();
        const int myproc = ParallelDescriptor::MyProc();
        const int ioproc = ParallelDescriptor::IOProcessorNumber();
        MPI_Comm comm = ParallelDescriptor::Communicator();

        const int n = 10;

        // rank p seeds entry i with (p + i), so the sum over all ranks is
        //   sum_{p=0}^{nprocs-1} (p + i)
        auto expected_sum = [=] (int i) -> Real {
            Real s = Real(0);
            for (int p = 0; p < nprocs; ++p) { s += Real(p + i); }
            return s;
        };

        // ParallelAllReduce::Sum -> result valid on every rank
        {
            Vector<Real> h(n);
            for (int i = 0; i < n; ++i) { h[i] = Real(myproc + i); }
            auto d = to_device(h);

            ParallelAllReduce::Sum(d, comm);

            auto const r = to_host(d);
            for (int i = 0; i < n; ++i) {
                AMREX_ALWAYS_ASSERT(r[i] == expected_sum(i));
            }
        }

        // ParallelReduce::Sum -> result valid on root only
        {
            Vector<Real> h(n);
            for (int i = 0; i < n; ++i) { h[i] = Real(myproc + i); }
            auto d = to_device(h);

            ParallelReduce::Sum(d, ioproc, comm);

            if (myproc == ioproc) {
                auto const r = to_host(d);
                for (int i = 0; i < n; ++i) {
                    AMREX_ALWAYS_ASSERT(r[i] == expected_sum(i));
                }
            }
        }

        // ParallelDescriptor::Bcast -> root's data to every rank.
        // Contract: every rank pre-allocates the receiver to the root's length.
        {
            Gpu::DeviceVector<Real> d(n, Real(0));
            if (myproc == ioproc) {
                Vector<Real> h(n);
                for (int i = 0; i < n; ++i) { h[i] = Real(100 + i); }
                Gpu::copy(Gpu::hostToDevice, h.begin(), h.end(), d.begin());
            }

            ParallelDescriptor::Bcast(d, ioproc, comm);

            auto const r = to_host(d);
            for (int i = 0; i < n; ++i) {
                AMREX_ALWAYS_ASSERT(r[i] == Real(100 + i));
            }
        }

        amrex::Print() << "GpuParallelReduce: all tests passed on " << nprocs
                       << " rank(s).\n";
    }
    amrex::Finalize();
}
