#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Morton.H>
#include <AMReX_Print.H>

#include <cmath>
#include <cstdint>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // Domain [-1,1): (xhi - xmin)/(xmax - xmin) rounds to exactly 1.
        const Real xmin = Real(-1.0);
        const Real xmax = Real( 1.0);
        const Real xmid = Real( 0.5);
        const Real xhi  = std::nextafter(xmax, xmin);
        AMREX_ALWAYS_ASSERT(xhi < xmax);
        amrex::Print() << "(xhi-xmin)/(xmax-xmin) == 1: "
                       << ((xhi-xmin)/(xmax-xmin) == Real(1.0)) << '\n';

        // Each encoder must map xhi to the top of its range, not wrap to 0.
        AMREX_ALWAYS_ASSERT(Morton::toUInt10(xmin, xmin, xmax) == 0U);
        AMREX_ALWAYS_ASSERT(Morton::toUInt16(xmin, xmin, xmax) == 0U);
        AMREX_ALWAYS_ASSERT(Morton::toUInt32(xmin, xmin, xmax) == 0U);
        AMREX_ALWAYS_ASSERT(Morton::toUInt10(xhi, xmin, xmax) == (1U << 10) - 1U);
        AMREX_ALWAYS_ASSERT(Morton::toUInt16(xhi, xmin, xmax) == (1U << 16) - 1U);
        AMREX_ALWAYS_ASSERT(Morton::toUInt32(xhi, xmin, xmax) >  Morton::toUInt32(xmid, xmin, xmax));
        AMREX_ALWAYS_ASSERT(Morton::toUInt32(xhi, xmin, xmax) >= 0xFFFFFE00U); // SP bound

        // Codes along the diagonal must be non-decreasing and end at the max.
        GpuArray<Real,AMREX_SPACEDIM> plo, phi;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { plo[d] = xmin; phi[d] = xmax; }
        const int n = 1000;
        Gpu::DeviceVector<std::uint32_t> dcode(n);
        auto* pc = dcode.data();
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            Real x = (i == n-1) ? xhi : xmin + (xmax-xmin) * Real(i) / Real(n);
            pc[i] = Morton::get32BitCode(AMREX_D_DECL(x,x,x), plo, phi);
        });
        Gpu::PinnedVector<std::uint32_t> hcode(n);
        Gpu::copy(Gpu::deviceToHost, dcode.begin(), dcode.end(), hcode.begin());
        for (int i = 1; i < n; ++i) {
            AMREX_ALWAYS_ASSERT(hcode[i] >= hcode[i-1]);
        }
        AMREX_ALWAYS_ASSERT(hcode[0] == 0U);
        AMREX_ALWAYS_ASSERT(hcode[n-1] > hcode[n/2]);
        AMREX_ALWAYS_ASSERT(hcode[n-1] == Morton::get32BitCode(AMREX_D_DECL(xhi,xhi,xhi), plo, phi));
        amrex::Print() << "Morton test passed\n";
    }
    amrex::Finalize();
}
