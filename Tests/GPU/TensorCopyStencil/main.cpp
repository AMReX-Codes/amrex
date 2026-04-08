#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Print.H>
#include <AMReX_StencilFor.H>

#include <cmath>

using namespace amrex;

namespace {

[[nodiscard]] Box make_domain (int nx, int ny = 1, int nz = 1) noexcept
{
    return Box(IntVect(0),
               IntVect(AMREX_D_DECL(nx-1, ny-1, nz-1)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real source_value (int i, int j, int k) noexcept
{
    return Real(0.125)*(i*i + 3*j + 5*k*k)
         + Real(0.5)*(2*i - j + 4*k)
         + Real(1.0);
}

template <typename Arr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real stencil_sum (Arr const& arr, IntVect const& halo, int i, int j, int k) noexcept
{
    Dim3 const h = halo.dim3(0);
    Real sum = Real(0.0);

    for (int kk = -h.z; kk <= h.z; ++kk) {
        for (int jj = -h.y; jj <= h.y; ++jj) {
            for (int ii = -h.x; ii <= h.x; ++ii) {
                int const ai = (ii < 0) ? -ii : ii;
                int const aj = (jj < 0) ? -jj : jj;
                int const ak = (kk < 0) ? -kk : kk;
                Real const weight = Real(1.0)
                    / (Real(1.0) + Real(ai) + Real(2*aj) + Real(3*ak));
                sum += weight * arr(i+ii, j+jj, k+kk);
            }
        }
    }

    return sum;
}

template <int MT = 128>
void check_case (Box const& bx,
                 IntVect const& tile,
                 IntVect const& halo,
                 IntVect const& shared_padding,
                 Gpu::TensorCopyPolicy policy,
                 char const* label)
{
    Box const src_box = amrex::grow(bx, halo);
    FArrayBox src(src_box, 1);
    FArrayBox baseline(bx, 1);
    FArrayBox result(bx, 1);

    auto const src_arr = src.array();
    ParallelFor(src_box,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        src_arr(i, j, k) = source_value(i, j, k);
    });

    auto const src_const = src.const_array();
    auto const baseline_arr = baseline.array();
    ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        baseline_arr(i, j, k) = stencil_sum(src_const, halo, i, j, k);
    });

    auto info = Gpu::StencilInfo{}
        .setTile(tile)
        .setHalo(halo)
        .setSharedPadding(shared_padding)
        .setTensorCopyPolicy(policy);

    auto const result_arr = result.array();
    StencilFor<MT>(bx, info, src.const_array(),
    [=] AMREX_GPU_DEVICE (auto const& tile_view, int i, int j, int k) noexcept
    {
        result_arr(i, j, k) = stencil_sum(tile_view, halo, i, j, k);
    });

    Gpu::streamSynchronize();

    Long const npts = bx.numPts();
    Vector<Real> h_baseline(npts);
    Vector<Real> h_result(npts);
    Gpu::copyAsync(Gpu::deviceToHost, baseline.dataPtr(), baseline.dataPtr()+npts, h_baseline.begin());
    Gpu::copyAsync(Gpu::deviceToHost, result.dataPtr(), result.dataPtr()+npts, h_result.begin());
    Gpu::streamSynchronize();

    for (Long i = 0; i < npts; ++i) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(h_baseline[i] - h_result[i]) < Real(1.e-12),
                                         label);
    }

    amrex::Print() << "TensorCopyStencil: " << label << " passed\n";
}

void run_tests ()
{
    check_case(make_domain(AMREX_D_DECL(24, 10, 8)),
               IntVect(AMREX_D_DECL(8, 5, 4)),
               IntVect(1),
               IntVect(0),
               Gpu::TensorCopyPolicy::Auto,
               "uniform-radius1-auto");

    if (Gpu::supportsTensorCopy()) {
        check_case(make_domain(AMREX_D_DECL(24, 10, 8)),
                   IntVect(AMREX_D_DECL(8, 5, 4)),
                   IntVect(1),
                   IntVect(0),
                   Gpu::TensorCopyPolicy::Always,
                   "uniform-radius1-always");
    }

    check_case(make_domain(AMREX_D_DECL(24, 10, 8)),
               IntVect(AMREX_D_DECL(8, 5, 4)),
               IntVect(2),
               IntVect(0),
               Gpu::TensorCopyPolicy::Never,
               "uniform-radius2-never");

    check_case(make_domain(AMREX_D_DECL(29, 19, 13)),
               IntVect(AMREX_D_DECL(7, 6, 4)),
               IntVect(4),
               IntVect(0),
               Gpu::TensorCopyPolicy::Auto,
               "partial-radius4-auto");

    check_case<100>(make_domain(AMREX_D_DECL(24, 10, 8)),
                    IntVect(AMREX_D_DECL(8, 5, 4)),
                    IntVect(1),
                    IntVect(0),
                    Gpu::TensorCopyPolicy::Auto,
                    "uniform-radius1-auto-mt100");

    check_case(make_domain(AMREX_D_DECL(24, 10, 8)),
               IntVect(AMREX_D_DECL(8, 5, 4)),
               IntVect(1),
               IntVect(AMREX_D_DECL(1, 0, 0)),
               Gpu::TensorCopyPolicy::Auto,
               "uniform-radius1-auto-paddingx1");

    check_case(make_domain(AMREX_D_DECL(29, 19, 13)),
               IntVect(AMREX_D_DECL(7, 6, 4)),
               IntVect(AMREX_D_DECL(5, 2, 1)),
               IntVect(0),
               Gpu::TensorCopyPolicy::Never,
               "partial-anisotropic-never");
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        run_tests();
    }
    amrex::Finalize();
    return 0;
}
